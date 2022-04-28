
#ifndef INVERTER_FULL_C_
#define INVERTER_FULL_C_

#include "../Include/common_defines.h"
#include "./struct_c_def.h"
#include "./fermionic_utilities.h"
#include "./fermion_matrix.h"
#include "./inverter_full.h"
#include "../Include/inverter_tricks.h"
#include "../Mpi/multidev.h"


#ifndef __GNUC__
#include "openacc.h"
#endif

#define SAFETY_MARGIN 0.95

int ker_invert_openacc(__restrict const su3_soa * u, // non viene aggiornata mai qui dentro
        ferm_param *pars,
        __restrict vec3_soa * solution,
        __restrict const vec3_soa * in, // non viene aggiornato mai qui dentro
        double res,
        __restrict vec3_soa * loc_r,
        __restrict vec3_soa * loc_h,
        __restrict vec3_soa * loc_s,
        __restrict vec3_soa * loc_p,
        const int  max_cg,
        double shift,
        int * cg_return )
{

    int cg = 0;
    long int i;
    double delta, alpha, lambda, omega, gammag;

    double source_norm = l2norm2_global(in);
    if(verbosity_lv>5){
        printf("source nor_:%f\n",source_norm);
        printf("MPI%02d - ker_invert_openacc - pointers: u=%p,sol=%p,in=%p\n",devinfo.myrank,
                u,solution,in);
    }



    if(verbosity_lv > 4 ){
        printf("MPI%02d: max_cg  %d, restarting every %d, shift %f\n", devinfo.myrank,
                max_cg,inverter_tricks.restartingEvery,shift);
#pragma acc update host(u[0:1])
        printf("MPI%02d: u[sizeh/2]: %f\n", devinfo.myrank, creal(u[0].r0.c0[sizeh/2]));
#pragma acc update host(solution[0:1])
        printf("MPI%02d: solution[sizeh/2]: %f\n", devinfo.myrank, creal(solution[0].c0[sizeh/2]));
#pragma acc update host(in[0:1])
        printf("MPI%02d: in[sizeh/2]: %f\n", devinfo.myrank, creal(in[0].c0[sizeh/2]));
    }

    do {

        fermion_matrix_multiplication_shifted(u,loc_s,solution,loc_h,pars,shift);

        combine_in1_minus_in2(in,loc_s,loc_r);
 
        assign_in_to_out(loc_r,loc_p);


        delta=l2norm2_global(loc_r);
        if(verbosity_lv>5) printf("Delta:%f",delta);

        // loop over cg iterations
        if (verbosity_lv > 3 && 0==devinfo.myrank )
            printf("STARTING CG:\nCG\tR\n");

        int cg_restarted = 0 ;
        do {
            cg++;    
            cg_restarted++;
            // s=(M^dag M)p    alpha=(p,s)

            fermion_matrix_multiplication_shifted(u,loc_s,loc_p,loc_h,pars,shift);
            alpha = real_scal_prod_global(loc_p,loc_s);
            if(verbosity_lv>5) printf("alpha:%e\n",alpha);

            omega=delta/alpha;     
            if(verbosity_lv>5) printf("omega:%e\n",omega);
            // solution+=omega*p  r-=omega*s
            // lambda=(r,r);

            combine_in1xfactor_plus_in2(loc_p,omega,solution,solution);
            combine_in1xfactor_plus_in2(loc_s,-omega,loc_r,loc_r);

            lambda = l2norm2_global(loc_r);
            if(verbosity_lv>5) printf("lambda:%e\n",lambda);
            gammag=lambda/delta;
            delta=lambda;
            if(verbosity_lv>5) printf("Delta:%e\n",delta);
            // p=r+gammag*p
            combine_in1xfactor_plus_in2(loc_p,gammag,loc_r,loc_p);

            if (verbosity_lv > 3 && cg%100==0 && 0==devinfo.myrank  ){

                printf("%d\t%1.1e\n",cg, sqrt(lambda/source_norm)/res);fflush(stdout);

            }

        } while( (sqrt(lambda/source_norm)>res*SAFETY_MARGIN) && 
                cg_restarted<inverter_tricks.restartingEvery);
        if(0==devinfo.myrank && verbosity_lv >4)
            printf("Exited inner cycle in inverter\n");

    } while( (sqrt(lambda/source_norm)>res) && cg<max_cg);


    if (verbosity_lv > 3  && 0==devinfo.myrank ) printf("\n");

    fermion_matrix_multiplication_shifted(u,loc_s,solution,loc_h,pars,shift);
    combine_in1_minus_in2(in,loc_s,loc_h); // r = s - y  
    double  giustoono=l2norm2_global(loc_h)/source_norm;
    if(verbosity_lv > 1 && 0==devinfo.myrank  ){
        printf("Terminated invert after   %d    iterations", cg);
        printf("[res/stop_res=  %e , stop_res=%e ]\n",
                sqrt(giustoono)/res,res);
    }
    if(cg==max_cg  && 0==devinfo.myrank )
    {
        printf("WARNING: maximum number of iterations reached in invert\n");
    }

    *cg_return = cg;
    if (sqrt(giustoono) <= res)
      return INVERTER_SUCCESS;
    else return INVERTER_FAILURE;

}


#endif

