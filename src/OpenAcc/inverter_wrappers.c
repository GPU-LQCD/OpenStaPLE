#ifndef INVERTER_WRAPPERS_C_
#define INVERTER_WRAPPERS_C_

#include "../Include/fermion_parameters.h"
#include "../Include/inverter_tricks.h"
#include "../RationalApprox/rationalapprox.h"
#include "./alloc_vars.h"
#include "./geometry.h"
#include "./inverter_full.h"
#include "./inverter_mixedp.h"
#include "./inverter_multishift_full.h"
#include "./inverter_package.h"
#include "./inverter_wrappers.h"
#include "./sp_alloc_vars.h"
#include "./sp_inverter_full.h"
#include "./sp_inverter_multishift_full.h"
#include "./struct_c_def.h"
#include "./float_double_conv.h"
#include "../Mpi/multidev.h"




void convergence_messages(int conv_importance, int inverter_status){

        if(INVERTER_FAILURE == inverter_status ){
            if(CONVERGENCE_CRITICAL == conv_importance){
                if(0==devinfo.myrank){
                    printf("\n\n\t\tERROR : inverter failed to converge in a critical region of the code. Exiting now!!\n\n");
                }
                exit(1);
            }
            if(0==devinfo.myrank){
                printf("\n\t\tWARNING : inverter failed to converge.\n");
                }
        }


}



int multishift_invert_iterations ; // global count of CG iterations

int inverter_multishift_wrapper(inverter_package ip,
        ferm_param *pars,
        RationalApprox * approx,
        vec3_soa * out,
        const vec3_soa * in,
        double res,
        int max_cg,
        int convergence_importance )
{
    int total_iterations = 0;
    int cg_return;
    int temp_conv_check;

    if(inverter_tricks.singlePInvAccelMultiInv){

        convert_double_to_float_vec3_soa(in,aux1_f);
        // multishift inverter in single precision
        float singlePMultiInvTargetRes = 8e-7f*sqrtf(sizeh);
        if(singlePMultiInvTargetRes < res)
          singlePMultiInvTargetRes = res;
        //ferm_shiftmulti_acc_f,aux1_f: global variables
        if(0 == devinfo.myrank && verbosity_lv > 3)
           printf("Multishift inverter, single precision, target res %e\n",
                   singlePMultiInvTargetRes);
        temp_conv_check = multishift_invert_f(ip.u_f,pars,approx, 
                ferm_shiftmulti_acc_f,aux1_f, singlePMultiInvTargetRes, 
                ip.loc_r_f, ip.loc_h_f, ip.loc_s_f, ip.loc_p_f, ip.ferm_shift_temp_f, 
                max_cg, &cg_return);

        convergence_messages(convergence_importance, temp_conv_check);
        total_iterations += cg_return;



        // single inversions  
        int ishift;
        for(ishift=0;ishift< approx->approx_order;
                ishift++){

            double bshift = approx->RA_b[ishift];
            printf("Shift %d, %f\n", ishift,bshift);

            convert_float_to_double_vec3_soa(&(ferm_shiftmulti_acc_f[ishift]),
                    &(out[ishift]));


            temp_conv_check = inverter_wrapper(ip,pars,&(out[ishift]),
                        in,res,max_cg,bshift,convergence_importance);

            convergence_messages(convergence_importance, temp_conv_check);
            total_iterations += cg_return;
        }

    }
    else{ 
    
        if(0 == devinfo.myrank && verbosity_lv > 3)
           printf("Multishift inverter, DOUBLE precision, target res %e\n",res);

        temp_conv_check = multishift_invert(ip.u,pars,approx,out,in,res,
                ip.loc_r,ip.loc_h,ip.loc_s,ip.loc_p,ip.ferm_shift_temp,max_cg,
                &cg_return);

        convergence_messages(convergence_importance, temp_conv_check);
        total_iterations += cg_return;
        
        
    }
    return total_iterations;

}

int inverter_wrapper(inverter_package ip,
        ferm_param *pars,
        vec3_soa * out,
        const vec3_soa * in,
        double res,
        int max_cg,
        double shift,
        int convergence_importance)
{

    int total_iterations = 0;
    int cg_return;
    int temp_conv_check;

    if(verbosity_lv > 4 ){
        printf("MPI%02d: (wrapper) max_cg  %d, restarting every %d, shift %f\n", devinfo.myrank,
                max_cg,inverter_tricks.restartingEvery,shift);

#pragma acc update host(ip.u[0:1])
        printf("MPI%02d: (wrapper) u[sizeh/2]: %f\n", devinfo.myrank, creal(ip.u[0].r0.c0[sizeh/2]));
#pragma acc update host(out[0:1])
        printf("MPI%02d: (wrapper) out[sizeh/2]: %f\n", devinfo.myrank, creal(out[0].c0[sizeh/2]));
#pragma acc update host(in[0:1])
        printf("MPI%02d: (wrapper) in[sizeh/2]: %f\n", devinfo.myrank, creal(in[0].c0[sizeh/2]));
    }






    if(inverter_tricks.useMixedPrecision)
        temp_conv_check = inverter_mixed_precision(ip,pars,out,
                    in,res,max_cg,shift,&cg_return);
    else temp_conv_check = ker_invert_openacc(ip.u, pars,out,in,res,
            ip.loc_r,ip.loc_h,ip.loc_s,ip.loc_p,max_cg,shift,&cg_return);

    convergence_messages(convergence_importance, temp_conv_check);

    total_iterations += cg_return;

   return total_iterations; 
}



#endif


