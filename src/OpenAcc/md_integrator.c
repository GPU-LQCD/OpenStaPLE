// MULTISTEP VERSION of minimum_norm2_B

// 2nd order Minimum Norm integrator (2MN)
// See reference hep-lat/0505020 Takaishi, De Forcrand
//
// Scheme (needs three force calculations per step):
// 2MN(dt) = exp (-l * dt * dS/dq) exp (dt/2 * p) * ...
//           exp (-(1- 2l) * dt * dS/dq) exp (dt/2 * p)  exp (-l * dt * dS/dq)
// p       : momenta
// - dS/sq : force contibution
// dt      : integration step
// l       : lambda parameter (1/6 Sexton Weingarten, 0.193... Omelyan et al)
//    
// Total scheme:
// [2MN]^N
// N       : number of Molecular Dynamics steps
//
// See reference hep-lat/0506011 Urbach et al. for the multiple time scale 
//
// p->p+a dS/dq=p+ia(-i dS/dq)=p+ia*ipdot

// the various dt involved are stored into the delta array which is allocated into UPDATE_ACC(args)

#ifndef MD_INTEGRATOR_C
#define MD_INTEGRATOR_C


#include "../DbgTools/dbgtools.h"
#include "../Include/common_defines.h"
#include "../Include/debug.h"
#include "../Include/fermion_parameters.h"
#include "../Mpi/multidev.h"
#include "./action.h"
#include "./alloc_vars.h"
#include "./fermion_force.h"
#include "./fermionic_utilities.h"
#include "./ipdot_gauge.h"
#include "./md_integrator.h"
#include "./md_parameters.h"
#include "./struct_c_def.h"
#include "./su3_utilities.h"
#include "./inverter_package.h"
#include "../tests_and_benchmarks/test_and_benchmarks.h"
#include "./alloc_settings.h"
#include "./topological_force.h"

#include <sys/time.h>

#define TOPO_MACRO 1
#define TOPO_MICRO 0

#ifdef MULTIDEVICE
#include "../Mpi/communications.h"
#endif

#ifdef MULTIDEVICE

#if defined(USE_MPI_CUDA_AWARE) || defined(__GNUC__)
void multistep_2MN_gauge_async_bloc(su3_soa *tconf_acc_old, su3_soa *tconf_acc_new,
        su3_soa *local_staples, tamat_soa *tipdot,thmat_soa *tmomenta, int omelyan_index)
{

    if(verbosity_lv > 3) printf("DOUBLE PRECISION VERSION OF MULTISTEP_2MN_GAUGE_ASYNC_BLOc\n");

    if(md_dbg_print_count<debug_settings.md_dbg_print_max_count){
        char conffilename[50];
        char momfilename[50];
        char ipdotfilename[50];
        sprintf(conffilename,"conf_md_%d_%d",devinfo.myrank, md_dbg_print_count);
        sprintf(momfilename,"tmomenta_%d_%d",devinfo.myrank, md_dbg_print_count);
        sprintf(ipdotfilename,"tipdot_%d_%d",devinfo.myrank, md_dbg_print_count);

        dbg_print_su3_soa(tconf_acc_old,conffilename, 1);
        print_thmat_soa(tmomenta,momfilename);
        print_tamat_soa(tipdot,ipdotfilename);

        char global_conffilename[50];
        char global_momfilename[50];
        char global_ipdotfilename[50];

        sprintf(global_conffilename ,"global_conf_md_%d" , md_dbg_print_count);
        sprintf(global_momfilename  ,"global_tmomenta_%d", md_dbg_print_count);
        sprintf(global_ipdotfilename,"global_tipdot_%d"  , md_dbg_print_count);

        print_gl3_soa_wrapper(tconf_acc_old,global_conffilename);
        print_thmat_soa_wrapper(tmomenta,global_momfilename);
        print_tamat_soa_wrapper(tipdot,global_ipdotfilename);

        md_dbg_print_count++;
    }

    struct timeval t[9];
    double dt[7], dtcom;
    MPI_Request send_border_requests[96]; 
    MPI_Request recv_border_requests[96];
    if(verbosity_lv > 2 && 0 == devinfo.myrank) printf("\tMPI%02d - In async bloc - Index %d\n",
            devinfo.myrank, omelyan_index);

    if(verbosity_lv > 3 && 0 == devinfo.myrank) printf("\t\tMPI%02d - calc_ipdot_gauge_d3c()\n", devinfo.myrank);

    gettimeofday ( &t[0], NULL );
    // NOT embarrassingly parallel
    calc_ipdot_gauge_soloopenacc_d3c(tconf_acc_old,local_staples,tipdot,
            HALO_WIDTH,GAUGE_HALO); 
    calc_ipdot_gauge_soloopenacc_d3c(tconf_acc_old,local_staples,tipdot,
            nd3-HALO_WIDTH-GAUGE_HALO,GAUGE_HALO); 

    if(verbosity_lv > 3 && 0 == devinfo.myrank) printf("\t\tMPI%02d - mom_sum_mult_d3c()\n", devinfo.myrank);
    gettimeofday ( &t[1], NULL );
    // embarrassingly parallel
    mom_sum_mult_d3c(tmomenta,tipdot,deltas_Omelyan,omelyan_index,
            HALO_WIDTH,GAUGE_HALO);
    mom_sum_mult_d3c(tmomenta,tipdot,deltas_Omelyan,omelyan_index,
            nd3-HALO_WIDTH-GAUGE_HALO,GAUGE_HALO); 


    if(verbosity_lv > 3 && 0 == devinfo.myrank) printf("\t\tMPI%02d - mom_exp_times_conf_soloopenacc_d3c()\n",devinfo.myrank);
    gettimeofday ( &t[2], NULL );
    // this function should have differen in and out for the gauge conf
    // embarrassingly parallel
    mom_exp_times_conf_soloopenacc_d3c(
            tconf_acc_old, tconf_acc_new, tmomenta,
            deltas_Omelyan,4,
            HALO_WIDTH,GAUGE_HALO);
    // this function should have differen in and out for the gauge conf
    // embarrassingly parallel
    mom_exp_times_conf_soloopenacc_d3c(
            tconf_acc_old, tconf_acc_new, tmomenta,
            deltas_Omelyan,4,
            nd3-HALO_WIDTH-GAUGE_HALO,GAUGE_HALO); 
    gettimeofday ( &t[3], NULL );

    if(verbosity_lv > 2 && 0 == devinfo.myrank) printf("\tMPI%02d -communicate_su3_borders_async()\n",devinfo.myrank);

    communicate_su3_borders_async(tconf_acc_new,GAUGE_HALO,
            send_border_requests,recv_border_requests);


    gettimeofday ( &t[4], NULL );
    if(verbosity_lv > 3 && 0 == devinfo.myrank) printf("\t\tMPI%02d - calc_ipdot_gauge_bulk()\n", devinfo.myrank);
    // NOT embarrassingly parallel
    calc_ipdot_gauge_soloopenacc_bulk(tconf_acc_old,local_staples,tipdot);
    gettimeofday ( &t[5], NULL );

    if(verbosity_lv > 3 && 0 == devinfo.myrank) printf("\t\tMPI%02d - mom_sum_mult_bulk()\n", devinfo.myrank);
    // embarrassingly parallel
    mom_sum_mult_bulk(tmomenta,tipdot,deltas_Omelyan,omelyan_index);
    gettimeofday ( &t[6], NULL );

    if(verbosity_lv > 3 && 0 == devinfo.myrank) printf("\t\tMPI%02d - mom_exp_times_conf_soloopenacc_bulk()\n",devinfo.myrank);
    // this function should have differen in and out for the gauge conf
    // embarrassingly parallel
    mom_exp_times_conf_soloopenacc_bulk(
            tconf_acc_old, tconf_acc_new, tmomenta,
            deltas_Omelyan,4);
    gettimeofday ( &t[7], NULL );


    MPI_Waitall(96,send_border_requests,MPI_STATUSES_IGNORE);
    MPI_Waitall(96,recv_border_requests,MPI_STATUSES_IGNORE);
    gettimeofday ( &t[8], NULL );

    if(verbosity_lv > 2 && 0 == devinfo.myrank) printf("\tMPI%02d - End of async bloc, index %d\n",
            devinfo.myrank, omelyan_index);




    // DIAGNOSTICS - PERFORMANCE MEASUREMENTS
    int di;
    for(di = 0 ; di<7; di++)
        dt[di] = (double)(t[di+1].tv_sec - t[di].tv_sec) +
            ((double)(t[di+1].tv_usec - t[di].tv_usec)/1.0e6);

    dtcom =  (double)(t[8].tv_sec - t[3].tv_sec) +
        ((double)(t[8].tv_usec - t[3].tv_usec)/1.0e6);

    gauge_mdtimes.calcIpdotTimeBorder      += dt[0];
    gauge_mdtimes.calcIpdotTimeBulk        += dt[4];
    gauge_mdtimes.momSumMultTimeBorder     += dt[1];
    gauge_mdtimes.momSumMultTimeBulk       += dt[5];
    gauge_mdtimes.momExpTimesConfTimeBorder+= dt[2];
    gauge_mdtimes.momExpTimesConfTimeBulk  += dt[6];
    gauge_mdtimes.communicationsStartTime  += dt[3];
    gauge_mdtimes.communicationsTime       += dtcom;
    gauge_mdtimes.count++ ;



    if(verbosity_lv > 2 && 0 == devinfo.myrank){
        printf("\t|          Function              \t|Border\t|Bulk\t|\n");
        printf("\t| Calc ipdot                     \t|%e|%e|\n",dt[0],dt[4] );
        printf("\t| Mom sum mult                   \t|%e|%e|\n",dt[1],dt[5] );
        printf("\t| mom_exp_times_conf_soloopenacc \t|%e|%e|\n",dt[2],dt[6] );
        printf("\t| Communications starting        \t|%e|\n",dt[3] );
        printf("\t| Communications total           \t|%e|\n",dtcom );

    }

    if(md_dbg_print_count<debug_settings.md_dbg_print_max_count
            && 1 == debug_settings.md_dbg_be_verbose ){
            
        char genericfilename[50];
        // staples
        sprintf(genericfilename,"staples_%d_%d",
                devinfo.myrank, md_dbg_print_count);
        dbgprint_gl3_soa(local_staples,genericfilename,1000);

    }

}
#endif 


void multistep_2MN_gauge_async(su3_soa *tconf_acc,su3_soa *local_staples,tamat_soa *tipdot,thmat_soa *tmomenta)
{
    if(verbosity_lv > 3) printf("DOUBLE PRECISION VERSION OF MULTISTEP_2MN_GAUGE_ASYNC\n");

#if defined(USE_MPI_CUDA_AWARE) || defined(__GNUC__)
    int md;
    if(verbosity_lv>1) 
        printf("MPI%02d - Performing Gauge substeps - Async communications\n",
                devinfo.myrank);

    // tconf_acc[0:8]   --> old conf
    // tconf_acc[8:16]  --> new conf

    multistep_2MN_gauge_async_bloc(tconf_acc,&tconf_acc[8],// for async we need 2 confs
            local_staples, tipdot,tmomenta,3);

    for(md=1; md<md_parameters.gauge_scale; md++){
        if(verbosity_lv > 2) printf("MPI%02d - Gauge step %d of %d...\n",
                devinfo.myrank,md,md_parameters.gauge_scale);

        multistep_2MN_gauge_async_bloc(&tconf_acc[8],tconf_acc,//for async we need 2 confs
                local_staples, tipdot,tmomenta,5);

        multistep_2MN_gauge_async_bloc(tconf_acc,&tconf_acc[8],//for async we need 2 confs
                local_staples, tipdot,tmomenta,6);

    }
    if(verbosity_lv > 2) printf("MPI%02d - Last Gauge step of %d...\n",
            devinfo.myrank,md_parameters.gauge_scale);

    multistep_2MN_gauge_async_bloc(&tconf_acc[8],tconf_acc, 
            local_staples, tipdot,tmomenta,5);

    // LAST STEP IS NOT ASYNCED
    // Step for the P
    // P' = P - l*dt*dS/dq
    // deltas_Omelyan[3]=-cimag(ieps_acc)*lambda*scale;
    calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);
    mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,3);

#else
    printf("ERROR, Async gauge evolution cannot be performed on accelerators,\n");
    printf("       if USE_MPI_CUDA_AWARE is not #defined. Exiting now.\n");
    MPI_Finalize();
    exit(1);
#endif 


}
#endif

void multistep_2MN_gauge_bloc(su3_soa *tconf_acc,
        su3_soa *local_staples, tamat_soa *tipdot,thmat_soa *tmomenta,
        int omelyan_index)
{
    if(verbosity_lv > 3) printf("DOUBLE PRECISION VERSION OF MULTISTEP_2MN_GAUGE_BLOC\n");

    if(verbosity_lv > 2 && 0 == devinfo.myrank) printf("\tMPI%02d - In bloc - Index %d\n",
            devinfo.myrank, omelyan_index);

    if(md_dbg_print_count<debug_settings.md_dbg_print_max_count){
        char conffilename[50];
        char momfilename[50];
        char ipdotfilename[50];
        sprintf(conffilename,"conf_md_%d_%d",devinfo.myrank, md_dbg_print_count);
        sprintf(momfilename,"tmomenta_%d_%d",devinfo.myrank, md_dbg_print_count);
        sprintf(ipdotfilename,"tipdot_%d_%d",devinfo.myrank, md_dbg_print_count);

        dbg_print_su3_soa(tconf_acc,conffilename, 1);
        print_thmat_soa(tmomenta,momfilename);
        print_tamat_soa(tipdot,ipdotfilename);

        char global_conffilename[50];
        char global_momfilename[50];
        char global_ipdotfilename[50];

        sprintf(global_conffilename ,"global_conf_md_%d" , md_dbg_print_count);
        sprintf(global_momfilename  ,"global_tmomenta_%d", md_dbg_print_count);
        sprintf(global_ipdotfilename,"global_tipdot_%d"  , md_dbg_print_count);

        print_gl3_soa_wrapper(tconf_acc,global_conffilename);
        print_thmat_soa_wrapper(tmomenta,global_momfilename);
        print_tamat_soa_wrapper(tipdot,global_ipdotfilename);

        md_dbg_print_count++;
    }


    // Step for the P
    // P' = P - l*dt*dS/dq
    // deltas_Omelyan[3]=-cimag(ieps_acc)*scale*lambda;
    // deltas_Omelyan[5]=-cimag(ieps_acc)*(1.0-2.0*lambda)*scale;
    // deltas_Omelyan[6]=-cimag(ieps_acc)*2.0*lambda*scale;    
    struct timeval t[5];
    double dt[3], dtcom;

    if(verbosity_lv > 3 && 0 == devinfo.myrank) printf("\t\tMPI%02d - calc_ipdot_gauge()\n", devinfo.myrank);
    gettimeofday ( &t[0], NULL );
    calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);
    if(verbosity_lv > 3 && 0 == devinfo.myrank) printf("\t\tMPI%02d - mom_sum_mult()\n", devinfo.myrank);
    gettimeofday ( &t[1], NULL );
    mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,omelyan_index);
    if(verbosity_lv > 3 && 0 == devinfo.myrank) printf("\t\tMPI%02d - mom_exp_times_conf_soloopenacc()\n",devinfo.myrank);

    // Step for the Q
    // Q' = exp[dt/2 *i P] Q
    // deltas_Omelyan[4]=cimag(iepsh_acc)*scale;
    gettimeofday ( &t[2], NULL );
    mom_exp_times_conf_soloopenacc( tconf_acc, tmomenta,
            deltas_Omelyan,4);
    gettimeofday ( &t[3], NULL );

#ifdef MULTIDEVICE
    communicate_su3_borders(tconf_acc,GAUGE_HALO);
#endif
    gettimeofday ( &t[4], NULL );
    if(verbosity_lv > 3) printf("\tMPI%02d - End of bloc, index %d\n",
            devinfo.myrank, omelyan_index);



    // DIAGNOSTICS - PERFORMANCE MEASUREMENTS
    int di;
    for(di = 0 ; di<3; di++)
        dt[di] = (double)(t[di+1].tv_sec - t[di].tv_sec) +
            ((double)(t[di+1].tv_usec - t[di].tv_usec)/1.0e6);

    dtcom =  (double)(t[4].tv_sec - t[3].tv_sec) +
        ((double)(t[4].tv_usec - t[3].tv_usec)/1.0e6);

    gauge_mdtimes.calcIpdotTimeBulk        += dt[0];
    gauge_mdtimes.momSumMultTimeBulk       += dt[1];
    gauge_mdtimes.momExpTimesConfTimeBulk  += dt[2];
    gauge_mdtimes.communicationsTime       += dtcom;
    gauge_mdtimes.count++ ;

    if(verbosity_lv > 2 && 0 == devinfo.myrank){
        printf("\t|          Function              \t|Bulk\t|\n");
        printf("\t| Calc ipdot                     \t|%e|\n",dt[0] );
        printf("\t| Mom sum mult                   \t|%e|\n",dt[1] );
        printf("\t| mom_exp_times_conf_soloopenacc \t|%e|\n",dt[2] );
        printf("\t| Communications total           \t|%e|\n",dtcom );

    }

}

void multistep_2MN_gauge(su3_soa *tconf_acc,su3_soa *local_staples,tamat_soa *tipdot,thmat_soa *tmomenta)
{
    if(verbosity_lv > 3) printf("DOUBLE PRECISION VERSION OF MULTISTEP_2MN_GAUGE\n");
    int md;
    if(verbosity_lv>1) 
        printf("MPI%02d - Performing Gauge substeps\n",
                devinfo.myrank);

    multistep_2MN_gauge_bloc(tconf_acc,local_staples, tipdot,tmomenta,3);

    for(md=1; md<md_parameters.gauge_scale; md++){
        if(verbosity_lv > 2) printf("MPI%02d - Gauge step %d of %d...\n",
                devinfo.myrank,md,md_parameters.gauge_scale);

        multistep_2MN_gauge_bloc(tconf_acc,local_staples, tipdot,tmomenta,5);

        multistep_2MN_gauge_bloc(tconf_acc,local_staples, tipdot,tmomenta,6);

    }
    if(verbosity_lv > 2) printf("MPI%02d - Last Gauge step of %d...\n",
            devinfo.myrank,md_parameters.gauge_scale);

    multistep_2MN_gauge_bloc(tconf_acc, local_staples, tipdot,tmomenta,5);

    // Step for the P
    // P' = P - l*dt*dS/dq
    // deltas_Omelyan[3]=-cimag(ieps_acc)*lambda*scale;
    calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);
    mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,3);

}

/*
   void multistep_2MN_gauge(su3_soa *tconf_acc,su3_soa *local_staples,tamat_soa *tipdot,thmat_soa *tmomenta)
   {
   if(verbosity_lv>1) 
   printf("MPI%02d - Performing Gauge substeps\n",
   devinfo.myrank);
   int md;
// Step for the P
// P' = P - l*dt*dS/dq
// deltas_Omelyan[3]=-cimag(ieps_acc)*scale*lambda;

calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot); 
#ifdef DEBUG_MD
char conffilename[50];
char momfilename[50];
char ipdotfilename[50];
if(!already_printed_debug){
sprintf(conffilename,"conf_md_0_%s",devinfo.myrankstr);
sprintf(momfilename,"tmomenta_0_%s",devinfo.myrankstr);
sprintf(ipdotfilename,"tipdot_0_%s",devinfo.myrankstr);
dbg_print_su3_soa(tconf_acc,conffilename, 1);
print_thmat_soa(tmomenta,momfilename);
print_tamat_soa(tipdot,ipdotfilename);
}
#endif

mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,3);
#ifdef DEBUG_MD
if(!already_printed_debug){
sprintf(momfilename,"tmomenta_1_%s",devinfo.myrankstr);
print_thmat_soa(tmomenta,momfilename);
}
#endif
//---------------------------
// Step for the Q
// Q' = exp[dt/2 *i P] Q
// deltas_Omelyan[4]=cimag(iepsh_acc)*scale;
mom_exp_times_conf_soloopenacc(tconf_acc,tmomenta,deltas_Omelyan,4);
#ifdef MULTIDEVICE
communicate_su3_borders(tconf_acc,GAUGE_HALO);  
#endif
//-------------------------




for(md=1; md<md_parameters.gauge_scale; md++){
if(verbosity_lv > 2) printf("MPI%02d - Gauge step %d of %d...\n",
devinfo.myrank,md,md_parameters.gauge_scale);
// Step for the P
// P' = P - (1-2l)*dt*dS/dq
// deltas_Omelyan[5]=-cimag(ieps_acc)*(1.0-2.0*lambda)*scale;
calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);

mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,5);
// Step for the Q
// Q' = exp[dt/2 *i P] Q
// deltas_Omelyan[4]=cimag(iepsh_acc)*scale;
mom_exp_times_conf_soloopenacc(tconf_acc,tmomenta,deltas_Omelyan,4);
#ifdef MULTIDEVICE
communicate_su3_borders(tconf_acc, GAUGE_HALO);  
#endif
// Step for the P
// P' = P - 2l*dt*dS/dq
// deltas_Omelyan[6]=-cimag(ieps_acc)*2.0*lambda*scale;
calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);

mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,6);
//--------------
// Step for the Q
// Q' = exp[dt/2 *i P] Q
// deltas_Omelyan[4]=cimag(iepsh_acc)*scale;
mom_exp_times_conf_soloopenacc(tconf_acc,tmomenta,deltas_Omelyan,4);
#ifdef MULTIDEVICE
communicate_su3_borders(tconf_acc, GAUGE_HALO);  
#endif
//---------------------



}
if(verbosity_lv > 2) printf("MPI%02d - Last Gauge step of %d...\n",
        devinfo.myrank,md_parameters.gauge_scale);


// Step for the P
// P' = P - (1-2l)*dt*dS/dq
calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);

// calc_ipdot_gauge();
// deltas_Omelyan[5]=-cimag(ieps_acc)*(1.0-2.0*lambda)*scale;
mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,5);

#ifdef DEBUG_MD
if(!already_printed_debug){
    sprintf(conffilename,"conf_md_1_%s",devinfo.myrankstr);
    sprintf(momfilename,"tmomenta_3_%s",devinfo.myrankstr);
    sprintf(ipdotfilename,"tipdot_2_%s",devinfo.myrankstr);
    dbg_print_su3_soa(tconf_acc,conffilename, 2);
    print_thmat_soa(tmomenta,momfilename);
    print_tamat_soa(tipdot,ipdotfilename);
}
already_printed_debug = 1;
#endif

// Step for the Q
// Q' = exp[dt/2 *i P] Q
// deltas_Omelyan[4]=cimag(iepsh_acc)*scale;
mom_exp_times_conf_soloopenacc(tconf_acc,tmomenta,deltas_Omelyan,4);
#ifdef MULTIDEVICE
communicate_su3_borders(tconf_acc, GAUGE_HALO);  
#endif
// Step for the P
// P' = P - l*dt*dS/dq
// deltas_Omelyan[3]=-cimag(ieps_acc)*lambda*scale;
calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);
mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,3);


}
*/


void multistep_2MN_SOLOOPENACC( tamat_soa * tipdot_acc,
        su3_soa  * tconf_acc,
#if (defined STOUT_FERMIONS) || (defined STOUT_TOPO)
        su3_soa  * tstout_conf_acc_arr, // huge parking for stouting
#endif
        su3_soa  * tauxbis_conf_acc, 
        su3_soa  * taux_conf_acc,
        ferm_param * tfermions_parameters,// [nflavs]
        int tNDiffFlavs,
        vec3_soa * ferm_in_acc, //[NPS_tot], will be ferm_chi_acc
        vec3_soa * tferm_shiftmulti_acc,// parking variable [maxNeededShifts]
        inverter_package ip,
        thmat_soa * tmomenta,
        dcomplex_soa * local_sums,
        double res,
        const int max_cg)
{
    if(verbosity_lv > 3) printf("DOUBLE PRECISION VERSION OF MULTISTEP_2MN_SOLOOPENACC\n");

    ipdot_f_reset = 1;
    ipdot_g_reset = 1;

    vec3_soa * invOuts = tferm_shiftmulti_acc;
    nMdInversionPerformed = 0; // used to recycle inversion results
    // after first force calculation is done
    int ishift;
    for(ishift =0; ishift < alloc_info.maxNeededShifts; ishift++)
        set_vec3_soa_to_zero(&tferm_shiftmulti_acc[ishift]);

    int md;

    // Step for the P
    // P' = P - l*dt*dS/dq
    //    deltas_Omelyan[0]=-cimag(ieps_acc)*lambda;
    fermion_force_soloopenacc(tconf_acc, 
#ifdef STOUT_FERMIONS
            tstout_conf_acc_arr,
#endif
            tauxbis_conf_acc, // parkeggio
            tipdot_acc, tfermions_parameters, tNDiffFlavs, 
            ferm_in_acc, res, taux_conf_acc, invOuts, ip,max_cg);

    if(verbosity_lv > 4) printf("MPI%02d - Calculated first fermion force\n", 
            devinfo.myrank);
    
    if(TOPO_MACRO == 1 && act_params.topo_action == 1)
      {
    	calc_ipdot_topo(tconf_acc,  
#ifdef STOUT_TOPO
			tstout_conf_acc_arr,
#endif
			tauxbis_conf_acc,
			taux_conf_acc,
			tipdot_acc);
	
	if(verbosity_lv > 4) printf("MPI%02d - Calculated first topological force\n", 
				    devinfo.myrank);
	
      }

    mom_sum_mult(tmomenta,tipdot_acc,deltas_Omelyan,0);
    for(md=1; md<md_parameters.no_md; md++){

        if(md_parameters.extrapolateInvsForce){
            printf("ERROR, not implemented correctly! %s : %d",__FILE__,__LINE__); exit(1);
            invOuts = &tferm_shiftmulti_acc[totalMdShifts];
        }

        printf("\n\nMPI%02d\t\tRUNNING MD STEP %d OF %d...\n",
                devinfo.myrank, md, md_parameters.no_md);
        // Step for the Q
        // Q' = exp[dt/2 *i P] Q
#ifdef MULTIDEVICE
        if(devinfo.async_comm_gauge)
            multistep_2MN_gauge_async(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);
        else 
#endif
	multistep_2MN_gauge(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);
        // Step for the P
        // P' = P - (1-2l)*dt*dS/dq
        // deltas_Omelyan[1]=-cimag(ieps_acc)*(1.0-2.0*lambda);
        fermion_force_soloopenacc(tconf_acc, 
#ifdef STOUT_FERMIONS
                tstout_conf_acc_arr,
#endif
                tauxbis_conf_acc, // parkeggio
                tipdot_acc, tfermions_parameters, tNDiffFlavs,
                ferm_in_acc, res, taux_conf_acc, invOuts,
                ip,max_cg);

	if(TOPO_MACRO == 1 && act_params.topo_action == 1)
	  {
	    calc_ipdot_topo(tconf_acc,  
#ifdef STOUT_TOPO
			    tstout_conf_acc_arr,
#endif
			    tauxbis_conf_acc,
			    taux_conf_acc,
			    tipdot_acc);
	  }

	mom_sum_mult(tmomenta,tipdot_acc,deltas_Omelyan,1);

        // Step for the Q
        // Q' = exp[dt/2 *i P] Q
#ifdef MULTIDEVICE
        if(devinfo.async_comm_gauge)
            multistep_2MN_gauge_async(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);
        else
#endif           
            multistep_2MN_gauge(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);

        if(md_parameters.extrapolateInvsForce){
            printf("ERROR, not implemented correctly! %s : %d",__FILE__,__LINE__); exit(1);
            calc_new_trialsol_for_inversion_in_force(totalMdShifts,tferm_shiftmulti_acc,
                    nMdInversionPerformed); // nMdInversionPerformed even - trial sol 
            // in the first half of the vector
            invOuts = tferm_shiftmulti_acc;
        }        
        // Step for the P
        // P' = P - 2l*dt*dS/dq
        // deltas_Omelyan[2]=-cimag(ieps_acc)*(2.0*lambda);
        fermion_force_soloopenacc(tconf_acc,
#ifdef STOUT_FERMIONS
                tstout_conf_acc_arr, 
#endif
                tauxbis_conf_acc, // parkeggio
                tipdot_acc, tfermions_parameters, tNDiffFlavs,
                ferm_in_acc, res, taux_conf_acc, invOuts, ip,max_cg);


	if(TOPO_MACRO == 1 && act_params.topo_action == 1)
	  {

	    calc_ipdot_topo(tconf_acc,  
#ifdef STOUT_TOPO
			    tstout_conf_acc_arr,
#endif
			    tauxbis_conf_acc,
			    taux_conf_acc,
			    tipdot_acc);
	
	  }

	    mom_sum_mult(tmomenta,tipdot_acc,deltas_Omelyan,2);

        if(md_parameters.extrapolateInvsForce){
            printf("ERROR, not implemented correctly! %s : %d",__FILE__,__LINE__); exit(1);
            calc_new_trialsol_for_inversion_in_force(totalMdShifts,tferm_shiftmulti_acc,
                    nMdInversionPerformed); // nMdInversionPerformed odd - trial sol 
            // in the second half of the vector
            invOuts = &tferm_shiftmulti_acc[totalMdShifts];
        }        


    }  

    printf("\n\nMPI%02d\t\tRUNNING LAST MD STEP OF %d...\n",
            devinfo.myrank,md_parameters.no_md);
    // Step for the Q
    // Q' = exp[dt/2 *i P] Q
#ifdef MULTIDEVICE
    if(devinfo.async_comm_gauge)
        multistep_2MN_gauge_async(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);
    else
#endif
        multistep_2MN_gauge(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);
    
    
    // Step for the P
    // P' = P - (1-2l)*dt*dS/dq
    // deltas_Omelyan[1]=-cimag(ieps_acc)*(1.0-2.0*lambda);
    fermion_force_soloopenacc(tconf_acc,
#ifdef STOUT_FERMIONS
            tstout_conf_acc_arr,
#endif
            tauxbis_conf_acc, // parkeggio
            tipdot_acc, tfermions_parameters, tNDiffFlavs,
            ferm_in_acc, res, taux_conf_acc, invOuts, 
            ip,max_cg);

    
    if(TOPO_MACRO == 1 && act_params.topo_action == 1)
      {

	calc_ipdot_topo(tconf_acc,  
#ifdef STOUT_TOPO
			tstout_conf_acc_arr,
#endif
			tauxbis_conf_acc,
			taux_conf_acc,
			tipdot_acc);
	
      }
    
	mom_sum_mult(tmomenta,tipdot_acc,deltas_Omelyan,1);
    
    // Step for the Q
    // Q' = exp[dt/2 *i P] Q
#ifdef MULTIDEVICE
    if(devinfo.async_comm_gauge)
        multistep_2MN_gauge_async(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);
    else 
#endif
        multistep_2MN_gauge(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);

    if(md_parameters.extrapolateInvsForce){ 
        printf("ERROR, not implemented correctly! %s : %d",__FILE__,__LINE__); exit(1);
        calc_new_trialsol_for_inversion_in_force(totalMdShifts,tferm_shiftmulti_acc,
                nMdInversionPerformed); // nMdInversionPerformed even - trial sol 
        // in the first half of the vector
        invOuts = tferm_shiftmulti_acc;
    }        

    // Step for the P
    // P' = P - l*dt*dS/dq
    // deltas_Omelyan[0]=-cimag(ieps_acc)*lambda;
    fermion_force_soloopenacc(tconf_acc,
#ifdef STOUT_FERMIONS
            tstout_conf_acc_arr, 
#endif
            tauxbis_conf_acc, // parkeggio
            tipdot_acc, tfermions_parameters, tNDiffFlavs,
            ferm_in_acc, res, taux_conf_acc, invOuts, ip,max_cg);

    
    if(TOPO_MACRO == 1 && act_params.topo_action == 1)
      {

	calc_ipdot_topo(tconf_acc,  
#ifdef STOUT_TOPO
			tstout_conf_acc_arr,
#endif
			tauxbis_conf_acc,
			taux_conf_acc,
			tipdot_acc);
	
      }
    
	mom_sum_mult(tmomenta,tipdot_acc,deltas_Omelyan,0);

}// end multistep_2MN_SOLOOPENACC()

#endif

