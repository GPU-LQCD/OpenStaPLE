#include <stdio.h>
#include <stdlib.h>

#ifndef __GNUC__
#include "openacc.h"
#endif
#include <sys/time.h>

#include "../Include/common_defines.h"
#include "../Include/debug.h"
#include "../Include/fermion_parameters.h"
#include "../Include/montecarlo_parameters.h"
#include "../Include/setting_file_parser.h"
#include "../Include/tell_geom_defines.h"
#include "../Mpi/multidev.h"
#include "../OpenAcc/action.h"
#include "../OpenAcc/alloc_vars.h"
#include "../OpenAcc/deviceinit.h" 
#include "../OpenAcc/find_min_max.h"
#include "../OpenAcc/float_double_conv.h"
#include "../OpenAcc/inverter_multishift_full.h"
#include "../OpenAcc/io.h"
#include "../OpenAcc/md_integrator.h"
#include "../OpenAcc/random_assignement.h"
#include "../OpenAcc/sp_alloc_vars.h"
#include "../OpenAcc/sp_fermion_matrix.h"
#include "../OpenAcc/sp_inverter_multishift_full.h"
#include "../Rand/random.h"
#include "../RationalApprox/rationalapprox.h"
#include "../DbgTools/dbgtools.h"
#include "../DbgTools/sp_dbgtools.h"
#include "./test_and_benchmarks.h"
#include "../OpenAcc/alloc_settings.h"

#ifdef MULTIDEVICE
#include "../Mpi/communications.h"
#include "../Mpi/sp_communications.h"
#include <mpi.h>
#endif


#include "../Include/stringify.h"

int id_iter;
int verbosity_lv;

void init_fake_rational_approx(RationalApprox* rational_approx,
        double a,double b, int order){
    // a0=0, an=a, bn=b  
    printf("Initializing fake rational approximation, constant term = 0, all terms \
            %e/(x+%e)\n",a,b);

    rational_approx->exponent_den = 0;
    rational_approx->exponent_num = 0;
    rational_approx->approx_order = order;
    rational_approx->lambda_min = 0.1 ; 
    rational_approx->lambda_max = 1.0 ;
    rational_approx->gmp_remez_precision = 0;
    rational_approx->error = 1.0;
    rational_approx->RA_a0 = 0; 
    int i;
    for(i = 0; i < order; i++){
        rational_approx->RA_a[i]= a;
        rational_approx->RA_b[i]= b;
    }
}



int main(int argc, char* argv[]){

    const char confname[50] = "test_conf";
    const char fermionname[50] = "test_fermion";
    char myconfname             [50];
    char myfermionname          [50];


    act_params.stout_steps = 0;
    alloc_info.NPS_tot = 1;


    // INITIALIZATION
#ifdef MULTIDEVICE
    pre_init_multidev1D(&devinfo);
#endif
    fflush(stdout);
    printf("Multishift inverter test\n");
    printf("commit: %s\n", xstr(COMMIT_HASH) );
    // INIT FERM PARAMS AND READ RATIONAL APPROX COEFFS
    printf("WELCOME! \n");
    // READ input file.
    set_global_vars_and_fermions_from_input_file(argv[1]);

    if(! test_settings.parametersAreSet){
        printf("ERROR : Test parameters are not set in %s. Is a proper TestSetting section present?\n", argv[1]);
        printf("Exiting now!\n");
        exit(1);
    }

    if(test_settings.benchmarkMode)
        printf("ENTERING BENCHMARK MODE: using fake rational approximation.\n");
    else printf("NOT ENTERING BENCHMARK MODE\n");

    initrand((unsigned int) mc_params.seed+devinfo.myrank);
    verbosity_lv = debug_settings.input_vbl;
    // INIT FERM PARAMS AND READ RATIONAL APPROX COEFFS
    if(init_ferm_params(fermions_parameters)){
        if( test_settings.benchmarkMode ){
            printf("Issues in init_ferm_params, even if\
                    this is a multishift inverter benchmark.\n");
            exit(1);
        }
        else {
            printf("Issues in init_ferm_params(), exiting now.\n");
            printf("[Maybe a rational appriximation file is missing?]\n");
            exit(1);
        }
    }
#pragma acc enter data copyin(fermions_parameters[0:alloc_info.NDiffFlavs])

    if(alloc_info.NPS_tot == 0) 
    {
        printf("ERROR: For technical reasons, at least a fermion in the input file is needed\n");
        exit(1);
    }
    else {
        printf("N fermions found: %d\n", alloc_info.NPS_tot);
    }
    //

    printf("Setting stout level to zero - no stout is needed for this benchmark.\n");
    alloc_info.stoutAllocations = 0;
    act_params.stout_steps = 0;



#ifdef MULTIDEVICE
    init_multidev1D(&devinfo);
#else
    devinfo.myrank = 0;
    devinfo.nranks = 1;
#endif
    if(0==devinfo.myrank) print_geom_defines();

#ifndef __GNUC__
    //////  OPENACC CONTEXT INITIALIZATION    //////////////////////////////////////
    // NVIDIA GPUs
    acc_device_t my_device_type = acc_device_nvidia;
    // AMD GPUs
    // acc_device_t my_device_type = acc_device_radeon;
    // Intel XeonPhi
    //acc_device_t my_device_type = acc_device_xeonphi;
    // Select device ID
#ifdef MULTIDEVICE
    printf("MPI%02d: Selecting device.\n",devinfo.myrank, devinfo.myrank%devinfo.proc_per_node );
    select_init_acc_device(my_device_type, devinfo.myrank%devinfo.proc_per_node);
#else
    printf("MPI%02d: Selecting device %d.\n",devinfo.myrank, devinfo.single_dev_choice );
    select_init_acc_device(my_device_type, devinfo.single_dev_choice);
#endif
    printf("Device Selected : OK \n");
#endif



    mem_alloc_core();
    mem_alloc_extended();
    mem_alloc_core_f();
    mem_alloc_extended_f();
    printf("Allocazione della memoria : OK \n");
    compute_nnp_and_nnm_openacc();
#pragma acc enter data copyin(nnp_openacc)
#pragma acc enter data copyin(nnm_openacc)
    printf("nn computation : OK \n");
    init_all_u1_phases(backfield_parameters,fermions_parameters);
#pragma acc update device(u1_back_phases[0:8*alloc_info.NDiffFlavs])
#pragma acc update device(u1_back_phases_f[0:8*alloc_info.NDiffFlavs])
    printf("u1_backfield initialization : OK \n");
    sprintf(myconfname             ,"%s_MPI%02d",confname             ,devinfo.myrank);
    sprintf(myfermionname          ,"%s_MPI%02d",fermionname          ,devinfo.myrank);




    // INITIALIZING GAUGE CONFIGURATION
    int conf_id_iter = 0;
    if(!read_conf_wrapper(conf_acc,mc_params.save_conf_name,
                &conf_id_iter,debug_settings.use_ildg)){
        // READS ALSO THE conf_id_iter
        printf("MPI%02d - Stored Gauge Conf \"%s\" Read : OK \n",
                devinfo.myrank, mc_params.save_conf_name);

    }
    else{
        generate_Conf_cold(conf_acc,mc_params.eps_gen);
        printf("MPI%02d - Cold Gauge Conf Generated : OK \n",
                devinfo.myrank);
        if(test_settings.saveResults)
            save_conf_wrapper(conf_acc,mc_params.save_conf_name,0,debug_settings.use_ildg);
        if(debug_settings.use_ildg)
            printf("MPI%02d: You're using ILDG format.\n", devinfo.myrank);
        conf_id_iter=0;
    }
#pragma acc update device(conf_acc[0:8])


    // init fermion
    if(!read_vec3_soa_wrapper(ferm_chi_acc,fermionname )){
        printf("MPI%02d - Fermion READ : OK \n",devinfo.myrank);
    }else{
        printf("MPI%02d - Fermion READ : FAIL, generating\n",devinfo.myrank);
        generate_vec3_soa_gauss(ferm_chi_acc);
        printf("MPI%02d - Fermion GENERATED : OK\n",devinfo.myrank);
        printf("MPI%02d - Writing file %s.\n",devinfo.myrank, fermionname);
        print_vec3_soa_wrapper(ferm_chi_acc,fermionname);
    }


#pragma acc update device(ferm_chi_acc[0:1])

    print_vec3_soa(ferm_chi_acc,myfermionname);




    // fake rational approx

    RationalApprox * rationalApproxToUse;
    /*
#pragma acc data  copy(conf_acc[0:8]) copy(ferm_chi_acc[0:1])\
copy(ferm_phi_acc[0:1])  copy(u1_back_phases[0:8*alloc_info.NDiffFlavs]) \
create(kloc_r[0:1]) create(kloc_h[0:1]) create(kloc_s[0:1]) create(kloc_p[0:1]) \
create(ferm_shiftmulti_acc[0:alloc_info.maxNeededShifts]) \
copy(conf_acc_f[0:8]) copy(ferm_chi_acc_f[0:1])\
copy(ferm_phi_acc_f[0:1])  copy(u1_back_phases_f[0:8*alloc_info.NDiffFlavs]) \
create(kloc_r_f[0:1]) create(kloc_h_f[0:1])\
create(kloc_s_f[0:1]) create(kloc_p_f[0:1]) \
create(ferm_shiftmulti_acc_f[0:alloc_info.maxNeededShifts]) \
create(k_p_shiftferm[0:alloc_info.maxApproxOrder] )\
create(k_p_shiftferm_f[0:alloc_info.maxApproxOrder] )
{
*/
    if(test_settings.benchmarkMode){
        rationalApproxToUse = (RationalApprox*) malloc(sizeof(RationalApprox));
        init_fake_rational_approx(rationalApproxToUse, 1,
                test_settings.fakeShift, 15);
    }
    else{
        printf("Using rational approximation for molecular dynamics\n");
        printf("For the first quark\n");
        //just choosing one, but rescaling it, but rescaling it first.
        double minmaxeig[2];
        generate_vec3_soa_gauss(kloc_p);
#pragma acc update device(kloc_p[0:1])
        find_min_max_eigenvalue_soloopenacc(conf_acc,fermions_parameters,kloc_r,kloc_h,kloc_p,kloc_s,minmaxeig);
        printf("Found eigenvalues of dirac operator: %e,  %e\n",
                minmaxeig[0],minmaxeig[1]);


        printf("Rescaling rational approximation...\n");
        rescale_rational_approximation(
                &(fermions_parameters[0].approx_md_mother),
                &(fermions_parameters[0].approx_md),
                minmaxeig);
        rationalApproxToUse = &(fermions_parameters[0].approx_md);//just choosing one

    }




    struct timeval t0,t1,t2,t3,t4,t5;
    int r;
    int cg_return;
    if(0 == devinfo.myrank){
        printf("Multishift Inversion, %d times, with residue %e, shift %e\n",
                test_settings.multiShiftInverterRepetitions,md_parameters.residue_metro,
                test_settings.fakeShift );
        printf("max_cg_iterations: %d\n", md_parameters.max_cg_iterations);

    }

    // double precision
    if(0 == devinfo.myrank){
        printf("####################\n");
        printf("# DOUBLE PRECISION #\n");
        printf("####################\n");
    }
    for(r=0; r<test_settings.multiShiftInverterRepetitions; r++){
        gettimeofday(&t0,NULL);
        multishift_invert(conf_acc,&fermions_parameters[0],
                rationalApproxToUse,
                ferm_shiftmulti_acc,
                ferm_chi_acc,
                2e-144,
                kloc_r,
                kloc_h,
                kloc_s,
                kloc_p,
                k_p_shiftferm,
                md_parameters.max_cg_iterations,
                &cg_return);
        gettimeofday(&t1,NULL);
        if(0==devinfo.myrank){
            double dt_cgm = (double)(t1.tv_sec - t0.tv_sec) + 
                ((double)(t1.tv_usec - t0.tv_usec)/1.0e6);
            if(0 == devinfo.myrank)
                printf("Time for 1 step of multishift inversion   : %e\n",
                        dt_cgm/cg_return);
        }
    }

    if(test_settings.saveResults){
#pragma acc update host(ferm_shiftmulti_acc[0:rationalApproxToUse->approx_order]) // update on host
        for(r=0; r<rationalApproxToUse->approx_order; r++){

            char fermionname_shift[50];
            sprintf(fermionname_shift,"fermion_shift_%d.dat",r);

            // shift fermio names
            if(0 == devinfo.myrank)
                printf("Writing file %s.\n", fermionname_shift);

            print_vec3_soa_wrapper(&ferm_shiftmulti_acc[r],fermionname_shift);
        }
    }


    // single precision
    if(0 == devinfo.myrank){
        printf("####################\n");
        printf("# SINGLE PRECISION #\n");
        printf("####################\n");
    }
    // conversion to single precision
    convert_double_to_float_su3_soa(conf_acc,conf_acc_f);
    convert_double_to_float_vec3_soa(ferm_chi_acc,ferm_chi_acc_f);

    struct timeval t0_f,t1_f,t2_f,t3_f,t4_f,t5_f;
    for(r=0; r<test_settings.multiShiftInverterRepetitions; r++){
        gettimeofday(&t0_f,NULL);
        multishift_invert_f(conf_acc_f,&fermions_parameters[0],
                rationalApproxToUse,
                ferm_shiftmulti_acc_f,
                ferm_chi_acc_f,
                1e-33,
                kloc_r_f,
                kloc_h_f,
                kloc_s_f,
                kloc_p_f,
                k_p_shiftferm_f,
                md_parameters.max_cg_iterations,
                &cg_return);
        gettimeofday(&t1_f,NULL);
        if(0==devinfo.myrank){
            double dt_cgm_f = (double)(t1_f.tv_sec - t0_f.tv_sec) + 
                ((double)(t1_f.tv_usec - t0_f.tv_usec)/1.0e6);
            if(0 == devinfo.myrank)
                printf("Time for 1 step of multishift inversion   : %e\n",
                        dt_cgm_f/cg_return);
        }
    }


    if(test_settings.saveResults){
#pragma acc update host(ferm_shiftmulti_acc_f[0:rationalApproxToUse->approx_order]) // update on host
        for(r=0; r<rationalApproxToUse->approx_order; r++){

            char fermionname_shift[50];
            sprintf(fermionname_shift,"sp_fermion_shift_%d.dat",r);

            // shift fermio names
            if(0 == devinfo.myrank)
                printf("Writing file %s.\n", fermionname_shift);

            print_vec3_soa_wrapper_f(&ferm_shiftmulti_acc_f[r],fermionname_shift);
        }
    }


    if(test_settings.benchmarkMode) free(rationalApproxToUse);

    mem_free_extended_f();
    mem_free_core_f();
    mem_free_extended();
    mem_free_core();

    printf("MPI%02d: Test completed.\n",devinfo.myrank);

#ifndef __GNUC__
    shutdown_acc_device(my_device_type);
#endif
    printf("MPI%02d: Before MPI_Finalize()...\n",devinfo.myrank);
#ifdef MULTIDEVICE
    MPI_Finalize();
#endif 





    return 0; 


}
