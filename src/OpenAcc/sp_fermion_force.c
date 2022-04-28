#ifndef SP_FERMION_FORCE_C
#define SP_FERMION_FORCE_C


#include "../DbgTools/sp_dbgtools.h"
#include "../Include/common_defines.h"
#include "../Include/debug.h"
#include "../Include/fermion_parameters.h"
#include "../Mpi/multidev.h"
#include "./action.h"
#include "./sp_alloc_vars.h"
#include "./sp_backfield.h"
#include "./sp_fermion_force.h"
#include "./sp_fermion_force_utilities.h"
#include "./sp_inverter_multishift_full.h"
#include "../Include/inverter_tricks.h"
#include "./inverter_package.h"
#include "./inverter_wrappers.h"
#include "./md_parameters.h"
#include "./sp_plaquettes.h"
#include "./sp_stouting.h"
#include "./sp_struct_c_def.h"
#include "./sp_su3_measurements.h"
#include "./sp_su3_utilities.h"
#include "./sp_inverter_multishift_full.h"
#include "./sp_inverter_full.h"




#ifndef __GNUC__
#define TIMING_FERMION_FORCE
#endif

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
#define __restrict
#endif

#ifdef MULTIDEVICE
#include "../Mpi/sp_communications.h"
#endif 


extern int verbosity_lv;

void compute_sigma_from_sigma_prime_backinto_sigma_prime_f(  __restrict su3_soa_f    * Sigma, // la var globale e' auxbis_conf_acc_f [sia input che ouptput]
        __restrict thmat_soa_f  * Lambda, // la var globale e' aux_th_f
        __restrict tamat_soa_f  * QA, // la var globale e' aux_ta_f
        __restrict const su3_soa_f * const U,// la var globale e' .... per adesso conf_acc_f
				__restrict su3_soa_f * const TMP, //la var globale e' aux_conf_acc_f //PARCHEGGIO??
        const int istopo //istopo = {0,1} -> rho={fermrho,toporho}
        ){
    
    if(verbosity_lv > 3) printf("SINGLE PRECISION VERSION OF COMPUTE_SIGMA_FROM_SIGMA_PRIME_BACKINTO_SIGMA_PRIME\n");




    if(verbosity_lv > 2) printf("MPI%02d:\t\tSIGMA_PRIME --> SIGMA\n",
            devinfo.myrank);
    if(verbosity_lv > 5){// printing stuff
#pragma acc update host(Sigma[0:8])
        printf("-------------Sigma[old]------------------\n");                                                                                             
        printf("Sigma[old]00 = %f + (%f)*I\n",crealf(Sigma[0].r0.c0[0]),cimagf(Sigma[0].r0.c0[0]));                                               
        printf("Sigma[old]01 = %f + (%f)*I\n",crealf(Sigma[0].r0.c1[0]),cimagf(Sigma[0].r0.c1[0]));                                               
        printf("Sigma[old]02 = %f + (%f)*I\n",crealf(Sigma[0].r0.c2[0]),cimagf(Sigma[0].r0.c2[0]));                                               
        printf("Sigma[old]10 = %f + (%f)*I\n",crealf(Sigma[0].r1.c0[0]),cimagf(Sigma[0].r1.c0[0]));                                               
        printf("Sigma[old]11 = %f + (%f)*I\n",crealf(Sigma[0].r1.c1[0]),cimagf(Sigma[0].r1.c1[0]));                                               
        printf("Sigma[old]12 = %f + (%f)*I\n",crealf(Sigma[0].r1.c2[0]),cimagf(Sigma[0].r1.c2[0]));                                               
        printf("Sigma[old]20 = %f + (%f)*I\n",crealf(Sigma[0].r2.c0[0]),cimagf(Sigma[0].r2.c0[0]));                                               
        printf("Sigma[old]21 = %f + (%f)*I\n",crealf(Sigma[0].r2.c1[0]),cimagf(Sigma[0].r2.c1[0]));                                               
        printf("Sigma[old]22 = %f + (%f)*I\n\n",crealf(Sigma[0].r2.c2[0]),cimagf(Sigma[0].r2.c2[0]));                
    
#pragma acc update self(U[0:8])
        printf("-------------U------------------\n");                                                                                             
        printf("U00 = %.18lf + (%.18lf)*I\n",creal(U[0].r0.c0[0]),cimag(U[0].r0.c0[0]));                                               
        printf("U01 = %.18lf + (%.18lf)*I\n",creal(U[0].r0.c1[0]),cimag(U[0].r0.c1[0]));                                               
        printf("U02 = %.18lf + (%.18lf)*I\n",creal(U[0].r0.c2[0]),cimag(U[0].r0.c2[0]));                                               
        printf("U10 = %.18lf + (%.18lf)*I\n",creal(U[0].r1.c0[0]),cimag(U[0].r1.c0[0]));                                               
        printf("U11 = %.18lf + (%.18lf)*I\n",creal(U[0].r1.c1[0]),cimag(U[0].r1.c1[0]));                                               
        printf("U12 = %.18lf + (%.18lf)*I\n",creal(U[0].r1.c2[0]),cimag(U[0].r1.c2[0]));                                               
        printf("U20 = %.18lf + (%.18lf)*I\n",creal(U[0].r2.c0[0]),cimag(U[0].r2.c0[0]));                                               
        printf("U21 = %.18lf + (%.18lf)*I\n",creal(U[0].r2.c1[0]),cimag(U[0].r2.c1[0]));                                               
        printf("U22 = %.18lf + (%.18lf)*I\n\n",creal(U[0].r2.c2[0]),cimag(U[0].r2.c2[0]));                
    }

    set_su3_soa_to_zero_f(TMP);

    calc_loc_staples_nnptrick_all_f(U,TMP);
    if(verbosity_lv > 4)printf("MPI%02d:\t\tcomputed staples  \n",
            devinfo.myrank);
#ifdef MULTIDEVICE
    communicate_gl3_borders_f(TMP,1);
#endif
    RHO_times_conf_times_staples_ta_part_f(U,TMP,QA,istopo);

#ifdef MULTIDEVICE
    communicate_tamat_soa_borders_f(QA,1);
#endif

    // check: TMP = local staples.
    if(verbosity_lv > 4) printf("MPI%02d:\t\tcomputed Q  \n",
            devinfo.myrank);
    if(verbosity_lv > 5) {// printing stuff
#pragma acc update host(QA[ 0:8])
        printf("-------------Q------------------\n");
        printf("Q00 = %f\n",QA[0].ic00[0]);
        printf("Q00 = %f\n",QA[0].ic11[0]);
        printf("Q01 = %f + (%f)*I\n",crealf(QA[0].c01[0]),cimagf(QA[0].c01[0]));
        printf("Q02 = %f + (%f)*I\n",crealf(QA[0].c02[0]),cimagf(QA[0].c02[0]));
        printf("Q12 = %f + (%f)*I\n\n",crealf(QA[0].c12[0]),cimagf(QA[0].c12[0]));
    }
    compute_lambda_f(Lambda,Sigma,U,QA,TMP);

#ifdef MULTIDEVICE
    communicate_thmat_soa_borders_f(Lambda,1);
#endif

    if(verbosity_lv > 4)   printf("MPI%02d:\t\tcomputed Lambda  \n",
            devinfo.myrank);

    if(verbosity_lv > 5) {// printing stuff
#pragma acc update host(Lambda[0:8])
        printf("-------------LAMBDA------------------\n");
        printf("Lambda00 = %f\n",Lambda[0].rc00[0]);
        printf("Lambda00 = %f\n",Lambda[0].rc11[0]);
        printf("Lambda01 = %f + (%f)*I\n",crealf(Lambda[0].c01[0]),cimagf(Lambda[0].c01[0]));
        printf("Lambda02 = %f + (%f)*I\n",crealf(Lambda[0].c02[0]),cimagf(Lambda[0].c02[0]));
        printf("Lambda12 = %f + (%f)*I\n\n",crealf(Lambda[0].c12[0]),cimagf(Lambda[0].c12[0]));
    }
    compute_sigma_f(Lambda,U,Sigma,QA,TMP,istopo);
    if(verbosity_lv > 4)   printf("MPI%02d:\t\tcomputed Sigma  \n",
            devinfo.myrank);

    if(verbosity_lv > 5) {// printing stuff
#pragma acc update host(Sigma[0:8])
        printf("-------------Sigma[new]------------------\n");                                                                                             
        printf("Sigma[new]00 = %f + (%f)*I\n",crealf(Sigma[0].r0.c0[0]),cimagf(Sigma[0].r0.c0[0]));                                               
        printf("Sigma[new]01 = %f + (%f)*I\n",crealf(Sigma[0].r0.c1[0]),cimagf(Sigma[0].r0.c1[0]));                                               
        printf("Sigma[new]02 = %f + (%f)*I\n",crealf(Sigma[0].r0.c2[0]),cimagf(Sigma[0].r0.c2[0]));                                               
        printf("Sigma[new]10 = %f + (%f)*I\n",crealf(Sigma[0].r1.c0[0]),cimagf(Sigma[0].r1.c0[0]));                                               
        printf("Sigma[new]11 = %f + (%f)*I\n",crealf(Sigma[0].r1.c1[0]),cimagf(Sigma[0].r1.c1[0]));                                               
        printf("Sigma[new]12 = %f + (%f)*I\n",crealf(Sigma[0].r1.c2[0]),cimagf(Sigma[0].r1.c2[0]));                                               
        printf("Sigma[new]20 = %f + (%f)*I\n",crealf(Sigma[0].r2.c0[0]),cimagf(Sigma[0].r2.c0[0]));                                               
        printf("Sigma[new]21 = %f + (%f)*I\n",crealf(Sigma[0].r2.c1[0]),cimagf(Sigma[0].r2.c1[0]));                                               
        printf("Sigma[new]22 = %f + (%f)*I\n\n",crealf(Sigma[0].r2.c2[0]),cimagf(Sigma[0].r2.c2[0]));                
    }

#ifdef MULTIDEVICE
        communicate_gl3_borders_f(Sigma,1);
#endif

}



void fermion_force_soloopenacc_f(__restrict su3_soa_f    * tconf_acc,
#ifdef STOUT_FERMIONS        
			       __restrict su3_soa_f * tstout_conf_acc_arr,// parking
#endif
			       __restrict su3_soa_f * gl3_aux, // gl(3) parking
			       __restrict tamat_soa_f  * tipdot_acc,
			       __restrict ferm_param * tfermion_parameters,// [nflavs] 
			       int tNDiffFlavs,
			       __restrict const vec3_soa_f * ferm_in_acc, // [NPS_tot]         
			       float res,
			       __restrict su3_soa_f  * taux_conf_acc,
			       __restrict vec3_soa_f * tferm_shiftmulti_acc,//parking variable [max_ps*max_approx_order]           
                   inverter_package ipt,
                   const int max_cg )
{
    if(verbosity_lv > 3) printf("SINGLE PRECISION VERSION OF FERMION_FORCE_SOLOOPENACC\n");

    if(verbosity_lv > 2){
        printf("MPI%02d:\tCalculation of fermion force...\n", 
                devinfo.myrank);
    }

#ifdef TIMING_FERMION_FORCE
    struct timeval t1,t2;
    gettimeofday ( &t1, NULL );
#endif

    __restrict su3_soa_f * conf_to_use; // CONF TO USE IN CALCULATION OF 
    // FERMION FORCE
#ifdef STOUT_FERMIONS
    stout_wrapper_f(tconf_acc,tstout_conf_acc_arr,0);// calcolo 
    if(act_params.stout_steps > 0) 
        conf_to_use =  
            &(tstout_conf_acc_arr[8*(act_params.stout_steps-1)]);
    else conf_to_use = tconf_acc;
#else
    conf_to_use = tconf_acc;
#endif
    set_su3_soa_to_zero_f(gl3_aux); // pseudo ipdot
    set_tamat_soa_to_zero_f(tipdot_acc);

    for(int iflav = 0; iflav < tNDiffFlavs; iflav++) {
        set_su3_soa_to_zero_f(taux_conf_acc);
        int ifps = tfermion_parameters[iflav].index_of_the_first_ps;
        for(int ips = 0 ; ips < tfermion_parameters[iflav].number_of_ps ; ips++){


            // modified from converted 
            int converged, cg_return;

            if(1==md_parameters.recycleInvsForce && nMdInversionPerformed >= 2){
                printf("ERROR, not implemented correctly! %s : %d",__FILE__,__LINE__); exit(1);
                int fshift_index = tfermion_parameters[iflav].index_of_the_first_shift;
                int md_approx_order = tfermion_parameters[iflav].approx_md.approx_order;
                
                int ishift;
                for(ishift =0; ishift < md_approx_order; ishift++){
                    double shift = tfermion_parameters[iflav].approx_md.RA_b[ishift];
                    int shiftindex = fshift_index + ips * md_approx_order + ishift;
                    ker_invert_openacc_f(conf_to_use,&tfermion_parameters[iflav],
                            &tferm_shiftmulti_acc[shiftindex],
                            &ferm_in_acc[ifps+ips], res,
                            ipt.loc_r_f, ipt.loc_h_f, ipt.loc_s_f, ipt.loc_p_f,
                            max_cg, shift, &cg_return);
                }
            }else converged = multishift_invert_f(conf_to_use,&tfermion_parameters[iflav],
                     &(tfermion_parameters[iflav].approx_md),
                     tferm_shiftmulti_acc, &(ferm_in_acc[ifps+ips]), res,
                     ipt.loc_r_f, ipt.loc_h_f, ipt.loc_s_f, ipt.loc_p_f,
                     ipt.ferm_shift_temp_f, max_cg, &cg_return);
            convergence_messages(CONVERGENCE_NONCRITICAL, converged);

            ker_openacc_compute_fermion_force_f(conf_to_use, taux_conf_acc, tferm_shiftmulti_acc, ipt.loc_s_f, ipt.loc_h_f, &(tfermion_parameters[iflav]));

        }

        // JUST MULTIPLY BY STAGGERED PHASES,
        // BACK FIELD AND/OR CHEMICAL POTENTIAL 
        multiply_backfield_times_force_f(&(tfermion_parameters[iflav]),taux_conf_acc,gl3_aux);
        if(md_dbg_print_count<debug_settings.md_dbg_print_max_count){
            char taux_conf_acc_name[50];
            sprintf(taux_conf_acc_name,
                    "taux_conf_acc_%s_%d_%d",tfermion_parameters[iflav].name,
                    devinfo.myrank, md_dbg_print_count);
            dbg_print_su3_soa_f(taux_conf_acc,taux_conf_acc_name, 1);
        }


    }
    nMdInversionPerformed++;
#ifdef STOUT_FERMIONS

    for(int stout_level = act_params.stout_steps ; stout_level > 1 ; 
            stout_level--){
        if(verbosity_lv > 1) 
            printf("MPI%02d:\t\tSigma' to Sigma [lvl %d to lvl %d]\n",
                    devinfo.myrank, stout_level,stout_level-1);
        conf_to_use = &(tstout_conf_acc_arr[8*(stout_level-2)]);
        compute_sigma_from_sigma_prime_backinto_sigma_prime_f(gl3_aux,
					           aux_th_f,aux_ta_f,conf_to_use, taux_conf_acc, 0 );
        if(md_dbg_print_count<debug_settings.md_dbg_print_max_count){
            char gl3_aux_name[50];
            sprintf(gl3_aux_name,
                    "gl3_aux_name_%dstout_%d_%d",stout_level,
                    devinfo.myrank, md_dbg_print_count);
            dbg_print_su3_soa_f(gl3_aux,gl3_aux_name,1);
        }

    }
    if(act_params.stout_steps > 0 ){
    if(verbosity_lv > 1) 
        printf("MPI%02d:\t\tSigma' to Sigma [lvl 1 to lvl 0]\n",
                devinfo.myrank);
    compute_sigma_from_sigma_prime_backinto_sigma_prime_f(gl3_aux,
		           		aux_th_f,aux_ta_f,tconf_acc, taux_conf_acc, 0 );
    }
#endif




    multiply_conf_times_force_and_take_ta_nophase_f(tconf_acc, gl3_aux,
            tipdot_acc);

    /*
#pragma acc update host(tipdot_acc[0:8])
printf("-------------FFORCE------------------\n");
printf("F00 = %f\n",tipdot_acc[0].rc00[0]);
printf("F11 = %f\n",tipdot_acc[0].rc11[0]);
printf("F01 = %f + (%f)*I\n",crealf(tipdot_acc[0].c01[0]),cimagf(tipdot_acc[0].c01[0]));
printf("F02 = %f + (%f)*I\n",crealf(tipdot_acc[0].c02[0]),cimagf(tipdot_acc[0].c02[0]));
printf("F12 = %f + (%f)*I\n\n",crealf(tipdot_acc[0].c12[0]),cimagf(tipdot_acc[0].c12[0]));

*/

#ifdef TIMING_FERMION_FORCE
    gettimeofday ( &t2, NULL );
    float dt_preker_to_postker = (float)(t2.tv_sec - t1.tv_sec) + ((float)(t2.tv_usec - t1.tv_usec)/1.0e6f);
    printf("MPI%02d\t\t\
FULL FERMION FORCE COMPUTATION  PreKer->PostKer :%f sec  \n",
dt_preker_to_postker,devinfo.myrank);
#endif
    if(verbosity_lv > 0){
        printf("MPI%02d:\t\tCompleted fermion force openacc\n",
                devinfo.myrank);
    }

    if(debug_settings.save_diagnostics == 1 ){


        float  force_norm, diff_force_norm;
        if((md_diag_count_fermion % debug_settings.md_diag_print_every) == 0){
            ipdot_f_reset = 0;
            copy_ipdot_into_old_f(tipdot_acc,ipdot_f_old_f);
        }


        if((md_diag_count_fermion % debug_settings.md_diag_print_every) == 1 &&
                ipdot_f_reset == 0){
            force_norm = calc_force_norm_f(tipdot_acc);
            diff_force_norm = calc_diff_force_norm_f(tipdot_acc,ipdot_f_old_f);

            if(0 == devinfo.myrank){

                FILE *foutfile = 
                    fopen(debug_settings.diagnostics_filename,"at");
                fprintf(foutfile,"%d\tFFHN %e\n%d\tDFFHN %e\n",
                        md_diag_count_fermion, force_norm,
                        md_diag_count_fermion, diff_force_norm);
                fclose(foutfile);

                if(verbosity_lv > 1)
                    printf("MPI%02d:\
                            \t\t\tFermion Force Half Norm: %e, Diff with previous:%e\n",
                            devinfo.myrank, force_norm, diff_force_norm);
            }
        }

        md_diag_count_fermion++;
    } 


}


#endif
