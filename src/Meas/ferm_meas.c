// use z2 noise instead of gaussian noise (see hep-lat/9308015)
// use the global defined fermions loc_chi, loc_phi, rnd_o, rnd_e, chi_o and loc_h

#ifndef FERM_MEAS_C 
#define FERM_MEAS_C

#define ALIGN 128

#ifdef __GNUC__
#define _POSIX_C_SOURCE 200809L   // not to have warning on posix memalign
#endif



#include <stdio.h>
#include <stdlib.h>

#include "../Include/common_defines.h"
#include "../Include/fermion_parameters.h"
#include "../Include/montecarlo_parameters.h"
#include "../OpenAcc/alloc_vars.h"
#include "../OpenAcc/sp_alloc_vars.h"
#include "../OpenAcc/fermion_matrix.h"
#include "../OpenAcc/field_times_fermion_matrix.h"
#include "../OpenAcc/fermionic_utilities.h"
#include "../OpenAcc/inverter_package.h"
#include "../OpenAcc/inverter_wrappers.h"
#include "../OpenAcc/float_double_conv.h"
#include "../OpenAcc/random_assignement.h"
#include "../OpenAcc/struct_c_def.h"
#include "./baryon_number_utilities.h"
#include "../OpenAcc/alloc_settings.h"
#include "./ferm_meas.h"
#include "./magnetic_susceptibility_utilities.h"
#include "../Include/inverter_tricks.h"

#ifdef STOUT_FERMIONS
#include "../OpenAcc/stouting.h"
#endif
#include "../OpenAcc/action.h"

#ifdef MULTIDEVICE
#include "../Mpi/communications.h"
#endif
#include "../Mpi/multidev.h"

ferm_meas_params fm_par;

// vedi tesi LS F.Negro per ragguagli (Appendici)
void eo_inversion(inverter_package ip,
        ferm_param * tfermions_parameters,
        double res, int max_cg,
        vec3_soa * in_e,     // z2 noise
        vec3_soa * in_o,     // z2 noise
        vec3_soa * out_e,
        vec3_soa * out_o,
        vec3_soa * phi_e,    // parking variable
        vec3_soa * phi_o){   // parking variable for the inverter

    acc_Deo(ip.u, phi_e, in_o,tfermions_parameters->phases);


    combine_in1_x_fact1_minus_in2_back_into_in2(in_e, tfermions_parameters->ferm_mass , phi_e);
    set_vec3_soa_to_zero(out_e);
    inverter_wrapper(ip,tfermions_parameters,
            out_e,phi_e,res,max_cg,0,CONVERGENCE_CRITICAL);
    acc_Doe(ip.u, phi_o, out_e,tfermions_parameters->phases);
    combine_in1_minus_in2_allxfact(in_o,phi_o,(double)1/tfermions_parameters->ferm_mass,out_o);


}// end eo_inversion




void set_fermion_file_header(ferm_meas_params * fmpar, ferm_param * tferm_par){

    strcpy(fmpar->fermionic_outfile_header,"#0.conf\t1.icopy\t");
    int col_count=2;
    if(1 == fmpar->printPlaqAndRect){
        char strtocat[200];
        sprintf(strtocat, "%02d.Plaquette\t%02d.Rectangle\t",col_count,col_count+1);
        strcat(fmpar->fermionic_outfile_header,strtocat);
        col_count +=2;
    }

    for(int iflv=0;iflv<alloc_info.NDiffFlavs;iflv++){
        char strtocat[200];
        sprintf(strtocat, "%02d.Reff_%-18s%02d.Imff_%-18s",col_count,
                tferm_par[iflv].name,col_count+1,tferm_par[iflv].name);
        strcat(fmpar->fermionic_outfile_header,strtocat);
        col_count +=2;

        sprintf(strtocat, "%02d.ReN_%-19s%02d.ImN_%-19s",col_count,
                tferm_par[iflv].name,col_count+1,tferm_par[iflv].name);
        strcat(fmpar->fermionic_outfile_header,strtocat);
        col_count +=2;

        sprintf(strtocat, "%02d.ReMag_%-18s%02d.ImMag_%-18s",col_count,
                tferm_par[iflv].name,col_count+1,tferm_par[iflv].name);
        strcat(fmpar->fermionic_outfile_header,strtocat);
        col_count +=2;

        if (fmpar->DoubleInvNVectorsChiral>0){
            sprintf(strtocat, 
                    "%02d.ReChSuscConn_%-9s%02d.ImChSuscConn_%-9s",
                    col_count,tferm_par[iflv].name,
                    col_count+1,tferm_par[iflv].name);
            strcat(fmpar->fermionic_outfile_header,strtocat);
            col_count +=2;
        }
        if (fmpar->DoubleInvNVectorsQuarkNumber>0){
            // Actually, the first connected piece does not need a second inversion
            sprintf(strtocat, 
                    "%02d.ReQNSuscConn1_%-8s%02d.ImQNSuscConn1_%-8s",
                    col_count,tferm_par[iflv].name,
                    col_count+1,tferm_par[iflv].name);
            strcat(fmpar->fermionic_outfile_header,strtocat);
            col_count +=2;
            sprintf(strtocat, 
                    "%02d.ReQNSuscConn2_%-9s%02d.ImQNSuscConn2_%-9s",
                    col_count,tferm_par[iflv].name,
                    col_count+1,tferm_par[iflv].name);
            strcat(fmpar->fermionic_outfile_header,strtocat);
            col_count +=2;

        }

    }
    strcat(fmpar->fermionic_outfile_header,"\n");

}

void fermion_measures( su3_soa * tconf_acc,
        ferm_param * tfermions_parameters,
        ferm_meas_params * tfm_par,
        double res, int max_cg, int conf_id_iter,
        double plaq, double rect )
{
    vec3_soa * rnd_e,* rnd_o;
    vec3_soa * chi_e,* chi_o; //results of eo_inversion
    vec3_soa * magchi_e,* magchi_o; //results of eo_inversion
    vec3_soa *bnchi_e,*bnchi_o; //for baryon number calculation 
    vec3_soa *chi2_e,*chi2_o; //results of second eo_inversion, 
    // for susceptibilities

    vec3_soa * phi_e,* phi_o; // parking variables for eo_inverter
    vec3_soa * trial_sol;
    su3_soa * conf_to_use;
    su3_soa_f * conf_to_use_f;

#ifdef STOUT_FERMIONS


    if(act_params.stout_steps > 0){
			  stout_wrapper(tconf_acc ,gstout_conf_acc_arr,0);
        conf_to_use = &gstout_conf_acc_arr[8*(act_params.stout_steps-1)];
    }
    else conf_to_use = tconf_acc;
#else
    conf_to_use = tconf_acc;
#endif
    if(inverter_tricks.useMixedPrecision){ 
        conf_to_use_f = conf_acc_f[0];// global variable
        convert_double_to_float_su3_soa(conf_to_use,conf_to_use_f); 
    }


    int allocation_check;
#define ALLOCCHECK(control_int,var)  if(control_int != 0 ) \
    printf("\tMPI%02d: Error in  allocation of %s . \n",devinfo.myrank, #var);\
    else if(verbosity_lv > 2) printf("\tMPI%02d: Allocation of %s : OK , %p\n",\
            devinfo.myrank, #var, var );\
    fflush(stdout);

    allocation_check =  posix_memalign((void **)&rnd_e, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,rnd_e);     
    allocation_check =  posix_memalign((void **)&rnd_o, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,rnd_o);     
    allocation_check =  posix_memalign((void **)&phi_e, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,phi_e);     
    allocation_check =  posix_memalign((void **)&phi_o, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,phi_o);     
    allocation_check =  posix_memalign((void **)&chi_e, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,chi_e);     
    allocation_check =  posix_memalign((void **)&chi_o, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,chi_o);     
    allocation_check =  posix_memalign((void **)&magchi_e, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,magchi_e);     
    allocation_check =  posix_memalign((void **)&magchi_o, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,magchi_o);     
    allocation_check =  posix_memalign((void **)&bnchi_e, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,bnchi_e);     
    allocation_check =  posix_memalign((void **)&bnchi_o, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,bnchi_o);     
    allocation_check =  posix_memalign((void **)&chi2_e, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,chi2_e);     
    allocation_check =  posix_memalign((void **)&chi2_o, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,chi2_o);    

    allocation_check =  posix_memalign((void **)&trial_sol, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,trial_sol);

#undef ALLOCCHECK

    // FILE checks
    // find the actual file size
    FILE *foutfile;
    int fsize;

    if(devinfo.myrank == 0 ){
        foutfile = fopen(tfm_par->fermionic_outfilename,"ab");

        if(foutfile){
            fseek(foutfile, 0L, SEEK_END);
            fsize = ftell(foutfile);
            fseek(foutfile, 0L, SEEK_SET);
            fsize -= ftell(foutfile);

        }else {
            printf("File %s can't be opened for writing. Exiting.\n", fm_par.fermionic_outfilename);
            exit(1);
        }
        fclose(foutfile);// found file size
    }


#pragma acc data create(phi_e[0:1]) create(phi_o[0:1])\
    create(chi_e[0:1]) create(chi_o[0:1])   \
    create(magchi_e[0:1]) create(magchi_o[0:1])   \
    create(bnchi_e[0:1]) create(bnchi_o[0:1])   \
    create(chi2_e[0:1]) create(chi2_o[0:1])   \
    create(rnd_e[0:1]) create(rnd_o[0:1]) create(trial_sol[0:1]) 
    {



        // cycle on copies
        int icopy = mc_params.measures_done;
        while( icopy < tfm_par->SingleInvNVectors &&
                RUN_CONDITION_TERMINATE != mc_params.run_condition ) {
            FILE *foutfile;

            if(devinfo.myrank == 0 ){
                foutfile = fopen(tfm_par->fermionic_outfilename,"at");
                if(fsize == 0){
                    set_fermion_file_header(tfm_par, tfermions_parameters);
                    fprintf(foutfile,"%s",fm_par.fermionic_outfile_header);
                    fsize++;
                }

                if(!foutfile) {
                    printf("File %s can't be opened for writing. Exiting.\n", fm_par.fermionic_outfilename);
                    exit(1);
                }

                fprintf(foutfile,"%d\t%d\t",conf_id_iter,icopy);
                if(1==tfm_par->printPlaqAndRect){
                    fprintf(foutfile,"%.18lf\t%.18lf\t",plaq,rect);
                }
            }        

            struct timeval pre_flavour_cycle, post_flavour_cycle;
            gettimeofday(&pre_flavour_cycle,NULL);
            // cycle on flavours
            for(int iflv=0; iflv < alloc_info.NDiffFlavs ; iflv++){

                if(verbosity_lv > 1 && devinfo.myrank == 0){
                    printf("MPI%02d: Performing %d of %d chiral measures for quark %s",
                            devinfo.myrank,icopy+1,tfm_par->SingleInvNVectors, 
                            tfermions_parameters[iflv].name);

                    printf("(%d for chiral susc, %d for quark number susc).\n",
                            tfm_par->DoubleInvNVectorsChiral,
                            tfm_par->DoubleInvNVectorsQuarkNumber);

                }


                generate_vec3_soa_z2noise(rnd_e);
                generate_vec3_soa_z2noise(rnd_o);
                generate_vec3_soa_gauss(trial_sol);
#pragma acc update device(rnd_e[0:1]) 
#pragma acc update device(rnd_o[0:1]) 
#pragma acc update device(trial_sol[0:1]) 
                d_complex chircond_size = 0.0 + 0.0*I;
                d_complex magnetization_size = 0.0 + 0.0*I;
                d_complex barnum_size = 0.0 + 0.0*I; // https://en.wikipedia.org/wiki/P._T._Barnum
                d_complex trMinvSq_size = 0.0 + 0.0*I; // for connected chiral susc
                d_complex trd2M_dmu2_Minv_size = 0.0 + 0.0*I; //for connected baryon susc,1
                d_complex trdM_dmuMinv_sq_size = 0.0 + 0.0*I; //for connected baryon susc,2

                double factor = tfermions_parameters[iflv].degeneracy*0.25/GL_SIZE;
                // preparing inverter_package with global variables
                inverter_package ip;
                setup_inverter_package_dp(&ip,conf_to_use,ferm_shiftmulti_acc,1,kloc_r,kloc_h,
                        kloc_s,kloc_p);

                if(inverter_tricks.useMixedPrecision) 
                    setup_inverter_package_sp(&ip,conf_to_use_f,ferm_shiftmulti_acc_f,1,
                            kloc_r_f,kloc_h_f,kloc_s_f,kloc_p_f,aux1_f); 


                // FIRST INVERSION
                // (chi_e,chi_o) = M^{-1} (rnd_e,rnd_o)
                eo_inversion(ip,&tfermions_parameters[iflv],res, max_cg,   
                        rnd_e,rnd_o,chi_e,chi_o,phi_e,phi_o);



                // CHIRAL CONDENSATE
                chircond_size = scal_prod_global(rnd_e,chi_e)+
                    scal_prod_global(rnd_o,chi_o);

                if(devinfo.myrank == 0)

                    fprintf(foutfile,"%.16lf\t%.16lf\t",
                            creal(chircond_size*factor),
                            cimag(chircond_size*factor));

                // QUARK NUMBER
                // (bnchi_e, * ) = dM/dmu (* , chi_o)
                dM_dmu_eo[geom_par.tmap](conf_to_use,bnchi_e,chi_o,
                        tfermions_parameters[iflv].phases);
                // ( * ,bnchi_o) = dM/dmu (chi_e, * )
                dM_dmu_oe[geom_par.tmap](conf_to_use,bnchi_o,chi_e,
                        tfermions_parameters[iflv].phases);


                // (bnchi_e,bnchi_o) = dM/dmu (chi_e,chi_o) 
                barnum_size = scal_prod_global(rnd_e,bnchi_e)+
                    scal_prod_global(rnd_o,bnchi_o);


                if(devinfo.myrank == 0)
                    fprintf(foutfile,"%.16lf\t%.16lf\t",
                            creal(barnum_size*factor),
                            cimag(barnum_size*factor));

                // MAGNETIZATION


                acc_Deo_wf(conf_to_use,magchi_e,chi_o,tfermions_parameters[iflv].phases,
                        tfermions_parameters[iflv].mag_re,
                        tfermions_parameters[iflv].mag_im);
                acc_Doe_wf(conf_to_use,magchi_o,chi_e,tfermions_parameters[iflv].phases,
                        tfermions_parameters[iflv].mag_re,
                        tfermions_parameters[iflv].mag_im);

                magnetization_size = - scal_prod_global(rnd_e,magchi_e) -
                    scal_prod_global(rnd_o,magchi_o);

                // Notice: 
                // the minus sign comes from the fact that G = U - TS - MH/beta = 1/beta ln Z 
                // (right) 
                // and the results of the scalar product are dG/dH

                if(devinfo.myrank == 0)
                    fprintf(foutfile,"%.16lf\t%.16lf\t",
                            creal(magnetization_size*factor),
                            cimag(magnetization_size*factor));

                // SUSCEPTIBILITIES
                if(icopy < tfm_par->DoubleInvNVectorsChiral ){

                    // CHIRAL SUSCEPTIBILITY
                    // (chi2_e,chi2_o) = M^{-1} (chi_e,chi_o) = 
                    // = M^{-2} (rnd_e, rnd_o)
                    eo_inversion(ip,&tfermions_parameters[iflv],
                            res,max_cg,chi_e,chi_o,chi2_e,chi2_o,phi_e,phi_o);

                    trMinvSq_size = -scal_prod_global(rnd_e,chi2_e)-
                        scal_prod_global(rnd_o,chi2_o); 

                    if(devinfo.myrank == 0)
                        fprintf(foutfile,"%.16lf\t%.16lf\t",
                                creal(trMinvSq_size*factor),
                                cimag(trMinvSq_size*factor));

                } else if (tfm_par->DoubleInvNVectorsChiral>0)                  
                    if(devinfo.myrank == 0)
                        fprintf(foutfile,"%-24s%-24s","none","none" );

                if(icopy < tfm_par->DoubleInvNVectorsQuarkNumber ){
                    // Actually, the first connected piece 
                    // does not need a second inversion

                    // CONNECTED QUARK NUMBER SUSCEPTIBILITY (1),
                    // (dM/dmu)^2 M^{-1}

                    // (chi2_e, * ) = d^2M/dmu^2 (* , chi_o)
                    d2M_dmu2_eo[geom_par.tmap](conf_to_use,chi2_e,chi_o,
                            tfermions_parameters[iflv].phases); 
                    // ( * ,chi2_o) = d^2M/dmu^2 (chi_e, * )      
                    d2M_dmu2_oe[geom_par.tmap](conf_to_use,chi2_o,chi_e,
                            tfermions_parameters[iflv].phases);

                    // (chi2_e,chi2_o) = d^2M/dmu^2 M^{-1} (rnd_e,rnd_o)
                    trd2M_dmu2_Minv_size = scal_prod_global(rnd_e,chi2_e)
                        + scal_prod_global(rnd_o,chi2_o);

                    if(devinfo.myrank == 0)
                        fprintf(foutfile,"%.16lf\t%.16lf\t",
                                creal(trd2M_dmu2_Minv_size*factor),
                                cimag(trd2M_dmu2_Minv_size*factor));

                    // CONNECTED QUARK NUMBER SUSCEPTIBILITY (2),
                    // (dM/dmu M^{-1})^2

                    // (chi2_e,chi2_o) = M^{-1} (bnchi_e,bnchi_o) =
                    // = M^{-1} dM/dmu M^{-1} (rnd_e,rnd_o)
                    eo_inversion(ip,&tfermions_parameters[iflv],
                            res,max_cg,bnchi_e,bnchi_o,chi2_e,chi2_o,
                            phi_e,phi_o);

                    // (bnchi_e, * ) = dM/dmu (* , chi2_o)
                    dM_dmu_eo[geom_par.tmap](conf_to_use,bnchi_e,chi2_o,
                            tfermions_parameters[iflv].phases);
                    // ( * ,bnchi_o) = dM/dmu (chi2_e, * )
                    dM_dmu_oe[geom_par.tmap](conf_to_use,bnchi_o,chi2_e,
                            tfermions_parameters[iflv].phases);
                    // (bnchi_e,bnchi_o) = (dM/dmu M^{-1})^2 (rnd_e,rnd_o)

                    trdM_dmuMinv_sq_size = 
                        -scal_prod_global(rnd_e,bnchi_e)
                        -scal_prod_global(rnd_o,bnchi_o); 

                    if(devinfo.myrank == 0)
                        fprintf(foutfile,"%.16lf\t%.16lf\t",
                                creal(trdM_dmuMinv_sq_size*factor),
                                cimag(trdM_dmuMinv_sq_size*factor));

                } else if (tfm_par->DoubleInvNVectorsQuarkNumber>0)                  

                    if(devinfo.myrank == 0)
                        fprintf(foutfile,"%-24s%-24s%-24s%-24s","none",
                                "none","none","none" );

            }// end of cycle over flavours
            if(devinfo.myrank == 0){
                fprintf(foutfile,"\n");
                fclose(foutfile);
            }
            gettimeofday(&post_flavour_cycle,NULL);
            
            icopy++;
            

            double flavour_cycle_time = 
            (post_flavour_cycle.tv_sec - pre_flavour_cycle.tv_sec)+
                (double)(post_flavour_cycle.tv_usec - pre_flavour_cycle.tv_usec)/1.0e6;

            mc_params.max_flavour_cycle_time = 
                (flavour_cycle_time > mc_params.max_flavour_cycle_time)?
                flavour_cycle_time : mc_params.max_flavour_cycle_time;

            double total_duration = (double) 
                (post_flavour_cycle.tv_sec - mc_params.start_time.tv_sec)+
                (double)(post_flavour_cycle.tv_usec - mc_params.start_time.tv_usec)/1.0e6;

            double max_expected_duration_with_another_cycle = 
                total_duration + 2*mc_params.max_flavour_cycle_time;

            if( 0 == devinfo.myrank  && 
                    max_expected_duration_with_another_cycle > mc_params.MaxRunTimeS){
                printf("Time is running out (%d of %d seconds elapsed),",
                        (int) total_duration, (int) mc_params.MaxRunTimeS);
                printf(" shutting down now.\n");
                printf("Total max expected duration: %d seconds",
                        (int) max_expected_duration_with_another_cycle);
                printf("(%d elapsed now, 2 * %d maximum expected)\n",(int) total_duration,
                        (int) mc_params.max_flavour_cycle_time);
                //https://www.youtube.com/watch?v=MfGhlVcrc8U
                // but without that much pathos
                mc_params.run_condition = RUN_CONDITION_TERMINATE;
            
            }

            if(0 == devinfo.myrank && RUN_CONDITION_TERMINATE != mc_params.run_condition){
                FILE * test_stop = fopen("stop","r");
                if(test_stop){
                    fclose(test_stop);
                    printf("File  \'stop\' found, stopping cycle now.\n");
                    mc_params.run_condition = RUN_CONDITION_TERMINATE;
                }
            }



#ifdef MULTIDEVICE
	    MPI_Bcast((void*)&(mc_params.run_condition),1,MPI_INT,0,devinfo.mpi_comm);
	    printf("MPI%02d - Broadcast of run condition %d from master...\n",
			    devinfo.myrank, mc_params.run_condition);
#endif
            
	    mc_params.measures_done = icopy;

        }// end of cycle over copies


    }

    free(rnd_e);
    free(rnd_o);
    free(phi_e);
    free(phi_o);
    free(bnchi_e);
    free(bnchi_o);
    free(magchi_e);
    free(magchi_o);
    free(chi_e);
    free(chi_o);
    free(chi2_e);
    free(chi2_o);
    free(trial_sol);


}

#endif
