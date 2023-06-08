#ifndef FERMION_FORCE_C
#define FERMION_FORCE_C


#include "../DbgTools/dbgtools.h"
#include "../Include/common_defines.h"
#include "../Include/debug.h"
#include "../Include/fermion_parameters.h"
#include "../Mpi/multidev.h"
#include "./action.h"
#include "./alloc_vars.h"
#include "./backfield.h"
#include "./fermion_force.h"
#include "./fermion_force_utilities.h"
#include "./inverter_multishift_full.h"
#include "../Include/inverter_tricks.h"
#include "./inverter_package.h"
#include "./md_parameters.h"
#include "./plaquettes.h"
#include "./stouting.h"
#include "./struct_c_def.h"
#include "./su3_measurements.h"
#include "./su3_utilities.h"
#include "./inverter_wrappers.h"

// includes related to single precision acceleration
// global variable used for convenience

#include "./alloc_vars.h"
#include "./float_double_conv.h"
#include "./sp_alloc_vars.h"
#include "./sp_struct_c_def.h"



#ifndef __GNUC__
#define TIMING_FERMION_FORCE
#endif

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
#define __restrict
#endif

#ifdef MULTIDEVICE
#include "../Mpi/communications.h"
#endif 

extern int verbosity_lv;


void compute_sigma_from_sigma_prime_backinto_sigma_prime(  __restrict su3_soa    * Sigma, // corresponding global variable: auxbis_conf_acc
																													 __restrict thmat_soa  * Lambda, // corresponding global variable: aux_th
																													 __restrict tamat_soa  * QA, // corresponding global variable: aux_ta
																													 __restrict const su3_soa * const U,// smeared conf
																													 __restrict su3_soa * const TMP,// corresponding global variable: aux_conf_acc
																													 const int istopo // istopo = {0,1} -> rho = {fermrho,toporho}
																													 ){
	if(verbosity_lv > 3) printf("DOUBLE PRECISION VERSION OF COMPUTE_SIGMA_FROM_SIGMA_PRIME_BACKINTO_SIGMA_PRIME\n");
    



	if(verbosity_lv > 2) printf("MPI%02d:\t\tSIGMA_PRIME --> SIGMA\n",
															devinfo.myrank);
	if(verbosity_lv > 5){ // printing stuff
		#pragma acc update host(Sigma[0:8])
		printf("-------------Sigma[old]------------------\n");                                                                                             
		printf("Sigma[old]00 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r0.c0[0]),cimag(Sigma[0].r0.c0[0]));                                               
		printf("Sigma[old]01 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r0.c1[0]),cimag(Sigma[0].r0.c1[0]));                                               
		printf("Sigma[old]02 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r0.c2[0]),cimag(Sigma[0].r0.c2[0]));                                               
		printf("Sigma[old]10 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r1.c0[0]),cimag(Sigma[0].r1.c0[0]));                                               
		printf("Sigma[old]11 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r1.c1[0]),cimag(Sigma[0].r1.c1[0]));                                               
		printf("Sigma[old]12 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r1.c2[0]),cimag(Sigma[0].r1.c2[0]));                                               
		printf("Sigma[old]20 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r2.c0[0]),cimag(Sigma[0].r2.c0[0]));                                               
		printf("Sigma[old]21 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r2.c1[0]),cimag(Sigma[0].r2.c1[0]));                                               
		printf("Sigma[old]22 = %.18lf + (%.18lf)*I\n\n",creal(Sigma[0].r2.c2[0]),cimag(Sigma[0].r2.c2[0]));                
    
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



	set_su3_soa_to_zero(TMP);

	calc_loc_staples_nnptrick_all(U,TMP);
	if(verbosity_lv > 4)printf("MPI%02d:\t\tcomputed staples  \n",
														 devinfo.myrank);
#ifdef MULTIDEVICE
	communicate_gl3_borders(TMP,1);
#endif
	RHO_times_conf_times_staples_ta_part(U,TMP,QA,istopo);

#ifdef MULTIDEVICE
	communicate_tamat_soa_borders(QA,1);
#endif

	// check: TMP = local staples.
	if(verbosity_lv > 4) printf("MPI%02d:\t\tcomputed Q  \n",
															devinfo.myrank);
	if(verbosity_lv > 5) { // printing stuff
		#pragma acc update host(QA[ 0:8])
		printf("-------------Q------------------\n");
		printf("Q00 = %.18lf\n",QA[0].ic00[0]);
		printf("Q00 = %.18lf\n",QA[0].ic11[0]);
		printf("Q01 = %.18lf + (%.18lf)*I\n",creal(QA[0].c01[0]),cimag(QA[0].c01[0]));
		printf("Q02 = %.18lf + (%.18lf)*I\n",creal(QA[0].c02[0]),cimag(QA[0].c02[0]));
		printf("Q12 = %.18lf + (%.18lf)*I\n\n",creal(QA[0].c12[0]),cimag(QA[0].c12[0]));
	}
	compute_lambda(Lambda,Sigma,U,QA,TMP);

#ifdef MULTIDEVICE
	communicate_thmat_soa_borders(Lambda,1);
#endif

	if(verbosity_lv > 4)   printf("MPI%02d:\t\tcomputed Lambda  \n",
																devinfo.myrank);

	if(verbosity_lv > 5) { // printing stuff
		#pragma acc update host(Lambda[0:8])
		printf("-------------LAMBDA------------------\n");
		printf("Lambda00 = %.18lf\n",Lambda[0].rc00[0]);
		printf("Lambda00 = %.18lf\n",Lambda[0].rc11[0]);
		printf("Lambda01 = %.18lf + (%.18lf)*I\n",creal(Lambda[0].c01[0]),cimag(Lambda[0].c01[0]));
		printf("Lambda02 = %.18lf + (%.18lf)*I\n",creal(Lambda[0].c02[0]),cimag(Lambda[0].c02[0]));
		printf("Lambda12 = %.18lf + (%.18lf)*I\n\n",creal(Lambda[0].c12[0]),cimag(Lambda[0].c12[0]));
	}
	compute_sigma(Lambda,U,Sigma,QA,TMP,istopo);
	if(verbosity_lv > 4)   printf("MPI%02d:\t\tcomputed Sigma  \n",
																devinfo.myrank);

	if(verbosity_lv > 5) { // printing stuff
		#pragma acc update host(Sigma[0:8])
		printf("-------------Sigma[new]------------------\n");                                                                                             
		printf("Sigma[new]00 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r0.c0[0]),cimag(Sigma[0].r0.c0[0]));                                               
		printf("Sigma[new]01 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r0.c1[0]),cimag(Sigma[0].r0.c1[0]));                                               
		printf("Sigma[new]02 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r0.c2[0]),cimag(Sigma[0].r0.c2[0]));                                               
		printf("Sigma[new]10 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r1.c0[0]),cimag(Sigma[0].r1.c0[0]));                                               
		printf("Sigma[new]11 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r1.c1[0]),cimag(Sigma[0].r1.c1[0]));                                               
		printf("Sigma[new]12 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r1.c2[0]),cimag(Sigma[0].r1.c2[0]));                                               
		printf("Sigma[new]20 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r2.c0[0]),cimag(Sigma[0].r2.c0[0]));                                               
		printf("Sigma[new]21 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r2.c1[0]),cimag(Sigma[0].r2.c1[0]));                                               
		printf("Sigma[new]22 = %.18lf + (%.18lf)*I\n\n",creal(Sigma[0].r2.c2[0]),cimag(Sigma[0].r2.c2[0]));                
	}



#ifdef MULTIDEVICE
	communicate_gl3_borders(Sigma,1);
#endif

}



void fermion_force_soloopenacc(__restrict su3_soa    * tconf_acc,
#ifdef STOUT_FERMIONS        
															 __restrict su3_soa * tstout_conf_acc_arr, // container for global variable stout_conf_acc_arr
#endif
															 __restrict su3_soa * gl3_aux, // gl(3) matrix container
															 __restrict tamat_soa  * tipdot_acc,
															 __restrict ferm_param * tfermion_parameters, // [nflavs] 
															 int tNDiffFlavs,
															 __restrict const vec3_soa * ferm_in_acc, // [NPS_tot]         
															 double res,
															 __restrict su3_soa  * taux_conf_acc,
															 __restrict vec3_soa * tferm_shiftmulti_acc, // container for global variable ferm_shiftmulti_acc [maxNeededShifts]           
															 inverter_package ipt,
															 const int max_cg )
{
	if(verbosity_lv > 3) printf("DOUBLE PRECISION VERSION OF FERMION_FORCE_SOLOOPENACC\n");

	if(verbosity_lv > 2){
		printf("MPI%02d:\tCalculation of fermion force...\n", 
					 devinfo.myrank);
	}

#ifdef TIMING_FERMION_FORCE
	struct timeval t1,t2;
	gettimeofday ( &t1, NULL );
#endif

	__restrict su3_soa * conf_to_use; // conf to use in calculation
	__restrict su3_soa_f * conf_to_use_f; // conf to use in calculation
	// FERMION FORCE
#ifdef STOUT_FERMIONS
	stout_wrapper(tconf_acc,tstout_conf_acc_arr,0); // computation
	if(act_params.stout_steps > 0) 
		conf_to_use =  
			&(tstout_conf_acc_arr[8*(act_params.stout_steps-1)]);
	else conf_to_use = tconf_acc;
#else
	conf_to_use = tconf_acc;
#endif
	set_su3_soa_to_zero(gl3_aux); // pseudo ipdot
	set_tamat_soa_to_zero(tipdot_acc);

	ipt.u = conf_to_use;

	if(1 == inverter_tricks.singlePInvAccelMultiInv || 
		 1 == md_parameters.recycleInvsForce){
		if(0==devinfo.myrank && verbosity_lv >2) 
			printf("Converting gauge conf to single precision...\n");
		conf_to_use_f = conf_acc_f; // using global variable for convenience.
		convert_double_to_float_su3_soa(conf_to_use,conf_to_use_f);
		ipt.u_f = conf_to_use_f;
	}
	else setup_inverter_package_sp(&ipt,0,0,0,0,0,0,0,0); // passed to this function by copy
	// this instruction has effect only here
	// in this scope  

	for(int iflav = 0; iflav < tNDiffFlavs; iflav++) {
		set_su3_soa_to_zero(taux_conf_acc);
		int ifps = tfermion_parameters[iflav].index_of_the_first_ps;
		for(int ips = 0 ; ips < tfermion_parameters[iflav].number_of_ps ; ips++){

			if(1==md_parameters.recycleInvsForce && nMdInversionPerformed >= 2 ){
				printf("ERROR, not implemented correctly! %s : %d",__FILE__,__LINE__); exit(1);
				if(0==devinfo.myrank && verbosity_lv > 3) 
					printf("Recycling old results for force calculation.\n");
				int fshift_index = tfermion_parameters[iflav].index_of_the_first_shift;
				int md_approx_order = tfermion_parameters[iflav].approx_md.approx_order;
                
				int ishift;
				for(ishift =0; ishift < md_approx_order; ishift++){
					double shift = tfermion_parameters[iflav].approx_md.RA_b[ishift];
					int shiftindex = fshift_index + ips * md_approx_order + ishift;
					inverter_wrapper(ipt, &tfermion_parameters[iflav],
													 &tferm_shiftmulti_acc[shiftindex],&ferm_in_acc[ifps+ips],
													 res, max_cg, shift, CONVERGENCE_NONCRITICAL);
				}
			} else inverter_multishift_wrapper(ipt,&tfermion_parameters[iflav],
																				 &tfermion_parameters[iflav].approx_md,
																				 tferm_shiftmulti_acc, &ferm_in_acc[ifps+ips], res, max_cg,
																				 CONVERGENCE_NONCRITICAL);
			ker_openacc_compute_fermion_force(ipt.u, taux_conf_acc, tferm_shiftmulti_acc,
																				ipt.loc_s, ipt.loc_h, &(tfermion_parameters[iflav]));

		}


		// just multiply by staggered phases,
		// back field and/or chemical potential 
		multiply_backfield_times_force(&(tfermion_parameters[iflav]),taux_conf_acc,gl3_aux);
		if(md_dbg_print_count<debug_settings.md_dbg_print_max_count){
			char taux_conf_acc_name[50];
			sprintf(taux_conf_acc_name,
							"taux_conf_acc_%s_%d_%d",tfermion_parameters[iflav].name,
							devinfo.myrank, md_dbg_print_count);
			dbg_print_su3_soa(taux_conf_acc,taux_conf_acc_name, 1);
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
		compute_sigma_from_sigma_prime_backinto_sigma_prime(gl3_aux,
																												aux_th,aux_ta,conf_to_use, taux_conf_acc, 0);
		if(md_dbg_print_count<debug_settings.md_dbg_print_max_count){
			char gl3_aux_name[50];
			sprintf(gl3_aux_name,
							"gl3_aux_name_%dstout_%d_%d",stout_level,
							devinfo.myrank, md_dbg_print_count);
			dbg_print_su3_soa(gl3_aux,gl3_aux_name,1);
		}

	}

	if(act_params.stout_steps > 0 ){
    if(verbosity_lv > 1) 
			printf("MPI%02d:\t\tSigma' to Sigma [lvl 1 to lvl 0]\n",
						 devinfo.myrank);
    compute_sigma_from_sigma_prime_backinto_sigma_prime(gl3_aux,
																												aux_th,aux_ta,tconf_acc, taux_conf_acc, 0);
	}

#endif




	multiply_conf_times_force_and_take_ta_nophase(tconf_acc, gl3_aux,
																								tipdot_acc);

	/*
		#pragma acc update host(tipdot_acc[0:8])
		printf("-------------FFORCE------------------\n");
		printf("F00 = %.18lf\n",tipdot_acc[0].rc00[0]);
		printf("F11 = %.18lf\n",tipdot_acc[0].rc11[0]);
		printf("F01 = %.18lf + (%.18lf)*I\n",creal(tipdot_acc[0].c01[0]),cimag(tipdot_acc[0].c01[0]));
		printf("F02 = %.18lf + (%.18lf)*I\n",creal(tipdot_acc[0].c02[0]),cimag(tipdot_acc[0].c02[0]));
		printf("F12 = %.18lf + (%.18lf)*I\n\n",creal(tipdot_acc[0].c12[0]),cimag(tipdot_acc[0].c12[0]));

	*/

#ifdef TIMING_FERMION_FORCE
	gettimeofday ( &t2, NULL );
	double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
	printf("MPI%02d\t\t\
FULL FERMION FORCE COMPUTATION  PreKer->PostKer :%f sec  \n",
				 devinfo.myrank,dt_preker_to_postker);
#endif
	if(verbosity_lv > 0){
		printf("MPI%02d:\t\tCompleted fermion force openacc\n",
					 devinfo.myrank);
	}


	if(debug_settings.save_diagnostics == 1 ){

		double  force_norm, diff_force_norm;
		if((md_diag_count_fermion % debug_settings.md_diag_print_every) == 0){
			ipdot_f_reset = 0;
			copy_ipdot_into_old(tipdot_acc,ipdot_f_old);
		}
		if((md_diag_count_fermion % debug_settings.md_diag_print_every) == 1 && 
			 ipdot_f_reset == 0){
			force_norm = calc_force_norm(tipdot_acc);
			diff_force_norm = calc_diff_force_norm(tipdot_acc,ipdot_f_old);


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
