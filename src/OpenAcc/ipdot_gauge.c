#ifndef IPDOT_GAUGE_C
#define IPDOT_GAUGE_C

#include "./md_parameters.h"
#include "./struct_c_def.h"
#include "./su3_utilities.h"
#include "./su3_measurements.h"
#include "./plaquettes.h"
#include "./rectangles.h"
#include "./ipdot_gauge.h"
#include "./geometry.h"
#include "../Include/debug.h"
#include "./alloc_vars.h"

#include "sys/time.h"
// #define TIMING_STAPLES

extern int verbosity_lv;

#include "../Mpi/multidev.h"
#include "../DbgTools/dbgtools.h"
#include "./action.h"

// #define DEBUG_FORCE

#ifdef DEBUG_FORCE
#ifdef MULTIDEVICE
#include "../Mpi/communications.h"
#endif
#endif

void calc_ipdot_gauge_soloopenacc_std(__restrict const su3_soa * const tconf_acc, 
																			__restrict su3_soa * const local_staples,
																			__restrict tamat_soa * const tipdot)
{

#ifdef TIMING_STAPLES
	struct timeval t1,t2;
	gettimeofday ( &t1, NULL );
#endif

	set_su3_soa_to_zero(local_staples); // staples = 0
	calc_loc_staples_nnptrick_all(tconf_acc,local_staples); // compute staples = dS/dU
#ifdef PAR_TEMP
	add_defect_coeffs_to_staple(tconf_acc, local_staples); // staple_mu(x) *= k_mu(x) for every link (x,mu)
#endif
	conf_times_staples_ta_part(tconf_acc,local_staples,tipdot); // U * (dS/dU - dS/dU^dagger)/2

#ifdef DEBUG_FORCE
	// NOTA BENE: this debug only works for a 16^4 lattice with the following mapping: x->1, y->2, z->3, t->0
	int d0=0, d1=15, d2=0, d3=D3_HALO;
	int mu=1;
	int idxh = snum_acc(d0, d1, d2, d3);
	int parity = (d0 + d1 + d2 + d3) % 2;
	int dir_link = 2*mu + parity;

	#pragma acc update self(tconf_acc[0:8]) self(local_staples[0:8]) self(tipdot[0:8])
	if(devinfo.myrank == 0){
		printf("DEBUG FORCE WILSON GREP HERE\n");
		printf("check multiplication into tamat:\n");
		STAMPA_DEBUG_SU3_SOA(tconf_acc,dir_link,idxh);
		STAMPA_DEBUG_SU3_SOA(local_staples,dir_link,idxh);
		STAMPA_DEBUG_TAMAT_SOA(tipdot,dir_link,idxh);
	}

#define SQRT_3 1.732050807568877
	// i*Gell-mann matrices as from eq.A.10 of Gattringer - note that T=lambda/2
	single_tamat i_gell_mann_matr[8]={ 
																		{ 0+1*I, 0+0*I, 0+0*I, 0, 0 },
																		{ 1+0*I, 0+0*I, 0+0*I, 0, 0 },
																		{ 0+0*I, 0+0*I, 0+0*I, 1, -1 },
																		{ 0+0*I, 0+1*I, 0+0*I, 0, 0 },
																		{ 0+0*I, 1+0*I, 0+0*I, 0, 0 },
																		{ 0+0*I, 0+0*I, 0+1*I, 0, 0 },
																		{ 0+0*I, 0+0*I, 1+0*I, 0, 0 },
																		{ 0+0*I, 0+0*I, 0+0*I, 1/SQRT_3, 1/SQRT_3 }
	};
	double eps=1.0e-5;

	// store initial link and comp action
	single_su3 sto;
	#pragma acc update self(tconf_acc[0:8])
	if(devinfo.myrank==0){
		single_su3_from_su3_soa(&tconf_acc[dir_link],idxh,&sto);
		rebuild3row(&sto);
	}
	double ori_act = - calc_plaquette_soloopenacc(tconf_acc, aux_conf_acc, local_sums);

	// store derivative
	single_tamat posi={ 0+0*I, 0+0*I, 0+0*I, 0, 0 };
	single_tamat nega={ 0+0*I, 0+0*I, 0+0*I, 0, 0 };

	for(int igen=0;igen<8;igen++){	    
		// prepare increment and change
		single_tamat ba;
		if(devinfo.myrank==0)
			single_tamat_times_scalar_into_tamat(&ba,&i_gell_mann_matr[igen],eps/2+0*I);

		single_su3 exp_mod;
		if(devinfo.myrank==0){
			CH_exponential_antihermitian_nissalike(&exp_mod,&ba); // exp( -i eps lambda/2 )
			rebuild3row(&exp_mod);
		}
		// change -, compute action
		single_su3 shilink;
		if(devinfo.myrank==0){
			single_su3xsu3(&shilink, &exp_mod, &sto);
			single_su3_into_su3_soa(&tconf_acc[dir_link], idxh, &shilink);
		}

		#pragma acc update device(tconf_acc[0:8])
#if NRANKS_D3 > 1
    communicate_su3_borders(tconf_acc, GAUGE_HALO);  
#endif

		double act_minus = - calc_plaquette_soloopenacc(tconf_acc, aux_conf_acc, local_sums);
		// change +, compute action
		if(devinfo.myrank==0){
			gl3_dagger(&exp_mod); // exp( i eps lambda/2 )
			single_su3xsu3(&shilink, &exp_mod, &sto);
			single_su3_into_su3_soa(&tconf_acc[dir_link], idxh, &shilink);
		}

		#pragma acc update device(tconf_acc[0:8])      
#if NRANKS_D3 > 1
    communicate_su3_borders(tconf_acc, GAUGE_HALO);  
#endif
	
		double act_plus = - calc_plaquette_soloopenacc(tconf_acc, aux_conf_acc, local_sums);

		// set back everything
		if(devinfo.myrank==0)
			single_su3_into_su3_soa(&tconf_acc[dir_link], idxh, &sto);

		#pragma acc update device(tconf_acc[0:8])
#if NRANKS_D3 > 1
    communicate_su3_borders(tconf_acc, GAUGE_HALO);  
#endif

		double check_ori_act = - calc_plaquette_soloopenacc(tconf_acc, aux_conf_acc, local_sums);

		if(devinfo.myrank==0){
			printf("ori = %.15lg\n",ori_act);
			printf("ori_check = %.15lg\n",check_ori_act);
			printf("plus = %.15lg\n",act_plus);
			printf("minus = %.15lg\n",act_minus);
		}
		double gr_plus  =  (act_plus  - ori_act)/eps; //   [ S(U+eps) - S(U) ] / eps = dS/dU + O(eps)
		double gr_minus = -(act_minus - ori_act)/eps; // - [ S(U-eps) - S(U) ] / eps = dS/dU + O(eps)
		if(devinfo.myrank==0){		
			single_tamat_times_scalar_add_to_tamat(&posi,&i_gell_mann_matr[igen],gr_plus+0*I);
			single_tamat_times_scalar_add_to_tamat(&nega,&i_gell_mann_matr[igen],gr_minus+0*I);
		}
	}
	// take the average
	single_tamat Numerical_derivative;
	if(devinfo.myrank==0){
		summ_single_tamats_times_scalar(&Numerical_derivative, &posi, &nega, 0.5+0*I);
		STAMPA_DEBUG_SINGLE_TAMAT(Numerical_derivative);
		printf("Are those similar? If yes, the routine works (unless they are all NaN).\n");
	}
	mem_free_core();
	mem_free_extended();
#ifdef MULTIDEVICE
	shutdown_multidev();
#endif
	exit(1);
#endif

	if(md_dbg_print_count < debug_settings.md_dbg_print_max_count 
		 && 1 == debug_settings.md_dbg_be_verbose ){
		char genericfilename[50];
		// staples
		sprintf(genericfilename,"std_staples_%d_%d",
						devinfo.myrank, md_dbg_print_count);
		dbgprint_gl3_soa(local_staples,genericfilename,1000);

		// tipdot
		sprintf(genericfilename,"std_tipdot_staples_%d_%d",
						devinfo.myrank, md_dbg_print_count);
		print_tamat_soa(tipdot,genericfilename);
	}

#ifdef TIMING_STAPLES
	gettimeofday ( &t2, NULL );
	double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
	printf("FULL STAPLES CALC OPENACC                       PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
#endif

}

void calc_ipdot_gauge_soloopenacc_tlsm(__restrict const su3_soa * const tconf_acc,  
																			 __restrict su3_soa * const local_staples,
																			 __restrict tamat_soa * const tipdot)
{

#ifdef TIMING_STAPLES
	struct timeval t1,t2;
	gettimeofday ( &t1, NULL );
#endif

	set_su3_soa_to_zero(local_staples);
	calc_loc_staples_nnptrick_all(tconf_acc,local_staples);

	// this all-toghether version seems to be quite slower
	// calc_loc_improved_staples_typeABC_nnptrick_all(tconf_acc,local_staples);

	calc_loc_improved_staples_typeA_nnptrick_all(tconf_acc,local_staples);
	calc_loc_improved_staples_typeB_nnptrick_all(tconf_acc,local_staples);
	calc_loc_improved_staples_typeC_nnptrick_all(tconf_acc,local_staples);

#ifdef PAR_TEMP
	add_defect_coeffs_to_staple(tconf_acc, local_staples); // staple_mu(x) *= k_mu(x) for every link (x,mu)
#endif
		
	conf_times_staples_ta_part(tconf_acc,local_staples,tipdot);
 
#ifdef DEBUG_FORCE
	// NOTA BENE: this debug only works for a 16^4 lattice with the following mapping: x->1, y->2, z->3, t->0
	int d0=0, d1=15, d2=0, d3=D3_HALO;
	int mu=1;
	int idxh = snum_acc(d0, d1, d2, d3);
	int parity = (d0 + d1 + d2 + d3) % 2;
	int dir_link = 2*mu + parity;

	#pragma acc update self(tconf_acc[0:8]) self(local_staples[0:8]) self(tipdot[0:8])
	if(devinfo.myrank == 0){
		printf("DEBUG FORCE WILSON GREP HERE\n");
		printf("check multiplication into tamat:\n");
		STAMPA_DEBUG_SU3_SOA(tconf_acc,dir_link,idxh);
		STAMPA_DEBUG_SU3_SOA(local_staples,dir_link,idxh);
		STAMPA_DEBUG_TAMAT_SOA(tipdot,dir_link,idxh);
	}

#define SQRT_3 1.732050807568877
	// i*Gell-mann matrices as from eq.A.10 of Gattringer - note that T=lambda/2
	single_tamat i_gell_mann_matr[8]={ 
																		{ 0+1*I, 0+0*I, 0+0*I, 0, 0 },
																		{ 1+0*I, 0+0*I, 0+0*I, 0, 0 },
																		{ 0+0*I, 0+0*I, 0+0*I, 1, -1 },
																		{ 0+0*I, 0+1*I, 0+0*I, 0, 0 },
																		{ 0+0*I, 1+0*I, 0+0*I, 0, 0 },
																		{ 0+0*I, 0+0*I, 0+1*I, 0, 0 },
																		{ 0+0*I, 0+0*I, 1+0*I, 0, 0 },
																		{ 0+0*I, 0+0*I, 0+0*I, 1/SQRT_3, 1/SQRT_3 }
	};
	double eps=1.0e-5;

	// store initial link and comp action
	single_su3 sto;
	#pragma acc update self(tconf_acc[0:8])
	if(devinfo.myrank==0){
		single_su3_from_su3_soa(&tconf_acc[dir_link],idxh,&sto);
		rebuild3row(&sto);
	}
	double ori_act = - C_ZERO * calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
	ori_act -= C_ONE * calc_rettangolo_soloopenacc(tconf_acc,aux_conf_acc,local_sums);

	// store derivative
	single_tamat posi={ 0+0*I, 0+0*I, 0+0*I, 0, 0 };
	single_tamat nega={ 0+0*I, 0+0*I, 0+0*I, 0, 0 };

	for(int igen=0;igen<8;igen++){	    
		// prepare increment and change
		single_tamat ba;
		if(devinfo.myrank==0)
			single_tamat_times_scalar_into_tamat(&ba,&i_gell_mann_matr[igen],eps/2+0*I);

		single_su3 exp_mod;
		if(devinfo.myrank==0){
			CH_exponential_antihermitian_nissalike(&exp_mod,&ba); // exp( -i eps lambda/2 )
			rebuild3row(&exp_mod);
		}
		// change -, compute action
		single_su3 shilink;
		if(devinfo.myrank==0){
			single_su3xsu3(&shilink, &exp_mod, &sto);
			single_su3_into_su3_soa(&tconf_acc[dir_link], idxh, &shilink);
		}

		#pragma acc update device(tconf_acc[0:8])
#if NRANKS_D3 > 1
    communicate_su3_borders(tconf_acc, GAUGE_HALO);  
#endif
      
		double act_minus = - C_ZERO * calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
		act_minus -= C_ONE * calc_rettangolo_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
		// change +, compute action
		if(devinfo.myrank==0){
			gl3_dagger(&exp_mod); // exp( i eps lambda/2 )
			single_su3xsu3(&shilink, &exp_mod, &sto);
			single_su3_into_su3_soa(&tconf_acc[dir_link], idxh, &shilink);
		}

		#pragma acc update device(tconf_acc[0:8])      
#if NRANKS_D3 > 1
    communicate_su3_borders(tconf_acc, GAUGE_HALO);  
#endif
	
		double act_plus = - C_ZERO * calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
		act_plus -= C_ONE * calc_rettangolo_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
		// set back everything
		if(devinfo.myrank==0)
			single_su3_into_su3_soa(&tconf_acc[dir_link], idxh, &sto);

		#pragma acc update device(tconf_acc[0:8])
#if NRANKS_D3 > 1
    communicate_su3_borders(tconf_acc, GAUGE_HALO);  
#endif

		double check_ori_act = - C_ZERO * calc_plaquette_soloopenacc(tconf_acc, aux_conf_acc, local_sums);
		check_ori_act -= C_ONE * calc_rettangolo_soloopenacc(tconf_acc,aux_conf_acc,local_sums);

		if(devinfo.myrank==0){
			printf("ori = %.15lg\n",ori_act);
			printf("ori_check = %.15lg\n",check_ori_act);
			printf("plus = %.15lg\n",act_plus);
			printf("minus = %.15lg\n",act_minus);
		}
		double gr_plus  =  (act_plus  - ori_act)/eps; //   [ S(U+eps) - S(U) ] / eps = dS/dU + O(eps)
		double gr_minus = -(act_minus - ori_act)/eps; // - [ S(U-eps) - S(U) ] / eps = dS/dU + O(eps)
		if(devinfo.myrank==0){		
			single_tamat_times_scalar_add_to_tamat(&posi,&i_gell_mann_matr[igen],gr_plus+0*I);
			single_tamat_times_scalar_add_to_tamat(&nega,&i_gell_mann_matr[igen],gr_minus+0*I);
		}
	}
	// take the average
	single_tamat Numerical_derivative;
	if(devinfo.myrank==0){
		summ_single_tamats_times_scalar(&Numerical_derivative, &posi, &nega, 0.5+0*I);
		STAMPA_DEBUG_SINGLE_TAMAT(Numerical_derivative);
		printf("Are those similar? If yes, the routine works (unless they are all NaN)\n");
	}
	mem_free_core();
	mem_free_extended();
#ifdef MULTIDEVICE
	shutdown_multidev();
#endif
	exit(1);
#endif
  
	if(md_dbg_print_count<debug_settings.md_dbg_print_max_count
		 && 1 == debug_settings.md_dbg_be_verbose ){
		char genericfilename[50];
		sprintf(genericfilename,"impr_staples_%d_%d",
						devinfo.myrank, md_dbg_print_count);
		dbgprint_gl3_soa(local_staples,genericfilename,1000);
		sprintf(genericfilename,"impr_tipdot_staples_%d_%d",
						devinfo.myrank, md_dbg_print_count);
		print_tamat_soa(tipdot,genericfilename);
	}

#ifdef TIMING_STAPLES
	gettimeofday ( &t2, NULL );
	double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
	printf("FULL STAPLES CALC OPENACC                       PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
#endif

}

void calc_ipdot_gauge_soloopenacc(__restrict const su3_soa * const tconf_acc,  
																	__restrict su3_soa * const local_staples,
																	__restrict tamat_soa * const tipdot)
{
	if(GAUGE_ACTION==0){
		calc_ipdot_gauge_soloopenacc_std(tconf_acc,local_staples,tipdot);
	}
	if(GAUGE_ACTION==1){
		calc_ipdot_gauge_soloopenacc_tlsm(tconf_acc,local_staples,tipdot);
	}
	
	if(debug_settings.save_diagnostics == 1){
		double  force_norm, diff_force_norm;
		int printEvery = debug_settings.md_diag_print_every*md_parameters.gauge_scale;
#if NRANKS_D3 > 1
		if(devinfo.nranks != 1 && devinfo.async_comm_gauge)
			printEvery /= md_parameters.gauge_scale;
#endif
		if((md_diag_count_gauge % printEvery) == 0){
			ipdot_g_reset = 0;
			copy_ipdot_into_old(tipdot,ipdot_g_old);
		}

		if((md_diag_count_gauge % printEvery) == 1 && ipdot_g_reset == 0 ){

			force_norm = calc_force_norm(tipdot)*BETA_BY_THREE;
			diff_force_norm = calc_diff_force_norm(tipdot,ipdot_g_old)*BETA_BY_THREE;

			double diff_force_norm_corrected = diff_force_norm;
#if NRANKS_D3 > 1
			if(devinfo.nranks != 1 && devinfo.async_comm_gauge)
				diff_force_norm_corrected /= (2*md_parameters.gauge_scale+1); 
#endif

			if(0 == devinfo.myrank){
				FILE *foutfile = 
					fopen(debug_settings.diagnostics_filename,"at");
				fprintf(foutfile,"%d\tGFHN %e\n%d\tDGFHN %e",
								md_diag_count_gauge, force_norm,
								md_diag_count_gauge, diff_force_norm_corrected);
#if NRANKS_D3 > 1
				if(devinfo.nranks != 1 && devinfo.async_comm_gauge)
					fprintf(foutfile,"  !!! [Not divided by (2ngs+1): %e ]", diff_force_norm);
#endif
				fprintf(foutfile,"\n");

				fclose(foutfile);

				if(verbosity_lv > 1){
					printf("\t\t\tGauge Force Half Norm: %e, Diff with previous: %e \n", 
								 force_norm, diff_force_norm);
#if NRANKS_D3 > 1
					if(devinfo.nranks != 1 && devinfo.async_comm_gauge)
						printf("\t\t\t (Diff with end of previous gauge cycle)");
#endif
				}

			}
		}

		md_diag_count_gauge++;
	} 

}

#if NRANKS_D3 > 1
void calc_ipdot_gauge_soloopenacc_std_bulk(__restrict const su3_soa * const tconf_acc, 
																					 __restrict su3_soa * const local_staples,
																					 __restrict tamat_soa * const tipdot)
{

#ifdef TIMING_STAPLES
	struct timeval t1,t2;
	gettimeofday ( &t1, NULL );
#endif

	set_su3_soa_to_zero_bulk(local_staples);
	calc_loc_staples_nnptrick_all_bulk(tconf_acc,local_staples);

#ifdef PAR_TEMP
	add_defect_coeffs_to_staple_bulk(tconf_acc, local_staples); // staple_mu(x) *= k_mu(x) for every link (x,mu) on the bulk
#endif

	conf_times_staples_ta_part_bulk(tconf_acc,local_staples,tipdot);

	if(md_dbg_print_count<debug_settings.md_dbg_print_max_count
		 && 1 == debug_settings.md_dbg_be_verbose ){
		char genericfilename[50];
		sprintf(genericfilename,"std_staples_%d_%d_bulk",
						devinfo.myrank, md_dbg_print_count);
		dbgprint_gl3_soa(local_staples,genericfilename,1000);
		sprintf(genericfilename,"std_tipdot_staples_%d_%d_bulk",
						devinfo.myrank, md_dbg_print_count);
		print_tamat_soa(tipdot,genericfilename);
	}


#ifdef TIMING_STAPLES
	gettimeofday ( &t2, NULL );
	double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
	printf("BULK STAPLES CALC OPENACC PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
#endif

}

void calc_ipdot_gauge_soloopenacc_tlsm_bulk(__restrict const su3_soa * const tconf_acc,  
																						__restrict su3_soa * const local_staples,
																						__restrict tamat_soa * const tipdot)
{

#ifdef TIMING_STAPLES
	struct timeval t1,t2;
	gettimeofday ( &t1, NULL );
#endif

	set_su3_soa_to_zero_bulk(local_staples);
	calc_loc_staples_nnptrick_all_bulk(tconf_acc,local_staples);

	// this all-toghether version seems to be quite slower
	// calc_loc_improved_staples_typeABC_nnptrick_all(tconf_acc,local_staples);

	calc_loc_improved_staples_typeA_nnptrick_all_bulk(tconf_acc,local_staples);
	calc_loc_improved_staples_typeB_nnptrick_all_bulk(tconf_acc,local_staples);
	calc_loc_improved_staples_typeC_nnptrick_all_bulk(tconf_acc,local_staples);

#ifdef PAR_TEMP
	add_defect_coeffs_to_staple_bulk(tconf_acc, local_staples); // staple_mu(x) *= k_mu(x) for every link (x,mu) on the bulk
#endif

	conf_times_staples_ta_part_bulk(tconf_acc,local_staples,tipdot);

	if(md_dbg_print_count<debug_settings.md_dbg_print_max_count
		 && 1 == debug_settings.md_dbg_be_verbose ){
		char genericfilename[50];
		sprintf(genericfilename,"impr_staples_%d_%d",
						devinfo.myrank, md_dbg_print_count);
		dbgprint_gl3_soa(local_staples,genericfilename,1000);
		sprintf(genericfilename,"impr_tipdot_staples_%d_%d",
						devinfo.myrank, md_dbg_print_count);
		print_tamat_soa(tipdot,genericfilename);
	}
#ifdef TIMING_STAPLES
	gettimeofday ( &t2, NULL );
	double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
	printf("BULK STAPLES CALC OPENACC, PreKer->PostKer   : %f sec  \n",
				 dt_preker_to_postker);
#endif

}

void calc_ipdot_gauge_soloopenacc_bulk(__restrict const su3_soa * const tconf_acc,  
																			 __restrict su3_soa * const local_staples,
																			 __restrict tamat_soa * const tipdot)
{
	if(GAUGE_ACTION==0){
		calc_ipdot_gauge_soloopenacc_std_bulk(tconf_acc,local_staples,tipdot);
	}
	if(GAUGE_ACTION==1){
		calc_ipdot_gauge_soloopenacc_tlsm_bulk(tconf_acc,local_staples,tipdot);
	}

}

void calc_ipdot_gauge_soloopenacc_std_d3c(__restrict const su3_soa * const tconf_acc, 
																					__restrict su3_soa * const local_staples,
																					__restrict tamat_soa * const tipdot,
																					int offset3, int thickness3)
{

#ifdef TIMING_STAPLES
	struct timeval t1,t2;
	gettimeofday ( &t1, NULL );
#endif

	set_su3_soa_to_zero_d3c(local_staples, offset3, thickness3);
	calc_loc_staples_nnptrick_all_d3c(tconf_acc, local_staples, offset3, thickness3);

#ifdef PAR_TEMP
	add_defect_coeffs_to_staple_d3c(tconf_acc, local_staples, offset3, thickness3); // staple_mu(x) *= k_mu(x) for every link (x,mu) on the border
#endif

	conf_times_staples_ta_part_d3c(tconf_acc,local_staples,tipdot, offset3,thickness3);
    
	if(md_dbg_print_count<debug_settings.md_dbg_print_max_count
		 && 1 == debug_settings.md_dbg_be_verbose ){
		char genericfilename[50];
		sprintf(genericfilename,"std_staples_%d_%d_d3c",
						devinfo.myrank, md_dbg_print_count);
		dbgprint_gl3_soa(local_staples,genericfilename,1000);
		sprintf(genericfilename,"std_tipdot_staples_%d_%d_d3c",
						devinfo.myrank, md_dbg_print_count);
		print_tamat_soa(tipdot,genericfilename);
	}
    

#ifdef TIMING_STAPLES
	gettimeofday ( &t2, NULL );
	double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
	printf("D3C (%d-%d) STAPLES CALC OPENACC PreKer->PostKer :%f sec\n",
				 offset3,offset3+thickness3,dt_preker_to_postker);
#endif

}

void calc_ipdot_gauge_soloopenacc_tlsm_d3c(__restrict const su3_soa * const tconf_acc,  
																					 __restrict su3_soa * const local_staples,
																					 __restrict tamat_soa * const tipdot,
																					 int offset3, int thickness3)
{

#ifdef TIMING_STAPLES
	struct timeval t1,t2;
	gettimeofday ( &t1, NULL );
#endif

	set_su3_soa_to_zero_d3c(local_staples,offset3,thickness3);
	calc_loc_staples_nnptrick_all_d3c(tconf_acc,local_staples,offset3,thickness3);

	// this all-toghether version seems to be quite slower
	// calc_loc_improved_staples_typeABC_nnptrick_all(tconf_acc,local_staples);

	calc_loc_improved_staples_typeA_nnptrick_all_d3c(tconf_acc,local_staples,offset3,thickness3);
	calc_loc_improved_staples_typeB_nnptrick_all_d3c(tconf_acc,local_staples,offset3,thickness3);
	calc_loc_improved_staples_typeC_nnptrick_all_d3c(tconf_acc,local_staples,offset3,thickness3);

#ifdef PAR_TEMP
	add_defect_coeffs_to_staple_d3c(tconf_acc, local_staples, offset3, thickness3); // staple_mu(x) *= k_mu(x) for every link (x,mu) on the border
#endif

	conf_times_staples_ta_part_d3c(tconf_acc,local_staples,tipdot,offset3,thickness3);

	if(md_dbg_print_count<debug_settings.md_dbg_print_max_count
		 && 1 == debug_settings.md_dbg_be_verbose ){
		char genericfilename[50];
		sprintf(genericfilename,"impr_staples_%d_%d_d3c",
						devinfo.myrank, md_dbg_print_count);
		dbgprint_gl3_soa(local_staples,genericfilename,1000);
		sprintf(genericfilename,"impr_tipdot_staples_%d_%d_d3c",
						devinfo.myrank, md_dbg_print_count);
		print_tamat_soa(tipdot,genericfilename);
	}



#ifdef TIMING_STAPLES
	gettimeofday ( &t2, NULL );
	double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
	printf("D3C (%d-%d) STAPLES CALC OPENACC,PreKer->PostKer:%f sec\n",
				 offset3,offset3+thickness3,dt_preker_to_postker);
#endif

}

void calc_ipdot_gauge_soloopenacc_d3c( 
																			__restrict const su3_soa * const tconf_acc,  
																			__restrict su3_soa * const local_staples,
																			__restrict tamat_soa * const tipdot,
																			int offset3, int thickness3)
{
	if(GAUGE_ACTION==0){
		calc_ipdot_gauge_soloopenacc_std_d3c(tconf_acc,
																				 local_staples,tipdot,
																				 offset3,thickness3);
	}
	if(GAUGE_ACTION==1){
		calc_ipdot_gauge_soloopenacc_tlsm_d3c(tconf_acc,
																					local_staples,tipdot,
																					offset3,thickness3);
	}

}
#endif


#endif
