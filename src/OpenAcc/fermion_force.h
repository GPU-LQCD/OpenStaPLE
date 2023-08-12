#ifndef FERMION_FORCE_H
#define FERMION_FORCE_H


#include "./struct_c_def.h"
#include "../Include/fermion_parameters.h"
#include "./inverter_package.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
#define __restrict
#endif


void compute_sigma_from_sigma_prime_backinto_sigma_prime(__restrict su3_soa    * Sigma, // corresponding global variable: auxbis_conf_acc
																												 __restrict thmat_soa  * Lambda, // corresponding global variable: aux_th
																												 __restrict tamat_soa  * QA, // corresponding global variable: aux_ta
																												 __restrict const su3_soa * const U,// smeared conf
																												 __restrict su3_soa * const TMP,// corresponding global variable: aux_conf_acc
																												 int nnp_openacc[sizeh][4][2], int nnm_openacc[sizeh][4][2],
																												 const int istopo // istopo = {0,1} -> rho = {fermrho,toporho}
																												 );

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
															 int nnp_openacc[sizeh][4][2], int nnm_openacc[sizeh][4][2],
															 inverter_package ipt,
															 const int max_cg );



#endif
