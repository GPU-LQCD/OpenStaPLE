#ifndef STOUTING_H
#define STOUTING_H

#include "./struct_c_def.h"
#include "../Include/common_defines.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
#define __restrict
#endif


void exp_minus_QA_times_conf(__restrict const su3_soa * const tu,
														 __restrict const tamat_soa * QA,
														 __restrict su3_soa * const tu_out,
														 __restrict su3_soa * const exp_aux);

void stout_isotropic(__restrict const su3_soa * const u, // input conf
										 __restrict su3_soa * const uprime, // output conf [stouted]
										 __restrict su3_soa * const local_staples, // parking variable
										 __restrict su3_soa * const auxiliary, // parking variable
										 __restrict tamat_soa * const tipdot, // parking variable
										 int nnp_openacc[sizeh][4][2], int nnm_openacc[sizeh][4][2],
										 const int istopo); // istopo = {0,1} -> rho = {fermrho,toporho}


void compute_lambda(__restrict thmat_soa * const L, // Lambda --> ouput (this is for next step fermion force computation)
										__restrict const su3_soa   * const SP, // Sigma primo --> input (previous step fermion force)
										__restrict const su3_soa   * const U, // gauge configuration --> input
										__restrict const tamat_soa * const QA, // Cayley Hamilton Qs --> input (rho*ta(staples))
										__restrict su3_soa   * const TMP // parking variable
										);

#if (defined STOUT_FERMIONS) || (defined STOUT_TOPO)
void stout_wrapper(__restrict const su3_soa * const tconf_acc,
									 __restrict su3_soa * tstout_conf_acc_arr, int nnp_openacc[sizeh][4][2], int nnm_openacc[sizeh][4][2],const int istopo); // istopo = {0,1} ==> stout steps = {fermionic,topological} stout steps
#endif

void compute_sigma(__restrict const thmat_soa * const L, // Lambda --> ouput  (una cosa che serve per calcolare la forza fermionica successiva)
									 __restrict const su3_soa   * const U, // gauge configuration --> input
									 __restrict su3_soa   * const S,  // in input it is Sigma prime (input: previous step fermforce);in output it is Sigma
									 __restrict const tamat_soa * const QA, // Cayley hamilton Qs --> input (rho*ta(staples))
									 __restrict su3_soa   * const TMP, // parking variable
									 int nnp_openacc[sizeh][4][2], int nnm_openacc[sizeh][4][2],
									 const int istopo // istopo = {0,1} -> rho = {fermrho,toporho}
									 );


#endif
