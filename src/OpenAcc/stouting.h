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

void stout_isotropic(
        __restrict const su3_soa * const u,               // --> input conf
        __restrict su3_soa * const uprime,          // --> output conf [stouted]
        __restrict su3_soa * const local_staples,   // --> parking variable
        __restrict su3_soa * const auxiliary,       // --> parking variable
        __restrict tamat_soa * const tipdot,       // --> parking variable
				const int istopo); //istopo = {0,1} -> rho={fermrho,toporho}


void compute_lambda(__restrict thmat_soa * const L, // la Lambda --> ouput  (una cosa che serve per calcolare la forza fermionica successiva)
__restrict const su3_soa   * const SP, // Sigma primo --> input (forza fermionica del passo precedente)
__restrict const su3_soa   * const U,    // la configurazione di gauge --> input
__restrict const tamat_soa * const QA, // gli stessi Q che arrivano a Cayley hamilton --> input (sostanzialmente sono rho*ta(staples))
__restrict su3_soa   * const TMP  // variabile di parcheggio
		    );

#if (defined STOUT_FERMIONS) || (defined STOUT_TOPO)
void stout_wrapper(__restrict const su3_soa * const tconf_acc,
									 __restrict su3_soa * tstout_conf_acc_arr, const int istopo);//istopo=1 -> stout steps = topological stout steps; istopo=0 -> stout steps = fermionic stout steps;
#endif

void compute_sigma(__restrict const thmat_soa * const L,  // la Lambda --> ouput  (una cosa che serve per calcolare la forza fermionica successiva)
 __restrict const su3_soa   * const U,  // la configurazione di gauge --> input
 __restrict su3_soa   * const S,  // entra Sigma primo (input: fermforce del passo precedente) ED esce Sigma --> sia input che ouput
 __restrict const tamat_soa * const QA, // gli stessi Q che arrivano a Cayley hamilton --> input (sostanzialmente sono rho*ta(staples))
__restrict su3_soa   * const TMP, // variabile di parcheggio
const int istopo //istopo = {0,1} -> rho={fermrho,toporho}
);


#endif
