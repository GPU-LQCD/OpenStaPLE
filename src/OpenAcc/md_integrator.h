#ifndef MD_INTEGRATOR_H
#define MD_INTEGRATOR_H

#include "./struct_c_def.h"
#include "./md_parameters.h"
#include "../Include/fermion_parameters.h"
#include "./inverter_package.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
#define __restrict
#endif

void multistep_2MN_gauge(su3_soa *tconf_acc,su3_soa *local_staples,tamat_soa *tipdot,thmat_soa *tmomenta,int nnp_openacc[sizeh][4][2],int nnm_openacc[sizeh][4][2]);


void multistep_2MN_SOLOOPENACC(tamat_soa * tipdot_acc,
															 su3_soa  * tconf_acc,
#if (defined STOUT_FERMIONS) || (defined STOUT_TOPO)
															 su3_soa  * tstout_conf_acc_arr, // huge parking for stouting
#endif
															 su3_soa  * tauxbis_conf_acc, 
															 su3_soa  * taux_conf_acc,
															 ferm_param * tfermions_parameters, // [nflavs]
															 int tNDiffFlavs,
															 vec3_soa * ferm_in_acc, // [NPS_tot], will be ferm_chi_acc
															 vec3_soa * tferm_shiftmulti_acc, // parking variable [max_ps*max_approx_order]
															 int nnp_openacc[sizeh][4][2],
															 int nnm_openacc[sizeh][4][2],
															 inverter_package ip,
															 thmat_soa * tmomenta,
															 dcomplex_soa * local_sums,
															 double res,
															 const int max_cg);




#endif

