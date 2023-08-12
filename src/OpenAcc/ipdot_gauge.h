#ifndef IPDOT_GAUGE_H
#define IPDOT_GAUGE_H

#include "./struct_c_def.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
#define __restrict
#endif


void calc_ipdot_gauge_soloopenacc_std(__restrict const  su3_soa * const tconf_acc,
																			__restrict su3_soa * const local_staples,
																			__restrict tamat_soa * const tipdot,
																			int nnp_openacc[sizeh][4][2], int nnm_openacc[sizeh][4][2]);

void calc_ipdot_gauge_soloopenacc_tlsm(__restrict const  su3_soa * const tconf_acc,  
																			 __restrict su3_soa * const local_staples,
																			 __restrict tamat_soa * const tipdot,
																			 int nnp_openacc[sizeh][4][2], int nnm_openacc[sizeh][4][2]);


// versatile wrapper which chooses between std gauge action or tlsm gauge action
void calc_ipdot_gauge_soloopenacc(__restrict const  su3_soa * const tconf_acc,  
																	__restrict su3_soa * const local_staples,
																	__restrict tamat_soa * const tipdot,
																	int nnp_openacc[sizeh][4][2], int nnm_openacc[sizeh][4][2]);

#if NRANKS_D3 > 1
void calc_ipdot_gauge_soloopenacc_std_bulk(__restrict const  su3_soa * const tconf_acc, 
																					 __restrict su3_soa * const local_staples,
																					 __restrict tamat_soa * const tipdot,
																					 int nnp_openacc[sizeh][4][2],int nnm_openacc[sizeh][4][2]);

void calc_ipdot_gauge_soloopenacc_tlsm_bulk(__restrict const  su3_soa * const tconf_acc,  
																						__restrict su3_soa * const local_staples,
																						__restrict tamat_soa * const tipdot,
																						int nnp_openacc[sizeh][4][2],int nnm_openacc[sizeh][4][2]);

void calc_ipdot_gauge_soloopenacc_bulk(__restrict const  su3_soa * const tconf_acc,  
																			 __restrict su3_soa * const local_staples,
																			 __restrict tamat_soa * const tipdot);

void calc_ipdot_gauge_soloopenacc_std_d3c(__restrict const  su3_soa * const tconf_acc, 
																					__restrict su3_soa * const local_staples,
																					__restrict tamat_soa * const tipdot,
																					int nnp_openacc[sizeh][4][2],int nnm_openacc[sizeh][4][2],
																					int offset3, int thickness3);

void calc_ipdot_gauge_soloopenacc_tlsm_d3c(__restrict const  su3_soa * const tconf_acc,  
																					 __restrict su3_soa * const local_staples,
																					 __restrict tamat_soa * const tipdot,
																					 int nnp_openacc[sizeh][4][2],int nnm_openacc[sizeh][4][2],
																					 int offset3, int thickness3);

void calc_ipdot_gauge_soloopenacc_d3c(__restrict const  su3_soa * const tconf_acc,  
																			__restrict su3_soa * const local_staples,
																			__restrict tamat_soa * const tipdot,
																			int nnp_openacc[sizeh][4][2],int nnm_openacc[sizeh][4][2],
																			int offset3, int thickness3);

#endif


#endif
