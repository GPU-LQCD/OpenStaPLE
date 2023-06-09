#ifndef PLAQUETTES_H_
#define PLAQUETTES_H_


#include "../OpenAcc/struct_c_def.h"


// routine for the computation of the average of the plaquettes computed on the plane mu-nu
// 1) all the plaquettes on the plane mu-nu are computed and saved locally
// 2) finally the reduction of the traces is performed
// routine to compute the staples for each site on a given plane mu-nu and sum the result to the local stored staples
double calc_loc_plaquettes_nnptrick(__restrict const su3_soa * const u,
																		__restrict su3_soa * const loc_plaq,
																		dcomplex_soa * const tr_local_plaqs,
																		const int mu, const int nu);
void calc_loc_staples_nnptrick_all(__restrict const su3_soa * const u,
																	 __restrict su3_soa * const loc_stap);
void calc_loc_staples_nnptrick_all_onlyferms(__restrict const su3_soa * const u,
																						 __restrict su3_soa * const loc_stap);
void calc_loc_staples_nnptrick_all_only_even(__restrict const su3_soa * const u,
																						 __restrict su3_soa * const loc_stap);
void calc_loc_staples_nnptrick_all_only_odd(__restrict const su3_soa * const u,
																						__restrict su3_soa * const loc_stap);

#if NRANKS_D3 > 1 // #ifdef MULTIDEVICE
void calc_loc_staples_nnptrick_all_bulk(__restrict const su3_soa * const u,
																				__restrict su3_soa * const loc_stap);

void calc_loc_staples_nnptrick_all_d3c(__restrict const su3_soa * const u,
																			 __restrict su3_soa * const loc_stap,
																			 int offset, int thickness);


#endif


#endif
