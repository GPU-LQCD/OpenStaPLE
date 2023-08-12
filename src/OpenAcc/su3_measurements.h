#ifndef SU3_MEASUREMENTS_H_
#define SU3_MEASUREMENTS_H_

#include "../Include/common_defines.h"
#include "../OpenAcc/struct_c_def.h"


void check_unitarity_device(__restrict const su3_soa * const u, double * max_unitarity_deviation, double * avg_unitarity_deviation);
void check_unitarity_host(__restrict su3_soa * const u, double * max_unitarity_deviation, double * avg_unitarity_deviation);

double calc_momenta_action(const __restrict thmat_soa * const mom,
													 double_soa * tr_local, const int mu);

double  calc_plaquette_soloopenacc(__restrict  su3_soa * const tconf_acc, __restrict su3_soa * const local_plaqs, dcomplex_soa * const tr_local_plaqs,int nnp_openacc[sizeh][4][2]);

double calc_force_norm(const __restrict tamat_soa * tipdot);


double calc_diff_force_norm(const __restrict tamat_soa * tipdot,
                            const __restrict tamat_soa * tipdot_old);
void copy_ipdot_into_old(const __restrict tamat_soa * tipdot,  
												 __restrict tamat_soa * tipdot_old);


void invert_momenta(__restrict thmat_soa * mom);

double calc_sum_momenta_norm(const __restrict thmat_soa * tmomenta,  
														 const __restrict thmat_soa * tmomenta_old);

double calc_diff_su3_soa_norm(const __restrict su3_soa * tconf,  
															const __restrict su3_soa * tconf_old);

void copy_momenta_into_old(const __restrict thmat_soa * tmomenta,  
													 __restrict thmat_soa * tmomenta_old);


#endif
