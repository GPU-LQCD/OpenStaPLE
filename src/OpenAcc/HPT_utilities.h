#ifndef HPT_UTILITIES_H
#define HPT_UTILITIES_H

#include "../Include/common_defines.h"
#ifdef PAR_TEMP

#include <stdio.h>
#include "./struct_c_def.h"
#include "./alloc_settings.h"
#include "./alloc_vars.h"
#include "./action.h"
#include "../Include/rep_info.h"

typedef struct defect_info_t{
	// min and max values of the logical coordinates where Delta_S_SWAP =/= 0
		
	// Wilson case
	// def_min/max[nu][dirs]. Delta_S_SWAP=/=0 for the 3 nu=/=def_axis. For each nu, store min/max for each dir (NB: 1 of the 4 nu index is not used, defined just for convenience)
	int defect_swap_min[4][4];
	int defect_swap_max[4][4];
		
	// Symanzik case
	// def_min/max[j][nu][dirs]: j specifies which case is considered ( 0 = 1x2 rect, 1 = 2x1 rect), nu and dirs same as Wilson case
#ifdef GAUGE_ACT_TLSM
	int defect_swap_min_TLSM[2][4][4];
	int defect_swap_max_TLSM[2][4][4];
#endif

	// defect location and extension
	int def_axis_mapped; // boundary where the defect is put MAPPED
	int def_mapped_perp_dir[3]; // the 3 logical directions orthogonal to def_axis_mapped

} defect_info;

extern defect_info *def;

void counter_size_function(int d0,int d1,int d2,int d3);
void init_k(su3_soa * conf,double c_r,int def_axis,int * def_vec,defect_info * def, int defect_info_config);
void printing_k_mu(su3_soa * conf);
void trasl_conf( __restrict const su3_soa *  const tconf_acc,
                 __restrict const su3_soa *  const taux_conf);

//void replicas_swap(su3_soa * conf1,su3_soa * conf2, int lab1, int lab2, rep_info * hpt_params);
void label_print(rep_info * hpt_params, FILE *file, int step_number);

//#ifdef GAUGE_ACT_TLSM
//double calc_Delta_S_Symanzik_SWAP( __restrict const su3_soa * const u,
//                                   __restrict const su3_soa * const w,
//                                   __restrict su3_soa * const loc_rects,
//                                   dcomplex_soa * const tr_local_rects,
//                                   const int mu, const int nu, defect_info * def);
//#endif //GAUGE_ACT_TLSM

#ifdef GAUGE_ACT_TLSM
double calc_S_Symanzik_defect(__restrict const su3_soa * const u,
                              __restrict su3_soa * const loc_rects,
                              dcomplex_soa * const tr_local_rects,
                              const int mu, const int nu, defect_info * def);
#endif //GAUGE_ACT_TLSM

//double calc_Delta_S_Wilson_SWAP(__restrict const su3_soa * const u,
//                                __restrict const su3_soa * const w,
//                                __restrict su3_soa * const loc_plaq,
//                                dcomplex_soa * const tr_local_plaqs, 
//                                const int mu, const int nu, defect_info * def);

double calc_S_Wilson_defect(__restrict const su3_soa * const u,
                              __restrict su3_soa * const loc_plaq,
                              dcomplex_soa * const tr_local_plaqs,
                              const int mu, const int nu, defect_info * def);


//double calc_Delta_S_soloopenacc_SWAP( __restrict  su3_soa * const tconf_acc,
//																			__restrict  su3_soa * const tconf_acc2,
//																			__restrict su3_soa * const local_plaqs,
//																			dcomplex_soa * const tr_local_plaqs,defect_info * def);

double calc_S_soloopenacc_defect(__restrict  su3_soa * const tconf_acc,
																 __restrict su3_soa * const local_plaqs,
																 dcomplex_soa * const tr_local_plaqs,defect_info * def);


//int metro_SWAP(su3_soa ** conf_acc,
//               __restrict su3_soa * const loc_plaq,
//               dcomplex_soa * const tr_local_plaqs,
//               int rep_indx1, int rep_indx2,defect_info * def, rep_info *hpt_params);

int metro_SWAP_worldmaster(double Delta_S_SWAP);


//void All_Conf_SWAP(su3_soa ** conf_acc,
//                   __restrict su3_soa * const loc_plaq,
//                   dcomplex_soa * const tr_local_plaqs,                   
//                   defect_info * def, 
//									 int *swap_num,int * all_swap_vet, int * acceptance_vet, rep_info *hpt_params);
void compute_S_of_replicas(
                        __restrict su3_soa * const tconf_acc, 
                        __restrict su3_soa * const local_plaq, 
                        dcomplex_soa * const tr_local_plaqs,
                        defect_info * def,
                        double * S_arr);


#endif // PAR_TEMP
#endif // HPT_UTILITIES_H
