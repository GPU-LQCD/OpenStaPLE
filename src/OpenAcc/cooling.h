#ifndef COOLING_C
#define COOLING_C

#include "./struct_c_def.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
#define __restrict
#endif

void compute_cooled_even_links(__restrict su3_soa   * const U,__restrict su3_soa   * const Ucool,  __restrict su3_soa   * const STAP);

void compute_cooled_odd_links(__restrict su3_soa   * const U,__restrict su3_soa   * const Ucool, 
															__restrict su3_soa   * const STAP);

void cool_conf(__restrict su3_soa   * const U,__restrict su3_soa   * const Ucool, 
							 __restrict su3_soa   * const TMP, int nnp_openacc[sizeh][4][2], int nnm_openacc[sizeh][4][2]);


#endif
