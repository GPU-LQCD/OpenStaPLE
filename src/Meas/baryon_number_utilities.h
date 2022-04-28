#ifndef BARYON_NUMBER_UTILITIES_H_
#define BARYON_NUMBER_UTILITIES_H_

#include "../OpenAcc/struct_c_def.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif

// FIRST DERIVATIVES

void (*dM_dmu_eo[4])( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield);
void (*dM_dmu_oe[4]) ( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield);
   
// SECOND DERIVATIVES

void (*d2M_dmu2_eo[4])( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield);
void (*d2M_dmu2_oe[4]) ( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield);
    

#endif
