#ifndef FIELD_TIMES_FERMION_MATRIX_H
#define FIELD_TIMES_FERMION_MATRIX_H


#include "./struct_c_def.h"
// #ifndef __GNUC__
#include "openacc.h"
// #endif
#include "./fermionic_utilities.h"
#include "../Include/fermion_parameters.h"


// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
#define __restrict
#endif


// 'wf = with a field (for magnetic susceptibility)'
void acc_Deo_wf_unsafe(__restrict const su3_soa * const u, 
											 __restrict vec3_soa * const out, 
											 __restrict const vec3_soa * const in,
											 const double_soa * phases,
											 __restrict const double_soa* field_re,    // e.g. e^(i phi_1 ) - e^(i phi_2),
											 __restrict const double_soa* field_im);   // or i dphi/dbz  

// 'wf = with a field (for magnetic susceptibility)'
void acc_Doe_wf_unsafe(__restrict const su3_soa * const u,
											 __restrict vec3_soa * const out,
											 __restrict const vec3_soa * const in,
											 const double_soa * phases,
											 __restrict const double_soa* field_re,    // e.g. e^(i phi_1 ) - e^(i phi_2),
											 __restrict const double_soa* field_im);   // or i dphi/dbz                   

// 'wf = with a field (for magnetic susceptibility)'
void acc_Deo_wf(__restrict const su3_soa * const u, 
								__restrict vec3_soa * const out, 
								__restrict const vec3_soa * const in,
								const double_soa * phases,
								__restrict const double_soa* field_re,    // e.g. e^(i phi_1 ) - e^(i phi_2),
								__restrict const double_soa* field_im);   // or i dphi/dbz                    

// 'wf = with a field (for magnetic susceptibility)'
void acc_Doe_wf(__restrict const su3_soa * const u,
								__restrict vec3_soa * const out,
								__restrict const vec3_soa * const in,
								const double_soa * phases,
								__restrict const double_soa* field_re,    // e.g. e^(i phi_1 ) - e^(i phi_2),
								__restrict const double_soa* field_im);   // or i dphi/dbz   


#endif
