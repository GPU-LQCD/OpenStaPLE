#ifndef FIELD_TIMES_FERMION_MATRIX_C
#define FIELD_TIMES_FERMION_MATRIX_C

#include "./field_times_fermion_matrix.h"
#include "../Include/common_defines.h"
#include "./geometry.h"
#include "./struct_c_def.h"
#ifndef __GNUC__
#include "openacc.h"
#endif
#include "./fermionic_utilities.h"
#include "./fermion_matrix.h"
#include "./matvecmul.h"

#ifdef MULTIDEVICE
#include "../Mpi/communications.h"
#endif 
#include "../Mpi/multidev.h"


#define DEO_DOE_PREAMBLE																								\
	int d0, d0m, d1m, d2m, d3m, d0p, d1p, d2p, d3p, idxh, matdir,dirindex;\
	vec3 aux=(vec3) {0,0,0};																							\
	d1m = d1 - 1;																													\
	d1m = d1m + (((d1m >> 31) & 0x1) * nd1);															\
	d2m = d2 - 1;																													\
	d2m = d2m + (((d2m >> 31) & 0x1) * nd2);															\
	d3m = d3 - 1;																													\
	d3m = d3m + (((d3m >> 31) & 0x1) * nd3);															\
	d1p = d1 + 1;																													\
	d1p *= (((d1p-nd1) >> 31) & 0x1);																			\
	d2p = d2 + 1;																													\
	d2p *= (((d2p-nd2) >> 31) & 0x1);																			\
	d3p = d3 + 1;																													\
	d3p *= (((d3p-nd3) >> 31) & 0x1);																			\



#define CLOSE4CYCLES }}}}


// 'wf = with a field' (for magnetic susceptibility)
/*****************************************************
 * NOTE: here is treated the case of an Abelian field.
 * The Abelian field will be present in the action in the
 * usual way, that is, as a phase which multiplies the 
 * SU(3) links:
 *
 * e^{i_phase(x)_\mu} U(x)_\mu \delta_{x+\hat{\mu},y} - h.c.
 * 
 * (where h.c. stands for the term where also the delta is 
 * "transposed).
 * The following represent the version of the Doe and Deo 
 * operators, composed of the links (su3 * phase), where the links
 * are multiplied by an external, complex, field.
 *
 * IMPORTANT EXAMPLE
 * They CAN be used to represent, e.g. \partial M / \partial b,
 * BUT the external field added must be IMAGINARY PURE, that is
 * the derivative of I*PHASE, not the derivative of PHASE only.
 * Notice that indeed the derivative would be 
 *
 * i (dphase(x)_\mu / db) e^{i_phase(x)_\mu} U(x)_\mu \delta_{x+\hat{\mu},y} 
 * + i dphase(x)_\mu / db e^{-i_phase(x)_\mu} U^\dagger(x)_\mu \delta_{y+\hat{\mu},x} 
 *
 * THE SIGN is PLUS. The last expression can be obviously rewritten as 
 *
 * i (dphase(x)_\mu / db) e^{i_phase(x)_\mu} U(x)_\mu \delta_{x+\hat{\mu},y} - h.c.
 *
 * The following functions take as input :
 * - the SU(3) links
 * - the abelian phase
 * - the derivative of I*phase with respect to the field parameters (e.g. b_z): this 
 *   is a complex field.
 *
 * *********************************************************************************/
void acc_Deo_wf_unsafe( __restrict const su3_soa * const u, 
												__restrict vec3_soa * const out, 
												__restrict const vec3_soa * const in,
												const double_soa * phases,
												__restrict const double_soa* field_re,   // e.g. e^(i phi_1 ) - e^(i phi_2),
												__restrict const double_soa* field_im)   // or I dphi/dbz : notice the I
{
	int hd0, d1, d2, d3;
	#pragma acc kernels present(u) present(out) present(in) present(phases)	\
		present(field_re) present(field_im)
	#pragma acc loop independent gang(DEODOEGANG3)
	for(d3=D3_HALO; d3<D3_HALO+LOC_N3;d3++) {
		#pragma acc loop independent tile(DEODOETILE0,DEODOETILE1,DEODOETILE2)
		for(d2=0; d2<nd2; d2++) {
			for(d1=0; d1<nd1; d1++) {
				for(hd0=0; hd0 < nd0h; hd0++) {

					DEO_DOE_PREAMBLE;
					// the following depends on the function
					// (d0+d1+d2+d3) even
					d0 = 2*hd0 + ((d1+d2+d3) & 0x1);

					d0m = d0 - 1;
					d0m = d0m + (((d0m >> 31) & 0x1) * nd0);
					d0p = d0 + 1;
					d0p *= (((d0p-nd0) >> 31) & 0x1);

#define SUB_RESULT_DEO(matdir, indexm)																	\
					aux = subResult(aux,conjmat_vec_mul_arg_wf( &u[matdir],indexm, \
																											in,indexm,&phases[matdir],&field_re[matdir], \
																											&field_im[matdir]));  
					SUB_RESULT_DEO(1,snum_acc(d0m,d1,d2,d3));
					SUB_RESULT_DEO(3,snum_acc(d0,d1m,d2,d3));
					SUB_RESULT_DEO(5,snum_acc(d0,d1,d2m,d3));
					SUB_RESULT_DEO(7,snum_acc(d0,d1,d2,d3m));
#undef SUB_RESULT_DEO

					//////////////////////////////////////////////////////////////
					idxh = snum_acc(d0,d1,d2,d3);

#define SUM_RESULT_DEO(matdir, indexp)																	\
					aux   = sumResult(aux, mat_vec_mul_arg_wf(&u[matdir],idxh,		\
																										in,indexp,&phases[matdir],&field_re[matdir], \
																										&field_im[matdir])); 
					SUM_RESULT_DEO(0,snum_acc(d0p,d1,d2,d3));
					SUM_RESULT_DEO(2,snum_acc(d0,d1p,d2,d3));
					SUM_RESULT_DEO(4,snum_acc(d0,d1,d2p,d3));
					SUM_RESULT_DEO(6,snum_acc(d0,d1,d2,d3p));
#undef SUM_RESULT_DEO

					///////////////////////////////////////////////////////////// 

					out->c0[idxh] = (aux.c0)*0.5;
					out->c1[idxh] = (aux.c1)*0.5;
					out->c2[idxh] = (aux.c2)*0.5;

				}}}}
}

// 'wf = with a field (for magnetic susceptibility)'
void acc_Doe_wf_unsafe( __restrict const su3_soa * const u,
												__restrict vec3_soa * const out,
												__restrict const vec3_soa * const in,
												const double_soa * phases,
												__restrict const double_soa* field_re,   // e.g. e^(i phi_1 ) - e^(i phi_2),
												__restrict const double_soa* field_im)   // or I * dphi/dbz  : notice the I
{
	int hd0, d1, d2, d3;
	#pragma acc kernels present(u) present(out) present(in) present(phases)	\
		present(field_re) present(field_im)
	#pragma acc loop independent gang(DEODOEGANG3)
	for(d3=D3_HALO; d3<D3_HALO+LOC_N3;d3++) {
		#pragma acc loop independent tile(DEODOETILE0,DEODOETILE1,DEODOETILE2)
		for(d2=0; d2<nd2; d2++) {
			for(d1=0; d1<nd1; d1++) {
				for(hd0=0; hd0 < nd0h; hd0++) {


					DEO_DOE_PREAMBLE;

					// the following depends on the function
					// (d0+d1+d2+d3) odd
					d0 = 2*hd0 + ((d1+d2+d3+1) & 0x1);

					d0m = d0 - 1;
					d0m = d0m + (((d0m >> 31) & 0x1) * nd0);
					d0p = d0 + 1;
					d0p *= (((d0p-nd0) >> 31) & 0x1);

#define SUB_RESULT_DOE(matdir,index)																		\
					aux = subResult(aux, conjmat_vec_mul_arg_wf( &u[matdir],index, \
																											 in,index,&phases[matdir],&field_re[matdir], \
																											 &field_im[matdir]));
					SUB_RESULT_DOE(0,snum_acc(d0m,d1,d2,d3));
					SUB_RESULT_DOE(2,snum_acc(d0,d1m,d2,d3));
					SUB_RESULT_DOE(4,snum_acc(d0,d1,d2m,d3));
					SUB_RESULT_DOE(6,snum_acc(d0,d1,d2,d3m));
#undef SUB_RESULT_DOE        

					///////////////////////////////////////////////////////////////////

					idxh = snum_acc(d0,d1,d2,d3);

#define SUM_RESULT_DOE(matdir,index)																		\
					aux   = sumResult(aux, mat_vec_mul_arg_wf(&u[matdir],idxh,		\
																										in,index,&phases[matdir],&field_re[matdir],	\
																										&field_im[matdir]));

					SUM_RESULT_DOE(1,snum_acc(d0p,d1,d2,d3));
					SUM_RESULT_DOE(3,snum_acc(d0,d1p,d2,d3));
					SUM_RESULT_DOE(5,snum_acc(d0,d1,d2p,d3));
					SUM_RESULT_DOE(7,snum_acc(d0,d1,d2,d3p));
#undef SUM_RESULT_DOE

					/////////////////////////////////////////////////////////////////

					out->c0[idxh] = aux.c0*0.5;
					out->c1[idxh] = aux.c1*0.5;
					out->c2[idxh] = aux.c2*0.5;

				}}}}
}

// 'wf = with a field (for magnetic susceptibility)'
void acc_Deo_wf( __restrict const su3_soa * const u, 
								 __restrict vec3_soa * const out, 
								 __restrict const vec3_soa * const in,
								 const double_soa * phases,
								 __restrict const double_soa* field_re,   // e.g. e^(i phi_1 ) - e^(i phi_2),
								 __restrict const double_soa* field_im)   // or i dphi/dbz 
{

#if NRANKS_D3 > 1	//#ifdef MULTIDEVICE
	acc_Deo_wf_unsafe(u, out, in, phases, field_re,field_im);
	communicate_fermion_borders(out); // contains host-device communications
#else 
	acc_Deo_wf_unsafe(u, out, in, phases, field_re,field_im);
#endif

}

// 'wf = with a field (for magnetic susceptibility)'
void acc_Doe_wf( __restrict const su3_soa * const u,
								 __restrict vec3_soa * const out,
								 __restrict const vec3_soa * const in,
								 const double_soa * phases,
								 __restrict const double_soa* field_re,   // e.g. e^(i phi_1 ) - e^(i phi_2),
								 __restrict const double_soa* field_im)   // or i dphi/dbz 
{

#if NRANKS_D3 > 1	//#ifdef MULTIDEVICE
	acc_Doe_wf_unsafe(u, out, in, phases,field_re,field_im);
	communicate_fermion_borders(out); // contains host-device communications
#else 
	acc_Doe_wf_unsafe(u, out, in, phases,field_re,field_im);
#endif

}


#endif
