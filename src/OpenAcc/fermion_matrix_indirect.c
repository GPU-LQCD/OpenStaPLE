#ifndef FERMION_MATRIX_C
#define FERMION_MATRIX_C

#include "../Include/common_defines.h"
#include "./geometry.h"
#include "./struct_c_def.h"
// #ifndef __GNUC__
#include "openacc.h"
// #endif
#include "./fermionic_utilities.h"
#include "./fermion_matrix.h"
#include "./matvecmul.h"

#ifdef MULTIDEVICE
#include "../Mpi/communications.h"
#endif 
#include "../Mpi/multidev.h"



void acc_Deo_unsafe_indaddr(__restrict const su3_soa * const u, 
														__restrict vec3_soa * const out, 
														__restrict const vec3_soa * const in,
														__restrict const double_soa * const backfield) {

	int hd0, d1, d2, d3;
	#pragma acc kernels present(u) present(out) present(in) present(backfield)
	#pragma acc loop independent gang(DEODOEGANG3)
	for(d3=D3_HALO; d3<D3_HALO+LOC_N3;d3++) {
		#pragma acc loop independent tile(DEODOETILE0,DEODOETILE1,DEODOETILE2)
		for(d2=0; d2<nd2; d2++) {
			for(d1=0; d1<nd1; d1++) {
				for(hd0=0; hd0 < nd0h; hd0++) {

					// the following depends on the function
					// (d0+d1+d2+d3) even
					int d0 = 2*hd0 + ((d1+d2+d3) & 0x1);
					int idxh = snum_acc(d0,d1,d2,d3);

#define SUB_RESULT_DEO(matdir, indexm)																	\
					aux = subResult(aux,conjmat_vec_mul_arg( &u[matdir],indexm,		\
																									 in,indexm,&backfield[matdir]));  
					SUB_RESULT_DEO(1,nnm_openacc[idxh][0][0]);
					SUB_RESULT_DEO(3,nnm_openacc[idxh][1][0]);
					SUB_RESULT_DEO(5,nnm_openacc[idxh][2][0]);
					SUB_RESULT_DEO(7,nnm_openacc[idxh][3][0]);
#undef SUB_RESULT_DEO

					//////////////////////////////////////////////////////////////

#define SUM_RESULT_DEO(matdir, indexp)																	\
					aux   = sumResult(aux, mat_vec_mul_arg(&u[matdir],idxh,				\
																								 in,indexp,&backfield[matdir])); 
					SUM_RESULT_DEO(0,nnp_openacc[idxh][0][0]);
					SUM_RESULT_DEO(2,nnp_openacc[idxh][1][0]);
					SUM_RESULT_DEO(4,nnp_openacc[idxh][2][0]);
					SUM_RESULT_DEO(6,nnp_openacc[idxh][3][0]);
#undef SUM_RESULT_DEO

					///////////////////////////////////////////////////////////// 

					out->c0[idxh] = (aux.c0)*0.5;
					out->c1[idxh] = (aux.c1)*0.5;
					out->c2[idxh] = (aux.c2)*0.5;

				}}}}
}

void acc_Doe_unsafe_indaddr(__restrict const su3_soa * const u,
														__restrict vec3_soa * const out,
														__restrict const vec3_soa * const in,
														__restrict const double_soa * const backfield) {

	int hd0, d1, d2, d3;
	#pragma acc kernels present(u) present(out) present(in) present(backfield)
	#pragma acc loop independent gang(DEODOEGANG3)
	for(d3=D3_HALO; d3<D3_HALO+LOC_N3;d3++) {
		#pragma acc loop independent tile(DEODOETILE0,DEODOETILE1,DEODOETILE2)
		for(d2=0; d2<nd2; d2++) {
			for(d1=0; d1<nd1; d1++) {
				for(hd0=0; hd0 < nd0h; hd0++) {


					// the following depends on the function
					// (d0+d1+d2+d3) odd
					int d0 = 2*hd0 + ((d1+d2+d3+1) & 0x1);
					idxh = snum_acc(d0,d1,d2,d3);

#define SUB_RESULT_DOE(matdir,index)																		\
					aux = subResult(aux, conjmat_vec_mul_arg( &u[matdir],index,		\
																										in,index,&backfield[matdir]));
					SUB_RESULT_DOE(0,nnm_openacc[idxh][0][1]);
					SUB_RESULT_DOE(2,nnm_openacc[idxh][1][1]);
					SUB_RESULT_DOE(4,nnm_openacc[idxh][2][1]);
					SUB_RESULT_DOE(6,nnm_openacc[idxh][3][1]);
#undef SUB_RESULT_DOE        

					///////////////////////////////////////////////////////////////////


#define SUM_RESULT_DOE(matdir,index)																		\
					aux   = sumResult(aux, mat_vec_mul_arg(&u[matdir],idxh,				\
																								 in,index,&backfield[matdir]));

					SUM_RESULT_DOE(1,nnm_openacc[idxh][0][1]);
					SUM_RESULT_DOE(3,nnm_openacc[idxh][1][1]);
					SUM_RESULT_DOE(5,nnm_openacc[idxh][2][1]);
					SUM_RESULT_DOE(7,nnm_openacc[idxh][3][1]);
#undef SUM_RESULT_DOE

					/////////////////////////////////////////////////////////////////

					out->c0[idxh] = aux.c0*0.5;
					out->c1[idxh] = aux.c1*0.5;
					out->c2[idxh] = aux.c2*0.5;

				}}}}
}

#endif
