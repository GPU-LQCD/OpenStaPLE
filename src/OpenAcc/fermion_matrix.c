#ifndef FERMION_MATRIX_C
#define FERMION_MATRIX_C

#include "../Include/common_defines.h"
#include "./geometry.h"
#include "./struct_c_def.h"
#ifndef __GNUC__
#include "openacc.h"
#endif

#ifdef __GNUC__
#include "sys/time.h"
#endif
#include "./fermionic_utilities.h"
#include "./fermion_matrix.h"
#include "./matvecmul.h"
#include "../tests_and_benchmarks/test_and_benchmarks.h"

//TODO: check if this can be changed to NRANKS_D3 > 1
#ifdef MULTIDEVICE
#include "../Mpi/communications.h"
#endif 
#include "../Mpi/multidev.h"


#define DEO_DOE_PREAMBLE																\
	int d0, d0m, d1m, d2m, d3m, d0p, d1p, d2p, d3p, idxh;	\
	vec3 aux=(vec3) {0,0,0};															\
	d1m = d1 - 1;																					\
	d1m = d1m + (((d1m >> 31) & 0x1) * nd1);							\
	d2m = d2 - 1;																					\
	d2m = d2m + (((d2m >> 31) & 0x1) * nd2);							\
	d3m = d3 - 1;																					\
	d3m = d3m + (((d3m >> 31) & 0x1) * nd3);							\
	d1p = d1 + 1;																					\
	d1p *= (((d1p-nd1) >> 31) & 0x1);											\
	d2p = d2 + 1;																					\
	d2p *= (((d2p-nd2) >> 31) & 0x1);											\
	d3p = d3 + 1;																					\
	d3p *= (((d3p-nd3) >> 31) & 0x1);											\



#define CLOSE4CYCLES }}}}


void acc_Deo_unsafe(__restrict const su3_soa * const u, 
										__restrict vec3_soa * const out, 
										__restrict const vec3_soa * const in,
										__restrict const double_soa * const backfield)
{
	int hd0, d1, d2, d3;
	#pragma acc kernels present(u) present(out) present(in) present(backfield)
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
					aux = subResult(aux,conjmat_vec_mul_arg( &u[matdir],indexm,		\
																									 in,indexm,&backfield[matdir]));  
					SUB_RESULT_DEO(1,snum_acc(d0m,d1,d2,d3));
					SUB_RESULT_DEO(3,snum_acc(d0,d1m,d2,d3));
					SUB_RESULT_DEO(5,snum_acc(d0,d1,d2m,d3));
					SUB_RESULT_DEO(7,snum_acc(d0,d1,d2,d3m));
#undef SUB_RESULT_DEO

					//////////////////////////////////////////////////////////////
					idxh = snum_acc(d0,d1,d2,d3);

#define SUM_RESULT_DEO(matdir, indexp)																	\
					aux   = sumResult(aux, mat_vec_mul_arg(&u[matdir],idxh,				\
																								 in,indexp,&backfield[matdir])); 
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

void acc_Doe_unsafe(__restrict const su3_soa * const u,
										__restrict vec3_soa * const out,
										__restrict const vec3_soa * const in,
										__restrict const double_soa * const backfield)
{
	int hd0, d1, d2, d3;
	#pragma acc kernels present(u) present(out) present(in) present(backfield)
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
					aux = subResult(aux, conjmat_vec_mul_arg( &u[matdir],index,		\
																										in,index,&backfield[matdir]));
					SUB_RESULT_DOE(0,snum_acc(d0m,d1,d2,d3));
					SUB_RESULT_DOE(2,snum_acc(d0,d1m,d2,d3));
					SUB_RESULT_DOE(4,snum_acc(d0,d1,d2m,d3));
					SUB_RESULT_DOE(6,snum_acc(d0,d1,d2,d3m));
#undef SUB_RESULT_DOE        

					///////////////////////////////////////////////////////////////////

					idxh = snum_acc(d0,d1,d2,d3);

#define SUM_RESULT_DOE(matdir,index)																		\
					aux   = sumResult(aux, mat_vec_mul_arg(&u[matdir],idxh,				\
																								 in,index,&backfield[matdir]));

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

void acc_Deo(__restrict const su3_soa * const u, 
						 __restrict vec3_soa * const out, 
						 __restrict const vec3_soa * const in,
						 __restrict const double_soa * const backfield)
{
#if NRANKS_D3 > 1
	if(devinfo.async_comm_fermion){

#if defined(USE_MPI_CUDA_AWARE) || defined(__GNUC__)
		MPI_Request send_border_requests[6], recv_border_requests[6];

		//        acc_Deo_d3c(u, out, in, backfield,HALO_WIDTH,1);
		//        acc_Deo_d3c(u, out, in, backfield,nd3-HALO_WIDTH-1,1);
		acc_Deo_d3p(u, out, in, backfield); // async(2)
		acc_Deo_d3m(u, out, in, backfield); // async(3)

		acc_Deo_bulk(u, out, in, backfield); // async(1)
		#pragma acc wait(2)
		#pragma acc wait(3)
		communicate_fermion_borders_async(out,send_border_requests,
																			recv_border_requests);
		//        acc_Deo_bulk(u, out, in, backfield); //async(1)
		MPI_Waitall(6,recv_border_requests,MPI_STATUSES_IGNORE);
		MPI_Waitall(6,send_border_requests,MPI_STATUSES_IGNORE);
		#pragma acc wait(1)
#else 
		printf("ERROR:Async transfer involving accelerators require USE_MPI_CUDA_AWARE!\n");
		printf("EXITING NOW\n");
		MPI_Finalize();
		exit(1);
#endif
	}else{ 
		struct timeval t_transfer_start;
		struct timeval t_transfer_end;

		acc_Deo_unsafe(u, out, in, backfield);

		if (0 == devinfo.myrank ) gettimeofday(&t_transfer_start,NULL);
		communicate_fermion_borders(out); // contains host-device communications
		if (0 == devinfo.myrank ){
			gettimeofday(&t_transfer_end,NULL);
			dirac_times.totTransferTime += (double) 
				(t_transfer_end.tv_sec -t_transfer_start.tv_sec)+
				(double)(t_transfer_end.tv_usec -t_transfer_start.tv_usec)/1.0e6;
			++dirac_times.count;

		}
	}
	//    MPI_Barrier(devinfo.mpi_comm);
#else 
	acc_Deo_unsafe(u, out, in, backfield);
#endif

}

void acc_Doe(__restrict const su3_soa * const u,
						 __restrict vec3_soa * const out,
						 __restrict const vec3_soa * const in,
						 __restrict const double_soa * const backfield)
{
#if NRANKS_D3 > 1
	if(devinfo.async_comm_fermion){

#if defined(USE_MPI_CUDA_AWARE) || defined(__GNUC__)
		MPI_Request send_border_requests[6], recv_border_requests[6];

		//        acc_Doe_d3c(u, out, in, backfield,HALO_WIDTH,1);
		//        acc_Doe_d3c(u, out, in, backfield,nd3-HALO_WIDTH-1,1);
		acc_Doe_d3p(u, out, in, backfield); // async(2)
		acc_Doe_d3m(u, out, in, backfield); // async(3)

		acc_Doe_bulk(u, out, in, backfield); // async(1)
		#pragma acc wait(2)
		#pragma acc wait(3)
		communicate_fermion_borders_async(out,send_border_requests,
																			recv_border_requests);
		//        acc_Doe_bulk(u, out, in, backfield); //async(1)
		MPI_Waitall(6,recv_border_requests,MPI_STATUSES_IGNORE);
		MPI_Waitall(6,send_border_requests,MPI_STATUSES_IGNORE);
		#pragma acc wait(1)
#else 
		printf("ERROR:Async transfer involving accelerators require USE_MPI_CUDA_AWARE!\n");
		printf("EXITING NOW\n");
		MPI_Finalize();
		exit(1);
#endif

	}else{ 
		struct timeval t_transfer_start;
		struct timeval t_transfer_end;

		acc_Doe_unsafe(u, out, in, backfield);

		if (0 == devinfo.myrank ) gettimeofday(&t_transfer_start,NULL);
		communicate_fermion_borders(out); // contains host-device communications
		if (0 == devinfo.myrank ){
			gettimeofday(&t_transfer_end,NULL);
			dirac_times.totTransferTime += (double) 
				(t_transfer_end.tv_sec -t_transfer_start.tv_sec)+
				(double)(t_transfer_end.tv_usec -t_transfer_start.tv_usec)/1.0e6;
			++dirac_times.count;

		}
	}
	//    MPI_Barrier(devinfo.mpi_comm);
#else 
	acc_Doe_unsafe(u, out, in, backfield);
#endif

}


void acc_Deo_bulk(__restrict const su3_soa * const u, 
									__restrict vec3_soa * const out, 
									__restrict const vec3_soa * const in,
									__restrict const double_soa * const backfield)
{
	int hd0, d1, d2, d3;
	#pragma acc kernels present(u) present(out) present(in) present(backfield) async(1)
	#pragma acc loop independent gang(DEODOEGANG3)
	for(d3=D3_HALO+1; d3<D3_HALO+1+LOC_N3-2;d3++) {
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
					aux = subResult(aux,conjmat_vec_mul_arg( &u[matdir],indexm,		\
																									 in,indexm,&backfield[matdir]));  
					SUB_RESULT_DEO(1,snum_acc(d0m,d1,d2,d3));
					SUB_RESULT_DEO(3,snum_acc(d0,d1m,d2,d3));
					SUB_RESULT_DEO(5,snum_acc(d0,d1,d2m,d3));
					SUB_RESULT_DEO(7,snum_acc(d0,d1,d2,d3m));
#undef SUB_RESULT_DEO

					//////////////////////////////////////////////////////////////
					idxh = snum_acc(d0,d1,d2,d3);

#define SUM_RESULT_DEO(matdir, indexp)																	\
					aux   = sumResult(aux, mat_vec_mul_arg(&u[matdir],idxh,				\
																								 in,indexp,&backfield[matdir])); 
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

void acc_Doe_bulk(__restrict const su3_soa * const u,
									__restrict vec3_soa * const out,
									__restrict const vec3_soa * const in,
									__restrict const double_soa * const backfield)
{
	int hd0, d1, d2, d3;
	#pragma acc kernels present(u) present(out) present(in) present(backfield) async(1)
	#pragma acc loop independent gang(DEODOEGANG3)
	for(d3=D3_HALO+1; d3<D3_HALO+1+LOC_N3-2;d3++) {
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
					aux = subResult(aux, conjmat_vec_mul_arg( &u[matdir],index,		\
																										in,index,&backfield[matdir]));
					SUB_RESULT_DOE(0,snum_acc(d0m,d1,d2,d3));
					SUB_RESULT_DOE(2,snum_acc(d0,d1m,d2,d3));
					SUB_RESULT_DOE(4,snum_acc(d0,d1,d2m,d3));
					SUB_RESULT_DOE(6,snum_acc(d0,d1,d2,d3m));
#undef SUB_RESULT_DOE        

					///////////////////////////////////////////////////////////////////

					idxh = snum_acc(d0,d1,d2,d3);

#define SUM_RESULT_DOE(matdir,index)																		\
					aux   = sumResult(aux, mat_vec_mul_arg(&u[matdir],idxh,				\
																								 in,index,&backfield[matdir]));

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



void acc_Deo_d3c(__restrict const su3_soa * const u, 
								 __restrict vec3_soa * const out, 
								 __restrict const vec3_soa * const in,
								 __restrict const double_soa * const backfield, int off3, int thick3)
{
	int hd0, d1, d2, d3;
	for(d3=off3; d3<off3+thick3;d3++) {
		#pragma acc kernels present(u) present(out) present(in) present(backfield) async(2)
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
					aux = subResult(aux,conjmat_vec_mul_arg( &u[matdir],indexm,		\
																									 in,indexm,&backfield[matdir]));  
					SUB_RESULT_DEO(1,snum_acc(d0m,d1,d2,d3));
					SUB_RESULT_DEO(3,snum_acc(d0,d1m,d2,d3));
					SUB_RESULT_DEO(5,snum_acc(d0,d1,d2m,d3));
					SUB_RESULT_DEO(7,snum_acc(d0,d1,d2,d3m));
#undef SUB_RESULT_DEO

					//////////////////////////////////////////////////////////////
					idxh = snum_acc(d0,d1,d2,d3);

#define SUM_RESULT_DEO(matdir, indexp)																	\
					aux   = sumResult(aux, mat_vec_mul_arg(&u[matdir],idxh,				\
																								 in,indexp,&backfield[matdir])); 
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

void acc_Doe_d3c(__restrict const su3_soa * const u,
								 __restrict vec3_soa * const out,
								 __restrict const vec3_soa * const in,
								 __restrict const double_soa * const backfield, int off3, int thick3)
{
	int hd0, d1, d2, d3;
	for(d3=off3; d3<off3+thick3;d3++) {
		#pragma acc kernels present(u) present(out) present(in) present(backfield) async(2)
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
					aux = subResult(aux, conjmat_vec_mul_arg( &u[matdir],index,		\
																										in,index,&backfield[matdir]));
					SUB_RESULT_DOE(0,snum_acc(d0m,d1,d2,d3));
					SUB_RESULT_DOE(2,snum_acc(d0,d1m,d2,d3));
					SUB_RESULT_DOE(4,snum_acc(d0,d1,d2m,d3));
					SUB_RESULT_DOE(6,snum_acc(d0,d1,d2,d3m));
#undef SUB_RESULT_DOE        

					///////////////////////////////////////////////////////////////////

					idxh = snum_acc(d0,d1,d2,d3);

#define SUM_RESULT_DOE(matdir,index)																		\
					aux   = sumResult(aux, mat_vec_mul_arg(&u[matdir],idxh,				\
																								 in,index,&backfield[matdir]));

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

void acc_Deo_d3p(__restrict const su3_soa * const u, 
								 __restrict vec3_soa * const out, 
								 __restrict const vec3_soa * const in,
								 __restrict const double_soa * const backfield)
{
	int hd0, d1, d2;
	const int d3 = nd3-D3_HALO-1;
    
	#pragma acc kernels present(u) present(out) present(in) present(backfield) async(2)
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
				aux = subResult(aux,conjmat_vec_mul_arg( &u[matdir],indexm,			\
																								 in,indexm,&backfield[matdir]));  
				SUB_RESULT_DEO(1,snum_acc(d0m,d1,d2,d3));
				SUB_RESULT_DEO(3,snum_acc(d0,d1m,d2,d3));
				SUB_RESULT_DEO(5,snum_acc(d0,d1,d2m,d3));
				SUB_RESULT_DEO(7,snum_acc(d0,d1,d2,d3m));
#undef SUB_RESULT_DEO

				//////////////////////////////////////////////////////////////
				idxh = snum_acc(d0,d1,d2,d3);

#define SUM_RESULT_DEO(matdir, indexp)																	\
				aux   = sumResult(aux, mat_vec_mul_arg(&u[matdir],idxh,					\
																							 in,indexp,&backfield[matdir])); 
				SUM_RESULT_DEO(0,snum_acc(d0p,d1,d2,d3));
				SUM_RESULT_DEO(2,snum_acc(d0,d1p,d2,d3));
				SUM_RESULT_DEO(4,snum_acc(d0,d1,d2p,d3));
				SUM_RESULT_DEO(6,snum_acc(d0,d1,d2,d3p));
#undef SUM_RESULT_DEO

				///////////////////////////////////////////////////////////// 

				out->c0[idxh] = (aux.c0)*0.5;
				out->c1[idxh] = (aux.c1)*0.5;
				out->c2[idxh] = (aux.c2)*0.5;

			}}}
}

void acc_Doe_d3p(__restrict const su3_soa * const u,
								 __restrict vec3_soa * const out,
								 __restrict const vec3_soa * const in,
								 __restrict const double_soa * const backfield)
{
	int hd0, d1, d2;
	const int d3 = nd3-D3_HALO-1;

	#pragma acc kernels present(u) present(out) present(in) present(backfield) async(2)
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
				aux = subResult(aux, conjmat_vec_mul_arg( &u[matdir],index,			\
																									in,index,&backfield[matdir]));
				SUB_RESULT_DOE(0,snum_acc(d0m,d1,d2,d3));
				SUB_RESULT_DOE(2,snum_acc(d0,d1m,d2,d3));
				SUB_RESULT_DOE(4,snum_acc(d0,d1,d2m,d3));
				SUB_RESULT_DOE(6,snum_acc(d0,d1,d2,d3m));
#undef SUB_RESULT_DOE        

				///////////////////////////////////////////////////////////////////

				idxh = snum_acc(d0,d1,d2,d3);

#define SUM_RESULT_DOE(matdir,index)																		\
				aux   = sumResult(aux, mat_vec_mul_arg(&u[matdir],idxh,					\
																							 in,index,&backfield[matdir]));

				SUM_RESULT_DOE(1,snum_acc(d0p,d1,d2,d3));
				SUM_RESULT_DOE(3,snum_acc(d0,d1p,d2,d3));
				SUM_RESULT_DOE(5,snum_acc(d0,d1,d2p,d3));
				SUM_RESULT_DOE(7,snum_acc(d0,d1,d2,d3p));
#undef SUM_RESULT_DOE

				/////////////////////////////////////////////////////////////////

				out->c0[idxh] = aux.c0*0.5;
				out->c1[idxh] = aux.c1*0.5;
				out->c2[idxh] = aux.c2*0.5;

			}}}
}

void acc_Deo_d3m(__restrict const su3_soa * const u, 
								 __restrict vec3_soa * const out, 
								 __restrict const vec3_soa * const in,
								 __restrict const double_soa * const backfield)
{
	int hd0, d1, d2;
	const int  d3 = D3_HALO;
	#pragma acc kernels present(u) present(out) present(in) present(backfield) async(3)
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
				aux = subResult(aux,conjmat_vec_mul_arg( &u[matdir],indexm,			\
																								 in,indexm,&backfield[matdir]));  
				SUB_RESULT_DEO(1,snum_acc(d0m,d1,d2,d3));
				SUB_RESULT_DEO(3,snum_acc(d0,d1m,d2,d3));
				SUB_RESULT_DEO(5,snum_acc(d0,d1,d2m,d3));
				SUB_RESULT_DEO(7,snum_acc(d0,d1,d2,d3m));
#undef SUB_RESULT_DEO

				//////////////////////////////////////////////////////////////
				idxh = snum_acc(d0,d1,d2,d3);

#define SUM_RESULT_DEO(matdir, indexp)																	\
				aux   = sumResult(aux, mat_vec_mul_arg(&u[matdir],idxh,					\
																							 in,indexp,&backfield[matdir])); 
				SUM_RESULT_DEO(0,snum_acc(d0p,d1,d2,d3));
				SUM_RESULT_DEO(2,snum_acc(d0,d1p,d2,d3));
				SUM_RESULT_DEO(4,snum_acc(d0,d1,d2p,d3));
				SUM_RESULT_DEO(6,snum_acc(d0,d1,d2,d3p));
#undef SUM_RESULT_DEO

				///////////////////////////////////////////////////////////// 

				out->c0[idxh] = (aux.c0)*0.5;
				out->c1[idxh] = (aux.c1)*0.5;
				out->c2[idxh] = (aux.c2)*0.5;

			}}}
}

void acc_Doe_d3m(__restrict const su3_soa * const u,
								 __restrict vec3_soa * const out,
								 __restrict const vec3_soa * const in,
								 __restrict const double_soa * const backfield)
{
	int hd0, d1, d2;
	const int  d3 = D3_HALO;
	#pragma acc kernels present(u) present(out) present(in) present(backfield) async(3)
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
				aux = subResult(aux, conjmat_vec_mul_arg( &u[matdir],index,			\
																									in,index,&backfield[matdir]));
				SUB_RESULT_DOE(0,snum_acc(d0m,d1,d2,d3));
				SUB_RESULT_DOE(2,snum_acc(d0,d1m,d2,d3));
				SUB_RESULT_DOE(4,snum_acc(d0,d1,d2m,d3));
				SUB_RESULT_DOE(6,snum_acc(d0,d1,d2,d3m));
#undef SUB_RESULT_DOE        

				///////////////////////////////////////////////////////////////////

				idxh = snum_acc(d0,d1,d2,d3);

#define SUM_RESULT_DOE(matdir,index)																		\
				aux   = sumResult(aux, mat_vec_mul_arg(&u[matdir],idxh,					\
																							 in,index,&backfield[matdir]));

				SUM_RESULT_DOE(1,snum_acc(d0p,d1,d2,d3));
				SUM_RESULT_DOE(3,snum_acc(d0,d1p,d2,d3));
				SUM_RESULT_DOE(5,snum_acc(d0,d1,d2p,d3));
				SUM_RESULT_DOE(7,snum_acc(d0,d1,d2,d3p));
#undef SUM_RESULT_DOE

				/////////////////////////////////////////////////////////////////

				out->c0[idxh] = aux.c0*0.5;
				out->c1[idxh] = aux.c1*0.5;
				out->c2[idxh] = aux.c2*0.5;

			}}}
}




void fermion_matrix_multiplication(__restrict const su3_soa * const u, 
																	 __restrict vec3_soa * const out,  
																	 __restrict const vec3_soa * const in, 
																	 __restrict vec3_soa * const temp1, 
																	 ferm_param *pars) {

	acc_Doe(u, temp1, in, pars->phases);
	acc_Deo(u, out, temp1, pars->phases);
	combine_in1xferm_mass_minus_in2(in,pars->ferm_mass*pars->ferm_mass,out);// function added in OpenAcc/fermionic_utilities.c

}

void fermion_matrix_multiplication_shifted(__restrict const su3_soa * const u, 
																					 __restrict vec3_soa * const out, 
																					 __restrict const vec3_soa * const in, 
																					 __restrict vec3_soa * const temp1, 
																					 ferm_param *pars, 
																					 double shift) {

	acc_Doe(u, temp1, in, pars->phases);
	acc_Deo(u, out, temp1, pars->phases);
	combine_in1xferm_mass_minus_in2(in,pars->ferm_mass*pars->ferm_mass+shift,out);// function added in OpenAcc/fermionic_utilities.c

}

#endif
