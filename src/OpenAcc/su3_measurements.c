#ifndef SU3_MEASUREMENTS_C_
#define SU3_MEASUREMENTS_C_

#include "../Include/common_defines.h"
#include "../OpenAcc/struct_c_def.h"
#include "./plaquettes.h"
#include "./single_types.h"
#include "./su3_measurements.h"
#include "./su3_utilities.h"
#include "math.h"

// #ifndef __GNUC__
#include <openacc.h>
// #endif

#ifdef MULTIDEVICE
#include <mpi.h>
#include "../Mpi/multidev.h"
#endif


void check_unitarity_device(__restrict const su3_soa * const u, double * max_unitarity_deviation, double *avg_unitarity_deviation)
{
	// removing stag phases
	double r = 0;
	double rmax = 0;

#pragma acc kernels present(u)
#pragma acc loop reduction(+:r) reduction(max:rmax)
	for(int idx = 0; idx < sizeh ; idx++){
		for(int dir = 0; dir < 8 ; dir++){
			single_su3 m;
			single_su3_from_su3_soa(&u[dir],idx,&m);
			rebuild3row(&m);

			d_complex err = 1 - detSu3(&m);
			r += creal(err * conj(err));
			rmax = fmax(rmax,creal(err * conj(err)));

		}
	}
	printf("Warning: check max unitarity deviation calculation (%s:%d)\n",
				 __FILE__,__LINE__);

	// adding them again
	*avg_unitarity_deviation = r/(sizeh*8);
	*max_unitarity_deviation = rmax;


}
void check_unitarity_host( __restrict su3_soa * const u, double * max_unitarity_deviation, double *avg_unitarity_deviation){

	// removing stag phases


	double r = 0;
	double rmax = 0;

	for(int dir = 0; dir < 8 ; dir++){
		for(int idx = 0; idx < sizeh ; idx++){
			single_su3 m;
			single_su3_from_su3_soa(&u[dir],idx,&m);
			rebuild3row(&m);

			d_complex err = 1 - detSu3(&m);
			r += creal(err * conj(err));
			rmax = fmax(rmax,creal(err * conj(err)));

		}
	}

	// adding them again
	*avg_unitarity_deviation = r/(sizeh*8);
	*max_unitarity_deviation = rmax;


}

void check_unitarity_device_loc( __restrict su3_soa * const u, double * max_unitarity_deviation, double *avg_unitarity_deviation){

	// removing stag phases

	double r = 0;
	double rmax = 0;

#pragma acc kernels present(u)
#pragma acc loop reduction(+:r) reduction(max:rmax)
	for(int idx = (LNH_SIZEH-LOC_SIZEH)/2; idx < (LNH_SIZEH+LOC_SIZEH)/2 ; idx++){
		for(int dir = 0; dir < 8 ; dir++){
			single_su3 m;
			single_su3_from_su3_soa(&u[dir],idx,&m);
			rebuild3row(&m);

			d_complex err = 1 - detSu3(&m);
			r += creal(err * conj(err));
			rmax = fmax(rmax,creal(err * conj(err)));

		}
	}
	printf("Warning: check max unitarity deviation calculation (%s:%d)\n",
				 __FILE__,__LINE__);

	// adding them again
	*avg_unitarity_deviation = r/(sizeh*8);
	*max_unitarity_deviation = rmax;


}
void check_unitarity_host_loc( __restrict su3_soa * const u, double * max_unitarity_deviation, double *avg_unitarity_deviation){

	// removing stag phases

	double r = 0;
	double rmax = 0;

	for(int dir = 0; dir < 8 ; dir++){
		for(int idx = (LNH_SIZEH-LOC_SIZEH)/2; idx < (LNH_SIZEH+LOC_SIZEH)/2 ; idx++){
			single_su3 m;
			single_su3_from_su3_soa(&u[dir],idx,&m);
			rebuild3row(&m);

			d_complex err = 1 - detSu3(&m);
			r += creal(err * conj(err));
			rmax = fmax(rmax,creal(err * conj(err)));

		}
	}
	// adding them again

	*avg_unitarity_deviation = r/(sizeh*8);
	*max_unitarity_deviation = rmax;


}



double calc_momenta_action( const __restrict thmat_soa * const mom,
														double_soa * tr_local,
														const int mu){
	int t;

#pragma acc kernels present(mom) present(tr_local)
#pragma acc loop independent // gang(nt)
	for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
		tr_local[0].d[t] = half_tr_thmat_squared(&mom[mu],t);
	} // t


	double result=0.0;
	double total_result=0.0;
#pragma acc kernels present(tr_local)
#pragma acc loop reduction(+:result)
	for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
		result += tr_local[0].d[t];
	}


#if NRANKS_D3 > 1
	MPI_Allreduce((void*)&result,(void*)&total_result,1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);
#else
	total_result = result;
	#endif
	return total_result;

} // closes routine


void invert_momenta( __restrict thmat_soa * mom)
{
	int t,mu;

#pragma acc kernels present(mom)
	#pragma acc loop independent
	for(mu=0; mu < 8 ; mu++){
		#pragma acc loop independent
		for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
			mom[mu].c01[t] = -mom[mu].c01[t] ;
			mom[mu].c02[t] = -mom[mu].c02[t] ;
			mom[mu].c12[t] = -mom[mu].c12[t] ;
			mom[mu].rc00[t]= -mom[mu].rc00[t];
			mom[mu].rc11[t]= -mom[mu].rc11[t];
		} // t
	}
} // closes routine




double  calc_plaquette_soloopenacc(__restrict  su3_soa * const tconf_acc,
																	 __restrict su3_soa * const local_plaqs,
																	 dcomplex_soa * const tr_local_plaqs)
{

	double result=0.0;
	double total_result=0.0;

	int i_counter=0;
	for(int mu=0;mu<3;mu++){
		for(int nu=mu+1;nu<4;nu++){
			result  += calc_loc_plaquettes_nnptrick(tconf_acc,local_plaqs,tr_local_plaqs,mu,nu); // here all the plaquettes of a specific plane's choice are computed.
		}
	}

#if NRANKS_D3 > 1
	MPI_Allreduce((void*)&result,(void*)&total_result,
								1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);
#else
	total_result = result;
	#endif
	return total_result;
}

double calc_force_norm(const __restrict tamat_soa * tipdot)
{

	int t,mu;
	double result=0.0;
	double total_result=0.0;



#pragma acc kernels present(tipdot)
#pragma acc loop reduction(+:result)
	for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
#pragma acc loop reduction(+:result)
		for(mu=0; mu < 8; mu++) // NOTE: check if it works
			{
				result += half_tr_tamat_squared(&tipdot[mu],t);
			}
	}  
#if NRANKS_D3 > 1
	MPI_Allreduce((void*)&result,(void*)&total_result,
								1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);
#else
	total_result = result;
	#endif
	return sqrt(total_result);

};

double calc_diff_force_norm(const __restrict tamat_soa * tipdot,
														const __restrict tamat_soa * tipdot_old)
{

	int t,mu;
	double result=0.0;
	double total_result=0.0;


	for(mu=0; mu < 8; mu++) // NOTE: check if it works
		{
#pragma acc kernels present(tipdot) present(tipdot_old)
#pragma acc loop reduction(+:result)
			for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
				double A = tipdot[mu].ic00[t]-tipdot_old[mu].ic00[t];
				double B = tipdot[mu].ic11[t]-tipdot_old[mu].ic11[t];
				d_complex C = tipdot[mu].c01[t]-tipdot_old[mu].c01[t];
				d_complex D = tipdot[mu].c02[t]-tipdot_old[mu].c02[t];
				d_complex E = tipdot[mu].c12[t]-tipdot_old[mu].c12[t];
				result +=  A*A + B*B + A*B
					+ creal(C)*creal(C) + cimag(C)*cimag(C)
					+ creal(D)*creal(D) + cimag(D)*cimag(D)
					+ creal(E)*creal(E) + cimag(E)*cimag(E);

			}
		}

#if NRANKS_D3 > 1
	MPI_Allreduce((void*)&result,(void*)&total_result,
								1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);
#else
	total_result = result;
	#endif
	return sqrt(total_result);

};


double calc_sum_momenta_norm(const __restrict thmat_soa * tmomenta,
														 const __restrict thmat_soa * tmomenta_old)
{

	int t,mu;
	double result=0.0;
	double total_result=0.0;


	for(mu=0; mu < 8; mu++) // NOTE: check if it works
		{
#pragma acc kernels present(tmomenta) present(tmomenta_old)
#pragma acc loop reduction(+:result)
			for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
				double A = tmomenta[mu].rc00[t]+tmomenta_old[mu].rc00[t];
				double B = tmomenta[mu].rc11[t]+tmomenta_old[mu].rc11[t];
				d_complex C = tmomenta[mu].c01[t]+tmomenta_old[mu].c01[t];
				d_complex D = tmomenta[mu].c02[t]+tmomenta_old[mu].c02[t];
				d_complex E = tmomenta[mu].c12[t]+tmomenta_old[mu].c12[t];
				result +=  A*A + B*B + A*B
					+ creal(C)*creal(C) + cimag(C)*cimag(C)
					+ creal(D)*creal(D) + cimag(D)*cimag(D)
					+ creal(E)*creal(E) + cimag(E)*cimag(E);
			}
		}

#if NRANKS_D3 > 1
	MPI_Allreduce((void*)&result,(void*)&total_result,
								1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);
#else
	total_result = result;
	#endif
	return sqrt(total_result);

};


double calc_diff_su3_soa_norm(const __restrict su3_soa * tconf,
															const __restrict su3_soa * tconf_old)
{

	int t,mu;
	double result=0.0;
	double total_result=0.0;

	for(mu=0; mu < 8; mu++) // NOTE: check if it works
		{
#pragma acc kernels present(tconf) present(tconf_old)
#pragma acc loop reduction(+:result)
			for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
				d_complex A = tconf[mu].r0.c0[t]-tconf_old[mu].r0.c0[t];
				d_complex B = tconf[mu].r0.c1[t]-tconf_old[mu].r0.c1[t];
				d_complex C = tconf[mu].r0.c2[t]-tconf_old[mu].r0.c2[t];
				d_complex D = tconf[mu].r1.c0[t]-tconf_old[mu].r1.c0[t];
				d_complex E = tconf[mu].r1.c1[t]-tconf_old[mu].r1.c1[t];
				d_complex F = tconf[mu].r1.c2[t]-tconf_old[mu].r1.c2[t];
				result += creal(A)*creal(A) + cimag(A)*cimag(A)
					+ creal(B)*creal(B) + cimag(B)*cimag(B)
					+ creal(C)*creal(C) + cimag(C)*cimag(C)
					+ creal(D)*creal(D) + cimag(D)*cimag(D)
					+ creal(E)*creal(E) + cimag(E)*cimag(E)
					+ creal(F)*creal(F) + cimag(F)*cimag(F);

			}
		}

#if NRANKS_D3 > 1
	MPI_Allreduce((void*)&result,(void*)&total_result,
								1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);
#else
	total_result = result;
	#endif
	return sqrt(total_result);

};



void copy_ipdot_into_old(const __restrict tamat_soa * tipdot,
												 __restrict tamat_soa * tipdot_old)
{

	// NOTE: possibly some problems here
	int t,mu;

#pragma acc kernels present(tipdot) present(tipdot_old)
	#pragma acc loop independent
	for(mu=0; mu < 8; mu++) // NOTE: check if it works
		{
			#pragma acc loop independent
			for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
				tipdot_old[mu].ic00[t] = tipdot[mu].ic00[t];
				tipdot_old[mu].ic11[t] = tipdot[mu].ic11[t];
				tipdot_old[mu].c01[t]  = tipdot[mu].c01[t] ;
				tipdot_old[mu].c02[t]  = tipdot[mu].c02[t] ;
				tipdot_old[mu].c12[t]  = tipdot[mu].c12[t] ;

			}
		}

}

void copy_momenta_into_old(
													 const __restrict thmat_soa * tmomenta,
													 __restrict thmat_soa * tmomenta_old)
{

	// NOTE: possibly some problems here
	int t,mu;

#pragma acc kernels present(tmomenta) present(tmomenta_old)
	#pragma acc loop independent
	for(mu=0; mu < 8; mu++) // NOTE: check if it works
		{
			#pragma acc loop independent
			for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
				tmomenta_old[mu].rc00[t] = tmomenta[mu].rc00[t];
				tmomenta_old[mu].rc11[t] = tmomenta[mu].rc11[t];
				tmomenta_old[mu].c01[t]  = tmomenta[mu].c01[t] ;
				tmomenta_old[mu].c02[t]  = tmomenta[mu].c02[t] ;
				tmomenta_old[mu].c12[t]  = tmomenta[mu].c12[t] ;
			}
		}
}


#endif




















































































































































































































































































































































































































