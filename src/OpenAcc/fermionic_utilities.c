// MEMO: functions defined in this file:
/*
	- scalar product utilities
	+ scal_prod_loc_2double(in1,in2,idx,d_re,d_im)         --> compute the scalar product locally
	+ scal_prod_global(in1,in2,result)                     --> compute the global scalar product
	- inverter (and other stuff) utilities
	+ combine_in1xm2_minus_in2_minus_in3(in1,in2,in3,out)  --> compute out = in1*m2-in2-in3
	+ combine_in1xm2_minus_in2(in1,in2,out)                --> compute out = in1*m2-in2
	+ combine_in1_minus_in2(in1,in2,out)                   --> compute out = in1-in2
	+ assign_in_to_out(in,out)                             --> set     out = in
	+ combine_in1xfactor_plus_in2(in1,fattore,in2,out)     --> compute out = in1*fattore+in2
*/


#ifndef FERMIONIC_UTILITIES_C_
#define FERMIONIC_UTILITIES_C_

#ifdef __GNUC__
#define __restrict
#define _POSIX_C_SOURCE 200809L // not to have warning on posix memalign
#endif

#include "./geometry.h"
#include "./fermionic_utilities.h"

#include <stdlib.h>

#if NRANKS_D3 > 1
#include <mpi.h>
#include "../Mpi/multidev.h"

d_complex scal_prod_loc_1Dcut(const __restrict vec3_soa * const in_vect1,
															const __restrict vec3_soa * const in_vect2)
{
	int t;
	double resR = 0.0;
	double resI = 0.0;

	#pragma acc kernels present(in_vect1) present(in_vect2)
	#pragma acc loop reduction(+:resR) reduction(+:resI) 
	for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
		d_complex color_sum  =  conj(in_vect1->c0[t]) *  in_vect2->c0[t] ;
		color_sum +=  conj(in_vect1->c1[t]) *  in_vect2->c1[t] ;
		color_sum +=  conj(in_vect1->c2[t]) *  in_vect2->c2[t] ;

		resR+=creal(color_sum);
		resI+=cimag(color_sum);

	}
	return resR+resI*I;

}
double real_scal_prod_loc_1Dcut(const __restrict vec3_soa * in_vect1,
																const __restrict vec3_soa * in_vect2)
{
	int t;
	double res_R_p;
	double resR = 0.0;

	#pragma acc kernels present(in_vect1) present(in_vect2)
	#pragma acc loop reduction(+:resR) // debug
	for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
		res_R_p=scal_prod_loc_1double(in_vect1,in_vect2,t);
		resR+=res_R_p;
	}
	return resR;
}
double l2norm2_loc_1Dcut(const __restrict vec3_soa * const in_vect1)
{
	int t;
	double res_R_p;
	double resR = 0.0;

	#pragma acc kernels present(in_vect1)
	#pragma acc loop reduction(+:resR) // debug
	for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
		res_R_p=l2norm2_loc(in_vect1,t);
		resR+=res_R_p;
	}
	return resR;
}


// mpi requiring functions
d_complex scal_prod_global(const __restrict vec3_soa * const in_vect1,
													 const __restrict vec3_soa * const in_vect2)
{
	// returns the total result using mpi

	d_complex local_res;

	local_res = scal_prod_loc_1Dcut(in_vect1,in_vect2); 
	double lrr = creal(local_res);
	double lri = cimag(local_res);
	double trr,tri;

	MPI_Allreduce((void*)&lrr,(void*)&trr,1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);
	MPI_Allreduce((void*)&lri,(void*)&tri,1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);

	return trr + tri * I;
}
double real_scal_prod_global(const __restrict vec3_soa * const in_vect1,
														 const __restrict vec3_soa * const in_vect2)
{
	// returns the total result using mpi
	double total_res,local_res;
	// local_res = real_scal_prod_loc(in_vect1,in_vect2); 
	local_res = real_scal_prod_loc_1Dcut(in_vect1,in_vect2); 
	MPI_Allreduce((void*)&local_res,(void*)&total_res,1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);
	return total_res;
}
double l2norm2_global(const __restrict vec3_soa * const in_vect1)
{

	double total_res,local_res;
	// local_res = l2norm2_loc(in_vect1); 
	local_res = l2norm2_loc_1Dcut(in_vect1); 
	MPI_Allreduce((void*)&local_res,(void*)&total_res,1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);
	return total_res;
}


#else // no NRANKS_D3 > 1


d_complex scal_prod_global(const __restrict vec3_soa * in_vect1,
													 const __restrict vec3_soa * in_vect2)
{
	int t;

	double resR = 0.0;
	double resI = 0.0;

	#pragma acc kernels present(in_vect1) present(in_vect2)
	#pragma acc loop reduction(+:resR) reduction(+:resI)
	for(t=0; t<sizeh; t++) {
		d_complex color_sum  =  conj(in_vect1->c0[t]) *  in_vect2->c0[t] ;
		color_sum +=  conj(in_vect1->c1[t]) *  in_vect2->c1[t] ;
		color_sum +=  conj(in_vect1->c2[t]) *  in_vect2->c2[t] ;

		resR+=creal(color_sum);
		resI+=cimag(color_sum);

	}
	return resR+resI*I;
}
double real_scal_prod_global(__restrict const  vec3_soa * in_vect1,
														 __restrict const vec3_soa * in_vect2)
{
	int t;
	double res_R_p;
	double resR = 0.0;

	#pragma acc kernels present(in_vect1) present(in_vect2)
	#pragma acc loop reduction(+:resR)
	for(t=0; t<sizeh; t++) {
		res_R_p=scal_prod_loc_1double(in_vect1,in_vect2,t);
		resR+=res_R_p;
	}
	return resR;
}
double l2norm2_global( __restrict  const vec3_soa * in_vect1 )
{
	int t;
	double res_R_p;
	double resR = 0.0;

	#pragma acc kernels present(in_vect1)
	#pragma acc loop reduction(+:resR)
	for(t=0; t<sizeh; t++) {
		res_R_p=l2norm2_loc(in_vect1,t);
		resR+=res_R_p;
	}
	return resR;
}

#endif // NRANKS_D3 > 1

// starting from here, only "embarassingly parallel" functions, not requiring any reduction.
void combine_in1xfactor_plus_in2(__restrict const vec3_soa * in_vect1,
																 const double factor,
																 __restrict const vec3_soa * in_vect2,
																 __restrict vec3_soa * const out)
{

	#pragma acc kernels present(in_vect1) present(in_vect2) present(out)
	#pragma acc loop independent 
	for(int ih=LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih < LNH_SIZEH-LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih++) {
		out->c0[ih]=(in_vect1->c0[ih]*factor)+in_vect2->c0[ih];
		out->c1[ih]=(in_vect1->c1[ih]*factor)+in_vect2->c1[ih];
		out->c2[ih]=(in_vect1->c2[ih]*factor)+in_vect2->c2[ih];
	}
}


void multiply_fermion_x_doublefactor(__restrict vec3_soa * const in1,
																		 const double factor)
{

	#pragma acc kernels present(in1)
	#pragma acc loop independent
	for(int ih=LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih < LNH_SIZEH-LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih++) {
		in1->c0[ih] = factor * (in1->c0[ih]);
		in1->c1[ih] = factor * (in1->c1[ih]);
		in1->c2[ih] = factor * (in1->c2[ih]);
	}
}

void combine_add_factor_x_in2_to_in1(__restrict vec3_soa * const in1,
																		 __restrict const vec3_soa * const in2, 
																		 double factor)
{
	#pragma acc kernels present(in1) present(in2)
	#pragma acc loop independent
	for(int ih=LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih < LNH_SIZEH-LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih++) {
		in1->c0[ih] += (factor) *  (in2->c0[ih]);
		in1->c1[ih] += (factor) *  (in2->c1[ih]);
		in1->c2[ih] += (factor) *  (in2->c2[ih]);
	}
}




void combine_in1xferm_mass2_minus_in2_minus_in3(__restrict const vec3_soa * const in_vect1,
																								double ferm_mass,
																								__restrict const vec3_soa * const in_vect2,
																								__restrict const vec3_soa * const in_vect3,
																								__restrict vec3_soa * const out)
{
	#pragma acc kernels present(in_vect1) present(in_vect2) present(in_vect3) present(out)
	#pragma acc loop independent 
	for(int ih=LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih < LNH_SIZEH-LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih++) {
		out->c0[ih]=(in_vect1->c0[ih]*ferm_mass)-in_vect2->c0[ih]-in_vect3->c0[ih];
		out->c1[ih]=(in_vect1->c1[ih]*ferm_mass)-in_vect2->c1[ih]-in_vect3->c1[ih];
		out->c2[ih]=(in_vect1->c2[ih]*ferm_mass)-in_vect2->c2[ih]-in_vect3->c2[ih];
	}
}


void combine_inside_loop(__restrict vec3_soa * const vect_out,
												 __restrict vec3_soa * const vect_r,
												 const __restrict vec3_soa * const vect_s,
												 const __restrict vec3_soa * const vect_p,
												 const double omega)
{
	#pragma acc kernels present(vect_out) present(vect_r) present(vect_s) present(vect_p)
	#pragma acc loop independent 
	for(int ih=LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih < LNH_SIZEH-LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih++) {
		// out+=omega*p
		vect_out->c0[ih] += (vect_p->c0[ih]*omega);
		vect_out->c1[ih] += (vect_p->c1[ih]*omega);
		vect_out->c2[ih] += (vect_p->c2[ih]*omega);
		// r-=omega*s
		vect_r->c0[ih]   -= (vect_s->c0[ih]*omega);
		vect_r->c1[ih]   -= (vect_s->c1[ih]*omega);
		vect_r->c2[ih]   -= (vect_s->c2[ih]*omega);
	}
}

void combine_in1xferm_mass_minus_in2(__restrict const vec3_soa * const in_vect1,
																		 double ferm_mass2,
																		 __restrict vec3_soa * const in_vect2)
{
	#pragma acc kernels present(in_vect1) present(in_vect2)
	#pragma acc loop independent
	for(int ih=LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih < LNH_SIZEH-LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih++) {
		in_vect2->c0[ih]=(in_vect1->c0[ih]*ferm_mass2)-in_vect2->c0[ih];
		in_vect2->c1[ih]=(in_vect1->c1[ih]*ferm_mass2)-in_vect2->c1[ih];
		in_vect2->c2[ih]=(in_vect1->c2[ih]*ferm_mass2)-in_vect2->c2[ih];
	}
}

void combine_in1_minus_in2(__restrict const vec3_soa * in_vect1,
													 __restrict const vec3_soa * in_vect2,
													 __restrict vec3_soa * out)
{
	#pragma acc kernels present(in_vect1) present(in_vect2) present(out)
	#pragma acc loop independent 
	for(int ih=LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih < LNH_SIZEH-LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih++) {
		out->c0[ih]=(in_vect1->c0[ih])-in_vect2->c0[ih];
		out->c1[ih]=(in_vect1->c1[ih])-in_vect2->c1[ih];
		out->c2[ih]=(in_vect1->c2[ih])-in_vect2->c2[ih];
	}
}

void assign_in_to_out(__restrict const vec3_soa * in_vect1,
											__restrict vec3_soa * out)
{
	#pragma acc kernels present(in_vect1)  present(out)
	#pragma acc loop independent 
	for(int ih=LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih < LNH_SIZEH-LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih++) {
		out->c0[ih]=(in_vect1->c0[ih]);
		out->c1[ih]=(in_vect1->c1[ih]);
		out->c2[ih]=(in_vect1->c2[ih]);
	}

}

// functions added for linear algebra in multishift "versatilizzato"
void set_vec3_soa_to_zero(__restrict vec3_soa* const fermion)
{
	//  printf(" indicted pointer (%p) \n",fermion);
	#pragma acc kernels present(fermion)
	#pragma acc loop independent 
	// printf("Setting to zero memory from %p to %p ", fermion->c0 , &(fermion->c0[sizeh-1])+1);
	// printf(", %p to %p  ", fermion->c1 , &(fermion->c1[sizeh-1])+1);
	// printf(", %p to %p\n", fermion->c2 , &(fermion->c2[sizeh-1])+1);
	for(int i= 0; i < sizeh ; i++){
		fermion->c0[i]=0;
		fermion->c1[i]=0;
		fermion->c2[i]=0;
	}
}
void multiple_combine_in1_minus_in2x_factor_back_into_in1(__restrict vec3_soa * const out, 
																													__restrict const vec3_soa * const in, 
																													const int maxiter, 
																													__restrict const int * const flag, 
																													__restrict const double * const omegas)
{
	int ia, ih;

	#pragma acc kernels present(in) present(out) copyin(omegas[0:maxiter]) copyin(flag[0:maxiter]) async(4) 
	{
		#pragma acc cache(omegas[0:maxiter])
		#pragma acc cache(flag[0:maxiter])
		#pragma acc loop independent gang
		for (ia=0; ia<maxiter; ia++) {
			#pragma acc loop independent gang vector(512)
			for (ih=LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih < LNH_SIZEH-LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih++) {
				if (flag[ia] == 1) {
					double factor = omegas[ia];
					out[ia].c0[ih] -= (factor)*(in[ia].c0[ih]);
					out[ia].c1[ih] -= (factor)*(in[ia].c1[ih]);
					out[ia].c2[ih] -= (factor)*(in[ia].c2[ih]);
				} // if
			} // ih
		} // ia
	} // acc kernels
}

void multiple1_combine_in1_x_fact1_plus_in2_x_fact2_back_into_in1(__restrict vec3_soa * const in1,
																																	int maxiter, 
																																	__restrict const int * const flag, 
																																	__restrict const double * const gammas, 
																																	__restrict const vec3_soa * const in2, 
																																	__restrict const double * const zeta_iii ) 
{

	int ia, ih;

	#pragma acc kernels present(in1) present(in2) copyin(gammas[0:maxiter]) copyin(flag[0:maxiter]) copyin(zeta_iii[0:maxiter])
	{

		#pragma acc cache(gammas[0:maxiter]) 
		#pragma acc cache(zeta_iii[0:maxiter])
		#pragma acc cache(flag[0:maxiter]) 

		#pragma acc loop independent gang
		for (ia=0; ia<maxiter; ia++) {
			#pragma acc loop independent gang vector(512)
			for(ih=LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih < LNH_SIZEH-LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih++) {

				if (flag[ia] == 1) {

					double fact1 = gammas[ia];
					double fact2 = zeta_iii[ia];

					in1[ia].c0[ih] = fact1 * (in1[ia].c0[ih]) + fact2 * (in2->c0[ih]);
					in1[ia].c1[ih] = fact1 * (in1[ia].c1[ih]) + fact2 * (in2->c1[ih]);
					in1[ia].c2[ih] = fact1 * (in1[ia].c2[ih]) + fact2 * (in2->c2[ih]);

				} // if flag
			} // ih
		} // ia
	} // acc data

}
void combine_in1_x_fact1_minus_in2_back_into_in2(__restrict const vec3_soa * const in1, 
																								 double fact1, 
																								 __restrict vec3_soa * const in2 )
{

	int ih;


	#pragma acc kernels present(in1) present(in2) 
	{


		#pragma acc loop independent gang vector(512)
		for(ih=LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih < LNH_SIZEH-LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih++) {


			in2->c0[ih] = fact1 * (in1->c0[ih]) - (in2->c0[ih]);
			in2->c1[ih] = fact1 * (in1->c1[ih]) - (in2->c1[ih]);
			in2->c2[ih] = fact1 * (in1->c2[ih]) - (in2->c2[ih]);

		} // ih
	} // acc data

}
void combine_in1_minus_in2_allxfact(__restrict const vec3_soa * const in1, 
																		__restrict const vec3_soa * const in2, 
																		double fact,
																		__restrict vec3_soa * const out )
{
	int ih;
	#pragma acc kernels present(in1) present(in2) present(out) 
	{
		#pragma acc loop independent gang vector(512)
		for(ih=LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih < LNH_SIZEH-LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih++) {
			out->c0[ih] = fact * (in1->c0[ih] - in2->c0[ih]);
			out->c1[ih] = fact * (in1->c1[ih] - in2->c1[ih]);
			out->c2[ih] = fact * (in1->c2[ih] - in2->c2[ih]);
		} // ih
	} // acc data
}

void calc_new_trialsol_for_inversion_in_force(int halfLen,
																							__restrict vec3_soa * inout, int nPrecCalculations){
	int ih,iv;
	#pragma acc kernels present(inout)
	{


		if(1 == nPrecCalculations%2){    
			#pragma acc loop independent gang            
			for(iv=0; iv < halfLen; iv++){
				#pragma acc loop independent gang vector(512)
				for(ih=LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih < LNH_SIZEH-LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih++) {
					inout[iv+halfLen].c0[ih]=2*inout[iv].c0[ih]-inout[iv+halfLen].c0[ih];
					inout[iv+halfLen].c1[ih]=2*inout[iv].c1[ih]-inout[iv+halfLen].c1[ih];
					inout[iv+halfLen].c2[ih]=2*inout[iv].c2[ih]-inout[iv+halfLen].c2[ih];
					// ^^ next trial         ^^ last solution    ^^ old solution
				} // ih
			}// iv
		}// if


		if(0 == nPrecCalculations%2){   // nPrecCalculations even : next trial 
			// written in the first half of the vector     
			#pragma acc loop independent gang
			for(iv=0; iv < halfLen; iv++){
				#pragma acc loop independent gang vector(512)
				for(ih=LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih < LNH_SIZEH-LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih++) {
					inout[iv].c0[ih]=2*inout[iv+halfLen].c0[ih]-inout[iv].c0[ih];
					inout[iv].c1[ih]=2*inout[iv+halfLen].c1[ih]-inout[iv].c1[ih];
					inout[iv].c2[ih]=2*inout[iv+halfLen].c2[ih]-inout[iv].c2[ih];
					// ^^ next trial         ^^ last solution    ^^ old solution
				} // ih
			} // iv
		} // if
	} // acc data
}


#endif
