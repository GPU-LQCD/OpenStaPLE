#ifndef COOLING_C
#define COOLING_C

#include "../Include/common_defines.h"
#include "./struct_c_def.h"
#include "../OpenAcc/su3_utilities.h"
#include "../OpenAcc/plaquettes.h"
#include "./cooling.h"
#ifdef __GNUC__
#include <math.h>
#else // assuming PGI is used for compilation  on accelerators
#include <accelmath.h>
#endif

#ifdef MULTIDEVICE
#include "../Mpi/communications.h"
#endif

#pragma acc routine seq
static inline void cabibbo_marinari_cooling(__restrict su3_soa   * const U, // unsmeared input conf
																						__restrict su3_soa   * const Ucool, // smeared output conf
																						__restrict su3_soa   * const STAP, // staples
																						int idx){

  // link initialization
  d_complex U00 = U->r0.c0[idx];
  d_complex U01 = U->r0.c1[idx];
  d_complex U02 = U->r0.c2[idx];
  d_complex U10 = U->r1.c0[idx];
  d_complex U11 = U->r1.c1[idx];
  d_complex U12 = U->r1.c2[idx];

  // Compute 3rd U row from the first two
  d_complex U20 = conj( ( U01 * U12 ) - ( U02 * U11) ) ;
  d_complex U21 = conj( ( U02 * U10 ) - ( U00 * U12) ) ;
  d_complex U22 = conj( ( U00 * U11 ) - ( U01 * U10) ) ;

  // staple initialization
  d_complex S00 = STAP->r0.c0[idx]/C_ZERO;
  d_complex S01 = STAP->r0.c1[idx]/C_ZERO;
  d_complex S02 = STAP->r0.c2[idx]/C_ZERO;
  d_complex S10 = STAP->r1.c0[idx]/C_ZERO;
  d_complex S11 = STAP->r1.c1[idx]/C_ZERO;
  d_complex S12 = STAP->r1.c2[idx]/C_ZERO;
  d_complex S20 = STAP->r2.c0[idx]/C_ZERO;
  d_complex S21 = STAP->r2.c1[idx]/C_ZERO;
  d_complex S22 = STAP->r2.c2[idx]/C_ZERO;


  //U*P complete product
  /*
		d_complex P00 = U00*S00 + U01*S10 + U02*S20;
		d_complex P01 = U00*S01 + U01*S11 + U02*S21;
		d_complex P02 = U00*S02 + U01*S12 + U02*S22;
		d_complex P10 = U10*S00 + U11*S10 + U12*S20;
		d_complex P11 = U10*S01 + U11*S11 + U12*S21;
		d_complex P12 = U10*S02 + U11*S12 + U12*S22;
		d_complex P20 = U20*S00 + U21*S10 + U22*S20;
		d_complex P21 = U20*S01 + U21*S11 + U22*S21;
		d_complex P22 = U20*S02 + U21*S12 + U22*S22;
  */
  /*
		double c0,c1,c2,c3;
		c0 = creal(Pii) + creal(Pjj);
		c1 = cimag(Pij) + cimag(Pji);
		c2 = creal(Pij) - creal(Pji);
		c3 = cimag(Pii) - cimag(Pjj);
  */


  // first subgroup (i=0,j=1)
  d_complex Pii = U00*S00 + U01*S10 + U02*S20;
  d_complex Pij = U00*S01 + U01*S11 + U02*S21;
  d_complex Pji = U10*S00 + U11*S10 + U12*S20;
  d_complex Pjj = U10*S01 + U11*S11 + U12*S21;
  d_complex A = (creal(Pii) + creal(Pjj)) - (cimag(Pii) - cimag(Pjj))*I;
  d_complex B = -(creal(Pij) - creal(Pji)) - (cimag(Pij) + cimag(Pji))*I;
  double norm = 1.0/sqrt(creal(A)*creal(A)+cimag(A)*cimag(A)+creal(B)*creal(B)+cimag(B)*cimag(B));
  A = A * norm;
  B = B * norm;
  // T matrix is U link cooled along the first subgroup
  d_complex T00 = A * U00 + B * U10;
  d_complex T01 = A * U01 + B * U11;
  d_complex T02 = A * U02 + B * U12;
  d_complex T10 = -conj(B) * U00 + conj(A) * U10;
  d_complex T11 = -conj(B) * U01 + conj(A) * U11;
  d_complex T12 = -conj(B) * U02 + conj(A) * U12;
  d_complex T20 = U20;
  d_complex T21 = U21;
  d_complex T22 = U22;

  // second subgroup (i=0,j=2)
  Pii = T00*S00 + T01*S10 + T02*S20;
  Pij = T00*S02 + T01*S12 + T02*S22;
  Pji = T20*S00 + T21*S10 + T22*S20;
  Pjj = T20*S02 + T21*S12 + T22*S22;
  A = (creal(Pii) + creal(Pjj)) - (cimag(Pii) - cimag(Pjj))*I;
  B = -(creal(Pij) - creal(Pji)) - (cimag(Pij) + cimag(Pji))*I;
  norm = 1.0/sqrt(creal(A)*creal(A)+cimag(A)*cimag(A)+creal(B)*creal(B)+cimag(B)*cimag(B));
  A = A * norm;
  B = B * norm;
  // U matrix is T link cooled along the second subgroup
  U00 = A * T00 + B * T20;
  U01 = A * T01 + B * T21;
  U02 = A * T02 + B * T22;
  U10 = T10;
  U11 = T11;
  U12 = T12;
  U20 = -conj(B) * T00 + conj(A) * T20;
  U21 = -conj(B) * T01 + conj(A) * T21;
  U22 = -conj(B) * T02 + conj(A) * T22;


  // third subgroup (i=1,j=2)
  Pii = U10*S01 + U11*S11 + U12*S21;
  Pij = U10*S02 + U11*S12 + U12*S22;
  Pji = U20*S01 + U21*S11 + U22*S21;
  Pjj = U20*S02 + U21*S12 + U22*S22;
  A = (creal(Pii) + creal(Pjj)) - (cimag(Pii) - cimag(Pjj))*I;
  B = -(creal(Pij) - creal(Pji)) - (cimag(Pij) + cimag(Pji))*I;
  norm = 1.0/sqrt(creal(A)*creal(A)+cimag(A)*cimag(A)+creal(B)*creal(B)+cimag(B)*cimag(B));
  A = A * norm;
  B = B * norm;
  // T matrix is U link cooled along the third subgroup
  Ucool->r0.c0[idx] = U00;
  Ucool->r0.c1[idx] = U01;
  Ucool->r0.c2[idx] = U02;
  Ucool->r1.c0[idx] = A * U10 + B * U20;
  Ucool->r1.c1[idx] = A * U11 + B * U21;
  Ucool->r1.c2[idx] = A * U12 + B * U22;
  Ucool->r2.c0[idx] = -conj(B) * U10 + conj(A) * U20;
  Ucool->r2.c1[idx] = -conj(B) * U11 + conj(A) * U21;
  Ucool->r2.c2[idx] = -conj(B) * U12 + conj(A) * U22;
  
}



void compute_cooled_even_links(__restrict su3_soa   * const U, // unsmeared input conf
															 __restrict su3_soa   * const Ucool, // smeared output conf
															 __restrict su3_soa   * const STAP){

  int dirindex,idxh;
	#pragma acc kernels present(U) present(STAP)
	#pragma acc loop independent
  for(dirindex = 0 ; dirindex < 4 ; dirindex++){
    const int dir_link = 2*dirindex; // only even sites
		#pragma acc loop independent
    for( idxh = 0 ; idxh < sizeh; idxh++){
      cabibbo_marinari_cooling(&U[dir_link],&Ucool[dir_link],&STAP[dir_link],idxh);
    }
  }
}

void compute_cooled_odd_links(__restrict su3_soa   * const U, // unsmeared input conf
															__restrict su3_soa   * const Ucool, // smeared output conf
															__restrict su3_soa   * const STAP){
  
  int dirindex,idxh;
	#pragma acc kernels present(U) present(STAP)
	#pragma acc loop independent
  for(dirindex = 0 ; dirindex < 4 ; dirindex++){
    const int dir_link = 2*dirindex+1; // only odd sites
		#pragma acc loop independent
    for( idxh = 0 ; idxh < sizeh; idxh++){
      cabibbo_marinari_cooling(&U[dir_link],&Ucool[dir_link],&STAP[dir_link],idxh);
    }
  }
}


void cool_conf(__restrict su3_soa   * const U, // unsmeared input conf
							 __restrict su3_soa   * const Ucool, // smeared output conf
							 __restrict su3_soa   * const TMP){ // auxiliary staple container

  set_su3_soa_to_zero(TMP);
#if NRANKS_D3 > 1	//#ifdef MULTIDEVICE
  communicate_su3_borders(U, GAUGE_HALO);  
#endif
 
  calc_loc_staples_nnptrick_all_only_even(U,TMP);

#if NRANKS_D3 > 1	//#ifdef MULTIDEVICE
  communicate_gl3_borders(TMP, GAUGE_HALO);  
#endif

  compute_cooled_even_links(U,Ucool,TMP);
  calc_loc_staples_nnptrick_all_only_odd(U,TMP);

#if NRANKS_D3 > 1	//#ifdef MULTIDEVICE
  communicate_gl3_borders(TMP, GAUGE_HALO);  
#endif  

  compute_cooled_odd_links(U,Ucool,TMP);
  unitarize_conf(Ucool);

#if NRANKS_D3 > 1	//#ifdef MULTIDEVICE
  communicate_su3_borders(Ucool, GAUGE_HALO);  
#endif


}


#endif
