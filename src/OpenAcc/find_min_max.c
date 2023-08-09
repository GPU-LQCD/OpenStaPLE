#ifndef FIND_MIN_MAX_C_
#define FIND_MIN_MAX_C_


#include "./struct_c_def.h"
#include "./fermionic_utilities.h"
#include "./fermion_matrix.h"
#include <stdio.h>
#include "./find_min_max.h"
#include "../Mpi/multidev.h"

// #ifndef __GNUC__
#include "openacc.h"
// #endif


extern int verbosity_lv;

// find the maximum eigenvalue of the fermion matrix
// use loc_h, loc_p, loc_r
double ker_find_max_eigenvalue_openacc(  __restrict su3_soa * const u,
																				 ferm_param *pars,
																				 __restrict vec3_soa * const loc_r,
																				 __restrict vec3_soa * const loc_h,
																				 __restrict vec3_soa * const loc_p
																				 ){

  int loop_count;
  double norm, inorm, old_norm;
  double temp;
  // starting gauss vector p (it has to be gaussian yet)
  norm=sqrt(l2norm2_global(loc_p));
  if(verbosity_lv>4) printf("MPI%02d: (ker_find_max_eigenvalue_openacc) Norm:%lf\n",devinfo.myrank,norm );
  loop_count=0;
  // loop start
  do{
		// normalize  p
		// and r=p
		inorm=1.0/norm;
		multiply_fermion_x_doublefactor(loc_p,inorm);
		assign_in_to_out(loc_p,loc_r);

		old_norm=norm;
		// p=DeoDoe r
      
      
		fermion_matrix_multiplication(u,loc_p,loc_r,loc_h,pars);


		norm=sqrt(l2norm2_global(loc_p));
		old_norm=fabs(old_norm-norm);
		old_norm/=norm;
		loop_count++;
		if(loop_count %100 == 0 && verbosity_lv > 2 && 0 == devinfo.myrank )
			printf("        iterations of max computat  = %d, norm = %lf, old_norm = %lf\n",
						 loop_count, norm, old_norm);
  } while(old_norm>1.0e-5); // loop end
  double max=norm;
  return max;
}

// find the minimum eigenvalue of the fermion matrix
// use loc_h, loc_p, loc_r
double ker_find_min_eigenvalue_openacc(  __restrict su3_soa * const u,
																				 __restrict ferm_param * const pars,
																				 __restrict vec3_soa * const loc_r,
																				 __restrict vec3_soa * const loc_h,
																				 __restrict vec3_soa * const loc_p,
																				 double max
																				 ){
  int loop_count;
  double norm, inorm, old_norm, delta;
  double m2 =(pars->ferm_mass)*(pars->ferm_mass);
  delta=max-m2;
  norm=sqrt(l2norm2_global(loc_p));
  loop_count=0;
  // loop start
  do{
    // normalize  p 
    // and r=p 
    inorm=1.0/norm;
    multiply_fermion_x_doublefactor(loc_p,inorm);
    // r=p
    assign_in_to_out(loc_p,loc_r);
    old_norm=norm;
    // p=DeoDoe r 
    // p=max r - (M^dag M)r 
    fermion_matrix_multiplication_shifted(u,loc_p,loc_r,loc_h,pars,delta-m2);
    norm=sqrt(l2norm2_global(loc_p));
    old_norm=fabs(old_norm-norm);
    old_norm/=norm;
    loop_count++;
  }  while(old_norm>1.0e-5); // loop end
  double  min=max-norm;
  return min;
}


void find_min_max_eigenvalue_soloopenacc(  __restrict su3_soa * const u,
																					 ferm_param *pars,
																					 __restrict vec3_soa * const loc_r,
																					 __restrict vec3_soa * const loc_h,
																					 __restrict vec3_soa * const loc_p1,
																					 __restrict vec3_soa * const loc_p2,
																					 double *minmax
																					 ){
  // minmax[0] --> minimum
  // minmax[1] --> maximum

  minmax[0] = pars->ferm_mass * pars->ferm_mass;
  minmax[1] = ker_find_max_eigenvalue_openacc(u,pars,loc_r,loc_h,loc_p1);

  if(verbosity_lv > 3 && 0 == devinfo.myrank) printf("    Computed the min and max eigs : [ min= %.18lf; max= %.18lf ]  \n",minmax[0],minmax[1]);
  //  now minimum is m*m, it could be computed using the following.
  //  minmax[0] = ker_find_min_eigenvalue_openacc(u,backfield,pars,loc_r,loc_h,loc_p2,minmax[1]); //--> it could be used mass2

}



#endif

