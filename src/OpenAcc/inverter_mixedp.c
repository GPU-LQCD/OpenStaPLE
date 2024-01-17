#ifndef INVERTER_MIXEDP_C_
#define INVERTER_MIXEDP_C_

#include "../Include/common_defines.h"
#include "../Include/inverter_tricks.h"
#include "./fermion_matrix.h"
#include "./fermionic_utilities.h"
#include "./float_double_conv.h"
#include "./sp_fermion_matrix.h"
#include "./sp_fermionic_utilities.h"
#include "./sp_struct_c_def.h"
#include "./struct_c_def.h"
#include "./inverter_package.h"

#ifndef __GNUC__
#include "openacc.h"
#define __restrict
#endif
#include "./sp_inverter_full.h"

#include "../Mpi/multidev.h"


void combine_add_in2_into_in1_mixed_precision(__restrict vec3_soa * in1,__restrict const  vec3_soa_f * in2)
{
	#pragma acc kernels present(in1) present(in2)
	#pragma acc loop independent
  for(int ih=LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih < LNH_SIZEH-LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih++) {
    in1->c0[ih] += (double)creal(in2->c0[ih])+ I* (double) cimag(in2->c0[ih]);
    in1->c1[ih] += (double)creal(in2->c1[ih])+ I* (double) cimag(in2->c1[ih]);
    in1->c2[ih] += (double)creal(in2->c2[ih])+ I* (double) cimag(in2->c2[ih]);
  }

}


// based on https://arxiv.org/pdf/0911.3191v2.pdf, pg.20

#define SAFETY_MARGIN 0.9

int inverter_mixed_precision(inverter_package ip,
														 ferm_param *pars,
														 __restrict vec3_soa * solution, // double precision output
														 __restrict const vec3_soa * in, // constant
														 double res,
														 const int  max_cg,
														 double shift,
														 int * cg_return)
{

  int cg;
  int magicTouchCount = 0 ;
  long int i;
  double delta, alpha, lambda, omega, gammag;
  double lastMaxResNorm=0;
  double source_norm = l2norm2_global(in);
  
  int printevery = 10000;
  int qmax = (verbosity_lv < 13) ? verbosity_lv : 13;
  for(int q=0;q<qmax;q++) printevery /= 2;
  if (verbosity_lv > 3 && 0==devinfo.myrank )
		{ 
      printf("source_norm: %e, ",source_norm);
      printf("res = %e, ",res);
      printf("target residue = %e\n", res*res*source_norm);
		}

  const su3_soa_f * u = ip.u_f;
  __restrict vec3_soa_f * loc_r = ip.loc_r_f ;
  __restrict vec3_soa_f * loc_h = ip.loc_h_f ;
  __restrict vec3_soa_f * loc_s = ip.loc_s_f ;
  __restrict vec3_soa_f * loc_p = ip.loc_p_f ;
  __restrict vec3_soa_f * out = ip.out_f ;

  // r = in  - M * solution
  // s = M * solution
  fermion_matrix_multiplication_shifted(ip.u,ip.loc_s,solution,ip.loc_h,pars,shift);
  // r = in - s
  combine_in1_minus_in2(in,ip.loc_s,ip.loc_r);
  
  // \hat{r} = (float)r
  convert_double_to_float_vec3_soa(ip.loc_r,loc_r);

  // starting with CG, Single precision
  assign_in_to_out_f(loc_r,loc_p);
  
  delta=l2norm2_global_f(loc_r);

  // \hat{x} = 0 // 
  set_vec3_soa_to_zero_f(out);

  // loop over cg iterations
  cg=0;
	if (verbosity_lv > 3 && 0==devinfo.myrank )
		printf("STARTING CG:\nCG\tR - mixed precision\n");

  do {
    
    cg++;    
    // s=(M^dag M)p    alpha=(p,s)
    fermion_matrix_multiplication_shifted_f(u,loc_s,loc_p,loc_h,pars,shift);
    alpha = real_scal_prod_global_f(loc_p,loc_s);

    omega=delta/alpha;     
    // solution+=omega*p  r-=omega*s
    // lambda=(r,r);
    
    combine_in1xfactor_plus_in2_f(loc_p,omega,out,out);
  
    // calculation of loc_r_f
    // EITHER THIS WAY 
    if(lastMaxResNorm < delta) lastMaxResNorm = delta;
    // the double precision magic touch
    if(delta < inverter_tricks.mixedPrecisionDelta * lastMaxResNorm){
			// adding the partial solution back to the full solution
			combine_add_in2_into_in1_mixed_precision(solution,out);
			// r = in  - M * solution
			// s = M * solution
			fermion_matrix_multiplication_shifted(ip.u,ip.loc_s,solution,ip.loc_h,pars,shift);
			// r = in - s
			combine_in1_minus_in2(in,ip.loc_s,ip.loc_r);

			// \hat{r} = (float)r
			convert_double_to_float_vec3_soa(ip.loc_r,loc_r);

			// \hat{x} = 0 // 
			set_vec3_soa_to_zero_f(out);

			lastMaxResNorm=0;
			magicTouchCount++;


			// OR THIS WAY
    }else combine_in1xfactor_plus_in2_f(loc_s,-omega,loc_r,loc_r);



    lambda = l2norm2_global_f(loc_r);
    gammag=lambda/delta;
    delta=lambda;

    if(0==devinfo.myrank && cg%printevery==0 )
			{
        printf("Iteration %d, Inverter Mixed Precision:", cg);
        printf("residue norm: %1.3e ",lambda);
        printf("(%d \"magic touches\" until now)\n", magicTouchCount);
			}




    // p=r+gammag*p
    combine_in1xfactor_plus_in2_f(loc_p,gammag,loc_r,loc_p);

  } while( (sqrt(lambda/source_norm)>res*SAFETY_MARGIN) && cg<max_cg);

  combine_add_in2_into_in1_mixed_precision(solution,out);

  if (verbosity_lv > 3  && 0==devinfo.myrank ) printf("\n");
#if ((defined DEBUG_MODE) || (defined DEBUG_INVERTER_FULL_OPENACC))
  // s = M x
  fermion_matrix_multiplication_shifted(ip.u,ip.loc_s,solution,ip.loc_h,pars,shift);
  combine_in1_minus_in2(in,ip.loc_s,ip.loc_h); // r = s - y  
  double  giustoono=l2norm2_global(ip.loc_h)/source_norm;
  if(verbosity_lv > 1 && 0==devinfo.myrank  ){
		printf("Terminated invert after   %d    iterations", cg);
		printf("[res/stop_res=  %e , stop_res=%e ]\n",
					 sqrt(giustoono)/res,res);
  }
#endif
  if(cg==max_cg  && 0==devinfo.myrank )
    {
      printf("WARNING: maximum number of iterations reached in invert\n");
    }
  
	*cg_return = cg;
	if (sqrt(giustoono) <= res)
		return INVERTER_SUCCESS;
	else return INVERTER_FAILURE;

}


#endif
