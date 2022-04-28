#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "./rationalapprox.h"
#include"../src/RationalApprox/Remez/bigfloat.h"
#include"../src/RationalApprox/Remez/alg_remez.h"


//using namespace std;

const int gmp_remez_precision=256; // The precision that gmp uses

//Rational approximations for Metropolis
const double lambda_min=4.0e-7;  // rational approx valid on [lambda_min_metro, 1.0]
const double residue=1.0e-8;    // stopping residual for CG

int verbosity_lv = 5 ; 

int main(int argc, char **argv){
  if(argc!=5){
    printf("Use as arguments:    approx_order exponent_num exponent_den lambda_min\n");
    return 0;
  }

  RationalApprox* approx = new RationalApprox;

  approx->approx_order = atoi(argv[1]);
  approx->exponent_num = atoi(argv[2]);
  approx->exponent_den = atoi(argv[3]);
  approx->gmp_remez_precision = gmp_remez_precision; // The precision that gmp uses
  approx->lambda_min = atof(argv[4]);
  approx->lambda_max = 1.0;              // The bounds of the approximation

  // The error from the approximation (the relative error is minimised
  // - if another error minimisation is requried, then line 398 in
  // alg_remez.C is where to change it)

  // The partial fraction expansion takes the form 
  // r(x) = norm + sum_{k=1}^{n} res[k] / (x + pole[k])
 
 
  printf("\nApproximation to f(x) = (x)^(%d/%d) on [%e, %e]\n\n",approx->exponent_num, approx->exponent_den,approx->lambda_min, approx->lambda_max);
  fflush(stdout);
  // Instantiate the Remez class
  AlgRemez remez1(approx->lambda_min,approx->lambda_max,approx->gmp_remez_precision);
  // Generate the required approximation
  approx->error = remez1.generateApprox(approx->approx_order,approx->approx_order,abs(approx->exponent_num),abs(approx->exponent_den));
  printf("approximation error = %e\n\n", approx->error);

  // Find the partial fraction expansion of the approximation 
  // to the function x^{y/z} (this only works currently for 
  // the special case that n = d)

  if(approx->exponent_num*approx->exponent_den > 0)
      remez1.getPFE(approx->RA_a,approx->RA_b,&(approx->RA_a0));
  else
      remez1.getIPFE(approx->RA_a,approx->RA_b,&(approx->RA_a0));

  char * nomefile = rational_approx_filename(approx->approx_order, approx->exponent_num, approx->exponent_den, approx->lambda_min);

  rationalapprox_save(nomefile, approx);

}
