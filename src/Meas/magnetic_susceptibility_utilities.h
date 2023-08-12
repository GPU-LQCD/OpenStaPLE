#ifndef MAGNETIC_SUSCEPTIBILITY_UTILITIES_H_
#define MAGNETIC_SUSCEPTIBILITY_UTILITIES_H_

#include "../OpenAcc/struct_c_def.h"
#include "../Include/fermion_parameters.h"

// calculation of the derivative of the phase with respect to parameter bz.
// basically, just I*( phases(b=(0,0,1)) - phases(b=(0,0,0)))
// linearity is assumed in the field.
void idphase_dbz(double_soa * dphi_dbz_re, double_soa * dphi_dbz_im,  // dphi_dbz_im = 0
								 ferm_param *tfermion_parameters,int nnp_openacc[sizeh][4][2]);


// FINITE differences in e^(i phase_1) - e^(i phase_0)
void delta_u_delta_bz(double_soa * delta_u_re,double_soa * delta_u_im,
        double_soa * tp1, double_soa * tp2,
        ferm_param *tfermion_parameters, bf_param ext_field);


#endif
