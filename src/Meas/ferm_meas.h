// use z2 noise instead of gaussian noise (see hep-lat/9308015)
// use the global defined fermions loc_chi, loc_phi, rnd_o, rnd_e, chi_o and loc_h
#ifndef FERM_MEAS_H 
#define FERM_MEAS_H

#include <stdio.h>
#include "../OpenAcc/struct_c_def.h"
#include "../OpenAcc/inverter_package.h"
#include "../Include/fermion_parameters.h"


typedef struct ferm_meas_param_t{

    int measEvery;
    char fermionic_outfilename[50];
    char fermionic_outfile_header[1000];
    int SingleInvNVectors, DoubleInvNVectorsChiral;
    int DoubleInvNVectorsQuarkNumber;
    int printPlaqAndRect ;
} ferm_meas_params;

extern ferm_meas_params fm_par;




void eo_inversion(inverter_package ip,
		  ferm_param * tfermions_parameters,
          double res, int max_cg,
		  vec3_soa * in_e,     // z2 noise
		  vec3_soa * in_o,     // z2 noise
		  vec3_soa * out_e,
		  vec3_soa * out_o,
		  vec3_soa * phi_e,    // parking variable
		  vec3_soa * phi_o);   // parking variable


void fermion_measures(su3_soa * tconf_acc,
									ferm_param * tfermions_parameters, 
                  ferm_meas_params * tfm_par,
                  double res, int max_cg, int conf_id_iter,
									int nnp_openacc[sizeh][4][2],int nnm_openacc[sizeh][4][2],
                  double plaq, double rect);


#endif
