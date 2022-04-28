#ifndef FERMION_PARAMETERS_H
#define FERMION_PARAMETERS_H

#include "../RationalApprox/rationalapprox.h"
#include "../OpenAcc/struct_c_def.h"
#include "../OpenAcc/sp_struct_c_def.h"
#include "../OpenAcc/backfield.h"

typedef struct ferm_param_t{
  double ferm_mass;
  int degeneracy;
  int number_of_ps;
  char name[30];
  // chem_pot or backfield related things
  double ferm_charge;
  double ferm_im_chem_pot;
  // automatic from here on
  int index_of_the_first_ps;
  int index_of_the_first_shift;
  double_soa * phases; //this incorporates staggered phases,
                       // external u(1) fields and 
                       // imaginary chemical potential
                      
  double_soa * mag_re;  // Real and impaginary  parts of the 'algebra-prefix'
  double_soa * mag_im;  // of magnetization observable 
                            // 8 for each flavour

  int printed_bf_dbg_info;

  float_soa * phases_f;//this incorporates staggered phases,
                       // external u(1) fields and 
                       // imaginary chemical potential




  RationalApprox approx_fi_mother; // first inv   -> mother
  RationalApprox approx_md_mother; // md approx   -> mother
  RationalApprox approx_li_mother; // last inv    -> mother
  RationalApprox approx_fi;        // first inv
  RationalApprox approx_md;        // md approx
  RationalApprox approx_li;        // last inv
} ferm_param;


extern int max_ps;
extern int totalMdShifts;
extern ferm_param *fermions_parameters;


//FERMION PARAMETERS


int init_ferm_params(ferm_param * fermion_settings); 
// returns != 0 if errors are encountered.

// to be performed after all allocations (back fields should already be 
// allocated )
void init_all_u1_phases(bf_param bfpars, ferm_param *fpar);


void init_fermion_backfield(bf_param bf_pars, ferm_param *fermion_parameters);

int rat_approx_file_or_script_create(RationalApprox* rational_approx);

#endif

