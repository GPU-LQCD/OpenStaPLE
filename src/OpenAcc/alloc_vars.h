#ifndef ALLOC_VARS_H_
#define ALLOC_VARS_H_


#include "../Include/common_defines.h"
#include "./struct_c_def.h"

// used in the dynamical allocation of structures


extern global_su3_soa  * conf_rw; // the gauge configuration, only for read-write

// ued in debugging/testing
extern global_vec3_soa  * ferm_rw; // a global fermion, only for read-write
extern global_tamat_soa  * tamat_rw; // a global tamat, only for read-write
extern global_thmat_soa  * thmat_rw; // a global thmat, only for read-write
extern global_dcomplex_soa *dcomplex_rw; // a global dcomplex_soa, only for read-write
extern global_double_soa *double_rw; // a global double_soa, only for read-write

extern su3_soa  * conf_acc_bkp; // the old stored conf that will be recovered if the metro test fails.

extern su3_soa  * aux_conf_acc; // auxiliary 
extern su3_soa  * auxbis_conf_acc; // auxiliary

extern double_soa * u1_back_phases; // background,staggered,chempot phases
                                    // 8 for each flavour

extern double_soa * mag_obs_re;     // real part of the 'algebra-prefix'
                                    // of magnetization observable 
                                    // 8 for each flavour

extern double_soa * mag_obs_im;     // imaginary part of the 'algebra-prefix'
                                    // of magnetization observable 
                                    // 8 for each flavour

extern double_soa * topo_loc; // topological charge auxiliary

extern thmat_soa * momenta; // gauge field evolution
extern int momenta_backupped;
extern thmat_soa * momenta_backup; // gauge field evolution - reversibility test
extern tamat_soa * ipdot_acc; // gauge field evolution
extern tamat_soa * ipdot_g_old; // for HMC diagnostics
extern tamat_soa * ipdot_f_old; // for HMC diagnostics



extern su3_soa * gconf_as_fermionmatrix; // (only a pointer) conf to use in either cases 
                                         // in fermion related computation (with or without stouting)

extern su3_soa * conf_acc; // gauge configuration(s)

// STOUTING 
extern su3_soa * gstout_conf_acc_arr; // all stouting steps except the zeroth
extern su3_soa * glocal_staples;
extern tamat_soa * gipdot;
extern tamat_soa * aux_ta;
extern thmat_soa * aux_th;

// FERMIONS
extern vec3_soa * ferm_chi_acc; // chi [NPS_tot]
extern vec3_soa * ferm_phi_acc; // phi [NPS_tot]
extern vec3_soa * ferm_out_acc; // auxiliary [NPS_tot]
extern vec3_soa * ferm_shiftmulti_acc; // multishift inverter auxiliary [maxNeededShifts]
extern vec3_soa * kloc_r;  // auxiliary vector
extern vec3_soa * kloc_h;  // auxiliary vector
extern vec3_soa * kloc_s;  // auxiliary vector
extern vec3_soa * kloc_p;  // auxiliary vector
extern vec3_soa * k_p_shiftferm; // auxiliary [maxApproxOrder]

extern vec3_soa * aux1; // used in fermion force calculation, 
                        // for single precision acceleration


// LOCAL SUMS
extern dcomplex_soa * local_sums;
extern double_soa * d_local_sums;

void mem_alloc_core();
void mem_alloc_extended();

void mem_free_core();
void mem_free_extended();

#endif
