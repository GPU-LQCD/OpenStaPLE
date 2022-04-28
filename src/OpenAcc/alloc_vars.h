#ifndef ALLOC_VARS_H_
#define ALLOC_VARS_H_


#include "../Include/common_defines.h"
#include "./struct_c_def.h"

// used in the dynamical allocation of structures




// the gauge configuration, only for read-write
extern global_su3_soa  * conf_rw;

// ued in debugging/testing
extern global_vec3_soa  * ferm_rw; // a global fermion, only for read-write
extern global_tamat_soa  * tamat_rw; // a global tamat, only for read-write
extern global_thmat_soa  * thmat_rw; // a global thmat, only for read-write
// a global dcomplex_soa, only for read-write
extern global_dcomplex_soa *dcomplex_rw;
// a global double_soa, only for read-write
extern global_double_soa *double_rw; 


extern su3_soa  * conf_acc; // the gauge configuration.
extern su3_soa  * conf_acc_bkp; // the old stored conf that will be recovered 
                                // if the metro test fails.
extern su3_soa  * aux_conf_acc; // auxiliary 
extern su3_soa  * auxbis_conf_acc; // auxiliary 
extern double_soa * u1_back_phases; //Background,staggered,chempot phases
                                    // 8 for each flavour
extern double_soa * mag_obs_re;     // Real part of the 'algebra-prefix'
                                    // of magnetization observable 
                                    // 8 for each flavour

extern double_soa * mag_obs_im;     // Imaginary part of the 'algebra-prefix'
                                    // of magnetization observable 
                                    // 8 for each flavour

extern double_soa * topo_loc; //topological charge auxiliary

extern thmat_soa * momenta;// GAUGE FIELD EVOLUTION
extern int momenta_backupped;
extern thmat_soa * momenta_backup;// GAUGE FIELD EVOLUTION - REVERSIBILITY TEST
extern tamat_soa * ipdot_acc;// GAUGE FIELD EVOLUTION
extern tamat_soa * ipdot_g_old;// for HMC diagnostics
extern tamat_soa * ipdot_f_old;// for HMC diagnostics



extern su3_soa * gconf_as_fermionmatrix; //(only a pointer) conf to use in either cases 
                      // in fermion related computation (with or without stouting)

// STOUTING 
extern su3_soa * gstout_conf_acc_arr; // all stouting steps except the zeroth
extern su3_soa * glocal_staples;
extern tamat_soa * gipdot;
extern tamat_soa * aux_ta; // aggiunta per il calcolo della forza stoutata
extern thmat_soa * aux_th; // aggiunta per il calcolo della forza stoutata

// FERMIONS

extern vec3_soa * ferm_chi_acc; // questo e' il chi [NPS_tot]
extern vec3_soa * ferm_phi_acc; // questo e' il phi [NPS_tot]
extern vec3_soa * ferm_out_acc; // questo e' uno ausiliario [NPS_tot]
extern vec3_soa * ferm_shiftmulti_acc; // ausiliario per l'invertitore multishift [maxNeededShifts]
extern vec3_soa * kloc_r;  // vettore ausiliario
extern vec3_soa * kloc_h;  // vettore ausiliario
extern vec3_soa * kloc_s;  // vettore ausiliario
extern vec3_soa * kloc_p;  // vettore ausiliario
extern vec3_soa * k_p_shiftferm; // ausiliario [maxApproxOrder]

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
