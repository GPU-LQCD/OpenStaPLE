#ifndef MAGNETIC_SUSCEPTIBILITY_UTILITIES_C_
#define MAGNETIC_SUSCEPTIBILITY_UTILITIES_C_

#include "./magnetic_susceptibility_utilities.h"
#include "../OpenAcc/struct_c_def.h"
#include "../OpenAcc/backfield.h"
#include "../DbgTools/dbgtools.h"
#include "../Mpi/multidev.h"
#include "../Include/debug.h"

#define acc_twopi 2*3.14159265358979323846

void idphase_dbz(double_soa * idphi_dbz_re, double_soa * idphi_dbz_im, // dphi_dbz_re = 0 
								 ferm_param *fpar,int nnp_openacc[sizeh][4][2]){

    bf_param zero_point = (bf_param){0,0,0,0,0,0};
    bf_param bz1_point  = zero_point;
    bz1_point.bz += 1;

    // We exploit linearity of the UNBOUNDED phases in the field parameters.
    // BOUNDING to -pi +pi BREAKS LINEARITY.

    // NOTICE: the u1 phases so calculated contain the staggere phases and the 
    // chemical potential phases  [NOT BOUNDED, NOT MULTIPLIED BY 2PI] 
    calc_u1_phases_unb_no2pi(idphi_dbz_im,bz1_point ,0,fpar->ferm_charge);// mu = 0 
    // using idphi_dbz_re as temp
    calc_u1_phases_unb_no2pi(idphi_dbz_re,zero_point,0,fpar->ferm_charge);// mu = 0 
   
    // the difference removes the staggered phases and the chemical 
    // potential phases
    phase_diff_in_place(idphi_dbz_im,idphi_dbz_re); // in place

    mult_u1_phases(idphi_dbz_im, acc_twopi);

    set_double_soa_to_zero(idphi_dbz_re); 

    if(debug_settings.print_bfield_dbginfo){
        char tempname[50];                           
        // phases
        sprintf(tempname,"du_dbz_%s_c%d",fpar->name,fpar->printed_bf_dbg_info);
#ifdef MULTIDEVICE      
        strcat(tempname,devinfo.myrankstr);          
#endif
        print_double_soa(idphi_dbz_im,tempname);   

        // plaquettes
        sprintf(tempname,"du_dz_plz_%s_c%d",fpar->name,fpar->printed_bf_dbg_info);
#ifdef MULTIDEVICE      
        strcat(tempname,devinfo.myrankstr);          
#endif

        print_all_abelian_plaquettes(idphi_dbz_im,tempname,nnp_openacc);
        fpar->printed_bf_dbg_info++; 
    }

}


// FINITE differences in e^(i phase_1) - e^(i phase_0)
void delta_u_delta_bz(double_soa * delta_u_re,double_soa * delta_u_im,
        double_soa *tp1, double_soa *tp2,
        ferm_param *fpar, bf_param ext_field){

    bf_param bz1_point  = ext_field;
    bz1_point.bz += 1;
    
    calc_u1_phases(tp1,bz1_point,0,fpar->ferm_charge);// mu = 0 
    calc_u1_phases(tp2,ext_field,0,fpar->ferm_charge);// mu = 0 

    u1_diff( delta_u_re, delta_u_im, tp1, tp2);
    
}




#endif
