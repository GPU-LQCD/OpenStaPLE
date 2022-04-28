double casuale(void);
void initrand(unsigned long s);
void su2_rand(double *pp);
#define PRINT_DETAILS_INSIDE_UPDATE


#include "openacc.h"
#include "../DbgTools/debug_macros_glvarcheck.c"
#include "../RationalApprox/rationalapprox.c"
#include "./struct_c_def.c"
#include "./alloc_vars.c"
#include "./dbgtools.c"
#include "./fermionic_utilities.c"
#include "./single_types.c"
#include "./su3_utilities.c"
#include "./random_assignement.c"
#include "./fermion_matrix.c"
#include "./inverter_full.c"
#include "./find_min_max.c"
#include "./inverter_multishift_full.c"
#include "./rettangoli.c"
#include "./ipdot_gauge.c"
#include "../Meas/ferm_meas.c"
#include "./homebrew_acos.c"
#include "./stouting.c"
#include "./stouting_deottimizzato.c"
#include "./fermion_force.c"
#include "./md_integrator.c"
#include "./update_versatile.c"


int main(){

  initrand(111);
  fflush(stdout);
  printf("INIZIO DEL PROGRAMMA \n");
  su3_soa  * conf_acc;
  su3_soa  * sigma_primo;
  su3_soa  * sigma;
  int  allocation_check =  posix_memalign((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa));
  int  allocation_check =  posix_memalign((void **)&sigma_primo, ALIGN, 8*sizeof(su3_soa));
  int  allocation_check =  posix_memalign((void **)&sigma, ALIGN, 8*sizeof(su3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di conf_acc \n");
  conf_acc->status = IN_USE;
  printf("Allocazione della configurazione : OK \n");

  generate_Conf_cold(sigma_primo,0.3);// 
  generate_Conf_cold(conf_acc,0.1);
  //  generate_Conf_cold(sigma,0.3);


  // INIT FERM PARAMS AND READ RATIONAL APPROX COEFFS
  init_ferm_params();



  mem_alloc();
  printf("Allocazione della memoria : OK \n");
  initialize_global_variables();
  printf("init vars : OK \n");
  compute_nnp_and_nnm_openacc();
  printf("nn computation : OK \n");
#ifdef BACKFIELD
  init_backfield();
  print_double_soa(u1_back_field_phases,"backfield");
  printf("u1_backfield initialization : OK \n");
#endif

  //////  OPENACC CONTEXT INITIALIZATION    //////////////////////////////////////////////////////
  // NVIDIA GPUs
  acc_device_t my_device_type = acc_device_nvidia;
  // AMD GPUs
  // acc_device_t my_device_type = acc_device_radeon;
  // Intel XeonPhi
  //acc_device_t my_device_type = acc_device_xeonphi;
  // Select device ID
  int dev_index = 0;
  SELECT_INIT_ACC_DEVICE(my_device_type, dev_index);
  printf("Device Selected : OK \n");


  //###################### INIZIALIZZAZIONE DELLA CONFIGURAZIONE #################################
  // cold start
  if(start_opt==0){
    generate_Conf_cold(conf_acc,0.01);
    printf("Cold Gauge Conf Generated : OK \n");
    conf_id_iter=0;
  }
 // start from saved conf
  if(start_opt==1){
    read_su3_soa(conf_acc,"stored_config"); // READS ALSO THE conf_id_iter
    printf("Stored Gauge Conf Read : OK \n");
  }
  //###############################################################################################


#pragma acc data   copy(conf_acc[0:8]) copyin(u1_back_field_phases[0:8]) create(ipdot_acc[0:8]) create(aux_conf_acc[0:8]) create(auxbis_conf_acc[0:8]) create(ferm_chi_acc[0:NPS_tot]) create(ferm_phi_acc[0:NPS_tot])  create(ferm_out_acc[0:NPS_tot]) create(ferm_shiftmulti_acc[0:max_ps*max_approx_order]) create(kloc_r[0:1])  create(kloc_h[0:1])  create(kloc_s[0:1])  create(kloc_p[0:1])  create(k_p_shiftferm[0:max_approx_order]) create(momenta[0:8]) copyin(nnp_openacc) copyin(nnm_openacc) create(local_sums[0:2]) create(d_local_sums[0:1])  copyin(fermions_parameters[0:NDiffFlavs])
  {
#ifdef STOUT_FERMIONS
#pragma acc data create(aux_th[0:8]) create(aux_ta[0:8]) create(gstout_conf_acc_arr[0:(8*STOUT_STEPS)]) create(glocal_staples[0:8]) create(gipdot[0:8]) 
      {
#endif
	
	double plq,plqbis,rect;

    int accettate_therm=0;
    int accettate_metro=0;
    int accettate_therm_old=0;
    int accettate_metro_old=0;
    int id_iter_offset=conf_id_iter;
    //################### THERMALIZATION & METRO    ----   UPDATES ####################//
    for(int id_iter=id_iter_offset;id_iter<(ITERATIONS+id_iter_offset);id_iter++){
               accettate_therm_old = accettate_therm;
               accettate_metro_old = accettate_metro;
               conf_id_iter++;
	       printf("\n#################################################\n");
               printf(  "            GENERATING CONF %d of %d\n",conf_id_iter,ITERATIONS+id_iter_offset);
	       printf(  "#################################################\n\n");
	       //--------- CONF UPDATE ----------------//
	       if(id_iter<therm_ITERATIONS){
		 accettate_therm = UPDATE_SOLOACC_UNOSTEP_VERSATILE(conf_acc,
                 residue_metro,residue_md,id_iter-id_iter_offset,
                 accettate_therm,0);
	       }else{
		 accettate_metro = UPDATE_SOLOACC_UNOSTEP_VERSATILE(conf_acc,residue_metro,residue_md,id_iter-id_iter_offset-accettate_therm,accettate_metro,1);
	       }
#pragma acc update host(conf_acc[0:8])
	       //---------------------------------------//

	       //--------- MISURA ROBA FERMIONICA ----------------//
	       FILE *foutfile = fopen(nome_file_ferm_output,"at");
	       if(!foutfile) foutfile = fopen(nome_file_ferm_output,"wt");
	       if(foutfile){
		 fprintf(foutfile,"%d\t",conf_id_iter);
		 for(int iflv=0;iflv<NDiffFlavs;iflv++) perform_chiral_measures(conf_acc,u1_back_field_phases,&(fermions_parameters[iflv]),residue_metro,foutfile);
		 fprintf(foutfile,"\n");
	       }
	       fclose(foutfile);
	       //-------------------------------------------------// 
	       //--------- MISURA ROBA DI GAUGE ------------------//
	       plq = calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
	       rect = calc_rettangolo_soloopenacc(conf_acc,aux_conf_acc,local_sums);
	       
               FILE *goutfile = fopen(nome_file_gauge_output,"at");
               if(!goutfile) goutfile = fopen(nome_file_gauge_output,"wt");
               if(goutfile){
		 if(id_iter<therm_ITERATIONS){
		   printf("Therm_iter %d   Placchetta= %.18lf \n",conf_id_iter,plq/NSITES/6.0/3.0);
		   printf("Therm_iter %d   Placchetta_bis= %.18lf \n",conf_id_iter,plqbis/NSITES/6.0/3.0);
		   printf("Therm_iter %d   Rettangolo= %.18lf \n",conf_id_iter,rect/NSITES/6.0/3.0/2.0);
		   
		   fprintf(goutfile,"%d\t%d\t",conf_id_iter,accettate_therm-accettate_therm_old);
		   fprintf(goutfile,"%.18lf\t%.18lf\n",plq/NSITES/6.0/3.0,rect/NSITES/6.0/3.0/2.0);
		   
		 }else{
		   printf("Metro_iter %d   Placchetta= %.18lf \n",conf_id_iter,plq/NSITES/6.0/3.0);
		   printf("Metro_iter %d   Placchetta_bis= %.18lf \n",conf_id_iter,plqbis/NSITES/6.0/3.0);
		   printf("Metro_iter %d   Rettangolo= %.18lf \n",conf_id_iter,rect/NSITES/6.0/3.0/2.0);
		   
		   fprintf(goutfile,"%d\t%d\t",conf_id_iter,accettate_metro-accettate_metro_old);
		   fprintf(goutfile,"%.18lf\t%.18lf\n",plq/NSITES/6.0/3.0,rect/NSITES/6.0/3.0/2.0);
		   
		 }
	       }
               fclose(goutfile);
	       //-------------------------------------------------//
	       
	       //--------- SALVA LA CONF SU FILE ------------------//
	       if(conf_id_iter%save_conf_every==0)	print_su3_soa(conf_acc,"stored_config");
	       //-------------------------------------------------//
	       
    }// id_iter loop ends here
    


    //--------- SALVA LA CONF SU FILE ------------------//
    print_su3_soa(conf_acc,"stored_config");
    //-------------------------------------------------//

#ifdef STOUT_FERMIONS
    }
#endif
 
  }// end pragma acc data

  CHECKSTATUS(conf_acc);

  //////  OPENACC CONTEXT CLOSING    //////////////////////////////////////////////////////////////
  SHUTDOWN_ACC_DEVICE(my_device_type);
  /////////////////////////////////////////////////////////////////////////////////////////////////


  free(conf_acc);
  mem_free();

  return 0;
}

