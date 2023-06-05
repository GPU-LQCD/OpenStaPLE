// here macros are defined
#define PRINT_DETAILS_INSIDE_UPDATE
#define ALIGN 128

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
#define __restrict
#endif

#ifndef __GNUC__
#include "openacc.h"
#endif

#ifdef ONE_FILE_COMPILATION
#include "../Include/all_include.h"
#endif

#include "../DbgTools/dbgtools.h"
#include "../DbgTools/debugger_hook.h"
#include "../Include/debug.h"
#include "../Include/fermion_parameters.h"
#include "../Include/montecarlo_parameters.h"
#include "../Include/inverter_tricks.h"
#include "../Include/memory_wrapper.h"
#include "../Include/setting_file_parser.h"
#include "../Include/tell_geom_defines.h"
#include "../Include/rep_info.h"
#include "../Include/acceptances_info.h"
#include "../Meas/ferm_meas.h"
#include "../Meas/gauge_meas.h"
#include "../Meas/polyakov.h"
#include "../Meas/measure_topo.h"
#include "../Mpi/communications.h"
#include "../Mpi/multidev.h"
#include "../Rand/random.h"
#include "../RationalApprox/rationalapprox.h"
#include "./action.h"
#include "./alloc_vars.h"
#include "./backfield_parameters.h"
#include "./deviceinit.h"
#include "./fermion_matrix.h"
#include "./fermionic_utilities.h"
#include "./find_min_max.h"
#include "./float_double_conv.h"
#include "./inverter_full.h"
#include "./inverter_multishift_full.h"
#include "./io.h"
#include "./ipdot_gauge.h"
#include "./md_integrator.h"
#include "./md_parameters.h"
#include "./random_assignement.h"
#include "./rectangles.h"
#include "./sp_alloc_vars.h"
#include "./stouting.h"
#include "./struct_c_def.h"
#include "./su3_measurements.h"
#include "./su3_utilities.h"
#include "./update_versatile.h"
#include "./alloc_settings.h"

#ifdef __GNUC__
#include "sys/time.h"
#endif
#include "./cooling.h"
#include <unistd.h>
#include <mpi.h>
#include "../Include/stringify.h"

#include <errno.h>
#include "./HPT_utilities.h"
#include <time.h>

// double level macro, necessary to stringify
// https://gcc.gnu.org/onlinedocs/cpp/Stringification.html
#define xstr(s) str(s) 
#define str(s) #s

// definitions outside the main.

int conf_id_iter;
int verbosity_lv;

int main(int argc, char* argv[]){
 
  gettimeofday ( &(mc_params.start_time), NULL );
  
  if(0==devinfo.myrank_world){
    printf("****************************************************\n");
		if (argc!=2) printf("          COMPILATION INFO                        \n");
    if (argc==2) printf("          PRE INIT - READING SETTING  FILE          \n");
    if (argc==2) printf("     check which parameter corresponds to what! \n");
    printf("commit: %s\n", xstr(COMMIT_HASH) );
    printf("****************************************************\n");
  }
  
  if (GAUGE_ACTION == 1 ){
    if(0==devinfo.myrank_world) printf("\nCOMPILED WITH TREE-LEVEL SYMANZIK IMPROVED GAUGE ACTION\n\n");
  }
  else{
    if(0==devinfo.myrank_world) printf("COMPILED WITH WILSON GAUGE ACTION\n\n");
  }

#ifdef PAR_TEMP
  if(0==devinfo.myrank_world) printf("COMPILED FOR PARALLEL TEMPERING ON BOUNDARY CONDITIONS\n\n");
#else
  if(0==devinfo.myrank_world) printf("COMPILED WITHOUT PARALLEL TEMPERING (1 REPLICA RUN)\n\n");
#endif

  if(argc!=2){
    if(0==devinfo.myrank_world) print_geom_defines();
    if(0==devinfo.myrank_world) printf("\n\nERROR! Use mpirun -n <num_tasks> %s input_file to execute the code!\n\n", argv[0]);
    exit(EXIT_FAILURE);			
  }
  
  FILE *hmc_acc_file;
  FILE *swap_acc_file;
  FILE *file_label;
 
  srand(time(NULL));
    
  // read input file.
#ifdef MULTIDEVICE
	pre_init_multidev1D(&devinfo);
	gdbhook();
#endif

  int input_file_read_check = set_global_vars_and_fermions_from_input_file(argv[1]);

//TODO: abort replicas communication group
    
#ifdef MULTIDEVICE
  if(input_file_read_check){
    printf("MPI%02d: input file reading failed, Aborting...\n",devinfo.myrank);
    MPI_Abort(MPI_COMM_WORLD,1);
  }else init_multidev1D(&devinfo);
#else
  devinfo.myrank = 0;
  devinfo.nranks = 1;
#endif

	if(input_file_read_check){
		printf("MPI%02d: input file reading failed, aborting...\n",devinfo.myrank);
		exit(1);
	}

	if(0==devinfo.myrank_world) print_geom_defines();
	verbosity_lv = debug_settings.input_vbl;

	if(0==devinfo.myrank_world){
		if(0 != mc_params.JarzynskiMode){
			printf("****************************************************\n");
			printf("                   JARZYNSKI MODE              \n");
			printf("     check which parameter corresponds to what! \n");
			printf("****************************************************\n");
		}
	}

	// initialization
#ifdef PAR_TEMP
	int *all_swap_vector;
  int *acceptance_vector;
  double mean_acceptance;
  int *acceptance_vector_old;
    
  all_swap_vector=malloc(sizeof(int)*rep->replicas_total_number-1);
  acceptance_vector=malloc(sizeof(int)*rep->replicas_total_number-1);
  acceptance_vector_old=malloc(sizeof(int)*rep->replicas_total_number-1);
    
  for(int i=0;i<rep->replicas_total_number-1;i++){
    acceptance_vector[i]=0;
    acceptance_vector_old[i]=0;
    all_swap_vector[i]=0;
  }
#endif

  // just printing headtitles
  if(0==devinfo.myrank_world){
    if(0 != mc_params.JarzynskiMode){
      printf("****************************************************\n");
      printf("                   JARZYNSKI MODE              \n");
      printf("     check which parameter corresponds to what! \n");
      printf("****************************************************\n");


    }else {
      printf("****************************************************\n");
      printf("                    NORMAL MODE                \n");
      printf("****************************************************\n");
    }
    if(debug_settings.do_norandom_test){
      printf("****************************************************\n");
      printf("         WELCOME. This is a NORANDOM test.    \n");
      printf("     MOST things will not be random generated,\n");
      printf("            but read from memory instead.     \n");
      printf("                  CHECK THE CODE!!            \n");
      printf("   ALSO: setting the number of trajectories to 1.\n");
      printf("****************************************************\n");
      mc_params.ntraj = 1;

    }
  }

  if(verbosity_lv > 2) 
    printf("MPI%02d, Input file read and initialized multidev1D...\n",
					 devinfo.myrank);

#ifndef __GNUC__
  // OpenAcc context initialization
  // NVIDIA GPUs
  acc_device_t my_device_type = acc_device_nvidia;
  // AMD GPUs
  // acc_device_t my_device_type = acc_device_radeon;
  // Intel XeonPhi
  // acc_device_t my_device_type = acc_device_xeonphi;
  // Select device ID
  printf("MPI%02d: Selecting device.\n", devinfo.myrank);
#ifdef MULTIDEVICE
  select_init_acc_device(my_device_type, (devinfo.single_dev_choice + devinfo.myrank)%devinfo.proc_per_node);
	// select_init_acc_device(my_device_type, devinfo.myrank%devinfo.proc_per_node);
#else
  select_init_acc_device(my_device_type, devinfo.single_dev_choice);
#endif
  printf("Device Selected : OK \n"); //checking printing.
#endif

  unsigned int myseed_default =  (unsigned int) mc_params.seed; 

#ifdef MULTIDEVICE
  myseed_default =  (unsigned int) (myseed_default + devinfo.myrank) ;
  char myrank_string[6];
  sprintf(myrank_string,".R%d",devinfo.myrank);
  strcat(mc_params.RandGenStatusFilename,myrank_string);
#endif

  initrand_fromfile(mc_params.RandGenStatusFilename,myseed_default);

  // init ferm params and read rational approx coeffs
  if(init_ferm_params(fermions_parameters)){
    printf("MPI%02d - Finalizing...\n",devinfo.myrank); //cp
#ifdef MULTIDEVICE
    MPI_Finalize();
#endif
    exit(1);
  }
	#pragma acc enter data copyin(fermions_parameters[0:alloc_info.NDiffFlavs])
    
  mem_alloc_core();
  mem_alloc_extended();
 
  // single/double precision allocation
  //TODO: extend print to all strings
  printf("\n   MPI%02d:%02d - Memory allocation (double) : OK \n\n\n",devinfo.replica_idx,devinfo.myrank);
  if(inverter_tricks.useMixedPrecision || md_parameters.singlePrecMD){
    mem_alloc_core_f();
    printf("\n  MPI%02d - Memory allocation (float) [CORE]: OK \n\n\n",devinfo.myrank);
  }

  if( md_parameters.singlePrecMD){
    mem_alloc_extended_f();
    printf("\n  MPI%02d - Memory allocation (float) [EXTENDED]: OK \n\n\n",devinfo.myrank);
  }
   
  printf("\n  MPI%02d - Total allocated memory: %zu \n\n\n",devinfo.myrank,max_memory_used);
  
	gl_stout_rho=act_params.stout_rho;
	gl_topo_rho=act_params.topo_rho;
	#pragma acc enter data copyin(gl_stout_rho)
	#pragma acc enter data copyin(gl_topo_rho)
  compute_nnp_and_nnm_openacc();

	#pragma acc enter data copyin(nnp_openacc)
	#pragma acc enter data copyin(nnm_openacc)
  printf("MPI%02d - nn computation : OK\n",devinfo.myrank);
  init_all_u1_phases(backfield_parameters,fermions_parameters);
	#pragma acc update device(u1_back_phases[0:8*alloc_info.NDiffFlavs])
	#pragma acc update device(mag_obs_re[0:8*alloc_info.NDiffFlavs])
	#pragma acc update device(mag_obs_im[0:8*alloc_info.NDiffFlavs])

  if(inverter_tricks.useMixedPrecision || md_parameters.singlePrecMD){
		#pragma acc update device(u1_back_phases_f[0:8*alloc_info.NDiffFlavs])
  }

  printf("MPI%02d - u1_backfield initialization (float & double): OK \n",devinfo.myrank);

  initialize_md_global_variables(md_parameters);
  printf("MPI%02d - init md vars : OK \n",devinfo.myrank);

	
#ifdef PAR_TEMP
	defect_info def;
  char rep_str [20];
  char aux_name_file[200];
  strcpy(aux_name_file,mc_params.save_conf_name);
#endif
	
  for(int r=0;r<rep->replicas_total_number;r++){
#ifdef PAR_TEMP
    snprintf(rep_str,20,"replica_%d",r);
		strcat(mc_params.save_conf_name,rep_str);
#endif
		
    if(debug_settings.do_norandom_test){
      if(!read_conf_wrapper(conf_acc[r],"conf_norndtest",&conf_id_iter,debug_settings.use_ildg)){
				printf("MPI%02d - Stored Gauge Conf conf_norndtest Read : OK\n",devinfo.myrank);
      }
      else{
				// cold start
				printf("MPI%02d - COMPILED IN NORANDOM MODE. A CONFIGURATION FILE NAMED\
                    \"conf_norndtest\" MUST BE PRESENT\n",devinfo.myrank);
				exit(1);
      }
    }
    else{
      if(!read_conf_wrapper(conf_acc[r],mc_params.save_conf_name,
														&conf_id_iter,debug_settings.use_ildg)){
				printf("MPI%02d - Stored Gauge Conf \"%s\" Read : OK \n",
							 devinfo.myrank, mc_params.save_conf_name);
      }
      else{
            
				generate_Conf_cold(conf_acc[r],mc_params.eps_gen);
				printf("MPI%02d - Cold Gauge Conf Generated : OK \n",
							 devinfo.myrank);
				conf_id_iter=0;
      }
    }
#ifdef PAR_TEMP
    strcpy(mc_params.save_conf_name,aux_name_file);

    //TODO: check consistency with rep->replicas_total_number and nranks_world/NRANKS_D3
    
		if (0==devinfo.myrank_world) printf("%d/%d Defect initialization\n",r,rep->replicas_total_number); 
		rep->label[r]=r;
		init_k(conf_acc[r],rep->cr_vec[r],rep->defect_boundary,rep->defect_coordinates,&def,r);
		
#ifdef MULTIDEVICE
    if(devinfo.async_comm_gauge) init_k(&conf_acc[r][8],rep->cr_vec[r],rep->defect_boundary,rep->defect_coordinates,&def,1);
#endif
		
		#pragma acc update device(conf_acc[r:1][0:alloc_info.conf_acc_size])
		if(md_parameters.singlePrecMD){
			convert_double_to_float_su3_soa(conf_acc[r],conf_acc_f[r]);
			//^^ NOTE: doing this because a K initialization for su3_soa_f doesn't exist. Please create it.
#ifdef MULTIDEVICE
			if(devinfo.async_comm_gauge){
				convert_double_to_float_su3_soa(&conf_acc[r][8],&conf_acc_f[r][8]);
				//^^ NOTE: doing this because a K initialization for su3_soa_f doesn't exist. Please create it.
			}
#endif
			
			#pragma acc update host(conf_acc_f[r:1][0:alloc_info.conf_acc_size])
		}
#else
		#pragma acc update device(conf_acc[r:1][0:alloc_info.conf_acc_size])
#endif
  }

#ifdef PAR_TEMP

  int vec_aux_bound[3]={1,1,1};

	if (0==devinfo.myrank_world) printf("Auxiliary confs defect initialization\n");
  init_k(aux_conf_acc,1,0,vec_aux_bound,&def,1);
  init_k(auxbis_conf_acc,1,0,vec_aux_bound,&def,1);
	#pragma acc update device(aux_conf_acc[0:8])
	#pragma acc update device(auxbis_conf_acc[0:8])

	if(md_parameters.singlePrecMD){
		convert_double_to_float_su3_soa(aux_conf_acc,aux_conf_acc_f);
		convert_double_to_float_su3_soa(auxbis_conf_acc,auxbis_conf_acc_f);
		#pragma acc update host(aux_conf_acc_f[0:8])
		#pragma acc update host(auxbis_conf_acc_f[0:8])
	}

	if(alloc_info.stoutAllocations){
		int stout_steps = ((act_params.topo_stout_steps>act_params.stout_steps) & (act_params.topo_action==1)?
											  act_params.topo_stout_steps:act_params.stout_steps );
		for (int i = 0; i < stout_steps; i++)
			init_k(&gstout_conf_acc_arr[8*i],1,0,vec_aux_bound,&def,1);
		#pragma acc update device(gstout_conf_acc_arr[0:8*stout_steps])
		if(md_parameters.singlePrecMD){
			for (int i = 0; i < stout_steps; i++)
				convert_double_to_float_su3_soa(&gstout_conf_acc_arr[8*i],&gstout_conf_acc_arr_f[8*i]);
			#pragma acc update host(gstout_conf_acc_arr_f[0:8*stout_steps])
		}
	}
#endif
	
  double max_unitarity_deviation,avg_unitarity_deviation;
    
  for(int r=0;r<rep->replicas_total_number;r++){
    check_unitarity_host(conf_acc[r],&max_unitarity_deviation,&avg_unitarity_deviation);
    printf("\tMPI%02d: Avg_unitarity_deviation on host: %e\n", devinfo.myrank, 
					 avg_unitarity_deviation);
    printf("\tMPI%02d: Max_unitarity_deviation on host: %e\n", devinfo.myrank,
					 max_unitarity_deviation);
  }

	
  // measures
    
  double plq,rect;
  double cool_topo_ch[meastopo_params.coolmeasstep/meastopo_params.cool_measinterval+1];
  double stout_topo_ch[meastopo_params.stoutmeasstep/meastopo_params.stout_measinterval+1];
  d_complex poly;

  int *accettate_therm;
  int *accettate_metro;
  int *iterations;
    
  int *accettate_therm_old;
  int *accettate_metro_old;
  int id_iter_offset=conf_id_iter;

  accettate_therm=malloc(sizeof(int)*rep->replicas_total_number);
  accettate_metro=malloc(sizeof(int)*rep->replicas_total_number);
  accettate_therm_old=malloc(sizeof(int)*rep->replicas_total_number);
  accettate_metro_old=malloc(sizeof(int)*rep->replicas_total_number);
  iterations=malloc(sizeof(int)*rep->replicas_total_number);
    
  // inizialization to 0
  for(int i=0;i<rep->replicas_total_number;i++){
    accettate_metro[i]=0;
    accettate_therm[i]=0;
    iterations[i]=0;
    accettate_therm_old[i]=0;
    accettate_metro_old[i]=0;
  }
	int swap_number=0;

  // plaquette measures and polyakov loop measures.
  printf("PLAQUETTE START\n");
    
  plq = calc_plaquette_soloopenacc(conf_acc[0],aux_conf_acc,local_sums);
  printf("\tMPI%02d: Therm_iter %d Placchetta    = %.18lf \n",
				 devinfo.myrank, conf_id_iter,plq/GL_SIZE/6.0/3.0);
    
  printf("PLAQUETTE END\n");

#if !defined(GAUGE_ACT_WILSON) || !defined(MULTIDEVICE)
  rect = calc_rettangolo_soloopenacc(conf_acc[0],aux_conf_acc,local_sums);
  printf("\tMPI%02d: Therm_iter %d Rettangolo    = %.18lf \n",
				 devinfo.myrank, conf_id_iter,rect/GL_SIZE/6.0/3.0/2.0);
#else
  printf("\tMPI%02d: multidevice rectangle computation with Wilson action not implemented\n",devinfo.myrank);
#endif

  poly =  (*polyakov_loop[geom_par.tmap])(conf_acc[0]);
  printf("\tMPI%02d: Therm_iter %d Polyakov Loop = (%.18lf, %.18lf)  \n",
				 devinfo.myrank, conf_id_iter,creal(poly),cimag(poly));

	// conf swap debug routine

	/*printf("SWAP_DEBOIAG_GREPPA_QUI\n");
		for(int i=0; i<(rep->replicas_total_number-1); i++) {
		int rep_indx1=i;
		int rep_indx2=i+1;
		double delta_S_STUPID=0.0;

		// action old
		double S_1 = - C_ZERO * BETA_BY_THREE * calc_plaquette_soloopenacc(conf_acc[rep_indx1], aux_conf_acc, local_sums);
		double S_2 = - C_ZERO * BETA_BY_THREE * calc_plaquette_soloopenacc(conf_acc[rep_indx2], aux_conf_acc, local_sums);
		#ifdef GAUGE_ACT_TLSM
		S_1 += - C_ONE * BETA_BY_THREE * calc_rettangolo_soloopenacc(conf_acc[rep_indx1], aux_conf_acc, local_sums);
		S_2 += - C_ONE * BETA_BY_THREE * calc_rettangolo_soloopenacc(conf_acc[rep_indx2], aux_conf_acc, local_sums);
		#endif
		delta_S_STUPID -= (S_1 + S_2);  

		// swap conf
		replicas_swap(conf_acc[rep_indx1], conf_acc[rep_indx2], rep_indx1, rep_indx2, rep);
		#pragma acc update device(conf_acc[0:rep->replicas_total_number][0:8])

		// action new
		S_1 = - C_ZERO * BETA_BY_THREE * calc_plaquette_soloopenacc(conf_acc[rep_indx1], aux_conf_acc, local_sums);
		S_2 = - C_ZERO * BETA_BY_THREE * calc_plaquette_soloopenacc(conf_acc[rep_indx2], aux_conf_acc, local_sums);
		#ifdef GAUGE_ACT_TLSM
		S_1 += - C_ONE * BETA_BY_THREE * calc_rettangolo_soloopenacc(conf_acc[rep_indx1], aux_conf_acc, local_sums);
		S_2 += - C_ONE * BETA_BY_THREE * calc_rettangolo_soloopenacc(conf_acc[rep_indx2], aux_conf_acc, local_sums);
		#endif
		delta_S_STUPID += (S_1 + S_2);

		// go back
		replicas_swap(conf_acc[rep_indx1], conf_acc[rep_indx2], rep_indx1, rep_indx2, rep);
		#pragma acc update device(conf_acc[0:rep->replicas_total_number][0:8])

		// compute delta_S_SWAP smart
		double delta_S_SWAP = calc_Delta_S_soloopenacc_SWAP(conf_acc[rep_indx1], conf_acc[rep_indx2], aux_conf_acc, local_sums, &def);

		// compare
		if(0==devinfo.myrank) printf("(%d<->%d) %.15lg %.15lg\n", rep_indx1, rep_indx1+1, delta_S_SWAP, delta_S_STUPID);
		}

		// shut down program
		mem_free_core();
		mem_free_extended();
		#ifdef MULTIDEVICE
		shutdown_multidev();
		#endif
		exit(1);*/

	// end debug

  //Here we are in Jarzynski mode.
    
  if(0 == mc_params.ntraj && 0 == mc_params.JarzynskiMode ){ // measures only
      
    printf("\n#################################################\n");
    printf("\tMEASUREMENTS ONLY ON FILE %s\n", mc_params.save_conf_name);
    printf("\n#################################################\n");

    // gauge stuff measures

		if(0 == devinfo.myrank ) printf("Gauge Measures:\n");
        
    plq = calc_plaquette_soloopenacc(conf_acc[0],aux_conf_acc,local_sums);
#if !defined(GAUGE_ACT_WILSON) || !defined(MULTIDEVICE)
    rect = calc_rettangolo_soloopenacc(conf_acc[0],aux_conf_acc,local_sums);
#else
    printf("\tMPI%02d: multidevice rectangle computation with Wilson action not implemented\n",devinfo.myrank);
#endif 
    poly =  (*polyakov_loop[geom_par.tmap])(conf_acc[0]);//misura polyakov loop
        
    printf("Plaquette     : %.18lf\n" ,plq/GL_SIZE/3.0/6.0);
    printf("Rectangle     : %.18lf\n" ,rect/GL_SIZE/3.0/6.0/2.0);
    printf("Polyakov Loop : (%.18lf,%.18lf) \n",creal(poly),cimag(poly));

    // fermionic stuff measures

    if(0 == devinfo.myrank)  printf("Fermion Measurements: see file %s\n",
																		fm_par.fermionic_outfilename);
    fermion_measures(conf_acc[0],fermions_parameters,
										 &fm_par, md_parameters.residue_metro, 
										 md_parameters.max_cg_iterations, id_iter_offset,
										 plq/GL_SIZE/3.0/6.0,
										 rect/GL_SIZE/3.0/6.0/2.0);   

  }else printf("MPI%02d: Starting generation of Configurations.\n",
							 devinfo.myrank);
     
  // thermalization & metropolis updates

  int id_iter=id_iter_offset;
    
  init_global_program_status(); 

  printf("run_condition: %d\n",mc_params.run_condition) ;
  if ( 0 != mc_params.ntraj ) {
    while ( RUN_CONDITION_TERMINATE != mc_params.run_condition)
      {
        if(GPSTATUS_UPDATE == mc_params.next_gps){
					struct timeval tstart_cycle,tend_cycle;
					gettimeofday(&tstart_cycle, NULL);

					if(0 != mc_params.JarzynskiMode ){

						bf_param new_backfield_parameters = backfield_parameters;

						// direct mode 
						if(1 == mc_params.JarzynskiMode)
							new_backfield_parameters.bz = backfield_parameters.bz + 
								(double) id_iter/mc_params.MaxConfIdIter;
						// reverse mode
						if(-1 == mc_params.JarzynskiMode)
							new_backfield_parameters.bz = backfield_parameters.bz -
								(double) id_iter/mc_params.MaxConfIdIter;

						if(0==devinfo.myrank_world){

							if(1 == mc_params.JarzynskiMode)
								printf("\n\nJarzynskiMode - DIRECT - From bz=%f to bz=%f+1 in %d steps.\n",
											 backfield_parameters.bz , backfield_parameters.bz, 
											 mc_params.MaxConfIdIter);
							if(-1 == mc_params.JarzynskiMode)
								printf("\n\nJarzynskiMode - REVERSE - From bz=%f to bz=%f-1 in %d steps.\n",
											 backfield_parameters.bz , backfield_parameters.bz, 
											 mc_params.MaxConfIdIter);

							if(1 == mc_params.JarzynskiMode)
								printf("\n\nJarzynskiMode - DIRECT - From bz=%f to bz=%f+1 in %d steps.\n",
											 backfield_parameters.bz , backfield_parameters.bz, 
											 mc_params.MaxConfIdIter);
							if(-1 == mc_params.JarzynskiMode)
								printf("\n\nJarzynskiMode - REVERSE - From bz=%f to bz=%f-1 in %d steps.\n",
											 backfield_parameters.bz , backfield_parameters.bz, 
											 mc_params.MaxConfIdIter);

							printf("JarzynskiMode, iteration %d/%d (%d max for this run)\n",
										 id_iter,mc_params.MaxConfIdIter,mc_params.ntraj);
							printf("JarzynskiMode - current bz value : %f\n", new_backfield_parameters.bz);
						}

						init_all_u1_phases(new_backfield_parameters,fermions_parameters);
						#pragma acc update device(u1_back_phases[0:8*alloc_info.NDiffFlavs])
						#pragma acc update device(u1_back_phases_f[0:8*alloc_info.NDiffFlavs])
						printf("Jarzynski mode's end\n");
					}


					check_unitarity_device(conf_acc[0],&max_unitarity_deviation,
																 &avg_unitarity_deviation);
					printf("\tMPI%02d: Avg/Max unitarity deviation on device: %e / %e\n", 
								 devinfo.myrank,avg_unitarity_deviation,max_unitarity_deviation);
            
					for (int i=0;i<rep->replicas_total_number;i++){
            accettate_therm_old [i]= accettate_therm[i];
            accettate_metro_old [i]= accettate_metro[i];
					}

					if(devinfo.myrank ==0 ){
						printf("\n#################################################\n");
						printf(  "   GENERATING CONF %d of %d, %dx%dx%dx%d,%1.3f \n",
										 conf_id_iter,mc_params.ntraj+id_iter_offset,
										 geom_par.gnx,geom_par.gny,
										 geom_par.gnz,geom_par.gnt,
										 act_params.beta);
						printf(  "#################################################\n\n");
					}

#ifdef PAR_TEMP
					for(int i=0;i<rep->replicas_total_number-1;i++){
						acceptance_vector_old[i]=acceptance_vector[i];
					}
#endif
		
					// replicas update - hpt step
        
					for(int r=0;r<rep->replicas_total_number;r++){
						printf("REPLICA %d:\n",r);

						// initial action
			
						if (verbosity_lv>10){
							double action;
							action  = - C_ZERO * BETA_BY_THREE * calc_plaquette_soloopenacc(conf_acc[r], aux_conf_acc, local_sums);
#ifdef GAUGE_ACT_TLSM
							action += - C_ONE  * BETA_BY_THREE * calc_rettangolo_soloopenacc(conf_acc[r], aux_conf_acc, local_sums);
#endif
							printf("ACTION BEFORE HMC STEP REPLICA %d: %.15lg\n", r, action);
						}

						// HMC step
			
						if(id_iter<mc_params.therm_ntraj){
							accettate_therm[r] = UPDATE_SOLOACC_UNOSTEP_VERSATILE(conf_acc[r],
																																		md_parameters.residue_metro,md_parameters.residue_md,
																																		id_iter-id_iter_offset,
																																		accettate_therm[r],0,md_parameters.max_cg_iterations);
						}
						else{
							accettate_metro[r] = UPDATE_SOLOACC_UNOSTEP_VERSATILE(conf_acc[r],
																																		md_parameters.residue_metro,md_parameters.residue_md,
																																		id_iter-id_iter_offset-accettate_therm[r],
																																		accettate_metro[r],1,md_parameters.max_cg_iterations);
							if(0==devinfo.myrank_world){
								iterations[r] = id_iter-id_iter_offset-accettate_therm[r]+1;
								double acceptance = (double) accettate_metro[r] / iterations[r];
								double acc_err = sqrt((double)accettate_metro[r]*(iterations[r]-accettate_metro[r])/iterations[r])/iterations[r];
								printf("Estimated HMC acceptance for this run [replica %d]: %f +- %f\n. Iterations: %d",r,acceptance, acc_err, iterations[r]);
							}
						}
						#pragma acc update host(conf_acc[0:rep->replicas_total_number][0:8])

						// final action

						if (verbosity_lv>10){
							double action;
							action  = - C_ZERO * BETA_BY_THREE * calc_plaquette_soloopenacc(conf_acc[r], aux_conf_acc, local_sums);
#ifdef GAUGE_ACT_TLSM
							action += - C_ONE  * BETA_BY_THREE * calc_rettangolo_soloopenacc(conf_acc[r], aux_conf_acc, local_sums);
#endif
							printf("ACTION AFTER HMC STEP REPLICA %d: %.15lg\n", r, action);
						}

#ifdef PAR_TEMP
						if(rep->replicas_total_number>1){
							// conf swap
							if (0==devinfo.myrank_world) {printf("CONF SWAP PROPOSED\n");}
							All_Conf_SWAP(conf_acc,aux_conf_acc,local_sums, &def, &swap_number,all_swap_vector,acceptance_vector, rep);
							if (0==devinfo.myrank_world) {printf("Number of accepted swaps: %d\n", swap_number);}       
							#pragma acc update host(conf_acc[0:rep->replicas_total_number][0:8])
                
							// periodic conf translation
							trasl_conf(conf_acc[0],auxbis_conf_acc);
						}
						#pragma acc update host(conf_acc[0:rep->replicas_total_number][0:8])
#endif
					}

					id_iter++;
					conf_id_iter++;
            
#ifdef PAR_TEMP
					printf("MPI%02d - Printing acceptances - only by master rank...\n",
								 devinfo.myrank);
					if(devinfo.myrank ==0){
                
						if(rep->replicas_total_number>1){
							file_label=fopen(acc_info->file_label_name,"at");
							if(!file_label){file_label=fopen(acc_info->file_label_name,"wt");}

							hmc_acc_file=fopen(acc_info->hmc_file_name,"at");
							if(!hmc_acc_file){hmc_acc_file=fopen(acc_info->hmc_file_name,"wt");}
                    
							swap_acc_file=fopen(acc_info->swap_file_name,"at");
							if(!swap_acc_file){swap_acc_file=fopen(acc_info->swap_file_name,"wt");}

							fprintf(hmc_acc_file,"%d\t",conf_id_iter);
							fprintf(swap_acc_file,"%d\t",conf_id_iter);
							label_print(rep, file_label, conf_id_iter);
						}
            // print acceptances
						for(int i=0;i<rep->replicas_total_number;i++){
							if(i<rep->replicas_total_number-1){
								mean_acceptance=(double)acceptance_vector[i]/all_swap_vector[i];
								if(0==devinfo.myrank_world){ printf("replica couple [%d/%d]: proposed %d, accepted %d, mean_acceptance %f\n",i,i+1,all_swap_vector[i],acceptance_vector[i],mean_acceptance);}
								if(rep->replicas_total_number>1){
									fprintf(swap_acc_file,"%d\t",acceptance_vector[i]-acceptance_vector_old[i]);}
							}
            
							if(rep->replicas_total_number>1){
                fprintf(hmc_acc_file,"%d\t", accettate_therm[i]+accettate_metro[i]
                        -accettate_therm_old[i]-accettate_metro_old[i]);
							}
              
            }
            
						if(rep->replicas_total_number>1){
            
							fprintf(hmc_acc_file,"\n");
							fprintf(swap_acc_file,"\n");
            
							fclose(hmc_acc_file);
							fclose(swap_acc_file);
							fclose(file_label);
            
						}
        
					}
#endif
          
					// gauge stuff measures
					printf("===========GAUGE MEASURING============\n");
            
					plq  = calc_plaquette_soloopenacc(conf_acc[0],aux_conf_acc,local_sums);
#if !defined(GAUGE_ACT_WILSON) || !defined(MULTIDEVICE)
					rect = calc_rettangolo_soloopenacc(conf_acc[0],aux_conf_acc,local_sums);
#else
					printf("\tMPI%02d: multidevice rectangle computation with Wilson action not implemented\n",devinfo.myrank);
#endif
					poly =  (*polyakov_loop[geom_par.tmap])(conf_acc[0]);
            
					if(meastopo_params.meascool && conf_id_iter%meastopo_params.cooleach==0){
						su3_soa *conf_to_use;
						cool_topo_ch[0]=compute_topological_charge(conf_acc[0],auxbis_conf_acc,topo_loc);
						for(int cs = 1; cs <= meastopo_params.coolmeasstep; cs++){
							if(cs==1)
								conf_to_use=(su3_soa*)conf_acc[0];
							else
								conf_to_use=(su3_soa*)aux_conf_acc;
							cool_conf(conf_to_use,aux_conf_acc,auxbis_conf_acc);
							if(cs%meastopo_params.cool_measinterval==0)
								cool_topo_ch[cs/meastopo_params.cool_measinterval]=compute_topological_charge(aux_conf_acc,auxbis_conf_acc,topo_loc);
						}
						printf("MPI%02d - Printing cooled charge - only by master rank...\n",
									 devinfo.myrank);
						if(devinfo.myrank ==0){
							FILE *cooloutfile = fopen(meastopo_params.pathcool,"at");
							if(!cooloutfile){
								cooloutfile = fopen(meastopo_params.pathcool,"wt");
								char coolheader[35];
								strcpy(coolheader,"#conf_id\tCoolStp\tTopoChCool\n");
								fprintf(cooloutfile,"%s",coolheader);
							}
							if(cooloutfile){
								for(int i = 0; i <= meastopo_params.coolmeasstep/meastopo_params.cool_measinterval;i++)
									fprintf(cooloutfile,"%d\t%d\t%18.18lf\n",conf_id_iter,
													i*meastopo_params.cool_measinterval,
													cool_topo_ch[i]);
							}
							fclose(cooloutfile);
						}
					}

					if(meastopo_params.measstout && conf_id_iter%meastopo_params.stouteach==0){
						stout_wrapper(conf_acc[0],gstout_conf_acc_arr,1);
						stout_topo_ch[0]=compute_topological_charge(conf_acc[0],auxbis_conf_acc,topo_loc);
						for(int ss = 0; ss < meastopo_params.stoutmeasstep; ss+=meastopo_params.stout_measinterval){
							int topoindx =1+ss/meastopo_params.stout_measinterval; 
							stout_topo_ch[topoindx]=compute_topological_charge(&gstout_conf_acc_arr[8*ss],auxbis_conf_acc,topo_loc);
						}
	    
						printf("MPI%02d - Printing stouted charge - only by master rank...\n",
									 devinfo.myrank);
						if(devinfo.myrank ==0){
							FILE *stoutoutfile = fopen(meastopo_params.pathstout,"at");
							if(!stoutoutfile){
								stoutoutfile = fopen(meastopo_params.pathstout,"wt");
								char stoutheader[35];
								strcpy(stoutheader,"#conf_id\tStoutStp\tTopoChStout\n");
								fprintf(stoutoutfile,"%s",stoutheader);
							}
							if(stoutoutfile){
								for(int i = 0; i <= meastopo_params.stoutmeasstep/meastopo_params.stout_measinterval;i++)
									fprintf(stoutoutfile,"%d\t%d\t%18.18lf\n",conf_id_iter,
													i*meastopo_params.stout_measinterval,
													stout_topo_ch[i]);
							}
							fclose(stoutoutfile);
						}
					}//if stout end

					printf("MPI%02d - Printing gauge obs - only by master rank...\n",
								 devinfo.myrank);
					if(devinfo.myrank ==0){
						FILE *goutfile = fopen(gauge_outfilename,"at");
						if(!goutfile){
							goutfile = fopen(gauge_outfilename,"wt");
							strcpy(gauge_outfile_header,"#conf_id\tacc\tplq\trect\tReP\tImP\n");
							fprintf(goutfile,"%s",gauge_outfile_header);
						}
						if(goutfile){
							if(id_iter<mc_params.therm_ntraj){
								printf("Therm_iter %d",conf_id_iter );
								printf("Plaquette = %.18lf    ", plq/GL_SIZE/6.0/3.0);
								printf("Rectangle = %.18lf\n",rect/GL_SIZE/6.0/3.0/2.0);
							}else printf("Metro_iter %d   Plaquette= %.18lf    Rectangle = %.18lf\n",conf_id_iter,plq/GL_SIZE/6.0/3.0,rect/GL_SIZE/6.0/3.0/2.0);

							fprintf(goutfile,"%d\t%d\t",conf_id_iter,
											accettate_therm[0]+accettate_metro[0]
											-accettate_therm_old[0]-accettate_metro_old[0]);
                    
							fprintf(goutfile,"%.18lf\t%.18lf\t%.18lf\t%.18lf\n",
											plq/GL_SIZE/6.0/3.0,
											rect/GL_SIZE/6.0/3.0/2.0, 
											creal(poly), cimag(poly));
						}
						fclose(goutfile);
					}
       
					// saving conf_store

					// saves gauge conf and rng status to file
					if(conf_id_iter%mc_params.storeconfinterval==0){ 
						char tempname[50];
						char serial[10];
						strcpy(tempname,mc_params.store_conf_name);
						sprintf(serial,".%05d",conf_id_iter);
						strcat(tempname,serial);
						printf("MPI%02d - Storing conf %s.\n",
									 devinfo.myrank, tempname);
						save_conf_wrapper(conf_acc[0],tempname,conf_id_iter,
															debug_settings.use_ildg);
						strcpy(tempname,mc_params.RandGenStatusFilename);
						sprintf(serial,".%05d",conf_id_iter);
						strcat(tempname,serial);
						printf("MPI%02d - Storing rng status in %s.\n", 
									 devinfo.myrank , tempname);
						saverand_tofile(tempname);
					}
  
					if(conf_id_iter%mc_params.saveconfinterval==0){
						for(int r=0;r<rep->replicas_total_number;r++){
#ifdef PAR_TEMP        
							snprintf(rep_str,20,"replica_%d",r);
							strcat(mc_params.save_conf_name,rep_str);
#endif
							if (debug_settings.SaveAllAtEnd){
								printf("MPI%02d - Saving conf %s.\n", devinfo.myrank,
											 mc_params.save_conf_name);
								save_conf_wrapper(conf_acc[r],mc_params.save_conf_name, conf_id_iter,
																	debug_settings.use_ildg);
							}else
								printf("\n\nMPI%02d: WARNING, \'SaveAllAtEnd\'=0,NOT SAVING/OVERWRITING CONF AND RNG STATUS.\n\n\n", devinfo.myrank);
#ifdef PAR_TEMP        
							strcpy(mc_params.save_conf_name,aux_name_file);
#endif
						}
						if (debug_settings.SaveAllAtEnd){
							printf("MPI%02d - Saving rng status in %s.\n", devinfo.myrank, 
										 mc_params.RandGenStatusFilename);
							saverand_tofile(mc_params.RandGenStatusFilename);
						}
					}

					gettimeofday(&tend_cycle, NULL);

					double update_time = (double) 
						(tend_cycle.tv_sec - tstart_cycle.tv_sec)+
						(double)(tend_cycle.tv_usec - tstart_cycle.tv_usec)/1.0e6;
            
					mc_params.max_update_time = (update_time > mc_params.max_update_time)?
						update_time :mc_params.max_update_time;


					if(0==devinfo.myrank_world){
						printf("Tot time : %f sec (with measurements)\n", update_time);
						if(debug_settings.save_diagnostics == 1){
							FILE *foutfile = fopen(debug_settings.diagnostics_filename,"at");
							fprintf(foutfile,"TOTTIME  %f \n",update_time);
							fclose(foutfile);
						}

					}
  
        }


        if (GPSTATUS_FERMION_MEASURES == mc_params.next_gps){
    
            
					// fermionic stuff measures
            
					if(0 != mc_params.JarzynskiMode ){ // halfway measurements for Jarzynski

						bf_param new_backfield_parameters = backfield_parameters;

						// direct mode 
						if(1 == mc_params.JarzynskiMode)
							new_backfield_parameters.bz = backfield_parameters.bz + 
								(double) (id_iter+0.5)/mc_params.MaxConfIdIter;
						// reverse mode
						if(-1 == mc_params.JarzynskiMode)
							new_backfield_parameters.bz = backfield_parameters.bz -
								(double) (id_iter+0.5)/mc_params.MaxConfIdIter;


						if(0==devinfo.myrank_world){

							printf("JarzynskiMode, iteration %d/%d (%d max for this run) - MEASUREMENTS AT HALFWAY \n",
										 id_iter,mc_params.MaxConfIdIter,mc_params.ntraj);
							printf("JarzynskiMode - current bz value : %f (HALFWAY)\n", new_backfield_parameters.bz);
						}

						init_all_u1_phases(new_backfield_parameters,fermions_parameters);
						#pragma acc update device(u1_back_phases[0:8*alloc_info.NDiffFlavs])
						#pragma acc update device(u1_back_phases_f[0:8*alloc_info.NDiffFlavs])

					}

					check_unitarity_device(conf_acc[0],&max_unitarity_deviation,
																 &avg_unitarity_deviation);
					printf("\tMPI%02d: Avg/Max unitarity deviation on device: %e / %e\n", 
								 devinfo.myrank,avg_unitarity_deviation,max_unitarity_deviation);

					if(conf_id_iter % fm_par.measEvery == 0 )
						mc_params.next_gps = GPSTATUS_FERMION_MEASURES;

					struct timeval tf0, tf1;
					gettimeofday(&tf0, NULL);
					fermion_measures(conf_acc[0],fermions_parameters,
													 &fm_par, md_parameters.residue_metro,
													 md_parameters.max_cg_iterations,conf_id_iter,
													 plq/GL_SIZE/3.0/6.0,
													 rect/GL_SIZE/3.0/6.0/2.0);   

					gettimeofday(&tf1, NULL);

					double fermionMeasureTiming =
						(double) (tf1.tv_sec - tf0.tv_sec)+
						(double)(tf1.tv_usec - tf0.tv_usec)/1.0e6;

					if(debug_settings.save_diagnostics == 1){
						FILE *foutfile = 
							fopen(debug_settings.diagnostics_filename,"at");

						if(conf_id_iter % fm_par.measEvery == 0 )
							fprintf(foutfile,"FERMMEASTIME  %f \n",fermionMeasureTiming);
						fclose(foutfile);
					}
                  
					// save RNG status
					if(conf_id_iter%mc_params.storeconfinterval==0){
						char tempname[50];
						char serial[10];
						strcpy(tempname,mc_params.RandGenStatusFilename);
						sprintf(serial,".%05d",conf_id_iter);
						strcat(tempname,serial);
						printf("MPI%02d - Storing rng status in %s.\n", 
									 devinfo.myrank , tempname);
						saverand_tofile(tempname);
					} 

					if(conf_id_iter%mc_params.saveconfinterval==0){
						if( debug_settings.SaveAllAtEnd){
							printf("MPI%02d - Saving rng status in %s.\n", devinfo.myrank, 
										 mc_params.RandGenStatusFilename);
							saverand_tofile(mc_params.RandGenStatusFilename);
						}
						else printf(
												"\n\nMPI%02d: WARNING, \'SaveAllAtEnd\'=0,NOT SAVING/OVERWRITING RNG STATUS.\n\n\n", devinfo.myrank);
					}

        }
        // determining next thing to do
        if(0 == conf_id_iter % fm_par.measEvery)
					mc_params.next_gps = GPSTATUS_FERMION_MEASURES;
        if(mc_params.measures_done == fm_par.SingleInvNVectors){
					mc_params.next_gps = GPSTATUS_UPDATE;
					mc_params.measures_done = 0;
        }

        // determining run condition
        if(0 == devinfo.myrank && RUN_CONDITION_TERMINATE != mc_params.run_condition){
        
					// program exits if it finds a file called "stop"

					FILE * test_stop = fopen("stop","r");
					if(test_stop){
						fclose(test_stop);
						printf("File  \'stop\' found, stopping cycle now.\n");
						mc_params.run_condition = RUN_CONDITION_TERMINATE;
					}

					// program exits if time is running out
            
					struct timeval now;
					gettimeofday(&now,NULL);
					double total_duration = (double) 
						(now.tv_sec - mc_params.start_time.tv_sec)+
						(double)(now.tv_usec - mc_params.start_time.tv_usec)/1.0e6;

					double max_expected_duration_with_another_cycle;
					if(GPSTATUS_UPDATE == mc_params.next_gps){
						max_expected_duration_with_another_cycle = 
							total_duration + 1.3*mc_params.max_update_time;
						printf("Next step, update : %ds\n",(int) mc_params.max_update_time);
					}
					if(GPSTATUS_FERMION_MEASURES == mc_params.next_gps){
						max_expected_duration_with_another_cycle = 
							total_duration + 2*mc_params.max_flavour_cycle_time;
						printf("Next step, flavour measure cycle : %ds\n",
									 (int) mc_params.max_flavour_cycle_time);
					}

					if(max_expected_duration_with_another_cycle > mc_params.MaxRunTimeS){
						printf("Time is running out (%d of %d seconds elapsed),",
									 (int) total_duration, (int) mc_params.MaxRunTimeS);
						printf(" shutting down now.\n");
						printf("Total max expected duration: %d seconds",
									 (int) max_expected_duration_with_another_cycle);
						printf("(%d elapsed now)\n",(int) total_duration);
						// https://www.youtube.com/watch?v=MfGhlVcrc8U
						// but without that much pathos
						mc_params.run_condition = RUN_CONDITION_TERMINATE;
					}

					// program exits if MaxConfIdIter is reached
					if(conf_id_iter >= mc_params.MaxConfIdIter ){

						printf("%s - MaxConfIdIter=%d reached, job done!",
									 devinfo.myrankstr, mc_params.MaxConfIdIter);
						printf("%s - shutting down now.\n", devinfo.myrankstr);
						mc_params.run_condition = RUN_CONDITION_TERMINATE;
					}
					// program exits if MTraj is reached
					if( id_iter >= (mc_params.ntraj+id_iter_offset)){
						printf("%s - NTraj=%d reached, job done!",
									 devinfo.myrankstr, mc_params.ntraj);
						printf("%s - shutting down now.\n", devinfo.myrankstr);
						mc_params.run_condition = RUN_CONDITION_TERMINATE;
					}
					if (0==mc_params.ntraj) {
						printf("%s - NTraj=%d reached, job done!",
									 devinfo.myrankstr, mc_params.ntraj);
						printf("%s - shutting down now.\n", devinfo.myrankstr);
						mc_params.run_condition = RUN_CONDITION_TERMINATE;
					}
        }
#ifdef MULTIDEVICE

        MPI_Bcast((void*)&(mc_params.run_condition),1,MPI_INT,0,devinfo.mpi_comm);
        printf("MPI%02d - Broadcast of run condition %d from master...\n",
							 devinfo.myrank, mc_params.run_condition);
        MPI_Bcast((void*)&(mc_params.next_gps),1,MPI_INT,0,devinfo.mpi_comm);
        printf("MPI%02d - Broadcast of next global program status %d from master...\n",
							 devinfo.myrank, mc_params.next_gps);

#endif
      } // while id_iter loop ends here             
  } // closes if (0 != mc_params.ntraj)
    
  // saving gauge conf and RNG status to file
  for(int r=0;r<rep->replicas_total_number;r++){
#ifdef PAR_TEMP
    snprintf(rep_str,20,"replica_%d",r); // initialize rep_str
    strcat(mc_params.save_conf_name,rep_str); // append rep_str
#endif
		
    if (debug_settings.SaveAllAtEnd){
      printf("MPI%02d - Saving conf %s.\n", devinfo.myrank,
						 mc_params.save_conf_name);
      save_conf_wrapper(conf_acc[r],mc_params.save_conf_name, conf_id_iter,
												debug_settings.use_ildg);
    }else printf("\n\nMPI%02d: WARNING, \'SaveAllAtEnd\'=0,NOT SAVING/OVERWRITING CONF AND RNG STATUS.\n\n\n", devinfo.myrank);
#ifdef PAR_TEMP
    strcpy(mc_params.save_conf_name,aux_name_file);
#endif
	} // end replicas loop.
	
	if (debug_settings.SaveAllAtEnd){
		printf("MPI%02d - Saving rng status in %s.\n", devinfo.myrank, 
					 mc_params.RandGenStatusFilename);
		saverand_tofile(mc_params.RandGenStatusFilename);
	}

  if(0 == devinfo.myrank && debug_settings.SaveAllAtEnd){
    save_global_program_status(mc_params); // WARNING: this function in some cases does not work
  }

  printf("MPI%02d: Double precision free [CORE]\n", devinfo.myrank);
  mem_free_core();
    
  printf("MPI%02d: Double precision free [EXTENDED]\n", devinfo.myrank);
  mem_free_extended();

  if(inverter_tricks.useMixedPrecision || md_parameters.singlePrecMD){
    printf("MPI%02d: Single precision free [CORE]\n", devinfo.myrank);
    mem_free_core_f();
  }
  if( md_parameters.singlePrecMD){
    printf("MPI%02d: Signle precision free [EXTENDED]\n", devinfo.myrank);
    mem_free_extended_f();
  }

#ifdef PAR_TEMP
  free(all_swap_vector);
  free(acceptance_vector);

  // freeing rep_info vectors
  free(rep->cr_vec);
  free(rep->label);
#endif
	
  printf("MPI%02d: freeing device nnp and nnm\n", devinfo.myrank);
	#pragma acc exit data delete(nnp_openacc)
	#pragma acc exit data delete(nnm_openacc)
	#pragma acc exit data delete(gl_stout_rho)
	#pragma acc exit data delete(gl_topo_rho)

  printf("\n  MPI%02d - Allocated memory before the shutdown: %zu \n\n\n",devinfo.myrank,memory_used);
  struct memory_allocated_t *all=memory_allocated_base;
    
  while(all!=NULL)
    {
      printf("\n  MPI%02d - To be deallocated: %s having size %zu (or maybe is to be counted)\n\n\n",devinfo.myrank,all->varname,all->size);
      // free_wrapper(all->ptr);
      all=all->next;
    };

#ifndef __GNUC__
  // OpenAcc context closing
  shutdown_acc_device(my_device_type);
#endif

#ifdef MULTIDEVICE
  shutdown_multidev();
#endif
    
  if(0==devinfo.myrank_world){printf("The End\n");}
  return(EXIT_SUCCESS);
}
