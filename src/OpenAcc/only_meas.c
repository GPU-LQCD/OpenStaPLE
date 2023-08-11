#define PRINT_DETAILS_INSIDE_UPDATE
#define ALIGN 128
#include "../Include/stringify.h"
#include <errno.h>
#include <time.h>

#include "../Include/memory_wrapper.h"
#include "../Include/tell_geom_defines.h"
#include "../Meas/gauge_meas.h"
#include "../Meas/polyakov.h"
#include "../Meas/measure_topo.h"
#include "./HPT_utilities.h"
#include "./plaquettes.h"
#include "./field_corr.h"
#include "./geometry.h"
#include "./deviceinit.h"
#include "./io.h"
#include "./struct_c_def.h"
#include "./su3_measurements.h"
#include "./su3_utilities.h"
#include "./alloc_settings.h"
#include "./cooling.h"
#include <unistd.h>

#ifdef __GNUC__
#include "sys/time.h"
#endif

#include "../Mpi/multidev.h"
#if NRANKS_D3 > 1
#include <mpi.h>
#include "../Mpi/communications.h"
#endif

#define MAXLINES 300
#define MAXLINELENGTH 500

// global variables
char input_file_str[MAXLINES*MAXLINELENGTH];
int conf_id_iter;
int verbosity_lv;

// fool defintion to avoid include fermion_parameters.h and
// insert dependence on action.{c,h}
action_param act_params;
typedef struct ferm_param_t{
  double ferm_mass;
} ferm_param;
ferm_param *fermions_parameters;
///// end fool definitions /////

#define ALLOCCHECK(control_int,var)																			\
  if(control_int != 0 )	printf("MPI%02d: \tError in  allocation of %s . \n",devinfo.myrank, #var); \
	else printf("MPI%02d: \tAllocation of %s : OK , %p\n",								\
							devinfo.myrank, #var, var );															\
	
#define FREECHECK(var)														\
	printf("\tFreed %s ...", #var); fflush(stdout);	\
	free_wrapper(var); printf(" done.\n");

int main(int argc, char **argv){
	if(argc!=3)
		{
			printf("ERROR: wrong usage.\nCorrect syntax: mpirun -n <nranks> %s <conf_list_file_name> <output_file_name>\n",argv[0]);
			exit(1);
		}

#if NRANKS_D3 > 1
	pre_init_multidev1D(&devinfo);
	gdbhook();
#endif		

#if NRANKS_D3 > 1       
	devinfo.single_dev_choice  = 0;
	devinfo.async_comm_fermion = 1;
	devinfo.async_comm_gauge   = 1;
	devinfo.proc_per_node      = NRANKS_D3; // It only works on one node
#endif
	devinfo.nranks_read=devinfo.nranks;
	printf("devinfo.nranks: %d\n",devinfo.nranks);
	
#if NRANKS_D3 > 1        
	init_multidev1D(&devinfo);
#else
  devinfo.myrank = 0;
  devinfo.nranks = 1;
#endif
	
	//////  OPENACC CONTEXT INITIALIZATION    ////
	acc_device_t my_device_type = acc_device_nvidia;
	printf("MPI%02d: Selecting device.\n", devinfo.myrank);
#if NRANKS_D3 > 1
	printf("devinfo.single_dev_choice: %d\ndevinfo.myrank: %d\n devinfo.proc_per_node: %d\n", devinfo.single_dev_choice, devinfo.myrank, devinfo.proc_per_node);
	select_init_acc_device(my_device_type, (devinfo.single_dev_choice + devinfo.myrank)%devinfo.proc_per_node);
#else
	select_init_acc_device(my_device_type, devinfo.single_dev_choice);
#endif
	printf("Device Selected : OK \n");


	su3_soa * conf_acc;
	su3_soa * aux_conf_acc;
	su3_soa * field_corr;
	su3_soa * field_corr_aux;
	single_su3 * closed_corr;
	single_su3 * loc_plaq_aux;
	global_su3_soa * conf;
	d_complex * corr ;
	dcomplex_soa * local_sums;
	int allocation_check;

	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&conf, ALIGN,8*sizeof(global_su3_soa));
	ALLOCCHECK(allocation_check, conf);

	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa)); 
	ALLOCCHECK(allocation_check, conf_acc);
#pragma acc enter data create(conf_acc[0:8])

	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&aux_conf_acc, ALIGN, 8*sizeof(su3_soa));
	ALLOCCHECK(allocation_check, aux_conf_acc);
#pragma acc enter data create(aux_conf_acc[0:8])

	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&local_sums, ALIGN, 2*sizeof(dcomplex_soa));
	ALLOCCHECK(allocation_check, local_sums) ;
#pragma acc enter data create(local_sums[0:2])

	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&field_corr, ALIGN, 8*sizeof(su3_soa)); 
	ALLOCCHECK(allocation_check, field_corr );
#pragma acc enter data create(field_corr[0:8])

	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&field_corr_aux, ALIGN, 8*sizeof(su3_soa));
	ALLOCCHECK(allocation_check, field_corr_aux );
#pragma acc enter data create(field_corr_aux[0:8])

	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&loc_plaq_aux, ALIGN, 8*sizeof(su3_soa));
	ALLOCCHECK(allocation_check, loc_plaq_aux);
#pragma acc enter data create(loc_plaq_aux[0:8])

	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&corr, ALIGN, nd0*sizeof(d_complex));
	ALLOCCHECK(allocation_check, corr);
#pragma acc enter data create(corr[0:nd0])

	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&closed_corr, ALIGN, sizeof(single_su3));
	ALLOCCHECK(allocation_check, closed_corr);
#pragma acc enter data create(closed_corr[0:1])

	if(0==devinfo.myrank) print_geom_defines(); 
	compute_nnp_and_nnm_openacc();
#pragma acc enter data copyin(nnp_openacc) 
#pragma acc enter data copyin(nnm_openacc) 

	if(0==devinfo.myrank)
		{
			printf("HARDCODED LATTICE DIMENSIONS:\n");
			printf("GL_N0: %d\n", GL_N0) ;
			printf("GL_N1: %d\n", GL_N1) ;
			printf("GL_N2: %d\n", GL_N2) ;
			printf("GL_N3: %d\n", GL_N3) ;
		}
	
	geom_par.gnx = GL_N0 ;
	geom_par.gny = GL_N1 ;
	geom_par.gnz = GL_N2 ;
	geom_par.gnt = GL_N3 ;
		
	geom_par.xmap = 0 ;
	geom_par.ymap = 1 ;
	geom_par.zmap = 2 ;
	geom_par.tmap = 3 ;

	set_geom_glv(&geom_par);

	// read conf paths from input and count number of configurations.
	if(0==devinfo.myrank)
		printf("Reading input...\n");
	long int confmax;
	char ** confs = read_list_of_confs(argv[1], &confmax);

	
	FILE *fplq;
	char filename[100];
	sprintf(filename,"%s_plaq",argv[2]);
	// open output file only by master rank
	if(0==devinfo.myrank)
		{
			printf("Opening file %s\n",filename);
			fplq=fopen(filename, "w");
			if (fplq == NULL)
				{
					printf("Could not open file\n");
					exit(1);
				}
			// printing header of measurements file.
			fprintf(fplq,"#conf_id ReTrPlq\n");
			fflush(fplq);			
		}

#ifdef MEASURE_PLAQ_CORR
	FILE *fcor;
	char* filename;
	sprintf(filename,"%s_plaq",argv[2]);
	// open output file only by master rank
	if(0==devinfo.myrank)
		{
			fcor=fopen(filename, "w");
			if (fcor == NULL)
				{
					printf("Could not open file\n");
					exit(1);
				}
			// printing header of measurements file.
			fprintf(fcor,"#conf_id cooling_step L D_para D_perp\n");
			fflush(fcor);
		}
	
	// hardcoded number of cooling steps. Could become a command line argument.
	int maxstep=100;
#endif

	int conf_id;	
	for(int conf_num=0; conf_num<confmax; conf_num++){

		// reading gauge conf
		if(0==devinfo.myrank)
			printf("Reading file %s\n",confs[conf_num]);
		if(! local_confrw_read_conf_wrapper(conf, conf_acc, confs[conf_num], &conf_id, 1))
			{
				printf("MPI%02d - Gauge Conf \"%s\" Read : OK \n",
							 devinfo.myrank, confs[conf_num]);
      }
		else
			{
				printf("MPI%02d - Gauge Conf \"%s\" Read : FAILED \n\nABORTING",
							 devinfo.myrank, confs[conf_num]);
#if NRANKS_D3 > 1
				MPI_Abort(MPI_COMM_WORLD,0);
#endif
			}
#pragma acc update device(conf_acc[0:8])

		// performing measurements
		if(0==devinfo.myrank)
			printf("Computing plaquette\n");

		double plq;
		plq = calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);

		if(0==devinfo.myrank)
			fprintf(fplq,"%d\t%.18lf\n", conf_id, plq/GL_SIZE/6.0/3.0);

#ifdef MEASURE_PLAQ_CORR
		if(0==devinfo.myrank)
			printf("Computing correlators\n");
		double  D_paral[nd0], D_perp[nd0], D_temp_paral[nd0], D_temp_perp[nd0];
		calc_field_corr(conf_acc, field_corr, field_corr_aux, aux_conf_acc, local_sums, corr, closed_corr, D_paral, D_perp, D_temp_paral, D_temp_perp); 

		if(0==devinfo.myrank)
			for(int L=0; L<nd0/2; L++)
			  fprintf(fcor,"%d\t%d\t%d\t%.18lf\t%.18lf\t%.18lf\t%.18lf\\n", conf_id, 0, L+1,
								D_paral[L]/((double)18*GL_SIZE),  D_perp[L]/((double)18*GL_SIZE),
								D_temp_paral[L]/((double)6*GL_SIZE),  D_temp_perp[L]/((double)6*GL_SIZE));			

		for(int coolstep=1; coolstep<=maxstep; coolstep++){
			cool_conf(conf_acc, conf_acc, aux_conf_acc);

			calc_field_corr(conf_acc, field_corr, field_corr_aux, aux_conf_acc, local_sums, corr, closed_corr, D_paral, D_perp, D_temp_paral, D_temp_perp);
			
			if(0==devinfo.myrank)
				for(int L=0; L<nd0/2; L++)
					fprintf(fcor,"%d\t%d\t%d\t%.18lf\t%.18lf\t%.18lf\t%.18lf\\n", conf_id, coolstep, L+1,
									D_paral[L]/((double)18*GL_SIZE),  D_perp[L]/((double)18*GL_SIZE),
									D_temp_paral[L]/((double)6*GL_SIZE),  D_temp_perp[L]/((double)6*GL_SIZE));			
			
		}// close coolstep
#endif
		// end of measurements
		
		if(0==devinfo.myrank)
			printf("%s analysis completed. Progress: %d/%d\n", confs[conf_num],conf_num+1,confmax);
	}// close conf_num
	if(0==devinfo.myrank){
		fclose(fplq);
#ifdef MEASURE_PLAQ_CORR
		fclose(fcor);
#endif
	}

#pragma acc exit data delete(nnp_openacc)
#pragma acc exit data delete(nnm_openacc)					 

	
	FREECHECK(closed_corr);
#pragma acc exit data delete(closed_corr)
	
	FREECHECK(corr);
#pragma acc exit data delete(corr)

	FREECHECK(loc_plaq_aux);
#pragma acc exit data delete(loc_plaq_aux)

	FREECHECK(field_corr_aux );
#pragma acc exit data delete(field_corr_aux)

	FREECHECK(field_corr );
#pragma acc exit data delete(field_corr)

	FREECHECK(local_sums);
#pragma acc exit data delete(local_sums)

	FREECHECK(aux_conf_acc);
#pragma acc exit data delete(aux_conf_acc)

	FREECHECK(conf_acc);
#pragma acc exit data delete(conf_acc)

	FREECHECK(conf);

#ifndef __GNUC__
  // OpenAcc context closing
  shutdown_acc_device(my_device_type);
#endif

#ifdef MULTIDEVICE
  shutdown_multidev();
#endif
    
  if(0==devinfo.myrank){printf("The End\n");}
  return(EXIT_SUCCESS);
}
