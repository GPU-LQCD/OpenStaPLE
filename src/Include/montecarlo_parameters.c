#ifndef _MONTECARLO_PARAMS_C_
#define _MONTECARLO_PARAMS_C_

#include <stdlib.h> // this apparently useless include must be here.
                    // exactly here, before 'include "./montecarlo_parameters.h'.
                    // Otherwise, printf, fscanf and fprintf DO NOT 
                    // work correctly while operating on the MC_PARAM_T structure,
                    // while compiling with PGCC.
                    // Not moving it into 'montecarlo_parameters.h' itself because
                    // I don't completely understand the issue here 
                    // trying to be as minimal as possible
#include "./montecarlo_parameters.h" 
#include <stdio.h>                   
#include "./common_defines.h"        
#include "../Mpi/multidev.h"         



mc_params_t mc_params;

void init_global_program_status(){
    
	MPI_PRINTF1("reading global program status  from file %s\n", mc_params.statusFileName);

	FILE * gps_file = fopen(mc_params.statusFileName, "r");  
	if(gps_file){
		int reads = fscanf(gps_file,"%d %f %f %d",
											 &(mc_params.next_gps),
											 &(mc_params.max_flavour_cycle_time),
											 &(mc_params.max_update_time),
											 &(mc_params.measures_done));
		if (reads != 4 ){
			if(0 == devinfo.myrank)
				printf("ERROR: %s:%d: montecarlo status file %s not readable\n",
							 __FILE__,__LINE__, mc_params.statusFileName);

		}

		fclose(gps_file);
	}
	else 
    {
			MPI_PRINTF1("file %s not readable\n", mc_params.statusFileName);

			mc_params.next_gps  = GPSTATUS_UPDATE;
			mc_params.max_flavour_cycle_time = 1.0f;
			mc_params.max_update_time = 1.0f;
			mc_params.measures_done = 0 ;
    }

	mc_params.run_condition = RUN_CONDITION_GO; 
	printf("run_condition %d ",mc_params.run_condition);

	printf("%d %f %f %d\n",
				 mc_params.next_gps,
				 mc_params.max_flavour_cycle_time,
				 mc_params.max_update_time,
				 mc_params.measures_done);

	printf("#mc_params.next_gps,mc_params.max_flavour_cycle_time,\n\
#mc_params.max_update_time,mc_params.measures_done\n");


}

void save_global_program_status(mc_params_t mcp, int max_update_times, int max_flavour_cycle_times){

	printf("Saving global program status...\n");
	printf("%d %f %f %d\n",
				 mcp.next_gps,
				 max_flavour_cycle_times,
				 max_update_times,
				 mcp.measures_done);

	printf("#mc_params.next_gps,mc_params.max_flavour_cycle_time,\n\
#mc_params.max_update_time,mc_params.measures_done\n");


	FILE * gps_file = fopen(mc_params.statusFileName, "w");  
	fprintf(gps_file,"%d %f %f %d\n",
					mcp.next_gps,
					max_flavour_cycle_times,
					max_update_times,
					mcp.measures_done);

	fprintf(gps_file,"#mc_params.next_gps,mc_params.max_flavour_cycle_time,\n\
#mc_params.max_update_time,mc_params.measures_done\n");

	fclose(gps_file);


}




#endif
