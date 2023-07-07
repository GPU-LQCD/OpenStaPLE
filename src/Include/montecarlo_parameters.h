#ifndef _MONTECARLO_PARAMS_H_
#define _MONTECARLO_PARAMS_H_


#ifdef __GNUC__
#include "sys/time.h"
#endif

#define GPSTATUS_UPDATE 1 
#define RUN_CONDITION_TERMINATE 0
#define RUN_CONDITION_GO 1
#define GPSTATUS_FERMION_MEASURES 2

typedef struct MC_PARAM_T{
	int ntraj;
	int MaxConfIdIter; 
	double MaxRunTimeS;
	int therm_ntraj;
	int storeconfinterval;
	int saveconfinterval;
	char RandGenStatusFilename[200];
	char store_conf_name[200];
	char save_conf_name[200];
	int seed;
	double eps_gen;
	int JarzynskiMode;
	char statusFileName[50];

	// these are not read from input file
	int next_gps; //next_global_program_status
	struct timeval start_time;
	float max_flavour_cycle_time;
	float max_update_time;
	int measures_done;

	int run_condition;

} mc_params_t;

extern mc_params_t mc_params;

void init_global_program_status();
void save_global_program_status(mc_params_t mcp, int max_update_times, int max_flavour_cycle_times);

#endif
