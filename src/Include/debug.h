#ifndef DEBUG_H_
#define DEBUG_H_

typedef struct DEBUG_SETTINGS_T{

    int use_ildg;
    int input_vbl; // verbosity level
    int save_diagnostics;
    char diagnostics_filename[200];
    int SaveAllAtEnd;
    int do_reversibility_test;
    int do_norandom_test;
    int rng_fakeness_level ; 
    int md_dbg_print_max_count;
    int print_bfield_dbginfo;


}debug_settings_t;

extern debug_settings_t debug_settings; 
extern int md_dbg_print_count;


#endif
