
#include "../OpenAcc/io.h"
#include "../OpenAcc/struct_c_def.h"
#include "../OpenAcc/geometry.h"
#include "../Mpi/geometry_multidev.h"
#ifdef __GNUC__
#include <stdlib.h>
#endif
#include <stdio.h>
#include <string.h>

#include "../Include/stringify.h"

int verbosity_lv;

#define TOILDG 0
#define TOASCII 1

#include "../OpenAcc/action.h"
#include "../Include/fermion_parameters.h"
#include "../OpenAcc/alloc_settings.h"
#include "../OpenAcc/alloc_vars.h"
#include "../OpenAcc/backfield.h"
#include "../OpenAcc/backfield_parameters.h"

int main(int argc, char ** argv){

    if(argc != 4){
        printf("USAGE: %s conf_input_name conf_output_name MODE\n",argv[0]);
        printf("Where MODE can be:\n");
        printf("TOILDG if you want to convert from ASCII to ILDG\n");
        printf("TOASCII if you want to convert from ILDG to ASCII\n");
        exit(1);
    }

    global_su3_soa * conf = (global_su3_soa * ) malloc(8*sizeof(global_su3_soa));
    int conf_id_iter;
    int mode; // TOILDG or TOASCII

    if(0 == strcmp(argv[3],"TOILDG")) mode = TOILDG ;
    else if(0 == strcmp(argv[3],"TOASCII"))  mode = TOASCII;
    else {
        printf("Option %s not recognized: please choose either TOILDG or TOASCII\n",
            argv[3]);
        exit(1);
        }



    printf("HARDCODED LATTICE DIMENSIONS:\n");
    printf("GL_N0: %d\n", GL_N0) ; 
    printf("GL_N1: %d\n", GL_N1) ; 
    printf("GL_N2: %d\n", GL_N2) ; 
    printf("GL_N3: %d\n", GL_N3) ; 


    geom_par.gnx = GL_N0 ;
    geom_par.gny = GL_N1 ;
    geom_par.gnz = GL_N2 ;
    geom_par.gnt = GL_N3 ;

    geom_par.xmap = 0 ;
    geom_par.ymap = 1 ;
    geom_par.zmap = 2 ;
    geom_par.tmap = 3 ;

    set_geom_glv(&geom_par);


    printf("Reading conf %s...\n", argv[1]);
    int readcheck;

    if(TOILDG == mode)
        readcheck = read_su3_soa_ASCII(conf,argv[1],&conf_id_iter );
    else if (TOASCII == mode )
        readcheck = read_su3_soa_ildg_binary(conf,argv[1],&conf_id_iter);

    if(readcheck){
        printf("There are issues, closing now...\n");
        exit(1);
    }
          


    
    act_params.beta = 0.0;
    // allocating one fermion for convenience
    fermions_parameters = malloc(sizeof(ferm_param));
    fermions_parameters[0].ferm_mass = 0.0;
    fermions_parameters[0].ferm_charge = 0.0;
    alloc_info.NDiffFlavs = 1;

    // setting EM field to zero
    backfield_parameters.ex = backfield_parameters.ey = backfield_parameters.ez = 0;
    backfield_parameters.bx = backfield_parameters.by = backfield_parameters.bz = 0;

    // allocating and preparing staggered phases field
    u1_back_phases = malloc(8*sizeof(double_soa));

    calc_u1_phases(u1_back_phases,backfield_parameters,0,0);


    
    printf("conf_id_iter is %d.\n", conf_id_iter);

    printf("Writing conf in file %s \n", argv[2]);\

    if(TOILDG == mode)
        print_su3_soa_ildg_binary(conf,argv[2],conf_id_iter);
    else if (TOASCII == mode )
        print_su3_soa_ASCII(conf,argv[2],conf_id_iter,u1_back_phases);

    free(conf);
    free(fermions_parameters);
    printf("(commit: %s)\n", xstr(COMMIT_HASH) );

}
