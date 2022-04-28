#ifndef IO_H_
#define IO_H_

#include <stdio.h>
#include <stdlib.h>
#include "./struct_c_def.h"
#include "./alloc_vars.h"

#ifdef MULTIDEVICE
#include "mpi.h"
#endif 

#include "../Mpi/multidev.h"
#include "../Mpi/communications.h"


int print_su3_soa_ASCII(global_su3_soa * const conf, const char* nomefile,
        int conf_id_iter,double_soa *back_phases);


int print_su3_soa_ildg_binary(global_su3_soa * conf, const char* nomefile,
        int conf_id_iter);


//WRAPPER
static inline int save_conf(global_su3_soa * const conf_rw, const char* nomefile,
        int conf_id_iter, int use_ildg)
{

#ifdef MULTIDEVICE
    if(devinfo.myrank != 0){
        printf("MPI%02d: Rank is not allowed to use this function!\n",devinfo.myrank );
        printf("ERROR: %s:%d\n",__FILE__, __LINE__);
        MPI_Abort(MPI_COMM_WORLD,1);
    }
#endif 

    if(use_ildg){
        printf("MPI%02d: Using ILDG format.\n",devinfo.myrank);
        return print_su3_soa_ildg_binary(conf_rw,nomefile,conf_id_iter);
    }
    else {
        printf("MPI%02d: Using ASCII format.\n", devinfo.myrank);
        return print_su3_soa_ASCII(conf_rw,nomefile,conf_id_iter,NULL);
    }

}


// WRAPPER OF WRAPPER
static inline void save_conf_wrapper(su3_soa* conf, const char* nomefile,
        int conf_id_iter, int use_ildg)
 {

    printf("MPI%02d - Saving whole configuration...\n", devinfo.myrank);
    int writeOutcome;

#ifdef MULTIDEVICE

    if(devinfo.myrank == 0){
        int irank;
        for(irank = 1 ; irank < devinfo.nranks; irank++)
            recv_loc_subconf_from_rank(conf_rw,irank,irank);

        recv_loc_subconf_from_buffer(conf_rw,conf,0);
        writeOutcome = save_conf(conf_rw, nomefile,conf_id_iter, use_ildg);
    }
    else send_lnh_subconf_to_master(conf,devinfo.myrank);
    
    MPI_Bcast((void*) &writeOutcome,1,MPI_INT,0,MPI_COMM_WORLD);

    if(0 != writeOutcome){
        MPI_Finalize();
    }
#else 
    recv_loc_subconf_from_buffer(conf_rw,conf,0);
    writeOutcome = save_conf(conf_rw, nomefile,conf_id_iter, use_ildg);
#endif
    if(0 != writeOutcome){ 
        printf("MPI%02d: Something went wrong in saving conf. Terminating.\n",devinfo.myrank );
        exit(1);
    }


}


int read_su3_soa_ASCII(global_su3_soa * conf, const char* nomefile,int * conf_id_iter );
int read_su3_soa_ildg_binary(global_su3_soa * conf, const char* nomefile,
        int * conf_id_iter );


//WRAPPER
static inline int read_conf(global_su3_soa * conf, const char* nomefile,
        int * conf_id_iter, int use_ildg )
{

#ifdef MULTIDEVICE
    if(devinfo.myrank != 0){
        printf("MPI%02d: Rank is not allowed to use this function!\n",devinfo.myrank );
        printf("ERROR: %s:%d\n",__FILE__, __LINE__);
        exit(1);
    }
#endif 

    printf("Reading configuration. ");
    if(use_ildg){ 
        printf("Using ILDG format.\n");
        return read_su3_soa_ildg_binary(conf,nomefile,conf_id_iter);
    }
    else{
        printf("Using ascii format.\n");
        return read_su3_soa_ASCII(conf,nomefile,conf_id_iter );
    }


}


// WRAPPER OF WRAPPER
static inline int read_conf_wrapper(su3_soa* conf, const char* nomefile,
        int * conf_id_iter, int use_ildg)
{

    int error =  0;
#ifdef MULTIDEVICE

    if(devinfo.myrank == 0){
        if(verbosity_lv > 2)
            printf("MPI%02d - reading global conf \n",devinfo.myrank );
        error = read_conf(conf_rw, nomefile,conf_id_iter, use_ildg);
        MPI_Bcast((void*) &error,1,MPI_INT,0,MPI_COMM_WORLD);

        if(!error) {
            send_lnh_subconf_to_buffer(conf_rw,conf,0);
            int irank;
            for(irank = 1 ; irank < devinfo.nranks; irank++)
                send_lnh_subconf_to_rank(conf_rw,irank);
        
            MPI_Bcast((void*) conf_id_iter,1,MPI_INT,0,MPI_COMM_WORLD);

        }
        else 
        if(verbosity_lv > 2)
            printf("MPI%02d - no conf sent!\n",devinfo.myrank );

    }
    else{
        if(verbosity_lv > 2)
            printf("MPI%02d - receiving conf \n",devinfo.myrank );
        
        MPI_Bcast((void*) &error,1,MPI_INT,0,MPI_COMM_WORLD);
        if(!error){ 
            receive_lnh_subconf_from_master(conf);
            MPI_Bcast((void*) conf_id_iter,1,MPI_INT,0,MPI_COMM_WORLD);
        }
        else 
        if(verbosity_lv > 2)
            printf("MPI%02d - no conf received!\n",devinfo.myrank );
    }


    if(2 == error) MPI_Finalize();

#else 
    error = read_conf(conf_rw, nomefile,conf_id_iter, use_ildg);
    if(0 == error) send_lnh_subconf_to_buffer(conf_rw,conf,0);
#endif
    if(2 == error){
        printf("MPI%02d: Configuration exists, but errors in reading it. Terminating now.\n",
                devinfo.myrank);
        exit(1);
    }

    return error;

}

void rw_iterate_on_global_sites_lx_xyzt_axis_ordering( 
        void (*single_element_rw)(
            const int /*idxh*/, const int /*parity*/, 
            const int /*direction*/, const void* /*data*/,
            const int /*conf_machine_endianness_disagreement*/,
            FILE * /*fp*/), 
        const void* datastruct, FILE * fp,
        const int scalar_even_mode);

void binaryread_single_su3_into_su3_soa( // machine big endian
        const int idx, const int parity, const int direction, 
        const void* datastruct,
        const int conf_machine_endianness_disagreement, FILE* fp);

void binarywrite_single_su3_into_su3_soa(
        const int idx, const int parity, const int direction, 
        const void* datastruct, 
        const int conf_machine_endianness_disagreement, FILE* fp);


// READ AND WRITE FOR THE POOR MAN
//void write_conf_binary_chunks(su3_soa * conf, const char *rootname,int nchunks);
//void read_conf_binary_chunks(su3_soa * conf, const char *rootname,int chunk);

#endif 
