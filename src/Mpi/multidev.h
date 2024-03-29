#ifndef MULTIDEV_H_
#define MULTIDEV_H_

#include "./geometry_multidev.h"

#ifdef MULTIDEVICE
#include <mpi.h>
#endif

typedef struct dev_info_t{

    // FROM MPI INIT AND SIMILAR
    
    int single_dev_choice; // from input file
    int myrank_world, nranks_world;
    int myrank, nranks;
    int replica_idx, num_replicas;
    char myrankstr[50];
    int nranks_read;

#ifdef MULTIDEVICE
    MPI_Comm mpi_comm;
    int async_comm_fermion;
    int async_comm_gauge;

    int myrank_L,myrank_R; // nearest neighbours, salamino
    int node_subrank; // which card it will use?
    char processor_name[MPI_MAX_PROCESSOR_NAME]; int namelen;

    // GEOMETRIC QUANTITIES
    vec4int myrank4int;
    int nnranks[4][2];// ranks of nearest neighbour sublattices
    vec4int gl_loc_origin4int;
    // FROM INPUT FILE
    int proc_per_node;
    int halo_widths0123[4];
    int origin_0123[4];
#endif


} dev_info ; 


extern dev_info devinfo;

#define MPI_PRINTF0(fmt) { \
  if(devinfo.num_replicas>1){ \
    printf("MPI%02d:%02d - " fmt ,devinfo.replica_idx, devinfo.myrank); \
  }else{ \
    printf("MPI%02d - " fmt ,devinfo.myrank); \
  }}
#define MPI_PRINTF1(fmt,...) { \
  if(devinfo.num_replicas>1){ \
    printf("MPI%02d:%02d - " fmt ,devinfo.replica_idx, devinfo.myrank, __VA_ARGS__); \
  }else{ \
    printf("MPI%02d - " fmt ,devinfo.myrank, __VA_ARGS__); \
  }}

//TODO: this doesn't work for some reason (it works with gcc), find out why
// #define MPI_PRINTF(fmt,...) printf("MPI%02d:%02d - " fmt ,devinfo.replica_idx, devinfo.myrank, ##__VA_ARGS__);  


#ifdef MULTIDEVICE

void pre_init_multidev1D(dev_info * mdi);
void init_multidev1D(dev_info * mdi);

void shutdown_multidev();

#endif

#endif
