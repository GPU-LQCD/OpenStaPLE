#ifndef MULTIDEV_C_
#define MULTIDEV_C_

#include <stdio.h>
#include <stdlib.h>
#include "./geometry_multidev.h"
#include "../OpenAcc/geometry.h"
#include "../Include/rep_info.h"

#include "./multidev.h"

dev_info devinfo;

#ifdef MULTIDEVICE
#include <mpi.h>

extern int verbosity_lv;


void pre_init_multidev1D(dev_info * mdi)
{
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &(mdi->myrank_world));
    MPI_Comm_size(MPI_COMM_WORLD, &(mdi->nranks_world));
    MPI_Get_processor_name(mdi->processor_name,&(mdi->namelen));

    // associate to specific replica communication group
    if(mdi->nranks_world/NRANKS_D3>1){
      mdi->replica_idx = mdi->myrank_world/NRANKS_D3;
      MPI_Comm_split(MPI_COMM_WORLD, mdi->replica_idx, mdi->myrank_world, &(mdi->mpi_comm));
      MPI_Comm_rank(mdi->mpi_comm,&(mdi->myrank));
      MPI_Comm_size(mdi->mpi_comm,&(mdi->nranks));

      int color = (mdi->myrank==0)? 0: MPI_UNDEFINED;
      MPI_Comm_split(MPI_COMM_WORLD, color, mdi->replica_idx, &(mdi->mpi_comm_salamino));
    }else{
      mdi->replica_idx=0;
      mdi->mpi_comm=MPI_COMM_WORLD;
      mdi->myrank=mdi->myrank_world;
      mdi->nranks=mdi->nranks_world;
    }


    if(mdi->nranks != NRANKS_D3){
        MPI_PRINTF0("NRANKS_D3 is different from nranks: no salamino? Exiting now\n");
        MPI_PRINTF1("NRANKS_D3 = %d, nranks = %d\n",mdi->myrank, NRANKS_D3);
        exit(1);
    }

    if(mdi->nranks_world != NREPLICAS*NRANKS_D3){
        MPI_PRINTF0("NREPLICAS is different from nranks_world/nranks. Exiting now\n");
        MPI_PRINTF1("NREPLICAS = %d, nranks_world/nranks = %d\n",NREPLICAS, mdi->nranks_world/mdi->nranks);
        exit(1);
    }
    
    if(verbosity_lv > 2){
        MPI_PRINTF0("- Called MPI_Init\n");
    }



}

void init_multidev1D(dev_info * mdi)
{
    mdi->myrank_L = (mdi->myrank + (mdi->nranks-1))%mdi->nranks;//SALAMINO
    mdi->myrank_R = (mdi->myrank + 1 ) % mdi->nranks;       //SALAMINO
    mdi->node_subrank =  mdi->myrank % mdi->proc_per_node;   //SALAMINO

    sprintf(mdi->myrankstr,"MPI%02d",mdi->myrank);

    MPI_PRINTF1("of \"%02d\" tasks running on host \"%s\", local rank: %d, rankL: %d, rankR: %d\n", mdi->processor_name, 
    mdi->node_subrank, mdi->myrank_L, mdi->myrank_R); //SALAMINO

    mdi->myrank4int = xyzt_rank(mdi->myrank); // ALL RANKS


    mdi->myrank4int = xyzt_rank(mdi->myrank);
    int myrank_vect[4];
    myrank_vect[0] = mdi->myrank4int.d0 ;
    myrank_vect[1] = mdi->myrank4int.d1 ; 
    myrank_vect[2] = mdi->myrank4int.d2 ; 
    myrank_vect[3] = mdi->myrank4int.d3 ; 

    int dir;
    for(dir=0;dir<4;dir++){
        int lrank_dir = (myrank_vect[dir]-1);
        if(lrank_dir == -1) lrank_dir = geom_par.nranks[dir]-1;
        int rrank_dir = (myrank_vect[dir]+1)% geom_par.nranks[dir];
        mdi->nnranks[dir][0] = lrank_dir;
        mdi->nnranks[dir][1] = rrank_dir;
    }

    mdi->gl_loc_origin4int = gl_loc_origin_from_rank(mdi->myrank);

    mdi->halo_widths0123[0] = D0_HALO ;  
    mdi->halo_widths0123[1] = D1_HALO ;  
    mdi->halo_widths0123[2] = D2_HALO ;  
    mdi->halo_widths0123[3] = D3_HALO ;  
    mdi->origin_0123[0]     = mdi->gl_loc_origin4int.d0;
    mdi->origin_0123[1]     = mdi->gl_loc_origin4int.d1;
    mdi->origin_0123[2]     = mdi->gl_loc_origin4int.d2;
    mdi->origin_0123[3]     = mdi->gl_loc_origin4int.d3;
    if(verbosity_lv > 2){
        MPI_PRINTF0("- Finished init_multidev1D\n");   
        MPI_PRINTF1("- Origin(%d,%d,%d,%d)",
                mdi->origin_0123[0],mdi->origin_0123[1],
                mdi->origin_0123[2], mdi->origin_0123[3]);
   
    }


}






void shutdown_multidev()
{

    MPI_PRINTF0("Finalizing...\n" );
    MPI_Finalize();     

}





#endif
#endif
