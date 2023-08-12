#ifndef GEOMETRY_C_
#define GEOMETRY_C_

#include "./geometry.h"
#include "stdio.h"

#include "../Mpi/multidev.h"


geom_parameters geom_par = {.initialized_check = 0};


// this function well defines the boundaries conditions.
void compute_nnp_and_nnm_openacc(int nnp_openacc[sizeh][4][2],int nnm_openacc[sizeh][4][2])
{
  int d0, d1, d2, d3,parity;
	/*
		char nnfilename[50];
		#ifdef MULTIDEVICE
		sprintf(nnfilename,"nnfile_%s",devinfo.myrankstr);
		#else
		sprintf(nnfilename,"nnfile");
		#endif

		FILE * nnfile = fopen(nnfilename,"w");
	*/
  for(d3=0; d3<nd3; d3++) {
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
        for(d0=0; d0 < nd0; d0++) {
          int  d0m, d1m, d2m, d3m, d0p, d1p, d2p, d3p, idxh;
          idxh = snum_acc(d0,d1,d2,d3);
          parity = (d0+d1+d2+d3) % 2; // 0 even, 1 odd.

          d0m = d0 - 1; // the negative neighbours
          d1m = d1 - 1;
          d2m = d2 - 1;
          d3m = d3 - 1;
          d0m = d0m + (((d0m >> 31) & 0x1) * nd0);
          d1m = d1m + (((d1m >> 31) & 0x1) * nd1);
          d2m = d2m + (((d2m >> 31) & 0x1) * nd2);
          d3m = d3m + (((d3m >> 31) & 0x1) * nd3);

          nnm_openacc[idxh][0][parity] = snum_acc(d0m,d1,d2,d3);
          nnm_openacc[idxh][1][parity] = snum_acc(d0,d1m,d2,d3);
          nnm_openacc[idxh][2][parity] = snum_acc(d0,d1,d2m,d3);
          nnm_openacc[idxh][3][parity] = snum_acc(d0,d1,d2,d3m);

          d0p = d0 + 1; // the positive neighbours
          d1p = d1 + 1;
          d2p = d2 + 1;
          d3p = d3 + 1;
          d0p *= (((d0p-nd0) >> 31) & 0x1);
          d1p *= (((d1p-nd1) >> 31) & 0x1);
          d2p *= (((d2p-nd2) >> 31) & 0x1);
          d3p *= (((d3p-nd3) >> 31) & 0x1);
          
          nnp_openacc[idxh][0][parity] = snum_acc(d0p,d1,d2,d3);
          nnp_openacc[idxh][1][parity] = snum_acc(d0,d1p,d2,d3);
          nnp_openacc[idxh][2][parity] = snum_acc(d0,d1,d2p,d3);
          nnp_openacc[idxh][3][parity] = snum_acc(d0,d1,d2,d3p);



          if(NRANKS_D0 != 1){
						if(d0m == nd0-1 )  nnm_openacc[idxh][0][parity] = -1; 
						if(d0p == 0 )  nnp_openacc[idxh][0][parity] = -1;
          }  
          if(NRANKS_D1 != 1){
						if(d1m == nd1-1 )  nnm_openacc[idxh][1][parity] = -1; 
						if(d1p == 0 )  nnp_openacc[idxh][1][parity] = -1;
          }  
          if(NRANKS_D2 != 1){
						if(d2m == nd2-1 )  nnm_openacc[idxh][2][parity] = -1; 
						if(d2p == 0 )  nnp_openacc[idxh][2][parity] = -1;
          }  
          if(NRANKS_D3 != 1){
						if(d3m == nd3-1 )  nnm_openacc[idxh][3][parity] = -1;
						if(d3p == 0 )  nnp_openacc[idxh][3][parity] = -1;
          }  

					/*          
											fprintf(nnfile,"%d %d %d %d    %d    %d   ",d0,d1,d2,d3,parity,idxh);
											int dir;
											for(dir = 0; dir < 4 ;dir ++)
											fprintf(nnfile," %d ", nnm_openacc[idxh][dir][parity]);
											for(dir = 0; dir < 4 ;dir ++)
											fprintf(nnfile," %d ", nnp_openacc[idxh][dir][parity]);
											fprintf(nnfile,"\n");
					*/


        }
      }
    }
  }


	//  fclose(nnfile);
}


int set_geom_glv(geom_parameters* gp)
{
	gp->initialized_check = 1;

	int res = 0;

	// hardcoded
	gp->nd[0]= nd0; gp->nd[1] = nd1; gp->nd[2] = nd2; gp->nd[3] = nd3;

	gp->nloc[0]=  LOC_N0; gp->nloc[1] = LOC_N1; 
	gp->nloc[2] = LOC_N2; gp->nloc[3] = LOC_N3;
    
	gp->vol3s[0] = vol4/nd0; gp->vol3s[1] = vol4/nd1;
	gp->vol3s[2] = vol4/nd2; gp->vol3s[3] = vol4/nd3;

	gp->halos[0] = D0_HALO; 
	gp->halos[1] = D1_HALO; 
	gp->halos[2] = D2_HALO; 
	gp->halos[3] = D3_HALO;

	gp->nranks[0] = NRANKS_D0; gp->nranks[1] = NRANKS_D1; 
	gp->nranks[2] = NRANKS_D2; gp->nranks[3] = NRANKS_D3;
 
	// read from input file
	gp->xyztmap[0] = gp->xmap;  gp->xyztmap[1] = gp->ymap;
	gp->xyztmap[2] = gp->zmap;  gp->xyztmap[3] = gp->tmap;

	gp->d0123map[gp->xmap] = 0;   gp->d0123map[gp->ymap] = 1;
	gp->d0123map[gp->zmap] = 2;   gp->d0123map[gp->tmap] = 3;

	// consistency checks


	int expnx =(gp->nd[gp->xmap]-2*gp->halos[gp->xmap])*
		gp->nranks[gp->xmap]; 
	int expny =(gp->nd[gp->ymap]-2*gp->halos[gp->ymap])*
		gp->nranks[gp->ymap]; 
	int expnz =(gp->nd[gp->zmap]-2*gp->halos[gp->zmap])*
		gp->nranks[gp->zmap]; 
	int expnt =(gp->nd[gp->tmap]-2*gp->halos[gp->tmap])*
		gp->nranks[gp->tmap]; 

	if(gp->gnx != expnx || gp->gny != expny ||
		 gp->gnz != expnz  || gp->gnt != expnt ){ 

		if(0==devinfo.myrank){
			printf("Error, input file lattice dimensions are not compatible\n");
			printf("       with the lattice dimensions written in geometry.h.\n");
			printf("       Either modify the input file, or recompile,\n");
			printf("(input) nx=%d\tny=%d\tnz=%d\tnt=%d\n",
						 gp->gnx,gp->gny,gp->gnz,gp->gnt);
			printf("With this mapping, the following lattice is expected:\n");
			printf(" nx %d ny %d nz %d nt %d ",
						 (gp->nd[gp->xmap]-(gp->nranks[gp->xmap]>1?2*HALO_WIDTH:0))*gp->nranks[gp->xmap],
						 (gp->nd[gp->ymap]-(gp->nranks[gp->ymap]>1?2*HALO_WIDTH:0))*gp->nranks[gp->ymap],
						 (gp->nd[gp->zmap]-(gp->nranks[gp->zmap]>1?2*HALO_WIDTH:0))*gp->nranks[gp->zmap],
						 (gp->nd[gp->tmap]-(gp->nranks[gp->tmap]>1?2*HALO_WIDTH:0))*gp->nranks[gp->tmap]);
		}
		res = 1;
	}
	int maps[4] = {gp->xmap,gp->ymap,gp->zmap,gp->tmap};
	int stop = 0;
	int imap,jmap;
	for(imap = 0 ; imap<3; imap++) for(jmap = imap+1 ; jmap<4; jmap++)
																	 stop = stop || (maps[imap] == maps[jmap]);

	if(stop){

		if(0==devinfo.myrank)
			printf("ERROR: found two equal direction mappings (%s:%d)\n",
						 __FILE__,__LINE__);
		res = 1;

	}

	return res;

}


#endif
