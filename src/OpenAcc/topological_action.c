#ifndef TOPOLOGICAL_ACTION_C
#define TOPOLOGICAL_ACTION_C

#include "./topological_action.h"
#include "../Mpi/multidev.h"
#include "../Meas/gauge_meas.h"
#include "./alloc_vars.h"
#include "./action.h"
#ifndef STOUTING_H
#include "./stouting.h"
#endif
#ifdef MULTIDEVICE
#include <mpi.h>
#endif


void load_topo(const char *path,const double barr,const double width, double *grid,const int ngrid)
{

	if(verbosity_lv>4)
		printf("Searching file %s\n",path);
	double widthh = width/2;
	FILE *fin=fopen(path,"r");
	if(fin==NULL)
		{
			printf("ERROR: topological file not found");
			exit(1);
		}
	for(int igrid=0; igrid<=ngrid; igrid++)
		{
			double xread;
			int rc=fscanf(fin,"%lf %lf\n",&xread,&grid[igrid]);
			if(rc!=2)
				{
					printf("ERROR: topological file ill-formatted\n");
					exit(1);
				}
			int jgrid=(int)floor((xread+barr+widthh)/width);
			if(igrid!=jgrid)
				{
					printf("ERROR: found %d (%lf) when expecting %d - topological file ill-formatted",jgrid,xread,igrid);
					exit(1);
				}
		}
	int check_closure=fclose(fin);
	if(check_closure!=0)
		{
			printf("ERROR: file %s can't be closed",path);
			exit(1);
		}

}



double topodynamical_pot(double Q)
{
	double barr=act_params.barrier, width=act_params.width;
	double widthh = width/2;
	double igridfloat = (Q+barr)/width;
	int ngrid= (int) floor((2*barr+widthh)/width);
	int igrid= (int) floor(igridfloat);
	double grid[ngrid];
	if(verbosity_lv>4)
		printf("\t\t\tMPI%02d - load_topo(path,barr,width,grid,ngrid)\n",devinfo.myrank);
	load_topo(act_params.topo_file_path,barr,width,grid,ngrid);
	if(igrid>=0 && igrid<=ngrid)
		{
			// interpolate
			double x0=igrid*width-barr;
			double m=(grid[igrid+1]-grid[igrid])/width;
			double q=grid[igrid]-m*x0;
			printf("topotential: q+mx = %f; q = %f m = %f, x0 = %f, igrid = %d, ngrid = %d, Q = %f\n", q+m*x0,q,m,x0,igrid,ngrid,Q);
			printf("barr = %f, width = %f, widthh = %f, \n",barr,width,widthh);
			return q+m*Q;
		}
	else if(igrid<0){
		printf("topotential: grid[0] = %f\n, igrid = %d, ngrid = %d, Q = %f\n",grid[0],igrid,ngrid,Q);
		printf("barr = %f, width = %f, widthh = %f\n",barr,width,widthh);
		return grid[0];
	}
	else{
		printf("topotential: grid[ngrid] = %f, igrid = %d, ngrid = %d, Q = %f\n",grid[ngrid],igrid,ngrid,Q);
		printf("barr = %f, width = %f, widthh = %f\n",barr,width,widthh);
		return grid[ngrid];
	}
}




double compute_topo_action(__restrict su3_soa * u
#ifdef STOUT_TOPO
													 ,__restrict su3_soa * const tstout_conf_acc_arr
#endif
													 ,int nnp_openacc[sizeh][4][2],int nnm_openacc[sizeh][4][2])
{
  __restrict su3_soa * quadri;
  __restrict su3_soa * conf_to_use;
  __restrict double_soa *loc_q;
  

#ifdef STOUT_TOPO
	#pragma acc update host(tstout_conf_acc_arr[0:8])
  if(act_params.topo_stout_steps > 0){
    stout_wrapper(u,tstout_conf_acc_arr,nnp_openacc,nnm_openacc,1);
    conf_to_use = &(tstout_conf_acc_arr[8*(act_params.topo_stout_steps-1)]);
  }
  else conf_to_use=u;
#else
  conf_to_use=u;
#endif
  
  if(verbosity_lv>4)
    printf("\t\t\tMPI%02d - compute_topological_charge(u,quadri,loc_q, nnp_openacc, nnm_openacc)\n",devinfo.myrank);
  
  posix_memalign((void **)&quadri,128,8*sizeof(su3_soa));
	#pragma acc enter data create(quadri[0:8])
  
  posix_memalign((void **)&loc_q,128,2*sizeof(double_soa));
	#pragma acc enter data create(loc_q[0:2])  

  double Q = compute_topological_charge(conf_to_use, quadri, loc_q, nnp_openacc, nnm_openacc);
  
	#pragma acc exit data delete(quadri)  
  free(quadri);

	#pragma acc exit data delete(loc_q)  
  free(loc_q);
  
  if(verbosity_lv>4)
    printf("Topological Charge: %lf\n", Q);
  
  if(verbosity_lv>4)
    printf("\t\t\tMPI%02d - topodynamical_pot(Q)\n",devinfo.myrank);

  double topo_action = topodynamical_pot(Q);

  return topo_action;
}


#endif
