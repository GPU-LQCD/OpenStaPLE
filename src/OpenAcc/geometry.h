#ifndef GEOMETRY_H_
#define GEOMETRY_H_


// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
#define __restrict
#endif

#include "../Mpi/geometry_multidev.h"

#define nd0h (nd0 >> 1) // nx/2
#define nd1h (nd1 >> 1)
#define nd2h (nd2 >> 1)
#define nd3h (nd3 >> 1)

#define sizehh nd0*nd1*nd2*nd3/2 

#define ANTIPERIODIC_T_BC // else periodic time bc are taken

#define vol1 nd0
#define vol2 (nd1 * vol1)
#define vol3 (nd2 * vol2)
#define vol4 (nd3 * vol3)
#define NSITES vol4
#define size2 (2*NSITES)
#define size3 (3*NSITES)
#define sizeh (NSITES / 2)
#define no_links (4 * vol4)

typedef struct geom_parameters_t{

	int initialized_check;
	int gnx,gny,gnz,gnt; // equal to the dimensions of the configuration
	// map of physical directions onto logical directions
	int xmap,ymap,zmap,tmap;// tmap is the "antiperiodic direction" for fermions

	// The following will be set by set_geom_glv()
	int nd[4];
	int nloc[4];
	int vol3s[4];
	int xyztmap[4];
	int d0123map[4];

	int nranks[4]; // ranks number
	int halos[4];

} geom_parameters;

extern geom_parameters geom_par;

// tables for the nearest neighbours
// nnp[site_half][dir][par] = nearest neighbour in the positive direction "dir"            
//                            starting from the site "site_half" (from 0 to sizeh) of parity "par"
// nnm[site_half][dir][par] = nearest neighbour in the negative direction "dir"
//                            starting from the site "site_half" (from 0 to sizeh) of parity "par"

int set_geom_glv(geom_parameters* gp);
void compute_nnp_and_nnm_openacc(int nnp_openacc[sizeh][4][2],int nnm_openacc[sizeh][4][2]);

#endif
