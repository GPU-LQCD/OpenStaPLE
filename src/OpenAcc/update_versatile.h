#ifndef UPDATE_VERSATILE_H_
#define UPDATE_VERSATILE_H_

#include "./struct_c_def.h"

// metro==1 --> then metropolis test is performed
// metro==0 --> thermalization


int UPDATE_SOLOACC_UNOSTEP_VERSATILE(su3_soa *tconf_acc,
																		 double res_metro, double res_md, int id_iter,int acc,int metro,
																		 int nnp_openacc[sizeh][4][2], int nnm_openacc[sizeh][4][2], int max_cg);


#endif
