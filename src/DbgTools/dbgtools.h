#ifndef DBGTOOLS_H_
#define DBGTOOLS_H_

#include "../OpenAcc/struct_c_def.h"

extern int ipdot_f_reset; // flag to be set to 1 when md starts, 
                   // and to be set to zero when first force calculation happens
                   // if it is 1 force difference is not calculated.
extern int ipdot_g_reset; // same


// multi rank data structure read/write functions - for debug/testing
// multi rank wrappers

// gl3 / su3
void print_gl3_soa_wrapper(const su3_soa * gl3, const char* nomefile);
int read_gl3_soa_wrapper(su3_soa * gl3, const char* nomefile);

// fermions
void print_vec3_soa_wrapper(const vec3_soa * fermion, const char* nomefile);
int read_vec3_soa_wrapper(vec3_soa * fermion, const char* nomefile);

// tamat
void print_tamat_soa_wrapper(const tamat_soa * tamat, const char* nomefile);
int read_tamat_soa_wrapper(tamat_soa * tamat, const char* nomefile);

// thmat
void print_thmat_soa_wrapper(const thmat_soa * thmat, const char* nomefile);
int read_thmat_soa_wrapper(thmat_soa * thmat, const char* nomefile);

// dcomplex
void print_dcomplex_soa_wrapper(const dcomplex_soa *dcarr,
        const char* nomefile);
int read_dcomplex_soa_wrapper(dcomplex_soa * dcarr, const char* nomefile);


// double
void print_double_soa_wrapper(const double_soa * darr, const char* nomefile);
int read_double_soa_wrapper(double_soa * darr, const char* nomefile);


// single rank functions
void print_vec3_soa(vec3_soa * const fermion, const char* nomefile);
int read_vec3_soa(vec3_soa * fermion, const char* nomefile);
void print_su3_soa(su3_soa * const conf, const char* nomefile,int conf_id_iter);
void print_1su3_soa(su3_soa * const conf, const char* nomefile);
void read_su3_soa(su3_soa * conf, const char* nomefile,int * conf_id_iter );
void print_tamat_soa(tamat_soa * const ipdot, const char* nomefile);
int read_tamat_soa(tamat_soa * ipdot, const char* nomefile);
void print_thmat_soa(thmat_soa * const ipdot, const char* nomefile);
void print_1thmat_soa(thmat_soa * const ipdot, const char* nomefile);
int read_thmat_soa(thmat_soa * ipdot, const char* nomefile);
void print_double_soa(double_soa * const backfield, const char* nomefile);
void print_1double_soa(double_soa * const vettore, const char* nomefile);
void dbg_print_su3_soa(su3_soa * const conf, const char* nomefile,int conf_id_iter);
int dbgread_su3_soa(su3_soa * conf, const char* nomefile,int * conf_id_iter );
void calc_loc_abelian_plaquettes(const double_soa* phases,
        double_soa * loc_abelian_plaquettes, const int mu, const int nu );

void print_all_abelian_plaquettes(const double_soa* phases, const char * filename);


void dbgprint_gl3_soa(su3_soa * const conf, const char* nomefile,int conf_id_iter);
int dbgread_gl3_soa(su3_soa * conf, const char* nomefile,int * conf_id_iter );

#define STAMPA_DEBUG_SU3_SOA(var,dir,idx)					\
  printf("%s[%d], idx %d:\n(%le,%le)    (%le,%le)    (%le,%le)\n(%le,%le)    (%le,%le)    (%le,%le)\n(%le,%le)    (%le,%le)    (%le,%le)\n\n", #var, dir, idx, \
	   creal(var[dir].r0.c0[idx]),cimag(var[dir].r0.c0[idx]),creal(var[dir].r0.c1[idx]),cimag(var[dir].r0.c1[idx]),creal(var[dir].r0.c2[idx]),cimag(var[dir].r0.c2[idx]), \
	   creal(var[dir].r1.c0[idx]),cimag(var[dir].r1.c0[idx]),creal(var[dir].r1.c1[idx]),cimag(var[dir].r1.c1[idx]),creal(var[dir].r1.c2[idx]),cimag(var[dir].r1.c2[idx]), \
	   creal(var[dir].r2.c0[idx]),cimag(var[dir].r2.c0[idx]),creal(var[dir].r2.c1[idx]),cimag(var[dir].r2.c1[idx]),creal(var[dir].r2.c2[idx]),cimag(var[dir].r2.c2[idx]));

#define STAMPA_DEBUG_SU3_SOA_DAG(var,dir,idx)					\
  printf("%s[%d], idx %d:\n(%le,%le)    (%le,%le)    (%le,%le)\n(%le,%le)    (%le,%le)    (%le,%le)\n(%le,%le)    (%le,%le)    (%le,%le)\n\n", #var, dir, idx, \
	   creal(var[dir].r0.c0[idx]),-cimag(var[dir].r0.c0[idx]),creal(var[dir].r1.c0[idx]),-cimag(var[dir].r1.c0[idx]),creal(var[dir].r2.c0[idx]),-cimag(var[dir].r2.c0[idx]), \
	   creal(var[dir].r0.c1[idx]),-cimag(var[dir].r0.c1[idx]),creal(var[dir].r1.c1[idx]),-cimag(var[dir].r1.c1[idx]),creal(var[dir].r2.c1[idx]),-cimag(var[dir].r2.c1[idx]), \
	   creal(var[dir].r0.c2[idx]),-cimag(var[dir].r0.c2[idx]),creal(var[dir].r1.c2[idx]),-cimag(var[dir].r1.c2[idx]),creal(var[dir].r2.c2[idx]),-cimag(var[dir].r2.c2[idx]));

#define STAMPA_DEBUG_TAMAT_SOA(var,dir,idx)				\
  printf("%s[%d], idx %d :\n(%le,%le)\t(%le,%le)\t(%le,%le)\n(%le,%le)\t(%le,%le)\t(%le,%le)\n(%le,%le)\t(%le,%le)\t(%le,%le)\n\n", #var, dir, idx, \
	 0.0,var[dir].ic00[idx],creal(var[dir].c01[idx]),cimag(var[dir].c01[idx]),creal(var[dir].c02[idx]),cimag(var[dir].c02[idx]), \
	 -creal(var[dir].c01[idx]),cimag(var[dir].c01[idx]),0.0,var[dir].ic11[idx],creal(var[dir].c12[idx]),cimag(var[dir].c12[idx]), \
	 -creal(var[dir].c02[idx]),cimag(var[dir].c02[idx]),-creal(var[dir].c12[idx]),cimag(var[dir].c12[idx]),0.0,-(var[dir].ic00[idx]+var[dir].ic11[idx]));

#define STAMPA_DEBUG_SINGLE_SU3(var)					\
  printf("%s:\n(%le,%le)    (%le,%le)    (%le,%le)\n(%le,%le)    (%le,%le)    (%le,%le)\n(%le,%le)    (%le,%le)    (%le,%le)\n\n", #var, \
	 creal(var.comp[0][0]),cimag(var.comp[0][0]),creal(var.comp[0][1]),cimag(var.comp[0][1]),creal(var.comp[0][2]),cimag(var.comp[0][2]), \
	 creal(var.comp[1][0]),cimag(var.comp[1][0]),creal(var.comp[1][1]),cimag(var.comp[1][1]),creal(var.comp[1][2]),cimag(var.comp[1][2]), \
	 creal(var.comp[2][0]),cimag(var.comp[2][0]),creal(var.comp[2][1]),cimag(var.comp[2][1]),creal(var.comp[2][2]),cimag(var.comp[2][2]));

#define STAMPA_DEBUG_SINGLE_TAMAT(var)					\
  printf("%s:\n(%le,%le)\t(%le,%le)\t(%le,%le)\n(%le,%le)\t(%le,%le)\t(%le,%le)\n(%le,%le)\t(%le,%le)\t(%le,%le)\n\n", #var, \
	 0.0,var.ic00,creal(var.c01),cimag(var.c01),creal(var.c02),cimag(var.c02), \
	 -creal(var.c01),cimag(var.c01),0.0,var.ic11,creal(var.c12),cimag(var.c12), \
	 -creal(var.c02),cimag(var.c02),-creal(var.c12),cimag(var.c12),0.0,-(var.ic00+var.ic11));


#endif
