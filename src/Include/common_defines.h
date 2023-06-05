#ifndef COMMON_DEFINES_H_
#define COMMON_DEFINES_H_

// #define BACKFIELD
// #define IMCHEMPOT

// if BACKFIELD or IMCHEMPOT are defined, then the routine which multiplies by
// the suitable phases is used in the Dirac matrix application. Otherwise, it
// uses the other routine (in an hard coded way), which multiplies links by
// fermion and no more.
// This is what PHASE_MAT_VEC_MULT does.

#if defined(BACKFIELD) || defined (IMCHEMPOT)
#define PHASE_MAT_VEC_MULT 
#endif

// the following definitions are necessary with some compilers, for example the
// one in the nvidia hpc-sdk package.
// When using the PGI compiler in M100, they give errors. Keep them commented
// in that case
/*
	#define conj __builtin_conj
	#define creal __builtin_creal
	#define cimag __builtin_cimag
	#define conjf __builtin_conjf
	#define crealf __builtin_crealf
	#define cimagf __builtin_cimagf
*/

#define DIM_BLOCK_X 8 // This should divide (nx/2)
#define DIM_BLOCK_Y 8 // This should divide ny
#define DIM_BLOCK_Z 8 // This should divide nz*nt

// #define TIMING_ALL // if defined many computation times are printed in the output

// WARNING!!!!!! comment the one you don't wanna use and decoment the other one
// ACHTUNG! Please, use script_compile_Open_StaPLE.sh to comment/uncomment these macros

#define GAUGE_ACT_TLSM
//#define GAUGE_ACT_WILSON

// WARNING! uncomment only if you want to perform parallel tempering
// ACHTUNG! Please, use script_compile_Open_StaPLE.sh to comment/uncomment this macro

#ifdef NREPLICAS 
#if NREPLICAS > 1
#define PAR_TEMP
#endif
#endif

#define STOUT_FERMIONS
#ifdef STOUT_FERMIONS
#define STOUT_STEPS 2
#endif

#define STOUT_TOPO

#define ONE_BY_THREE   0.33333333333333333333333
#define ONE_BY_THREEF   0.33333333333333333333333f
#define ONE_BY_SIX     0.16666666666666666666666
#define ONE_BY_SIXF     0.16666666666666666666666f


#ifdef GAUGE_ACT_TLSM
#define GAUGE_ACTION   1     // 0 --> Wilson; 1 --> tree level Symanzik
#define C_ZERO         (5.0*ONE_BY_THREE)
#define C_ZEROF         (5.0f*ONE_BY_THREEF)
#define C_ONE          (-0.25*ONE_BY_THREE)
#define C_ONEF          (-0.25f*ONE_BY_THREEF)
#define PLAQ_EXTENT    2
#endif

#ifdef GAUGE_ACT_WILSON
#define GAUGE_ACTION   0     // 0 --> Wilson; 1 --> tree level Symanzik
#define C_ZERO         1.0
#define C_ZEROF         1.0f
#define C_ONE          0.0
#define C_ONEF          0.0f
#define PLAQ_EXTENT    1
#endif


extern int verbosity_lv;


// macro to check the outcome of 'fscanf', reporting the file and line where the problem
// happened.
 
#define CHECKREAD(expr,should_read)																			\
	{int read = expr ;if(read != should_read) {														\
			printf("%s:%d, Error, not read expected number of entries : %d read vs %d should_read\n.", __FILE__, __LINE__ , read,should_read); \
			exit(1);}}																												\


#endif
