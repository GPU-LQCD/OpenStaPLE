#ifndef GEOMETRY_MULTIDEV_H_
#define GEOMETRY_MULTIDEV_H_

    // GEOMETRIC DEFINES
#include "../Include/common_defines.h"
#ifdef GAUGE_ACT_TLSM
 #define HALO_WIDTH 2  // for rectangles
#elif defined(GAUGE_ACT_WILSON)
 #define HALO_WIDTH 1
#endif
#define FERMION_HALO 1
#define GAUGE_HALO HALO_WIDTH
//LOCAL lattice dimensions


// lattice dimensions
//#define LOC_N0  // #define these with a -D when invoking compiler
//#define LOC_N1  //  see makefile and generate_makefile.py
//#define LOC_N2
//#define LOC_N3

// MULTIDEVICE
#define NRANKS_D0 1   // Keep 1 - only "salamino" allowed
#define NRANKS_D1 1   // Keep 1 - only "salamino" allowed
#define NRANKS_D2 1   // Keep 1 - only "salamino" allowed
//#define NRANKS_D3 // #define this with a -D when invoking compiler, 
                    // see makefile and generate_makefile.py


///HALO STRUCTURE (border widths)
// AUTOMATIC DEFINITION OF HALOS
#if NRANKS_D0 == 1 
#define D0_HALO 0
#define D0_FERMION_HALO 0
#else 
#define D0_HALO HALO_WIDTH
#define D0_FERMION_HALO FERMION_HALO
#define MULTIDEVICE
#endif

#if NRANKS_D1 == 1 
#define D1_HALO 0
#define D1_FERMION_HALO 0
#else 
#define D1_HALO HALO_WIDTH
#define D1_FERMION_HALO FERMION_HALO
#define MULTIDEVICE
#endif

#if NRANKS_D2 == 1 
#define D2_HALO 0
#define D2_FERMION_HALO 0
#else 
#define D2_HALO HALO_WIDTH
#define D2_FERMION_HALO FERMION_HALO
#define MULTIDEVICE
#endif

#if NRANKS_D3 == 1 
#define D3_HALO 0
#define D3_FERMION_HALO 0
#else 
#define D3_HALO HALO_WIDTH
#define D3_FERMION_HALO FERMION_HALO
#define MULTIDEVICE
#endif

#ifdef NREPLICAS 
#if NREPLICAS > 1
#define MULTIDEVICE
#endif
#endif

// Local aNd Halo lattice dimensions 
// Hopefully automatically calculated by the compiler during optimization
//AUTOMATIC DEFINITIONS OF LNH STUFF

#define LNH_N0 (LOC_N0+2*D0_HALO)
#define LNH_N1 (LOC_N1+2*D1_HALO)
#define LNH_N2 (LOC_N2+2*D2_HALO)
#define LNH_N3 (LOC_N3+2*D3_HALO)

#define nd0 LNH_N0
#define nd1 LNH_N1
#define nd2 LNH_N2
#define nd3 LNH_N3

// GLOBAL quantities
//AUTOMATIC DEFINITIONS OF GLOBAL STUFF
#define GL_N0 (LOC_N0*NRANKS_D0)
#define GL_N1 (LOC_N1*NRANKS_D1)
#define GL_N2 (LOC_N2*NRANKS_D2)
#define GL_N3 (LOC_N3*NRANKS_D3)

#define GL_VOL1  GL_N0
#define GL_VOL2 (GL_N1 * GL_VOL1)
#define GL_VOL3 (GL_N2 * GL_VOL2)
#define GL_VOL4 (GL_N3 * GL_VOL3)
#define GL_N0H ( GL_N0 >> 1) // nx/2
#define GL_N1H ( GL_N1 >> 1)
#define GL_N2H ( GL_N2 >> 1)
#define GL_N3H ( GL_N3 >> 1)

#define GL_SIZE GL_VOL4
#define GL_SIZE2 (2*GL_SIZE)
#define GL_SIZE3 (3*GL_SIZE)
#define GL_SIZEH (GL_SIZE / 2)
#define GL_NO_LINKS (4 * GL_VOL4)

// LOCAL STUFF
#define LOC_VOL1  LOC_N0
#define LOC_VOL2 (LOC_N1 * LOC_VOL1)
#define LOC_VOL3 (LOC_N2 * LOC_VOL2)
#define LOC_VOL4 (LOC_N3 * LOC_VOL3)
#define LOC_N0H ( LOC_N0 >> 1) // nx/2
#define LOC_N1H ( LOC_N1 >> 1)
#define LOC_N2H ( LOC_N2 >> 1)
#define LOC_N3H ( LOC_N3 >> 1)

#define LOC_SIZE LOC_VOL4
#define LOC_SIZE2 (2*LOC_SIZE)
#define LOC_SIZE3 (3*LOC_SIZE)
#define LOC_SIZEH ((int)(LOC_SIZE / 2) + LOC_SIZE % 2)
#define LOC_SIZEH_A ((int)(LOC_SIZE / 2) + LOC_SIZE % 2) //no of sites having the same parity as the local origin
#define LOC_SIZEH_B ((int)(LOC_SIZE / 2)) // no of sites having opposite parity to the local origin
#define LOC_NO_LINKS (4 * LOC_VOL4)

#ifdef MULTIDEVICE

#define NRANKS NRANKS_D0*NRANKS_D1*NRANKS_D2*NRANKS_D3


#define LNH_VOL1  LNH_N0
#define LNH_VOL2 (LNH_N1 * LNH_VOL1)
#define LNH_VOL3 (LNH_N2 * LNH_VOL2)
#define LNH_VOL4 (LNH_N3 * LNH_VOL3)
#define LNH_N0H ( LNH_N0 >> 1) // nx/2
#define LNH_N1H ( LNH_N1 >> 1)
#define LNH_N2H ( LNH_N2 >> 1)
#define LNH_N3H ( LNH_N3 >> 1)

#define LNH_SIZE LNH_VOL4
#define LNH_SIZE2 (2*LNH_SIZE)
#define LNH_SIZE3 (3*LNH_SIZE)
#define LNH_SIZEH ((int)(LNH_SIZE / 2)+ LNH_SIZE % 2 )
#define LNH_SIZEH_A ((int)(LNH_SIZE / 2)+ LNH_SIZE % 2 ) //no of sites having the same parity as the local origin
#define LNH_SIZEH_B ((int)(LNH_SIZE / 2)) // no of sites having opposite parity to the local origin
#define LNH_NO_LINKS (4 * LNH_VOL4)

//SURFACE STRUCTURE (surface widths)
//AUTOMATIC DEFINITION OF BULK STUFF
#define D0_SURF D0_HALO 
#define D1_SURF D1_HALO
#define D2_SURF D2_HALO
#define D3_SURF D3_HALO

//Local Bulk lattice dimensions

#define LB_N0 (LOC_N0-2*D0_SURF)
#define LB_N1 (LOC_N1-2*D1_SURF)
#define LB_N2 (LOC_N2-2*D2_SURF)
#define LB_N3 (LOC_N3-2*D3_SURF)

#define LB_VOL1  LB_N0
#define LB_VOL2 (LB_N1 * LB_VOL1)
#define LB_VOL3 (LB_N2 * LB_VOL2)
#define LB_VOL4 (LB_N3 * LB_VOL3)
#define LB_N0H ( LB_N0 >> 1) // nx/2
#define LB_N1H ( LB_N1 >> 1)
#define LB_N2H ( LB_N2 >> 1)
#define LB_N3H ( LB_N3 >> 1)

#define LB_SIZE LB_VOL4
#define LB_SIZE2 (2*LB_SIZE)
#define LB_SIZE3 (3*LB_SIZE)
#define LB_SIZEH (LB_SIZE / 2)
#define LB_NO_LINKS (4 * LB_VOL4)


struct dev_info;

#else

#define LNH_VOL1      LOC_VOL1
#define LNH_VOL2      LOC_VOL2
#define LNH_VOL3      LOC_VOL3
#define LNH_VOL4      LOC_VOL4
#define LNH_N0H       LOC_N0H 
#define LNH_N1H       LOC_N1H 
#define LNH_N2H       LOC_N2H 
#define LNH_N3H       LOC_N3H 
                                  
#define LNH_SIZE      LOC_SIZE 
#define LNH_SIZE2     LOC_SIZE2
#define LNH_SIZE3     LOC_SIZE3
#define LNH_SIZEH     LOC_SIZEH
#define LNH_NO_LINKS  LOC_NO_LINKS
    
                       
#endif                 
typedef struct vec4int_t{
int d0,d1,d2,d3;

} vec4int;



#pragma acc routine seq
static inline int gl_to_gl_snum(int gl_0, int gl_1, int gl_2, int gl_3){
//global coordinates to global 'snum' index

    int ris = gl_0 + gl_1 * GL_VOL1 + gl_2 * GL_VOL2 + gl_3 * GL_VOL3 ;
    return ris/2;// <---  /2 Pay attention to even/odd  (see init_geo) 

}

#pragma acc routine seq
//accesso al sito del reticolo. Lattice site access function.
static inline int snum_acc(int lnh_0, int lnh_1, int lnh_2, int lnh_3){
// local'n'halo coordinates to lnh 'snum' index
 
    //NOTE: Since lnh_0 goes from 0 to LNH_N0-1,
    //      here loc_0 goes from -D0_HALO to LOC_N0+D0_HALO-1

    int ris = lnh_0 + lnh_1 * LNH_VOL1 + lnh_2 * LNH_VOL2 + lnh_3 * LNH_VOL3;

    return ris/2;// <---  /2 Pay attention to even/odd  (see init_geo) 
//obviously with this uncast division we get snum_acc(2n,k1,k2,k4)=snum(2n+1,k1,k2,k3)
}

#pragma acc routine seq
static inline int loc_to_lnh_snum(int loc_0, int loc_1, int loc_2, int loc_3){
//local coordinates to loc'n'halo 'snum'index
    // Actually, the real memory layout is of the LNH type.
    // So, at the moment I don't see any real use for this function.
    
    loc_0+=D0_HALO;
    loc_1+=D1_HALO;
    loc_2+=D2_HALO;
    loc_3+=D3_HALO;

    return snum_acc(loc_0,loc_1,loc_2,loc_3);
    // ^^ /2 Pay attention to even/odd  (see init_geo) 
}

#pragma acc routine seq
static inline int lnh_to_gl_snum(int lnh_0, int lnh_1, int lnh_2, int lnh_3, vec4int myrank4int){

    lnh_0 += LOC_N0 * myrank4int.d0 - D0_HALO; // to global ref frame
    lnh_1 += LOC_N1 * myrank4int.d1 - D1_HALO;
    lnh_2 += LOC_N2 * myrank4int.d2 - D2_HALO;
    lnh_3 += LOC_N3 * myrank4int.d3 - D3_HALO;

    // in case lnh_* is out of bonds
    lnh_0 = lnh_0 % GL_N0; lnh_0-= GL_N0 * ( lnh_0 >> 31 );
    lnh_1 = lnh_1 % GL_N1; lnh_1-= GL_N1 * ( lnh_1 >> 31 );
    lnh_2 = lnh_2 % GL_N2; lnh_2-= GL_N2 * ( lnh_2 >> 31 );
    lnh_3 = lnh_3 % GL_N3; lnh_3-= GL_N3 * ( lnh_3 >> 31 );

    return gl_to_gl_snum(lnh_0,lnh_1,lnh_2,lnh_3);
    // pay attention to even/odd

}

#pragma acc routine seq
static inline int target_lnh_to_gl_snum(int lnh_0, int lnh_1, int lnh_2,
        int lnh_3,
        vec4int target_gl_loc_origin4int ){

    lnh_0 += target_gl_loc_origin4int.d0 - D0_HALO;
    lnh_1 += target_gl_loc_origin4int.d1 - D1_HALO;
    lnh_2 += target_gl_loc_origin4int.d2 - D2_HALO;
    lnh_3 += target_gl_loc_origin4int.d3 - D3_HALO;

    lnh_0 = lnh_0 % GL_N0; lnh_0-= GL_N0 * ( lnh_0 >> 31 );
    lnh_1 = lnh_1 % GL_N1; lnh_1-= GL_N1 * ( lnh_1 >> 31 );
    lnh_2 = lnh_2 % GL_N2; lnh_2-= GL_N2 * ( lnh_2 >> 31 );
    lnh_3 = lnh_3 % GL_N3; lnh_3-= GL_N3 * ( lnh_3 >> 31 );

    return gl_to_gl_snum(lnh_0,lnh_1,lnh_2,lnh_3);
    // pay attention to even/odd

}
/***************************************************
links used according to this scheme

0           LNH_SIZE       2*LNH_SIZE      3*LNH_SIZE       LNH_NO_LINKS
|-------|-------|-------|-------|-------|-------|-------|-------|
    e       o        e      o       e       o       e       o
      0-dir           1-dir           2-dir           3-dir
NOTE:the number of even and odd sites on the local lattice can be 
different for the local sublattice if the local sublattice dimensions are all odd.
initialize geometry
periodic spatial bc are always assumed (this is relevant only if
the considered direction is not paralelized))
*/

#pragma acc routine seq
static inline vec4int xyzt_rank(int rank){

    vec4int rank4int;

    rank4int.d0 = rank % NRANKS_D0;
    rank4int.d1 =( rank / NRANKS_D0) % NRANKS_D1;
    rank4int.d2 =( rank / (NRANKS_D0*NRANKS_D1)) % NRANKS_D2;
    rank4int.d3 =( rank / (NRANKS_D0*NRANKS_D1*NRANKS_D2)) % NRANKS_D3;

    return rank4int;

}
#pragma acc routine seq
static inline int rank_from_0123_rank(int rank_0, int rank_1, 
                               int rank_2, int rank_3){

    return rank_0 + ( NRANKS_D0 * (rank_1 + NRANKS_D1 * (rank_2  + NRANKS_D2 * rank_3)));

}

#pragma acc routine seq
static inline vec4int gl_loc_origin_from_rank(int rank){
    vec4int res, rank4int = xyzt_rank(rank);
    res.d0 = rank4int.d0 * LOC_N0;
    res.d1 = rank4int.d1 * LOC_N1;
    res.d2 = rank4int.d2 * LOC_N2;
    res.d3 = rank4int.d3 * LOC_N3;
    return res;

}




#endif
