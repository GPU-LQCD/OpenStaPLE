#ifndef POLYAKOV_C_
#define POLYAKOV_C_

#ifdef __GNUC__
#define _POSIX_C_SOURCE 200809L // not to have warning on posix memalign
#endif

#include <stdlib.h>
#include <stdio.h>
#include "../OpenAcc/geometry.h"
#include "../OpenAcc/su3_utilities.h"
#include "./polyakov.h"
#include "../OpenAcc/su3_measurements.h"


#ifdef MULTIDEVICE
#include <mpi.h>
#endif
#include "../Mpi/multidev.h"


// For Polyakov loop calculations
#define ALIGN 128
#ifdef NRANKS_D3 > 1
#define vol30h LOC_VOL4/LOC_N0/2 
#define vol31h LOC_VOL4/LOC_N1/2 
#define vol32h LOC_VOL4/LOC_N2/2 
#define vol33h LOC_VOL4/LOC_N3/2 
#else
#define vol30h vol4/nd0/2
#define vol31h vol4/nd1/2
#define vol32h vol4/nd2/2
#define vol33h vol4/nd3/2
#endif


#define ALLOCCHECK(control_int,var)  if(control_int != 0 ) \
                                                       printf("\tError in  allocation of %s . \n", #var);\
else if(verbosity_lv > 2) printf("\tAllocation of %s : OK , %p\n", #var, var );\


typedef struct vec3_plk0_t {
    d_complex c0[vol30h];    d_complex c1[vol30h];    d_complex c2[vol30h];
} vec3_plk0;
typedef struct su3_plk0_t {
    vec3_plk0 r0;  vec3_plk0 r1;
} su3_plk0;


typedef struct vec3_plk1_t {
    d_complex c0[vol31h];    d_complex c1[vol31h];    d_complex c2[vol31h];
} vec3_plk1;
typedef struct su3_plk1_t {
    vec3_plk1 r0;  vec3_plk1 r1;
} su3_plk1;


typedef struct vec3_plk2_t {
    d_complex c0[vol32h];     d_complex c1[vol32h];     d_complex c2[vol32h];
} vec3_plk2;
typedef struct su3_plk2_t {
    vec3_plk2 r0;  vec3_plk2 r1;
} su3_plk2;

typedef struct vec3_plk3_t {
    d_complex c0[vol33h];     d_complex c1[vol33h];     d_complex c2[vol33h];
} vec3_plk3;
typedef struct su3_plk3_t {
    vec3_plk3 r0;  vec3_plk3 r1;
} su3_plk3;


#define DEF_SU3_PLKN_SU3_SOA_MULTINPLACE(FUNCNAME, TYPE)\
    static inline void    FUNCNAME ( __restrict TYPE * const mat1,\
            const int idx_mat1,\
            __restrict const su3_soa * const mat2,\
            const int idx_mat2) {\
        d_complex mat1_00 = mat1->r0.c0[idx_mat1];\
        d_complex mat1_01 = mat1->r0.c1[idx_mat1];\
        d_complex mat1_02 = mat1->r0.c2[idx_mat1];\
        \
        d_complex mat1_10 = mat1->r1.c0[idx_mat1];\
        d_complex mat1_11 = mat1->r1.c1[idx_mat1];\
        d_complex mat1_12 = mat1->r1.c2[idx_mat1];\
        \
        d_complex mat2_00 = mat2->r0.c0[idx_mat2];\
        d_complex mat2_01 = mat2->r0.c1[idx_mat2];\
        d_complex mat2_02 = mat2->r0.c2[idx_mat2];\
        \
        d_complex mat2_10 = mat2->r1.c0[idx_mat2];\
        d_complex mat2_11 = mat2->r1.c1[idx_mat2];\
        d_complex mat2_12 = mat2->r1.c2[idx_mat2];\
        \
        d_complex mat2_20 = conj( ( mat2_01 * mat2_12 ) - ( mat2_02 * mat2_11) ) ;\
        d_complex mat2_21 = conj( ( mat2_02 * mat2_10 ) - ( mat2_00 * mat2_12) ) ;\
        d_complex mat2_22 = conj( ( mat2_00 * mat2_11 ) - ( mat2_01 * mat2_10) ) ;\
        \
        mat1->r0.c0[idx_mat1] = mat1_00 * mat2_00 + mat1_01 * mat2_10 + mat1_02 * mat2_20 ;\
        mat1->r0.c1[idx_mat1] = mat1_00 * mat2_01 + mat1_01 * mat2_11 + mat1_02 * mat2_21 ;\
        mat1->r0.c2[idx_mat1] = mat1_00 * mat2_02 + mat1_01 * mat2_12 + mat1_02 * mat2_22 ;\
        \
        mat1->r1.c0[idx_mat1] = mat1_10 * mat2_00 + mat1_11 * mat2_10 + mat1_12 * mat2_20 ;\
        mat1->r1.c1[idx_mat1] = mat1_10 * mat2_01 + mat1_11 * mat2_11 + mat1_12 * mat2_21 ;\
        mat1->r1.c2[idx_mat1] = mat1_10 * mat2_02 + mat1_11 * mat2_12 + mat1_12 * mat2_22 ;\
    }


#pragma acc routine seq
DEF_SU3_PLKN_SU3_SOA_MULTINPLACE(su3_plk0_su3_soa_multinplace, su3_plk0)
#pragma acc routine seq
DEF_SU3_PLKN_SU3_SOA_MULTINPLACE(su3_plk1_su3_soa_multinplace, su3_plk1)
#pragma acc routine seq
DEF_SU3_PLKN_SU3_SOA_MULTINPLACE(su3_plk2_su3_soa_multinplace, su3_plk2)
#pragma acc routine seq
DEF_SU3_PLKN_SU3_SOA_MULTINPLACE(su3_plk3_su3_soa_multinplace, su3_plk3)

#define DEF_SU3_PLKN_TRACE(FUNCNAME,TYPE)\
        static inline d_complex FUNCNAME (__restrict const TYPE * const loc_plaq, const int idx)\
{\
    d_complex loc_plaq_00 = loc_plaq->r0.c0[idx];\
    d_complex loc_plaq_01 = loc_plaq->r0.c1[idx];\
    d_complex loc_plaq_10 = loc_plaq->r1.c0[idx];\
    d_complex loc_plaq_11 = loc_plaq->r1.c1[idx];\
    d_complex loc_plaq_22 =  conj( ( loc_plaq_00 * loc_plaq_11 ) - ( loc_plaq_01 * loc_plaq_10) ) ;\
    return (loc_plaq_00 + loc_plaq_11 + loc_plaq_22);\
}

#pragma acc routine seq
DEF_SU3_PLKN_TRACE(su3_plk0_trace,su3_plk0)
#pragma acc routine seq
DEF_SU3_PLKN_TRACE(su3_plk1_trace,su3_plk1)
#pragma acc routine seq
DEF_SU3_PLKN_TRACE(su3_plk2_trace,su3_plk2)
#pragma acc routine seq
DEF_SU3_PLKN_TRACE(su3_plk3_trace,su3_plk3)


d_complex polyakov_loop0(__restrict const su3_soa * const u)
{

    if(geom_par.nranks[0] !=1 ){
        printf("MPI%02d, POLYAKOV LOOP CALCULATION NOT IMPLEMENTED!\n", devinfo.myrank);
        return 0+0*I;
    }
    su3_plk0 *loopplk0;
    double rel,iml;
    int allocation_check;
    allocation_check = posix_memalign((void **)&loopplk0, ALIGN, 2*sizeof(su3_plk0));
    ALLOCCHECK(allocation_check,loopplk0);

#pragma acc data  create(loopplk0[0:2])
    {
#pragma acc kernels present(loopplk0) 
#pragma acc loop independent
        for(int i=0;i<vol30h;i++){
            loopplk0[0].r0.c0[i] = 1;loopplk0[0].r0.c1[i] = 0;loopplk0[0].r0.c2[i] = 0;
            loopplk0[0].r1.c0[i] = 0;loopplk0[0].r1.c1[i] = 1;loopplk0[0].r1.c2[i] = 0;

            loopplk0[1].r0.c0[i] = 1;loopplk0[1].r0.c1[i] = 0;loopplk0[1].r0.c2[i] = 0;
            loopplk0[1].r1.c0[i] = 0;loopplk0[1].r1.c1[i] = 1;loopplk0[1].r1.c2[i] = 0;
        }
        int d0, d1, d2, d3,h,idxh,parity;
        double r = 0;
#pragma acc kernels present(u) present(loopplk0)
#pragma acc loop gang // NOT INDEPENDENT
        for(d0=0; d0 < nd0; d0++) {	    
#pragma acc loop independent gang vector //gang(nd2/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
            for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent gang vector //gang(nd1/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
                for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent vector //vector(DIM_BLOCK_X)
                    for(d1=0; d1<nd1; d1++) {
                        int ld0,ld1,ld2,ld3;
                        ld0 = d0 - D0_HALO;
                        ld1 = d1 - D1_HALO;
                        ld2 = d2 - D2_HALO;
                        ld3 = d3 - D3_HALO;

                        int parity0 = ( 0+ld1+ld2+ld3) % 2;
                        int id0h = (ld1+LOC_N1*(ld2+LOC_N2*ld3))/2;

                        parity = (d0+d1+d2+d3) % 2;
                        idxh = snum_acc(d0,d1,d2,d3);  	

                        su3_plk0_su3_soa_multinplace(&loopplk0[parity0],id0h,
                                &u[0*2+parity],idxh);

                    }//end d0
                }//end d1
            }//end d2
        }//end d3

        rel=0,iml=0;
#pragma acc kernels present(loopplk0) 
#pragma acc loop independent reduction(+:rel) reduction(+:iml)
        for(int i=0;i<vol30h;i++){
            d_complex r = su3_plk0_trace(&loopplk0[0],i);
            r += su3_plk0_trace(&loopplk0[1],i);
            rel += creal(r);
            iml += cimag(r);
        }

    }
    free(loopplk0);

    double trr,tri;
#if NRANKS_D3 > 1
     MPI_Allreduce((void*)&rel,(void*)&trr,
             1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);
     MPI_Allreduce((void*)&iml,(void*)&tri,
             1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);
#else
     trr = rel;
     tri = iml;
#endif 

    return (trr + tri * I)/(vol30h*2*3)/devinfo.nranks;

}

d_complex polyakov_loop1(__restrict const su3_soa * const u)
{

    if(geom_par.nranks[1] !=1 ){
        printf("MPI%02d, POLYAKOV LOOP CALCULATION NOT IMPLEMENTED!\n", devinfo.myrank);
        return 0+0*I;
    }

    su3_plk1 *loopplk1;
    double rel,iml;
    int allocation_check =  posix_memalign((void **)&loopplk1, ALIGN, 2*sizeof(su3_plk1));
    ALLOCCHECK(allocation_check,loopplk1);

#pragma acc data  create(loopplk1[0:2])
    {
#pragma acc kernels present(loopplk1) 
#pragma acc loop independent
        for(int i=0;i<vol31h;i++){
            loopplk1[0].r0.c0[i] = 1;loopplk1[0].r0.c1[i] = 0;loopplk1[0].r0.c2[i] = 0;
            loopplk1[0].r1.c0[i] = 0;loopplk1[0].r1.c1[i] = 1;loopplk1[0].r1.c2[i] = 0;

            loopplk1[1].r0.c0[i] = 1;loopplk1[1].r0.c1[i] = 0;loopplk1[1].r0.c2[i] = 0;
            loopplk1[1].r1.c0[i] = 0;loopplk1[1].r1.c1[i] = 1;loopplk1[1].r1.c2[i] = 0;
        }
        int d0, d1, d2, d3,h,idxh,parity;

#pragma acc kernels present(u) present(loopplk1)
#pragma acc loop gang // NOT INDEPENDENT 
        for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent gang vector //gang(nd2/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
            for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent gang vector //gang(nd1/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
                for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent vector //vector(DIM_BLOCK_X)
                    for(d0=0; d0 < nd0; d0++) {	    
                        int ld0,ld1,ld2,ld3;
                        ld0 = d0 - D0_HALO;
                        ld1 = d1 - D1_HALO;
                        ld2 = d2 - D2_HALO;
                        ld3 = d3 - D3_HALO;

                        int parity1 = (ld0+ 0+ld2+ld3) % 2;
                        int id1h = (ld0+LOC_N0*(ld2+LOC_N2*ld3))/2;


                        parity = (d0+d1+d2+d3) % 2;
                        idxh = snum_acc(d0,d1,d2,d3);  	

                        su3_plk1_su3_soa_multinplace(&loopplk1[parity1],id1h,
                                &u[1*2+parity],idxh);

                    }//end d0
                }//end d1
            }//end d2
        }//end d3

        rel=0,iml=0;
#pragma acc kernels present(loopplk1) 
#pragma acc loop independent reduction(+:rel) reduction(+:iml)
        for(int i=0;i<vol31h;i++){
            d_complex r = su3_plk1_trace(&loopplk1[0],i);
            r += su3_plk1_trace(&loopplk1[1],i);
            rel += creal(r);
            iml += cimag(r);
        }

    }
    free(loopplk1);

    double trr,tri;
#if NRANKS_D3 > 1
     MPI_Allreduce((void*)&rel,(void*)&trr,
             1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);
     MPI_Allreduce((void*)&iml,(void*)&tri,
             1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);
#else
     trr = rel;
     tri = iml;
#endif 

    return (trr + tri * I)/(vol31h*2*3)/devinfo.nranks;

}

d_complex polyakov_loop2(__restrict const su3_soa * const u)
{

    if(geom_par.nranks[2] !=1 ){
        printf("MPI%02d, POLYAKOV LOOP CALCULATION NOT IMPLEMENTED!\n", devinfo.myrank);
        return 0+0*I;
    }

    su3_plk2 *loopplk2;
    double rel,iml;
    int allocation_check = posix_memalign((void **)&loopplk2, ALIGN, 2*sizeof(su3_plk2));
    ALLOCCHECK(allocation_check,loopplk2);

#pragma acc data  create(loopplk2[0:2])
    {
#pragma acc kernels present(loopplk2) 
#pragma acc loop independent
        for(int i=0;i<vol32h;i++){
            loopplk2[0].r0.c0[i] = 1;loopplk2[0].r0.c1[i] = 0;loopplk2[0].r0.c2[i] = 0;
            loopplk2[0].r1.c0[i] = 0;loopplk2[0].r1.c1[i] = 1;loopplk2[0].r1.c2[i] = 0;

            loopplk2[1].r0.c0[i] = 1;loopplk2[1].r0.c1[i] = 0;loopplk2[1].r0.c2[i] = 0;
            loopplk2[1].r1.c0[i] = 0;loopplk2[1].r1.c1[i] = 1;loopplk2[1].r1.c2[i] = 0;
        }
        int d0, d1, d2, d3,h,idxh,parity;
#pragma acc kernels present(u) present(loopplk2)
#pragma acc loop gang  // NOT INDEPENDENT
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent gang vector
            for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent gang vector
                for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent vector //vector(DIM_BLOCK_X)
                    for(d0=0; d0 < nd0; d0++) {	    
                        int ld0,ld1,ld2,ld3;
                        ld0 = d0 - D0_HALO;
                        ld1 = d1 - D1_HALO;
                        ld2 = d2 - D2_HALO;
                        ld3 = d3 - D3_HALO;
                        int parity2 = (ld0+ld1+ 0+ld3) % 2;
                        int id2h = (ld0+LOC_N0*(ld1+LOC_N1*ld3))/2;

                        parity = (d0+d1+d2+d3) % 2;
                        idxh = snum_acc(d0,d1,d2,d3);  	

                        su3_plk2_su3_soa_multinplace(&loopplk2[parity2],id2h,
                                &u[2*2+parity],idxh);

                    }//end d0
                }//end d1
            }//end d2
        }//end d3

        rel=0,iml=0;
#pragma acc kernels present(loopplk2) 
#pragma acc loop independent reduction(+:rel) reduction(+:iml)
        for(int i=0;i<vol32h;i++){
            d_complex r = su3_plk2_trace(&loopplk2[0],i);
            r += su3_plk2_trace(&loopplk2[1],i);
            rel += creal(r);
            iml += cimag(r);
        }

    }
    free(loopplk2);


    double trr,tri;
#if NRANKS_D3 > 1
     MPI_Allreduce((void*)&rel,(void*)&trr,
             1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);
     MPI_Allreduce((void*)&iml,(void*)&tri,
             1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);
#else
     trr = rel;
     tri = iml;
#endif 

    return (trr + tri * I)/(vol32h*2*3)/devinfo.nranks;

}

d_complex polyakov_loop3(__restrict const su3_soa * const u)
{

    if(geom_par.nranks[3] !=1 ){
        printf("MPI%02d, POLYAKOV LOOP CALCULATION NOT IMPLEMENTED!\n", devinfo.myrank);
        return 0+0*I;
    }

    su3_plk3 *loopplk3;
    double rel,iml;
    int allocation_check = posix_memalign((void **)&loopplk3, ALIGN, 2*sizeof(su3_plk3));
    ALLOCCHECK(allocation_check,loopplk3);

#pragma acc data  create(loopplk3[0:2])
    {
#pragma acc kernels present(loopplk3) 
#pragma acc loop independent
        for(int i=0;i<vol33h;i++){
            loopplk3[0].r0.c0[i] = 1;loopplk3[0].r0.c1[i] = 0;loopplk3[0].r0.c2[i] = 0;
            loopplk3[0].r1.c0[i] = 0;loopplk3[0].r1.c1[i] = 1;loopplk3[0].r1.c2[i] = 0;

            loopplk3[1].r0.c0[i] = 1;loopplk3[1].r0.c1[i] = 0;loopplk3[1].r0.c2[i] = 0;
            loopplk3[1].r1.c0[i] = 0;loopplk3[1].r1.c1[i] = 1;loopplk3[1].r1.c2[i] = 0;
        }
        int d0, d1, d2, d3,h,idxh,parity;
#pragma acc kernels present(u) present(loopplk3)
#pragma acc loop gang // NOT INDEPENDENT
        for(d3=D3_HALO; d3<nd3-D3_HALO; d3++){
#pragma acc loop independent gang vector
            for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent gang vector
                for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent vector //vector(DIM_BLOCK_X)
                    for(d0=0; d0 < nd0; d0++) {	    
                        int ld0,ld1,ld2,ld3;
                        ld0 = d0 - D0_HALO;
                        ld1 = d1 - D1_HALO;
                        ld2 = d2 - D2_HALO;
                        ld3 = d3 - D3_HALO;

                        int parity3 = (ld0+ld1+ld2+0) % 2;
                        int id3h = (ld0+LOC_N0*(ld1+LOC_N1*ld2))/2;

                        parity = (d0+d1+d2+d3) % 2;
                        idxh = snum_acc(d0,d1,d2,d3);  	
 

                        su3_plk3_su3_soa_multinplace(&loopplk3[parity3],id3h,
                                &u[3*2+parity],idxh);

                    }//end d0
                }//end d1
            }//end d2
        }//end d3

        rel=0,iml=0;
#pragma acc kernels present(loopplk3) 
#pragma acc loop independent reduction(+:rel) reduction(+:iml)
        for(int i=0;i<vol33h;i++){
            d_complex r = su3_plk3_trace(&loopplk3[0],i);
            r += su3_plk3_trace(&loopplk3[1],i);
            rel += creal(r);
            iml += cimag(r);
        }

    }
    free(loopplk3);

    double trr,tri;
#if NRANKS_D3 > 1
     MPI_Allreduce((void*)&rel,(void*)&trr,
             1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);
     MPI_Allreduce((void*)&iml,(void*)&tri,
             1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);
#else
     trr = rel;
     tri = iml;
#endif 

    return (trr + tri * I)/(vol33h*2*3)/devinfo.nranks;

}


d_complex (*polyakov_loop[4])(__restrict const su3_soa * const u) =
{ polyakov_loop0, polyakov_loop1, polyakov_loop2, polyakov_loop3};




#endif
