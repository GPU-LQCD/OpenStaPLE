#ifndef BARYON_NUMBER_UTILITIES_C_
#define BARYON_NUMBER_UTILITIES_C_
#include "../OpenAcc/struct_c_def.h"
#include "../OpenAcc/matvecmul.h"
#include "./baryon_number_utilities.h"
#include "../OpenAcc/geometry.h"

#ifdef MULTIDEVICE
#include "../Mpi/communications.h"
#endif

// FIRST DERIVATIVES

void dM_dmu_eo0( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)
{
    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent
    for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent
                for(d0h=0; d0h < nd0h; d0h++) {


                    int d0, d0m, d0p, idxh;
                    vec3 aux;

                    // (d0+d1+d2+d3) even
                    d0 = 2*d0h + ((d1+d2+d3) & 0x1);

                    idxh = snum_acc(d0,d1,d2,d3);

                    // DIRECTION DEPENDENT START
                    d0m = d0 - 1;
                    d0m = d0m + (((d0m >> 31) & 0x1) * nd0);

                    d0p = d0 + 1;
                    d0p *= (((d0p-nd0) >> 31) & 0x1);

                    //direzione tempo positiva,
                    aux = mat_vec_mul_arg( &u[0], idxh, in, snum_acc(d0p,d1,d2,d3),
                            &backfield[0] );

                    // direzione tempo negativa
                    aux = sumResult(aux, conjmat_vec_mul_arg( &u[1],
                                snum_acc(d0m,d1,d2,d3),in,snum_acc(d0m,d1,d2,d3),&backfield[1]));

                    // DIRECTION DEPENDENT END
                    /////////////////////////////////////////////////////////////////////     

                    out->c0[idxh] = (aux.c0)*0.5;
                    out->c1[idxh] = (aux.c1)*0.5;
                    out->c2[idxh] = (aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3

#ifdef MULTIDEVICE
    communicate_fermion_borders(out);
#endif

}
void dM_dmu_eo1( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)
{
    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent
    for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent
                for(d0h=0; d0h < nd0h; d0h++) {


                    int d0, d1m, d1p, idxh;
                    vec3 aux;

                    // (d0+d1+d2+d3) even
                    d0 = 2*d0h + ((d1+d2+d3) & 0x1);

                    idxh = snum_acc(d0,d1,d2,d3);

                    // DIRECTION DEPENDENT START
                    d1m = d1 - 1;
                    d1m = d1m + (((d1m >> 31) & 0x1) * nd1);

                    d1p = d1 + 1;
                    d1p *= (((d1p-nd1) >> 31) & 0x1);

                    //direzione tempo positiva,
                    aux = mat_vec_mul_arg( &u[2], idxh, in, snum_acc(d0,d1p,d2,d3),
                            &backfield[2] );

                    // direzione tempo negativa
                    aux = sumResult(aux, conjmat_vec_mul_arg( &u[3],
                                snum_acc(d0,d1m,d2,d3),in,snum_acc(d0,d1m,d2,d3),&backfield[3]));

                    // DIRECTION DEPENDENT END
                    /////////////////////////////////////////////////////////////////////     

                    out->c0[idxh] = (aux.c0)*0.5;
                    out->c1[idxh] = (aux.c1)*0.5;
                    out->c2[idxh] = (aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3
#ifdef MULTIDEVICE
    communicate_fermion_borders(out);
#endif
}
void dM_dmu_eo2( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)
{
    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent
    for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent
                for(d0h=0; d0h < nd0h; d0h++) {


                    int d0, d2m, d2p, idxh;
                    vec3 aux;

                    // (d0+d1+d2+d3) even
                    d0 = 2*d0h + ((d1+d2+d3) & 0x1);

                    idxh = snum_acc(d0,d1,d2,d3);

                    // DIRECTION DEPENDENT START
                    d2m = d2 - 1;
                    d2m = d2m + (((d2m >> 31) & 0x1) * nd2);

                    d2p = d2 + 1;
                    d2p *= (((d2p-nd2) >> 31) & 0x1);

                    //direzione tempo positiva,
                    aux = mat_vec_mul_arg( &u[4], idxh, in, snum_acc(d0,d1,d2p,d3),
                            &backfield[4] );

                    // direzione tempo negativa
                    aux = sumResult(aux, conjmat_vec_mul_arg( &u[5],
                                snum_acc(d0,d1,d2m,d3),in,snum_acc(d0,d1,d2m,d3),&backfield[5]));

                    // DIRECTION DEPENDENT END
                    /////////////////////////////////////////////////////////////////////     

                    out->c0[idxh] = (aux.c0)*0.5;
                    out->c1[idxh] = (aux.c1)*0.5;
                    out->c2[idxh] = (aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3
#ifdef MULTIDEVICE
    communicate_fermion_borders(out);
#endif
}
void dM_dmu_eo3( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)
{
    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent
    for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent
                for(d0h=0; d0h < nd0h; d0h++) {


                    int d0, d3m, d3p, idxh;
                    vec3 aux;

                    // (d0+d1+d2+d3) even
                    d0 = 2*d0h + ((d1+d2+d3) & 0x1);

                    idxh = snum_acc(d0,d1,d2,d3);

                    // DIRECTION DEPENDENT START
                    d3m = d3 - 1;
                    d3m = d3m + (((d3m >> 31) & 0x1) * nd3);

                    d3p = d3 + 1;
                    d3p *= (((d3p-nd3) >> 31) & 0x1);

                    //direzione tempo positiva,
                    aux = mat_vec_mul_arg( &u[6], idxh, in, snum_acc(d0,d1,d2,d3p),
                            &backfield[6] );

                    // direzione tempo negativa
                    aux = sumResult(aux, conjmat_vec_mul_arg( &u[7],
                                snum_acc(d0,d1,d2,d3m),in,snum_acc(d0,d1,d2,d3m),&backfield[7]));

                    // DIRECTION DEPENDENT END
                    /////////////////////////////////////////////////////////////////////     

                    out->c0[idxh] = (aux.c0)*0.5;
                    out->c1[idxh] = (aux.c1)*0.5;
                    out->c2[idxh] = (aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3
#ifdef MULTIDEVICE
    communicate_fermion_borders(out);
#endif
}

void (*dM_dmu_eo[4])( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)=
{dM_dmu_eo0,dM_dmu_eo1,dM_dmu_eo2,dM_dmu_eo3};

void dM_dmu_oe0( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)
{
    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent
    for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent
                for(d0h=0; d0h < nd0h; d0h++) {

                    int d0, d0m, d0p, idxh;
                    vec3 aux;

                    // (d0+d1+d2+d3) odd
                    d0 = 2*d0h + ((d1+d2+d3+1) & 0x1);

                    idxh = snum_acc(d0,d1,d2,d3);

                    // DIRECTION DEPENDENT START
                    d0m = d0 - 1;
                    d0m = d0m + (((d0m >> 31) & 0x1) * nd0);

                    d0p = d0 + 1;
                    d0p *= (((d0p-nd0) >> 31) & 0x1);

                    //direzione tempo negativa
                    aux = conjmat_vec_mul_arg( &u[0],snum_acc(d0m,d1,d2,d3), in,
                            snum_acc(d0m,d1,d2,d3), &backfield[0]);

                    // direzione tempo positiva
                    aux = sumResult(aux, mat_vec_mul_arg( &u[1],idxh,
                                in,snum_acc(d0p,d1,d2,d3),&backfield[1]));
                    // DIRECTION DEPENDENT END
                    
                    /////////////////////////////////////////////////////////////////////     

                    out->c0[idxh] = (aux.c0)*0.5;
                    out->c1[idxh] = (aux.c1)*0.5;
                    out->c2[idxh] = (aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3
#ifdef MULTIDEVICE
    communicate_fermion_borders(out);
#endif
}
void dM_dmu_oe1( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)
{
    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent
    for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent
                for(d0h=0; d0h < nd0h; d0h++) {

                    int d0, d1m, d1p, idxh;
                    vec3 aux;

                    // (d0+d1+d2+d3) odd
                    d0 = 2*d0h + ((d1+d2+d3+1) & 0x1);

                    idxh = snum_acc(d0,d1,d2,d3);

                    // DIRECTION DEPENDENT START
                    d1m = d1 - 1;
                    d1m = d1m + (((d1m >> 31) & 0x1) * nd1);

                    d1p = d1 + 1;
                    d1p *= (((d1p-nd1) >> 31) & 0x1);

                    //direzione tempo negativa
                    aux = conjmat_vec_mul_arg( &u[2],snum_acc(d0,d1m,d2,d3), in,
                            snum_acc(d0,d1m,d2,d3), &backfield[2]);

                    // direzione tempo positiva
                    aux = sumResult(aux, mat_vec_mul_arg( &u[3],idxh,
                                in,snum_acc(d0,d1p,d2,d3),&backfield[3]));
                    // DIRECTION DEPENDENT END
                    
                    /////////////////////////////////////////////////////////////////////     

                    out->c0[idxh] = (aux.c0)*0.5;
                    out->c1[idxh] = (aux.c1)*0.5;
                    out->c2[idxh] = (aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3
#ifdef MULTIDEVICE
    communicate_fermion_borders(out);
#endif
}
void dM_dmu_oe2( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)
{
    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent
    for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent
                for(d0h=0; d0h < nd0h; d0h++) {

                    int d0, d2m, d2p, idxh;
                    vec3 aux;

                    // (d0+d1+d2+d3) odd
                    d0 = 2*d0h + ((d1+d2+d3+1) & 0x1);

                    idxh = snum_acc(d0,d1,d2,d3);

                    // DIRECTION DEPENDENT START
                    d2m = d2 - 1;
                    d2m = d2m + (((d2m >> 31) & 0x1) * nd2);

                    d2p = d2 + 1;
                    d2p *= (((d2p-nd2) >> 31) & 0x1);

                    //direzione tempo negativa
                    aux = conjmat_vec_mul_arg( &u[4],snum_acc(d0,d1,d2m,d3), in,
                            snum_acc(d0,d1,d2m,d3), &backfield[4]);

                    // direzione tempo positiva
                    aux = sumResult(aux, mat_vec_mul_arg( &u[5],idxh,
                                in,snum_acc(d0,d1,d2p,d3),&backfield[5]));
                    // DIRECTION DEPENDENT END
                    
                    /////////////////////////////////////////////////////////////////////     

                    out->c0[idxh] = (aux.c0)*0.5;
                    out->c1[idxh] = (aux.c1)*0.5;
                    out->c2[idxh] = (aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3
#ifdef MULTIDEVICE
    communicate_fermion_borders(out);
#endif
}
void dM_dmu_oe3( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)
{
    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent
    for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent
                for(d0h=0; d0h < nd0h; d0h++) {

                    int d0, d3m, d3p, idxh;
                    vec3 aux;

                    // (d0+d1+d2+d3) odd
                    d0 = 2*d0h + ((d1+d2+d3+1) & 0x1);

                    idxh = snum_acc(d0,d1,d2,d3);

                    // DIRECTION DEPENDENT START
                    d3m = d3 - 1;
                    d3m = d3m + (((d3m >> 31) & 0x1) * nd3);

                    d3p = d3 + 1;
                    d3p *= (((d3p-nd3) >> 31) & 0x1);

                    //direzione tempo negativa
                    aux = conjmat_vec_mul_arg( &u[6],snum_acc(d0,d1,d2,d3m), in,
                            snum_acc(d0,d1,d2,d3m), &backfield[6]);

                    // direzione tempo positiva
                    aux = sumResult(aux, mat_vec_mul_arg( &u[7],idxh,
                                in,snum_acc(d0,d1,d2,d3p),&backfield[7]));
                    // DIRECTION DEPENDENT END
                    
                    /////////////////////////////////////////////////////////////////////     

                    out->c0[idxh] = (aux.c0)*0.5;
                    out->c1[idxh] = (aux.c1)*0.5;
                    out->c2[idxh] = (aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3
#ifdef MULTIDEVICE
    communicate_fermion_borders(out);
#endif
}

void (*dM_dmu_oe[4]) ( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)=
{dM_dmu_oe0,dM_dmu_oe1,dM_dmu_oe2,dM_dmu_oe3};

// SECOND DERIVATIVES

void d2M_dmu2_eo0( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)
{
    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent
    for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent
                for(d0h=0; d0h < nd0h; d0h++) {


                    int d0, d0m, d0p, idxh;
                    vec3 aux;

                    // (d0+d1+d2+d3) even
                    d0 = 2*d0h + ((d1+d2+d3) & 0x1);

                    idxh = snum_acc(d0,d1,d2,d3);

                    // DIRECTION DEPENDENT START
                    d0m = d0 - 1;
                    d0m = d0m + (((d0m >> 31) & 0x1) * nd0);

                    d0p = d0 + 1;
                    d0p *= (((d0p-nd0) >> 31) & 0x1);

                    //direzione tempo positiva,
                    aux = mat_vec_mul_arg( &u[0], idxh, in, snum_acc(d0p,d1,d2,d3),
                            &backfield[0] );

                    // direzione tempo negativa
                    aux = subResult(aux, conjmat_vec_mul_arg( &u[1],
                                snum_acc(d0m,d1,d2,d3),in,snum_acc(d0m,d1,d2,d3),&backfield[1]));

                    // DIRECTION DEPENDENT END
                    /////////////////////////////////////////////////////////////////////     

                    out->c0[idxh] = (aux.c0)*0.5;
                    out->c1[idxh] = (aux.c1)*0.5;
                    out->c2[idxh] = (aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3
#ifdef MULTIDEVICE
    communicate_fermion_borders(out);
#endif
}
void d2M_dmu2_eo1( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)
{
    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent
    for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent
                for(d0h=0; d0h < nd0h; d0h++) {


                    int d0, d1m, d1p, idxh;
                    vec3 aux;

                    // (d0+d1+d2+d3) even
                    d0 = 2*d0h + ((d1+d2+d3) & 0x1);

                    idxh = snum_acc(d0,d1,d2,d3);

                    // DIRECTION DEPENDENT START
                    d1m = d1 - 1;
                    d1m = d1m + (((d1m >> 31) & 0x1) * nd1);

                    d1p = d1 + 1;
                    d1p *= (((d1p-nd1) >> 31) & 0x1);

                    //direzione tempo positiva,
                    aux = mat_vec_mul_arg( &u[2], idxh, in, snum_acc(d0,d1p,d2,d3),
                            &backfield[2] );

                    // direzione tempo negativa
                    aux = subResult(aux, conjmat_vec_mul_arg( &u[3],
                                snum_acc(d0,d1m,d2,d3),in,snum_acc(d0,d1m,d2,d3),&backfield[3]));

                    // DIRECTION DEPENDENT END
                    /////////////////////////////////////////////////////////////////////     

                    out->c0[idxh] = (aux.c0)*0.5;
                    out->c1[idxh] = (aux.c1)*0.5;
                    out->c2[idxh] = (aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3
#ifdef MULTIDEVICE
    communicate_fermion_borders(out);
#endif
}
void d2M_dmu2_eo2( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)
{
    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent
    for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent
                for(d0h=0; d0h < nd0h; d0h++) {


                    int d0, d2m, d2p, idxh;
                    vec3 aux;

                    // (d0+d1+d2+d3) even
                    d0 = 2*d0h + ((d1+d2+d3) & 0x1);

                    idxh = snum_acc(d0,d1,d2,d3);

                    // DIRECTION DEPENDENT START
                    d2m = d2 - 1;
                    d2m = d2m + (((d2m >> 31) & 0x1) * nd2);

                    d2p = d2 + 1;
                    d2p *= (((d2p-nd2) >> 31) & 0x1);

                    //direzione tempo positiva,
                    aux = mat_vec_mul_arg( &u[4], idxh, in, snum_acc(d0,d1,d2p,d3),
                            &backfield[4] );

                    // direzione tempo negativa
                    aux = subResult(aux, conjmat_vec_mul_arg( &u[5],
                                snum_acc(d0,d1,d2m,d3),in,snum_acc(d0,d1,d2m,d3),&backfield[5]));

                    // DIRECTION DEPENDENT END
                    /////////////////////////////////////////////////////////////////////     

                    out->c0[idxh] = (aux.c0)*0.5;
                    out->c1[idxh] = (aux.c1)*0.5;
                    out->c2[idxh] = (aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3
#ifdef MULTIDEVICE
    communicate_fermion_borders(out);
#endif
}
void d2M_dmu2_eo3( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)
{
    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent
    for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent
                for(d0h=0; d0h < nd0h; d0h++) {


                    int d0, d3m, d3p, idxh;
                    vec3 aux;

                    // (d0+d1+d2+d3) even
                    d0 = 2*d0h + ((d1+d2+d3) & 0x1);

                    idxh = snum_acc(d0,d1,d2,d3);

                    // DIRECTION DEPENDENT START
                    d3m = d3 - 1;
                    d3m = d3m + (((d3m >> 31) & 0x1) * nd3);

                    d3p = d3 + 1;
                    d3p *= (((d3p-nd3) >> 31) & 0x1);

                    //direzione tempo positiva,
                    aux = mat_vec_mul_arg( &u[6], idxh, in, snum_acc(d0,d1,d2,d3p),
                            &backfield[6] );

                    // direzione tempo negativa
                    aux = subResult(aux, conjmat_vec_mul_arg( &u[7],
                                snum_acc(d0,d1,d2,d3m),in,snum_acc(d0,d1,d2,d3m),&backfield[7]));

                    // DIRECTION DEPENDENT END
                    /////////////////////////////////////////////////////////////////////     

                    out->c0[idxh] = (aux.c0)*0.5;
                    out->c1[idxh] = (aux.c1)*0.5;
                    out->c2[idxh] = (aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3
#ifdef MULTIDEVICE
    communicate_fermion_borders(out);
#endif
}

void (*d2M_dmu2_eo[4])( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)=
{d2M_dmu2_eo0,d2M_dmu2_eo1,d2M_dmu2_eo2,d2M_dmu2_eo3};

void d2M_dmu2_oe0( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)
{
    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent
    for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent
                for(d0h=0; d0h < nd0h; d0h++) {

                    int d0, d0m, d0p, idxh;
                    vec3 aux;

                    // (d0+d1+d2+d3) odd
                    d0 = 2*d0h + ((d1+d2+d3+1) & 0x1);

                    idxh = snum_acc(d0,d1,d2,d3);

                    // DIRECTION DEPENDENT START
                    d0m = d0 - 1;
                    d0m = d0m + (((d0m >> 31) & 0x1) * nd0);

                    d0p = d0 + 1;
                    d0p *= (((d0p-nd0) >> 31) & 0x1);

                    //direzione tempo negativa
                    aux = conjmat_vec_mul_arg( &u[0],snum_acc(d0m,d1,d2,d3), in,
                            snum_acc(d0m,d1,d2,d3), &backfield[0]);

                    // direzione tempo positiva
                    aux = subResult(aux, mat_vec_mul_arg( &u[1],idxh,
                                in,snum_acc(d0p,d1,d2,d3),&backfield[1]));
                    // DIRECTION DEPENDENT END
                    
                    /////////////////////////////////////////////////////////////////////     

                    out->c0[idxh] =-(aux.c0)*0.5;
                    out->c1[idxh] =-(aux.c1)*0.5;
                    out->c2[idxh] =-(aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3
#ifdef MULTIDEVICE
    communicate_fermion_borders(out);
#endif
}
void d2M_dmu2_oe1( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)
{
    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent
    for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent
                for(d0h=0; d0h < nd0h; d0h++) {

                    int d0, d1m, d1p, idxh;
                    vec3 aux;

                    // (d0+d1+d2+d3) odd
                    d0 = 2*d0h + ((d1+d2+d3+1) & 0x1);

                    idxh = snum_acc(d0,d1,d2,d3);

                    // DIRECTION DEPENDENT START
                    d1m = d1 - 1;
                    d1m = d1m + (((d1m >> 31) & 0x1) * nd1);

                    d1p = d1 + 1;
                    d1p *= (((d1p-nd1) >> 31) & 0x1);

                    //direzione tempo negativa
                    aux = conjmat_vec_mul_arg( &u[2],snum_acc(d0,d1m,d2,d3), in,
                            snum_acc(d0,d1m,d2,d3), &backfield[2]);

                    // direzione tempo positiva
                    aux = subResult(aux, mat_vec_mul_arg( &u[3],idxh,
                                in,snum_acc(d0,d1p,d2,d3),&backfield[3]));
                    // DIRECTION DEPENDENT END
                    
                    /////////////////////////////////////////////////////////////////////     

                    out->c0[idxh] =-(aux.c0)*0.5;
                    out->c1[idxh] =-(aux.c1)*0.5;
                    out->c2[idxh] =-(aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3
#ifdef MULTIDEVICE
    communicate_fermion_borders(out);
#endif
}
void d2M_dmu2_oe2( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)
{
    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent
    for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent
                for(d0h=0; d0h < nd0h; d0h++) {

                    int d0, d2m, d2p, idxh;
                    vec3 aux;

                    // (d0+d1+d2+d3) odd
                    d0 = 2*d0h + ((d1+d2+d3+1) & 0x1);

                    idxh = snum_acc(d0,d1,d2,d3);

                    // DIRECTION DEPENDENT START
                    d2m = d2 - 1;
                    d2m = d2m + (((d2m >> 31) & 0x1) * nd2);

                    d2p = d2 + 1;
                    d2p *= (((d2p-nd2) >> 31) & 0x1);

                    //direzione tempo negativa
                    aux = conjmat_vec_mul_arg( &u[4],snum_acc(d0,d1,d2m,d3), in,
                            snum_acc(d0,d1,d2m,d3), &backfield[4]);

                    // direzione tempo positiva
                    aux = subResult(aux, mat_vec_mul_arg( &u[5],idxh,
                                in,snum_acc(d0,d1,d2p,d3),&backfield[5]));
                    // DIRECTION DEPENDENT END
                    
                    /////////////////////////////////////////////////////////////////////     

                    out->c0[idxh] =-(aux.c0)*0.5;
                    out->c1[idxh] =-(aux.c1)*0.5;
                    out->c2[idxh] =-(aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3
#ifdef MULTIDEVICE
    communicate_fermion_borders(out);
#endif
}
void d2M_dmu2_oe3( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)
{
    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent
    for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent
                for(d0h=0; d0h < nd0h; d0h++) {

                    int d0, d3m, d3p, idxh;
                    vec3 aux;

                    // (d0+d1+d2+d3) odd
                    d0 = 2*d0h + ((d1+d2+d3+1) & 0x1);

                    idxh = snum_acc(d0,d1,d2,d3);

                    // DIRECTION DEPENDENT START
                    d3m = d3 - 1;
                    d3m = d3m + (((d3m >> 31) & 0x1) * nd3);

                    d3p = d3 + 1;
                    d3p *= (((d3p-nd3) >> 31) & 0x1);

                    //direzione tempo negativa
                    aux = conjmat_vec_mul_arg( &u[6],snum_acc(d0,d1,d2,d3m), in,
                            snum_acc(d0,d1,d2,d3m), &backfield[6]);

                    // direzione tempo positiva
                    aux = subResult(aux, mat_vec_mul_arg( &u[7],idxh,
                                in,snum_acc(d0,d1,d2,d3p),&backfield[7]));
                    // DIRECTION DEPENDENT END
                    
                    /////////////////////////////////////////////////////////////////////     

                    out->c0[idxh] =-(aux.c0)*0.5;
                    out->c1[idxh] =-(aux.c1)*0.5;
                    out->c2[idxh] =-(aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3
#ifdef MULTIDEVICE
    communicate_fermion_borders(out);
#endif
}

void (*d2M_dmu2_oe[4]) ( __restrict const su3_soa * u,
        vec3_soa * out,
        __restrict const vec3_soa * in,
        __restrict const double_soa * backfield)=
{d2M_dmu2_oe0,d2M_dmu2_oe1,d2M_dmu2_oe2,d2M_dmu2_oe3};



#endif
