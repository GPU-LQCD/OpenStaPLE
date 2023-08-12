#ifndef FERMION_FORCE_UTILITIES_H
#define FERMION_FORCE_UTILITIES_H

#include "../Include/common_defines.h"
#include "./struct_c_def.h"
#include "../Include/fermion_parameters.h"

// if using GCC, there are some problems with __restrict.

#ifdef __GNUC__
#define __restrict
#endif


#pragma acc routine seq
static inline void vec1_directprod_conj_vec2_into_mat1( 
        __restrict su3_soa * const aux_u,   int idxh,
        __restrict const vec3_soa  * const fer_l, int idl,
        __restrict const vec3_soa  * const fer_r, int idr,
        double factor)
{
    // NOTE: maybe it is possible to avoid some computation using just real and imaginary parts,
    // and writing products explicitly
    // especially for diagonal entries
    d_complex r0 = factor * conj(fer_r->c0[idr]);
    d_complex r1 = factor * conj(fer_r->c1[idr]);
    d_complex r2 = factor * conj(fer_r->c2[idr]);
    d_complex l0 = fer_l->c0[idl];
    d_complex l1 = fer_l->c1[idl];
    d_complex l2 = fer_l->c2[idl];

    aux_u->r0.c0[idxh] += l0*r0;
    aux_u->r0.c1[idxh] += l0*r1;
    aux_u->r0.c2[idxh] += l0*r2;
    aux_u->r1.c0[idxh] += l1*r0;
    aux_u->r1.c1[idxh] += l1*r1;
    aux_u->r1.c2[idxh] += l1*r2;
    aux_u->r2.c0[idxh] += l2*r0;
    aux_u->r2.c1[idxh] += l2*r1;
    aux_u->r2.c2[idxh] += l2*r2;

}


#pragma acc routine seq
static inline  void mat1_times_auxmat_into_tamat( 
        __restrict const su3_soa * const mat1,  const  int idx,
        __restrict const su3_soa * const auxmat, const  int idx_aux,
        __restrict tamat_soa * const ipdot, const  int idipdot
        ,d_complex phase        )
{
    d_complex mat1_00 = mat1->r0.c0[idx];
    d_complex mat1_01 = mat1->r0.c1[idx];
    d_complex mat1_02 = mat1->r0.c2[idx];
    d_complex mat1_10 = mat1->r1.c0[idx];
    d_complex mat1_11 = mat1->r1.c1[idx];
    d_complex mat1_12 = mat1->r1.c2[idx];
    // compute 3rd matrix row from the first two
    d_complex mat1_20 = conj( ( mat1_01 * mat1_12 ) - ( mat1_02 * mat1_11) ) ;
    d_complex mat1_21 = conj( ( mat1_02 * mat1_10 ) - ( mat1_00 * mat1_12) ) ;
    d_complex mat1_22 = conj( ( mat1_00 * mat1_11 ) - ( mat1_01 * mat1_10) ) ;

    mat1_00 *= phase;
    mat1_01 *= phase;
    mat1_02 *= phase;
    mat1_10 *= phase;
    mat1_11 *= phase;
    mat1_12 *= phase;
    mat1_20 *= phase;
    mat1_21 *= phase;
    mat1_22 *= phase;

    d_complex auxmat_00 = auxmat->r0.c0[idx_aux];
    d_complex auxmat_01 = auxmat->r0.c1[idx_aux];
    d_complex auxmat_02 = auxmat->r0.c2[idx_aux];
    d_complex auxmat_10 = auxmat->r1.c0[idx_aux];
    d_complex auxmat_11 = auxmat->r1.c1[idx_aux];
    d_complex auxmat_12 = auxmat->r1.c2[idx_aux];
    d_complex auxmat_20 = auxmat->r2.c0[idx_aux];
    d_complex auxmat_21 = auxmat->r2.c1[idx_aux];
    d_complex auxmat_22 = auxmat->r2.c2[idx_aux];

    // product;
    d_complex a00 = mat1_00 * auxmat_00 + mat1_01 * auxmat_10 + mat1_02 * auxmat_20;
    d_complex a01 = mat1_00 * auxmat_01 + mat1_01 * auxmat_11 + mat1_02 * auxmat_21;
    d_complex a02 = mat1_00 * auxmat_02 + mat1_01 * auxmat_12 + mat1_02 * auxmat_22;

    mat1_00 = mat1_10 * auxmat_00 + mat1_11 * auxmat_10 + mat1_12 * auxmat_20;
    mat1_01 = mat1_10 * auxmat_01 + mat1_11 * auxmat_11 + mat1_12 * auxmat_21;
    mat1_02 = mat1_10 * auxmat_02 + mat1_11 * auxmat_12 + mat1_12 * auxmat_22;

    mat1_10 = mat1_20 * auxmat_00 + mat1_21 * auxmat_10 + mat1_22 * auxmat_20;
    mat1_11 = mat1_20 * auxmat_01 + mat1_21 * auxmat_11 + mat1_22 * auxmat_21; 
    mat1_12 = mat1_20 * auxmat_02 + mat1_21 * auxmat_12 + mat1_22 * auxmat_22;

    ipdot->c01[idipdot]  -= 0.5*((a01) - conj(mat1_00));
    ipdot->c02[idipdot]  -= 0.5*((a02) - conj(mat1_10));
    ipdot->c12[idipdot]  -= 0.5*((mat1_02) - conj(mat1_11));
    ipdot->ic00[idipdot] -= cimag(a00)-ONE_BY_THREE*(cimag(a00)+cimag(mat1_01)+cimag(mat1_12));
    ipdot->ic11[idipdot] -= cimag(mat1_01)-ONE_BY_THREE*(cimag(a00)+cimag(mat1_01)+cimag(mat1_12));

}


#pragma acc routine seq
static inline  void mat1_times_auxmat_into_tamat_nophase(  
        __restrict const su3_soa * const mat1, const  int idx,
        __restrict const su3_soa * const auxmat, const  int idx_aux,
        __restrict tamat_soa * const ipdot, const  int idipdot)
{
    d_complex mat1_00 = mat1->r0.c0[idx];
    d_complex mat1_01 = mat1->r0.c1[idx];
    d_complex mat1_02 = mat1->r0.c2[idx];
    d_complex mat1_10 = mat1->r1.c0[idx];
    d_complex mat1_11 = mat1->r1.c1[idx];
    d_complex mat1_12 = mat1->r1.c2[idx];
    // compute 3rd matrix row from the first two
    d_complex mat1_20 = conj( ( mat1_01 * mat1_12 ) - ( mat1_02 * mat1_11) ) ;
    d_complex mat1_21 = conj( ( mat1_02 * mat1_10 ) - ( mat1_00 * mat1_12) ) ;
    d_complex mat1_22 = conj( ( mat1_00 * mat1_11 ) - ( mat1_01 * mat1_10) ) ;
    d_complex auxmat_00 = auxmat->r0.c0[idx_aux];
    d_complex auxmat_01 = auxmat->r0.c1[idx_aux];
    d_complex auxmat_02 = auxmat->r0.c2[idx_aux];
    d_complex auxmat_10 = auxmat->r1.c0[idx_aux];
    d_complex auxmat_11 = auxmat->r1.c1[idx_aux];
    d_complex auxmat_12 = auxmat->r1.c2[idx_aux];
    d_complex auxmat_20 = auxmat->r2.c0[idx_aux];
    d_complex auxmat_21 = auxmat->r2.c1[idx_aux];
    d_complex auxmat_22 = auxmat->r2.c2[idx_aux];

    // product;
    d_complex a00 = mat1_00 * auxmat_00 + mat1_01 * auxmat_10 + mat1_02 * auxmat_20;
    d_complex a01 = mat1_00 * auxmat_01 + mat1_01 * auxmat_11 + mat1_02 * auxmat_21;
    d_complex a02 = mat1_00 * auxmat_02 + mat1_01 * auxmat_12 + mat1_02 * auxmat_22;

    mat1_00 = mat1_10 * auxmat_00 + mat1_11 * auxmat_10 + mat1_12 * auxmat_20;
    mat1_01 = mat1_10 * auxmat_01 + mat1_11 * auxmat_11 + mat1_12 * auxmat_21;
    mat1_02 = mat1_10 * auxmat_02 + mat1_11 * auxmat_12 + mat1_12 * auxmat_22;

    mat1_10 = mat1_20 * auxmat_00 + mat1_21 * auxmat_10 + mat1_22 * auxmat_20;
    mat1_11 = mat1_20 * auxmat_01 + mat1_21 * auxmat_11 + mat1_22 * auxmat_21; 
    mat1_12 = mat1_20 * auxmat_02 + mat1_21 * auxmat_12 + mat1_22 * auxmat_22;

    ipdot->c01[idipdot]  -= 0.5*((a01) - conj(mat1_00));
    ipdot->c02[idipdot]  -= 0.5*((a02) - conj(mat1_10));
    ipdot->c12[idipdot]  -= 0.5*((mat1_02) - conj(mat1_11));
    ipdot->ic00[idipdot] -= cimag(a00)-ONE_BY_THREE*(cimag(a00)+cimag(mat1_01)+cimag(mat1_12));
    ipdot->ic11[idipdot] -= cimag(mat1_01)-ONE_BY_THREE*(cimag(a00)+cimag(mat1_01)+cimag(mat1_12));

}


#pragma acc routine seq
static inline  void phase_times_auxmat_into_auxmat(
        __restrict const su3_soa * const auxmat,
        __restrict su3_soa * const pseudo_ipdot,
        const  int idx, d_complex phase )
{

    pseudo_ipdot->r0.c0[idx] += auxmat->r0.c0[idx]*phase ; 
    pseudo_ipdot->r0.c1[idx] += auxmat->r0.c1[idx]*phase ; 
    pseudo_ipdot->r0.c2[idx] += auxmat->r0.c2[idx]*phase ; 
    pseudo_ipdot->r1.c0[idx] += auxmat->r1.c0[idx]*phase ; 
    pseudo_ipdot->r1.c1[idx] += auxmat->r1.c1[idx]*phase ; 
    pseudo_ipdot->r1.c2[idx] += auxmat->r1.c2[idx]*phase ; 
    pseudo_ipdot->r2.c0[idx] += auxmat->r2.c0[idx]*phase ; 
    pseudo_ipdot->r2.c1[idx] += auxmat->r2.c1[idx]*phase ; 
    pseudo_ipdot->r2.c2[idx] += auxmat->r2.c2[idx]*phase ; 


}

#pragma acc routine seq
static inline void accumulate_auxmat1_into_auxmat2(
        __restrict const su3_soa * const auxmat1,
        __restrict su3_soa * const auxmat2,
        const  int idx )
{

    auxmat2->r0.c0[idx] += auxmat1->r0.c0[idx];
    auxmat2->r0.c1[idx] += auxmat1->r0.c1[idx];
    auxmat2->r0.c2[idx] += auxmat1->r0.c2[idx];
    auxmat2->r1.c0[idx] += auxmat1->r1.c0[idx];
    auxmat2->r1.c1[idx] += auxmat1->r1.c1[idx];
    auxmat2->r1.c2[idx] += auxmat1->r1.c2[idx];
    auxmat2->r2.c0[idx] += auxmat1->r2.c0[idx];
    auxmat2->r2.c1[idx] += auxmat1->r2.c1[idx];
    auxmat2->r2.c2[idx] += auxmat1->r2.c2[idx];

}


#pragma acc routine seq
static inline void assign_zero_to_tamat_soa_component(
        __restrict tamat_soa * const matrix_comp, int idx)
{
    matrix_comp->c01[idx]=0.0+I*0.0;
    matrix_comp->c02[idx]=0.0+I*0.0;
    matrix_comp->c12[idx]=0.0+I*0.0;
    matrix_comp->ic00[idx]=0.0;
    matrix_comp->ic11[idx]=0.0;
}

void set_tamat_soa_to_zero( __restrict tamat_soa * const matrix);



void direct_product_of_fermions_into_auxmat(
        __restrict const vec3_soa  * const loc_s,
        __restrict const vec3_soa  * const loc_h,
        __restrict su3_soa * const aux_u,
        const RationalApprox * const approx,
				int nnp_openacc[sizeh][4][2],
        int iter);

void multiply_conf_times_force_and_take_ta_nophase(
        __restrict const su3_soa * const u,
        __restrict const su3_soa * const auxmat,
        __restrict tamat_soa * const ipdot);


void multiply_backfield_times_force(
        __restrict ferm_param * const tpars, 
        __restrict const su3_soa * const auxmat,
        __restrict su3_soa * const pseudo_ipdot);

void accumulate_gl3soa_into_gl3soa(
        __restrict const su3_soa * const auxmat,
        __restrict su3_soa * const pseudo_ipdot);


void ker_openacc_compute_fermion_force( 
        __restrict const su3_soa * const u,
        __restrict su3_soa * const aux_u,
        __restrict const vec3_soa * const in_shiftmulti,
        __restrict vec3_soa  * const loc_s,
        __restrict vec3_soa  * const loc_h,
				int nnp_openacc[sizeh][4][2],
        __restrict ferm_param  *  tpars
        );

#endif
