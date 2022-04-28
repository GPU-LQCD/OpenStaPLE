#ifndef SINGLE_TYPES_H
#define SINGLE_TYPES_H

#include <complex.h>
#include "./struct_c_def.h"
#include <stdio.h>
#include <math.h>

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif


/****************************************************
 * (Perhaps) Check if adding                        *
 *                                                  *
 * #pragma acc routine seq                          *
 *                                                  *
 *before each function can improve the performances.*
 ****************************************************/

//Type definitions
typedef struct single_su3_t {
  d_complex comp[3][3];
} single_su3;
typedef struct single_tamat_t {
  d_complex c01; // comp_01
  d_complex c02; // comp_02
  d_complex c12; // comp_12
  double ic00;   // Im(comp_00)
  double ic11;   // Im(comp_11)
} single_tamat;
typedef struct single_thmat_t {
  d_complex c01; // comp_01
  d_complex c02; // comp_02
  d_complex c12; // comp_12
  double rc00;   // Re(comp_00)
  double rc11;   // Re(comp_11)
} single_thmat;


// 3rd row reconstruction on site
#pragma acc routine seq
static inline void rebuild3row(single_su3 *AUX){

  AUX->comp[2][0] = conj(AUX->comp[0][1] * AUX->comp[1][2] - AUX->comp[0][2] * AUX->comp[1][1]);
  AUX->comp[2][1] = conj(AUX->comp[0][2] * AUX->comp[1][0] - AUX->comp[0][0] * AUX->comp[1][2]);
  AUX->comp[2][2] = conj(AUX->comp[0][0] * AUX->comp[1][1] - AUX->comp[0][1] * AUX->comp[1][0]);

}
//Conversion functions
#pragma acc routine seq
static inline void thmat_to_su3(single_su3 * out, single_thmat * Q){

    out->comp[0][0] =      Q->rc00 ; 
    out->comp[0][1] =      Q->c01 ;
    out->comp[0][2] =      Q->c02 ;
    out->comp[1][0] = conj(Q->c01) ;
    out->comp[1][1] =      Q->rc11 ; 
    out->comp[1][2] =      Q->c12 ;
    out->comp[2][0] = conj(Q->c02);
    out->comp[2][1] = conj(Q->c12);
    out->comp[2][2] =    - Q->rc00 - Q->rc11 ; 


}
#pragma acc routine seq
static inline void tamat_to_su3(single_su3 * out, single_tamat * QA){

    out->comp[0][0] =     I* QA->ic00 ; 
    out->comp[0][1] =      QA->c01 ;
    out->comp[0][2] =      QA->c02 ;
    out->comp[1][0] = -conj(QA->c01) ;
    out->comp[1][1] =     I* QA->ic11 ; 
    out->comp[1][2] =      QA->c12 ;
    out->comp[2][0] = -conj(QA->c02);
    out->comp[2][1] = -conj(QA->c12);
    out->comp[2][2] =    (- QA->ic00 - QA->ic11)*I ; 


}
#pragma acc routine seq
static inline void i_times_tamat_to_su3(single_su3 * out, single_tamat * QA){
  //I*tamat is a thmat
    out->comp[0][0] =           - QA->ic00 ; 
    out->comp[0][1] =   (1.0*I) * QA->c01  ;
    out->comp[0][2] =   (1.0*I) * QA->c02  ;

    out->comp[1][0] =   (-1.0*I) * conj(QA->c01) ;
    out->comp[1][1] =           - QA->ic11 ; 
    out->comp[1][2] =   (1.0*I) * QA->c12  ;

    out->comp[2][0] =   (-1.0*I) * conj(QA->c02) ;
    out->comp[2][1] =   (-1.0*I) * conj(QA->c12) ;
    out->comp[2][2] =   QA->ic00 + QA->ic11; 
}
#pragma acc routine seq
static inline void gl3_to_tamat(single_su3 * in , single_tamat * out){

     double im_trace_third = (cimag(in->comp[0][0]) +  cimag(in->comp[1][1]) + cimag(in->comp[2][2]))/3;

     out->c01 = 0.5*(in->comp[0][1] - conj(in->comp[1][0]));
     out->c02 = 0.5*(in->comp[0][2] - conj(in->comp[2][0]));
     out->c12 = 0.5*(in->comp[1][2] - conj(in->comp[2][1]));

     out->ic00 = cimag(in->comp[0][0])-im_trace_third;
     out->ic11 = cimag(in->comp[1][1])-im_trace_third;

}
#pragma acc routine seq
static inline void gl3_to_thmat(single_su3 * in , single_thmat * out){

     double re_trace_third = (creal(in->comp[0][0]) +  creal(in->comp[1][1]) + creal(in->comp[2][2]))/3;

     out->c01 = 0.5*(in->comp[0][1]+ conj(in->comp[1][0]));
     out->c02 = 0.5*(in->comp[0][2]+ conj(in->comp[2][0]));
     out->c12 = 0.5*(in->comp[1][2]+ conj(in->comp[2][1]));

     out->rc00 = creal(in->comp[0][0])-re_trace_third;
     out->rc11 = creal(in->comp[1][1])-re_trace_third;

}

//extraction from soas
#pragma acc routine seq
static inline void single_su3_from_su3_soa( __restrict const su3_soa * const mat, const int idx_mat, single_su3 * Omat){
  Omat->comp[0][0] = mat->r0.c0[idx_mat];
  Omat->comp[0][1] = mat->r0.c1[idx_mat];
  Omat->comp[0][2] = mat->r0.c2[idx_mat];
  Omat->comp[1][0] = mat->r1.c0[idx_mat];
  Omat->comp[1][1] = mat->r1.c1[idx_mat];
  Omat->comp[1][2] = mat->r1.c2[idx_mat];

}
#pragma acc routine seq
static inline void single_su3_from_global_su3_soa( __restrict const global_su3_soa * const mat, const int idx_mat, single_su3 * Omat){
  Omat->comp[0][0] = mat->r0.c0[idx_mat];
  Omat->comp[0][1] = mat->r0.c1[idx_mat];
  Omat->comp[0][2] = mat->r0.c2[idx_mat];
  Omat->comp[1][0] = mat->r1.c0[idx_mat];
  Omat->comp[1][1] = mat->r1.c1[idx_mat];
  Omat->comp[1][2] = mat->r1.c2[idx_mat];

}
#pragma acc routine seq
static inline void single_gl3_from_su3_soa( __restrict const su3_soa * const mat, const int idx_mat, single_su3 * Omat){
    single_su3_from_su3_soa(mat,idx_mat,Omat);
    Omat->comp[2][0] = mat->r2.c0[idx_mat];
    Omat->comp[2][1] = mat->r2.c1[idx_mat];
    Omat->comp[2][2] = mat->r2.c2[idx_mat];

}
#pragma acc routine seq
static inline void single_gl3_from_global_su3_soa( __restrict const global_su3_soa * const mat, const int idx_mat, single_su3 * Omat){
    single_su3_from_global_su3_soa(mat,idx_mat,Omat);
    Omat->comp[2][0] = mat->r2.c0[idx_mat];
    Omat->comp[2][1] = mat->r2.c1[idx_mat];
    Omat->comp[2][2] = mat->r2.c2[idx_mat];

}
#pragma acc routine seq
static inline void single_tamat_from_tamat_soa(__restrict const tamat_soa * in, int idx, single_tamat * out){

    out->ic00 = in->ic00[idx];
    out->ic11 = in->ic11[idx];
    out->c01 = in->c01[idx];
    out->c02 = in->c02[idx];
    out->c12 = in->c12[idx];

}
#pragma acc routine seq
static inline void single_thmat_from_thmat_soa(__restrict const thmat_soa * in, int idx, single_thmat * out){

    out->rc00 = in->rc00[idx];
    out->rc11 = in->rc11[idx];
    out->c01 = in->c01[idx];
    out->c02 = in->c02[idx];
    out->c12 = in->c12[idx];

}



//insertion into soas
#pragma acc routine seq
static inline void single_su3_into_su3_soa( __restrict su3_soa * const mat, const int idx_mat,  single_su3 * Imat){
  mat->r0.c0[idx_mat] = Imat->comp[0][0];
  mat->r0.c1[idx_mat] = Imat->comp[0][1];
  mat->r0.c2[idx_mat] = Imat->comp[0][2];
  mat->r1.c0[idx_mat] = Imat->comp[1][0];
  mat->r1.c1[idx_mat] = Imat->comp[1][1];
  mat->r1.c2[idx_mat] = Imat->comp[1][2];
}
#pragma acc routine seq
static inline void single_su3_into_global_su3_soa( __restrict global_su3_soa * const mat, const int idx_mat,  single_su3 * Imat){
  mat->r0.c0[idx_mat] = Imat->comp[0][0];
  mat->r0.c1[idx_mat] = Imat->comp[0][1];
  mat->r0.c2[idx_mat] = Imat->comp[0][2];
  mat->r1.c0[idx_mat] = Imat->comp[1][0];
  mat->r1.c1[idx_mat] = Imat->comp[1][1];
  mat->r1.c2[idx_mat] = Imat->comp[1][2];
}
#pragma acc routine seq
static inline void single_su3_addinto_su3_soa( __restrict su3_soa * const mat, const int idx_mat,  single_su3 * Imat){
  mat->r0.c0[idx_mat] += Imat->comp[0][0];
  mat->r0.c1[idx_mat] += Imat->comp[0][1];
  mat->r0.c2[idx_mat] += Imat->comp[0][2];
  mat->r1.c0[idx_mat] += Imat->comp[1][0];
  mat->r1.c1[idx_mat] += Imat->comp[1][1];
  mat->r1.c2[idx_mat] += Imat->comp[1][2];
}
#pragma acc routine seq
static inline void single_su3_subtinto_su3_soa( __restrict su3_soa * const mat, const int idx_mat,  single_su3 * Imat){
  mat->r0.c0[idx_mat] -= Imat->comp[0][0];
  mat->r0.c1[idx_mat] -= Imat->comp[0][1];
  mat->r0.c2[idx_mat] -= Imat->comp[0][2];
  mat->r1.c0[idx_mat] -= Imat->comp[1][0];
  mat->r1.c1[idx_mat] -= Imat->comp[1][1];
  mat->r1.c2[idx_mat] -= Imat->comp[1][2];
}
#pragma acc routine seq
static inline void single_gl3_into_su3_soa( __restrict su3_soa * const mat, const int idx_mat, single_su3 * Imat){

 single_su3_into_su3_soa(mat,idx_mat,Imat);
  mat->r2.c0[idx_mat] = Imat->comp[2][0];
  mat->r2.c1[idx_mat] = Imat->comp[2][1];
  mat->r2.c2[idx_mat] = Imat->comp[2][2];
}
#pragma acc routine seq
static inline void single_gl3_into_global_su3_soa( __restrict global_su3_soa * const mat, const int idx_mat, single_su3 * Imat){

 single_su3_into_global_su3_soa(mat,idx_mat,Imat);
  mat->r2.c0[idx_mat] = Imat->comp[2][0];
  mat->r2.c1[idx_mat] = Imat->comp[2][1];
  mat->r2.c2[idx_mat] = Imat->comp[2][2];
}
#pragma acc routine seq
static inline void single_gl3_addinto_su3_soa( __restrict su3_soa * const mat, const int idx_mat, single_su3 * Imat){

  single_su3_addinto_su3_soa(mat,idx_mat,Imat);
  mat->r2.c0[idx_mat] += Imat->comp[2][0];
  mat->r2.c1[idx_mat] += Imat->comp[2][1];
  mat->r2.c2[idx_mat] += Imat->comp[2][2];
}
#pragma acc routine seq
static inline void single_gl3_subtinto_su3_soa( __restrict su3_soa * const mat, const int idx_mat, single_su3 * Imat){

  single_su3_subtinto_su3_soa(mat,idx_mat,Imat);
  mat->r2.c0[idx_mat] -= Imat->comp[2][0];
  mat->r2.c1[idx_mat] -= Imat->comp[2][1];
  mat->r2.c2[idx_mat] -= Imat->comp[2][2];
}
#pragma acc routine seq
static inline void single_tamat_into_tamat_soa(__restrict tamat_soa * out, int idx, single_tamat * in){

    out->ic00[idx] = in->ic00;
    out->ic11[idx] = in->ic11;
    out->c01[idx]  = in->c01  ;
    out->c02[idx]  = in->c02  ;
    out->c12[idx]  = in->c12  ;

}
#pragma acc routine seq
static inline void single_thmat_into_thmat_soa(__restrict thmat_soa * out, int idx, single_thmat * in){

    out->rc00[idx] = in->rc00;
    out->rc11[idx] = in->rc11;
    out->c01[idx]  = in->c01  ;
    out->c02[idx]  = in->c02  ;
    out->c12[idx]  = in->c12  ;

}


//misc. extraction and conversion
#pragma acc routine seq
static inline void i_times_tamat_soa_to_su3(single_su3 * out, __restrict tamat_soa * const QA,const int idx){
  //I*tamat is a thmat
    out->comp[0][0] =           - QA->ic00[idx] ; 
    out->comp[0][1] =   (1.0*I) * QA->c01[idx]  ;
    out->comp[0][2] =   (1.0*I) * QA->c02[idx]  ;

    out->comp[1][0] =   (-1.0*I) * conj(QA->c01[idx]) ;
    out->comp[1][1] =           - QA->ic11[idx] ; 
    out->comp[1][2] =   (1.0*I) * QA->c12[idx]  ;

    out->comp[2][0] =   (-1.0*I) * conj(QA->c02[idx]) ;
    out->comp[2][1] =   (-1.0*I) * conj(QA->c12[idx]) ;
    out->comp[2][2] =   QA->ic00[idx] + QA->ic11[idx]; 


}


// functions
#pragma acc routine seq
static inline void herm_conj_gl3(single_su3 * in, single_su3 *out){

    out->comp[0][0] = conj(in->comp[0][0]);
    out->comp[0][1] = conj(in->comp[1][0]);
    out->comp[0][2] = conj(in->comp[2][0]);
    out->comp[1][0] = conj(in->comp[0][1]);
    out->comp[1][1] = conj(in->comp[1][1]);
    out->comp[1][2] = conj(in->comp[2][1]);
    out->comp[2][0] = conj(in->comp[0][2]);
    out->comp[2][1] = conj(in->comp[1][2]);
    out->comp[2][2] = conj(in->comp[2][2]);

}
#pragma acc routine seq
static inline double detQ(single_thmat *Q){

    double rc22 = -Q->rc00-Q->rc11 ; 
    return creal(Q->rc00*Q->rc11*rc22 + 2*creal(Q->c01*Q->c12*conj(Q->c02)) -
        Q->rc00 * Q->c12 * conj(Q->c12) - rc22 * Q->c01 * conj(Q->c01) - 
        Q->rc11 * Q->c02 * conj(Q->c02));

}

#pragma acc routine seq
static inline double det_i_times_QA(__restrict single_tamat * const QA){

    double ic22 = -QA->ic00-QA->ic11 ; 

    return -creal(  QA->ic00 * QA->ic11 * ic22
		  + 2*cimag(QA->c01 * QA->c12 * conj(QA->c02))
		  - QA->ic00 * QA->c12 * conj(QA->c12)
		  - QA->ic11 * QA->c02 * conj(QA->c02)
      	          - ic22 * QA->c01 * conj(QA->c01) );

}
#pragma acc routine seq
static inline void print_su3_stdout(single_su3 *m){

    printf("\n");
   for(int r=0;r<3;r++){
    for(int c=0;c<3;c++) printf("%.18lf + %.18lf I   ",  creal(m->comp[r][c]), cimag(m->comp[r][c]));
    printf("\n");

   }
}

#pragma acc routine seq
static inline void summ_single_tamats_times_scalar(single_tamat * out, single_tamat * QA1 , single_tamat *QA2, d_complex scalar){
  out->c01 = (QA1->c01 + QA2->c01)*scalar; // comp_01
  out->c02 = (QA1->c02 + QA2->c02)*scalar; // comp_02
  out->c12 = (QA1->c12 + QA2->c12)*scalar; // comp_12
  out->ic00 = (QA1->ic00+ QA2->ic00)*creal(scalar);   // Im(comp_00)
  out->ic11 = (QA1->ic11+ QA1->ic11)*creal(scalar);   // Im(comp_11)
}

#pragma acc routine seq
static inline void single_tamat_times_scalar_into_tamat(single_tamat * out, single_tamat * QA , d_complex scalar){
  out->c01 = QA->c01 *scalar; // comp_01
  out->c02 = QA->c02 *scalar; // comp_02
  out->c12 = QA->c12 *scalar; // comp_12
  out->ic00 = QA->ic00*creal(scalar);   // Im(comp_00)
  out->ic11 = QA->ic11*creal(scalar);   // Im(comp_11)
}

#pragma acc routine seq
static inline void single_tamat_times_scalar_add_to_tamat(single_tamat * out, single_tamat * QA , d_complex scalar){
  out->c01 += QA->c01 *scalar; // comp_01
  out->c02 += QA->c02 *scalar; // comp_02
  out->c12 += QA->c12 *scalar; // comp_12
  out->ic00 += QA->ic00*creal(scalar);   // Im(comp_00)
  out->ic11 += QA->ic11*creal(scalar);   // Im(comp_11)
}


#pragma acc routine seq
static inline void single_su3_times_scalar(single_su3 * m , d_complex scalar){

   for(int r=0;r<3;r++)
    for(int c=0;c<3;c++)
     m->comp[r][c] *= scalar;
}

#pragma acc routine seq
static inline void single_su3xsu3(single_su3 * out , single_su3 *m1, single_su3 *m2){

   for(int r=0;r<3;r++)
    for(int c=0;c<3;c++){
        out->comp[r][c] = 0;
        for(int d=0;d<3;d++) out->comp[r][c] += m1->comp[r][d] * m2->comp[d][c] ;

    }
}
#pragma acc routine seq
static inline void single_su3xsu3_add_to_out(single_su3 * out , single_su3 *m1, single_su3 *m2){

   for(int r=0;r<3;r++)
    for(int c=0;c<3;c++){
        //   out->comp[r][c] = 0; // add to out!
        for(int d=0;d<3;d++) out->comp[r][c] += m1->comp[r][d] * m2->comp[d][c] ;

    }
}

#pragma acc routine seq
static inline void single_gl3xsu3_add_to_out(single_su3 * out , single_su3 *m1, single_su3 *m2){
  single_su3 *tmp=m2;
  rebuild3row(tmp);
  for(int r=0;r<3;r++)
    for(int c=0;c<3;c++){
      for(int d=0;d<3;d++) out->comp[r][c] += m1->comp[r][d] * tmp->comp[d][c] ;
    }

}

#pragma acc routine seq
static inline void single_su3xgl3_add_to_out(single_su3 * out , single_su3 *m1, single_su3 *m2){
  single_su3 *tmp=m1;
  rebuild3row(tmp);
   for(int r=0;r<3;r++)
    for(int c=0;c<3;c++){
        for(int d=0;d<3;d++) out->comp[r][c] += tmp->comp[r][d] * m2->comp[d][c] ;

    }
}

#pragma acc routine seq
static inline void single_su3add(single_su3 * out , single_su3 *m){

   for(int r=0;r<3;r++)// Magari fino alla seconda riga?
   //for(int r=0;r<2;r++) //??
    for(int c=0;c<3;c++)
        out->comp[r][c] += m->comp[r][c];

}
#pragma acc routine seq
static inline double TrQsq(single_thmat *Q){

  return 2 * creal( Q->rc00 * Q->rc00 +
                    Q->rc11 * Q->rc11 +
                    Q->c01 * conj(Q->c01) +
                    Q->c02 * conj(Q->c02) +
                    Q->c12 * conj(Q->c12) +
	            Q->rc00 * Q->rc11 );

}
#pragma acc routine seq
static inline d_complex detSu3(single_su3 *m){


    return  
    m->comp[0][0]* m->comp[1][1] * m->comp[2][2] +
    m->comp[0][1]* m->comp[1][2] * m->comp[2][0] +
    m->comp[0][2]* m->comp[1][0] * m->comp[2][1] -
    m->comp[0][0]* m->comp[1][2] * m->comp[2][1] -
    m->comp[0][1]* m->comp[1][0] * m->comp[2][2] -
    m->comp[0][2]* m->comp[1][1] * m->comp[2][0] ;

}
#pragma acc routine seq
static inline double Tr_i_times_QA_sq(__restrict single_tamat * const QA){
  // computes Tr( (i*QA)^2 )
    return 2 * creal( QA->ic00 * QA->ic00 +
		      QA->ic11 * QA->ic11 +
		      QA->ic00 * QA->ic11 +
		      QA->c01  * conj(QA->c01) +
		      QA->c02  * conj(QA->c02) +
		      QA->c12  * conj(QA->c12) );
}
#pragma acc routine seq
static inline void single_thmatAeqAmB(single_thmat *A, single_thmat* B){

    A->rc00 -= B->rc00;
    A->rc11 -= B->rc11;
    A->c01 -= B->c01;
    A->c02 -= B->c02;
    A->c12 -= B->c12;

}

#pragma acc routine seq
static inline void gl3_dagger(single_su3 * inout){

    inout->comp[0][0] = conj(inout->comp[0][0]);
    inout->comp[1][1] = conj(inout->comp[1][1]);
    inout->comp[2][2] = conj(inout->comp[2][2]);

    d_complex temp;

    temp = inout->comp[0][1]; inout->comp[0][1] = conj(inout->comp[1][0]);
    inout->comp[1][0] = conj(temp);

    temp = inout->comp[0][2]; inout->comp[0][2] = conj(inout->comp[2][0]);
    inout->comp[2][0] = conj(temp);

    temp = inout->comp[2][1]; inout->comp[2][1] = conj(inout->comp[1][2]);
    inout->comp[1][2] = conj(temp);


}
#pragma acc routine seq
static inline void single_su3_copy(const single_su3 * const in, single_su3 * out){

    for(int i =0;i<3;i++)    
        for(int j =0;j<3;j++)    
            out->comp[i][j] = in->comp[i][j];

}


//no 3rd row
#pragma acc routine seq
static inline void single_su3xsu3_no3rdrow(single_su3 * out , single_su3 *m1, single_su3 *m2){

   for(int r=0;r<2;r++)
    for(int c=0;c<3;c++){
        out->comp[r][c] = 0;
        for(int d=0;d<3;d++) out->comp[r][c] += m1->comp[r][d] * m2->comp[d][c] ;

    }
}
#pragma acc routine seq
static inline void single_su3_times_scalar_no3rdrow(single_su3 * m , d_complex scalar){

   for(int r=0;r<2;r++)
    for(int c=0;c<3;c++)
     m->comp[r][c] *= scalar;

}
#pragma acc routine seq
static inline void single_su3add_no3rdrow(single_su3 * out , single_su3 *m){

   for(int r=0;r<2;r++)
    for(int c=0;c<3;c++)
        out->comp[r][c] += m->comp[r][c];

}

#pragma acc routine seq
static inline void su3_soa_times_su3_soa_into_single_su3( const su3_soa * const mat1, const int idx_mat1, const su3_soa * const mat2, const int idx_mat2, single_su3 * Omat){

  //first single_su3
  single_su3 m1;
  m1.comp[0][0] = mat1->r0.c0[idx_mat1];
  m1.comp[0][1] = mat1->r0.c1[idx_mat1];
  m1.comp[0][2] = mat1->r0.c2[idx_mat1];
  m1.comp[1][0] = mat1->r1.c0[idx_mat1];
  m1.comp[1][1] = mat1->r1.c1[idx_mat1];
  m1.comp[1][2] = mat1->r1.c2[idx_mat1];
  //second single_su3
  single_su3 m2;
  m2.comp[0][0] = mat2->r0.c0[idx_mat2];
  m2.comp[0][1] = mat2->r0.c1[idx_mat2];
  m2.comp[0][2] = mat2->r0.c2[idx_mat2];
  m2.comp[1][0] = mat2->r1.c0[idx_mat2];
  m2.comp[1][1] = mat2->r1.c1[idx_mat2];
  m2.comp[1][2] = mat2->r1.c2[idx_mat2];
  rebuild3row(&m2);
  for(int r = 0; r < 2; r++)
    for(int c = 0; c < 3; c++){
      Omat->comp[r][c] = 0;
      for (int k = 0; k<3; k++)
	Omat->comp[r][c] += m1.comp[r][k]*m2.comp[k][c];
    }
}

#pragma acc routine seq
static inline void su3_soa_times_su3_soa_dag_into_single_su3( __restrict const su3_soa * const mat1, const int idx_mat1, __restrict const su3_soa * const mat2, const int idx_mat2, single_su3 * Omat){

  //first single_su3
  single_su3 m1;
  m1.comp[0][0] = mat1->r0.c0[idx_mat1];
  m1.comp[0][1] = mat1->r0.c1[idx_mat1];
  m1.comp[0][2] = mat1->r0.c2[idx_mat1];
  m1.comp[1][0] = mat1->r1.c0[idx_mat1];
  m1.comp[1][1] = mat1->r1.c1[idx_mat1];
  m1.comp[1][2] = mat1->r1.c2[idx_mat1];
  //second single_su3
  single_su3 m2;
  single_su3_from_su3_soa(mat2,idx_mat2,&m2);
  rebuild3row(&m2);
  gl3_dagger(&m2);

  for(int r = 0; r < 2; r++)
    for(int c = 0; c < 3; c++){
      Omat->comp[r][c] = 0;
      for (int k = 0; k<3; k++)
	Omat->comp[r][c] += m1.comp[r][k]*m2.comp[k][c];
    }
}

#pragma acc routine seq
static inline void su3_soa_dag_times_su3_soa_into_single_su3( __restrict const su3_soa * const mat1, const int idx_mat1, __restrict const su3_soa * const mat2, const int idx_mat2, single_su3 * Omat){

  //first single_su3
  single_su3 m1;
  single_su3_from_su3_soa(mat1,idx_mat1,&m1);
  rebuild3row(&m1);
  gl3_dagger(&m1);
  //second single_su3
  single_su3 m2;
  single_su3_from_su3_soa(mat2,idx_mat2,&m2);
  rebuild3row(&m2);
  
  for(int r = 0; r < 2; r++)
    for(int c = 0; c < 3; c++){
      Omat->comp[r][c] = 0;
      for (int k = 0; k<3; k++)
	Omat->comp[r][c] += m1.comp[r][k]*m2.comp[k][c];
    }
}

#pragma acc routine seq
static inline void su3_soa_times_gl3_soa_into_gl3(__restrict const su3_soa * const mat1, const int idx_mat1, __restrict const su3_soa * const mat2, const int idx_mat2, single_su3 * Omat)
{
  //single su3
  single_su3 su3;
  single_su3_from_su3_soa(mat1,idx_mat1,&su3);
  rebuild3row(&su3);
  //single gl3
  single_su3 gl3;
  single_su3_from_su3_soa(mat2,idx_mat2,&gl3);
  gl3.comp[2][0] = mat2->r2.c0[idx_mat2];
  gl3.comp[2][1] = mat2->r2.c1[idx_mat2];
  gl3.comp[2][2] = mat2->r2.c2[idx_mat2];

  single_su3xsu3(Omat,&su3,&gl3);
}

#pragma acc routine seq
static inline void gl3_soa_times_su3_soa_dag_into_gl3(__restrict const su3_soa * const mat1, const int idx_mat1, __restrict const su3_soa * const mat2, const int idx_mat2, single_su3 * Omat)
{
  //single gl3
  single_su3 gl3;
  single_su3_from_su3_soa(mat1,idx_mat1,&gl3);
  gl3.comp[2][0] = mat1->r2.c0[idx_mat1];
  gl3.comp[2][1] = mat1->r2.c1[idx_mat1];
  gl3.comp[2][2] = mat1->r2.c2[idx_mat1];

  //single su3
  single_su3 su3;
  single_su3_from_su3_soa(mat2,idx_mat2,&su3);
  rebuild3row(&su3);
  gl3_dagger(&su3);
  
  single_su3xsu3(Omat,&gl3,&su3);
}

#pragma acc routine seq
static inline void su3_soa_dag_times_gl3_soa_dag_into_gl3(__restrict const su3_soa * const mat1, const int idx_mat1, __restrict const su3_soa * const mat2, const int idx_mat2, single_su3 * Omat)
{
  //single su3
  single_su3 su3;
  single_su3_from_su3_soa(mat1,idx_mat1,&su3);
  rebuild3row(&su3);
  gl3_dagger(&su3);
  
  //single gl3
  single_su3 gl3;
  single_su3_from_su3_soa(mat2,idx_mat2,&gl3);
  gl3.comp[2][0] = mat2->r2.c0[idx_mat2];
  gl3.comp[2][1] = mat2->r2.c1[idx_mat2];
  gl3.comp[2][2] = mat2->r2.c2[idx_mat2];
  gl3_dagger(&gl3);
  
  single_su3xsu3(Omat,&su3,&gl3);
}

#pragma acc routine seq
static inline void gl3_soa_dag_times_su3_soa_into_gl3(__restrict const su3_soa * const mat1, const int idx_mat1, __restrict const su3_soa * const mat2, const int idx_mat2, single_su3 * Omat)
{
  //single gl3
  single_su3 gl3;
  single_su3_from_su3_soa(mat1,idx_mat1,&gl3);
  gl3.comp[2][0] = mat1->r2.c0[idx_mat1];
  gl3.comp[2][1] = mat1->r2.c1[idx_mat1];
  gl3.comp[2][2] = mat1->r2.c2[idx_mat1];
  gl3_dagger(&gl3);
  
  //single su3 - MAKE SURE IT IS AN SU3
  single_su3 su3;
  single_su3_from_su3_soa(mat2,idx_mat2,&su3);
  rebuild3row(&su3);
  
  single_su3xsu3(Omat,&gl3,&su3);
}

#pragma acc routine seq
static inline void gl3_soa_times_single_su3_addto_gl3(__restrict const su3_soa * const mat1, const int idx_mat1, single_su3 * mat2, single_su3 * Omat)
{
  //single gl3
  single_su3 gl3;
  single_su3_from_su3_soa(mat1,idx_mat1,&gl3);
  gl3.comp[2][0] = mat1->r2.c0[idx_mat1];
  gl3.comp[2][1] = mat1->r2.c1[idx_mat1];
  gl3.comp[2][2] = mat1->r2.c2[idx_mat1];
  
  //single su3 - MAKE SURE IT IS AN SU3
  single_su3 mat2_3row;
  mat2_3row.comp[0][0] = mat2->comp[0][0];
  mat2_3row.comp[0][1] = mat2->comp[0][1];
  mat2_3row.comp[0][2] = mat2->comp[0][2];

  mat2_3row.comp[1][0] = mat2->comp[1][0];
  mat2_3row.comp[1][1] = mat2->comp[1][1];
  mat2_3row.comp[1][2] = mat2->comp[1][2];
  rebuild3row(&mat2_3row);
  
  for(int r=0;r<3;r++)
    for(int c=0;c<3;c++)
      for(int d=0;d<3;d++)
	Omat->comp[r][c] += gl3.comp[r][d] * mat2_3row.comp[d][c] ;

}

#pragma acc routine seq
static inline void single_su3_times_gl3_soa_addto_gl3(single_su3 * mat1,__restrict const su3_soa * const mat2, const int idx_mat2,  single_su3 * Omat)
{
  //single su3 - MAKE SURE IT IS AN SU3
  single_su3 mat1_3row;
  mat1_3row.comp[0][0] = mat1->comp[0][0];
  mat1_3row.comp[0][1] = mat1->comp[0][1];
  mat1_3row.comp[0][2] = mat1->comp[0][2];

  mat1_3row.comp[1][0] = mat1->comp[1][0];
  mat1_3row.comp[1][1] = mat1->comp[1][1];
  mat1_3row.comp[1][2] = mat1->comp[1][2];
  rebuild3row(&mat1_3row);
  //single gl3
  single_su3 gl3;
  single_su3_from_su3_soa(mat2,idx_mat2,&gl3);
  gl3.comp[2][0] = mat2->r2.c0[idx_mat2];
  gl3.comp[2][1] = mat2->r2.c1[idx_mat2];
  gl3.comp[2][2] = mat2->r2.c2[idx_mat2];

  for(int r=0;r<3;r++)
    for(int c=0;c<3;c++)
      for(int d=0;d<3;d++)
	Omat->comp[r][c] += mat1_3row.comp[r][d] * gl3.comp[d][c] ;
  
}

#pragma acc routine seq
static inline void gl3_soa_dag_times_single_su3_addto_gl3(__restrict const su3_soa * const mat1, const int idx_mat1, single_su3 * mat2, single_su3 * Omat)
{
  //single gl3
  single_su3 gl3;
  single_su3_from_su3_soa(mat1,idx_mat1,&gl3);
  gl3.comp[2][0] = mat1->r2.c0[idx_mat1];
  gl3.comp[2][1] = mat1->r2.c1[idx_mat1];
  gl3.comp[2][2] = mat1->r2.c2[idx_mat1];
  gl3_dagger(&gl3);
  
  rebuild3row(mat2);
  
  single_su3xsu3_add_to_out(Omat,&gl3,mat2);
}

#pragma acc routine seq
static inline void single_su3_times_gl3_soa_dag_addto_gl3(single_su3 * mat1,__restrict const su3_soa * const mat2, const int idx_mat2,  single_su3 * Omat)
{
  //single su3 - MAKE SURE IT IS AN SU3
  single_su3 mat1_3row;
  mat1_3row.comp[0][0] = mat1->comp[0][0];
  mat1_3row.comp[0][1] = mat1->comp[0][1];
  mat1_3row.comp[0][2] = mat1->comp[0][2];

  mat1_3row.comp[1][0] = mat1->comp[1][0];
  mat1_3row.comp[1][1] = mat1->comp[1][1];
  mat1_3row.comp[1][2] = mat1->comp[1][2];
  rebuild3row(&mat1_3row);
  //single gl3
  single_su3 gl3;
  single_su3_from_su3_soa(mat2,idx_mat2,&gl3);
  gl3.comp[2][0] = mat2->r2.c0[idx_mat2];
  gl3.comp[2][1] = mat2->r2.c1[idx_mat2];
  gl3.comp[2][2] = mat2->r2.c2[idx_mat2];
  gl3_dagger(&gl3);
  
  single_su3xsu3_add_to_out(Omat,&mat1_3row,&gl3);
}

#pragma acc routine seq
static inline void su3_soa_times_single_su3_into_single_su3( __restrict const su3_soa * const mat1, const int idx_mat1, single_su3 * mat2, single_su3 * Omat){

  //first single_su3
  single_su3 m1;
  single_su3_from_su3_soa(mat1,idx_mat1,&m1);
  
  single_su3 m2;
  m2.comp[0][0] = mat2->comp[0][0];
  m2.comp[0][1] = mat2->comp[0][1];
  m2.comp[0][2] = mat2->comp[0][2];

  m2.comp[1][0] = mat2->comp[1][0];
  m2.comp[1][1] = mat2->comp[1][1];
  m2.comp[1][2] = mat2->comp[1][2];
  rebuild3row(&m2);
  
  single_su3xsu3_no3rdrow(Omat,&m1,&m2);
}

#pragma acc routine seq
static inline void single_su3_times_su3_soa_into_single_su3( single_su3 * mat1, __restrict const su3_soa * const mat2, const int idx_mat2, single_su3 * Omat){

  //second single_su3
  single_su3 m2;
  single_su3_from_su3_soa(mat2,idx_mat2,&m2);
  rebuild3row(&m2);

  single_su3xsu3_no3rdrow(Omat,mat1,&m2);
}

#pragma acc routine seq
static inline void set_to_zero_single_su3(single_su3 * Inout){
  Inout->comp[0][0] = 0+I*0;
  Inout->comp[0][1] = 0+I*0;
  Inout->comp[0][2] = 0+I*0;
  Inout->comp[1][0] = 0+I*0;
  Inout->comp[1][1] = 0+I*0;
  Inout->comp[1][2] = 0+I*0;
  Inout->comp[2][0] = 0+I*0;
  Inout->comp[2][1] = 0+I*0;
  Inout->comp[2][2] = 0+I*0;
}

/*
#pragma acc routine seq
static inline void su3_soa_times_su3_soa_into_single_gl3( __restrict const su3_soa * const mat1, const int idx_mat1, __restrict const su3_soa * const mat2, const int idx_mat2, single_su3 * Omat){

  //first matrix elements IT MUST BE SU3
  m1_00 = mat1->r0.c0[idx_mat1];
  m1_01 = mat1->r0.c1[idx_mat1];
  m1_02 = mat1->r0.c2[idx_mat1];

  m1_10 = mat1->r1.c0[idx_mat1];
  m1_11 = mat1->r1.c1[idx_mat1];
  m1_12 = mat1->r1.c2[idx_mat1];  

  //second matrix elements IT MUST BE SU3
  m2_00 = mat2->r0.c0[idx_mat2];
  m2_01 = mat2->r0.c1[idx_mat2];
  m2_02 = mat2->r0.c2[idx_mat2];  

  m2_10 = mat2->r1.c0[idx_mat2];
  m2_11 = mat2->r1.c1[idx_mat2];
  m2_12 = mat2->r1.c2[idx_mat2];  

  m2_20 = conj(mat2->r0.c1[idx_mat2] * mat2->r1.c2[idx_mat2] - mat2->r0.c2[idx_mat2] * mat2->r1.c1[idx_mat2]);
  m2_21 = conj(mat2->r0.c2[idx_mat2] * mat2->r1.c0[idx_mat2] - mat2->r0.c0[idx_mat2] * mat2->r1.c2[idx_mat2]);
  m2_22 = conj(mat2->r0.c0[idx_mat2] * mat2->r1.c1[idx_mat2] - mat2->r0.c1[idx_mat2] * mat2->r1.c0[idx_mat2]);  
  
  //computation
  Omat->comp[0][0] = m1_00*m2_00 + m1_01*m2_10 + m1_02*m2_20;
  Omat->comp[0][1] = m1_00*m2_01 + m1_01*m2_11 + m1_02*m2_21;
  Omat->comp[0][2] = m1_00*m2_02 + m1_01*m2_12 + m1_02*m2_22;

  Omat->comp[1][0] = m1_10*m2_00 + m1_11*m2_10 + m1_12*m2_20;
  Omat->comp[1][1] = m1_10*m2_01 + m1_11*m2_11 + m1_12*m2_21;
  Omat->comp[1][2] = m1_10*m2_02 + m1_11*m2_12 + m1_12*m2_22;

  Omat->comp[2][0] = m1_20*m2_00 + m1_21*m2_10 + m1_22*m2_20;
  Omat->comp[2][1] = m1_20*m2_01 + m1_21*m2_11 + m1_22*m2_21;
  Omat->comp[2][2] = m1_20*m2_02 + m1_21*m2_12 + m1_22*m2_22;

}
*/


#endif
