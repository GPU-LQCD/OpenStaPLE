#ifndef DEOTT_STOUTING_C
#define DEOTT_STOUTING_C

#include "../OpenAcc/cayley_hamilton.h"
#include "../OpenAcc/struct_c_def.h"
#include "../OpenAcc/single_types.h"


static inline void DEOTT_conf_left_exp_multiply_to_su3_soa(__restrict su3_soa * const cnf,
							   const int idx,
							   __restrict su3_soa * const  EXP,
							   __restrict su3_soa * const cnf_out){
  
  
  single_su3 AUX;
  //Multiply: U_new = EXP * U_old
  
  //Extraction
  single_gl3_from_su3_soa(cnf, idx, &AUX);
  
  rebuild3row(&AUX);
  
  
  // moltiplica
  cnf_out->r0.c0[idx] = EXP->r0.c0[idx] * AUX.comp[0][0] + EXP->r0.c1[idx] * AUX.comp[1][0] + EXP->r0.c2[idx] * AUX.comp[2][0];
  cnf_out->r0.c1[idx] = EXP->r0.c0[idx] * AUX.comp[0][1] + EXP->r0.c1[idx] * AUX.comp[1][1] + EXP->r0.c2[idx] * AUX.comp[2][1];
  cnf_out->r0.c2[idx] = EXP->r0.c0[idx] * AUX.comp[0][2] + EXP->r0.c1[idx] * AUX.comp[1][2] + EXP->r0.c2[idx] * AUX.comp[2][2];
  
  cnf_out->r1.c0[idx] = EXP->r1.c0[idx] * AUX.comp[0][0] + EXP->r1.c1[idx] * AUX.comp[1][0] + EXP->r1.c2[idx] * AUX.comp[2][0];
  cnf_out->r1.c1[idx] = EXP->r1.c0[idx] * AUX.comp[0][1] + EXP->r1.c1[idx] * AUX.comp[1][1] + EXP->r1.c2[idx] * AUX.comp[2][1];
  cnf_out->r1.c2[idx] = EXP->r1.c0[idx] * AUX.comp[0][2] + EXP->r1.c1[idx] * AUX.comp[1][2] + EXP->r1.c2[idx] * AUX.comp[2][2];

  
}


void DEOTT_exp_minus_QA_times_conf(__restrict su3_soa * const tu,
        __restrict tamat_soa * const QA,
        __restrict su3_soa * const tu_out,
        __restrict su3_soa * const exp_aux){

    SETINUSE(exp_aux);
    SETINUSE(tu_out);
    int x, y, z, t;
#pragma acc kernels present(tu) present(QA) present(tu_out) present(exp_aux)
#pragma acc loop independent gang
    for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector
        for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector
            for(y=0; y<ny; y++) {
#pragma acc loop independent vector
                for(x=0; x < nx; x++) {
                    const  int idxh   = snum_acc(x,y,z,t);  // r
                    const  int parity = (x+y+z+t) % 2;
                    int dir_link;
                    int mu;

                    for(mu=0;mu<4;mu++){

                        dir_link = 2*mu + parity;
                        CH_exponential_antihermitian_soa_nissalike(&exp_aux[dir_link],&QA[dir_link],idxh);

                        DEOTT_conf_left_exp_multiply_to_su3_soa(&tu[dir_link],idxh,&exp_aux[dir_link],&tu_out[dir_link]);
                    }

                }  // x
            }  // y
        }  // z
    }  // t

    SETFREE(exp_aux);
}// closes routine



//calcolo di lambda
#pragma acc routine seq
static inline void DEOTT_compute_loc_Lambda(__restrict thmat_soa * const L, // la Lambda --> ouput
        __restrict su3_soa   * const SP, // Sigma primo --> input
        __restrict su3_soa   * const U,    // la configurazione di gauge --> input
        __restrict tamat_soa * const QA, // gli stessi Q che arrivano a Cayley hamilton --> input
        __restrict su3_soa   * const TMP,  // variabile di parcheggio
        int idx
        ){

    double c0 = det_i_times_QA_soa(QA,idx);
    double c1  = 0.5 * Tr_i_times_QA_sq_soa(QA,idx);
    double c0max = 2*pow(c1/3,1.5);
    d_complex f0,f1,f2;
    d_complex b10,b11,b12,b20,b21,b22;
    d_complex  r0_1,r1_1,r2_1,r0_2,r1_2,r2_2;
    if(c1<4e-3)
    {
        f0 = (1-c0*c0/720) + (1.0*I)*(-c0*(1-c1*(1-c1/42)/20)/6);
        f1 = (c0*(1-c1*(1-3*c1/112)/15)/24) + (1.0*I)*(1-c1*(1-c1*(1-c1/42)/20)/6-c0*c0/5040);
        f2 = (0.5*(-1+c1*(1-c1*(1-c1/56)/30)/12+c0*c0/20160)) + (1.0*I)*(0.5*(c0*(1-c1*(1-c1/48)/21)/60));

        //from nissa
        /* 
           b[0][0][RE]=0;
           b[0][0][IM]=c0/120*(1-c1/21);
           b[0][1][RE]=-c0/360*(1-3*c1/56);
           b[0][1][IM]=-1.0/6*(1-c1/10*(1.0-c1/28));
           b[0][2][RE]=0.5*(1.0/12*(1-2*c1/30*(1-3*c1/112)));
           b[0][2][IM]=0.5*(-c0/1260*(1-c1/24));

           b[1][0][RE]=-c0/360;
           b[1][0][IM]=-1.0/6*(1-c1/20*(1-c1/42));
           b[1][1][RE]=1.0/24*(1-c1/15*(1-3*c1/112));
           b[1][1][IM]=-c0/2520;
           b[1][2][RE]=0.5*c0/10080;
           b[1][2][IM]=0.5*(1.0/60*(1-c1/21*(1-c1/48)));
           */

        b10 = 0.0 + (1.0*I)*(c0/120*(1-c1/21));
        b11 = (-c0/360*(1-3*c1/56)) + (1.0*I)*(-1.0/6*(1-c1/10*(1.0-c1/28)));
        b12 = (0.5*(1.0/12*(1-2*c1/30*(1-3*c1/112)))) + (1.0*I)*(0.5*(-c0/1260*(1-c1/24)));

        b20 = (-c0/360) + (1.0*I)*(-1.0/6*(1-c1/20*(1-c1/42)));
        b21 = (1.0/24*(1-c1/15*(1-3*c1/112))) + (1.0*I)*(-c0/2520);
        b22 = (0.5*c0/10080) + (1.0*I)*(0.5*(1.0/60*(1-c1/21*(1-c1/48))));

    }
    else
    {
        int segno=1;
        if(c0<0)
        {
            segno=-1;
            c0=-c0;
        }

        double eps=(c0max-c0)/c0max;
        double theta;
        if(eps<0) theta=0.0;
        else
            if(eps<1e-3) theta=sqrt(2*eps)*(1+(1.0/12+(3.0/160+(5.0/896+(35.0/18432+63.0/90112*eps)*eps)*eps)*eps)*eps);
            else theta=acos(c0/c0max);
        double u = sqrt(c1/3) * cos(theta/3) ;
        double w = sqrt(c1) * sin(theta/3) ;
        double u2=u*u,w2=w*w,u2mw2=u2-w2,w2p3u2=w2+3*u2,w2m3u2=w2-3*u2;
        double cu=cos(u),c2u=cos(2*u);
        double su=sin(u),s2u=sin(2*u);
        double cw=cos(w);
        double xi0w;
        if(fabs(w)<0.05)
        {
            double temp0=w*w,temp1=1-temp0/42,temp2=1.0-temp0/20*temp1;
            xi0w=1-temp0/6*temp2;
        }
        else xi0w=sin(w)/w;
        double xi1w; //eq. (67)
        if(fabs(w)<0.05) xi1w=-(1-w2*(1-w2*(1-w2/54)/28)/10)/3;
        else xi1w=cw/w2-sin(w)/(w2*w);

        double denom = 1/(9*u*u - w*w);

        f0=(u2mw2*c2u+ //(u2-w2)*cos(2u)
                cu*8*u2*cw+ //cos(u)*8*u2*cos(w)
                2*su*u*w2p3u2*xi0w) //sin(u)*2*mu*(3*u2+w2)*xi0(w)
                +(1.0*I)*(u2mw2*s2u+ //(u2-w2)*sin(2u)
                        -su*8*u2*cw+ //-sin(u)*8*u2*cos(w)
                        cu*2*u*w2p3u2*xi0w); //cos(u)*2*u*(3*u2+w2)*xi0(w)
        f0 *= denom;
        f1=(2*u*c2u+ //2*u*cos(2u)
                -cu*2*u*cw+ //cos(u)*2*u*cos(w)
                -su*w2m3u2*xi0w) //sin(u)*(u2-3*w2)*xi0(w)
                +(1.0*I)*(2*u*s2u+ //2*u*sin(2u)
                        su*2*u*cos(w)+ //sin(u)*2*u*cos(w)
                        -cu*w2m3u2*xi0w);//cos(u)*(3*u2-w2)*xi0(w)
        f1 *= denom;
        f2=(c2u+ //cos(2u)
                -cu*cw+ //-cos(u)*cos(w)
                -3*su*u*xi0w) //-3*sin(u)*u*xi0(w)
                +(1.0*I)*(s2u+ //sin(2u)
                        su*cw+ //sin(w)*cos(w)
                        -cu*3*u*xi0w);//-cos(u)*3*u*xi0(w)
        f2 *= denom;

        //from nissa
        /*
           complex r[2][3]=
           {{{2*c2u*u+s2u*(-2*u2+2*w2)+2*cu*u*(8*cw+3*u2*xi0w+w2*xi0w)+su*(-8*cw*u2+18*u2*xi0w+2*w2*xi0w),
           -8*cw*(2*su*u+cu*u2)+2*(s2u*u+c2u*u2-c2u*w2)+2*(9*cu*u2-3*su*u*u2+cu*w2 -su*u*w2)*xi0w},
           {2*c2u-4*s2u*u+su*(2*cw*u+6*u*xi0w)+cu*(-2*cw+3*u2*xi0w-w2*xi0w),
           2*s2u+4*c2u*u+2*cw*(su+cu*u)+(6*cu*u-3*su*u2+su*w2)*xi0w},
           {-2*s2u+cw*su-3*(su+cu*u)*xi0w,
           2*c2u+cu*cw+(-3*cu+3*su*u)*xi0w}},

           {{-2*c2u+2*cw*su*u+2*su*u*xi0w-8*cu*u2*xi0w+6*su*u*u2*xi1w,
           2*(-s2u+4*su*u2*xi0w+cu*u*(cw+xi0w+3*u2*xi1w))},
           {2*cu*u*xi0w+su*(-cw-xi0w+3*u2*xi1w),
           -2*su*u*xi0w-cu*(cw+xi0w-3*u2*xi1w)},
           {cu*xi0w-3*su*u*xi1w,
           -(su*xi0w)-3*cu*u*xi1w}}};
           */

        r0_1 = (2*c2u*u+s2u*(-2*u2+2*w2)+2*cu*u*(8*cw+3*u2*xi0w+w2*xi0w)+su*(-8*cw*u2+18*u2*xi0w+2*w2*xi0w))
            + (1.0*I)*(-8*cw*(2*su*u+cu*u2)+2*(s2u*u+c2u*u2-c2u*w2)+2*(9*cu*u2-3*su*u*u2+cu*w2 -su*u*w2)*xi0w);
        r1_1 = (2*c2u-4*s2u*u+su*(2*cw*u+6*u*xi0w)+cu*(-2*cw+3*u2*xi0w-w2*xi0w))
            + (1.0*I)*(2*s2u+4*c2u*u+2*cw*(su+cu*u)+(6*cu*u-3*su*u2+su*w2)*xi0w);
        r2_1 = (-2*s2u+cw*su-3*(su+cu*u)*xi0w)
            + (1.0*I)*(2*c2u+cu*cw+(-3*cu+3*su*u)*xi0w);

        r0_2 = (-2*c2u+2*cw*su*u+2*su*u*xi0w-8*cu*u2*xi0w+6*su*u*u2*xi1w)
            + (1.0*I)*(2*(-s2u+4*su*u2*xi0w+cu*u*(cw+xi0w+3*u2*xi1w)));
        r1_2 = (2*cu*u*xi0w+su*(-cw-xi0w+3*u2*xi1w))
            + (1.0*I)*(-2*su*u*xi0w-cu*(cw+xi0w-3*u2*xi1w));
        r2_2 = (cu*xi0w-3*su*u*xi1w)
            + (1.0*I)*(-(su*xi0w)-3*cu*u*xi1w);

        b10 = 0.5*denom*denom*(2.0*u*r0_1 + (3.0*u*u-w*w)*r0_2 -2.0*(15.0*u*u+w*w)*f0); // (57) 
        b11 = 0.5*denom*denom*(2.0*u*r1_1 + (3.0*u*u-w*w)*r1_2 -2.0*(15.0*u*u+w*w)*f1); // (57)
        b12 = 0.5*denom*denom*(2.0*u*r2_1 + (3.0*u*u-w*w)*r2_2 -2.0*(15.0*u*u+w*w)*f2); // (57)

        b20 = 0.5*denom*denom*(r0_1 - 3.0*u*r0_2 -24.0*u*f0); // (58)
        b21 = 0.5*denom*denom*(r1_1 - 3.0*u*r1_2 -24.0*u*f1); // (58)
        b22 = 0.5*denom*denom*(r2_1 - 3.0*u*r2_2 -24.0*u*f2); // (58)

        if(segno==-1){
            b10 = +conj(b10);
            b11 = -conj(b11);
            b12 = +conj(b12);
            b20 = -conj(b20);
            b21 = +conj(b21);
            b22 = -conj(b22);


            f0 =  conj(f0);
            f1 = -conj(f1);
            f2 =  conj(f2);
        }
    }


    single_su3 gl3_temp1,gl3_temp2;
    // estrazione di sQA,sU,sSP
    single_tamat sQA;
    single_tamat_from_tamat_soa(QA,idx,&sQA);

    single_su3 sU;
    single_su3_from_su3_soa(U,idx,&sU);
    rebuild3row(&sU);

    single_su3 sSP;
    single_gl3_from_su3_soa(SP,idx,&sSP);

    // Calcolo di B1   
    single_su3 B1; 
    Itamat_2ndDeg_poly(b10,b11,b12,&sQA,&B1);


    // CALCOLO DELLA TRACCIA => tr1 = Tr(Sigma' * B1 * U)

    single_su3xsu3(&gl3_temp1,&sU,&sSP);// gl3_temp1 = u*Sigma'
    single_su3xsu3(&gl3_temp2,&B1,&gl3_temp1);
    d_complex tr1= gl3_temp2.comp[2][2]+gl3_temp2.comp[1][1]+gl3_temp2.comp[0][0];

    // Calcolo di B2
    single_su3 B2; 
    Itamat_2ndDeg_poly(b20,b21,b22,&sQA,&B2);

    // CALCOLO DELLA TRACCIA => tr2 = Tr(Sigma' * B2 * U)
    single_su3xsu3(&gl3_temp2,&B2,&gl3_temp1);
    d_complex tr2= gl3_temp2.comp[2][2]+gl3_temp2.comp[1][1]+gl3_temp2.comp[0][0];

    ////////////////CALCOLO U*Sigma' e lo metto in TMP che non serve piu'
    // gia' fatto , e' in 'gl3_temp1'
    // single_su3xsu3(&gl3_temp1,&sU,&sSP);


    ///////////////// CALCOLO DI GAMMA = tr1*Q + tr2*Q^2 + f1* (U * Sigma') + f2 * (Q * U * Sigma' + U * Sigma' * Q ) = 
    /////////////////                  = tr1*Q + tr2*Q^2 + f1* TMP + f2 * (Q * TMP + TMP * Q)
    //                Soon to be:    --gl3_temp3---  -gl3_temp1-  ---gl3_temp2---
    single_su3 gl3_temp3; 
    i_times_tamat_to_su3(&gl3_temp3,&sQA);// getting Q= i QA from QA
    single_su3xsu3(&gl3_temp2,&gl3_temp3,&gl3_temp1);// gl3_temp2 = Q * U * S'
    single_su3xsu3_add_to_out(&gl3_temp2,&gl3_temp1,&gl3_temp3);
    //^^ gl3_temp2+=U*S'*Q

    single_su3_times_scalar(&gl3_temp2,f2);// f2 *(Q*U*S'+U*S'*Q)
    single_su3_times_scalar(&gl3_temp1,f1);// f1* (U * Sigma')
    Itamat_2ndDeg_poly(0,tr1,tr2,&sQA,&gl3_temp3);//gl3_temp3 = tr1*Q + tr2*Q^2 


    single_su3add(&gl3_temp1,&gl3_temp2);//f1*U*S' + f2*(Q*U*S'+U*S'*Q)
    single_su3add(&gl3_temp1,&gl3_temp3);// += tr1*Q + tr2*Q^2

    ///////////////// INFINE CALCOLO DI LAMBDA = 0.5*(GAMMA + GAMMA^CROCE) - (1/6)* Id * Tr(GAMMA + GAMMA^CROCE)
    single_thmat sLambda;
    gl3_to_thmat(&gl3_temp1,&sLambda);

    single_thmat_into_thmat_soa(L,idx,&sLambda);


       printf("c0 = %.18lf\n",c0);
       printf("c1 = %.18lf\n",c1);
       printf("f0 = %.18lf\n",f0);
       printf("f1 = %.18lf\n",f1);
       printf("f2 = %.18lf\n",f2);
//       printf("tr1 = (%.18lf) + (%.18lf)*I\n",creal(tr1),cimag(tr1));
//       printf("tr2 = (%.18lf) + (%.18lf)*I\n",creal(tr2),cimag(tr2));
       printf("b10 = (%.18lf) + (%.18lf)*I\n",creal(b10),cimag(b10));
       printf("b20 = (%.18lf) + (%.18lf)*I\n",creal(b20),cimag(b20));
       printf("b11 = (%.18lf) + (%.18lf)*I\n",creal(b11),cimag(b11));
       printf("b21 = (%.18lf) + (%.18lf)*I\n",creal(b21),cimag(b21));
       printf("b12 = (%.18lf) + (%.18lf)*I\n",creal(b12),cimag(b12));
       printf("b22 = (%.18lf) + (%.18lf)*I\n",creal(b22),cimag(b22));



/*
       printf("(U*SP)_00 = (%.18lf) + (%.18lf)*I \n",creal(TMP->r0.c0[idx]),cimag(TMP->r0.c0[idx]));
       printf("(U*SP)_01 = (%.18lf) + (%.18lf)*I \n",creal(TMP->r0.c1[idx]),cimag(TMP->r0.c1[idx]));
       printf("(U*SP)_02 = (%.18lf) + (%.18lf)*I \n",creal(TMP->r0.c2[idx]),cimag(TMP->r0.c2[idx]));
       printf("(U*SP)_10 = (%.18lf) + (%.18lf)*I \n",creal(TMP->r1.c0[idx]),cimag(TMP->r1.c0[idx]));
       printf("(U*SP)_11 = (%.18lf) + (%.18lf)*I \n",creal(TMP->r1.c1[idx]),cimag(TMP->r1.c1[idx]));
       printf("(U*SP)_12 = (%.18lf) + (%.18lf)*I \n",creal(TMP->r1.c2[idx]),cimag(TMP->r1.c2[idx]));
       printf("(U*SP)_20 = (%.18lf) + (%.18lf)*I \n",creal(TMP->r2.c0[idx]),cimag(TMP->r2.c0[idx]));
       printf("(U*SP)_21 = (%.18lf) + (%.18lf)*I \n",creal(TMP->r2.c1[idx]),cimag(TMP->r2.c1[idx]));
       printf("(U*SP)_22 = (%.18lf) + (%.18lf)*I \n",creal(TMP->r2.c2[idx]),cimag(TMP->r2.c2[idx]));

       printf("GAMMA_00  = (%.18lf) + (%.18lf)*I \n",creal(r0_1),cimag(r0_1));
       printf("GAMMA_01  = (%.18lf) + (%.18lf)*I \n",creal(r1_1),cimag(r1_1));
       printf("GAMMA_02  = (%.18lf) + (%.18lf)*I \n",creal(r2_1),cimag(r2_1));
       printf("GAMMA_10  = (%.18lf) + (%.18lf)*I \n",creal(r0_2),cimag(r0_2));
       printf("GAMMA_11  = (%.18lf) + (%.18lf)*I \n",creal(r1_2),cimag(r1_2));
       printf("GAMMA_12  = (%.18lf) + (%.18lf)*I \n",creal(r2_2),cimag(r2_2));
       printf("GAMMA_20  = (%.18lf) + (%.18lf)*I \n",creal(U20),cimag(U20));
       printf("GAMMA_21  = (%.18lf) + (%.18lf)*I \n",creal(U21),cimag(U21));
       printf("GAMMA_22  = (%.18lf) + (%.18lf)*I \n",creal(U22),cimag(U22));
*/

       printf("\n");


}



void DEOTT_compute_lambda(__restrict thmat_soa * const L, // la Lambda --> ouput  (una cosa che serve per calcolare la forza fermionica successiva)
        __restrict su3_soa   * const SP, // Sigma primo --> input (forza fermionica del passo precedente)
        __restrict su3_soa   * const U,    // la configurazione di gauge --> input
        __restrict tamat_soa * const QA, // gli stessi Q che arrivano a Cayley hamilton --> input (sostanzialmente sono rho*ta(staples))
        __restrict su3_soa   * const TMP  // variabile di parcheggio
        ){
  SETINUSE(L);
  SETINUSE(TMP);



    int x, y, z, t;
#pragma acc kernels present(L)  present(SP)  present(U)  present(QA)  present(TMP)
#pragma acc loop independent gang
    for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector
        for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector
            for(y=0; y<ny; y++) {
#pragma acc loop independent vector
                for(x=0; x < nx; x++) {
                    const  int idxh   = snum_acc(x,y,z,t);  // r
                    const  int parity = (x+y+z+t) % 2;
                    int dir_link;
                    int mu;
#pragma acc loop seq
                    for(mu=0;mu<4;mu++){
                        dir_link = 2*mu + parity;
                        DEOTT_compute_loc_Lambda(&L[dir_link],&SP[dir_link],&U[dir_link],&QA[dir_link],&TMP[dir_link],idxh);
                    }
                }  // x
            }  // y
        }  // z
    }  // t
}





#pragma acc routine seq
void DEOTT_compute_sigma_local_PEZZO1(__restrict thmat_soa * const L,  // la Lambda --> ouput  (una cosa che serve per calcolare la forza fermionica successiva)
        __restrict su3_soa   * const U,  // la configurazione di gauge --> input
        __restrict su3_soa   * const SP,  // entra Sigma primo (input: fermforce del passo precedente) ED esce Sigma --> sia input che ouput
        __restrict tamat_soa * const QA, // gli stessi Q che arrivano a Cayley hamilton --> input (sostanzialmente sono rho*ta(staples))
        __restrict su3_soa   * const TMP, // variabile di parcheggio
        int idx
        ){


    double c0 = det_i_times_QA_soa(QA,idx);
    double c1  = 0.5 * Tr_i_times_QA_sq_soa(QA,idx);
    double c0max = 2*pow(c1/3,1.5);
    d_complex f0,f1,f2;
    if(c1<4e-3)
    {
        f0 = (1-c0*c0/720) + (1.0*I)*(-c0*(1-c1*(1-c1/42)/20)/6);
        f1 = (c0*(1-c1*(1-3*c1/112)/15)/24) + (1.0*I)*(1-c1*(1-c1*(1-c1/42)/20)/6-c0*c0/5040);
        f2 = (0.5*(-1+c1*(1-c1*(1-c1/56)/30)/12+c0*c0/20160)) + (1.0*I)*(0.5*(c0*(1-c1*(1-c1/48)/21)/60));
    }
    else
    {
        int segno=1;
        if(c0<0)
        {
            segno=-1;
            c0=-c0;
        }

        double eps=(c0max-c0)/c0max;
        double theta;
        if(eps<0) theta=0.0;
        else
            if(eps<1e-3) theta=sqrt(2*eps)*(1+(1.0/12+(3.0/160+(5.0/896+(35.0/18432+63.0/90112*eps)*eps)*eps)*eps)*eps);
            else theta=acos(c0/c0max);
        double u = sqrt(c1/3) * cos(theta/3) ;
        double w = sqrt(c1) * sin(theta/3) ;
        double u2=u*u,w2=w*w,u2mw2=u2-w2,w2p3u2=w2+3*u2,w2m3u2=w2-3*u2;
        double cu=cos(u),c2u=cos(2*u);
        double su=sin(u),s2u=sin(2*u);
        double cw=cos(w);
        double xi0w;
        if(fabs(w)<0.05)
        {
            double temp0=w*w,temp1=1-temp0/42,temp2=1.0-temp0/20*temp1;
            xi0w=1-temp0/6*temp2;
        }
        else xi0w=sin(w)/w;

        double denom = 1/(9*u*u - w*w);

        f0=(u2mw2*c2u+ //(u2-w2)*cos(2u)
                cu*8*u2*cw+ //cos(u)*8*u2*cos(w)
                2*su*u*w2p3u2*xi0w) //sin(u)*2*mu*(3*u2+w2)*xi0(w)
                +(1.0*I)*(u2mw2*s2u+ //(u2-w2)*sin(2u)
                        -su*8*u2*cw+ //-sin(u)*8*u2*cos(w)
                        cu*2*u*w2p3u2*xi0w); //cos(u)*2*u*(3*u2+w2)*xi0(w)
        f0 *= denom;
        f1=(2*u*c2u+ //2*u*cos(2u)
                -cu*2*u*cw+ //cos(u)*2*u*cos(w)
                -su*w2m3u2*xi0w) //sin(u)*(u2-3*w2)*xi0(w)
                +(1.0*I)*(2*u*s2u+ //2*u*sin(2u)
                        su*2*u*cos(w)+ //sin(u)*2*u*cos(w)
                        -cu*w2m3u2*xi0w);//cos(u)*(3*u2-w2)*xi0(w)
        f1 *= denom;
        f2=(c2u+ //cos(2u)
                -cu*cw+ //-cos(u)*cos(w)
                -3*su*u*xi0w) //-3*sin(u)*u*xi0(w)
                +(1.0*I)*(s2u+ //sin(2u)
                        su*cw+ //sin(w)*cos(w)
                        -cu*3*u*xi0w);//-cos(u)*3*u*xi0(w)
        f2 *= denom;

        if(segno==-1){
            f0 =  conj(f0);
            f1 = -conj(f1);
            f2 =  conj(f2);
        }
    }

    ///   SIGMA = SIGMA'*EXP(iQ) + i * RHO * STAPLE(complete, non parte TA)*LAMBDA
    ///           - i * RHO * (  U Udag Udag LAMBDA + Udag Udag LAMBDA U + Udag LAMBDA Udag U
    ///                        - Udag Udag LAMBDA U - LAMBDA U Udag Udag - U Udag LAMBDA Udag ) =
    ///            PEZZO1 + PEZZO2 + PEZZO3

    //Calcolo il PEZZO 1 e lo metto in RES

    //////  CALCOLO TMP = EXP(iQ)

    single_tamat sQA;
    single_su3 sExpiQ,sSigmaP1,sSigmaPrimo;//sSigmaP1 : single Sigma Pezzo 1

    single_tamat_from_tamat_soa(QA,idx,&sQA);
    Itamat_2ndDeg_poly_no3rdrow(f0,f1,f2,&sQA,&sExpiQ);
    rebuild3row(&sExpiQ);

    // PEZZO2 e PEZZO3 sono gia' stati calcolati da un'altra routine (a meno del fatto di moltiplicare per RHO)
    // Adesso quindi faccio RES = RES * RHO + PEZZO1
    //////////  CALCOLO PEZZO1 ==>  SIGMA' * EXP(IQ) = SIGMA' * TMP

    single_gl3_from_su3_soa(SP,idx,&sSigmaPrimo);
    single_su3xsu3(&sSigmaP1,&sSigmaPrimo,&sExpiQ);
    single_gl3_into_su3_soa(SP,idx,&sSigmaP1);
}



#pragma acc routine seq
static inline void DEOTT_RIGHT_iABC_times_DminusE_absent_stag_phases(  __restrict su3_soa * const UA,
        const int idxA,
        __restrict su3_soa * const UB,
        const int idxB,
        __restrict su3_soa * const UC,
        const int idxC,
        __restrict thmat_soa * const LD,
        const int idxD,
        __restrict thmat_soa * const LE,
        const int idxE,
        __restrict su3_soa * const RES,
        const int idxRES){
    // Cosa calcoliamo in questa routine:
    //  RES += UA * dag(UB) * dag(UC) * ((RHO*I)*(LD - LE))

    single_su3 gl3_temp0,gl3_temp1,gl3_temp2;
    single_thmat thmat_temp0, thmat_temp1;

    single_thmat_from_thmat_soa(LD,idxD,&thmat_temp0);
    single_thmat_from_thmat_soa(LE,idxE,&thmat_temp1);
    single_thmatAeqAmB(&thmat_temp0,&thmat_temp1);// sLD = sLD-sLE

    thmat_to_su3(&gl3_temp0,&thmat_temp0);// to calculate products

    single_su3_from_su3_soa(UC,idxC,&gl3_temp1);//UC
    rebuild3row(&gl3_temp1);
    gl3_dagger(&gl3_temp1);//dag(UC)
    single_su3xsu3(&gl3_temp2,&gl3_temp1,&gl3_temp0); // dag(UC) * (LD-LE)
    // gl3_temp0  gl3_temp1  gl3_temp2
    //   LD-LE      dag(UC)  res   


    single_su3_from_su3_soa(UB,idxB,&gl3_temp0);//UB
    rebuild3row(&gl3_temp0);
    gl3_dagger(&gl3_temp0);//dag(UB)
    single_su3xsu3(&gl3_temp1,&gl3_temp0,&gl3_temp2);//dag(UB)*dag(UC)*(LD-LE)
    // gl3_temp0  gl3_temp1  gl3_temp2
    //  dag(UB)    res      dag(UC)*(LD-LE)

    single_su3_from_su3_soa(UA,idxA,&gl3_temp0);//UA
    rebuild3row(&gl3_temp0);

    //UA * dag(UB)*dag(UC)*(LD-LE)
    single_su3xsu3(&gl3_temp2,&gl3_temp0,&gl3_temp1);
    // gl3_temp0  gl3_temp1               gl3_temp2
    //  UA      dag(UB)*dag(UC)*(LD-LE)     res

    // multiplication by RHO*I
    single_su3_times_scalar(&gl3_temp2,RHO*I);

    //  RES += UA * dag(UB) * dag(UC) * ((RHO*I)*(LD - LE))
    single_gl3_addinto_su3_soa(RES,idxRES,&gl3_temp2);

}


// Questa e' di quelle della categoria RIGHT 
#pragma acc routine seq
static inline void DEOTT_RIGHT_iFABC_absent_stag_phases(  __restrict su3_soa * const UA,
        const int idxA,
        __restrict su3_soa * const UB,
        const int idxB,
        __restrict su3_soa * const UC,
        const int idxC,
        __restrict thmat_soa * const LF,
        const int idxF,
        __restrict su3_soa * const RES,
        const int idxRES){
    // Cosa calcoliamo in questa routine:
    //  RES += ((RHO*I)*LF) * UA * dag(UB) * dag(UC)

    single_su3 gl3_temp0,gl3_temp1,gl3_temp2;
    single_thmat thmat_temp0;

    single_su3_from_su3_soa(UC,idxC,&gl3_temp0);//UC
    rebuild3row(&gl3_temp0);
    gl3_dagger(&gl3_temp0);//dag(UC)

    single_su3_from_su3_soa(UB,idxB,&gl3_temp1);//UB
    rebuild3row(&gl3_temp1);
    gl3_dagger(&gl3_temp1);//dag(UB)


    single_su3xsu3(&gl3_temp2,&gl3_temp1,&gl3_temp0); // dag(UB) * dag(UC)
    // gl3_temp0  gl3_temp1  gl3_temp2
    //  dag(UC)    dag(UB)     res   

    single_su3_from_su3_soa(UA,idxA,&gl3_temp0);//UA
    rebuild3row(&gl3_temp0);

    //UA * dag(UB)*dag(UC)
    single_su3xsu3(&gl3_temp1,&gl3_temp0,&gl3_temp2);
    // gl3_temp0  gl3_temp1     gl3_temp2
    //  UA          res       dag(UB)*dag(UC)


    single_thmat_from_thmat_soa(LF,idxF,&thmat_temp0);// LF
    thmat_to_su3(&gl3_temp0,&thmat_temp0);// to calculate products

    single_su3xsu3(&gl3_temp2,&gl3_temp0,&gl3_temp1);
    // gl3_temp0    gl3_temp1        gl3_temp2
    //  LF      UA*dag(UB)*dag(UC)    res

    single_su3_times_scalar(&gl3_temp2, RHO*I);
    // ^^ ((RHO*I)* LF  * UA * dag(UB) * dag(UC)

    single_gl3_addinto_su3_soa(RES,idxRES,&gl3_temp2);
}



// Questa e' di quelle della categoria RIGHT 
#pragma acc routine seq
static inline void DEOTT_RIGHT_miABGC_absent_stag_phases(  __restrict su3_soa * const UA,
        const int idxA,
        __restrict su3_soa * const UB,
        const int idxB,
        __restrict su3_soa * const UC,
        const int idxC,
        __restrict thmat_soa * const LG,
        const int idxG,
        __restrict su3_soa * const RES,
        const int idxRES){
    // Cosa calcoliamo in questa routine:
    //  RES +=  UA * dag(UB) * ((-RHO*I)*LG) * dag(UC)


    single_su3 gl3_temp0,gl3_temp1,gl3_temp2;
    single_thmat thmat_temp0;

    single_su3_from_su3_soa(UC,idxC,&gl3_temp0);//UC
    rebuild3row(&gl3_temp0);
    gl3_dagger(&gl3_temp0);//dag(UC)

    single_thmat_from_thmat_soa(LG,idxG,&thmat_temp0);// LG
    thmat_to_su3(&gl3_temp1,&thmat_temp0);// to calculate products

    single_su3xsu3(&gl3_temp2,&gl3_temp1,&gl3_temp0);// LG * dag(UC)
    // gl3_temp0  gl3_temp1     gl3_temp2
    //  dag(UC)      LG              res


    single_su3_from_su3_soa(UB,idxB,&gl3_temp0);//UB
    rebuild3row(&gl3_temp0);
    gl3_dagger(&gl3_temp0);//dag(UB)

    single_su3xsu3(&gl3_temp1,&gl3_temp0,&gl3_temp2); // dag(UB) * LG * dag(UC)
    // gl3_temp0  gl3_temp1  gl3_temp2
    //  dag(UB)     res      LG*dag(UC)

    single_su3_from_su3_soa(UA,idxA,&gl3_temp0);//UA
    rebuild3row(&gl3_temp0);

    //UA * dag(UB)*LG*dag(UC)
    single_su3xsu3(&gl3_temp2,&gl3_temp0,&gl3_temp1);//UA*dag(UB)*LG*dag(UC)
    // gl3_temp0     gl3_temp1        gl3_temp2
    //  UA       dag(UB)*LG*dag(UC)     res


    single_su3_times_scalar(&gl3_temp2, -RHO*I);
    // ^^ (-RHO*I)*UA*dag(UB)*LG*dag(UC)


    single_gl3_addinto_su3_soa(RES,idxRES,&gl3_temp2);

}


#pragma acc routine seq
static inline void DEOTT_LEFT_iAB_times_GminusE_times_C_absent_stag_phases(  __restrict su3_soa * const UA,
        const int idxA,
        __restrict su3_soa * const UB,
        const int idxB,
        __restrict su3_soa * const UC,
        const int idxC,
        __restrict thmat_soa * const LE,
        const int idxE,
        __restrict thmat_soa * const LG,
        const int idxG,
        __restrict su3_soa * const RES,
        const int idxRES){
    // Cosa calcoliamo in questa routine:
    //  RES += dag(UA) * dag(UB) * ((RHO*I)*(LG - LE)) * UC


    single_su3 gl3_temp0,gl3_temp1,gl3_temp2;
    single_thmat thmat_temp0, thmat_temp1;

    single_su3_from_su3_soa(UC,idxC,&gl3_temp0);//UC
    rebuild3row(&gl3_temp0);

    single_thmat_from_thmat_soa(LG,idxG,&thmat_temp0);
    single_thmat_from_thmat_soa(LE,idxE,&thmat_temp1);
    single_thmatAeqAmB(&thmat_temp0,&thmat_temp1);// sLG = sLG-sLE

    thmat_to_su3(&gl3_temp1,&thmat_temp0);// to calculate products

    single_su3xsu3(&gl3_temp2,&gl3_temp1,&gl3_temp0); // (LG-LE)* UC
    // gl3_temp0  gl3_temp1  gl3_temp2
    //   UC        (LG-LE)    res   


    single_su3_from_su3_soa(UB,idxB,&gl3_temp0);//UB
    rebuild3row(&gl3_temp0);
    gl3_dagger(&gl3_temp0);//dag(UB)
    single_su3xsu3(&gl3_temp1,&gl3_temp0,&gl3_temp2);//dag(UB)*(LG-LE)*UC
    // gl3_temp0  gl3_temp1       gl3_temp2
    //  dag(UB)      res       dag(UC)*(LD-LE)

    single_su3_from_su3_soa(UA,idxA,&gl3_temp0);//UA
    rebuild3row(&gl3_temp0);
    gl3_dagger(&gl3_temp0);//dag(UA)

    //dag(UA) * dag(UB) * (LG-LE) * UC
    single_su3xsu3(&gl3_temp2,&gl3_temp0,&gl3_temp1);
    // gl3_temp0      gl3_temp1        gl3_temp2
    //  dag(UA)    dag(UB)*(LG-LE)*UC     res

    // multiplication by RHO*I
    single_su3_times_scalar(&gl3_temp2,RHO*I);

    //  RES += dag(UA) * dag(UB) * ((RHO*I)*(LG - LE)) * UC
    single_gl3_addinto_su3_soa(RES,idxRES,&gl3_temp2);
}

#pragma acc routine seq
static inline void DEOTT_LEFT_iABCD_absent_stag_phases(  __restrict su3_soa * const UA,
        const int idxA,
        __restrict su3_soa * const UB,
        const int idxB,
        __restrict su3_soa * const UC,
        const int idxC,
        __restrict thmat_soa * const LD,
        const int idxD,
        __restrict su3_soa * const RES,
        const int idxRES){
    // Cosa calcoliamo in questa routine:
    //  RES += dag(UA) * dag(UB) * UC * ((RHO*I)*LD


    single_su3 gl3_temp0,gl3_temp1,gl3_temp2;
    single_thmat thmat_temp0;

    single_thmat_from_thmat_soa(LD,idxD,&thmat_temp0);// LD
    thmat_to_su3(&gl3_temp0,&thmat_temp0);// to calculate products

    single_su3_from_su3_soa(UC,idxC,&gl3_temp1);//UC
    rebuild3row(&gl3_temp1);

    single_su3xsu3(&gl3_temp2,&gl3_temp1,&gl3_temp0);// UC * LD
    // gl3_temp0  gl3_temp1     gl3_temp2
    //  LD          UC              res

    single_su3_from_su3_soa(UB,idxB,&gl3_temp0);//UB
    rebuild3row(&gl3_temp0);
    gl3_dagger(&gl3_temp0);//dag(UB)

    single_su3xsu3(&gl3_temp1,&gl3_temp0,&gl3_temp2); // dag(UB) * UC * LD
    // gl3_temp0  gl3_temp1  gl3_temp2
    //   dag(UB)    res       UC * LD

    single_su3_from_su3_soa(UA,idxA,&gl3_temp0);//UA
    rebuild3row(&gl3_temp0);
    gl3_dagger(&gl3_temp0);//dag(UA)

    single_su3xsu3(&gl3_temp2,&gl3_temp0,&gl3_temp1);//dag(UA)*dag(UB)*UC*LD
    // gl3_temp0     gl3_temp1        gl3_temp2
    //  UA         dag(UB)*UC*LD          res


    single_su3_times_scalar(&gl3_temp2, RHO*I);
    // ^^ (RHO*I)*dag(UA)*dag(UB)*UC*LD

    single_gl3_addinto_su3_soa(RES,idxRES,&gl3_temp2);

}



#pragma acc routine seq
static inline void DEOTT_LEFT_miAFBC_absent_stag_phases(  __restrict su3_soa * const UA,
        const int idxA,
        __restrict su3_soa * const UB,
        const int idxB,
        __restrict su3_soa * const UC,
        const int idxC,
        __restrict thmat_soa * const LF,
        const int idxF,
        __restrict su3_soa * const RES,
        const int idxRES){
    // Cosa calcoliamo in questa routine:
    //  RES += dag(UA) * (-RHO*I)*LF * dag(UB) * UC

    single_su3 gl3_temp0,gl3_temp1,gl3_temp2;
    single_thmat thmat_temp0;

    single_su3_from_su3_soa(UC,idxC,&gl3_temp0);//UC
    rebuild3row(&gl3_temp0);

    single_su3_from_su3_soa(UB,idxB,&gl3_temp1);//UB
    rebuild3row(&gl3_temp1);
    gl3_dagger(&gl3_temp1);//dag(UB)


    single_su3xsu3(&gl3_temp2,&gl3_temp1,&gl3_temp0); // dag(UB) * UC
    // gl3_temp0  gl3_temp1  gl3_temp2
    //  UC         dag(UB)     res   

    single_thmat_from_thmat_soa(LF,idxF,&thmat_temp0);// LF
    thmat_to_su3(&gl3_temp0,&thmat_temp0);// to calculate products

    single_su3xsu3(&gl3_temp1,&gl3_temp0,&gl3_temp2);// LF * dag(UB) * UC
    // gl3_temp0   gl3_temp1     gl3_temp2
    //    LF          res        dag(UB)*UC


    single_su3_from_su3_soa(UA,idxA,&gl3_temp0);//UA
    rebuild3row(&gl3_temp0);
    gl3_dagger(&gl3_temp0);//dag(UA)

    single_su3xsu3(&gl3_temp2,&gl3_temp0,&gl3_temp1);// dag(UA)*LF*dag(UB)*UC
    // gl3_temp0   gl3_temp1     gl3_temp2
    //  dag(UA)  LF*dag(UB)*UC     res

    single_su3_times_scalar(&gl3_temp2, -RHO*I);
    // (-RHO*I) *dag(UA)*LF*dag(UB)*UC

    single_gl3_addinto_su3_soa(RES,idxRES,&gl3_temp2);
}


void DEOTT_compute_sigma(__restrict thmat_soa * const L,  // la Lambda --> ouput  (una cosa che serve per calcolare la forza fermionica successiva)
        __restrict su3_soa   * const U,  // la configurazione di gauge --> input
        __restrict su3_soa   * const S,  // entra Sigma primo (input: fermforce del passo precedente) ED esce Sigma --> sia input che ouput
        __restrict tamat_soa * const QA, // gli stessi Q che arrivano a Cayley hamilton --> input (sostanzialmente sono rho*ta(staples))
        __restrict su3_soa   * const TMP // variabile di parcheggio
        ){
  SETINUSE(TMP) ;
  SETINUSE(S);
    int x, y, z, t, mu, iter;

#pragma acc kernels present(L) present(U) present(nnp_openacc) present(nnm_openacc) present(S) present(QA) present(TMP)
#pragma acc loop independent gang
    for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector(4)
        for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector(4)
            for(y=0; y<ny; y++) {
#pragma acc loop independent vector(32)
                for(x=0; x < nx; x++) {
#pragma acc loop seq
                    for(mu=0; mu<4; mu++){

                        const int indice = snum_acc(x,y,z,t); // anche questa e' calcolata pure dopo ...
                        const int parita = (x+y+z+t) % 2;
                        const int dir_mu = 2*mu + parita; // definisco una variabile di tipo dirlink pure qui (ce ne e' anche una dentro il loop dopo)

                        // in questa routine faccio RES = PEZZO1
                        DEOTT_compute_sigma_local_PEZZO1(&L[dir_mu],&U[dir_mu],&S[dir_mu],&QA[dir_mu],&TMP[dir_mu],indice);

#pragma acc loop seq
                        for(iter=0; iter<3; iter++){
                            int nu;
                            if (mu==0) { nu = iter + 1; }
                            else if (mu==1) { nu = iter + (iter & 1) + (iter >> 1); }
                            else if (mu==2) { nu = iter + (iter >> 1); }
                            else if (mu==3) { nu = iter; }
                            else { //error
                            }

                            const int idxh = snum_acc(x,y,z,t);  // r
                            const int parity = (x+y+z+t) % 2;
#pragma acc cache (nnp_openacc[idxh:8])

                            const int dir_link = 2*mu + parity;
                            const int dir_mu_2R = 2*mu + (1-parity);
                            const int dir_mu_2L = 2*mu + (1-parity);
                            const int idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
#pragma acc cache (nnm_openacc[idx_pmu:8])
                            const int dir_nu_1R = 2*nu + (1-parity);
                            const int dir_nu_3R = 2*nu +  parity;
                            const int dir_nu_1L = 2*nu +  parity;
                            const int dir_nu_3L = 2*nu + (1-parity);

                            const int idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu

                            //  computation of the >>>Right<<< part of the staple with LAMBDAs
                            //  su3
                            //  A = U_nu(x+mu)
                            //  B = dagger U_mu(x+nu)
                            //  C = dagger U_nu(x)
                            //  thmat
                            //  D = LAMBDA_mu(x)
                            //  E = LAMBDA_nu(x)
                            //  F = LAMBDA_nu(x+mu)
                            //  G = LAMBDA_mu(x+nu)
                            //  iABCD and -iABCE
                            DEOTT_RIGHT_iABC_times_DminusE_absent_stag_phases(&U[dir_nu_1R],       idx_pmu, // A
                                    &U[dir_mu_2R],       idx_pnu, // B
                                    &U[dir_nu_3R],       idxh,    // C
                                    &L[dir_link],        idxh,    // D
                                    &L[dir_nu_3R],       idxh,    // E
                                    &S[dir_link],        idxh);
                            //  iFABC
                            DEOTT_RIGHT_iFABC_absent_stag_phases(&U[dir_nu_1R],       idx_pmu, // A
                                    &U[dir_mu_2R],       idx_pnu, // B
                                    &U[dir_nu_3R],       idxh,    // C
                                    &L[dir_nu_1R],       idx_pmu, // F
                                    &S[dir_link],        idxh);
                            // -iABGC
                            DEOTT_RIGHT_miABGC_absent_stag_phases(&U[dir_nu_1R],       idx_pmu, // A
                                    &U[dir_mu_2R],       idx_pnu, // B
                                    &U[dir_nu_3R],       idxh,    // C
                                    &L[dir_mu_2R],       idx_pnu, // G
                                    &S[dir_link],        idxh);

                            const int idx_mnu = nnm_openacc[idxh][nu][parity] ;         // r-nu
                            const int idx_pmu_mnu = nnm_openacc[idx_pmu][nu][(1-parity)];  // r+mu-nu

                            //  computation of the >>>Left<<< part of the staple with LAMBDAs
                            //  su3
                            //  A = dagger U_nu(x+mu-nu)
                            //  B = dagger U_mu(x-nu)
                            //  C = U_nu(x-nu)
                            //  thmat
                            //  D = LAMBDA_mu(x)
                            //  E = LAMBDA_mu(x-nu)
                            //  F = LAMBDA_nu(x+mu-nu)
                            //  G = LAMBDA_nu(x-nu)
                            //  iABGC and -iABEC
                            DEOTT_LEFT_iAB_times_GminusE_times_C_absent_stag_phases(&U[dir_nu_1L],       idx_pmu_mnu, // A
                                    &U[dir_mu_2L],       idx_mnu,     // B
                                    &U[dir_nu_3L],       idx_mnu,     // C
                                    &L[dir_mu_2L],       idx_mnu,     // E
                                    &L[dir_nu_3L],       idx_mnu,     // G
                                    &S[dir_link],        idxh);

                            //  iABCD

                            DEOTT_LEFT_iABCD_absent_stag_phases(&U[dir_nu_1L],       idx_pmu_mnu, // A
                                    &U[dir_mu_2L],       idx_mnu,     // B
                                    &U[dir_nu_3L],       idx_mnu,     // C
                                    &L[dir_link],        idxh,        // D
                                    &S[dir_link],        idxh);
                            // -iAFBC
                            DEOTT_LEFT_miAFBC_absent_stag_phases(&U[dir_nu_1L],       idx_pmu_mnu, // A
                                    &U[dir_mu_2L],       idx_mnu,     // B
                                    &U[dir_nu_3L],       idx_mnu,     // C
                                    &L[dir_nu_1L],       idx_pmu_mnu, // F
                                    &S[dir_link],        idxh);
                        }  // iter


                    } // mu

                }  // x
            }  // y
        }  // z
    }  // t

  SETFREE(TMP);
}// closes routine



#endif

