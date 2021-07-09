//
//  Mod_functions.c
//  
//
//  Created by Luca Parente on 02/06/21.
//

#include "Mod_functions.h"
#include "./struct_c_def.h"
#include "./alloc_settings.h"
#include "./alloc_vars.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "./geometry.h"
#include "../Mpi/multidev.h"
#include "./plaquettes.h"
#include "./su3_utilities.h"
#include "../Include/common_defines.h"
#include "./single_types.h"
#include "./rettangoli.h"





/////funzione aggiunta//////////////////////
//per ora la metto qui poi la sposto in struct_def.h

//this function does actually nothing useful, I wrote it down just for studying how bc have been implementated.
void counter_size_function(int d0,int d1, int d2, int d3){
    int id=0;
    for(d3=0; d3<nd3; d3++) {
        for(d2=0; d2<nd2; d2++) {
            for(d1=0; d1<nd1; d1++) {
                for(d0=0; d0 < nd0; d0++) { //i cicli sui siti del reticolo.
                    
                    id=snum_acc(d0,d1,d2,d3);
                    printf("(%d,%d,%d,%d):      id=%d\n",d0,d1,d2,d3,id);
                }
            }
        }
    }
    
    return;
}

/*
void init_k_values(su3_soa * conf,int c_r,int * pos_def){


    for(mu=0;mu<8;mu++){
        for(i=0; i<sizeh; i++ ){
            
            
            conf[mu].K.d[i]=c_r;
            }
    
    conf[mu].K.d[snum_acc(2*i,j,k,t)]
    return;
}
*/


int init_k(su3_soa * conf,double c_r,int def_axis,int * def_vet){
    int mu;
    int i,j,z,t;
    int res=0;
    int defect_volume;
    int parity;
    
    
    int counter=0;
   //il case considera le 4 possibili direzioni in cui viene fatto il defect.
    
    switch (def_axis) {
        case 0:
            printf("defect on x's boundary\n");
        for(mu=0;mu<4;mu++){
           for(t=0;t<nd3;t++) {
                for (z=0; z<nd2; z++){
                     for(j=0;j<nd1;j++){
                       for(i=0;i<nd0;i++){
                            parity = (i+j+z+t) % 2;
                           
                            if(j>=0 && j<def_vet[0] && z>=0 && z<def_vet[1] && t>=0 && t<def_vet[2] && i==((nd0)-1) && mu==0 )
                            {
                                
                                if (parity==0){conf[2*mu].K.d[snum_acc(i,j,z,t)]=c_r;
                                  /*  printf("(%d,%d,%d,%d):  k_mu[%d]=%f (%d)parity\n",i,j,z,t,snum_acc(i,j,z,t),conf[2*mu].K.d[snum_acc(i,j,z,t)],parity);*/
                                } //inizializza il vettore
                                if(parity!=0){conf[2*mu+1].K.d[snum_acc(i,j,z,t)]=c_r;
                                    /*printf("(%d,%d,%d,%d):  k_mu[%d]=%f (%d)parity\n",i,j,z,t,snum_acc(i,j,z,t),conf[2*mu+1].K.d[snum_acc(i,j,z,t)],parity);*/
                                } //inizializza il vettore
                                    
                                
                                counter=counter+1;
                                
                                ;
                     /* printf("(%d,%d,%d,%d):     k_mu[%d]=%f\n",2*i,j,z,t,snum_acc(2*i,j,z,t),conf[mu].K.d[snum_acc(2*i,j,z,t)]);*/
                                 }
                                
                           
                        
                            else{
                                if (parity==0){conf[2*mu].K.d[snum_acc(i,j,z,t)]=1;}
                                else{conf[2*mu+1].K.d[snum_acc(i,j,z,t)]=1; }
                                
                                
                                
                            } //else inizializza a 1.
                        
                     /*printf("(%d,%d,%d,%d):     k_mu[%d]=%f\n",2*i,j,z,t,snum_acc(2*i,j,z,t),conf[mu].K.d[snum_acc(2*i,j,z,t)]);*/
                       }
                    }
                }
            }
            
        }
        
            
            
            
            break;
            
        case 1:
            printf("defect on y's boundary\n");
            for(mu=0;mu<4;mu++){
                for(t=0;t<nd3;t++) {
                    for (z=0; z<nd2; z++){
                        for(j=0;j<nd1;j++){
                            for(i=0;i<nd0;i++){
                            if(i>=0 && i<def_vet[0] && z>=0 && z<def_vet[1] && t>=0 && t<def_vet[2] && j==nd1-1 && mu==1 )
                            {
                                 parity = (i+j+z+t) % 2;
                                
                                if (parity==0){conf[2*mu].K.d[snum_acc(i,j,z,t)]=c_r;} //inizializza il vettore}
                                if(parity!=0){conf[2*mu+1].K.d[snum_acc(i,j,z,t)]=c_r;} //inizializza il vettore}
                                
                                
                                counter=counter+1;
                                
                            }
                                
                                
                                
                            else{
                                if (parity==0){conf[2*mu].K.d[snum_acc(i,j,z,t)]=1;}
                                else{conf[2*mu+1].K.d[snum_acc(i,j,z,t)]=1;}
                                
                                
                                
                            } //else inizializza a 1.
                                
                                /*printf("(%d,%d,%d,%d):     k_mu[%d]=%f\n",2*i,j,z,t,snum_acc(2*i,j,z,t),conf[mu].K.d[snum_acc(2*i,j,z,t)]);*/
                                
                    }
              }
            }
          }
        }
            
            
            
            
            
            break;
            
        case 2:
            printf("defect on z's boundary\n");
            for(mu=0;mu<4;mu++){
                for(t=0;t<nd3;t++) {
                    for (z=0; z<nd2; z++){
                        for(j=0;j<nd1;j++){
                            for(i=0;i<nd0;i++){
                            if(i>=0 && i<def_vet[0] && j>=0 && j<def_vet[1] && t>=0 && t<def_vet[2] && z==nd2-1 && mu==2 )
                            {
                                 parity = (i+j+z+t) % 2;
                                if (parity==0){conf[2*mu].K.d[snum_acc(i,j,z,t)]=c_r;} //inizializza il vettore}
                                if(parity!=0){conf[2*mu+1].K.d[snum_acc(i,j,z,t)]=c_r;} //inizializza il vettore}
                                
                                
                                counter=counter+1;
                                /* printf("(%d,%d,%d,%d):     k_mu[%d]=%f\n",2*i,j,z,t,snum_acc(2*i,j,z,t),conf[mu].K.d[snum_acc(2*i,j,z,t)]);*/
                            }
                                
                                
                                
                            else{
                                if (parity==0){conf[2*mu].K.d[snum_acc(i,j,z,t)]=1;}
                                
                                else{conf[2*mu+1].K.d[snum_acc(i,j,z,t)]=1;}
                                
                                
                                
                            } //else inizializza a 1.
                                
                                /*printf("(%d,%d,%d,%d):     k_mu[%d]=%f\n",2*i,j,z,t,snum_acc(2*i,j,z,t),conf[mu].K.d[snum_acc(2*i,j,z,t)]);*/
                    }
                }
            }
            
           }
      }
            
            
            
            break;
            
        case 3:
            printf("defect on t's boundary\n");
            for(mu=0;mu<4;mu++){
                for(t=0;t<nd3;t++) {
                    for (z=0; z<nd2; z++){
                        for(j=0;j<nd1;j++){
                            for(i=0;i<nd0;i++){
                            if(i>=0 && i<def_vet[0] && j>=0 && j<def_vet[1] && z>=0 && z<def_vet[2] && t==nd3-1 && mu==3 )
                            {
                                 parity = (i+j+z+t) % 2;
                                
                                if (parity==0){conf[2*mu].K.d[snum_acc(i,j,z,t)]=c_r;} //inizializza il vettore}
                                if(parity!=0){conf[2*mu+1].K.d[snum_acc(i,j,z,t)]=c_r;} //inizializza il vettore}
                                
                                
                                counter=counter+1;
                                /* printf("(%d,%d,%d,%d):     k_mu[%d]=%f\n",2*i,j,z,t,snum_acc(2*i,j,z,t),conf[mu].K.d[snum_acc(2*i,j,z,t)]);*/
                            }
                                
                                
                                
                            else{
                                if (parity==0){conf[2*mu].K.d[snum_acc(i,j,z,t)]=1;}
                                else{conf[2*mu+1].K.d[snum_acc(i,j,z,t)]=1;}
                                
                                
                                
                            } //else inizializza a 1.
                                
                                /*printf("(%d,%d,%d,%d):     k_mu[%d]=%f\n",2*i,j,z,t,snum_acc(2*i,j,z,t),conf[mu].K.d[snum_acc(2*i,j,z,t)]);*/
                    }
                }
            }
            
         }
        }
            break;
            
            
            
        default:
            printf("ERROR WRONG AXIS CHOICE!\n");
            res=1;
            break;
    }
    
    defect_volume=(def_vet[0])*(def_vet[1])*(def_vet[2]);
    
   printf("counter %d\n",counter);
    
    if(counter!=defect_volume){printf("wrong defect initialization!\n"); res=1;}
    
    return res;
  
}






//function that tested the initialization of k_mu  by setting every value to one. //OUTADATED
int init_k_test(su3_soa *conf_acc,int c_r){
    int kk2=0;
    int mu1=0;
    for(mu1=0;mu1<8;mu1++){
        for(kk2=0;kk2<sizeh;kk2++){
            if(conf_acc[mu1].K.d[kk2]!=(1 && c_r)){
                return 1;
            }
        printf("%d:ku[%d]:%f\n",mu1,kk2,conf_acc[mu1].K.d[kk2]);
        }
    }
    return 0;
}
/* work in progress
int defect_shifter_function(su3_soa * conf,int def_axis,int * def_vet){
        int mu;
    double aux;
    
    switch (def_axis){
            
         case 0: for(mu=0;mu<8;mu++){
             for(t=def_vet[4];t<def_vet[5];t++) {
                 for (z=def_vet[2]; z<def_vet[3]; z++){
                     for(j=def_vet[0];j<def_vet[1];j++){
             
             
                         conf[mu].K.d[snum_acc((nd0/2)-1),j,z,t)];
                         
            
                    }
                 }
             }
         }
            break;
            
            
        case 1:
            
    }
    
    
    
    
    return 0;
}


*/




//Outdated and useless function. Once used before a better mod of the setting_parser_file.c
/*
int n_replicas_reader(const char* input_filename){
    int value_nr;
    int trovato=0;
    int i=0;
    char riga[20];
    char riga2[20]="Replicas number";
    FILE *input = fopen(input_filename,"r"); //questo ovviamente apre il file
    printf("LETTURA DEL NUMERO DI REPLICHE\n");
    while(trovato==0){
        fgets(riga,20,input);
      //  printf("ecc %d\n",i);
 
        //i=i+1;
        if(strncmp(riga,riga2,15)==0){
          fscanf(input,"%d",&value_nr );
            trovato=1;
        }
        
        
        
        
        
    }
    
    
    
    fclose(input);
    
    return value_nr;
}
*/


void printing_k_mu(su3_soa * conf){
    int mu,t,z,j,i;
    
    for(mu=0;mu<8;mu++){
        for(t=0;t<nd3;t++) {
            for (z=0; z<nd2; z++){
                for(j=0;j<nd1;j++){
                    for(i=0;i<(nd0/2);i++){
    
    printf("(%d,%d,%d,%d):     k_mu[%d]=%f\n",2*i,j,z,t,snum_acc(2*i,j,z,t),conf[mu].K.d[snum_acc(2*i,j,z,t)]);
                    }
                }
            }
            
        }
    }
    return;
}

int replicas_swap_1(su3_soa * conf1,su3_soa * conf2,int def_axis,int * def_vet ){
    double aux;
    int i,j,k,t,mu,parity;
    int res=0;
    //test variable
    int counter=0;

    switch (def_axis){
            
        case 0:
            i=nd0-1;
                mu=0;
                for(t=0;t<def_vet[2];t++){
                    for(k=0;k<def_vet[1];k++){
                        for(j=0;j<def_vet[0];j++){
                            
                            parity = (i+j+k+t) % 2;
                            
                            //test
                            
                            
                            //K_mu_values swap
                            if (parity==0){aux=conf1[mu].K.d[snum_acc(i,j,k,t)];
                           /* printf("beforrre (%d) %f %f",counter,conf1[mu].K.d[snum_acc(i,j,k,t)],conf2[mu].K.d[snum_acc(i,j,k,t)] );*/
                                conf1[mu].K.d[snum_acc(i,j,k,t)]=conf2[mu].K.d[snum_acc(i,j,k,t)];
                                conf2[mu].K.d[snum_acc(i,j,k,t)]=aux;
                         /*   printf("aftermath (%d) %f %f\n",counter,conf1[mu].K.d[snum_acc(i,j,k,t)],conf2[mu].K.d[snum_acc(i,j,k,t)] );*/
                            }
                            if(parity!=0){aux=conf1[mu+1].K.d[snum_acc(i,j,k,t)];
                         /*   printf("beforrre (%d) %f %f",counter,conf1[mu+1].K.d[snum_acc(i,j,k,t)],conf2[mu+1].K.d[snum_acc(i,j,k,t)] );*/
                            conf1[mu+1].K.d[snum_acc(i,j,k,t)]=conf2[mu+1].K.d[snum_acc(i,j,k,t)];
                            conf2[mu+1].K.d[snum_acc(i,j,k,t)]=aux;
                        /*    printf("aftermath (%d) %f %f\n",counter,conf1[mu+1].K.d[snum_acc(i,j,k,t)],conf2[mu+1].K.d[snum_acc(i,j,k,t)] );*/
                            }
                            //test
                            
                           /* counter=counter+1;*/
                        }
                    }
                }
            break;
            
        case 1:
            j=nd0-1;
            mu=1;
            for(t=0;t<def_vet[2];t++){
                for(k=0;k<def_vet[1];k++){
                    for(i=0;i<def_vet[0];i++){
                         parity = (i+j+k+t) % 2;
                        
                        //K_mu_values swap
                        if (parity==0){aux=conf1[mu].K.d[snum_acc(i,j,k,t)];
                            conf1[mu].K.d[snum_acc(i,j,k,t)]=conf2[mu].K.d[snum_acc(i,j,k,t)];
                            conf2[mu].K.d[snum_acc(i,j,k,t)]=aux;
                        }
                        if(parity!=0){aux=conf1[mu+1].K.d[snum_acc(i,j,k,t)];
                            conf1[mu+1].K.d[snum_acc(i,j,k,t)]=conf2[mu+1].K.d[snum_acc(i,j,k,t)];
                            conf2[mu+1].K.d[snum_acc(i,j,k,t)]=aux;
                        }
                        
                        
                        
            
                   }
                }
            }
            break;
            
        case 2:
            k=nd0-1;
            mu=2;
                    for(t=0;t<def_vet[2];t++){
                        for(k=0;k<def_vet[1];k++){
                            for(j=0;j<def_vet[0];j++){
                                
                                parity = (i+j+k+t) % 2;
                                
                                //K_mu_values swap
                                if (parity==0){aux=conf1[mu].K.d[snum_acc(i,j,k,t)];
                                    printf("beforrre (%d) %f %f",counter,conf1[mu].K.d[snum_acc(i,j,k,t)],conf2[mu].K.d[snum_acc(i,j,k,t)] );
                                    conf1[mu].K.d[snum_acc(i,j,k,t)]=conf2[mu].K.d[snum_acc(i,j,k,t)];
                                    conf2[mu].K.d[snum_acc(i,j,k,t)]=aux;
                                    printf("aftermath (%d) %f %f\n",counter,conf1[mu].K.d[snum_acc(i,j,k,t)],conf2[mu].K.d[snum_acc(i,j,k,t)] );
                                }
                                if(parity!=0){aux=conf1[mu+1].K.d[snum_acc(i,j,k,t)];
                                    printf("beforrre (%d) %f %f",counter,conf1[mu+1].K.d[snum_acc(i,j,k,t)],conf2[mu+1].K.d[snum_acc(i,j,k,t)] );
                                    conf1[mu+1].K.d[snum_acc(i,j,k,t)]=conf2[mu+1].K.d[snum_acc(i,j,k,t)];
                                    conf2[mu+1].K.d[snum_acc(i,j,k,t)]=aux;
                                     printf("aftermath (%d) %f %f\n",counter,conf1[mu+1].K.d[snum_acc(i,j,k,t)],conf2[mu+1].K.d[snum_acc(i,j,k,t)] );
                                }
                                
                                
                            }
                        }
                    }
                                
            break;
            
        case 3:
            t=nd0-1;
            mu=3;
            for(t=0;t<def_vet[2];t++){
                for(k=0;k<def_vet[1];k++){
                    for(j=0;j<def_vet[0];j++){
                        
                        parity = (i+j+k+t) % 2;
                        
                        //K_mu_values swap
                        if (parity==0){aux=conf1[mu].K.d[snum_acc(i,j,k,t)];
                            conf1[mu].K.d[snum_acc(i,j,k,t)]=conf2[mu].K.d[snum_acc(i,j,k,t)];
                            conf2[mu].K.d[snum_acc(i,j,k,t)]=aux;
                        }
                        if(parity!=0){aux=conf1[mu+1].K.d[snum_acc(i,j,k,t)];
                            conf1[mu+1].K.d[snum_acc(i,j,k,t)]=conf2[mu+1].K.d[snum_acc(i,j,k,t)];
                            conf2[mu+1].K.d[snum_acc(i,j,k,t)]=aux;
                        }
                        
                    }
                }
            }
                        
            break;
            

        default:
            printf("ERROR WRONG AXIS CHOICE! (SWAP)\n");
            res=1;
            
    }
    
    
    //test
   /* printf("counter:%d\n%",counter);*/
    
 
    return res;
}
//replicas_swap function: 2 confs defects are exchanged.
int replicas_swap(su3_soa * conf1,su3_soa * conf2,int def_axis,int * def_vet ){
    vec3_soa  aux;
    int aux_label;
        int res=0;
    int mu=0;
    
    for(mu=0;mu<8;mu++){
        printf("beforrre (%d) %.18lf %.18lf",mu,creal(conf1[mu].r1.c1[snum_acc(31,6,6,6)]),creal(conf2[mu].r1.c1[snum_acc(31,6,6,6)]));
        
        //label swap.
        aux_label=conf1[mu].label;
        conf1[mu].label=conf2[mu].label;
        conf2[mu].label=aux_label;
        
        
        //swap r0
        aux=conf1[mu].r0;
        conf1[mu].r0=conf2[mu].r0;
        conf2[mu].r0=aux;
        
        
        
        //swap r1
        aux=conf1[mu].r1;
        conf1[mu].r1=conf2[mu].r1;
        conf2[mu].r1=aux;
        
        //swap r2
        aux=conf1[mu].r2;
        conf1[mu].r2=conf2[mu].r2;
        conf2[mu].r2=aux;
        
    printf("aftermath (%d) %.18lf %.18lf\n",mu,creal(conf1[mu].r1.c1[snum_acc(31,6,6,6)]),creal(conf2[mu].r1.c1[snum_acc(31,6,6,6)]));

    }
    
    return res;
}


//function which print the confs'labels.
int label_print(su3_soa ** conf_hasen, int replicas_number,FILE *file,int step_number){
    int res=0;
    int i;
    
    fprintf(file,"#step: %d\n",step_number);
    
    for(i=0;i<replicas_number;i++){
        fprintf(file,"%d   %d\n",i,conf_hasen[i][0].label);
        
    }
  
    
    
    return res;
}


//Function which chooses the plane and iterate only on useful them.


double  calc_plaquette_soloopenacc_SWAP(
                                   __restrict  su3_soa * const tconf_acc,
                                   __restrict su3_soa * const local_plaqs,
                                   dcomplex_soa * const tr_local_plaqs,int def_axis, int improved )
{
    
    
    double result=0.0;
    double total_result=0.0;
    int mu;
    
    int i_counter=0;
    // calcolo il valore della plaquette sommata su tutti i siti a fissato piano mu-nu (6 possibili piani)//(the couple has to be chosen excluding same direction ones.(4 2)binomial coefficient.
    
    switch (def_axis) {
        case 0:
            mu=0;
            for(int nu=mu+1;nu<4;nu++){
                // sommo i 6 risultati in tempo
                 if(improved==0){
                result  += calc_loc_plaquettes_nnptrick_SWAP(tconf_acc,local_plaqs,tr_local_plaqs,mu,nu); //here ol the plaquettes of a specific plane's choice are computed.
                 }
                
                
                if(improved==1){
                    result  += calc_loc_plaquettes_rectangles_SWAP(tconf_acc,local_plaqs,tr_local_plaqs,mu,nu);
                
                    
                }
                
                
            }
        
            
            break;
            
        case 1:
            mu=1;
            for(int nu=0;nu<4;nu++){
                // sommo i 6 risultati in tempo
                if(nu!=mu){
                     if(improved==0){
                result  += calc_loc_plaquettes_nnptrick_SWAP(tconf_acc,local_plaqs,tr_local_plaqs,mu,nu); //here ol the plaquettes of a specific plane's choice are computed.
                     }
                    
                    if(improved==1){
                        result  += calc_loc_plaquettes_rectangles_SWAP(tconf_acc,local_plaqs,tr_local_plaqs,mu,nu);
                        
                        
                    }
                
                }
                
                
            }
            
            
            break;
      
        case 2:
            mu=2;
            for(int nu=0;nu<4;nu++){
                // sommo i 6 risultati in tempo
                if(nu!=mu){
                    if(improved==0){
                    result  += calc_loc_plaquettes_nnptrick_SWAP(tconf_acc,local_plaqs,tr_local_plaqs,mu,nu); //here all the plaquettes of a specific plane's choice are computed.
                    }
                    if(improved==1){
                        result  += calc_loc_plaquettes_rectangles_SWAP(tconf_acc,local_plaqs,tr_local_plaqs,mu,nu);
                        
                        
                    }
                }
                
                
            }
            
            
            break;
            
        case 3:
            mu=3;
            for(int nu=0;nu<3;nu++){
                // sommo i 6 risultati in tempo
                 if(improved==0){
                    result  += calc_loc_plaquettes_nnptrick_SWAP(tconf_acc,local_plaqs,tr_local_plaqs,mu,nu); //here all the plaquettes of a specific plane's choice are computed.
                 }
                if(improved==1){
                    result  += calc_loc_plaquettes_rectangles_SWAP(tconf_acc,local_plaqs,tr_local_plaqs,mu,nu);
                    
                }
                
            }
            
            
            break;
            
            
        default:
            printf("DELTA_S_SWAP ERROR!\n");
            break;
    }
 
    
    
    
#ifdef MULTIDEVICE
    MPI_Allreduce((void*)&result,(void*)&total_result,
                  1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
    total_result = result;
#endif
    return total_result;
    
}


double calc_loc_plaquettes_nnptrick_SWAP{
    __restrict const su3_soa * const u,//for an unknown reason the vet conf is called u. this is a vector odf su3_soa.
    __restrict su3_soa * const loc_plaq, //la placchetta locale.
    dcomplex_soa * const tr_local_plaqs, //complex number that states the value of the trace. Of course is a vector of the struct dcomplex_soa.
    const int mu, const int nu, int def_axis, int *def_vet)
    {
        double K_mu_nu; //MOD.
        
        int d0, d1, d2, d3;
        
        
#pragma acc kernels present(u) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang(STAPGANG3)
        switch(def_axis){
            case 0:
                d0=nd0-1;
                for(d3=D3_HALO; d3<def_vet[2]+1-D3_HALO; d3++) {//what?
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
                    for(d2=0; d2<def_vet[1]+1; d2++) {
                        for(d1=0; d1<def_vet[0]+1; d1++) {
                            
                            
                        
                        int idxh,idxpmu,idxpnu; //idxh is the half-lattice position, idxpmu and idxpnu the nearest neighbours.
                        int parity; //parity
                        int dir_muA,dir_nuB; //mu and nu directions.
                        int dir_muC,dir_nuD;
                        
                        idxh = snum_acc(d0,d1,d2,d3);// the site on the  half-lattice.
                        parity = (d0+d1+d2+d3) % 2; //obviously the parity_term
                        idxh=nnm_openacc[idxh][nu][parity]; // the previous one. //MOD
                        
                        dir_muA = 2*mu +  parity;
                        dir_muC = 2*mu + !parity;
                        idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
                        
                        dir_nuB = 2*nu + !parity;
                        dir_nuD = 2*nu +  parity;
                        idxpnu = nnp_openacc[idxh][nu][parity];// r+nu //the table that states which is the nearest neighbour.
                        //       r+nu (C)  r+mu+nu
                        //          +<---+
                        // nu       |    ^
                        // ^    (D) V    | (B)
                        // |        +--->+
                        // |       r  (A)  r+mu
                        // +---> mu
                        
                        //(&u[dir_muA] & &u[dir_nuB] States which part of the the conf will be used. It is important to pass them as pointer, cause loc_plaq has to be modified.
                        
                        mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                        mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                        mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                        
                        d_complex ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                        tr_local_plaqs[parity].c[idxh] = creal(ciao)+cimag(ciao)*I;
                        
                        /* printf("%f +i%f ||",creal(tr_local_plaqs[parity].c[idxh]),cimag(tr_local_plaqs[parity].c[idxh])*I);*/
                        //MOD****************************************//
                        
                        //K_mu_nu computation;
                        K_mu_nu=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_muC].K.d[idxpnu])*(u[dir_nuD].K.d[idxh]);
                        
                        
                        tr_local_plaqs[parity].c[idxh]=K_mu_nu*tr_local_plaqs[parity].c[idxh];
                        //*****************************************//
                        
                       
                        
                }  // d1
            }  // d2
        }  // d3
                break;
               
            case 1:
                d1=nd1-1;
                for(d3=D3_HALO; d3<def_vet[2]+1-D3_HALO; d3++) {//what?
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
                    for(d2=0; d2<def_vet[1]+1; d2++) {
                        for(d0=0; d0<def_vet[0]+1; d0++) {
                            
                            
                            
                            int idxh,idxpmu,idxpnu; //idxh is the half-lattice position, idxpmu and idxpnu the nearest neighbours.
                            int parity; //parity
                            int dir_muA,dir_nuB; //mu and nu directions.
                            int dir_muC,dir_nuD;
                            
                            idxh = snum_acc(d0,d1,d2,d3);// the site on the  half-lattice.
                            parity = (d0+d1+d2+d3) % 2; //obviously the parity_term
                            idxh=nnm_openacc[idxh][nu][parity]; // the previous one. //MOD
                            
                            dir_muA = 2*mu +  parity;
                            dir_muC = 2*mu + !parity;
                            idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
                            
                            dir_nuB = 2*nu + !parity;
                            dir_nuD = 2*nu +  parity;
                            idxpnu = nnp_openacc[idxh][nu][parity];// r+nu //the table that states which is the nearest neighbour.
                            //       r+nu (C)  r+mu+nu
                            //          +<---+
                            // nu       |    ^
                            // ^    (D) V    | (B)
                            // |        +--->+
                            // |       r  (A)  r+mu
                            // +---> mu
                            
                            //(&u[dir_muA] & &u[dir_nuB] States which part of the the conf will be used. It is important to pass them as pointer, cause loc_plaq has to be modified.
                            
                            mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                            
                            d_complex ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            tr_local_plaqs[parity].c[idxh] = creal(ciao)+cimag(ciao)*I;
                            
                            /* printf("%f +i%f ||",creal(tr_local_plaqs[parity].c[idxh]),cimag(tr_local_plaqs[parity].c[idxh])*I);*/
                            //MOD****************************************//
                            
                            //K_mu_nu computation;
                            K_mu_nu=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_muC].K.d[idxpnu])*(u[dir_nuD].K.d[idxh]);
                            
                            
                            tr_local_plaqs[parity].c[idxh]=K_mu_nu*tr_local_plaqs[parity].c[idxh];
                            //*****************************************//
                            
                            
                            
                        }  // d1
                    }  // d2
                }  // d3
                break;
                
            case 2:
                d2=nd2-1;
                for(d3=D3_HALO; d3<def_vet[2]+1-D3_HALO; d3++) {//what?
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
                    for(d1=0; d1<def_vet[1]+1; d1++) {
                        for(d0=0; d0<def_vet[0]+1; d0++) {
                            
                            
                            
                            int idxh,idxpmu,idxpnu; //idxh is the half-lattice position, idxpmu and idxpnu the nearest neighbours.
                            int parity; //parity
                            int dir_muA,dir_nuB; //mu and nu directions.
                            int dir_muC,dir_nuD;
                            
                            idxh = snum_acc(d0,d1,d2,d3);// the site on the  half-lattice.
                            parity = (d0+d1+d2+d3) % 2; //obviously the parity_term
                            idxh=nnm_openacc[idxh][nu][parity]; // the previous one. //MOD
                            
                            dir_muA = 2*mu +  parity;
                            dir_muC = 2*mu + !parity;
                            idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
                            
                            dir_nuB = 2*nu + !parity;
                            dir_nuD = 2*nu +  parity;
                            idxpnu = nnp_openacc[idxh][nu][parity];// r+nu //the table that states which is the nearest neighbour.
                            //       r+nu (C)  r+mu+nu
                            //          +<---+
                            // nu       |    ^
                            // ^    (D) V    | (B)
                            // |        +--->+
                            // |       r  (A)  r+mu
                            // +---> mu
                            
                            //(&u[dir_muA] & &u[dir_nuB] States which part of the the conf will be used. It is important to pass them as pointer, cause loc_plaq has to be modified.
                            
                            mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                            
                            d_complex ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            tr_local_plaqs[parity].c[idxh] = creal(ciao)+cimag(ciao)*I;
                            
                            /* printf("%f +i%f ||",creal(tr_local_plaqs[parity].c[idxh]),cimag(tr_local_plaqs[parity].c[idxh])*I);*/
                            //MOD****************************************//
                            
                            //K_mu_nu computation;
                            K_mu_nu=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_muC].K.d[idxpnu])*(u[dir_nuD].K.d[idxh]);
                            
                            
                            tr_local_plaqs[parity].c[idxh]=K_mu_nu*tr_local_plaqs[parity].c[idxh];
                            //*****************************************//
                            
                            
                            
                        }  // d1
                    }  // d2
                }  // d3
                break;
                
            case 3:
                d3=nd3-1-D3_HALO;
                for(d2=0; d2<def_vet[2]+1; d2++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
                    for(d1=0; d1<def_vet[1]+1; d1++) {
                        for(d0=0; d0<def_vet[0]+1; d0++) {
                            
                            
                            
                            int idxh,idxpmu,idxpnu; //idxh is the half-lattice position, idxpmu and idxpnu the nearest neighbours.
                            int parity; //parity
                            int dir_muA,dir_nuB; //mu and nu directions.
                            int dir_muC,dir_nuD;
                            
                            idxh = snum_acc(d0,d1,d2,d3);// the site on the  half-lattice.
                            parity = (d0+d1+d2+d3) % 2; //obviously the parity_term
                            idxh=nnm_openacc[idxh][nu][parity]; // the previous one. //MOD
                            
                            dir_muA = 2*mu +  parity;
                            dir_muC = 2*mu + !parity;
                            idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
                            
                            dir_nuB = 2*nu + !parity;
                            dir_nuD = 2*nu +  parity;
                            idxpnu = nnp_openacc[idxh][nu][parity];// r+nu //the table that states which is the nearest neighbour.
                            //       r+nu (C)  r+mu+nu
                            //          +<---+
                            // nu       |    ^
                            // ^    (D) V    | (B)
                            // |        +--->+
                            // |       r  (A)  r+mu
                            // +---> mu
                            
                            //(&u[dir_muA] & &u[dir_nuB] States which part of the the conf will be used. It is important to pass them as pointer, cause loc_plaq has to be modified.
                            
                            mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                            
                            d_complex ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            tr_local_plaqs[parity].c[idxh] = creal(ciao)+cimag(ciao)*I;
                            
                            /* printf("%f +i%f ||",creal(tr_local_plaqs[parity].c[idxh]),cimag(tr_local_plaqs[parity].c[idxh])*I);*/
                            //MOD****************************************//
                            
                            //K_mu_nu computation;
                            K_mu_nu=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_muC].K.d[idxpnu])*(u[dir_nuD].K.d[idxh]);
                            
                            
                            tr_local_plaqs[parity].c[idxh]=K_mu_nu*tr_local_plaqs[parity].c[idxh];
                            //*****************************************//
                            
                            
                            
                        }  // d1
                    }  // d2
                }  // d3
                break;


                
                
                
            default:
                printf("DELTA_S_SWAP 1x1 plaquette ERROR!\n");
                break;

                
                
    }
        double res_R_p = 0.0;
        double res_I_p = 0.0;
        double resR = 0.0;
        int t;
        
        
        
#pragma acc kernels present(tr_local_plaqs)
#pragma acc loop reduction(+:res_R_p) reduction(+:res_I_p)
        for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
            res_R_p += creal(tr_local_plaqs[0].c[t]); //even sites plaquettes
            
            res_R_p += creal(tr_local_plaqs[1].c[t]); //odd sites plaquettes
        }
        
        
        
        return res_R_p;
    }// closes routine
    


double calc_loc_plaquettes_rectangles_SWAP(tconf_acc,local_plaqs,tr_local_plaqs,mu,nu){
    
    return 0;
}
    
    


