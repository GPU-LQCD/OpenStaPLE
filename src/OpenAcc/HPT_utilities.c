#ifndef HPT_UTILITIES_C
#define HPT_UTILITIES_C

#include "../Include/common_defines.h"
#ifdef PAR_TEMP

// this library containes all Hasenbusch Parallel Tempering (HPT) routines

#include "HPT_utilities.h"
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
#include "./single_types.h"
#include "./rectangles.h"
#include "./plaquettes.h"
#include "../Include/debug.h"
#include "../Rand/random.h"
#include "../Include/rep_info.h"

#include <time.h>

void init_k(su3_soa * conf, double c_r, int def_axis, int * def_vec, defect_info * def, int defect_info_config){
  
  int mu, i, parity, condition, condition_2=1, condition_3=1;
  int def_vec_4d[4] = {0};
  int perp_dir[4][3] = { {1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2} };

  int map[4] = {geom_par.xmap, geom_par.ymap, geom_par.zmap, geom_par.tmap};
  int def_axis_mapped=map[def_axis];
  if(defect_info_config==0){ def->def_axis_mapped=def_axis_mapped; }

  // initialization
  int nd[4]={nd0-1,nd1-1,nd2-1,nd3-1};
  if(defect_info_config==0){
    for(int nu=0; nu<4;nu++){
      for (i=0; i<3; i++){
				def->defect_swap_max[nu][i]=0;
				def->defect_swap_min[nu][i]=nd[i];
#ifdef GAUGE_ACT_TLSM
				for(int j=0;j<2;j++){
					def->defect_swap_max_TLSM[j][nu][i]=0;
					def->defect_swap_min_TLSM[j][nu][i]=nd[i];
				}
#endif
      }
      def->defect_swap_max[nu][3]=D3_HALO;
      def->defect_swap_min[nu][3]=nd[3]-D3_HALO;
#ifdef GAUGE_ACT_TLSM
      for(int j=0;j<2;j++){
				def->defect_swap_max_TLSM[j][nu][3]=D3_HALO;
				def->defect_swap_min_TLSM[j][nu][3]=nd[3]-D3_HALO;
      }
#endif
    }
  }

  // NOTA BENE: def_axis and def_vec_4d are in PHYSICAL coordinates: 0=x, 1=y, 2=z, 3=t
  // def_vec_4d stores the defect extent along PHYSICAL coordinates, def_axis the boundary on which defect is put in PHYSICAL coords
  // def_vec represents the defect length along the PHYSICAL directions orthogonal to def_axis following the x-y-z-t order.
  // if def_axis=0 (x) => def_vec[0-1-2] = length_def_y-z-t,  if def_axis=1 (y) => def_vec[0-1-2] = length_def_x-z-t and so on...
  // WARNING: This is independent of the physical-logical axis map of the lattice specified in the input file.
  // thus, if in input defect_boundary = 0(1-2-3), def_axis will be put on the x(y-z-t) boundary REGARDLESS of the specified axis mapping

  def_vec_4d[def_axis]=1; // along def boundary, def extent is 1
  for(i=0;i<3;i++){
    def_vec_4d[perp_dir[def_axis][i]]=def_vec[i];
    if(defect_info_config==0){ def->def_mapped_perp_dir[i]=perp_dir[def_axis_mapped][i];}
  }

  // min defect coord along dir: = 0 if dir is orthogonal to def_axis_mapped, otherwise = N_sites_dir - 1
  int x_mind = (geom_par.xmap==def_axis_mapped)? geom_par.gnx-1 : 0;
  int y_mind = (geom_par.ymap==def_axis_mapped)? geom_par.gny-1 : 0;
  int z_mind = (geom_par.zmap==def_axis_mapped)? geom_par.gnz-1 : 0;
  int t_mind = (geom_par.tmap==def_axis_mapped)? geom_par.gnt-1 : 0;

  // max defect coord along every dir
  int x_maxd = x_mind + def_vec_4d[0];
  int y_maxd = y_mind + def_vec_4d[1];
  int z_maxd = z_mind + def_vec_4d[2];
  int t_maxd = t_mind + def_vec_4d[3];
  
  int tnx = geom_par.gnx;
  int tny = geom_par.gny;
  int tnz = geom_par.gnz;
  int tnt = geom_par.gnt;
  int count=0;
  
  int d[4], idxh, x, y, z, t;
	for(d[3]=0; d[3]<nd3; d[3]++)
		for(d[2]=0; d[2]<nd2; d[2]++)
			for(d[1]=0; d[1]<nd1; d[1]++)
				for(d[0]=0; d[0]<nd0; d[0]++){
	  
					idxh = snum_acc(d[0],d[1],d[2],d[3]);
	  
					x = d[geom_par.xmap];
					y = d[geom_par.ymap];
					z = d[geom_par.zmap];
					t = d[geom_par.tmap];

#ifdef MULTIDEVICE
					x+= devinfo.origin_0123[geom_par.xmap]
						- devinfo.halo_widths0123[geom_par.xmap]; // x is now physical x-coordinate for every MPI Rank, same for y,z and t
					y+= devinfo.origin_0123[geom_par.ymap]
						- devinfo.halo_widths0123[geom_par.ymap];
					z+= devinfo.origin_0123[geom_par.zmap]
						- devinfo.halo_widths0123[geom_par.zmap];
					t+= devinfo.origin_0123[geom_par.tmap]
						- devinfo.halo_widths0123[geom_par.tmap];

					if(x>tnx-1)  x-= tnx ; if(x<0)  x+= tnx ; // impose periodic boundary conditions
					if(y>tny-1)  y-= tny ; if(y<0)  y+= tny ;
					if(z>tnz-1)  z-= tnz ; if(z<0)  z+= tnz ;
					if(t>tnt-1)  t-= tnt ; if(t<0)  t+= tnt ;
#endif                  

					// int parity_fis = (x+y+z+t)%2; // physical parity - NOT USED
					int parity_log=(d[0]+d[1]+d[2]+d[3])%2; // logic parity -> USED
					parity = parity_log;

					// initialize K=1 on the whole lattice
					for(mu=0;mu<4;mu++) { conf[2*mu+parity].K.d[idxh] = 1.0; }

					// check if (x,y,z,t) is on defect
					condition = (x >= x_mind) && (y >= y_mind) && (z >= z_mind) && (t >= t_mind) && (x < x_maxd) && (y < y_maxd) && (z < z_maxd) && (t < t_maxd);

#ifdef MULTIDEVICE
					// check if site is NOT on the right halo (must use logic coords to this end)
					condition_2 = ( d[0]<=nd[0]-devinfo.halo_widths0123[0] && d[1]<=nd[1]-devinfo.halo_widths0123[1]
													&& d[2]<=nd[2]-devinfo.halo_widths0123[2] && d[3]<=nd[3]-devinfo.halo_widths0123[3]);
					// check if site is NOT on the left halo (must use logic coords to this end)
					condition_3 = ( d[0]>=devinfo.halo_widths0123[0] && d[1]>=devinfo.halo_widths0123[1]
													&& d[2]>=devinfo.halo_widths0123[2] && d[3]>=devinfo.halo_widths0123[3] );
#endif

					if(condition){
						// find min and max coordinates of sites beloning to the defect and not belonging to the halos
						if(defect_info_config==0){
							for(int nu=0; nu<4;nu++){
								for (i=0; i<4; i++){
									// find logic min coord of the defect along dir i only if site is not on halos
									if ( d[i] < def->defect_swap_min[nu][i] && condition_3 && condition_2 ){
										def->defect_swap_min[nu][i] = d[i];
#ifdef GAUGE_ACT_TLSM
										for(int j=0;j<2;j++) def->defect_swap_min_TLSM[j][nu][i] = d[i];
#endif
									}
									// find logic max coord of the defect along dir i only if site is not on halos
									if ( d[i] > def->defect_swap_max[nu][i] && condition_3 && condition_2 ){
										def->defect_swap_max[nu][i] = d[i];
#ifdef GAUGE_ACT_TLSM
										for(int j=0;j<2;j++) def->defect_swap_max_TLSM[j][nu][i] = d[i];
#endif
									}
								}
							}
						}
						if (condition_3 && condition_2) count++; // count number of sites belonging to the defect and not the halos sitting on the current MPI Rank
						conf[2*def_axis_mapped+parity].K.d[idxh] = c_r; // give the right value of c(r) also to links living on the halos (if the live on the defect)
						if(verbosity_lv>1)
							printf("This site lives on the defect || logic: %d %d %d %d || cart: %d %d %d %d || c(r)=%lf\n",d[0],d[1],d[2],d[3],x,y,z,t,c_r);
					}
				} // close loop on lattice sites

	printf("MPI%02d: num of defect sites on the bulk = %d\n", devinfo.myrank, count);

  if(defect_info_config==0){
    if(count==0){ // if no site (halos excluded) belong the defect for this MPI Rank ==> min=max=0
      for(int nu=0; nu<4;nu++){
				for (i=0; i<4; i++){
					def->defect_swap_min[nu][i] = 0;
					def->defect_swap_max[nu][i] = 0; 
#ifdef GAUGE_ACT_TLSM
					for(int j=0;j<2;j++){
						def->defect_swap_min_TLSM[j][nu][i] = 0;
						def->defect_swap_max_TLSM[j][nu][i] = 0;
					}
#endif
				}
      }
    }	
    else{ // (part of) the defect sits on this MPI Rank
    
      // adjust max extremes to correctly express the for loop guard  
      for(int nu=0; nu<4;nu++){
				for (i=0; i<4; i++){
					def->defect_swap_max[nu][i] += 1; // max must enter in the guard of a for loop as for(i=min; i<(max+1); i++) in order to have i going from min to max
#ifdef GAUGE_ACT_TLSM
					for(int j=0;j<2;j++){
						def->defect_swap_max_TLSM[j][nu][i] += 1; // same as above
					}
#endif
				}
      }

      // Wilson case
      // 1x1 plaq: when computing S_SWAP_Wilson on the plane (mu,nu) (with mu=def_axis), must include contribution from x+nu € def sites => min[nu][nu] -= 1
      def->defect_swap_min[0][0] -= 1;
      def->defect_swap_min[1][1] -= 1;
      def->defect_swap_min[2][2] -= 1;
      if(0==devinfo.myrank){ def->defect_swap_min[3][3] -= 1; } // dir 3 saliminized ==> if nu=3 go 1 step back only for master rank (automatically right for single-device too)
			
      // Symanzik case
#ifdef GAUGE_ACT_TLSM
      // 1x2 rect: include contributions from x+nu € def and x + 2nu € def sites
      def->defect_swap_min_TLSM[0][0][0] -= 2;
      def->defect_swap_min_TLSM[0][1][1] -= 2;
      def->defect_swap_min_TLSM[0][2][2] -= 2;
      if(0==devinfo.myrank){ def->defect_swap_min_TLSM[0][3][3] -= 2; } // dir 3 saliminized ==> if nu=3 go 2 steps back only for master rank (automatically right for single-device too)

      // 2x1 rect: include contributions from x+nu € def (first group of lines) and x + mu € def (second group of lines) sites
      def->defect_swap_min_TLSM[1][0][0] -= 1;
      def->defect_swap_min_TLSM[1][1][1] -= 1;
      def->defect_swap_min_TLSM[1][2][2] -= 1;
      if(0==devinfo.myrank){ def->defect_swap_min_TLSM[1][3][3] -= 1; }
     	
      def->defect_swap_min_TLSM[1][0][def_axis_mapped] -= 1;
      def->defect_swap_min_TLSM[1][1][def_axis_mapped] -= 1;
      def->defect_swap_min_TLSM[1][2][def_axis_mapped] -= 1;
      def->defect_swap_min_TLSM[1][3][def_axis_mapped] -= 1;
#endif
    }

    // print max and mins for each case
    if(verbosity_lv>8){
      printf("rank %d count: %d\n",devinfo.myrank,count);
      // 1x1 plaquettee
      // max
      printf("defect_swap_max: rank %d\n",devinfo.myrank);
      for(int nu=0;nu<4;nu++){
      	if(nu!=def->def_axis_mapped){
					printf("nu:%d||",nu);
					for(i=0;i<4;i++){
						printf("%d||",def->defect_swap_max[nu][i]);
					}
					printf("\n");
      	}
      }
      // min
      printf("defect_swap_min: rank %d\n",devinfo.myrank);
      for(int nu=0;nu<4;nu++){
      	if(nu!=def->def_axis_mapped){
					printf("nu:%d||",nu);
					for(i=0;i<4;i++){
						printf("%d||",def->defect_swap_min[nu][i]);
					}
					printf("\n");
      	}
      }

#ifdef GAUGE_ACT_TLSM
      // 1x2 rectangle
      // max
      printf("defect_swap_max_tlsm 1x2: rank %d\n",devinfo.myrank);
      for(int nu=0;nu<4;nu++){
				if(nu!=def_axis_mapped){
					printf("nu:%d||",nu);
					for(i=0;i<4;i++){
						printf("%d||",def->defect_swap_max_TLSM[0][nu][i]);
					}
					printf("\n");
				}
			}
      // min
      printf("defect_swap_min_tlsm 1x2: rank %d\n",devinfo.myrank);
      for(int nu=0;nu<4;nu++){
				if(nu!=def_axis_mapped){
					printf("nu:%d||",nu);
					for(i=0;i<4;i++){
						printf("%d||",def->defect_swap_min_TLSM[0][nu][i]);
					}
					printf("\n");
				}
			}
      // 2x1 rectangle
      // max
      printf("defect_swap_max_tlsm 2x1: rank %d\n",devinfo.myrank);
      for(int nu=0;nu<4;nu++){
				if(nu!=def_axis_mapped){
					printf("nu:%d||",nu);
					for(i=0;i<4;i++){
						printf("%d||",def->defect_swap_max_TLSM[1][nu][i]);
					}
					printf("\n");
				}
			}
      // min
      printf("defect_swap_min_tlsm 2x1: rank %d\n",devinfo.myrank);
      for(int nu=0;nu<4;nu++){
				if(nu!=def_axis_mapped){
					printf("nu:%d||",nu);
					for(i=0;i<4;i++){
						printf("%d||",def->defect_swap_min_TLSM[1][nu][i]);
					}
					printf("\n");
				}
			}
#endif
    } // close if(verbosity_lv > 8)
  } // close if(defect_info_config==0)
}

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
}

// print conf labels
void label_print(rep_info * hpt_params, FILE *file, int step_number){

  int i;
  fprintf(file,"%d    ",step_number);    
  for(i=0;i<hpt_params->replicas_total_number;i++){
    fprintf(file,"%d    ", hpt_params->label[i]);
  }
  fprintf(file,"\n");
}

double calc_S_soloopenacc_defect(__restrict  su3_soa * const tconf_acc,
																 __restrict su3_soa * const local_plaqs,
																 dcomplex_soa * const tr_local_plaqs,defect_info * def)
{
  double result=0.0;
  double total_result=0.0;
  int mu,nu;
  mu=def->def_axis_mapped;
  int counter=0;

  for(counter=0;counter<3;counter++){
    nu=def->def_mapped_perp_dir[counter];

		
		// result = C_0 * delta_S_plaq_1x1
    result += C_ZERO * calc_S_Wilson_defect(tconf_acc,local_plaqs,tr_local_plaqs,mu,nu,def);
#ifdef GAUGE_ACT_TLSM
		// result + C_1 * ( delta_S_rect_1x2 + delta_S_rect_2x1 )
    result += C_ONE * calc_S_Symanzik_defect(tconf_acc,local_plaqs,tr_local_plaqs,mu,nu,def);
#endif
  }
#ifdef MULTIDEVICE
  MPI_Allreduce((void*)&result,(void*)&total_result,
								1,MPI_DOUBLE,MPI_SUM,devinfo.mpi_comm);
#else
  total_result = result;
#endif
  return total_result;
}

void compute_S_of_replicas(
                        __restrict su3_soa * const tconf_acc, 
                        __restrict su3_soa * const local_plaq, 
                        dcomplex_soa * const tr_local_plaqs,
                        defect_info * def,
                        double * S_arr){
  
  // reduce on all world ranks
  double S_val[NREPLICAS];
  for(int r=0; r<NREPLICAS; ++r){
    S_val[r]=0.0;
  }
  //double S_prev_loc=(double)(devinfo.nranks*devinfo.replica_idx+devinfo.myrank); 
  double S_loc=calc_S_soloopenacc_defect(tconf_acc,local_plaq, tr_local_plaqs,def);

  // master decides swap attempts and send command of actions recomputation
  
  MPI_Barrier(MPI_COMM_WORLD);


  int rep_idx=devinfo.replica_idx;

  MPI_Reduce((void*)&S_loc,(void*)&(S_val[rep->label[rep_idx]]), 1,MPI_DOUBLE, MPI_SUM, 0, devinfo.mpi_comm);
  
  for(int lab=0; lab<NREPLICAS; ++lab){
    MPI_Reduce(&S_val[lab], &S_arr[lab], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
}


// compute S_plaq = (beta/3) * sum K_1x1 Tr(plaq_1x1)
double calc_S_Wilson_defect(__restrict const su3_soa * const u,
                              __restrict su3_soa * const loc_plaq,
                              dcomplex_soa * const tr_local_plaqs,
                              const int mu, const int nu, defect_info * def)
{
  double K_mu_nu;
  double K_mu_nu2;
  int d0,d1,d2,d3;
  int i;

  int D0_min,D1_min,D2_min,D3_min;
  int D0_max,D1_max,D2_max,D3_max;


  D0_min=def->defect_swap_min[nu][0];
  D1_min=def->defect_swap_min[nu][1];
  D2_min=def->defect_swap_min[nu][2];
  D3_min=def->defect_swap_min[nu][3];

  D0_max=def->defect_swap_max[nu][0];
  D1_max=def->defect_swap_max[nu][1];
  D2_max=def->defect_swap_max[nu][2];
  D3_max=def->defect_swap_max[nu][3];

  for(i=0; i<sizeh;i++){
    tr_local_plaqs[0].c[i]=0;
    tr_local_plaqs[1].c[i]=0;
  }
	#pragma acc update device(tr_local_plaqs[0:2])

	#pragma acc kernels present(u) present(loc_plaq) present(tr_local_plaqs)
	#pragma acc loop independent gang(STAPGANG3)
  for(d3=D3_min; d3<D3_max; d3++) {
		#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=D2_min; d2<D2_max; d2++) {
      for(d1=D1_min; d1<D1_max; d1++) {
				for(d0=D0_min; d0<D0_max; d0++) {
                
					int idxh,idxpmu,idxpnu;
					int parity;
					int dir_muA,dir_nuB;
					int dir_muC,dir_nuD;
                
					idxh = snum_acc(d0,d1,d2,d3); // site r
					parity = (d0+d1+d2+d3) % 2;

					// impose periodic boundary conditions
					if(d1==-1){parity = (d0+d2+d3) % 2; idxh=nnm_openacc[snum_acc(d0,0,d2,d3)][1][parity];parity=!parity;}
					if(d2==-1){parity = (d0+d1+d3) % 2; idxh=nnm_openacc[snum_acc(d0,d1,0,d3)][2][parity];parity=!parity;}
					if(d0==-1){parity = (d3+d1+d2) % 2; idxh=nnm_openacc[snum_acc(0,d1,d2,d3)][0][parity];parity=!parity;}
					if(d3==-1){parity = (d0+d1+d2) % 2; idxh=nnm_openacc[snum_acc(d0,d1,d2,0)][3][parity];parity=!parity;}
                
					dir_muA = 2*mu +  parity;
					dir_muC = 2*mu + !parity;
					idxpmu = nnp_openacc[idxh][mu][parity]; // r+mu
                
					dir_nuB = 2*nu + !parity;
					dir_nuD = 2*nu +  parity;
					idxpnu = nnp_openacc[idxh][nu][parity]; // r+nu
                
					//       r+nu (C)  r+mu+nu
					//          +<---+
					// nu       |    ^
					// ^    (D) V    | (B)
					// |        +--->+
					// |       r  (A)  r+mu
					// +---> mu
                
					// plaquette u
					mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);  // loc_plaq = A * B
					mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu);              // loc_plaq = loc_plaq * C
					mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuD],idxh);                // loc_plaq = loc_plaq * D
                
					d_complex ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh); // tr(plaq_u)
					tr_local_plaqs[parity].c[idxh] = creal(ciao)+cimag(ciao)*I;
                
					// K_{mu nu} u
					K_mu_nu=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_muC].K.d[idxpnu])*(u[dir_nuD].K.d[idxh]);

            tr_local_plaqs[parity].c[idxh] = K_mu_nu*tr_local_plaqs[parity].c[idxh];
				} // d0  
      } // d1
    } // d2
  } // d3

  double res_R_p = 0.0;
  double res_I_p = 0.0;
  double resR = 0.0;
  int t;

	#pragma acc kernels present(tr_local_plaqs)
	#pragma acc loop reduction(+:res_R_p) reduction(+:res_I_p)
  for(t=0; t<sizeh; t++) {
    res_R_p += creal(tr_local_plaqs[0].c[t]); // even sites plaquettes
    res_R_p += creal(tr_local_plaqs[1].c[t]); // odd sites plaquettes
  }
  res_R_p *= BETA_BY_THREE; // S = (beta/3) Sum K * Tr_plaq
  return res_R_p;
}


// compute S = (beta/3) Sum [ K_rect_1x2 * Tr(rect)_1x2 + 1x2 ---> 2x1 ]
#ifdef GAUGE_ACT_TLSM
double calc_S_Symanzik_defect(__restrict const su3_soa * const u,
                                  __restrict su3_soa * const loc_rects,
                                  dcomplex_soa * const tr_local_rects,
                                  const int mu, const int nu, defect_info * def)
{
  double K_mu_nu;
  double K_mu_nu2;
  int d0, d1, d2, d3;
  int D0_min,D1_min,D2_min,D3_min;
  int D0_max,D1_max,D2_max,D3_max;

  // min and max init
  D0_max=def->defect_swap_max_TLSM[1][nu][0];
  D1_max=def->defect_swap_max_TLSM[1][nu][1];
  D2_max=def->defect_swap_max_TLSM[1][nu][2];
  D3_max=def->defect_swap_max_TLSM[1][nu][3];

  D0_min=def->defect_swap_min_TLSM[1][nu][0];
  D1_min=def->defect_swap_min_TLSM[1][nu][1];
  D2_min=def->defect_swap_min_TLSM[1][nu][2];
  D3_min=def->defect_swap_min_TLSM[1][nu][3];

  int idxh,idxpmu,idxpnu;
  int parity;
  int dir_muA,dir_nuB;
  int dir_muC,dir_nuD;

  int dir_muB,dir_muD,dir_nuC;
  int dir_muE,dir_nuF,dir_nuE;
  int idxpmupmu,idxpmupnu; // 2x1
  int idxpnupnu; // 1x2

  int is;
  for(is=0; is<sizeh;is++){
    tr_local_rects[0].c[is]=0;
    tr_local_rects[1].c[is]=0;
  }
	#pragma acc update device(tr_local_rects[0:2])

  // 2x1
	#pragma acc kernels present(u) present(loc_rects) present(tr_local_rects)
	#pragma acc loop independent gang(STAPGANG3)
  for(d3=D3_min; d3< D3_max; d3++) {
		#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=D2_min; d2< D2_max; d2++) {
      for(d1=D1_min; d1< D1_max; d1++) {
				for(d0=D0_min; d0< D0_max; d0++) {

					idxh = snum_acc(d0,d1,d2,d3);
					parity = (d0+d1+d2+d3) % 2; 

					if(d1==-1){parity = (d0+d2+d3) % 2; idxh=nnm_openacc[snum_acc(d0,0,d2,d3)][1][parity];parity=!parity;}
					if(d2==-1){parity = (d0+d1+d3) % 2; idxh=nnm_openacc[snum_acc(d0,d1,0,d3)][2][parity];parity=!parity;}
					if(d0==-1){parity = (d3+d1+d2) % 2; idxh=nnm_openacc[snum_acc(0,d1,d2,d3)][0][parity];parity=!parity;}
					if(d3==-1){parity = (d0+d1+d2) % 2; idxh=nnm_openacc[snum_acc(d0,d1,d2,0)][3][parity];parity=!parity;}

					dir_muA = 2*mu +  parity;
					dir_muB = 2*mu + !parity;
					dir_nuC = 2*nu +  parity;
					dir_muD = 2*mu +  parity;
					dir_muE = 2*mu + !parity;
					dir_nuF = 2*nu +  parity;
					idxpmu = nnp_openacc[idxh][mu][parity]; // r+mu
					idxpmupmu = nnp_openacc[idxpmu][mu][!parity]; // r+2mu
					idxpmupnu = nnp_openacc[idxpmu][nu][!parity];// r+mu+nu
					idxpnu = nnp_openacc[idxh][nu][parity]; // r+nu

					//       r+nu r+mu+nu r+2mu+nu
					//          +<---+<---+
					// nu       | (E) (D) ^
					// ^    (F) V (A) (B) |  (C)
					// |        +--->+--->+
					// |       r   r+mu r+2mu
					// +---> mu

					// rect 2x1 u
					mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_muB],idxpmu,&loc_rects[parity],idxh);   // loc_rect = A * B
					mat1_times_mat2_into_mat1_absent_stag_phases(&loc_rects[parity],idxh,&u[dir_nuC],idxpmupmu);                 // loc_rect = loc_rect * C
					mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_rects[parity],idxh,&u[dir_muD],idxpmupnu);            // loc_rect = loc_rect * D
					mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_rects[parity],idxh,&u[dir_muE],idxpnu);               // loc_rect = loc_rect * E
					mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_rects[parity],idxh,&u[dir_nuF],idxh);                 // loc_rect = loc_rect * F
					d_complex ciao = matrix_trace_absent_stag_phase(&loc_rects[parity],idxh);

					// K_mu_nu_RECT 2x1 u
					double K_mu_nu_RECT;
					K_mu_nu_RECT=(u[dir_muA].K.d[idxh])*(u[dir_muB].K.d[idxpmu])*(u[dir_nuC].K.d[idxpmupmu])*(u[dir_muD].K.d[idxpmupnu])*(u[dir_muE].K.d[idxpnu])*(u[dir_nuF].K.d[idxh]);

					tr_local_rects[parity].c[idxh] += K_mu_nu_RECT*( creal(ciao)+cimag(ciao)*I); // K_2x1 * Tr_rect_2x1 = -S_2x1
				} //d0
      } // d1
    } // d2
  } // d3

  double res_R_p = 0.0;
  double res_I_p = 0.0;
  double resR = 0.0;
  int t;

	#pragma acc kernels present(tr_local_rects)
	#pragma acc loop reduction(+:res_R_p) reduction(+:res_I_p)
  for(t=0; t  < sizeh; t++) {
    res_R_p += creal(tr_local_rects[0].c[t]); // even sites rectangles
    res_R_p += creal(tr_local_rects[1].c[t]); // odd sites reactangles
  }

  // 1x2
  // min and max init
  int idxh1;
  D0_max=def->defect_swap_max_TLSM[0][nu][0];
  D1_max=def->defect_swap_max_TLSM[0][nu][1];
  D2_max=def->defect_swap_max_TLSM[0][nu][2];
  D3_max=def->defect_swap_max_TLSM[0][nu][3];

  D0_min=def->defect_swap_min_TLSM[0][nu][0];
  D1_min=def->defect_swap_min_TLSM[0][nu][1];
  D2_min=def->defect_swap_min_TLSM[0][nu][2];
  D3_min=def->defect_swap_min_TLSM[0][nu][3];

  for(is=0; is<sizeh;is++){
    tr_local_rects[0].c[is]=0;
    tr_local_rects[1].c[is]=0;
  }
	#pragma acc update device(tr_local_rects[0:2])

	#pragma acc kernels present(u) present(loc_rects) present(tr_local_rects)
	#pragma acc loop independent gang(STAPGANG3) 
  for(d3=D3_min; d3< D3_max; d3++) {
		#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=D2_min; d2< D2_max; d2++) {
      for(d1=D1_min; d1< D1_max; d1++) {
				for(d0=D0_min; d0< D0_max; d0++) {

					idxh = snum_acc(d0,d1,d2,d3);
					parity = (d0+d1+d2+d3) % 2;

					if(d1==-1){parity = (d0+d2+d3) % 2; idxh=nnm_openacc[snum_acc(d0,0,d2,d3)][1][parity];parity=!parity;}
					if(d2==-1){parity = (d0+d1+d3) % 2; idxh=nnm_openacc[snum_acc(d0,d1,0,d3)][2][parity];parity=!parity;}
					if(d0==-1){parity = (d3+d1+d2) % 2; idxh=nnm_openacc[snum_acc(0,d1,d2,d3)][0][parity];parity=!parity;}
					if(d3==-1){parity = (d0+d1+d2) % 2; idxh=nnm_openacc[snum_acc(d0,d1,d2,0)][3][parity];parity=!parity;}
					int idxh1;
					if(d1==-2){parity = (d0+d2+d3) % 2; idxh1=nnm_openacc[snum_acc(d0,0,d2,d3)][1][parity];idxh=nnm_openacc[idxh1][1][!parity]; }
					if(d2==-2){parity = (d0+d1+d3) % 2; idxh1=nnm_openacc[snum_acc(d0,d1,0,d3)][2][parity];idxh=nnm_openacc[idxh1][2][!parity];}
					if(d0==-2){parity = (d3+d1+d2) % 2; idxh1=nnm_openacc[snum_acc(0,d1,d2,d3)][0][parity];idxh=nnm_openacc[idxh1][0][!parity];}
					if(d3==-2){parity = (d0+d1+d2) % 2; idxh1=nnm_openacc[snum_acc(d0,d1,d2,0)][3][parity];idxh=nnm_openacc[idxh1][3][!parity]; }

					dir_muA = 2*mu +  parity;
					dir_nuB = 2*nu + !parity;
					dir_nuC = 2*nu +  parity;
					dir_muD = 2*mu +  parity;
					dir_nuE = 2*nu + !parity;
					dir_nuF = 2*nu +  parity;

					idxpmu = nnp_openacc[idxh][mu][parity];       // r+mu
					idxpmupnu = nnp_openacc[idxpmu][nu][!parity]; // r+mu+nu
					idxpnu = nnp_openacc[idxh][nu][parity];       // r+nu
					idxpnupnu = nnp_openacc[idxpnu][nu][!parity]; // r+nu+nu
                        
					//            (D)
					//    r+2nu +<---+ r+mu+2nu
					//          |    ^
					//      (E) V    | (C)
					//     r+nu +    + r+mu+nu
					// nu       |    ^
					// ^    (F) V    | (B)
					// |        +--->+
					// |       r  (A)  r+mu
					// +---> mu

					// rect 1x2 u
					mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_rects[parity],idxh);   // loc_rect = A * B
					mat1_times_mat2_into_mat1_absent_stag_phases(&loc_rects[parity],idxh,&u[dir_nuC],idxpmupnu);                 // loc_rect = loc_rect * C
					mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_rects[parity],idxh,&u[dir_muD],idxpnupnu);            // loc_rect = loc_rect * D
					mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_rects[parity],idxh,&u[dir_nuE],idxpnu);               // loc_rect = loc_rect * E
					mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_rects[parity],idxh,&u[dir_nuF],idxh);                 // loc_rect = loc_rect * F
					d_complex ciao = matrix_trace_absent_stag_phase(&loc_rects[parity],idxh);

					// K_mu_nu_RECT3 1x2 u
					double K_mu_nu_RECT3;
					K_mu_nu_RECT3=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_nuC].K.d[idxpmupnu])*(u[dir_muD].K.d[idxpnupnu])*(u[dir_nuE].K.d[idxpnu])*(u[dir_nuF].K.d[idxh]);
	  
					tr_local_rects[parity].c[idxh] += K_mu_nu_RECT3*(creal(ciao)+cimag(ciao)*I); // K_1x2 * Tr_rect_1x2 = - S_1x2
				} // d0
      } // d1
    } // d2
  } // d3

	#pragma acc kernels present(tr_local_rects)
	#pragma acc loop reduction(+:res_R_p) reduction(+:res_I_p)
  for(t=0; t  < sizeh; t++) {
    res_R_p += creal(tr_local_rects[0].c[t]); // even sites rectangles
    res_R_p += creal(tr_local_rects[1].c[t]); // odd sites rectangles
  }
  res_R_p *= BETA_BY_THREE; // S_Symanzik = (beta/3) sum ( K_1x2 * Tr_Rect_1x2 + 1x2 ---> 2x1 )
  return res_R_p;
}
#endif




int metro_SWAP_worldmaster(double Delta_S_SWAP){
  double p1,randoub;
  int accepted = 0;
  if(Delta_S_SWAP < 0) {
    accepted=1;
  } else {
		p1=exp(-Delta_S_SWAP);
		if(debug_settings.do_norandom_test) randoub=0.0; // NORANDOM
		else { // NORMAL, RANDOM
			randoub=casuale();
		}
		if(verbosity_lv>8) printf("p_metro = %.15lg, p_extracted= %.15lg\n", p1, randoub);
		if (randoub<p1) accepted=1;
	}
  return accepted;
}


void manage_replica_swaps(
                    su3_soa * tconf_acc,
										__restrict su3_soa * const loc_plaq,
										dcomplex_soa * const tr_local_plaqs, 
										defect_info * def, 
                    int* swap_num,
										int * all_swap_vet,
										int * acceptance_vet, rep_info * hpt_params){
    double S_arr_prev[NREPLICAS];
    double S_arr_next[NREPLICAS];

    MPI_Barrier(MPI_COMM_WORLD);

    // select relabeling and scatter coefficients to each mpirank
    double swap_order;
    int replicas_number = NREPLICAS;//hpt_params->replicas_total_number;

    if (0==devinfo.myrank_world){
      swap_order=casuale();
    }
    MPI_Bcast((void*)&swap_order,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    int accepted=0;
    int i_counter, j_counter; // indexes, not labels
    int rep_lab1,rep_lab2;
    for(i_counter=0;i_counter<replicas_number-1;i_counter++){
      // manages each pairs serially in a chained fashion

      for(int lab=0; lab<NREPLICAS; ++lab){
        S_arr_prev[lab]=0.0;
      }
      compute_S_of_replicas(tconf_acc, loc_plaq, tr_local_plaqs, def, &S_arr_prev[0]);


      if(swap_order<=0.5){ 
//        MPI_PRINTF0("Swap order: 0->N_r-1\n"); //
        rep_lab1=hpt_params->label[i_counter];
        rep_lab2=hpt_params->label[i_counter+1];
      }else{
//        MPI_PRINTF0("Swap order: N_r-1->0\n");
        rep_lab1=hpt_params->label[replicas_number-i_counter-1];
        rep_lab2=hpt_params->label[replicas_number-i_counter-2];
      }

//        if(verbosity_lv>4) printf("proposing swap %d %d\n",i_counter,i_counter+1);      
      if(devinfo.myrank_world==0){
//        for(int lab=0; lab<NREPLICAS; ++lab){
//          MPI_PRINTF1("i_counter: %d, S_arr_prev[%d]=%lg\n",i_counter, lab,S_arr_prev[lab]);
//        }
        printf("proposing swap %d %d\n",rep_lab1,rep_lab2);      
      }
      
      // rep_lab1 and rep_lab2 are tried for swap
      if(rep_lab1==hpt_params->label[devinfo.replica_idx]){
        // set defect as next
        MPI_PRINTF1("replica lab: %d gets coefficient %lf\n",rep_lab1,hpt_params->cr_vec[rep_lab2]);
        init_k(tconf_acc,hpt_params->cr_vec[rep_lab2],hpt_params->defect_boundary,hpt_params->defect_coordinates,&def,1);
      }
      if(rep_lab2==hpt_params->label[devinfo.replica_idx]){
        // set defect as prev
        MPI_PRINTF1("replica lab: %d gets coefficient %lf\n",rep_lab2,hpt_params->cr_vec[rep_lab1]);
        init_k(tconf_acc,hpt_params->cr_vec[rep_lab1],hpt_params->defect_boundary,hpt_params->defect_coordinates,&def,1);
      }
      //TODO: possibly optimize by updating only defect info
      #pragma acc update device(tconf_acc[0:alloc_info.conf_acc_size])

      MPI_Barrier(MPI_COMM_WORLD);

      for(int lab=0; lab<NREPLICAS; ++lab){
        S_arr_next[lab]=0.0;
      }
      compute_S_of_replicas(tconf_acc, loc_plaq, tr_local_plaqs, def, &S_arr_next[0]);

      if(devinfo.myrank_world==0){
//        for(int lab=0; lab<NREPLICAS; ++lab){
//          MPI_PRINTF1("icounter: %d, S_arr_next[%d]=%lg\n",i_counter, lab,S_arr_next[lab]);
//        }

        MPI_PRINTF1("All actions S_next1, S_next2, S_prev1, S_prev2: %lf, %lf, %lf, %lf\n",
            S_arr_next[rep_lab1], S_arr_next[rep_lab2], 
            S_arr_prev[rep_lab1], S_arr_prev[rep_lab2])
      }
//      // All mpiranks access this  
//      int do_perform_swap=1;

      if (0==devinfo.myrank_world){
        // compute acceptance:
        double Delta_S_SWAP =  S_arr_next[rep_lab1]
                              +S_arr_next[rep_lab2]
                              -S_arr_prev[rep_lab1]
                              -S_arr_prev[rep_lab2];
        accepted=metro_SWAP_worldmaster(Delta_S_SWAP);
        printf("DELTA_S_SWAP(lab1=%d,lab2=%d):%.15lg\n", rep_lab1, rep_lab2, Delta_S_SWAP);
        //if(verbosity_lv>8) printf("DELTA_S_SWAP(r1=%d,r2=%d):%.15lg\n",

        if(accepted){
          // search which indexes are associated the swapped labels
          int aux_label;
          int ii,jj; // indexes
          //TODO: maybe optimize with an inverse label map
          for(int idx=0; idx<replicas_number; ++idx){
            if(hpt_params->label[idx]==rep_lab1){
              ii=idx;
            }else if(hpt_params->label[idx]==rep_lab2){
              jj=idx;
            }
          }
          MPI_PRINTF1("Swap between %d and %d accepted\n",rep_lab1,rep_lab2);
          aux_label=hpt_params->label[ii];
          hpt_params->label[ii] = hpt_params->label[jj];
          hpt_params->label[jj] = aux_label;
        }
      }

      MPI_Bcast((void*)&(hpt_params->label[0]),NREPLICAS,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast((void*)&accepted,1,MPI_INT,0,MPI_COMM_WORLD); // each replica must know if pairs are swapped

      // if swap not accepted, the involved replicas must update their defect information 
      // back to their previous state (prev with prev and next with next)
      if(!accepted && rep_lab1==hpt_params->label[devinfo.replica_idx]){
        // set defect as next
        init_k(tconf_acc,hpt_params->cr_vec[rep_lab1],hpt_params->defect_boundary,hpt_params->defect_coordinates,&def,1);
      }
      if(!accepted && rep_lab2==hpt_params->label[devinfo.replica_idx]){
        // set defect as prev
        init_k(tconf_acc,hpt_params->cr_vec[rep_lab2],hpt_params->defect_boundary,hpt_params->defect_coordinates,&def,1);
      }
      //TODO: possibly optimize by updating only defect info
      #pragma acc update device(tconf_acc[0:alloc_info.conf_acc_size])

//      //XXX: debug!
//      for(int idx=0; idx<rep->replicas_total_number; ++idx){
//        MPI_PRINTF1("labels: label[%d]=%d\n",idx,hpt_params->label[idx]);
//      }
//      MPI_PRINTF1("before swapupdate: swap_num    =%d\n",*swap_num);
//      for(int labp=0; labp<rep->replicas_total_number; ++labp){
//        MPI_PRINTF1("labels: label[%d]=%d\n",labp,hpt_params->label[labp]);
//        MPI_PRINTF1("before swapupdate: all_swap_vet[%d]  =%d (prev: %d, next %d)\n",labp,all_swap_vet[labp],rep_lab1 ,  rep_lab2);
//        if(labp<rep->replicas_total_number-1){
//          MPI_PRINTF1("before swapupdate: acceptance_vet[%d]=%d (prev: %d, next %d)\n",labp,acceptance_vet[labp],rep_lab1 ,  rep_lab2);
//        }
//      }

      *swap_num=*swap_num+1;
      int replicas_num_fixed = (swap_order<=0.5)? rep_lab1 :  rep_lab2;
      all_swap_vet[replicas_num_fixed]++;
      if(accepted==1) acceptance_vet[replicas_num_fixed]++;
      //TODO: these arrays are managed by world master only, they could be deallocated for others

//      //XXX: debug!
//      MPI_PRINTF1("after  swapupdate: swap_num    =%d\n",*swap_num);
//      for(int labp=0; labp<rep->replicas_total_number; ++labp){
//        MPI_PRINTF1("after  swapupdate: all_swap_vet[%d]  =%d (prev: %d, next %d)\n",labp,all_swap_vet[labp],rep_lab1 ,  rep_lab2);
//        MPI_PRINTF1("after  swapupdate: acceptance_vet[%d]=%d (prev: %d, next %d)\n",labp,acceptance_vet[labp],rep_lab1 ,  rep_lab2);
//      }

      MPI_Barrier(MPI_COMM_WORLD); 
    }
}

    
void trasl_conf( __restrict const su3_soa *  const tconf_acc,
								 __restrict const su3_soa *  const taux_conf){
    
#ifdef MULTIDEVICE
  communicate_su3_borders(tconf_acc, GAUGE_HALO);
	#pragma acc update self(tconf_acc[0:8])
#endif
    
  set_su3_soa_to_su3_soa(tconf_acc,taux_conf);
    
  double dir0;
  int dir=0;
    
  if(0==devinfo.myrank){dir0=casuale();}
#ifdef MULTIDEVICE
  MPI_Bcast((void*) &dir0,1,MPI_DOUBLE,0,devinfo.mpi_comm);
  if(verbosity_lv>4)
    printf("MPI%02d dir0 : %f \n",devinfo.myrank,dir0);
#endif

  if(dir0<=0.25){dir=0;}
  if(dir0>0.25 && dir0<=0.5){dir=1;}
  if(dir0>0.5 && dir0<=0.75){dir=2;}
  if(dir0>0.75){dir=3;}
    
  set_su3_soa_to_su3_soa_trasl( taux_conf,tconf_acc, dir);
	#pragma acc update device(tconf_acc[0:8])  
    
#ifdef MULTIDEVICE
  communicate_su3_borders(tconf_acc, GAUGE_HALO);  
	#pragma acc update self(tconf_acc[0:8])
#endif
}

#endif


#endif
