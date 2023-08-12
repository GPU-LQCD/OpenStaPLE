#ifndef TOPOLOGICAL_FORCE_C
#define TOPOLOGICAL_FORCE_C
#include "../Include/common_defines.h"
#include "./topological_force.h"
#include <math.h>
#include "../Mpi/multidev.h"
#include "../Meas/gauge_meas.h"
#include "./fermion_force.h"
#include "./fermion_force_utilities.h"
#include "./alloc_vars.h"
#include "./action.h"
#include "./topological_action.h"
#include "./single_types.h"
#include "./stouting.h"
#include "../DbgTools/dbgtools.h"

#ifdef MULTIDEVICE
#include "../Mpi/communications.h"
#endif 

// #define DEBUG_LOLLO
#ifdef DEBUG_LOLLO
#include "./cayley_hamilton.h"
#include "./alloc_vars.h"

#define ALLOLLOCHECK(control_int,var)  if(control_int != 0 )						\
    MPI_PRINTF1("\tError in  allocation of %s . \n", #var); \
	else if(verbosity_lv > 2) MPI_PRINTF1("\tAllocation of %s : OK , %p\n", #var, var );

#define FREELOLLOCHECK(var) if(verbosity_lv >2)					\
    printf("\tFreed %s, %p ...", #var,var);							\
	free(var); if(verbosity_lv > 2)  printf(" done.\n");

#endif

double compute_topodynamical_potential_der(const double Q)
{
  double barr=act_params.barrier, width=act_params.width;
  int ngrid= (int) floor((2*barr+width/2)/width);
  int kgridfloat = (Q+barr)/width;
  int kgrid= (int) floor(kgridfloat);
  
  if(kgrid>=0 && kgrid<ngrid)
    {
      double grid[ngrid];
      if(verbosity_lv>3)
				printf("\t\t\tMPI%02d - load_topo(path,barr,width,grid,ngrid)\n",devinfo.myrank);
      
      load_topo(act_params.topo_file_path,barr,width,grid,ngrid);
      double dVdQ=(grid[kgrid+1]-grid[kgrid])/width;
      printf("dtopotential/dQ = %f, kgrid = %d, width = %f, Q = %f\n",dVdQ,kgrid,width,Q);
      return dVdQ;
    }
  else{
	  printf("dtopotential/dQ = %f, kgrid = %d, width = %f, Q = %f\n",0.0,kgrid,width,Q);
    return 0;}
}

void four_leaves(__restrict su3_soa * const leaves,__restrict const su3_soa * const u, int nnp_openacc[sizeh][4][2],int nnm_openacc[sizeh][4][2])
{
  int mu, nu;
  int d0, d1, d2, d3;
	#pragma acc kernels present(u) present(leaves) present(nnp_openacc) present(nnm_openacc)
  for(mu=0;mu<3;mu++)
    for(nu=mu+1;nu<4;nu++){
			#pragma acc loop independent gang
      for(d3 = D3_HALO; d3<nd3-D3_HALO; d3++){
				#pragma acc loop independent gang vector
				for(d2=0; d2<nd2; d2++){
					#pragma acc loop independent gang vector
					for(d1=0; d1<nd1; d1++){
						#pragma acc loop independent vector 
						for(d0=0; d0 < nd0; d0++){
							const int idxh = snum_acc(d0,d1,d2,d3);
							const int parity = (d0+d1+d2+d3)%2;
							const int idx_plane = 2*mu + nu - 1 - ((int)mu/2);// mu=0 nu=1 --> idx_plane=0;... mu=1 nu=2 --> idx_plane=3;... mu=2 nu=3 --> idx_plane=5.
							// leaves points
							int idx_p_mu = nnp_openacc[idxh][mu][parity];
							int idx_p_nu = nnp_openacc[idxh][nu][parity];
							int idx_m_mu = nnm_openacc[idxh][mu][parity];
							int idx_m_nu = nnm_openacc[idxh][nu][parity];
	      
							int idx_p_nu_m_mu = nnm_openacc[idx_p_nu][mu][!parity];
							int idx_m_mu_m_nu = nnm_openacc[idx_m_mu][nu][!parity];
							int idx_m_nu_p_mu = nnp_openacc[idx_m_nu][mu][!parity];
							
							// first leave
							comp_U_U_Udag_Udag(&u[2*mu+parity] , idxh,
																 &u[2*nu+!parity], idx_p_mu,
																 &u[2*mu+!parity], idx_p_nu,
																 &u[2*nu+parity] , idxh,
																 &leaves[2*idx_plane+parity], idxh);
							comp_and_add_U_Udag_Udag_U(&u[2*nu+parity] , idxh,
																				 &u[2*mu+parity] , idx_p_nu_m_mu,
																				 &u[2*nu+!parity], idx_m_mu,
																				 &u[2*mu+!parity], idx_m_mu,
																				 &leaves[2*idx_plane+parity], idxh);
							comp_and_add_Udag_Udag_U_U(&u[2*mu+!parity], idx_m_mu,           
																				 &u[2*nu+parity] , idx_m_mu_m_nu,  
																				 &u[2*mu+parity] , idx_m_mu_m_nu,      
																				 &u[2*nu+!parity], idx_m_nu,      
																				 &leaves[2*idx_plane+parity], idxh);
							comp_and_add_Udag_U_U_Udag(&u[2*nu+!parity], idx_m_nu,           
																				 &u[2*mu+!parity], idx_m_nu,  
																				 &u[2*nu+parity] , idx_m_nu_p_mu,      
																				 &u[2*mu+parity] , idxh,      
																				 &leaves[2*idx_plane+parity], idxh);
						} // d0
					} // d1
				} // d2
      } // d3
    } // nu
}

void antihermatize_unsafe(__restrict su3_soa * const leaves)
{
  int d0,d1,d2,d3,plane;

	#pragma acc kernels present(leaves)
	#pragma acc loop independent vector
  for(plane=0; plane<12; plane++)
		#pragma acc loop independent gang(STAPGANG3)
    for(d3 = D3_HALO; d3<nd3-D3_HALO; d3++)
			#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
      for(d2=0; d2<nd2; d2++)
				for(d1=0; d1<nd1; d1++)
					for(d0=0; d0<nd0; d0++)
						{
							const int idx = snum_acc(d0,d1,d2,d3);
	      
							leaves[plane].r0.c0[idx] = 2* cimag(leaves[plane].r0.c0[idx]) * _Complex_I;
							leaves[plane].r1.c1[idx] = 2* cimag(leaves[plane].r1.c1[idx]) * _Complex_I;
							leaves[plane].r2.c2[idx] = 2* cimag(leaves[plane].r2.c2[idx]) * _Complex_I;
	      
							d_complex tmp01 = conj(leaves[plane].r1.c0[idx]);
							d_complex tmp02 = conj(leaves[plane].r2.c0[idx]);
							d_complex tmp12 = conj(leaves[plane].r2.c1[idx]);
	      
							leaves[plane].r1.c0[idx] = (leaves[plane].r1.c0[idx] - conj(leaves[plane].r0.c1[idx]));
							leaves[plane].r2.c0[idx] = (leaves[plane].r2.c0[idx] - conj(leaves[plane].r0.c2[idx]));
							leaves[plane].r2.c1[idx] = (leaves[plane].r2.c1[idx] - conj(leaves[plane].r1.c2[idx]));
	      
							leaves[plane].r0.c1[idx] = (leaves[plane].r0.c1[idx] - tmp01);
							leaves[plane].r0.c2[idx] = (leaves[plane].r0.c2[idx] - tmp02);
							leaves[plane].r1.c2[idx] = (leaves[plane].r1.c2[idx] - tmp12);
	      
						} // closing all the loops at once
}

void topo_staples(__restrict su3_soa * u,__restrict su3_soa * const staples, double norm,int nnp_openacc[sizeh][4][2],int nnm_openacc[sizeh][4][2])
{
  // compute leaves
  su3_soa * leaves;
  posix_memalign((void **)&leaves,128,12*sizeof(su3_soa));
	#pragma acc enter data create(leaves[0:12]) copyin(norm)

  if(verbosity_lv>3)
    printf("\t\t\tMPI%02d - four_leaves(leaves,u)\n",devinfo.myrank);
  
  // compute leaves
#if NRANKS_D3 > 1
	communicate_su3_borders(u, GAUGE_HALO);  
#endif

  four_leaves(leaves,u, nnp_openacc,nnm_openacc);

  
  if(verbosity_lv>3)
    printf("\t\t\tMPI%02d - antihermatize_usafe(leaves)\n",devinfo.myrank);
  
  // antihermatize leaves
  
  antihermatize_unsafe(leaves);
  
#if NRANKS_D3 > 1
	communicate_gl3_borders(leaves, GAUGE_HALO);  
#endif

  int d0, d1, d2, d3;
  if(verbosity_lv>3)
    {
      printf("MPI%d - computing staples\n",devinfo.myrank);
    }  
  
  for(int mu=0; mu<4; mu++)
    for(int inu=0; inu<3; inu++){
			#pragma acc parallel loop independent gang present(leaves) present(u) present(nnp_openacc) present(nnm_openacc) present(staples)
      for(d3 = D3_HALO; d3<nd3-D3_HALO; d3++){
				#pragma acc loop independent gang vector
				for(d2=0; d2<nd2; d2++){
					#pragma acc loop independent gang vector
					for(d1=0; d1<nd1; d1++){
						#pragma acc loop independent gang vector
						for(d0=0; d0<nd0; d0++){
							single_su3 ABC, BCF, ABCF;
							single_su3 ADE, DEF, ADEF;
							single_su3 temp;
							single_su3 loc_stap;
							set_to_zero_single_su3(&loc_stap);
	      
							int perp_dir[4][3] = {{ 1, 2, 3}, { 0, 2, 3}, { 0, 1, 3}, { 0, 1, 2}};
							int plan_perp[4][3]= {{ 5, 4, 3}, { 5, 2, 1}, { 4, 2, 0}, { 3, 1, 0}};
							int plan_sign[4][3]= {{+1,-1,+1}, {-1,+1,-1}, {+1,-1,+1}, {-1,+1,-1}};
	      
							const int nu = perp_dir[mu][inu];       // E---F---C   
							const int idxA = snum_acc(d0,d1,d2,d3); // |   |   | mu
							const int parity = (d0+d1+d2+d3)%2;     // D---A---B
							const int iplan = plan_perp[mu][inu];   //       nu
	      
							const int idxF = nnp_openacc[idxA][mu][parity];
							const int idxB = nnp_openacc[idxA][nu][parity];
							const int idxC = nnp_openacc[idxB][mu][!parity];
							const int idxD = nnm_openacc[idxA][nu][parity];
							const int idxE = nnp_openacc[idxD][mu][!parity];
		  
							// compute ABC BCF ABCF
							su3_soa_times_su3_soa_into_single_su3(&u[2*nu+parity],  idxA,
																										&u[2*mu+!parity], idxB,
																										&ABC);
							su3_soa_times_su3_soa_dag_into_single_su3(&u[2*mu+!parity], idxB,
																												&u[2*nu+!parity], idxF,
																												&BCF);
							su3_soa_times_single_su3_into_single_su3(&u[2*nu+parity],idxA,
																											 &BCF, &ABCF);
							// compute ADE DEF ADEF
							su3_soa_dag_times_su3_soa_into_single_su3(&u[2*nu+!parity], idxD,
																												&u[2*mu+!parity], idxD,
																												&ADE);
							su3_soa_times_su3_soa_into_single_su3(&u[2*mu+!parity], idxD,
																										&u[2*nu+parity],  idxE,
																										&DEF);
							single_su3_times_su3_soa_into_single_su3(&ADE,
																											 &u[2*nu+parity], idxE,
																											 &ADEF);
							// compute local staples.
							gl3_soa_times_single_su3_addto_gl3(&leaves[2*iplan+parity],idxA,
																								 &ABCF,&loc_stap);
							su3_soa_times_gl3_soa_into_gl3(&u[2*nu+parity],idxA,
																						 &leaves[2*iplan+!parity],idxB,
																						 &temp);
							single_gl3xsu3_add_to_out(&loc_stap,&temp,&BCF);
							gl3_soa_times_su3_soa_dag_into_gl3(&leaves[2*iplan+parity],idxC,
																								 &u[2*nu+!parity],       idxF,
																								 &temp);
              single_su3xgl3_add_to_out(&loc_stap,&ABC,&temp);
							single_su3_times_gl3_soa_addto_gl3(&ABCF,
																								 &leaves[2*iplan+!parity],idxF,
																								 &loc_stap);
							// ^-- Right staple
							// v-- Left staple
							gl3_soa_dag_times_single_su3_addto_gl3(&leaves[2*iplan+parity],idxA,
																										 &ADEF,&loc_stap);
							su3_soa_dag_times_gl3_soa_dag_into_gl3(&u[2*nu+!parity],        idxD,
																										 &leaves[2*iplan+!parity],idxD,
																										 &temp);
							single_gl3xsu3_add_to_out(&loc_stap,&temp,&DEF);
							gl3_soa_dag_times_su3_soa_into_gl3(&leaves[2*iplan+parity],idxE,
																								 &u[2*nu+parity],        idxE,
																								 &temp);
							single_su3xgl3_add_to_out(&loc_stap,&ADE,&temp);
							single_su3_times_gl3_soa_dag_addto_gl3(&ADEF,
																										 &leaves[2*iplan+!parity],idxF,
																										 &loc_stap);

							// finalizing
							if(plan_sign[mu][inu] >= 0)
								single_gl3_addinto_su3_soa(&staples[2*mu+parity], idxA,&loc_stap);
							else
								single_gl3_subtinto_su3_soa(&staples[2*mu+parity], idxA,&loc_stap);		
	      
						} // d0
					} // d1
				} // d2
      } // d3
    } // inu

	#pragma acc exit data delete(leaves)
  free(leaves);

	#pragma acc parallel loop independent gang present(staples) present(norm)
  for(d3 = D3_HALO; d3<nd3-D3_HALO; d3++){
		#pragma acc loop independent gang vector
    for(d2=0; d2<nd2; d2++){
			#pragma acc loop independent gang vector
      for(d1=0; d1<nd1; d1++){
				#pragma acc loop independent gang vector  
				for(d0=0; d0<nd0; d0++){
	
					const int idxA = snum_acc(d0,d1,d2,d3);
					const int parity = (d0+d1+d2+d3)%2;
	
					for(int mu = 0; mu<4; mu++)
						gl3_dag_times_double_factor(&staples[2*mu+parity],idxA,norm);

      
				}
      }
    }	    
  }
	#pragma acc exit data delete(norm)
}



void calc_loc_topo_staples(__restrict su3_soa * u,
													 __restrict su3_soa * const staples,
													 int nnp_openacc[sizeh][4][2],int nnm_openacc[sizeh][4][2])
{    
  if(verbosity_lv>3)
    printf("\t\tMPI%02d - compute_topological_charge(u,quadri,loc_q)\n",devinfo.myrank);

  double Q = compute_topological_charge(u, staples, topo_loc, nnp_openacc, nnm_openacc);

  if(verbosity_lv>4)
    printf("Topological Charge: %lf\n",Q);
  
  
  if(verbosity_lv>3)
    printf("\t\tMPI%02d - compute_topodynamical_potential_der(Q)\n",devinfo.myrank);
	
  double pot_der=compute_topodynamical_potential_der(Q);
  double norm=pot_der/(M_PI*M_PI*128);
  
  
  if(verbosity_lv>3)
    printf("\t\tMPI%02d - topo_staples(u, staples, norm)\n",devinfo.myrank);
  

  set_su3_soa_to_zero(staples);
  topo_staples(u,staples,norm,nnp_openacc,nnm_openacc);
}


void calc_ipdot_topo(__restrict su3_soa const * const tconf_acc,  
#ifdef STOUT_TOPO
										 __restrict su3_soa * const tstout_conf_acc_arr,
#endif
										 __restrict su3_soa  * const taux_conf_acc,
										 __restrict su3_soa * const local_staples,
										 __restrict tamat_soa * const tipdot_acc,
										 int nnp_openacc[sizeh][4][2],int nnm_openacc[sizeh][4][2])
{
  set_su3_soa_to_zero(local_staples);
  __restrict su3_soa * conf_to_use; // configuration to use in calculation of topological force
  //  __restrict su3_soa_f * conf_to_use_f; // conf to use in calculation of topo force


#ifdef STOUT_TOPO
  if(act_params.topo_stout_steps > 0){
    stout_wrapper(tconf_acc,tstout_conf_acc_arr,nnp_openacc,nnm_openacc,1);
    conf_to_use = &(tstout_conf_acc_arr[8*(act_params.topo_stout_steps-1)]);
  }
  else conf_to_use=(su3_soa*)tconf_acc;
#else
  conf_to_use=(su3_soa*)tconf_acc;
#endif

  calc_loc_topo_staples(conf_to_use,local_staples,nnp_openacc,nnm_openacc);

#ifdef STOUT_TOPO
  for(int stout_level = act_params.topo_stout_steps; stout_level > 1; stout_level--){
    conf_to_use = &(tstout_conf_acc_arr[8*(stout_level-2)]);

    compute_sigma_from_sigma_prime_backinto_sigma_prime(local_staples,aux_th,aux_ta,
																												conf_to_use,taux_conf_acc,nnp_openacc,nnm_openacc,1);
  }    
  if(act_params.topo_stout_steps > 0 ){
    compute_sigma_from_sigma_prime_backinto_sigma_prime(local_staples,aux_th,aux_ta,
																												tconf_acc, taux_conf_acc,nnp_openacc,nnm_openacc,1);

  }
#endif

  conf_times_staples_ta_part_addto_tamat(tconf_acc,local_staples,tipdot_acc);  


	// DEBUG
#ifdef DEBUG_LOLLO
	if(devinfo.myrank == 0){
		#pragma acc update host(tconf_acc[0:8]) host(local_staples[0:8]) host(conf_to_use[0:8]) host(tipdot_acc[0:8])
		int O = snum_acc(0,0,0,D3_HALO); 
		printf("check multiplication into tamat:\n");
		STAMPA_DEBUG_SU3_SOA(tconf_acc,0,O);
		STAMPA_DEBUG_SU3_SOA(conf_to_use,0,O);
		STAMPA_DEBUG_SU3_SOA(local_staples,0,O);
		STAMPA_DEBUG_TAMAT_SOA(tipdot_acc,0,O);
	}

#define SQRT_3 1.732050807568877
  // i*Gell-mann matrices as from eq.A.10 of Gattringer - note that T=lambda/2
  single_tamat i_gell_mann_matr[8]={ 
																		{ 0+1*I, 0+0*I, 0+0*I, 0, 0 },
																		{ 1+0*I, 0+0*I, 0+0*I, 0, 0 },
																		{ 0+0*I, 0+0*I, 0+0*I, 1, -1 },
																		{ 0+0*I, 0+1*I, 0+0*I, 0, 0 },
																		{ 0+0*I, 1+0*I, 0+0*I, 0, 0 },
																		{ 0+0*I, 0+0*I, 0+1*I, 0, 0 },
																		{ 0+0*I, 0+0*I, 1+0*I, 0, 0 },
																		{ 0+0*I, 0+0*I, 0+0*I, 1/SQRT_3, 1/SQRT_3 }
  };
  double eps=1e-3;
  
  // store initial link and comp action
  single_su3 sto;
	#pragma acc update host(tconf_acc[0:8])
  single_su3_from_su3_soa(&tconf_acc[0],0,&sto);
  rebuild3row(&sto);
  
#ifdef STOUT_TOPO
	double ori_act = compute_topo_action(tconf_acc,tstout_conf_acc_arr,nnp_openacc,nnm_openacc);
#else
	double ori_act = compute_topo_action(tconf_acc,nnp_openacc,nnm_openacc);
#endif
  // store derivative
  single_tamat posi={ 0+0*I, 0+0*I, 0+0*I, 0, 0 };
  single_tamat nega={ 0+0*I, 0+0*I, 0+0*I, 0, 0 };
  
  for(int igen=0;igen<8;igen++)
    {	    
      // prepare increment and change
      single_tamat ba;
      single_tamat_times_scalar_into_tamat(&ba,&i_gell_mann_matr[igen],eps/2+0*I);
      
      single_su3 exp_mod;
      CH_exponential_antihermitian_nissalike(&exp_mod,&ba); // exp( -i eps Gell-Man/2 )
      rebuild3row(&exp_mod);
      
      // change -, compute action
      single_su3 shilink;
      single_su3xsu3(&shilink, &exp_mod, &sto);
      single_su3_into_su3_soa(&tconf_acc[0], 0, &shilink);
      
			#pragma acc update device(tconf_acc[0:8])      
#ifdef STOUT_TOPO
      double act_minus = compute_topo_action(tconf_acc,tstout_conf_acc_arr,nnp_openacc,nnm_openacc);
#else
      double act_minus = compute_topo_action(tconf_acc,nnp_openacc,nnm_openacc);
#endif
      // change +, compute action
      gl3_dagger(&exp_mod); // exp( i eps Gell-Mann/2 )
      single_su3xsu3(&shilink, &exp_mod, &sto);
      single_su3_into_su3_soa(&tconf_acc[0], 0, &shilink);
      
			#pragma acc update device(tconf_acc[0:8])      
#ifdef STOUT_TOPO
      double act_plus = compute_topo_action(tconf_acc,tstout_conf_acc_arr,nnp_openacc,nnm_openacc);
#else
      double act_plus = compute_topo_action(tconf_acc,nnp_openacc,nnm_openacc);
#endif

      // set back everything
      single_su3_into_su3_soa(&tconf_acc[0], 0, &sto);

			#pragma acc update device(tconf_acc[0:8])
#ifdef STOUT_TOPO
      double check_ori_act = compute_topo_action(tconf_acc,tstout_conf_acc_arr,nnp_openacc,nnm_openacc);
#else
      double check_ori_act = compute_topo_action(tconf_acc,nnp_openacc,nnm_openacc);
#endif
      /*
				printf("plus = %+016.016le ",act_plus);
				printf("ori = %+016.016le ",check_ori_act);
				printf("minus = %+016.016le\n",act_minus);
			*/
	    
      // double gr_plus = -(act_plus-ori_act)/eps;
      // double gr_minus = -(ori_act-act_minus)/eps;
      double gr_plus  =  (act_plus  - ori_act)/eps; // [S(U+eps)-S(U)]/eps = dS/dU + O(eps)
			double gr_minus = -(act_minus - ori_act)/eps; // -[S(U-eps)-S(U)]/eps = dS/dU + O(eps)
      single_tamat_times_scalar_add_to_tamat(&posi,&i_gell_mann_matr[igen],gr_plus+0*I);
      single_tamat_times_scalar_add_to_tamat(&nega,&i_gell_mann_matr[igen],gr_minus+0*I);

      /* printf("gr_plus = %+016.016le\ngr_minus = %+016.016le\n", gr_plus, gr_minus); */
    }
  // take the average
  single_tamat Numerical_derivative;
  summ_single_tamats_times_scalar(&Numerical_derivative, &posi, &nega, 0.5+0*I);
  STAMPA_DEBUG_SINGLE_TAMAT(Numerical_derivative);
	
  printf("If they are similar, be happy\n");
  mem_free_core();
  mem_free_extended();

#ifdef MULTIDEVICE
	shutdown_multidev();
#endif
  exit(0);
#endif

}

#endif
