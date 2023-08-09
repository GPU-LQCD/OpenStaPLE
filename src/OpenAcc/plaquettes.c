#ifndef PLAQUETTES_C_
#define PLAQUETTES_C_

#include "./geometry.h"
#include "./plaquettes.h"
#include "./su3_utilities.h"

// JUST FOR DEBUGGING PURPOSES -> TO BE REMOVED
// #include "../DbgTools/dbgtools.h"

// routine for the computation of the average of the plaquettes computed on the plane mu-nu
// 1) all the plaquettes on the plane mu-nu are computed and saved locally
// 2) finally the reduction of the traces is performed
double calc_loc_plaquettes_nnptrick(__restrict const su3_soa * const u,
																		__restrict su3_soa * const loc_plaq, 
																		dcomplex_soa * const tr_local_plaqs,
																		const int mu, const int nu)
{
	#pragma acc kernels present(u) present(loc_plaq) present(tr_local_plaqs)
	#pragma acc loop independent gang(STAPGANG3)
  for(int d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
		#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(int d2=0; d2<nd2; d2++) {
      for(int d1=0; d1<nd1; d1++) {
				for(int d0=0; d0 < nd0; d0++) {
              
					int idxh,idxpmu,idxpnu;
					int parity;
					int dir_muA,dir_nuB;
					int dir_muC,dir_nuD;
	    
					idxh = snum_acc(d0,d1,d2,d3);
					parity = (d0+d1+d2+d3) % 2;
	    
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
	    
					mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh); // loc_plaq = A * B
					mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu);             // loc_plaq = loc_plaq * C
					mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuD],idxh);               // loc_plaq = loc_plaq * D
	    
					d_complex tmp_aux = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
					tr_local_plaqs[parity].c[idxh] = creal(tmp_aux)+cimag(tmp_aux)*I;
		
#ifdef PAR_TEMP
					// K_mu_nu computation;
					double K_mu_nu=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_muC].K.d[idxpnu])*(u[dir_nuD].K.d[idxh]);
					tr_local_plaqs[parity].c[idxh] *= K_mu_nu;
#endif
		
				} // d0
      } // d1
    } // d2
  } // d3
  
  double res_R_p = 0.0;
  double res_I_p = 0.0;
  double resR = 0.0;
  int t; // WARNING: only good for 1D cut
	#pragma acc kernels present(tr_local_plaqs)
	#pragma acc loop reduction(+:res_R_p) reduction(+:res_I_p)
  for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
    res_R_p += creal(tr_local_plaqs[0].c[t]); // even sites plaquettes
    res_R_p += creal(tr_local_plaqs[1].c[t]); // odd sites plaquettes
  }
  return res_R_p;
} // closes routine

// routine to compute the staples for each site on a given plane mu-nu and sum the result to the local stored staples
void calc_loc_staples_nnptrick_all(__restrict const su3_soa * const u,
																	 __restrict su3_soa * const loc_stap)
{
  //       r+mu-nu  r+mu   r+mu+nu
  //          +<-----+----->+
  //          |  1L  ^  1R  |
  // mu    2L |      |      | 2R
  // ^        V  3L  |  3R  V
  // |        +----->+<-----+
  // |       r-nu    r     r+nu
  // +---> nu       
  //            r is idxh in the following      

	#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
	#pragma acc loop independent gang(STAPGANG3)
  for(int d3=D3_HALO; d3<nd3-D3_HALO; d3++){
		#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(int d2=0; d2<nd2; d2++){
      for(int d1=0; d1<nd1; d1++){
				for(int d0=0; d0<nd0; d0++){
					const int idxh = snum_acc(d0,d1,d2,d3); // r
					const int parity = (d0+d1+d2+d3) % 2;
					#pragma acc loop seq
					for(int mu=0; mu<4; mu++){
						const int dir_link = 2*mu + parity;
						#pragma acc loop seq
						for(int iter=0; iter<3; iter++){
							int perp_dirs[4][3] = { {1,2,3}, {0,2,3}, {0,1,3}, {0,1,2} };
							int nu = perp_dirs[mu][iter];
							/*if (mu==0) { nu = iter + 1; }
								else if (mu==1) { nu = iter + (iter & 1) + (iter >> 1); }
								else if (mu==2) { nu = iter + (iter >> 1); }
								else if (mu==3) { nu = iter; }
								else { printf("NU ERROR!\n"); }*/
            
							#pragma acc cache (nnp_openacc[idxh:8])

							const int dir_mu_2R = 2*mu + !parity;
							const int dir_mu_2L = 2*mu + !parity;
							const int idx_pmu = nnp_openacc[idxh][mu][parity]; // r+mu
							#pragma acc cache (nnm_openacc[idx_pmu:8])

							const int dir_nu_1R = 2*nu + !parity;
							const int dir_nu_3R = 2*nu +  parity;
							const int dir_nu_1L = 2*nu +  parity;
							const int dir_nu_3L = 2*nu + !parity;
							const int idx_pnu = nnp_openacc[idxh][nu][parity]; // r+nu

							// computation of the Right part of the staple
							// N.B.: inside this function, each link constituing the staple is multiplied by its K_mu(x) factor
							mat1_times_conj_mat2_times_conj_mat3_addto_mat4_absent_stag_phases(&u[dir_nu_1R], idx_pmu,
																																								 &u[dir_mu_2R], idx_pnu,
																																								 &u[dir_nu_3R], idxh,
																																								 &loc_stap[dir_link], idxh);

							const int idx_mnu = nnm_openacc[idxh][nu][parity]; // r-nu
							const int idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity]; // r+mu-nu

							// computation of the Left part of the staple
							// N.B.: also here, each link of the staple is multiplied by its K_mu(x) factor
							conj_mat1_times_conj_mat2_times_mat3_addto_mat4_absent_stag_phases(&u[dir_nu_1L], idx_pmu_mnu,
																																								 &u[dir_mu_2L], idx_mnu,
																																								 &u[dir_nu_3L], idx_mnu,
																																								 &loc_stap[dir_link], idxh);
						} // iter
						// in the end, each staple has to be multiplied by the K_mu(x) of the link emanating the staples
						// N.B.: This factor is common both for the left and for the right staples, and common for all values of nu
						// double K_mu = u[dir_link].K.d[idxh];
						// gl3_times_double_factor(&(loc_stap[dir_link]), idxh, K_mu);
					} // mu
				} // d0
      } // d1
    } // d2
  } // d3

  /*
  // DEBUG HMC WITH PARALLEL TEMPERING
  #pragma acc update host(loc_stap[0:8])
  int d0=0, d1=15, d2=0, d3=0;
  int mu=1;

  int idxh=snum_acc(d0,d1,d2,d3);
  int parity=(d0+d1+d2+d3)%2;
  int dir_link=2*mu+parity;

  printf("STAPLE WITH DEFECT\n");	
  STAMPA_DEBUG_SU3_SOA(loc_stap,dir_link,idxh);

  //#########################################

  d0=0, d1=13, d2=0, d3=0;
  mu=1;

  idxh=snum_acc(d0,d1,d2,d3);
  parity=(d0+d1+d2+d3)%2;
  dir_link=2*mu+parity;

  printf("STAPLE WITHOUT DEFECT\n");	
  STAMPA_DEBUG_SU3_SOA(loc_stap,dir_link,idxh);
  
  //##########################################

  printf("Ringrazia Zeb89 se si assomigliano almeno un po'\n");
  mem_free_core();
  mem_free_extended();
  #ifdef MULTIDEVICE
  shutdown_multidev();
  #endif
  exit(1);
  */
} // closes routine

// routine to compute the staples for each site on a given plane mu-nu and sum the result to the local stored staples
void calc_loc_staples_nnptrick_all_onlyferms(__restrict const su3_soa * const u,
																						 __restrict su3_soa * const loc_stap)
{
  //       r+mu-nu  r+mu   r+mu+nu
  //          +<-----+----->+
  //          |  1L  ^  1R  |
  // mu    2L |      |      | 2R
  // ^        V  3L  |  3R  V
  // |        +----->+<-----+
  // |       r-nu    r     r+nu
  // +---> nu       
  //            r is idxh in the following      

	#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
	#pragma acc loop independent gang(STAPGANG3)
  for(int d3=D3_HALO; d3<nd3-D3_HALO; d3++){
		#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(int d2=0; d2<nd2; d2++){
      for(int d1=0; d1<nd1; d1++){
				for(int d0=0; d0<nd0; d0++){
					const int idxh = snum_acc(d0,d1,d2,d3); // r
					const int parity = (d0+d1+d2+d3) % 2;
					#pragma acc loop seq
					for(int mu=0; mu<4; mu++){
						const int dir_link = 2*mu + parity;
						#pragma acc loop seq
						for(int iter=0; iter<3; iter++){
							int perp_dirs[4][3] = { {1,2,3}, {0,2,3}, {0,1,3}, {0,1,2} };
							int nu = perp_dirs[mu][iter];
            
							#pragma acc cache (nnp_openacc[idxh:8])

							const int dir_mu_2R = 2*mu + !parity;
							const int dir_mu_2L = 2*mu + !parity;
							const int idx_pmu = nnp_openacc[idxh][mu][parity]; // r+mu
							#pragma acc cache (nnm_openacc[idx_pmu:8])

							const int dir_nu_1R = 2*nu + !parity;
							const int dir_nu_3R = 2*nu +  parity;
							const int dir_nu_1L = 2*nu +  parity;
							const int dir_nu_3L = 2*nu + !parity;
							const int idx_pnu = nnp_openacc[idxh][nu][parity]; // r+nu

							// computation of the Right part of the staple
							mat1_times_conj_mat2_times_conj_mat3_addto_mat4_absent_stag_phases_onlyferms(&u[dir_nu_1R], idx_pmu,
																																													 &u[dir_mu_2R], idx_pnu,
																																													 &u[dir_nu_3R], idxh,
																																													 &loc_stap[dir_link], idxh);

							const int idx_mnu = nnm_openacc[idxh][nu][parity]; // r-nu
							const int idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity]; // r+mu-nu

							// computation of the Left part of the staple
							conj_mat1_times_conj_mat2_times_mat3_addto_mat4_absent_stag_phases_onlyferms(&u[dir_nu_1L], idx_pmu_mnu,
																																													 &u[dir_mu_2L], idx_mnu,
																																													 &u[dir_nu_3L], idx_mnu,
																																													 &loc_stap[dir_link], idxh);
						} // iter
					} // mu
				} // d0
      } // d1
    } // d2
  } // d3
} // closes routine


#if NRANKS_D3 > 1
void calc_loc_staples_nnptrick_all_bulk(__restrict const su3_soa * const u,
																				__restrict su3_soa * const loc_stap )
{
  //       r+mu-nu  r+mu   r+mu+nu
  //          +<-----+----->+
  //          |  1L  ^  1R  |
  // mu    2L |      |      | 2R
  // ^        V  3L  |  3R  V
  // |        +----->+<-----+
  // |       r-nu    r     r+nu
  // +---> nu       
  // r is idxh in the following      

	#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
	#pragma acc loop independent gang(STAPGANG3) 
  for(int d3=D3_HALO+GAUGE_HALO; d3<nd3-D3_HALO-GAUGE_HALO; d3++) {
		#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(int d2=0; d2<nd2; d2++) {
      for(int d1=0; d1<nd1; d1++) {
				for(int d0=0; d0 < nd0; d0++) {
					const int idxh = snum_acc(d0,d1,d2,d3); // r
					const int parity = (d0+d1+d2+d3) % 2;
					#pragma acc loop seq 
					for(int mu=0; mu<4; mu++){
						const int dir_link = 2*mu + parity;
						#pragma acc loop seq
						for(int iter=0; iter<3; iter++){

							int nu;
							if (mu==0) { nu = iter + 1; }
							else if (mu==1) { nu = iter + (iter & 1) + (iter >> 1); }
							else if (mu==2) { nu = iter + (iter >> 1); }
							else if (mu==3) { nu = iter; }
							else { //error 
							}

							#pragma acc cache (nnp_openacc[idxh:8])

							const int dir_mu_2R = 2*mu + !parity;
							const int dir_mu_2L = 2*mu + !parity;
							const int idx_pmu = nnp_openacc[idxh][mu][parity]; // r+mu
							#pragma acc cache (nnm_openacc[idx_pmu:8])

							const int dir_nu_1R = 2*nu + !parity;
							const int dir_nu_3R = 2*nu +  parity;
							const int dir_nu_1L = 2*nu +  parity;
							const int dir_nu_3L = 2*nu + !parity;

							const int idx_pnu = nnp_openacc[idxh][nu][parity]; // r+nu

							// computation of the Right part of the staple

							mat1_times_conj_mat2_times_conj_mat3_addto_mat4_absent_stag_phases(&u[dir_nu_1R],       idx_pmu,
																																								 &u[dir_mu_2R],       idx_pnu,
																																								 &u[dir_nu_3R],       idxh,
																																								 &loc_stap[dir_link], idxh);
            
							const int idx_mnu = nnm_openacc[idxh][nu][parity] ;        // r-nu
							const int idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity]; // r+mu-nu

							// computation of the Left  part of the staple

							conj_mat1_times_conj_mat2_times_mat3_addto_mat4_absent_stag_phases(&u[dir_nu_1L],       idx_pmu_mnu,
																																								 &u[dir_mu_2L],       idx_mnu,
																																								 &u[dir_nu_3L],       idx_mnu,
																																								 &loc_stap[dir_link], idxh);
            
						} // iter
						// in the end, each staple has to be multiplied by the K_mu(x) of the link emanating the staples
						// N.B.: this factor is common both for the left and for the right staples, and common for all values of nu
						// double K_mu=u[dir_link].K.d[idxh];
						// gl3_times_double_factor(&(loc_stap[dir_link]), idxh, K_mu);
					} // mu
				} // d0
      } // d1
    } // d2
  } // d3
} // closes routine

void calc_loc_staples_nnptrick_all_d3c(__restrict const su3_soa * const u,
																			 __restrict su3_soa * const loc_stap,
																			 int offset, int thickness)
{
  //       r+mu-nu  r+mu   r+mu+nu
  //          +<-----+----->+
  //          |  1L  ^  1R  |
  // mu    2L |      |      | 2R
  // ^        V  3L  |  3R  V
  // |        +----->+<-----+
  // |       r-nu    r     r+nu
  // +---> nu       
  // r is idxh in the following      


	#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
	#pragma acc loop independent gang
  for(int d3=offset; d3<offset+thickness; d3++) {
		#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(int d2=0; d2<nd2; d2++) {
      for(int d1=0; d1<nd1; d1++) {
				for(int d0=0; d0 < nd0; d0++) {
					const int idxh = snum_acc(d0,d1,d2,d3); // r
					const int parity = (d0+d1+d2+d3) % 2;
					#pragma acc loop seq 
					for(int mu=0; mu<4; mu++){
						const int dir_link = 2*mu + parity;
						#pragma acc loop seq
						for(int iter=0; iter<3; iter++){

							int nu;
							if (mu==0) { nu = iter + 1; }
							else if (mu==1) { nu = iter + (iter & 1) + (iter >> 1); }
							else if (mu==2) { nu = iter + (iter >> 1); }
							else if (mu==3) { nu = iter; }
							else { // error 
							}

							#pragma acc cache (nnp_openacc[idxh:8])

							const int dir_mu_2R = 2*mu + !parity;
							const int dir_mu_2L = 2*mu + !parity;
							const int idx_pmu = nnp_openacc[idxh][mu][parity]; // r+mu
							#pragma acc cache (nnm_openacc[idx_pmu:8])

							const int dir_nu_1R = 2*nu + !parity;
							const int dir_nu_3R = 2*nu +  parity;
							const int dir_nu_1L = 2*nu +  parity;
							const int dir_nu_3L = 2*nu + !parity;

							const int idx_pnu = nnp_openacc[idxh][nu][parity]; // r+nu

							// computation of the Right part of the staple
	      
							mat1_times_conj_mat2_times_conj_mat3_addto_mat4_absent_stag_phases(&u[dir_nu_1R],       idx_pmu,
																																								 &u[dir_mu_2R],       idx_pnu,
																																								 &u[dir_nu_3R],       idxh,
																																								 &loc_stap[dir_link], idxh);
            
           

							const int idx_mnu = nnm_openacc[idxh][nu][parity] ;        // r-nu
							const int idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity]; // r+mu-nu

							// computation of the Left  part of the staple

							conj_mat1_times_conj_mat2_times_mat3_addto_mat4_absent_stag_phases(&u[dir_nu_1L],       idx_pmu_mnu,
																																								 &u[dir_mu_2L],       idx_mnu,
																																								 &u[dir_nu_3L],       idx_mnu,
																																								 &loc_stap[dir_link], idxh);
            
						} // iter
						// in the end, each staple has to be multiplied by the K_mu(x) of the link emanating the staples
						// N.B.: this factor is common both for the left and for the right staples, and common for all values of nu
						// double K_mu = u[dir_link].K.d[idxh];
						// gl3_times_double_factor(&(loc_stap[dir_link]), idxh, K_mu);
					} // mu
				} // d0
      } // d1
    } // d2
  } // d3
} // closes routine


#endif


void calc_loc_staples_nnptrick_all_only_even(__restrict const su3_soa * const u,
																						 __restrict su3_soa * const loc_stap )
{
  //       r+mu-nu  r+mu   r+mu+nu
  //          +<-----+----->+
  //          |  1L  ^  1R  |
  // mu    2L |      |      | 2R
  // ^        V  3L  |  3R  V
  // |        +----->+<-----+
  // |       r-nu    r     r+nu
  // +---> nu       
  // r is idxh in the following      


	#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
	#pragma acc loop independent gang(STAPGANG3) 
  for(int d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
		#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(int d2=0; d2<nd2; d2++) {
      for(int d1=0; d1<nd1; d1++) {
				for(int hd0=0; hd0 < nd0h; hd0++) {
					const int d0 = 2*hd0 + ((d1+d2+d3) & 0x1); // d0 = 2*hd0 + ((d1+d2+d3+1) & 0x1); (for the odd case)
					const int idxh = snum_acc(d0,d1,d2,d3); // r 
					const int parity = (d0+d1+d2+d3) % 2;
					#pragma acc loop seq 
					for(int mu=0; mu<4; mu++){
						const int dir_link = 2*mu + parity;
						#pragma acc loop seq
						for(int iter=0; iter<3; iter++){

							int nu;
							if (mu==0) { nu = iter + 1; }
							else if (mu==1) { nu = iter + (iter & 1) + (iter >> 1); }
							else if (mu==2) { nu = iter + (iter >> 1); }
							else if (mu==3) { nu = iter; }
							else { // error 
							}
							#pragma acc cache (nnp_openacc[idxh:8])

							const int dir_mu_2R = 2*mu + !parity;
							const int dir_mu_2L = 2*mu + !parity;
							const int idx_pmu = nnp_openacc[idxh][mu][parity]; // r+mu
							#pragma acc cache (nnm_openacc[idx_pmu:8])

							const int dir_nu_1R = 2*nu + !parity;
							const int dir_nu_3R = 2*nu +  parity;
							const int dir_nu_1L = 2*nu +  parity;
							const int dir_nu_3L = 2*nu + !parity;

							const int idx_pnu = nnp_openacc[idxh][nu][parity]; // r+nu

							// computation of the Right part of the staple

							mat1_times_conj_mat2_times_conj_mat3_addto_mat4_absent_stag_phases(&u[dir_nu_1R],       idx_pmu,
																																								 &u[dir_mu_2R],       idx_pnu,
																																								 &u[dir_nu_3R],       idxh,
																																								 &loc_stap[dir_link], idxh);

							const int idx_mnu = nnm_openacc[idxh][nu][parity] ;        // r-nu
							const int idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity]; // r+mu-nu

							// computation of the Left  part of the staple

							conj_mat1_times_conj_mat2_times_mat3_addto_mat4_absent_stag_phases(&u[dir_nu_1L],       idx_pmu_mnu,
																																								 &u[dir_mu_2L],       idx_mnu,
																																								 &u[dir_nu_3L],       idx_mnu,
																																								 &loc_stap[dir_link], idxh);

						}  // iter
						// in the end, each staple has to be multiplied by the K_mu(x) of the link emanating the staples
						// N.B.: this factor is common both for the left and for the right staples, and common for all values of nu
						// double K_mu=u[dir_link].K.d[idxh];
						// gl3_times_double_factor(&(loc_stap[dir_link]), idxh, K_mu);
					} // mu
				} // d0
      } // d1
    } // d2
  } // d3
} // closes routine

void calc_loc_staples_nnptrick_all_only_odd(__restrict const su3_soa * const u,
																						__restrict su3_soa * const loc_stap )
{
  //       r+mu-nu  r+mu   r+mu+nu
  //          +<-----+----->+
  //          |  1L  ^  1R  |
  // mu    2L |      |      | 2R
  // ^        V  3L  |  3R  V
  // |        +----->+<-----+
  // |       r-nu    r     r+nu
  // +---> nu       
  // r is idxh in the following      


	#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
	#pragma acc loop independent gang(STAPGANG3) 
  for(int d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
		#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(int d2=0; d2<nd2; d2++) {
			for(int d1=0; d1<nd1; d1++) {
				for(int hd0=0; hd0 < nd0h; hd0++) {
					const int d0 = 2*hd0 + ((d1+d2+d3+1) & 0x1);
					const int idxh = snum_acc(d0,d1,d2,d3); // r 
					const int parity = (d0+d1+d2+d3) % 2;
	  
					#pragma acc loop seq 
					for(int mu=0; mu<4; mu++){
						const int dir_link = 2*mu + parity;
						#pragma acc loop seq
						for(int iter=0; iter<3; iter++){

							int nu;
							if (mu==0) { nu = iter + 1; }
							else if (mu==1) { nu = iter + (iter & 1) + (iter >> 1); }
							else if (mu==2) { nu = iter + (iter >> 1); }
							else if (mu==3) { nu = iter; }
							else { //error 
							}
							#pragma acc cache (nnp_openacc[idxh:8])

							const int dir_mu_2R = 2*mu + !parity;
							const int dir_mu_2L = 2*mu + !parity;
							const int idx_pmu = nnp_openacc[idxh][mu][parity]; // r+mu
							#pragma acc cache (nnm_openacc[idx_pmu:8])

							const int dir_nu_1R = 2*nu + !parity;
							const int dir_nu_3R = 2*nu +  parity;
							const int dir_nu_1L = 2*nu +  parity;
							const int dir_nu_3L = 2*nu + !parity;

							const int idx_pnu = nnp_openacc[idxh][nu][parity]; // r+nu

							// computation of the Right part of the staple

							mat1_times_conj_mat2_times_conj_mat3_addto_mat4_absent_stag_phases(&u[dir_nu_1R],       idx_pmu,
																																								 &u[dir_mu_2R],       idx_pnu,
																																								 &u[dir_nu_3R],       idxh,
																																								 &loc_stap[dir_link], idxh);

							const int idx_mnu = nnm_openacc[idxh][nu][parity] ;         // r-nu
							const int idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu

							// computation of the Left  part of the staple

							conj_mat1_times_conj_mat2_times_mat3_addto_mat4_absent_stag_phases(&u[dir_nu_1L],       idx_pmu_mnu,
																																								 &u[dir_mu_2L],       idx_mnu,
																																								 &u[dir_nu_3L],       idx_mnu,
																																								 &loc_stap[dir_link], idxh);
            
						} // iter
						// in the end, each staple has to be multiplied by the K_mu(x) of the link emanating the staples
						// N.B.: this factor is common both for the left and for the right staples, and common for all values of nu
						// double K_mu=u[dir_link].K.d[idxh];
						// gl3_times_double_factor(&(loc_stap[dir_link]), idxh, K_mu);
					} // mu
				} // d0
      } // d1
    } // d2
  } // d3
} // closes routine

#endif
