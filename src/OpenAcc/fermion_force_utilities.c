#ifndef FERMION_FORCE_UTILITIES_C
#define FERMION_FORCE_UTILITIES_C

#include "./geometry.h"
#include "./fermion_force_utilities.h"
#include "./struct_c_def.h"
#include "../Include/fermion_parameters.h"
#include "./fermionic_utilities.h"
#include "./fermion_matrix.h"
#include "./backfield.h"

#ifdef MULTIDEVICE
#include "../Mpi/communications.h"
#endif

void set_tamat_soa_to_zero( __restrict tamat_soa * const matrix)
{
	int mu, idxh;
	#pragma acc kernels present(matrix)
	#pragma acc loop independent gang(8)
	for(mu=0; mu<8; mu++) {
		#pragma acc loop independent vector(32)
		for(idxh=0; idxh<sizeh; idxh++) {
			assign_zero_to_tamat_soa_component(&matrix[mu],idxh);

		}
	}
}

void direct_product_of_fermions_into_auxmat(
																						__restrict const vec3_soa  * const loc_s,
																						__restrict const vec3_soa  * const loc_h, 
																						__restrict su3_soa * const aux_u,
																						const RationalApprox * const approx,
																						int nnp_openacc[sizeh][4][2],
																						int iter)
{

	// loop on even sites
	int hd0, d1, d2, d3;
  #pragma acc kernels present(loc_s) present(loc_h)  present(approx) present(aux_u) present(nnp_openacc)
	#pragma acc loop independent gang (STAPGANG3)
	for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
		#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
		for(d2=0; d2<nd2; d2++) {
			for(d1=0; d1<nd1; d1++) {
				for(hd0=0; hd0 < nd0h; hd0++) {
					int idxh,idxpmu,d0;
					int parity;
					int dir_mu;
					int mu;
					d0 = 2*hd0 + ((d1+d2+d3) & 0x1);
					idxh = snum_acc(d0,d1,d2,d3); // r
					//  parity = (d0+d1+d2+d3) % 2;
					parity = 0; // selecting even parity site

					for(mu=0;mu<4;mu++){
						idxpmu = nnp_openacc[idxh][mu][parity]; // r+mu        
						dir_mu = 2*mu +  parity;
						vec1_directprod_conj_vec2_into_mat1(&aux_u[dir_mu],idxh,loc_h,idxpmu,loc_s,idxh,approx->RA_a[iter]);
					} // mu

				}  // d0     
			}  // d1       
		}  // d2         
	}  // d3

	// loop on odd sites
	#pragma acc kernels present(loc_s) present(loc_h) present(approx) present(aux_u) present(nnp_openacc)
	#pragma acc loop independent gang(STAPGANG3)
	for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
		#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
		for(d2=0; d2<nd2; d2++) {
			for(d1=0; d1<nd1; d1++) {
				for(hd0=0; hd0 < nd0h; hd0++) {
					int idxh,idxpmu,d0;
					int parity;
					int dir_mu;
					int mu;
					d0 = 2*hd0 + ((d1+d2+d3+1) & 0x1);
					idxh = snum_acc(d0,d1,d2,d3); // r
					// parity = (d0+d1+d2+d3) % 2;
					parity = 1; // selecting odd parity site
					#pragma acc loop independent
					for(mu=0;mu<4;mu++){
						idxpmu = nnp_openacc[idxh][mu][parity]; // r+mu        
						dir_mu = 2*mu +  parity;
						vec1_directprod_conj_vec2_into_mat1(&aux_u[dir_mu],idxh,loc_s,idxpmu,loc_h,idxh,-approx->RA_a[iter]);
					} //mu
				}  // d0     
			}  // d1       
		}  // d2         
	}  // d3
}

void multiply_conf_times_force_and_take_ta_nophase(
																									 __restrict const su3_soa * const u,
																									 __restrict const su3_soa * const auxmat, 
																									 __restrict tamat_soa * const ipdot)
{
	int hd0,d1,d2,d3,idxh,dir;
	#pragma acc kernels present(u) present(auxmat) present(ipdot) 
	#pragma acc loop independent
	for(dir = 0; dir< 8 ; dir++){
		#pragma acc loop independent gang(STAPGANG3)
		for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
			#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
			for(d2=0; d2<nd2; d2++) {
				for(d1=0; d1<nd1; d1++) {
					for(hd0=0; hd0 < nd0h; hd0++) {

						// even sites
						idxh = snum_acc(2*hd0,d1,d2,d3);

						mat1_times_auxmat_into_tamat_nophase(&u[dir],idxh,&auxmat[dir],idxh,&ipdot[dir],idxh); 

					} // hd0
				} // d1
			} // d2
		} // d3
	} // dir
}

void multiply_backfield_times_force(
																		__restrict ferm_param * const tpars,
																		__restrict const su3_soa * const auxmat, 
																		__restrict su3_soa * const pseudo_ipdot)
{

	int idxh,dirindex;
	int hd0, d1, d2, d3;
	double_soa * args = tpars->phases;
	#pragma acc kernels present(args) present(auxmat) present(pseudo_ipdot) 
	#pragma acc loop independent
	for(dirindex = 0 ; dirindex < 8 ; dirindex++){
		#pragma acc loop independent gang(STAPGANG3)
		for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
			#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
			for(d2=0; d2<nd2; d2++) {
				for(d1=0; d1<nd1; d1++) {
					for(hd0=0; hd0 < nd0h; hd0++) {

						idxh = snum_acc(2*hd0,d1,d2,d3);
						double arg;                 
						d_complex phase;            
						arg = args[dirindex].d[idxh]; 
						phase = cos(arg) + I * sin(arg);                      
						phase_times_auxmat_into_auxmat(&auxmat[dirindex],&pseudo_ipdot[dirindex],idxh,phase);
					}
				}
			} 
		}
	}
}

void accumulate_gl3soa_into_gl3soa(
																	 __restrict const su3_soa * const auxmat,
																	 __restrict su3_soa * const pseudo_ipdot)
{


	int idxh, dirindex;
	int hd0, d1, d2, d3;
	#pragma acc kernels present(auxmat) present(pseudo_ipdot)
	#pragma acc loop independent
	for(dirindex = 0 ; dirindex < 8 ; dirindex++){
		#pragma acc loop independent gang(STAPGANG3)
		for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
			#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
			for(d2=0; d2<nd2; d2++) {
				for(d1=0; d1<nd1; d1++) {
					for(hd0=0; hd0 < nd0h; hd0++) {

						idxh = snum_acc(2*hd0,d1,d2,d3);
						accumulate_auxmat1_into_auxmat2(&auxmat[dirindex],&pseudo_ipdot[dirindex],idxh);
					}
				}
			}
		}
	}
}


void ker_openacc_compute_fermion_force( 
																			 __restrict const su3_soa * const u,
																			 __restrict su3_soa * const aux_u,
																			 __restrict const vec3_soa * const in_shiftmulti,
																			 __restrict vec3_soa  * const loc_s,
																			 __restrict vec3_soa  * const loc_h,
																			 int nnp_openacc[sizeh][4][2],
																			 ferm_param  *  tpars)
{
	int ih;
	int iter=0;

	for(iter=0; iter<tpars->approx_md.approx_order; iter++){
		assign_in_to_out(&in_shiftmulti[iter],loc_s);
		acc_Doe(u,loc_h,loc_s,tpars->phases);
		direct_product_of_fermions_into_auxmat(loc_s,loc_h,aux_u,&(tpars->approx_md),nnp_openacc,iter);


	}
}


#endif
