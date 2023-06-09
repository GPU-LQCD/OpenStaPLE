#ifndef RANDOM_ASSIGNEMENT_C_
#define RANDOM_ASSIGNEMENT_C_

#include "../Include/debug.h"
#include "../Mpi/multidev.h"
#include "../Rand/random.h"
#include "./alloc_vars.h"
#include "./geometry.h"
#include "./random_assignement.h"
#include "./single_types.h"
#include "./su3_utilities.h"

#define acc_twopi 2*3.14159265358979323846

#ifdef __GNUC__
#include "math.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#endif

#ifdef MULTIDEVICE
#include "../Mpi/communications.h"
#endif



// gaussian random number generator
double double_gauss(void)
{
	double phi, temp, radius, ris;

	phi=acc_twopi*casuale();
	temp=-1.0*log(casuale());
	radius=sqrt(temp);

	ris=radius*cos(phi);

	return ris;
}


// gaussian random number generator
void two_double_gauss(double *r)
{
	double phi, temp, radius;

	phi=acc_twopi*casuale();
	temp=-1.0*log(casuale());
	radius=sqrt(temp);

	r[0]=radius*cos(phi);
	r[1]=radius*sin(phi);
}


// complex random gaussian number
d_complex d_complex_gauss(void)
{
	double num[2];

	two_double_gauss(num);

	return num[0] + I*num[1];
}

void generate_vec3_soa_gauss(__restrict vec3_soa * const vect)
{

	int d0h, d1, d2, d3;

	for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) for(d2=0; d2<nd2; d2++)
																					for(d1=0; d1<nd1; d1++)	for(d0h=0; d0h < nd0h; d0h++) {

																							// I take the size to be even, but it's the same
																							int  d0 = 2*d0h + ((d1+d2+d3) & 0x1);
																							int t  = snum_acc(d0,d1,d2,d3);  

																							// for fake rng, fakeness level == 1
																							int glt1D = t + LOC_SIZEH * devinfo.myrank - LNH_VOL3/2 * D3_HALO;
																							rng_fake_gl_index = glt1D * 3+1 ;


																							vect->c0[t]=d_complex_gauss(); rng_fake_gl_index++; 
																							vect->c1[t]=d_complex_gauss(); rng_fake_gl_index++;
																							vect->c2[t]=d_complex_gauss();
																						}

#if NRANKS_D3 > 1	//#ifdef MULTIDEVICE
	communicate_fermion_borders_hostonly(vect);
#endif
}

static inline double norm3(d_complex a,d_complex b, d_complex c)
{

	double norm ;
	norm  = creal(a)*creal(a) +  cimag(a)*cimag(a);
	norm += creal(b)*creal(b) +  cimag(b)*cimag(b);
	norm += creal(c)*creal(c) +  cimag(c)*cimag(c);
	return norm;

}

// generate random SU(3) matrix
void generate_random_su3(single_su3 * m, double f)
{
	// only first two rows
	m->comp[0][0] =1+(1.0-2.0*casuale()+I*(1.0-2.0*casuale()))*f;
	m->comp[0][1] =  (1.0-2.0*casuale()+I*(1.0-2.0*casuale()))*f;
	m->comp[0][2] =  (1.0-2.0*casuale()+I*(1.0-2.0*casuale()))*f;
	m->comp[1][0] =  (1.0-2.0*casuale()+I*(1.0-2.0*casuale()))*f;
	m->comp[1][1] =1+(1.0-2.0*casuale()+I*(1.0-2.0*casuale()))*f;
	m->comp[1][2] =  (1.0-2.0*casuale()+I*(1.0-2.0*casuale()))*f;

	// unitarize it
	double norm = norm3(m->comp[0][0],m->comp[0][1],m->comp[0][2]);
	norm=1.0/sqrt(norm);

	m->comp[0][0] *= norm;
	m->comp[0][1] *= norm;
	m->comp[0][2] *= norm;

	d_complex prod  = conj(m->comp[0][0]) * m->comp[1][0] + conj(m->comp[0][1]) * m->comp[1][1] + conj(m->comp[0][2]) * m->comp[1][2];

	m->comp[1][0] -= prod * m->comp[0][0];
	m->comp[1][1] -= prod * m->comp[0][1];
	m->comp[1][2] -= prod * m->comp[0][2];

	norm = norm3(m->comp[1][0],m->comp[1][1],m->comp[1][2]);
	norm=1.0/sqrt(norm);



	m->comp[1][0] *= norm;
	m->comp[1][1] *= norm;
	m->comp[1][2] *= norm;


	rebuild3row(m);

}

// iterations of random assignements over all the 8 components
void generate_Momenta_gauss(__restrict thmat_soa * const mom8)
{
	int mu;
	for(mu=0; mu<8; mu++){
		thmat_soa * mom = &mom8[mu];

		double casuali[8], aux[2];
		double uno_su_radice_di_tre = 1.0/sqrt(3.0);
		int d0h, d1, d2, d3;

		for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) for(d2=0; d2<nd2; d2++)
																						for(d1=0; d1<nd1; d1++)	for(d0h=0; d0h < nd0h; d0h++) {

																								// I take the size to be even, but it's the same
																								int  d0 = 2*d0h + ((d1+d2+d3) & 0x1);
																								int t  = snum_acc(d0,d1,d2,d3); 

																								// for fake rng, fakeness level == 1
																								int glt1D = t + LOC_SIZEH * devinfo.myrank - LNH_VOL3/2 * D3_HALO;
																								rng_fake_gl_index = glt1D * 8+1 ;

																								int i;
																								for(i=0; i<4; i++)
																									{
																										two_double_gauss(aux); rng_fake_gl_index++;
																										casuali[2*i]   = aux[0];
																										casuali[2*i+1] = aux[1];
																									}
																								mom->rc00[t] =  casuali[2] + casuali[7] * uno_su_radice_di_tre;
																								mom->rc11[t] = -casuali[2] + casuali[7] * uno_su_radice_di_tre;
																								mom->c01[t]  =  casuali[0] - casuali[1] * I;
																								mom->c02[t]  =  casuali[3] - casuali[4] * I;
																								mom->c12[t]  =  casuali[5] - casuali[6] * I;
																							}   
	}

#ifdef MULTIDEVICE
	// should not be necessary to communicate momenta, they are not propagated.
	// communicate_thmat_soa_borders(mom8, HALO_WIDTH); 
#endif
}

// iterations of random assignements over all the 8 components
void generate_Conf_cold(__restrict su3_soa * const conf,double factor)
{
	printf("\t Generating random cold configuration with eps = %f .\n", factor);

	int d0h, d1, d2, d3;
	for(int mu=0; mu<8; mu++)
		for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) for(d2=0; d2<nd2; d2++)
																						for(d1=0; d1<nd1; d1++)	for(d0h=0; d0h < nd0h; d0h++) {
                
																								// I take the size to be even, but it's the same
																								int  d0 = 2*d0h + ((d1+d2+d3) & 0x1);
																								int t  = snum_acc(d0,d1,d2,d3);

																								single_su3 aux;
																								generate_random_su3(&aux,factor);
																								single_gl3_into_su3_soa(&conf[mu],t,&aux);

																							}
#if NRANKS_D3 > 1	//#ifdef MULTIDEVICE
	communicate_su3_borders_hostonly(conf, HALO_WIDTH);
#endif

}




void generate_vec3_soa_z2noise(__restrict vec3_soa * const vect)
{
	double p;
	int d0h, d1, d2, d3;
	for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) for(d2=0; d2<nd2; d2++)
																					for(d1=0; d1<nd1; d1++)	for(d0h=0; d0h < nd0h; d0h++) {

																							// I take the size to be even, but it's the same
																							int  d0 = 2*d0h + ((d1+d2+d3) & 0x1);
																							int t  = snum_acc(d0,d1,d2,d3); 


																							p = casuale();
																							if(p<0.5){
																								vect->c0[t]=  1.0 + 0.0*I;
																							}else{
																								vect->c0[t]= -1.0 + 0.0*I;
																							}
																							p = casuale();
																							if(p<0.5){
																								vect->c1[t]=  1.0 + 0.0*I;
																							}else{
																								vect->c1[t]= -1.0 + 0.0*I;
																							}
																							p = casuale();
																							if(p<0.5){
																								vect->c2[t]=  1.0 + 0.0*I;
																							}else{
																								vect->c2[t]= -1.0 + 0.0*I;
																							}
																						}
#if NRANKS_D3 > 1	//#ifdef MULTIDEVICE
	communicate_fermion_borders_hostonly(vect);
#endif

}


#endif
