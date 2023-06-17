#ifndef COMMUNICATIONS_H_
#define COMMUNICATIONS_H_

#include "./geometry_multidev.h"
#include "./geometry_multidev.h"
#include "../OpenAcc/struct_c_def.h"


#if NRANKS_D3 > 1
#include <mpi.h>
// fermion border communications (only FERMION HALO thick)
void communicate_fermion_borders(vec3_soa *lnh_fermion);
void communicate_fermion_borders_hostonly(vec3_soa *lnh_fermion);

#if defined(USE_MPI_CUDA_AWARE) || defined(__GNUC__)
void communicate_fermion_borders_async(vec3_soa *lnh_fermion, 
        MPI_Request* send_border_requests,// 6 element long
        MPI_Request* recv_border_requests);// 6 element long

void communicate_su3_borders_async(su3_soa* lnh_conf, int thickness,
        MPI_Request* send_border_requests, // 96 element long
        MPI_Request* recv_border_requests);// 96 element long

#endif

// gauge conf border communication (only GAUGE_HALO thick)
void communicate_su3_borders(su3_soa* lnh_conf, int thickness);

// gauge conf border communication (only GAUGE_HALO thick)
void communicate_su3_borders_hostonly(su3_soa* lnh_conf, int thickness);


// communications of gl3 quantities (e.g. Sigma for stouting)
void communicate_gl3_borders(su3_soa* lnh_conf, int thickness);

// used just at the beginning of MD trajectory, but also in stouting
void communicate_thmat_soa_borders(thmat_soa* lnh_momenta,int thickness);
// force communication
// this is performance critical and an async version should be produced
void communicate_tamat_soa_borders(tamat_soa* lnh_ipdot, int thickness);

// communication of lattice chunks, for file IO only on master
// chunks, conf
void send_lnh_subconf_to_rank(const global_su3_soa *gl_soa_conf, 
        int target_rank);
void recv_loc_subconf_from_rank(global_su3_soa *gl_soa_conf, 
        int target_rank, int tag);
void send_lnh_subconf_to_master(const su3_soa *lnh_soa_conf, int tag);
void receive_lnh_subconf_from_master(su3_soa* lnh_su3_conf);

// chunks, fermions
void send_lnh_subfermion_to_rank(const global_vec3_soa *gl_soa_ferm, 
        int target_rank);
void recv_loc_subfermion_from_rank(global_vec3_soa *gl_soa_ferm,
        int target_rank);
void send_lnh_subfermion_to_master(const vec3_soa *lnh_ferm, int tag);
void receive_lnh_subfermion_from_master(vec3_soa* lnh_ferm);

// chunks, tamats
void send_lnh_subtamat_to_rank(const global_tamat_soa *gl_soa_tamat, 
        int target_rank);
void recv_loc_subtamat_from_rank(global_tamat_soa *gl_soa_tamat,
        int target_rank);
void send_lnh_subtamat_to_master(const tamat_soa *lnh_tamat, int tag);
void receive_lnh_subtamat_from_master(tamat_soa* lnh_tamat);

// chunks, thmats
void send_lnh_subthmat_to_rank(const global_thmat_soa *gl_soa_thmat, 
        int target_rank);
void recv_loc_subthmat_from_rank(global_thmat_soa *gl_soa_thmat,
        int target_rank);
void send_lnh_subthmat_to_master(const thmat_soa *lnh_thmat, int tag);
void receive_lnh_subthmat_from_master(thmat_soa* lnh_thmat);

// chunks, dcomplex
void send_lnh_subdcomplex_to_rank(const global_dcomplex_soa *gl_soa_dcomplex, 
        int target_rank);
void recv_loc_subdcomplex_from_rank(global_dcomplex_soa *gl_soa_dcomplex,
        int target_rank);
void send_lnh_subdcomplex_to_master(const dcomplex_soa *lnh_dcomplex, int tag);
void receive_lnh_subdcomplex_from_master(dcomplex_soa* lnh_dcomplex);

// chunks, dcomplex
void send_lnh_subdouble_to_rank(const global_double_soa *gl_soa_double, 
        int target_rank);
void recv_loc_subdouble_from_rank(global_double_soa *gl_soa_double,
        int target_rank);
void send_lnh_subdouble_to_master(const double_soa *lnh_double, int tag);
void receive_lnh_subdouble_from_master(double_soa* lnh_double);




#endif



// internal communication on master rank
// chunks, conf
void send_lnh_subconf_to_buffer(const global_su3_soa *gl_soa_conf, 
        su3_soa *lnh_conf, int target_rank);
void recv_loc_subconf_from_buffer(global_su3_soa *gl_soa_conf,
        const su3_soa* lnh_conf, int target_rank);

// chunks, fermions
void send_lnh_subfermion_to_buffer(const global_vec3_soa *gl_soa_ferm,
        vec3_soa *lnh_ferm, int target_rank);
void recv_loc_subfermion_from_buffer(global_vec3_soa *gl_soa_ferm,
        const vec3_soa* lnh_ferm, int target_rank);

// chunks, tamat
void send_lnh_subtamat_to_buffer(const global_tamat_soa *gl_soa_tamat,
        tamat_soa *lnh_tamat, int target_rank);
void recv_loc_subtamat_from_buffer(global_tamat_soa *gl_soa_tamat,
        const tamat_soa* lnh_tamat, int target_rank);

// chunks, thmat
void send_lnh_subthmat_to_buffer(const global_thmat_soa *gl_soa_thmat,
        thmat_soa *lnh_thmat, int target_rank);
void recv_loc_subthmat_from_buffer(global_thmat_soa *gl_soa_thmat,
        const thmat_soa* lnh_thmat, int target_rank);

// chunks, dcomplex
void send_lnh_subdcomplex_to_buffer(const global_dcomplex_soa *gl_soa_dcomplex,
        dcomplex_soa *lnh_dcomplex, int target_rank);
void recv_loc_subdcomplex_from_buffer(global_dcomplex_soa *gl_soa_dcomplex,
        const dcomplex_soa* lnh_dcomplex, int target_rank);

// chunks, double
void send_lnh_subdouble_to_buffer(const global_double_soa *gl_soa_double,
        double_soa *lnh_double, int target_rank);
void recv_loc_subdouble_from_buffer(global_double_soa *gl_soa_double,
        const double_soa* lnh_double, int target_rank);





#endif
