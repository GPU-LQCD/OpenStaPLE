#ifndef DBGTOOLS_C_
#define DBGTOOLS_C_

#include <complex.h>

#include "../OpenAcc/geometry.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include "../Include/common_defines.h"
#include "../OpenAcc/struct_c_def.h"
#include "../OpenAcc/single_types.h"
#include "../OpenAcc/io.h"
#include "../OpenAcc/alloc_vars.h"
#include "../Mpi/communications.h"

#define acc_pi 3.14159265358979323846

int ipdot_f_reset; // flag to be set to 1 when md starts, 
// and to be set to zero when first force calculation happens
// if it is 1 force difference is not calculated.
int ipdot_g_reset; // same



// multi rank data structure read/write functions - for debug/testing
// GL3
void ascii_write_single_su3_into_su3_soa(
        const int idxh_machine, const int parity, const int dirmachine, 
        const void* datastruct,
        int conf_machine_endianness_disagreement_UNUSED, FILE* fp){

    global_su3_soa * conf = (global_su3_soa *) datastruct;

    single_su3 m;          
    single_su3_from_global_su3_soa(&conf[2*dirmachine+parity],idxh_machine,&m);
    rebuild3row(&m);
    for(int r=0; r < 3 ; r++) for(int c=0; c < 3 ; c++) 
        fprintf(fp, "%.18lf %.18lf\n",creal(m.comp[r][c]),cimag(m.comp[r][c])); 

}
void save_gl_gl3(const global_su3_soa * conf, const char* nomefile)
{

    FILE *fp;
    fp = fopen(nomefile,"w");

    rw_iterate_on_global_sites_lx_xyzt_axis_ordering(
            ascii_write_single_su3_into_su3_soa,
            (void*)conf, fp, 0);
    fclose(fp);
}

void ascii_read_single_su3_into_su3_soa(
        const int idxh_machine, const int parity, const int dirmachine, 
        const void* datastruct,
        const int conf_machine_endianness_disagreement_UNUSED, FILE* fp){

    global_su3_soa * conf = (global_su3_soa *) datastruct;
    single_su3 m;          

    double re, im;
    for(int r=0; r < 3 ; r++)
        for(int c=0; c < 3 ; c++){
            CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);
            m.comp[r][c] = re + im * I;
        }

    single_su3_into_global_su3_soa(&conf[2*dirmachine+parity],idxh_machine,&m);
}
int read_gl_gl3(global_su3_soa * conf, const char* nomefile)
{
    FILE *fp;
    fp = fopen(nomefile,"r");
    if(!fp){
        printf("gl3 file %s not readable.\n",nomefile);
        return 1;
    }

    rw_iterate_on_global_sites_lx_xyzt_axis_ordering(
            ascii_write_single_su3_into_su3_soa,
            (void*)conf,fp, 0);

    fclose(fp);
    return 0;
}


//FERMIONS

void ascii_write_single_vec3_into_vec3_soa(
        const int i, const int parity, const int dirmachine, 
        const void* datastruct,
        const int conf_machine_endianness_disagreement_UNUSED, FILE* fp){

    global_vec3_soa * fermion = (global_vec3_soa*) datastruct;
    fprintf(fp, "%.18lf\t%.18lf\n",creal(fermion->c0[i]),cimag(fermion->c0[i]));
    fprintf(fp, "%.18lf\t%.18lf\n",creal(fermion->c1[i]),cimag(fermion->c1[i]));
    fprintf(fp, "%.18lf\t%.18lf\n",creal(fermion->c2[i]),cimag(fermion->c2[i]));

}
void save_gl_fermion(const global_vec3_soa * fermion, 
        const char* nomefile)
{

    FILE *fp;
    fp = fopen(nomefile,"w");

    rw_iterate_on_global_sites_lx_xyzt_axis_ordering(
            ascii_write_single_vec3_into_vec3_soa,
            (void*) fermion, fp, 1);

    fclose(fp);

}

void ascii_read_single_vec3_into_vec3_soa(
        const int i, const int parity, const int dirmachine, 
        const void* datastruct,
        const int conf_machine_endianness_disagreement_UNUSED, FILE* fp){

    global_vec3_soa * fermion = (global_vec3_soa*) datastruct;
    double re,im;

    CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);fermion->c0[i] = re + im * I;
    CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);fermion->c1[i] = re + im * I;
    CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);fermion->c2[i] = re + im * I;

}
int read_gl_fermion(global_vec3_soa * fermion, const char* nomefile)
{


    FILE *fp;
    fp = fopen(nomefile,"r");
    if(!fp){
        printf("vec3_soa File %s not found.\n", nomefile );
        return 1;
    }
    if(verbosity_lv > 2) 
        printf("Reading vec3_soa %s\n", nomefile );

    rw_iterate_on_global_sites_lx_xyzt_axis_ordering(
            ascii_read_single_vec3_into_vec3_soa,
            (void* ) fermion, fp, 1);
    
    fclose(fp);
    return 0;

}

// TAMAT
void ascii_write_single_tamat_into_tamat_soa(
        const int i, const int parity, const int dirmachine, 
        const void* datastruct,
        const int conf_machine_endianness_disagreement_UNUSED, FILE* fp){

    global_tamat_soa * tamat = (global_tamat_soa*) datastruct;
    tamat += 2*dirmachine + parity;
    fprintf(fp, "%.18lf\t%.18lf\n",creal(tamat->c01[i]),cimag(tamat->c01[i]));
    fprintf(fp, "%.18lf\t%.18lf\n",creal(tamat->c02[i]),cimag(tamat->c02[i]));
    fprintf(fp, "%.18lf\t%.18lf\n",creal(tamat->c12[i]),cimag(tamat->c12[i]));
    fprintf(fp, "0\t%.18lf\n",tamat->ic00[i]);
    fprintf(fp, "0\t%.18lf\n",tamat->ic11[i]);

}
void save_gl_tamat(const global_tamat_soa * tamat, 
        const char* nomefile)
{
    FILE *fp;
    fp = fopen(nomefile,"w");

    rw_iterate_on_global_sites_lx_xyzt_axis_ordering(
            ascii_write_single_tamat_into_tamat_soa,
            (void*) tamat, fp, 0);

    fclose(fp);

}

void ascii_read_single_tamat_into_tamat_soa(
        const int i, const int parity, const int dirmachine, 
        const void* datastruct,
        const int conf_machine_endianness_disagreement_UNUSED, FILE* fp){

    global_tamat_soa * tamat = (global_tamat_soa*) datastruct;
    tamat += 2*dirmachine + parity;
    double re,im;
    CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);tamat->c01[i] = re + im * I;
    CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);tamat->c02[i] = re + im * I;
    CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);tamat->c12[i] = re + im * I;
    CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);tamat->ic00[i] =im ;
    CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);tamat->ic11[i] =im ;


}
int read_gl_tamat(global_tamat_soa * tamat, const char* nomefile)
{

    FILE *fp;
    fp = fopen(nomefile,"r");
    if(!fp){
        printf("tamat_soa File %s not found.\n", nomefile );
        return 1;
    }
    if(verbosity_lv > 2) 
        printf("Reading tamat_soa %s\n", nomefile );

    rw_iterate_on_global_sites_lx_xyzt_axis_ordering(
            ascii_read_single_tamat_into_tamat_soa,
            (void*) tamat, fp, 0);

    fclose(fp);
    return 0;

}

// THMAT
void ascii_write_single_thmat_into_thmat_soa(
        const int i, const int parity, const int dirmachine, 
        const void* datastruct,
        const int conf_machine_endianness_disagreement_UNUSED, FILE* fp){

    global_thmat_soa * thmat = (global_thmat_soa*) datastruct;
    thmat += 2*dirmachine + parity;
    fprintf(fp, "%.18lf\t%.18lf\n",creal(thmat->c01[i]),cimag(thmat->c01[i]));
    fprintf(fp, "%.18lf\t%.18lf\n",creal(thmat->c02[i]),cimag(thmat->c02[i]));
    fprintf(fp, "%.18lf\t%.18lf\n",creal(thmat->c12[i]),cimag(thmat->c12[i]));
    fprintf(fp, "%.18lf\t0\n",thmat->rc00[i]);
    fprintf(fp, "%.18lf\t0\n",thmat->rc11[i]);

}
void save_gl_thmat(const global_thmat_soa * thmat, 
        const char* nomefile)
{
    FILE *fp;
    fp = fopen(nomefile,"w");

    rw_iterate_on_global_sites_lx_xyzt_axis_ordering(
            ascii_write_single_thmat_into_thmat_soa,
            (void*) thmat, fp, 0 );

    fclose(fp);

}

void ascii_read_single_thmat_into_thmat_soa(
        const int i, const int parity, const int dirmachine, 
        const void* datastruct,
        const int conf_machine_endianness_disagreement_UNUSED, FILE* fp){

    global_thmat_soa * thmat = (global_thmat_soa*) datastruct;
    double re,im;
    thmat += 2*dirmachine + parity;
    CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);thmat->c01[i] = re + im * I;
    CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);thmat->c02[i] = re + im * I;
    CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);thmat->c12[i] = re + im * I;
    CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);thmat->rc00[i] =re ;
    CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);thmat->rc11[i] =re ;

}
int read_gl_thmat(global_thmat_soa * thmat, const char* nomefile)
{

    FILE *fp;
    fp = fopen(nomefile,"r");
    if(!fp){
        printf("thmat_soa File %s not found.\n", nomefile );
        return 1;
    }
    if(verbosity_lv > 2) 
        printf("Reading thmat_soa %s\n", nomefile );

    rw_iterate_on_global_sites_lx_xyzt_axis_ordering(
            ascii_read_single_thmat_into_thmat_soa,
            (void*) thmat, fp, 0);

    fclose(fp);
    return 0;

}

// DCOMPLEX 
void ascii_write_dcomplex_into_dcomplex_soa(
        const int i, const int parity, const int dirmachine, 
        const void* datastruct,
        const int conf_machine_endianness_disagreement_UNUSED,FILE *fp){

    global_dcomplex_soa * arr = (global_dcomplex_soa *) datastruct;
    arr += 2*dirmachine + parity;
    fprintf(fp, "%.18lf\t%.18lf\n",creal(arr->c[i]),cimag(arr->c[i]));

}

void save_gl_dcomplex(const global_dcomplex_soa * arr, 
        const char* nomefile){
    FILE *fp;
    fp = fopen(nomefile,"w");

    rw_iterate_on_global_sites_lx_xyzt_axis_ordering(
            ascii_write_dcomplex_into_dcomplex_soa,
            (void *) arr, fp , 0); 

    fclose(fp);
}

void ascii_read_dcomplex_into_dcomplex_soa(
        const int i, const int parity, const int dirmachine, 
        const void* datastruct,
        const int conf_machine_endianness_disagreement_UNUSED,FILE *fp){

    global_dcomplex_soa * arr = (global_dcomplex_soa *) datastruct;
    arr += 2*dirmachine + parity;
    double re,im;
    CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);arr->c[i] = re + im * I;

}
int read_gl_dcomplex(global_dcomplex_soa * arr, const char* nomefile)
{

    FILE *fp;
    fp = fopen(nomefile,"r");
    if(!fp){
        printf("dcomplex_soa File %s not found.\n", nomefile );
        return 1;
    }
    if(verbosity_lv > 2) 
        printf("Reading dcomplex_soa %s\n", nomefile );

    rw_iterate_on_global_sites_lx_xyzt_axis_ordering(
            ascii_read_dcomplex_into_dcomplex_soa,
            (void*) arr, fp, 0); 

    fclose(fp);
    return 0;

}


// double
void ascii_write_double_into_double_soa(
        const int i, const int parity, const int dirmachine, 
        const void* datastruct,
        const int conf_machine_endianness_disagreement_UNUSED,FILE *fp){

    global_double_soa * arr = (global_double_soa *) datastruct;
    arr += 2*dirmachine + parity;
    fprintf(fp, "%.18lf\n",arr->d[i]);

}

void save_gl_double(const global_double_soa * arr, 
        const char* nomefile){
    FILE *fp;
    fp = fopen(nomefile,"w");

    rw_iterate_on_global_sites_lx_xyzt_axis_ordering(
            ascii_write_double_into_double_soa,
            (void *) arr, fp , 0); 

    fclose(fp);
}

void ascii_read_double_into_double_soa(
        const int i, const int parity, const int dirmachine, 
        const void* datastruct,
        const int conf_machine_endianness_disagreement_UNUSED,FILE *fp){

    global_double_soa * arr = (global_double_soa *) datastruct;
    arr += 2*dirmachine + parity;
    CHECKREAD(fscanf(fp, "%lf\n",&(arr->d[i])),1);

}
int read_gl_double(global_double_soa * arr, const char* nomefile)
{

    FILE *fp;
    fp = fopen(nomefile,"r");
    if(!fp){
        printf("double_soa File %s not found.\n", nomefile );
        return 1;
    }
    if(verbosity_lv > 2) 
        printf("Reading double_soa %s\n", nomefile );

    rw_iterate_on_global_sites_lx_xyzt_axis_ordering(
            ascii_read_double_into_double_soa,
            (void*) arr, fp, 0); 

    fclose(fp);
    return 0;

}

// multi rank wrappers

// gl3 / su3
void print_gl3_soa_wrapper(const su3_soa * gl3, 
        const char* nomefile)
{

    printf("MPI%02d - Saving whole gl3...\n", devinfo.myrank);

#ifdef MULTIDEVICE

    if(devinfo.myrank == 0){
        int irank;
        for(irank = 1 ; irank < devinfo.nranks; irank++)
            recv_loc_subconf_from_rank(conf_rw,irank,irank);
        recv_loc_subconf_from_buffer(conf_rw, gl3,0); 
        save_gl_gl3(conf_rw, nomefile);
    }
    else  send_lnh_subconf_to_master(gl3,devinfo.myrank);

#else 
    recv_loc_subconf_from_buffer(conf_rw,gl3,0);
    save_gl_gl3(conf_rw, nomefile);
#endif
    printf("MPI%02d - Saved whole gl3...\n", devinfo.myrank);
}

int read_gl3_soa_wrapper(su3_soa * gl3, const char* nomefile)
{

    int error =  0;
#ifdef MULTIDEVICE

    if(devinfo.myrank == 0){
        if(verbosity_lv > 2)
            printf("MPI%02d - reading global gl3\n",devinfo.myrank );
        error = read_gl_gl3(conf_rw, nomefile);
        MPI_Bcast((void*) &error,1,MPI_INT,0,devinfo.mpi_comm);

        if(!error) {
            send_lnh_subconf_to_buffer(conf_rw,gl3,0);
            int irank;
            for(irank = 1 ; irank < devinfo.nranks; irank++)
                send_lnh_subconf_to_rank(conf_rw,irank);

        }
        else 
            if(verbosity_lv > 2)
                printf("MPI%02d - no gl3 sent!\n",devinfo.myrank );

    }
    else{
        if(verbosity_lv > 2)
            printf("MPI%02d - receiving gl3\n",devinfo.myrank );

        MPI_Bcast((void*) &error,1,MPI_INT,0,devinfo.mpi_comm);
        if(!error){ 
            receive_lnh_subconf_from_master(gl3);
        }
        else 
            if(verbosity_lv > 2)
                printf("MPI%02d - no fermion received!\n",devinfo.myrank );
    }

#else 
    error = read_gl_gl3(conf_rw, nomefile);
    if(!error)  send_lnh_subconf_to_buffer(conf_rw,gl3,0);
#endif
    return error;
}


// fermions
void print_vec3_soa_wrapper(const vec3_soa * fermion, 
        const char* nomefile)
{

    printf("MPI%02d - Saving whole fermion...\n", devinfo.myrank);

#ifdef MULTIDEVICE

    if(devinfo.myrank == 0){
        int irank;
        for(irank = 1 ; irank < devinfo.nranks; irank++)
            recv_loc_subfermion_from_rank(ferm_rw,irank);
        recv_loc_subfermion_from_buffer(ferm_rw,fermion,0);
        save_gl_fermion(ferm_rw, nomefile);
    }
    else  send_lnh_subfermion_to_master(fermion,devinfo.myrank);

#else 
    recv_loc_subfermion_from_buffer(ferm_rw,fermion,0);
    save_gl_fermion(ferm_rw, nomefile);
#endif

    printf("MPI%02d - Saved whole fermion...\n", devinfo.myrank);

}

int read_vec3_soa_wrapper(vec3_soa * fermion, const char* nomefile)
{

    int error =  0;
#ifdef MULTIDEVICE

    if(devinfo.myrank == 0){
        if(verbosity_lv > 2)
            printf("MPI%02d - reading global fermion \n",devinfo.myrank );
        error = read_gl_fermion(ferm_rw, nomefile);
        MPI_Bcast((void*) &error,1,MPI_INT,0,devinfo.mpi_comm);

        if(!error) {
            send_lnh_subfermion_to_buffer(ferm_rw,fermion,0);
            int irank;
            for(irank = 1 ; irank < devinfo.nranks; irank++)
                send_lnh_subfermion_to_rank(ferm_rw,irank);

        }
        else 
            if(verbosity_lv > 2)
                printf("MPI%02d - no fermion sent!\n",devinfo.myrank );

    }
    else{
        if(verbosity_lv > 2)
            printf("MPI%02d - receiving fermion \n",devinfo.myrank );

        MPI_Bcast((void*) &error,1,MPI_INT,0,devinfo.mpi_comm);
        if(!error){ 
            receive_lnh_subfermion_from_master(fermion);
        }
        else 
            if(verbosity_lv > 2)
                printf("MPI%02d - no fermion received!\n",devinfo.myrank );
    }

#else 
    error = read_gl_fermion(ferm_rw, nomefile);
    if(!error)  send_lnh_subfermion_to_buffer(ferm_rw,fermion,0);
#endif

    return error;

}

// tamat
void print_tamat_soa_wrapper(const tamat_soa * tamat, 
        const char* nomefile)
{

    printf("MPI%02d - Saving whole tamat...\n", devinfo.myrank);

#ifdef MULTIDEVICE

    if(devinfo.myrank == 0){
        int irank;
        for(irank = 1 ; irank < devinfo.nranks; irank++)
            recv_loc_subtamat_from_rank(tamat_rw,irank);
        recv_loc_subtamat_from_buffer(tamat_rw, tamat,0); 
        save_gl_tamat(tamat_rw, nomefile);
    }
    else  send_lnh_subtamat_to_master(tamat,devinfo.myrank);

#else 
    recv_loc_subtamat_from_buffer(tamat_rw,tamat,0);
    save_gl_tamat(tamat_rw, nomefile);
#endif
    printf("MPI%02d - Saved whole tamat...\n", devinfo.myrank);
}

int read_tamat_soa_wrapper(tamat_soa * tamat, const char* nomefile)
{

    int error =  0;
#ifdef MULTIDEVICE

    if(devinfo.myrank == 0){
        if(verbosity_lv > 2)
            printf("MPI%02d - reading global tamat\n",devinfo.myrank );
        error = read_gl_tamat(tamat_rw, nomefile);
        MPI_Bcast((void*) &error,1,MPI_INT,0,devinfo.mpi_comm);

        if(!error) {
            send_lnh_subtamat_to_buffer(tamat_rw,tamat,0);
            int irank;
            for(irank = 1 ; irank < devinfo.nranks; irank++)
                send_lnh_subtamat_to_rank(tamat_rw,irank);

        }
        else 
            if(verbosity_lv > 2)
                printf("MPI%02d - no tamat sent!\n",devinfo.myrank );

    }
    else{
        if(verbosity_lv > 2)
            printf("MPI%02d - receiving tamat\n",devinfo.myrank );

        MPI_Bcast((void*) &error,1,MPI_INT,0,devinfo.mpi_comm);
        if(!error){ 
            receive_lnh_subtamat_from_master(tamat);
        }
        else 
            if(verbosity_lv > 2)
                printf("MPI%02d - no fermion received!\n",devinfo.myrank );
    }

#else 
    error = read_gl_tamat(tamat_rw, nomefile);
    if(!error)  send_lnh_subtamat_to_buffer(tamat_rw,tamat,0);
#endif
    return error;
}


// thmat
void print_thmat_soa_wrapper(const thmat_soa * thmat, 
        const char* nomefile)
{

    printf("MPI%02d - Saving whole thmat...\n", devinfo.myrank);

#ifdef MULTIDEVICE

    if(devinfo.myrank == 0){
        int irank;
        for(irank = 1 ; irank < devinfo.nranks; irank++)
            recv_loc_subthmat_from_rank(thmat_rw,irank);
        recv_loc_subthmat_from_buffer(thmat_rw, thmat,0); 
        save_gl_thmat(thmat_rw, nomefile);
    }
    else  send_lnh_subthmat_to_master(thmat,devinfo.myrank);

#else 
    recv_loc_subthmat_from_buffer(thmat_rw,thmat,0);
    save_gl_thmat(thmat_rw, nomefile);
#endif
    printf("MPI%02d - Saved whole thmat...\n", devinfo.myrank);
}

int read_thmat_soa_wrapper(thmat_soa * thmat, const char* nomefile)
{

    int error =  0;
#ifdef MULTIDEVICE

    if(devinfo.myrank == 0){
        if(verbosity_lv > 2)
            printf("MPI%02d - reading global thmat\n",devinfo.myrank );
        error = read_gl_thmat(thmat_rw, nomefile);
        MPI_Bcast((void*) &error,1,MPI_INT,0,devinfo.mpi_comm);

        if(!error) {
            send_lnh_subthmat_to_buffer(thmat_rw,thmat,0);
            int irank;
            for(irank = 1 ; irank < devinfo.nranks; irank++)
                send_lnh_subthmat_to_rank(thmat_rw,irank);

        }
        else 
            if(verbosity_lv > 2)
                printf("MPI%02d - no thmat sent!\n",devinfo.myrank );

    }
    else{
        if(verbosity_lv > 2)
            printf("MPI%02d - receiving thmat\n",devinfo.myrank );

        MPI_Bcast((void*) &error,1,MPI_INT,0,devinfo.mpi_comm);
        if(!error){ 
            receive_lnh_subthmat_from_master(thmat);
        }
        else 
            if(verbosity_lv > 2)
                printf("MPI%02d - no fermion received!\n",devinfo.myrank );
    }

#else 
    error = read_gl_thmat(thmat_rw, nomefile);
    if(!error)  send_lnh_subthmat_to_buffer(thmat_rw,thmat,0);
#endif
    return error;
}

// dcomplex
void print_dcomplex_soa_wrapper(const dcomplex_soa * dcarr, 
        const char* nomefile)
{

    printf("MPI%02d - Saving whole dcarr...\n", devinfo.myrank);

#ifdef MULTIDEVICE

    if(devinfo.myrank == 0){
        int irank;
        for(irank = 1 ; irank < devinfo.nranks; irank++)
            recv_loc_subdcomplex_from_rank(dcomplex_rw,irank);
        recv_loc_subdcomplex_from_buffer(dcomplex_rw, dcarr,0); 
        save_gl_dcomplex(dcomplex_rw, nomefile);
    }
    else  send_lnh_subdcomplex_to_master(dcarr,devinfo.myrank);

#else 
    recv_loc_subdcomplex_from_buffer(dcomplex_rw,dcarr,0);
    save_gl_dcomplex(dcomplex_rw, nomefile);
#endif
    printf("MPI%02d - Saved whole dcarr...\n", devinfo.myrank);
}

int read_dcomplex_soa_wrapper(dcomplex_soa * dcarr, const char* nomefile)
{

    int error =  0;
#ifdef MULTIDEVICE

    if(devinfo.myrank == 0){
        if(verbosity_lv > 2)
            printf("MPI%02d - reading global dcarr\n",devinfo.myrank );
        error = read_gl_dcomplex(dcomplex_rw, nomefile);
        MPI_Bcast((void*) &error,1,MPI_INT,0,devinfo.mpi_comm);

        if(!error) {
            send_lnh_subdcomplex_to_buffer(dcomplex_rw,dcarr,0);
            int irank;
            for(irank = 1 ; irank < devinfo.nranks; irank++)
                send_lnh_subdcomplex_to_rank(dcomplex_rw,irank);

        }
        else 
            if(verbosity_lv > 2)
                printf("MPI%02d - no dcarr sent!\n",devinfo.myrank );

    }
    else{
        if(verbosity_lv > 2)
            printf("MPI%02d - receiving dcarr\n",devinfo.myrank );

        MPI_Bcast((void*) &error,1,MPI_INT,0,devinfo.mpi_comm);
        if(!error){ 
            receive_lnh_subdcomplex_from_master(dcarr);
        }
        else 
            if(verbosity_lv > 2)
                printf("MPI%02d - no fermion received!\n",devinfo.myrank );
    }

#else 
    error = read_gl_dcomplex(dcomplex_rw, nomefile);
    if(!error)  send_lnh_subdcomplex_to_buffer(dcomplex_rw,dcarr,0);
#endif
    return error;
}


// double
void print_double_soa_wrapper(const double_soa * darr, 
        const char* nomefile)
{

    printf("MPI%02d - Saving whole darr...\n", devinfo.myrank);

#ifdef MULTIDEVICE

    if(devinfo.myrank == 0){
        int irank;
        for(irank = 1 ; irank < devinfo.nranks; irank++)
            recv_loc_subdouble_from_rank(double_rw,irank);
        recv_loc_subdouble_from_buffer(double_rw, darr,0); 
        save_gl_double(double_rw, nomefile);
    }
    else  send_lnh_subdouble_to_master(darr,devinfo.myrank);

#else 
    recv_loc_subdouble_from_buffer(double_rw,darr,0);
    save_gl_double(double_rw, nomefile);
#endif
    printf("MPI%02d - Saved whole darr...\n", devinfo.myrank);
}

int read_double_soa_wrapper(double_soa * darr, const char* nomefile)
{

    int error =  0;
#ifdef MULTIDEVICE

    if(devinfo.myrank == 0){
        if(verbosity_lv > 2)
            printf("MPI%02d - reading global darr\n",devinfo.myrank );
        error = read_gl_double(double_rw, nomefile);
        MPI_Bcast((void*) &error,1,MPI_INT,0,devinfo.mpi_comm);

        if(!error) {
            send_lnh_subdouble_to_buffer(double_rw,darr,0);
            int irank;
            for(irank = 1 ; irank < devinfo.nranks; irank++)
                send_lnh_subdouble_to_rank(double_rw,irank);

        }
        else 
            if(verbosity_lv > 2)
                printf("MPI%02d - no darr sent!\n",devinfo.myrank );

    }
    else{
        if(verbosity_lv > 2)
            printf("MPI%02d - receiving darr\n",devinfo.myrank );

        MPI_Bcast((void*) &error,1,MPI_INT,0,devinfo.mpi_comm);
        if(!error){ 
            receive_lnh_subdouble_from_master(darr);
        }
        else 
            if(verbosity_lv > 2)
                printf("MPI%02d - no fermion received!\n",devinfo.myrank );
    }

#else 
    error = read_gl_double(double_rw, nomefile);
    if(!error)  send_lnh_subdouble_to_buffer(double_rw,darr,0);
#endif
    return error;
}






void calc_loc_abelian_plaquettes(const double_soa* phases, // 8*
        double_soa * loc_abelian_plaquettes,// 2*
				int nnp_openacc[sizeh][4][2],
        const int mu, const int nu )
{
    int d0, d1, d2, d3;
    for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
        for(d2=0; d2<nd2; d2++) {
            for(d1=0; d1<nd1; d1++) {
                for(d0=0; d0 < nd0; d0++) {
                    int idxh,idxpmu,idxpnu;
                    int parity;
                    int dir_muA,dir_nuB;
                    int dir_muC,dir_nuD;

                    idxh = snum_acc(d0,d1,d2,d3);  // r 
                    parity = (d0+d1+d2+d3) % 2;

                    dir_muA = 2*mu +  parity;
                    dir_muC = 2*mu + (1-parity);
                    idxpmu = nnp_openacc[idxh][mu][parity];// r+mu

                    dir_nuB = 2*nu + (1-parity);
                    dir_nuD = 2*nu +  parity;
                    idxpnu = nnp_openacc[idxh][nu][parity];// r+nu
                    //       r+nu (C)  r+mu+nu
                    //          +<---+
                    // nu       |    ^
                    // ^    (D) V    | (B)
                    // |        +--->+
                    // |       r  (A)  r+mu
                    // +---> mu

                    double result =
                        phases[dir_muA].d[idxh] + phases[dir_nuB].d[idxpmu]-
                        phases[dir_muC].d[idxpnu]-phases[dir_nuD].d[idxh];   
                    while (result >  acc_pi ) 
                        result -=  2*acc_pi;
                    while (result <= - acc_pi ) 
                        result +=  2*acc_pi;

                    loc_abelian_plaquettes[parity].d[idxh] = result;

                }  // d0
            }  // d1
        }  // d2
    }  // d3

}


void calc_loc_abelian_plaquettes_device(const double_soa* phases, // 8*
        double_soa * loc_abelian_plaquettes,// 2*
			  int nnp_openacc[sizeh][4][2],
        const int mu, const int nu )
{
    int d0, d1, d2, d3;
#pragma acc kernels present(phases) present(loc_abelian_plaquettes)
#pragma acc loop independent gang 
    for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent gang vector //gang(nd2/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent gang vector //gang(nd1/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent vector //vector(DIM_BLOCK_X)
                for(d0=0; d0 < nd0; d0++) {
                    int idxh,idxpmu,idxpnu;
                    int parity;
                    int dir_muA,dir_nuB;
                    int dir_muC,dir_nuD;

                    idxh = snum_acc(d0,d1,d2,d3);  // r 
                    parity = (d0+d1+d2+d3) % 2;

                    dir_muA = 2*mu +  parity;
                    dir_muC = 2*mu + (1-parity);
                    idxpmu = nnp_openacc[idxh][mu][parity];// r+mu

                    dir_nuB = 2*nu + (1-parity);
                    dir_nuD = 2*nu +  parity;
                    idxpnu = nnp_openacc[idxh][nu][parity];// r+nu
                    //       r+nu (C)  r+mu+nu
                    //          +<---+
                    // nu       |    ^
                    // ^    (D) V    | (B)
                    // |        +--->+
                    // |       r  (A)  r+mu
                    // +---> mu

                    double result =
                        phases[dir_muA].d[idxh] + phases[dir_nuB].d[idxpmu]-
                        phases[dir_muC].d[idxpnu]-phases[dir_nuD].d[idxh];   
                    while (result >  acc_pi ) 
                        result -=  2*acc_pi;
                    while (result <= - acc_pi ) 
                        result +=  2*acc_pi;

                    loc_abelian_plaquettes[parity].d[idxh] = result;

                }  // d0
            }  // d1
        }  // d2
    }  // d3

}



void print_all_abelian_plaquettes(const double_soa* phases,const char * filename,int nnp_openacc[sizeh][4][2]){


    printf("Printing abelian plaquettes in file %s ...\n", filename);
    double_soa * all_plaquettes = malloc(2*6*sizeof(double_soa));
    int mu,nu,i,j;
    int pcount = 0;


    FILE *fp = fopen(filename, "w");
    fprintf(fp,"#");
    for(mu=0;mu<3;mu++) 
        for(nu=mu+1;nu<4;nu++){
					calc_loc_abelian_plaquettes(phases, &all_plaquettes[2*pcount], nnp_openacc,mu, nu);
            pcount++;
            fprintf(fp,"%d-%d\t",mu,nu);
        }
    fprintf(fp,"\n");

    //EVEN
    for(i=0;i<sizeh;i++){
        for(j=0;j<6;j++)
            fprintf(fp, "%f\t",all_plaquettes[2*j].d[i]);
        fprintf(fp,"\n");
    }
    // ODD
    for(i=0;i<sizeh;i++){
        for(j=0;j<6;j++)
            fprintf(fp, "%f\t",all_plaquettes[2*j+1].d[i]);
        fprintf(fp,"\n");
    }
    fclose(fp);

}


// one-rank read/write functions
void dbg_print_su3_soa(su3_soa * const conf, const char* nomefile,int conf_id_iter)
{

    FILE *fp;
    fp = fopen(nomefile,"w");


    int nx = geom_par.gnx;
    int ny = geom_par.gny;
    int nz = geom_par.gnz;
    int nt = geom_par.gnt;



#pragma acc update host(conf[0:8])
    for(int q = 0 ; q < 8 ; q++){
        for(int i = 0 ; i < sizeh ; i++){
            // rebuilding the 3rd row
            single_su3 aux;
            single_su3_from_su3_soa(&conf[q],i,&aux);
            rebuild3row(&aux);


            for(int r=0; r < 3 ; r++)
                for(int c=0; c < 3 ; c++) 
                    fprintf(fp, "%.18lf %.18lf\n",
                            creal(aux.comp[r][c]),cimag(aux.comp[r][c])); 
        }
    }
    fclose(fp);
}
int dbgread_su3_soa(su3_soa * conf, const char* nomefile,int * conf_id_iter )
{


    FILE *fp;
    fp = fopen(nomefile,"r");
    if(!fp){
        printf("Gauge configuration file %s not readable.\n",nomefile);
        *conf_id_iter = -1;
        return 1;
    }
    else{

        int nx = geom_par.gnx;
        int ny = geom_par.gny;
        int nz = geom_par.gnz;
        int nt = geom_par.gnt;


        int nxt,nyt,nzt,ntt;
        int minus_det_count = 0;
        CHECKREAD(fscanf(fp,"%d\t%d\t%d\t%d\t%d\n",&nxt,&nyt,&nzt,&ntt,conf_id_iter),5);
        if((nx!=nxt)||(ny!=nyt)||(nz!=nzt)||(nz!=nzt)){
            printf("Error, configuration dimensions not compatible with code.\n");
            exit(1);
        }

        for(int q = 0 ; q < 8 ; q++){
            for(int i = 0 ; i < sizeh ; i++){
                double re,im;
                single_su3 m;double det;
                //      fscanf(fp, "%.18lf\t%.18lf\n",&re,&im);
                CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[0][0]=conf[q].r0.c0[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[0][1]=conf[q].r0.c1[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[0][2]=conf[q].r0.c2[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[1][0]=conf[q].r1.c0[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[1][1]=conf[q].r1.c1[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[1][2]=conf[q].r1.c2[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[2][0]=conf[q].r2.c0[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[2][1]=conf[q].r2.c1[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[2][2]=conf[q].r2.c2[i] = re + im * I;
                det = detSu3(&m);
                if(fabs(1+det) < 0.005 ){ // DEBUG, the limit should be FAR stricter.
                    if(verbosity_lv > 5)  printf("Warning in read_su3_soa_ASCII(), Det M = -1.\n");
                    minus_det_count ++;
                    conf[q].r0.c0[i] = -conf[q].r0.c0[i];
                    conf[q].r0.c1[i] = -conf[q].r0.c1[i];
                    conf[q].r0.c2[i] = -conf[q].r0.c2[i];
                    conf[q].r1.c0[i] = -conf[q].r1.c0[i];
                    conf[q].r1.c1[i] = -conf[q].r1.c1[i];
                    conf[q].r1.c2[i] = -conf[q].r1.c2[i];
                    conf[q].r2.c0[i] = -conf[q].r2.c0[i];
                    conf[q].r2.c1[i] = -conf[q].r2.c1[i];
                    conf[q].r2.c2[i] = -conf[q].r2.c2[i];

                }else if (fabs(fabs(det)-1.0)>0.005 ){  // DEBUG, SAME HERE
                    printf("Non unitary matrix read, Det M =  %.18lf\n",det);
                    double pippo =  fabs(det)-1.0;
                    printf("Diff %f,    ", pippo);
                    pippo = fabs(pippo);
                    printf("Diff %f\n", pippo);
                }

            }
        }

        if(minus_det_count !=0 ) {
            printf("Conf read has some matrices with determinant = -1.\n");
            printf("This might be due to the configuration being saved\n");
            printf("is multiplied by the staggered phase.\n");
            printf("In such case, we would expect %d matrices to have determinant = -1\n",
                    nx*ny*nz*nt*3/2 );
            printf("The count of matrices with determinant = -1 is %d. \n", minus_det_count);





        }
        fclose(fp);

#pragma acc update device(conf[0:8])
        return 0;
    }
}

void dbgprint_gl3_soa(su3_soa * const conf, const char* nomefile,int conf_id_iter)
{

    FILE *fp;
    fp = fopen(nomefile,"w");


    int nx = geom_par.gnx;
    int ny = geom_par.gny;
    int nz = geom_par.gnz;
    int nt = geom_par.gnt;



#pragma acc update host(conf[0:8])

    for(int q = 0 ; q < 8 ; q++){
        for(int i = 0 ; i < sizeh ; i++){
            // rebuilding the 3rd row
            single_su3 aux;
            single_gl3_from_su3_soa(&conf[q],i,&aux);


            for(int r=0; r < 3 ; r++)
                for(int c=0; c < 3 ; c++) 
                    fprintf(fp, "%.18lf %.18lf\n",
                            creal(aux.comp[r][c]),cimag(aux.comp[r][c])); 
        }
    }
    fclose(fp);
}
int dbgread_gl3_soa(su3_soa * conf, const char* nomefile,int * conf_id_iter )
{


    FILE *fp;
    fp = fopen(nomefile,"r");
    if(!fp){
        printf("Gauge configuration file %s not readable.\n",nomefile);
        *conf_id_iter = -1;
        return 1;
    }
    else{

        int nx = geom_par.gnx;
        int ny = geom_par.gny;
        int nz = geom_par.gnz;
        int nt = geom_par.gnt;


        int nxt,nyt,nzt,ntt;
        int minus_det_count = 0;

        for(int q = 0 ; q < 8 ; q++){
            for(int i = 0 ; i < sizeh ; i++){
                double re,im;
                single_su3 m;double det;
                //      fscanf(fp, "%.18lf\t%.18lf\n",&re,&im);
                CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[0][0]=conf[q].r0.c0[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[0][1]=conf[q].r0.c1[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[0][2]=conf[q].r0.c2[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[1][0]=conf[q].r1.c0[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[1][1]=conf[q].r1.c1[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[1][2]=conf[q].r1.c2[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[2][0]=conf[q].r2.c0[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[2][1]=conf[q].r2.c1[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[2][2]=conf[q].r2.c2[i] = re + im * I;
                det = detSu3(&m);
                if(fabs(1+det) < 0.005 ){ // DEBUG, the limit should be FAR stricter.
                    if(verbosity_lv > 5)  printf("Warning in read_su3_soa_ASCII(), Det M = -1.\n");
                    minus_det_count ++;
                    conf[q].r0.c0[i] = -conf[q].r0.c0[i];
                    conf[q].r0.c1[i] = -conf[q].r0.c1[i];
                    conf[q].r0.c2[i] = -conf[q].r0.c2[i];
                    conf[q].r1.c0[i] = -conf[q].r1.c0[i];
                    conf[q].r1.c1[i] = -conf[q].r1.c1[i];
                    conf[q].r1.c2[i] = -conf[q].r1.c2[i];
                    conf[q].r2.c0[i] = -conf[q].r2.c0[i];
                    conf[q].r2.c1[i] = -conf[q].r2.c1[i];
                    conf[q].r2.c2[i] = -conf[q].r2.c2[i];

                }else if (fabs(fabs(det)-1.0)>0.005 ){  // DEBUG, SAME HERE
                    printf("Non unitary matrix read, Det M =  %.18lf\n",det);
                    double pippo =  fabs(det)-1.0;
                    printf("Diff %f,    ", pippo);
                    pippo = fabs(pippo);
                    printf("Diff %f\n", pippo);
                }

            }
        }

        if(minus_det_count !=0 ) {
            printf("Conf read has some matrices with determinant = -1.\n");
            printf("This might be due to the configuration being saved\n");
            printf("is multiplied by the staggered phase.\n");
            printf("In such case, we would expect %d matrices to have determinant = -1\n",
                    nx*ny*nz*nt*3/2 );
            printf("The count of matrices with determinant = -1 is %d. \n", minus_det_count);





        }
        fclose(fp);
#pragma acc update device(conf[0:8])
        return 0;
    }
}


void print_vec3_soa(vec3_soa * const fermion, const char* nomefile)
{

    printf("MPI%02d: Writing vec3_soa %s\n", devinfo.myrank, nomefile);
    FILE *fp;
    fp = fopen(nomefile,"w");
    for(int i = 0 ; i < sizeh ; i++){
        fprintf(fp, "%.18lf\t%.18lf\n",creal(fermion->c0[i]),cimag(fermion->c0[i]));
        fprintf(fp, "%.18lf\t%.18lf\n",creal(fermion->c1[i]),cimag(fermion->c1[i]));
        fprintf(fp, "%.18lf\t%.18lf\n",creal(fermion->c2[i]),cimag(fermion->c2[i]));
    }
    fclose(fp);
    printf("MPI%02d: Written vec3_soa %s\n", devinfo.myrank, nomefile);

}

int read_vec3_soa(vec3_soa * fermion, const char* nomefile)
{


    FILE *fp;
    fp = fopen(nomefile,"r");
    if(!fp){
        printf("vec3_soa File %s not found.\n", nomefile );
        return 1;
    }
    else{
        if(verbosity_lv > 2) 
            printf("Reading vec3_soa %s\n", nomefile );

        for(int i = 0 ; i < sizeh ; i++){
            double re,im;
            CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);fermion->c0[i] = re + im * I;
            CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);fermion->c1[i] = re + im * I;
            CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);fermion->c2[i] = re + im * I;
        }
        fclose(fp);
        return 0;
    }

}

void print_1su3_soa(su3_soa * const conf, const char* nomefile)
{

    FILE *fp;
    fp = fopen(nomefile,"w");
    for(int q = 0 ; q < 1 ; q++){ // q = 0 !!!!!
        for(int i = 0 ; i < sizeh ; i++){
            fprintf(fp, "%.18lf\t%.18lf\n",creal(conf[q].r0.c0[i]),cimag(conf[q].r0.c0[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(conf[q].r0.c1[i]),cimag(conf[q].r0.c1[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(conf[q].r0.c2[i]),cimag(conf[q].r0.c2[i]));//
            fprintf(fp, "%.18lf\t%.18lf\n",creal(conf[q].r1.c0[i]),cimag(conf[q].r1.c0[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(conf[q].r1.c1[i]),cimag(conf[q].r1.c1[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(conf[q].r1.c2[i]),cimag(conf[q].r1.c2[i]));//
            fprintf(fp, "%.18lf\t%.18lf\n",creal(conf[q].r2.c0[i]),cimag(conf[q].r2.c0[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(conf[q].r2.c1[i]),cimag(conf[q].r2.c1[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(conf[q].r2.c2[i]),cimag(conf[q].r2.c2[i]));
        }
    }
    fclose(fp);
}

void print_tamat_soa(tamat_soa * const ipdot, const char* nomefile)
{
    FILE *fp;
#pragma acc update host(ipdot[0:8])
    fp = fopen(nomefile,"w");
    for(int q = 0 ; q < 8 ; q++){
        for(int i = 0 ; i < sizeh ; i++){
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c01[i]),cimag(ipdot[q].c01[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c02[i]),cimag(ipdot[q].c02[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c12[i]),cimag(ipdot[q].c12[i]));//
            fprintf(fp, "0\t%.18lf\n",ipdot[q].ic00[i]);
            fprintf(fp, "0\t%.18lf\n",ipdot[q].ic11[i]);
        }
    }
    fclose(fp);
}

int read_tamat_soa(tamat_soa * ipdot, const char* nomefile)
{

    FILE *fp;
    fp = fopen(nomefile,"r");
    if(!fp){
        printf("Tamat_soa File %s not found.\n", nomefile );
        return 1;
    }
    else{
        for(int q = 0 ; q < 8 ; q++){
            for(int i = 0 ; i < sizeh ; i++){
                double re,im;
                CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);ipdot[q].c01[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);ipdot[q].c02[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);ipdot[q].c12[i] = re + im * I;//
                CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&(ipdot[q].ic00[i])),2);
                CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&(ipdot[q].ic11[i])),2);
            }
        }
        fclose(fp);
#pragma acc update device(ipdot[0:8])
        return 0;
    }
}


void print_thmat_soa(thmat_soa * const ipdot, const char* nomefile)
{
    FILE *fp;
#pragma acc update host(ipdot[0:8])
    fp = fopen(nomefile,"w");
    for(int q = 0 ; q < 8 ; q++){
        for(int i = 0 ; i < sizeh ; i++){
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c01[i]),cimag(ipdot[q].c01[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c02[i]),cimag(ipdot[q].c02[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c12[i]),cimag(ipdot[q].c12[i]));//
            fprintf(fp, "%.18lf\t0\n",ipdot[q].rc00[i]);
            fprintf(fp, "%.18lf\t0\n",ipdot[q].rc11[i]);
        }
    }
    fclose(fp);
}
void print_1thmat_soa(thmat_soa * const ipdot, const char* nomefile)
{
    FILE *fp;
    fp = fopen(nomefile,"w");
    for(int q = 0 ; q < 1 ; q++){// q = 1 !!!!
        for(int i = 0 ; i < sizeh ; i++){
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c01[i]),cimag(ipdot[q].c01[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c02[i]),cimag(ipdot[q].c02[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c12[i]),cimag(ipdot[q].c12[i]));//
            fprintf(fp, "%.18lf\t0\n",ipdot[q].rc00[i]);
            fprintf(fp, "%.18lf\t0\n",ipdot[q].rc11[i]);
        }
    }
    fclose(fp);
}

int read_thmat_soa(thmat_soa * ipdot, const char* nomefile)
{

    FILE *fp;
    fp = fopen(nomefile,"r");
    if(!fp){
        printf("Thmat_soa File %s not found.\n", nomefile );
        return 1;
    }
    else{
        for(int q = 0 ; q < 8 ; q++){
            for(int i = 0 ; i < sizeh ; i++){
                double re,im;
                CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);ipdot[q].c01[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);ipdot[q].c02[i] = re + im * I;
                CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);ipdot[q].c12[i] = re + im * I;//
                CHECKREAD(fscanf(fp, "%lf\t%lf\n",&(ipdot[q].rc00[i]),&im),2);
                CHECKREAD(fscanf(fp, "%lf\t%lf\n",&(ipdot[q].rc11[i]),&im),2);
            }
        }
        fclose(fp);
#pragma acc update device(ipdot[0:8])
        return 0;
    }
}



void print_double_soa(double_soa * const backfield, const char* nomefile)
{


    printf("Saving file %s\n", nomefile);
    FILE *fp;
    fp = fopen(nomefile,"w");
    for(int dir = 0 ; dir < 4 ; dir++)
        for(int d3 = 0 ; d3 < nd3; d3++)
            for(int d2 = 0 ; d2 < nd2; d2++)
                for(int d1 = 0 ; d1 < nd1; d1++)
                    for(int d0 = 0 ; d0 < nd0; d0++){
                        int parity = (d0+d1+d2+d3)%2;

                        int idxh = snum_acc(d0,d1,d2,d3);

                        fprintf(fp,
                                "dir: %d (d0,d1,d2,d3): %d %d %d %d value: %.18lf\n",
                                dir, d0,d1,d2,d3,
                                backfield[geom_par.xyztmap[dir]+parity].d[idxh]);
                    }
    fclose(fp);
}

void print_1double_soa(double_soa * const vettore, const char* nomefile)
{

    FILE *fp;
    fp = fopen(nomefile,"w");
    for(int i = 0 ; i < sizeh ; i++){
        fprintf(fp, "%.18lf\n",vettore->d[i]);
    }
    fclose(fp);
}


void transpose_vec3_soa(vec3_soa * vecin, vec3_soa *vecout, int xmap, int ymap, int zmap, int tmap){

    // WARNING : NOT TESTED
    int x,y,z,t;
    int d[4], idxh, newidxh;

    for(d[3]=0; d[3] < nd3; d[3]++) for(d[2]=0; d[2] < nd2; d[2]++)
        for(d[1]=0; d[1] < nd1; d[1]++) for(d[0]=0; d[0] < nd0; d[0]++){

            idxh = snum_acc(d[0],d[1],d[2],d[3]);
            int tnd[4] = {nd0,nd1,nd2,nd3};
            x = d[xmap];int tnx = tnd[xmap] ;
            y = d[ymap];int tny = tnd[ymap] ;
            z = d[zmap];int tnz = tnd[zmap] ;
            t = d[tmap];int tnt = tnd[tmap] ;
            newidxh = (x+tnx*(y+tny*(z+tnz*t)))/2;

            vecout->c0[newidxh] = vecin->c0[idxh];
            vecout->c1[newidxh] = vecin->c1[idxh];
            vecout->c2[newidxh] = vecin->c2[idxh];

        }

}

void transpose_su3_soa(su3_soa * matin, su3_soa *matout, int xmap, int ymap, int zmap, int tmap){

    // WARNING : NOT TESTED
    int dir,parity;
    for(dir = 1; dir < 4 ; dir++) for(parity = 1; parity < 2 ; parity++)
    {
        int dirmod = geom_par.xyztmap[dir];
        transpose_vec3_soa(&(matin[2*dirmod+parity].r0),
                &(matout[2*dir+parity].r0),xmap,ymap,zmap,tmap);
        transpose_vec3_soa(&(matin[2*dirmod+parity].r1),
                &(matout[2*dir+parity].r1),xmap,ymap,zmap,tmap);
        transpose_vec3_soa(&(matin[2*dirmod+parity].r2),
                &(matout[2*dir+parity].r2),xmap,ymap,zmap,tmap);

    }

}


#endif
