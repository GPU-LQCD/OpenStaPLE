#ifndef IO_C_
#define IO_C_

// FOR ILDG IO
#define __USE_XOPEN2K
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <sys/types.h>
#define MAXILDGHEADERS 10


#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include "./io.h"
#include "./geometry.h"
#include "./struct_c_def.h"
#include "../Include/common_defines.h"
#include "./single_types.h"
#include "../Include/setting_file_parser.h"
#include "./action.h"
#include "./alloc_settings.h"
#include "../Include/fermion_parameters.h"

#define CHECKREAD_V2(expr,should_read) \
{int read = expr ;if(read != should_read) { \
    printf("%s:%d, Error, not read expected number of entries : %d read vs %d should_read\n.", __FILE__, __LINE__ , read,should_read);\
    return 2;}}\




int print_su3_soa_ASCII(global_su3_soa * const conf_rw, const char* nomefile,
        int conf_id_iter, double_soa *back_phases){

    FILE *fp;
    fp = fopen(nomefile,"w");
    if(! fp ){
        printf("ERROR, %s unreadable.\n",nomefile);
        return 2;
    }



    if (back_phases != NULL){
        printf("This function works correctly if the field passed as argument contains\n");
        printf("ONLY STAGGERED PHASES\n");
    }
    int nx = geom_par.gnx;
    int ny = geom_par.gny;
    int nz = geom_par.gnz;
    int nt = geom_par.gnt;


    double fm = 0;
    // if we're in pure gauge we can have problems
    if (alloc_info.NDiffFlavs != 0){ 
        fm = fermions_parameters[0].ferm_mass;
    }

    fprintf(fp,"%d %d %d %d %lf %lf %d %d\n",nx,ny,nz,nt, 
            act_params.beta,fm,alloc_info.NDiffFlavs,
            conf_id_iter);
    int expected_stag_phases_minuses_count = nx*ny*nz*nt*3/2;
    int stag_phases_minuses_count = 0;

    for(int q = 0 ; q < 8 ; q++){
        for(int i = 0 ; i < GL_SIZEH ; i++){
            // rebuilding the 3rd row
            single_su3 aux;
            single_su3_from_global_su3_soa(&conf_rw[q],i,&aux);
            rebuild3row(&aux);


            int stag_phase = 1 ;
            if (back_phases != NULL) if(fabs( back_phases[q].d[i]) > 0.001 ) {
                stag_phase = -1;
                stag_phases_minuses_count ++;
            }


            for(int r=0; r < 3 ; r++){
                for(int c=0; c < 3 ; c++) 
                    fprintf(fp, "(%.18lf, %.18lf)  ",
                            stag_phase*creal(aux.comp[r][c]),
                            stag_phase*cimag(aux.comp[r][c])); 

                fprintf(fp, "\n");

            }
            fprintf(fp, "\n");
        }
    }
    if (back_phases != NULL){
        printf("Written %d matrices with -1 determinant, of expected %d.\n",
                stag_phases_minuses_count, expected_stag_phases_minuses_count );

        if(stag_phases_minuses_count != expected_stag_phases_minuses_count) {
            printf("WRONG NUMBER OF DET=-1 MATRICES.\n");
            return 1;
        }

    }
    fclose(fp);
    return 0;
}


int read_su3_soa_ASCII(global_su3_soa * conf, const char* nomefile,int * conf_id_iter ){


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

        int boh;
        double tmpbeta, tmpmass;

        int minus_det_count = 0;
        CHECKREAD_V2(fscanf(fp,"%d %d %d %d %lf %lf %d %d\n",
                    &nxt,&nyt,&nzt,&ntt,
                    &tmpbeta,&tmpmass,&boh,conf_id_iter),8);
        if( nx!=nxt || ny!=nyt || nz!=nzt || (nt != ntt)){
            printf("Error, configuration dimensions not compatible with code.\n");
            printf("Code dimensions: d0 %d,d1 %d,d2 %d, d3 %d\n",nx,ny,nz,nt);
            printf("Conf dimensions: d0 %d,d1 %d,d2 %d, d3 %d\n",nxt,nyt,nzt,ntt);
            return 2;
        }else

            for(int q = 0 ; q < 8 ; q++){
                for(int i = 0 ; i < GL_SIZEH ; i++){
                    double re,im;
                    single_su3 m;double det;
                    //      fscanf(fp, "%.18lf\t%.18lf\n",&re,&im);
                    CHECKREAD_V2(fscanf(fp, "(%lf, %lf) ",&re,&im),2);m.comp[0][0]=conf[q].r0.c0[i] = re + im * I;
                    CHECKREAD_V2(fscanf(fp, "(%lf, %lf) ",&re,&im),2);m.comp[0][1]=conf[q].r0.c1[i] = re + im * I;
                    CHECKREAD_V2(fscanf(fp, "(%lf, %lf)\n",&re,&im),2);m.comp[0][2]=conf[q].r0.c2[i] = re + im * I;
                    CHECKREAD_V2(fscanf(fp, "(%lf, %lf) ",&re,&im),2);m.comp[1][0]=conf[q].r1.c0[i] = re + im * I;
                    CHECKREAD_V2(fscanf(fp, "(%lf, %lf) ",&re,&im),2);m.comp[1][1]=conf[q].r1.c1[i] = re + im * I;
                    CHECKREAD_V2(fscanf(fp, "(%lf, %lf)\n",&re,&im),2);m.comp[1][2]=conf[q].r1.c2[i] = re + im * I;
                    CHECKREAD_V2(fscanf(fp, "(%lf, %lf) ",&re,&im),2);m.comp[2][0]=conf[q].r2.c0[i] = re + im * I;
                    CHECKREAD_V2(fscanf(fp, "(%lf, %lf) ",&re,&im),2);m.comp[2][1]=conf[q].r2.c1[i] = re + im * I;
                    CHECKREAD_V2(fscanf(fp, "(%lf, %lf)\n",&re,&im),2);m.comp[2][2]=conf[q].r2.c2[i] = re + im * I;
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
        return 0;
    }
}


// STUFF FOR READING AND WRITING ILDG FILES

// ENDIANNESS SWAP

//check the endianness of the machine
// ILDG data is written BIG endian, if the machine is little endian -> disagreement true

int machine_is_little_endian()
{
    int little_endian=1;
    return (int) (*(char*)(&little_endian));
}

void doubleswe(double  * doubletoswap){

    char *p= (char*) doubletoswap; char t[4];
    t[0] = p[7]; p[7] = p[0] ; p[0] = t[0];
    t[1] = p[6]; p[6] = p[1] ; p[1] = t[1];
    t[2] = p[5]; p[5] = p[2] ; p[2] = t[2];
    t[3] = p[4]; p[4] = p[3] ; p[3] = t[3];

}
void uint64swe(uint64_t  * inttoswap){

    char *p= (char*) inttoswap; char t[4];
    t[0] = p[7]; p[7] = p[0] ; p[0] = t[0];
    t[1] = p[6]; p[6] = p[1] ; p[1] = t[1];
    t[2] = p[5]; p[5] = p[2] ; p[2] = t[2];
    t[3] = p[4]; p[4] = p[3] ; p[3] = t[3];

}
void uint32swe(uint32_t  * inttoswap){

    char *p= (char*) inttoswap; char t[2];
    t[0] = p[3]; p[3] = p[0] ; p[0] = t[0];
    t[1] = p[2]; p[2] = p[1] ; p[1] = t[1];

}
void uint16swe(uint16_t  * inttoswap){

    char *p= (char*) inttoswap; char t;
    t = p[1]; p[1] = p[0] ; p[0] = t;

}
typedef struct ILDG_header_t
{
    uint32_t magic_no; // LIME magic number 0x456789ab
    uint16_t version;
    uint16_t mbme_flag;
    uint64_t data_length;
    char type[128];
}ILDG_header;

int read_su3_soa_ildg_binary(
        global_su3_soa * conf,
        const char* nomefile,int * conf_id_iter )
{

    if(geom_par.initialized_check == 0){
        printf("ERROR: GEOMETRY MUST BE INITIALIZED!\n");
        printf("ERROR: in %s at line %d\n",__FILE__, __LINE__);
        return 2;
    } 


    FILE *fg;int reads;
    char header[1000];
    *conf_id_iter = 1000 ; // random number

    int nx = geom_par.gnx;
    int ny = geom_par.gny;
    int nz = geom_par.gnz;
    int nt = geom_par.gnt;

    ILDG_header ildg_headers[MAXILDGHEADERS]; // all ildg headers
    // note _FILE_OFFSET_BITS 64 
    off_t ildg_header_ends_positions[MAXILDGHEADERS]; // point in files 
    // where the message starts

    // where, in ildg_headers, is the ildg_format header
    int ildg_format_index =-1;   
    // where, in ildg_headers, is the ildg_binary_data header
    int ildg_binary_data_index =-1;
    // optionals
    int MD_traj_index = -1;
    int input_file_index = -1;// not used now 


    fg = fopen(nomefile,"r");
    if(!fg){
        printf("Gauge configuration ILDG file %s not readable.\n",nomefile);
        *conf_id_iter = -1;
        return 1;
    } 

    int alldimfound = 0;
    fseeko(fg,0,SEEK_SET);
    // Find all headers
    int i=0;
    int  conf_machine_endianness_disagreement = machine_is_little_endian();
    while((ildg_format_index==-1 ||  ildg_binary_data_index==-1) && fg){
        reads = fread(&ildg_headers[i],sizeof(ILDG_header),1,fg);
        if(reads!= 1 ){
            printf("Error in reading file: %s ,%d\n",__FILE__,__LINE__);
            return 2;
        }

        ildg_header_ends_positions[i] = ftello(fg);
        if(strcmp(ildg_headers[i].type,"ildg-format")==0) ildg_format_index = i;
        if(strcmp(ildg_headers[i].type,"ildg-binary-data")==0) ildg_binary_data_index = i;
        if(strcmp(ildg_headers[i].type,"MD_traj")==0) MD_traj_index = i;
        if(strcmp(ildg_headers[i].type,"input-file")==0) input_file_index = i;

        if(conf_machine_endianness_disagreement){
            uint32swe(&ildg_headers[i].magic_no);
            uint16swe(&ildg_headers[i].version);
            uint16swe(&ildg_headers[i].mbme_flag);
            uint64swe(&ildg_headers[i].data_length);
        }

        if(verbosity_lv>3){
            printf("header type: %s, data length:%"PRIu64", magic number: %"PRIu32"\n",
                    ildg_headers[i].type, ildg_headers[i].data_length,
                    ildg_headers[i].magic_no);
            printf("version: %d, mbme_flag:%d\n",
                    ildg_headers[i].version, ildg_headers[i].mbme_flag);

        }


        // padding to 8 bytes
        off_t missing_bytes = 
            (off_t)(ildg_headers[i].data_length%8==0?0:8-ildg_headers[i].data_length%8);
        // note _FILE_OFFSET_BITS 64 
        off_t mod_data_length =  ildg_headers[i].data_length + missing_bytes ;
        //        printf("%"PRIu64" vs %"PRIu64"\n", (uint64_t) mod_data_length, ildg_headers[i].data_length);

        // going to the next record  
        fseeko(fg,mod_data_length,SEEK_CUR);

        i++;
    }


    // read ildg-format  
    char ildg_format_str[1000];
    if(verbosity_lv>2)
        printf("Reading ildg-format...\n");


    fseeko(fg,ildg_header_ends_positions[ildg_format_index],SEEK_SET);
    reads = fread(ildg_format_str,1,ildg_headers[ildg_format_index].data_length,fg);
    if(reads!= ildg_headers[ildg_format_index].data_length){
        printf("Error in reading file: %s ,%d\n",__FILE__,__LINE__);
        return 2;
    }
    char * strfnx = strstr(ildg_format_str,"<lx>");
    char * strfny = strstr(ildg_format_str,"<ly>");
    char * strfnz = strstr(ildg_format_str,"<lz>");
    char * strfnt = strstr(ildg_format_str,"<lt>");
    int nx_r,ny_r,nz_r,nt_r;
    int allfound = 1;
    if(strfnx!= NULL){sscanf(strfnx," <lx> %d </lx> ",&nx_r);printf("nx: %d\n",nx_r);}
    else{allfound = 0;}
    if(strfny!= NULL){sscanf(strfny," <ly> %d </ly> ",&ny_r);printf("ny: %d\n",ny_r);}
    else{allfound = 0;}
    if(strfnz!= NULL){sscanf(strfnz," <lz> %d </lz> ",&nz_r);printf("nz: %d\n",nz_r);}
    else{allfound = 0;}
    if(strfnt!= NULL){sscanf(strfnt," <lt> %d </lt> ",&nt_r);printf("nt: %d\n",nt_r);}
    else{allfound = 0;}
    if(!allfound){
        printf("Error, %s:%d : lx,ly,lz or lt not found in \"ildg-format\"\n",
                __FILE__,__LINE__);
        printf("Please check. Exiting now.\n");
        return 2;
    }
    // check lattice dimensions
    if((geom_par.gnx!=nx_r)||(geom_par.gny!=ny_r)||
            (geom_par.gnz!=nz_r)||(geom_par.gnt!=nt_r)){
        printf("Error, Configuration dimensions not compatible with input file.\n");
        return 2;
    }


    // read MD_traj
    if(MD_traj_index != -1){
        fseeko(fg,ildg_header_ends_positions[MD_traj_index],SEEK_SET);
        reads = fscanf(fg,"%d",conf_id_iter);
        if(reads!= 1){
            printf("Error in reading file: %s ,%d\n",__FILE__,__LINE__);
            return 2;
        }
    }else *conf_id_iter = 1;

    // read ildg-binary-data (su3 gauge conf)
    if(verbosity_lv>2)
        printf("Reading ildg-binary-data...\n");



    off_t ibd_start = ildg_header_ends_positions[ildg_binary_data_index];

    fseeko(fg,ibd_start, SEEK_SET);
    rw_iterate_on_global_sites_lx_xyzt_axis_ordering(
            binaryread_single_su3_into_su3_soa,(void*)conf,fg,0);


    fclose(fg);

    return 0;

}

int print_su3_soa_ildg_binary(global_su3_soa * const conf, 
        const char* nomefile, int conf_id_iter)
{

    if(geom_par.initialized_check == 0){
        printf("ERROR: GEOMETRY MUST BE INITIALIZED!\n");
        printf("ERROR: in %s at line %d\n",__FILE__, __LINE__);
        return 2;
    } 


    printf("Writing ildg conf...\n");
    int writes; FILE *fp;
    fp = fopen(nomefile,"w");
    if(! fp ){
        printf("ERROR, %s unreadable.\n",nomefile);
        return 2;
    }

    int nx = geom_par.gnx;
    int ny = geom_par.gny;
    int nz = geom_par.gnz;
    int nt = geom_par.gnt;





    int  conf_machine_endianness_disagreement = machine_is_little_endian();
    char ildg_format_str[1000];  

    sprintf(ildg_format_str,
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
            <ildgFormat xmlns=\"http://www.lqcd.org/ildg\"\n\
            xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n\
            xsi:schemaLocation=\"http://www.lqcd.org/ildg filefmt.xsd\">\n\
            <version>1.0</version>\n\
            <field>su3gauge</field>\n\
            <precision>64</precision>\n\
            <lx>%d</lx>\n\
            <ly>%d</ly>\n\
            <lz>%d</lz>\n\
            <lt>%d</lt>\n\
            </ildgFormat>",nx,ny,nz,nt);
    uint64_t len = strlen(ildg_format_str) ;
    uint32_t magic_number =  0x456789ab;
    uint16_t version = 1;
    uint16_t mbme_flag = 0;
    off_t missing_bytes = (off_t)(len%8==0 ? 0:8-len%8);
    uint64_t len64 = len + missing_bytes;

    if(conf_machine_endianness_disagreement){
        uint32swe(&magic_number);
        uint16swe(&version);
        uint16swe(&mbme_flag);
        uint64swe(&len64);
    }
    ILDG_header ildg_format_header = 
        (ILDG_header){magic_number,version,mbme_flag,len64,"ildg-format"};
    fwrite(&ildg_format_header,sizeof(ILDG_header),1,fp);
    fwrite(ildg_format_str,1,len,fp);

    char pad[8] = {0,0,0,0,0,0,0,0};

    // padding to 8 bytes
    fwrite(pad,1,missing_bytes,fp);


    // writing trajectory number
    char MD_traj_str[15];
    sprintf(MD_traj_str,"%d",conf_id_iter);
    len = strlen(MD_traj_str);
    missing_bytes = (off_t)(len%8==0 ? 0:8-len%8);
    len64 = len + missing_bytes;
    if(verbosity_lv > 3) printf("Writing MD_traj, %" PRIu64 "\n", len64);
    if(conf_machine_endianness_disagreement) uint64swe(&len64);
    ILDG_header MD_traj_header = 
        (ILDG_header){magic_number,version,mbme_flag,len64,"MD_traj"};
    fwrite(&MD_traj_header,sizeof(ILDG_header),1,fp);
    fwrite(MD_traj_str,1,len,fp);

    // padding to 8 bytes
    fwrite(pad,1,missing_bytes,fp);

    // writing input file in the conf
    // the string 'input_file_str' should be initialized in 
    // ../Include/setting_file_parser.c
    len = strlen(input_file_str);
    missing_bytes = (off_t)(len%8==0 ? 0:8-len%8);
    len64 = len + missing_bytes;
    if(verbosity_lv > 3) printf("Writing input file into the conf: %" PRIu64 
            "\n", len64);
    if(conf_machine_endianness_disagreement) uint64swe(&len64);
    ILDG_header input_file_header = 
        (ILDG_header){magic_number,version,mbme_flag,len64,"input-file"};
    fwrite(&input_file_header,sizeof(ILDG_header),1,fp);
    fwrite(input_file_str,1,len,fp);

    // padding to 8 bytes
    fwrite(pad,1,missing_bytes,fp);


    // writing ildg binary data, that is conf
    len = (uint64_t) geom_par.gnx * geom_par.gny * geom_par.gnz * geom_par.gnt *4*3*3*2*8;
    missing_bytes = (off_t)(len%8==0?0:8-len%8);
    len64 = len + missing_bytes;
    if(verbosity_lv > 3) printf("Writing ildg-binary-data, %" PRIu64 
            "\n", len64);
    if(conf_machine_endianness_disagreement) uint64swe(&len64);

    ILDG_header ildg_binary_data_header = 
        (ILDG_header){magic_number,version,mbme_flag,len64,"ildg-binary-data"};
    fwrite(&ildg_binary_data_header,sizeof(ILDG_header),1,fp);


    rw_iterate_on_global_sites_lx_xyzt_axis_ordering(
            binarywrite_single_su3_into_su3_soa, conf, fp , 0);

    // padding to 8 bytes
    fwrite(pad,1,missing_bytes,fp);


    len64 = 0;
    if(conf_machine_endianness_disagreement) uint64swe(&len64);
    ILDG_header ildg_data_lfn_header = 
        (ILDG_header){magic_number,version,mbme_flag,len64,"ildg-data-lfn"};
    fwrite(&ildg_data_lfn_header,sizeof(ILDG_header),1,fp);

    fclose(fp);


    return  0;
}





void rw_iterate_on_global_sites_lx_xyzt_axis_ordering( 
        void (*single_element_rw)(
            const int /*idxh machine*/, const int /*parity*/, 
            const int /*dirmachine*/, const void* /*data*/, 
            const int /*conf_machine_endianness_disagreement*/,
            FILE * /*fp*/), 
        const void* datastruct, FILE * fp, const int scalar_even_mode){
    // fp must be already in the right position

    const int nx = geom_par.gnx;
    const int ny = geom_par.gny;
    const int nz = geom_par.gnz;
    const int nt = geom_par.gnt;
    const int conf_machine_endianness_disagreement = machine_is_little_endian();

    int x,y,z,t,dir;
    // iterating on the sites and directions in the order expected 
    // for the file 
    int xtmp,xlimit,dirlimit;
    if (scalar_even_mode){
        xlimit = nx/2;
        dirlimit = 1;

    }else {
        xlimit = nx;
        dirlimit = 4;
    }

    for(t=0;t<nt;t++) for(z=0;z<nz;z++)
        for(y=0;y<ny;y++) for(xtmp=0;xtmp<xlimit;xtmp++)
        {
            // ILDG format on disk is not transposed,
            //  but conf in machine memory could be.
            int parity;
            if(scalar_even_mode){
                parity = 0; // even
                x = 2*xtmp + ((y+z+t)&0x1);
            }
            else {
                x = xtmp;
                parity = (x+y+z+t)%2;// parity = (d0+d1+d2+d3)%2;
    						#if defined(MULTIDEVICE) && defined(GAUGE_ACT_WILSON)
								parity = !parity;
    						#endif
            }
            int xs[4] = {x,y,z,t};
            int d[4];
            for(dir = 0; dir<4;dir++) d[dir] = xs[geom_par.d0123map[dir]];

#ifdef MULTIDEVICE
            int idxh_machine = gl_to_gl_snum(d[0],d[1],d[2],d[3]);
#else
            int idxh_machine = snum_acc(d[0],d[1],d[2],d[3]);
#endif
            for(dir=0;dir<dirlimit;dir++){
							int dirmachine = geom_par.xyztmap[dir];
							single_element_rw(idxh_machine,parity,dirmachine,datastruct,
																conf_machine_endianness_disagreement,fp);
            }
        }
}

void binaryread_single_su3_into_su3_soa( // machine big endian
        const int idxh_machine, const int parity, const int dirmachine, 
        const void* datastruct,
        const int conf_machine_endianness_disagreement, FILE* fp){

    global_su3_soa * conf = (global_su3_soa *) datastruct;

    single_su3 m;          
    int reads = fread((void*)m.comp,sizeof(double),18,fp);
    if(reads!= 18){
        printf("Error in reading file: %s ,%d, idxh: %d, dir %d , reads %d\n",
                __FILE__,__LINE__, idxh_machine, dirmachine,reads );
        exit(2);
    }

    // fix enddiannes
    int irow,icol;
    for(irow=0;irow<3;irow++) for(icol=0;icol<3;icol++){
        double re = creal(m.comp[irow][icol]);
        double im = cimag(m.comp[irow][icol]);
        if(conf_machine_endianness_disagreement){
            doubleswe(&re);

            doubleswe(&im);

            m.comp[irow][icol] = re + im*I;
        }

    }
    rebuild3row(&m);
    single_su3_into_global_su3_soa(&conf[2*dirmachine+parity],idxh_machine,&m);
}

void binarywrite_single_su3_into_su3_soa(
        const int idxh_machine, const int parity, const int dirmachine, 
        const void* datastruct,
        const int conf_machine_endianness_disagreement, FILE* fp){

    global_su3_soa * conf = (global_su3_soa *) datastruct;

    single_su3 m;          
    single_su3_from_global_su3_soa(&conf[2*dirmachine+parity],idxh_machine,&m);
    rebuild3row(&m);

    // fix enddiannes
    int irow,icol;
    for(irow=0;irow<3;irow++) for(icol=0;icol<3;icol++){
        double re = creal(m.comp[irow][icol]);
        double im = cimag(m.comp[irow][icol]);
        if(conf_machine_endianness_disagreement){
            doubleswe(&re);
            doubleswe(&im);
        }
        int writes =  fwrite((void*)&re,sizeof(double),1,fp);
        if(writes!= 1){
            printf("Error in writing file: %s ,%d\n",__FILE__,__LINE__);
            exit(2);
        }
        writes =  fwrite((void*)&im,sizeof(double),1,fp);
        if(writes!= 1){
            printf("Error in writing file: %s ,%d\n",__FILE__,__LINE__);
            exit(2);
        }
    }
}




#endif
