#ifndef DBGTOOLS_MPI_H_
#define DBGTOOLS_MPI_H_


#include "../DbgTools/debug_tools.c"

THE FOLLOWING MUST BE FIXED // DEBUG

 // FERMIONS

// init-finalize
    void load_global_vec3_soa_from_file(global_vec3_soa * fermion, const char * filename){

        FILE *fp;
        double ar, ai, br, bi, cr, ci;
        int i = 0;
        int error = 0;

        fp = fopen(filename, "rt");

        if (fp == NULL) {
            printf("Could not open file %s \n", filename);
            exit(-1);
        }

        while ((i < GL_SIZEH) && (!error)) {

            if (fscanf(fp, "(%lf,%lf) (%lf,%lf) (%lf,%lf) \n", &ar, &ai, &br, &bi, &cr, &ci) == 6) {
                fermion->c0[i] = (ar + ai * I);
                fermion->c1[i] = (br + bi * I);
                fermion->c2[i] = (cr + ci * I);
            } else {
                printf("Read error... ");
                error = 1;
            }

            i++;

        }

        printf("Read %d vectors from file %s \n", i, filename);

        fclose(fp);

    }
void write_global_vec3_soa_to_file(const global_vec3_soa * fermion, const char * filename){

 FILE *fp;
  double ar, ai, br, bi, cr, ci;
  int idx;
  int i = 0;
  int error = 0;

  fp = fopen(filename, "wt");

  if (fp == NULL) {
    printf("Could not open file %s \n", filename);
    exit(-1);
  }


    while ( (i < GL_SIZEH) && (!error) ) {
    
      ar  = creal(fermion->c0[i]);
      br  = creal(fermion->c1[i]);
      cr  = creal(fermion->c2[i]);
      ai  = cimag(fermion->c0[i]);
      bi  = cimag(fermion->c1[i]);
      ci  = cimag(fermion->c2[i]);

      fprintf(fp, "(%lf,%lf) (%lf,%lf) (%lf,%lf)\n", ar, ai, br, bi, cr, ci);

     i++;

    }

  printf("Written %d vectors to file %s \n", i, filename);

  fclose(fp);

}
void read_lnh_vec3_soa_from_file(lnh_vec3_soa * u, const char *lnh_filename){

    // TO WRITE

}
void write_lnh_vec3_soa_to_file(const lnh_vec3_soa * lnh_fermion, const char * global_filename){

    FILE *fp;
    double ar, ai, br, bi, cr, ci;
    int i = 0;
    int error = 0;

    // defining a filename for the sublattice file 
    char filename[50] = "rank";
    char rank_str[5];
    snprintf(rank_str, 5*sizeof(char),"%d",myrank);
    char n_ranks_str[5];
    snprintf(n_ranks_str, 5*sizeof(char),"%d",n_ranks);

    strcat(filename, rank_str);
    strcat(filename, "of");
    strcat(filename, n_ranks_str);    
    strcat(filename, "_");
    strcat(filename,global_filename);//global filename appended

    fp = fopen(filename, "wt");
    
    if (fp == NULL) {
        printf("Could not open file %s \n", filename);
        exit(-1);
    }

    vec4int sl_glo = gl_loc_origin_from_rank(myrank);
  //Sub Lattice Global (coordinate of) Local Origin
    printf("#MPI%02d: Writing file %s, \n",myrank, filename);
    printf("%d %d %d %d\n", sl_glo.x, sl_glo.y, sl_glo.z, sl_glo.t); 
    
    while ( (i < LNH_SIZEH) && (!error) ) {

        ar  = creal(lnh_fermion->c0[i]);
        br  = creal(lnh_fermion->c1[i]);
        cr  = creal(lnh_fermion->c2[i]);
        ai  = cimag(lnh_fermion->c0[i]);
        bi  = cimag(lnh_fermion->c1[i]);
        ci  = cimag(lnh_fermion->c2[i]);

       fprintf(fp, "(%lf,%lf) (%lf,%lf) (%lf,%lf)\n", ar, ai, br, bi, cr, ci);
       i++;

    }

    printf("#MPI%02d: Written %d vectors to file %s \n",myrank, i, filename);

    fclose(fp);

}
void global_vec3_soa_init(int initmode, global_vec3_soa * fermion){
    if( initmode == 0) // Read from file
        load_global_vec3_soa_from_file(fermion, "fermion");
    else{
       
        int i = 0;
        while( i < GL_SIZEH ) {
            
            vec3 v;

            if (initmode == 2 ) v = zero_vec3;
            else if (initmode == 1 ) v = vec3_gauss();
            else if (initmode == 3 ) v = generator_vec3_debug();

            fermion->c0[i] = v.c0;
            fermion->c1[i] = v.c1;
            fermion->c2[i] = v.c2;

            i++;
        }
    }
}
void lnh_vec3_soa_init(int initmode, lnh_vec3_soa * fermion){
       
        int i = 0;
        while( i < LNH_SIZEH ) {
            
            vec3 v;

            if (initmode == 0 ) v = zero_vec3;
            else if (initmode == 1 ) v = vec3_gauss();
            else if (initmode == 3 ) v = generator_vec3_debug();

            fermion->c0[i] = v.c0;
            fermion->c1[i] = v.c1;
            fermion->c2[i] = v.c2;

            i++;
        }
}


// communications
void send_lnh_subfermion_to_rank(const global_vec3_soa * gl_soa_fermion, int target_rank, int tag){

    // USE ONLY FROM MASTER RANK 
    // Notice that this function cannot tell the difference between
    // even and odd fermions
    //target sublattice information
    vec4int target_gl_loc_origin_from_rank = gl_loc_origin_from_rank(target_rank);

    // building sublattice duplicate, target_conf
    lnh_vec3_soa* target_vec3_soa = malloc(sizeof(lnh_vec3_soa)); 

    int tg_lnh_xh,tg_lnh_y,tg_lnh_z,tg_lnh_t; //target-lnh coordinates

    // Copying all relevant vec3s into the sublattice
    for(tg_lnh_t=0;tg_lnh_t<LNH_NT; tg_lnh_t++)
    for(tg_lnh_z=0;tg_lnh_z<LNH_NZ; tg_lnh_z++)
    for(tg_lnh_y=0;tg_lnh_y<LNH_NY; tg_lnh_y++)
    for(tg_lnh_xh=0;tg_lnh_xh<LNH_NXH; tg_lnh_xh++){

        int tg_lnh_x = 2* tg_lnh_xh; // CHECK


        int target_gl_snum = target_lnh_to_gl_snum(tg_lnh_x, tg_lnh_y, tg_lnh_z, tg_lnh_t, target_gl_loc_origin_from_rank);
        int target_lnh_snum = snum_acc(tg_lnh_x, tg_lnh_y, tg_lnh_z, tg_lnh_t);

        vec3 aux = vec3_from_global_vec3_soa(gl_soa_fermion,target_gl_snum);
        vec3_into_vec3_soa(aux,target_vec3_soa,target_lnh_snum);

    }

    //sending the subfermion
    MPI_Send(target_vec3_soa, 3*2*LNH_SIZEH , MPI_DOUBLE, target_rank, tag, devinfo.mpi_comm);

}
static inline void recv_lnh_subfermion_from_master(lnh_vec3_soa * lnh_fermion, int tag){

    MPI_Recv(lnh_fermion, 3*2*LNH_SIZEH,MPI_DOUBLE,0,tag,devinfo.mpi_comm,MPI_STATUS_IGNORE );
    // ^^ CHECK
}
static inline void send_lnh_subfermion_to_master(lnh_vec3_soa * lnh_fermion, int tag){

    MPI_Send(lnh_fermion, 3*2*LNH_SIZEH,MPI_DOUBLE,0,tag,devinfo.mpi_comm );
    // ^^ CHECK
}
void recv_loc_subfermion_from_rank(const global_vec3_soa * gl_soa_fermion, int rank, int tag){
    // USE ONLY FROM MASTER RANK 
    // Notice that this function cannot tell the difference between
    // even and odd fermions

    //target sublattice information
    vec4int origin_gl_loc_origin_from_rank = gl_loc_origin_from_rank(rank);
    // building sublattice duplicate, target_conf
    lnh_vec3_soa* origin_vec3_soa = malloc(sizeof(lnh_vec3_soa)); 
    //receive the subfermion
    MPI_Recv(origin_vec3_soa, 3*2*LNH_SIZEH , MPI_DOUBLE, rank, tag, devinfo.mpi_comm, MPI_STATUS_IGNORE);
    // ^^ CHECK

    int or_loc_xh,or_loc_y,or_loc_z,or_loc_t; //target-lnh coordinates

    // NOTICE: only LOC sites are important, thus cycling only on terget local sites

    // Copying all relevant vec3s into the sublattice
    for(or_loc_t=0;or_loc_t<LOC_NT; or_loc_t++)
    for(or_loc_z=0;or_loc_z<LOC_NZ; or_loc_z++)
    for(or_loc_y=0;or_loc_y<LOC_NY; or_loc_y++)
    for(or_loc_xh=0;or_loc_xh<LOC_NXH; or_loc_xh++){


        int or_loc_x = 2* or_loc_xh; // CHECK
        int or_lnh_x,or_lnh_y,or_lnh_z,or_lnh_t; //target-lnh coordinates
        
        or_lnh_x = or_loc_x + X_HALO;
        or_lnh_y = or_loc_y + Y_HALO;
        or_lnh_z = or_loc_z + Z_HALO;
        or_lnh_t = or_loc_t + T_HALO;

        int origin_gl_snum = target_lnh_to_gl_snum(or_lnh_x, or_lnh_y, or_lnh_z, or_lnh_t, origin_gl_loc_origin_from_rank);
        int origin_lnh_snum = snum_acc(or_lnh_x, or_lnh_y, or_lnh_z, or_lnh_t);

        vec3 aux = vec3_from_vec3_soa(origin_vec3_soa,origin_lnh_snum);
        vec3_into_global_vec3_soa(aux,gl_soa_fermion,origin_gl_snum);
    }
}



// CONF

void load_global_su3_soa_from_file(global_su3_soa * u, const char *filename){

    FILE *fp;
    int nx_l, ny_l, nz_l, nt_l, update_iterations;
    double beta_l, mass_l, no_flavours_l;
    double ar, ai, br, bi, cr, ci;
    int idx;
    int i = 0;
    int j = 0;
    int error = 0;

    fp = fopen(filename, "rt");

    if (fp == NULL) {
        printf("Could not open file %s \n", filename);
        exit(-1);
    }


    fscanf(fp, "%d %d %d %d %lf %lf %lf %d \n", &nx_l, &ny_l, &nz_l, &nt_l, 
            &beta_l, &mass_l, &no_flavours_l, 
            &update_iterations);

    printf("Reading configuration file with header: \n");
    printf("nx_l: %d, ny_l: %d, nz_l: %d, nt_l: %d \n", nx_l, ny_l, nz_l, nt_l); 
    printf("beta_l: %lf, mass_l: %lf, no_flavours_l: %lf \n", beta_l, mass_l, no_flavours_l);
    printf("update_iterations: %d \n", update_iterations);

    while ( (i < GL_SIZEH * 8) && (!error) ) {

        j = i / GL_SIZEH;
        idx = i % GL_SIZEH;

        if (fscanf(fp, "(%lf,%lf) (%lf,%lf) (%lf,%lf) \n", &ar, &ai, &br, &bi, &cr, &ci) == 6) {
            u[j].r0.c0[idx] = (ar + ai * I);
            u[j].r0.c1[idx] = (br + bi * I);
            u[j].r0.c2[idx] = (cr + ci * I);
        } else {
            printf("Read error... ");
            error = 1;
        }

        if (fscanf(fp, "(%lf,%lf) (%lf,%lf) (%lf,%lf) \n", &ar, &ai, &br, &bi, &cr, &ci) == 6) {
            u[j].r1.c0[idx] = (ar + ai * I);
            u[j].r1.c1[idx] = (br + bi * I);
            u[j].r1.c2[idx] = (cr + ci * I);
        } else {
            printf("Read error... ");
            error = 1;
        }

        if (fscanf(fp, "(%lf,%lf) (%lf,%lf) (%lf,%lf) \n", &ar, &ai, &br, &bi, &cr, &ci) == 6) {
            u[j].r2.c0[idx] = (ar + ai * I);
            u[j].r2.c1[idx] = (br + bi * I);
            u[j].r2.c2[idx] = (cr + ci * I);
        } else {
            printf("Read error... ");
            error = 1;
        }

        i++;

    }

    printf("Read %d matrices from file %s \n", i, filename);

    fclose(fp);

}
void write_global_su3_soa_to_file(global_su3_soa * u, const char *filename){

    FILE *fp;
    double ar, ai, br, bi, cr, ci;
    int idx;
    int i = 0;
    int j = 0;
    int error = 0;

    fp = fopen(filename, "wt");

    if (fp == NULL) {
        printf("Could not open file %s \n", filename);
        exit(-1);
    }


    fprintf(fp, "%d %d %d %d %lf %lf %d %d \n", GL_NX, GL_NY, GL_NZ, GL_NT, 
            beta, mass, no_flavours, 
            update_iteration);

    printf("Writing configuration file with header: \n");
    printf("nx: %d, ny: %d, nz: %d, nt: %d \n",GL_NX, GL_NY, GL_NZ, GL_NT); 
    printf("beta: %lf, mass: %lf, no_flavours: %d \n", beta, mass, no_flavours);
    printf("update_iteration: %d \n", update_iteration);

    while ( (i < GL_SIZEH * 8) && (!error) ) {

        j = i / GL_SIZEH;
        idx = i % GL_SIZEH;

        ar  = creal(u[j].r0.c0[idx]);
        br  = creal(u[j].r0.c1[idx]);
        cr  = creal(u[j].r0.c2[idx]);
        ai = cimag( u[j].r0.c0[idx]);
        bi = cimag( u[j].r0.c1[idx]);
        ci = cimag( u[j].r0.c2[idx]);

        fprintf(fp, "(%lf,%lf) (%lf,%lf) (%lf,%lf) \n", ar, ai, br, bi, cr, ci);
        //      fprintf(fp, "(%d,%d) (%d,%d) (%d,%d) \n",(int) ar, (int)ai, (int)br, (int)bi, (int)cr, (int)ci);  // DEBUG
        ar  = creal(u[j].r1.c0[idx]);
        br  = creal(u[j].r1.c1[idx]);
        cr  = creal(u[j].r1.c2[idx]);
        ai = cimag( u[j].r1.c0[idx]);
        bi = cimag( u[j].r1.c1[idx]);
        ci = cimag( u[j].r1.c2[idx]);

        fprintf(fp, "(%lf,%lf) (%lf,%lf) (%lf,%lf) \n", ar, ai, br, bi, cr, ci);
        //      fprintf(fp, "(%d,%d) (%d,%d) (%d,%d) \n",(int) ar, (int)ai, (int)br, (int)bi, (int)cr, (int)ci);  // DEBUG
        ar  = creal(u[j].r2.c0[idx]);
        br  = creal(u[j].r2.c1[idx]);
        cr  = creal(u[j].r2.c2[idx]);
        ai = cimag( u[j].r2.c0[idx]);
        bi = cimag( u[j].r2.c1[idx]);
        ci = cimag( u[j].r2.c2[idx]);

        fprintf(fp, "(%lf,%lf) (%lf,%lf) (%lf,%lf) \n", ar, ai, br, bi, cr, ci);
        //      fprintf(fp, "(%d,%d) (%d,%d) (%d,%d) \n",(int) ar, (int)ai, (int)br, (int)bi, (int)cr, (int)ci);  // DEBUG
        i++;
    }

    printf("Written %d matrices to file %s \n", i, filename);
    fclose(fp);

}
void read_lnh_su3_soa_from_file(const lnh_su3_soa * u, const char *lnh_filename){

    // TO WRITE

}

void write_lnh_su3_soa_to_file(lnh_su3_soa * u, const char *global_filename, int rank){

    FILE *fp;
    double ar, ai, br, bi, cr, ci;
    int idx;
    int i = 0;
    int j = 0;
    int error = 0;

    // defining a filename for the sublattice file
    char filename[50] = "rank";
    char rank_str[5]; snprintf(rank_str, 5*sizeof(char), "%d", rank);
    char n_ranks_str[5]; snprintf(n_ranks_str, 5*sizeof(char), "%d",NRANKS);

    strcat(filename, rank_str);
    strcat(filename, "of");
    strcat(filename, n_ranks_str);
    strcat(filename, "_");
    strcat(filename,global_filename);//global filename appended
    fp = fopen(filename, "wt");

    if (fp == NULL) {
        printf("Could not open file %s \n", filename);
        exit(-1);
    }
    // HEADER STUFF 
    fprintf(fp, "%d %d %d %d %lf %lf %d %d \n", GL_NX, GL_NY, GL_NZ, GL_NT, beta, mass, no_flavours, update_iteration);
    fprintf(fp, "%d %d %d %d\n", LNH_NX, LNH_NY, LNH_NZ, LNH_NT); 
    vec4int sl_glo = gl_loc_origin_from_rank(rank);
    //Sub Lattice Global (coordinate of) Local Origin
    fprintf(fp, "%d %d %d %d\n", sl_glo.x, sl_glo.y, sl_glo.z, sl_glo.t); 
    /*  printf("Writing file %s, \n", filename);
        printf("Writing configuration file with header: \n");
        printf("nx: %d, ny: %d, nz: %d, nt: %d \n",GL_NX, GL_NY, GL_NZ, GL_NT); 
        printf("beta: %lf, mass: %lf, no_flavours: %lf \n", beta, mass, no_flavours);
        printf("update_iteration: %d \n", update_iteration);
        printf("LNH sub lattice dimensions:\n"); 
        printf("LNH_NX:%d,LNH_NY:%d,LNH_NZ:%d,LNH_NT:%d\n", LNH_NX, LNH_NY, LNH_NZ, LNH_NT);
        printf("SubLattice Global coordinates of Local Origin:\n");
        */
    printf("%d %d %d %d\n", sl_glo.x, sl_glo.y, sl_glo.z, sl_glo.t); 
    // BODY STUFF
    while ( (i < LNH_SIZEH * 8) && (!error) ){

        // HERE PROBABLY IS THE FUCKING BASTARD! 

        j = i / LNH_SIZEH;
        idx = i %  LNH_SIZEH ;

        ar = creal(u[j].r0.c0[idx]);
        br = creal(u[j].r0.c1[idx]);
        cr = creal(u[j].r0.c2[idx]);
        ai = cimag(u[j].r0.c0[idx]);
        bi = cimag(u[j].r0.c1[idx]);
        ci = cimag(u[j].r0.c2[idx]);

        fprintf(fp, "(%lf,%lf) (%lf,%lf) (%lf,%lf) \n", ar, ai, br, bi, cr, ci);
        //      fprintf(fp, "(%d,%d) (%d,%d) (%d,%d) \n",(int) ar, (int)ai, (int)br, (int)bi, (int)cr, (int)ci);  // DEBUG

        ar = creal(u[j].r1.c0[idx]);
        br = creal(u[j].r1.c1[idx]);
        cr = creal(u[j].r1.c2[idx]);
        ai = cimag(u[j].r1.c0[idx]);
        bi = cimag(u[j].r1.c1[idx]);
        ci = cimag(u[j].r1.c2[idx]);

        fprintf(fp, "(%lf,%lf) (%lf,%lf) (%lf,%lf) \n", ar, ai, br, bi, cr, ci);  
        //      fprintf(fp, "(%d,%d) (%d,%d) (%d,%d) \n",(int) ar, (int)ai, (int)br, (int)bi, (int)cr, (int)ci);  // DEBUG

        ar = creal(u[j].r2.c0[idx]);
        br = creal(u[j].r2.c1[idx]);
        cr = creal(u[j].r2.c2[idx]);
        ai = cimag(u[j].r2.c0[idx]);
        bi = cimag(u[j].r2.c1[idx]);
        ci = cimag(u[j].r2.c2[idx]);

        fprintf(fp, "(%lf,%lf) (%lf,%lf) (%lf,%lf) \n", ar, ai, br, bi, cr, ci);
        //      fprintf(fp, "(%d,%d) (%d,%d) (%d,%d) \n",(int) ar, (int)ai, (int)br, (int)bi, (int)cr, (int)ci);  // DEBUG

        i++;
    }
    printf("#MPI%02d: Written %d matrices to file %s \n",myrank, i, filename);

    fclose(fp);

}

void global_su3_soa_init(int initmode, global_su3_soa *global_soa_conf){

    if (initmode == 2 )
        load_global_su3_soa_from_file(global_soa_conf, "conf");
    else{

        int gl_x, gl_y,gl_z,gl_t,dir,parity,gl_snum_index;
        for(dir =0; dir < 4; dir++)
            for(gl_t=0; gl_t < GL_NT; gl_t++)for(gl_z=0; gl_z < GL_NZ; gl_z++)
                for(gl_y=0; gl_y < GL_NY; gl_y++)for(gl_x=0; gl_x < GL_NX; gl_x++){

                    //          printf("%d %d %d %d %d %d \n",parity,dir,gl_x, gl_y, gl_z, gl_t);

                    su3 aux;
                    // Producing a single su3 matrix and dispatching it to the right position
                    if(initmode==0) aux = su3_cold_random_matrix();//start from ordered conf 
                    else if(initmode==1) aux = su3_random_matrix();//start from random (hot) conf

#ifdef DEBUG            
                    //ONLY IN TEST/DEBUG BRANCH
                    else if(initmode==3) aux = generator_3x3_debug();
#endif

#ifndef IGNORE_ETA            
                    aux = su3_scal_mult(aux, gl_eta(gl_x,gl_y,gl_z,gl_t,dir));
#endif
                    parity = (gl_x + gl_y + gl_z + gl_t ) % 2;
                    gl_snum_index = gl_to_gl_snum(gl_x, gl_y, gl_z, gl_t);
                    su3_into_global_su3_soa(aux,global_soa_conf,dir,parity,gl_snum_index);
                }
    }
}    

void lnh_su3_soa_init(int initmode, lnh_su3_soa *lnh_soa_conf){

    if (initmode == 2 )
//        load_lnh_su3_soa_from_file(lnh_soa_conf, "conf");
        printf("MPI%02d: PARALLEL CONF READING NOT IMPLEMENTED YET\n", myrank);
    else{

        int dir,parity,lnh_snum_index;
        int loc_x, loc_y, loc_z, loc_t;
    
        for(dir =0; dir < 4; dir++)
        for(loc_t=0; loc_t<LOC_NT; loc_t++) {
            for(loc_z=0; loc_z<LOC_NZ; loc_z++) {
                for(loc_y=0; loc_y<LOC_NY; loc_y++) {
                    for(loc_x=0; loc_x <LOC_NX; loc_x++) {
    
                        int lnh_x , lnh_y, lnh_z, lnh_t, eta;
    
                        lnh_x = loc_x + X_HALO;
                        lnh_y = loc_y + Y_HALO;
                        lnh_z = loc_z + Z_HALO;
                        lnh_t = loc_t + T_HALO;
                        
                        su3 aux;
                        // Producing a single su3 matrix and dispatching it to the right position
                        if(initmode==0) aux = su3_cold_random_matrix();//start from ordered conf 
                        else if(initmode==1) aux = su3_random_matrix();//start from random (hot) conf
    
    #ifdef DEBUG            
                        //ONLY IN TEST/DEBUG BRANCH
                        else if(initmode==3) aux = generator_3x3_debug();
    #endif
    #ifndef IGNORE_ETA            
                        aux = su3_scal_mult(aux, gl_eta(loc_x,loc_y,loc_z,loc_t,dir)); // These coordinated should be gl_x etc, but we can exploit the 
         // periodicity of eta
    #endif
    
                        parity = (loc_x + loc_y + loc_z + loc_t ) % 2; // Same applies here as for eta
                        lnh_snum_index = snum_acc(lnh_x, lnh_y, lnh_z, lnh_t);
                        su3_into_su3_soa(aux,lnh_soa_conf,dir,parity,lnh_snum_index);
                    } //x
                } //y
            } //z
        }  //t
    } //else
}

void send_lnh_subconf_to_buffer(global_su3_soa *gl_soa_conf, lnh_su3_soa *lnh_conf, int target_rank){
    // USE ONLY FROM MASTER RANK


    //target sublattice information
    vec4int target_gl_loc_origin_from_rank = gl_loc_origin_from_rank(target_rank);
    /*
       int target_loc_origin_parity = (target_gl_loc_origin_from_rank.x +
       target_gl_loc_origin_from_rank.y +
       target_gl_loc_origin_from_rank.z +
       target_gl_loc_origin_from_rank.t)%2;// with the current setup, this should always be 0

       if(target_loc_origin_parity) printf("Problems\n");
       */
    // building sublattice duplicate, target_conf
    lnh_su3_soa* target_su3_soa = lnh_conf;

    int tg_lnh_x,tg_lnh_y,tg_lnh_z,tg_lnh_t,dir; //target-lnh coordinates
    // and link direction
    // Copying all relevant links into the sublattice
    for(dir =0; dir < 4; dir++)
        for(tg_lnh_t=0;tg_lnh_t<LNH_NT; tg_lnh_t++)
            for(tg_lnh_z=0;tg_lnh_z<LNH_NZ; tg_lnh_z++)
                for(tg_lnh_y=0;tg_lnh_y<LNH_NY; tg_lnh_y++)
                    for(tg_lnh_x=0;tg_lnh_x<LNH_NX; tg_lnh_x++){

                        //        int gtsp; // global target site parity
                        int tsprlo ; // target site parity respect (to his) local origin;

                        int target_gl_snum = target_lnh_to_gl_snum(tg_lnh_x, tg_lnh_y, tg_lnh_z, tg_lnh_t, target_gl_loc_origin_from_rank);
                        int target_lnh_snum = snum_acc(tg_lnh_x, tg_lnh_y, tg_lnh_z, tg_lnh_t);

                        tsprlo = (X_HALO+Y_HALO+Z_HALO+T_HALO+ tg_lnh_t+tg_lnh_z+tg_lnh_y+tg_lnh_x)%2;
                        //      gtsp = (target_loc_origin_parity + tsprlo )%2;

                        //      printf("%d %d  %d  %d  %d  %d  %d ",dir, tg_lnh_t, tg_lnh_z, tg_lnh_y, tg_lnh_x, target_gl_snum, tsprlo);
                        su3 aux = su3_from_global_su3_soa(gl_soa_conf,dir,tsprlo,target_gl_snum);

                        //    printf("ciao \n");
                        //        printf("ciao ");
                        //      print_su3(aux);
                        su3_into_su3_soa(aux,target_su3_soa,dir,tsprlo,target_lnh_snum);
                        //        printf("ciao \n");

                    }

    /*
    //sending the subconfiguration
    MPI_Send(target_su3_soa, 2*4*(6*3)*LNH_SIZEH,MPI_DOUBLE, rank, 0, devinfo.mpi_comm);
    // In case we remove the third line, this is possibly the thing to do
    //MPI_Send(target_su3_soa, 2*4*(6*2)*LNH_SIZEH,MPI_DOUBLE, rank, 0, devinfo.mpi_comm);
    //Or maybe do multiple sends in a more complicated way
    // ^^ CHECK
    */ // COMMENTED FOR TEST
}

void recv_loc_subconf_from_buffer(global_su3_soa *gl_soa_conf, lnh_su3_soa* lnh_conf, int target_rank){
    // USE ONLY FROM MASTER RANK
    //target sublattice information
    vec4int target_gl_loc_origin_from_rank = gl_loc_origin_from_rank(target_rank);
    /*
       int target_loc_origin_parity = (target_gl_loc_origin_from_rank.x +
       target_gl_loc_origin_from_rank.y +
       target_gl_loc_origin_from_rank.z +
       target_gl_loc_origin_from_rank.t)%2;// with the current setup, this should always be 0

       if(target_loc_origin_parity) printf("Problems\n");
       */
    // building sublattice duplicate, target_conf
    lnh_su3_soa* target_su3_soa = lnh_conf; 

    int tg_loc_x,tg_loc_y,tg_loc_z,tg_loc_t,dir; //target-loc coordinates
    // and link direction
    // Copying all relevant links from the sublattice to the global lattice
    for(dir =0; dir < 4; dir++)
        for(tg_loc_t=0;tg_loc_t<LOC_NT; tg_loc_t++)
            for(tg_loc_z=0;tg_loc_z<LOC_NZ; tg_loc_z++)
                for(tg_loc_y=0;tg_loc_y<LOC_NY; tg_loc_y++)
                    for(tg_loc_x=0;tg_loc_x<LOC_NX; tg_loc_x++){
    
                        int tg_lnh_x,tg_lnh_y,tg_lnh_z,tg_lnh_t; //target-lnh coordinates
                        tg_lnh_x = tg_loc_x + X_HALO;
                        tg_lnh_y = tg_loc_y + Y_HALO;
                        tg_lnh_z = tg_loc_z + Z_HALO;
                        tg_lnh_t = tg_loc_t + T_HALO;

                        //        int gtsp; // global target site parity
                        int tsprlo ; // target site parity respect (to his) local origin;

                        int target_gl_snum = target_lnh_to_gl_snum(tg_lnh_x, tg_lnh_y, tg_lnh_z, tg_lnh_t, target_gl_loc_origin_from_rank);
                        int target_lnh_snum = snum_acc(tg_lnh_x, tg_lnh_y, tg_lnh_z, tg_lnh_t);

                        tsprlo = (X_HALO+Y_HALO+Z_HALO+T_HALO+ tg_lnh_t+tg_lnh_z+tg_lnh_y+tg_lnh_x)%2;

                        su3 aux = su3_from_su3_soa(lnh_conf,dir,tsprlo,target_lnh_snum);

                        su3_into_global_su3_soa(aux,gl_soa_conf,dir,tsprlo,target_gl_snum);

                    }

}




#endif
