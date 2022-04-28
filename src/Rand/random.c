#ifndef RANDOM_C_
#define RANDOM_C_

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include"./RANDOM/dSFMT.c"
#include "../Include/debug.h"
#include "../OpenAcc/geometry.h"
#include"./random.h"

dsfmt_t dsfmt;

int rng_fake_gl_index;
const int rng_fake_gl_index_max = GL_SIZE * 4 * 18;

extern int verbosity_lv;

// random number generator in (0,1)
double casuale(void)
{ 
    switch(debug_settings.rng_fakeness_level){

        case 0:
            return dsfmt_genrand_open_open(&dsfmt);
            break;
        case 1:
            return (double) rng_fake_gl_index / rng_fake_gl_index_max;
            break;
        case 2:
            return 0.5;
            break;
        default: 
            return dsfmt_genrand_open_open(&dsfmt);
    }
}


// random number initialization
void initrand(unsigned long s)
{
    if(s==0)
    {
        time_t t;
        //      srand((unsigned) time(&t));

        dsfmt_init_gen_rand(&dsfmt, time(NULL));
    }
    else
    {
        //      srand(s);
        dsfmt_init_gen_rand(&dsfmt, s);
    }
}



void initrand_fromfile(const char * filename, unsigned long seed_default){

    FILE * seedfile = fopen(filename,"r");
    int reads;
    int problems_in_read = 0;
    if(!seedfile){
        problems_in_read = 1;
    }else{
        reads = fread(&dsfmt,1,sizeof(dsfmt_t), seedfile);
        if(reads!=sizeof(dsfmt_t))
            problems_in_read = 1;
        fclose(seedfile);
    }
    if(problems_in_read)
        initrand(seed_default);
    if(verbosity_lv > 0){
        if(problems_in_read){
            printf("Problems in reading file %s, ",filename );
            printf("Initialising RNG with seed %lu\n",seed_default);
        }
        else{
            printf("RNG initialised from file %s.\n",filename );

        }
    }
}

void saverand_tofile(const char * filename){

    FILE * seedfile = fopen(filename,"w");

    fwrite(&dsfmt, 1,sizeof(dsfmt_t), seedfile);
    fclose(seedfile);
    if(verbosity_lv > 0)
        printf("RNG status saved in file %s.\n",filename );

}






// 4 parameters for random SU(2) matrix
//void su2_rand(double &p0, double &p1, double &p2, double &p3)
void su2_rand(double *pp)
{ 
    double p=2.0;
    while(p>1.0)
    {
        pp[0]=1.0-2.0*casuale();
        pp[1]=1.0-2.0*casuale();
        pp[2]=1.0-2.0*casuale();
        pp[3]=1.0-2.0*casuale();
        p=sqrt(pp[0]*pp[0]+pp[1]*pp[1]+pp[2]*pp[2]+pp[3]*pp[3]);
    }

    pp[0]/=p;
    pp[1]/=p;
    pp[2]/=p;
    pp[3]/=p;
}





#endif
