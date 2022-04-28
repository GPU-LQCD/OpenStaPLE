#include "../src/RationalApprox/rationalapprox.h"
#include <stdio.h>
#include <stdlib.h>

int verbosity_lv= 0; 

int main(int argc, char **argv){

    if(argc<3){
        printf("usage: %s filename_approx x_where_to_evaluate\n",argv[0] );
        return 1;
    }

    double x = atof(argv[2]);
    RationalApprox in;
    rationalapprox_read_custom_nomefile(&in,argv[1]);

    printf("Value is: %.18e\n",rational_approx_evaluate(&in, x));

    return 0;

    






}




