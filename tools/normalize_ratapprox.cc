#include "./rationalapprox.h"
#include <stdio.h>

int verbosity_lv = 5 ; 

int main(int argc, char **argv){

    if(argc<2){
        printf("usage: %s filename_approx_to_normalize\n",argv[0] );
        return 1;
    }

    RationalApprox in,out;
    rationalapprox_read_custom_nomefile(&in,argv[1]);
    renormalize_rational_approximation(&in,&out);

    char* outnomefile = rational_approx_filename(out.approx_order, out.exponent_num, out.exponent_den, out.lambda_min);

    rationalapprox_save(outnomefile, &out);

    return 0;

    






}




