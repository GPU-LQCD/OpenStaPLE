#ifndef RATIONAL_APPROX_C_
#define RATIONAL_APPROX_C_

#include "./rationalapprox.h"
#include <stdlib.h>

#define CHECKREAD(expr,should_read) \
{int read = expr ;if(read != should_read) { \
        printf("%s:%d, Error, not read expected number of entries : %d read vs %d should_read\n.", __FILE__, __LINE__ , read,should_read);\
        exit(1);}}\



extern int verbosity_lv;

char* rational_approx_filename_old(int approx_order, int exponent_num, int exponent_den, double lambda_min)
{

    char Cn1[3];
    char Cy1[3];
    char Cz1[3];
    char mloglmin[5];

    sprintf(Cn1, "%d", approx_order);  // order
    sprintf(Cy1, "%d", exponent_num); // num
    sprintf(Cz1, "%d", exponent_den); // den
    sprintf(mloglmin, "%1.1f", -log(lambda_min)/log(10.0)); // lambda min

    char * nomefile = (char*)malloc(50*sizeof(char));
    strcpy(nomefile,"approx_");
    strcat(nomefile,Cy1);
    strcat(nomefile,"_over_");
    strcat(nomefile,Cz1);
    strcat(nomefile,"_order_");
    strcat(nomefile,Cn1);
    strcat(nomefile,"_mloglm_");
    strcat(nomefile,mloglmin);
    strcat(nomefile,".REMEZ");

    return nomefile;

}

char* rational_approx_filename(double error, int exponent_num, int exponent_den, double lambda_min)
{

    char mlogerr[5];
    char Cy1[3];
    char Cz1[3];
    char mloglmin[5];

    sprintf(mlogerr, "%1.1lf", -log(error)/log(10.0));  // error
    sprintf(Cy1, "%d", exponent_num); // num
    sprintf(Cz1, "%d", exponent_den); // den
    sprintf(mloglmin, "%1.1f", -log(lambda_min)/log(10.0)); // lambda min

    char * nomefile = (char*)malloc(50*sizeof(char));
    strcpy(nomefile,"approx_");
    strcat(nomefile,Cy1);
    strcat(nomefile,"_over_");
    strcat(nomefile,Cz1);
    strcat(nomefile,"_mlogerr_");
    strcat(nomefile,mlogerr);
    strcat(nomefile,"_mloglm_");
    strcat(nomefile,mloglmin);
    strcat(nomefile,".REMEZ");

    return nomefile;

}

int rationalapprox_read(RationalApprox* rational_approx)
{

    char * nomefile = rational_approx_filename(rational_approx->error,rational_approx->exponent_num,rational_approx->exponent_den,rational_approx->lambda_min);

    int error = rationalapprox_read_custom_nomefile(rational_approx,nomefile);
    free(nomefile);
    return error;
}


int rationalapprox_read_custom_nomefile(RationalApprox* rational_approx, char* nomefile)
{

    // CALCULATION OF COEFFICIENTS FOR FIRST_INV_APPROX_NORM_COEFF
    FILE *input = fopen(nomefile, "rt");
    printf("%s\n", nomefile );
    if (input == NULL) {
        printf("Rational Approximation File %s not found!\n", nomefile);
        return 1;
    }

    // The partial fraction expansion takes the form 
    // r(x) = norm + sum_{k=1}^{n} res[k] / (x + pole[k])
    CHECKREAD(fscanf(input,"\nApproximation to f(x) = (x)^(%i/%i)\n", &(rational_approx->exponent_num), &(rational_approx->exponent_den)),2) ;
    CHECKREAD(fscanf(input,"Order: %i\n", &(rational_approx->approx_order)),1);
    CHECKREAD(fscanf(input,"Lambda Min: %lf\n", &(rational_approx->lambda_min)),1);
    CHECKREAD(fscanf(input,"Lambda Max: %lf\n", &(rational_approx->lambda_max)),1);
    CHECKREAD(fscanf(input,"GMP Remez Precision: %i\n", &(rational_approx->gmp_remez_precision)),1);
    CHECKREAD(fscanf(input,"Error: %lf\n", &(rational_approx->error)),1);
    CHECKREAD(fscanf(input, "RA_a0 = %lf\n",&(rational_approx->RA_a0)),1);
    for(int i = 0; i < rational_approx->approx_order; i++)
    {
        CHECKREAD(fscanf(input, "RA_a[%i] = %lf, RA_b[%i] = %lf\n", &i, &(rational_approx->RA_a[i]), &i, &(rational_approx->RA_b[i])),4);
    }
    fclose(input);

    if(verbosity_lv > 4){
        printf("RA_a0 = %18.16e\n", rational_approx->RA_a0);
        for(int i = 0; i < rational_approx->approx_order; i++) 
        {
            printf("RA_a[%d] = %18.16e, RA_b[%d] = %18.16e\n", i, rational_approx->RA_a[i], i, rational_approx->RA_b[i]);
        }
    }
    return 0; 
}

void rationalapprox_save(const char* nomefile, RationalApprox* rational_approx){

    FILE *output = fopen(nomefile, "w");
    printf("%s\n", nomefile );
    if (output == NULL) {
        printf("Could not open file %s \n",nomefile );
        exit(-1);
    }

    // The partial fraction expansion takes the form 
    // r(x) = norm + sum_{k=1}^{n} res[k] / (x + pole[k])
    fprintf(output,"\nApproximation to f(x) = (x)^(%i/%i)\n", rational_approx->exponent_num, rational_approx->exponent_den);
    fprintf(output,"Order: %i\n", rational_approx->approx_order);
    fprintf(output,"Lambda Min: %18.16e\n", rational_approx->lambda_min);
    fprintf(output,"Lambda Max: %18.16e\n", rational_approx->lambda_max);
    fprintf(output,"GMP Remez Precision: %i\n", rational_approx->gmp_remez_precision);
    fprintf(output,"Error: %18.16e\n", rational_approx->error);
    fprintf(output, "RA_a0 = %18.16e\n",rational_approx->RA_a0);
    for(int i = 0; i < rational_approx->approx_order; i++)
    {
        fprintf(output,"RA_a[%d] = %18.16e, RA_b[%d] = %18.16e\n", i, rational_approx->RA_a[i], i, rational_approx->RA_b[i]);
    }
    fclose(output);

}

void rescale_rational_approximation(RationalApprox *in, RationalApprox *out, double *minmax){

    double power = (double) in->exponent_num/in->exponent_den;



    out->exponent_num        = in->exponent_num       ;
    out->exponent_den        = in->exponent_den       ;
    out->approx_order        = in->approx_order       ;
    out->gmp_remez_precision = in->gmp_remez_precision;              
    out->error               = in->error              ;

    // HERE THE ASSUMPTION IS THAT in->lambda_max  = 1

    double min =  minmax[0];
    double max =  minmax[1];
    min*=0.95;
    max*=1.05;
    double epsilon=pow(max, power);  
    out->RA_a0               = in->RA_a0       *     epsilon ;
    for(int order = 0; order < in->approx_order; order ++){
        out->RA_a[order] = in->RA_a[order]*max * epsilon;
        out->RA_b[order] = in->RA_b[order]*max ;
    }

    out->lambda_min  = in->lambda_min * max ;
    out->lambda_max  = max ;
    //pray
    if(out->lambda_min > minmax[0]){
        printf("ERROR: mother rational approx does not cover the range!\n");
        printf("out->lambda_max: %e \n", out->lambda_max);
        printf("out->lambda_min: %e , minmax[0]: %e\n", out->lambda_min, minmax[0] );
        printf("Consider modifying your input file, setting\n");
        printf("ExpMaxEigenvalue         %f     # OR LARGER!!\n" , max*1.2 );
        printf("Program will now terminate.\n");


        exit(1);
    }else if(out->lambda_min< (0.3* minmax[0])){
        printf("Warning, the range of your rational approximation is really large\n");
        printf("out->lambda_min: %e , minmax[0]: %e\n", out->lambda_min, minmax[0] );
        printf("You may consider reducing \'ExpMaxEigenvalue\' in the input file to %f.\n",
                max*1.2 );


    }



}


void renormalize_rational_approximation(RationalApprox *in, RationalApprox *out){
    double power = (double) in->exponent_num/in->exponent_den;



    out->exponent_num        = in->exponent_num       ;
    out->exponent_den        = in->exponent_den       ;
    out->approx_order        = in->approx_order       ;
    out->gmp_remez_precision = in->gmp_remez_precision;              
    out->error               = in->error              ;

    // HERE THE ASSUMPTION IS THAT in->lambda_max  = 1

    double rescale_ratio =  1/in->lambda_max;
    double epsilon=pow(rescale_ratio, power);  
    out->RA_a0               = in->RA_a0       *     epsilon ;
    for(int order = 0; order < in->approx_order; order ++){
        out->RA_a[order] = in->RA_a[order]*rescale_ratio * epsilon;
        out->RA_b[order] = in->RA_b[order]*rescale_ratio ;
    }

    out->lambda_min  = in->lambda_min * rescale_ratio ;
    out->lambda_max  = 1 ;

}


double rational_approx_evaluate(RationalApprox *in, double x){

    double res = in->RA_a0;

    int order;
    for(order = 0; order<in->approx_order;++order){

        res += in->RA_a[order]/(x+in->RA_b[order]);

    }
    return res;

}


#endif


