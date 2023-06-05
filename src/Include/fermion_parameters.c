#ifndef FERMION_PARAMETERS_C_
#define FERMION_PARAMETERS_C_

#include "../DbgTools/dbgtools.h"
#include "../Include/debug.h"
#include "../Include/inverter_tricks.h"
#include "../Meas/magnetic_susceptibility_utilities.h"
#include "../Mpi/multidev.h"
#include "../OpenAcc/alloc_vars.h"
#include "../OpenAcc/backfield.h"
#include "../OpenAcc/geometry.h"
#include "../OpenAcc/md_integrator.h"
#include "../OpenAcc/md_parameters.h"
#include "../OpenAcc/sp_alloc_vars.h"
#include "../OpenAcc/sp_backfield.h"
#include "../RationalApprox/rationalapprox.h"
#include "./fermion_parameters.h"
#include "./montecarlo_parameters.h"
#include "../OpenAcc/alloc_settings.h"
#include <math.h>
#include <string.h>
#include <stdio.h>


#define ALIGN 128
#define acc_twopi 2*3.14159265358979323846

int max_ps;
int totalMdShifts;
ferm_param *fermions_parameters; // set in init.c, from input file

int init_ferm_params(ferm_param *fermion_settings){
	alloc_info.maxNeededShifts = 0;
	alloc_info.maxApproxOrder = 0;

	int errorstatus = 0;


	printf("Initializing fermions...\n");

	alloc_info.NPS_tot = 0;
	max_ps = 0;

	// calculation of alloc_info.NPS_tot, max_ps,index_of_the_first_ps; 
	for(int i=0;i<alloc_info.NDiffFlavs;i++){
		fermion_settings[i].printed_bf_dbg_info = 0;
		// compute the total number of ps
		alloc_info.NPS_tot += fermion_settings[i].number_of_ps;
		// compute the max number of ps among the various flavs
		if(fermion_settings[i].number_of_ps>=max_ps) max_ps = fermion_settings[i].number_of_ps;
		// determine the offset (where does the ps of the flavour i starts?)
		if(i==0){
			fermion_settings[i].index_of_the_first_ps=0;
		}else{
			fermion_settings[i].index_of_the_first_ps = fermion_settings[i-1].index_of_the_first_ps + fermion_settings[i-1].number_of_ps;
		}
	}

	printf("alloc_info.NPS_tot = %d \n",alloc_info.NPS_tot);
	printf("max_ps = %d \n",max_ps);

	totalMdShifts = 0;
	// Rational Approximation related stuff
	for(int i=0;i<alloc_info.NDiffFlavs;i++){
		ferm_param *quark = &fermion_settings[i];
		quark->approx_fi_mother.exponent_num =  +quark->degeneracy;
		quark->approx_md_mother.exponent_num =  -quark->degeneracy;
		quark->approx_li_mother.exponent_num =  -quark->degeneracy;

		quark->approx_fi_mother.exponent_den =   quark->number_of_ps*8;
		quark->approx_md_mother.exponent_den =   quark->number_of_ps*4;
		quark->approx_li_mother.exponent_den =   quark->number_of_ps*4;

		quark->approx_fi_mother.lambda_min = quark->ferm_mass*quark->ferm_mass/md_parameters.expected_max_eigenvalue;
		quark->approx_md_mother.lambda_min = quark->ferm_mass*quark->ferm_mass/md_parameters.expected_max_eigenvalue;
		quark->approx_li_mother.lambda_min = quark->ferm_mass*quark->ferm_mass/md_parameters.expected_max_eigenvalue;

		quark->approx_fi_mother.lambda_max =  1.0;
		quark->approx_md_mother.lambda_max =  1.0;
		quark->approx_li_mother.lambda_max =  1.0;

		quark->approx_fi_mother.error =  md_parameters.residue_metro/
			pow(md_parameters.expected_max_eigenvalue,
					(double) quark->approx_fi_mother.exponent_num/
					quark->approx_fi_mother.exponent_den );
		quark->approx_md_mother.error =  md_parameters.residue_md/
			pow(md_parameters.expected_max_eigenvalue,
					(double) quark->approx_md_mother.exponent_num/
					quark->approx_md_mother.exponent_den );

		quark->approx_li_mother.error =  md_parameters.residue_metro/
			pow(md_parameters.expected_max_eigenvalue,
					(double) quark->approx_li_mother.exponent_num/
					quark->approx_li_mother.exponent_den );

		quark->approx_fi_mother.gmp_remez_precision = 100;
		quark->approx_md_mother.gmp_remez_precision = 100;
		quark->approx_li_mother.gmp_remez_precision = 100;

		// copy everything also in the daughter approxs
		quark->approx_fi.exponent_num =   quark->approx_fi_mother.exponent_num;
		quark->approx_md.exponent_num =   quark->approx_md_mother.exponent_num;
		quark->approx_li.exponent_num =   quark->approx_li_mother.exponent_num;
		quark->approx_fi.exponent_den =   quark->approx_fi_mother.exponent_den;
		quark->approx_md.exponent_den =   quark->approx_md_mother.exponent_den;
		quark->approx_li.exponent_den =   quark->approx_li_mother.exponent_den;
		quark->approx_fi.approx_order =   quark->approx_fi_mother.approx_order;
		quark->approx_md.approx_order =   quark->approx_md_mother.approx_order;
		quark->approx_li.approx_order =   quark->approx_li_mother.approx_order;
		quark->approx_fi.gmp_remez_precision =
			quark->approx_fi_mother.gmp_remez_precision;
		quark->approx_md.gmp_remez_precision =   
			quark->approx_md_mother.gmp_remez_precision;
		quark->approx_li.gmp_remez_precision =   
			quark->approx_li_mother.gmp_remez_precision;

		// read the rat approxs from the files
		int temp_errorstatus;
		// first inversion
		temp_errorstatus = rationalapprox_read(&(quark->approx_fi_mother));
		if(temp_errorstatus) 
			temp_errorstatus = rat_approx_file_or_script_create(&(quark->approx_fi_mother));
		errorstatus += temp_errorstatus;

		// molecular dynamics
		temp_errorstatus = rationalapprox_read(&(quark->approx_md_mother));
		if(temp_errorstatus) 
			temp_errorstatus = rat_approx_file_or_script_create(&(quark->approx_md_mother));
		errorstatus += temp_errorstatus;

		// last inversion
		temp_errorstatus = rationalapprox_read(&(quark->approx_li_mother));
		if(temp_errorstatus) 
			temp_errorstatus = rat_approx_file_or_script_create(&(quark->approx_li_mother));
		errorstatus += temp_errorstatus;

		// needed to reuse the results from the inversions
		quark->index_of_the_first_shift = totalMdShifts;
		totalMdShifts += quark->approx_md_mother.approx_order * quark->number_of_ps;
     
		if(alloc_info.maxNeededShifts < quark->approx_fi_mother.approx_order)
			alloc_info.maxNeededShifts = quark->approx_fi_mother.approx_order;
		if(alloc_info.maxNeededShifts < quark->approx_md_mother.approx_order)
			alloc_info.maxNeededShifts = quark->approx_md_mother.approx_order;
		if(alloc_info.maxNeededShifts < quark->approx_li_mother.approx_order)
			alloc_info.maxNeededShifts = quark->approx_li_mother.approx_order;

	}

	alloc_info.maxApproxOrder = alloc_info.maxNeededShifts;
	if(1 == md_parameters.recycleInvsForce && alloc_info.maxNeededShifts < totalMdShifts*2)
		alloc_info.maxNeededShifts = totalMdShifts*2;

	if(0==devinfo.myrank && verbosity_lv > 3){
		printf("MaxApproxOrder found/#define'd: %d / %d\n",alloc_info.maxApproxOrder, MAX_APPROX_ORDER);
		printf("MaxNeededShifts: %d ",alloc_info.maxNeededShifts);
		if(1 == md_parameters.recycleInvsForce)
			printf("(md_parameters.recycleInvsForce = 1)");
		printf("\n");
		printf("totalMdShifts: %d \n",totalMdShifts);
	}

	return errorstatus;

}


int rat_approx_file_or_script_create(RationalApprox* rational_approx){

	int error_status = 1;

	char * nomefile = rational_approx_filename(rational_approx->error,rational_approx->exponent_num,rational_approx->exponent_den,rational_approx->lambda_min);

#ifdef MULTIDEVICE
	printf("MPI%02d - Some error happened in reading %s ...\n", devinfo.myrank,nomefile );
	if(0==devinfo.myrank){
		FILE * bash_repair_commands = fopen("genappfiles.sh","a");

		printf("You may want to generate a rational approximation file using the tool \'rgen\' (look in the tools directory). Please try\n");
		printf("./rgen %e %d %d %e\n", rational_approx->error, 
					 rational_approx->exponent_num, rational_approx->exponent_den, 
					 rational_approx->lambda_min);
		printf("(see and modify \"genappfiles.sh\", check for doublers)\n");
		printf("(Or give command \n bash <(sort genappfiles.sh | uniq).\n");
		fprintf(bash_repair_commands,
						"echo \'./rgen %e %d %d %e >> rat_app_gen_log.txt &\'\n",
						rational_approx->error, rational_approx->exponent_num,
						rational_approx->exponent_den, rational_approx->lambda_min);
		fprintf(bash_repair_commands,"./rgen %e %d %d %e >> rat_app_gen_log.txt &\n",
						rational_approx->error, rational_approx->exponent_num,
						rational_approx->exponent_den, rational_approx->lambda_min);
		fclose(bash_repair_commands);
	}
	MPI_Barrier(devinfo.mpi_comm);
	error_status = 1; // error cannot be corrected here.
#else
	char command[100];
	sprintf(command,
					"./rgen %e %d %d %e\n",
					rational_approx->error, rational_approx->exponent_num,
					rational_approx->exponent_den, 
					rational_approx->lambda_min);

	printf("Creating (and caching) file %s, wait ...\n", nomefile);
	int status=system(command);

	error_status = rationalapprox_read_custom_nomefile(rational_approx,nomefile);
	free(nomefile);

#endif
	return error_status;

}


void init_all_u1_phases(bf_param bfpars, ferm_param *fpar  )
{

	for(int i=0;i<alloc_info.NDiffFlavs;i++){
		fpar[i].phases = &u1_back_phases[i*8];
		fpar[i].phases_f = &u1_back_phases_f[i*8];
		init_fermion_backfield(bfpars,&(fpar[i]));

		// printing debug info
		if(debug_settings.print_bfield_dbginfo){
			char tempname[50];                           
			// phases
			sprintf(tempname,"backfield_%s_c%d",fpar[i].name,fpar[i].printed_bf_dbg_info);
#ifdef MULTIDEVICE      
			strcat(tempname,devinfo.myrankstr);          
#endif
			print_double_soa(fpar[i].phases,tempname);   

			// plaquettes
			sprintf(tempname,"abelian_plq_%s_c%d",fpar[i].name,
							fpar[i].printed_bf_dbg_info);
#ifdef MULTIDEVICE      
			strcat(tempname,devinfo.myrankstr);          
#endif

			print_all_abelian_plaquettes(fpar[i].phases,tempname);
			fpar->printed_bf_dbg_info++; 
		}
    
        
		fpar[i].mag_re = &mag_obs_re[i*8];
		fpar[i].mag_im = &mag_obs_im[i*8];

		idphase_dbz(fpar[i].mag_re,fpar[i].mag_im,&fpar[i]);
    
	}


}

void init_fermion_backfield(bf_param bf_pars, ferm_param *fermion_parameters)
{

	if(verbosity_lv > 2 && 0 == devinfo.myrank ) { 
		printf("Generating external field (containing staggered phases) ");
		printf("for flavour %s\n",fermion_parameters->name);
	}


	calc_u1_phases(fermion_parameters->phases, bf_pars, 
								 fermion_parameters->ferm_im_chem_pot, fermion_parameters->ferm_charge);



	if(inverter_tricks.useMixedPrecision || md_parameters.singlePrecMD)
		calc_u1_phases_f(fermion_parameters->phases_f, bf_pars, 
										 fermion_parameters->ferm_im_chem_pot, fermion_parameters->ferm_charge);




}

#endif
