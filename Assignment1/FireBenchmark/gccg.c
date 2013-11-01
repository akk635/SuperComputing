#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <papi.h>
#include "xread.h"
#include "xwrite.h"
#include "binread.h"
#define NUMEVENTS 4
#define  NUM_FPEVENTS 1

int main(int argc, char *argv[])
{
	if( argc < 4 )
	{
		printf("Usage: %s format_file input_file output_file\n", argv[0]);
		return EXIT_FAILURE;
	}

	// For checking the library initialisation
	int retval;

	//EventSet for L2 & L3 cache misses and accesses
	int EventSet = PAPI_NULL;
	int EventSet1 = PAPI_NULL;

	//Eventset for calculating the mflops
	// int EventSet1[NUM_FEVENTS] = { PAPI_NULL };

	// Data pointer for getting the cpu info
	const PAPI_hw_info_t * hwinfo = NULL;

	PAPI_mh_info_t mem_hrch;

	// Initialising the library
	retval = PAPI_library_init( PAPI_VER_CURRENT );

	if (retval != PAPI_VER_CURRENT) {
		printf( "Initialisation of Papi failed \n" );
		exit(1);
	}

	if ( (hwinfo = PAPI_get_hardware_info()) == NULL){
		printf( "Unable to access hw info \n" );
		return 1;
	}

	/* Accessing the cpus per node, threads per core, memory, frequency */
	printf( "No. of cpus in one node : %d \n", hwinfo ->ncpu );
	printf( "Threads per core : %d \n", hwinfo ->threads );
	printf( "No. of cores per socket : %d \n", hwinfo ->cores );
	printf( "No. of sockets : %d \n", hwinfo ->sockets );
	printf( "Total CPUS in the entire system : %d \n",  hwinfo ->totalcpus );
	mem_hrch = hwinfo->mem_hierarchy;
	printf( "No. of cache levels : %d \n", mem_hrch.levels );
	int L1_size = mem_hrch.level[0].cache[0].size;
	printf( "Size of L1 cache : %d \n", L1_size );
	int L2_size = mem_hrch.level[1].cache[1].size;
	printf( "Size of L2 cache : %d \n", L2_size );

	/* Variables for reading counters of EventSet*/
	long long eventValues[ NUMEVENTS ] = {0};
	long long eventFpValue[ NUM_FPEVENTS ] = {0};

	char *format = argv[1];
	char *file_in = argv[2];
	char *file_out = argv[3];

	char delim[] = ".";
	char *cp = (char *) malloc( sizeof( char ) * 10 );
	cp = strcpy( cp, file_in );
	char *token = malloc( sizeof(char) * 10 );
	token = strtok( cp, delim );
	char * res_file = malloc( sizeof( char ) * 20 );
	res_file = strcpy( res_file, file_out );
	res_file = strcat( res_file, token );
	free (cp);

	FILE * res_fp = fopen( strcat( res_file, "_psdats.dat") , "w" );
	int status = 0;
	free( res_file );

	/** internal cells start and end index*/
	int nintci, nintcf;
	/** external cells start and end index. The external cells are only ghost cells. They are accessed only through internal cells*/
	int nextci, nextcf;
	/** link cell-to-cell array. Stores topology information*/
	int **lcc;
	/** red-black colouring of the cells*/
	int *nboard;

	/** boundary coefficients for each volume cell */
	double *bs, *be, *bn, *bw, *bl, *bh, *bp, *su;

	// Parameters for measuring the time
	long long startusec, endusec;

	//Initializes the library again
	//PAPI_flops(&real_time, &proc_time, &flpins, &mflops) ;

	// Creating the eventSets
	/*if ( PAPI_create_eventset( &EventSet ) != PAPI_OK ) {
		printf( "Problem in create eventset \n" );
		exit(1);
	}*/

	//Create the Flops eventSet
	if ( PAPI_create_eventset( &EventSet1 ) != PAPI_OK ) {
		printf( "Problem in creating the flops eventset \n" );
		exit(1);
	}

	//int EventCode[ NUMEVENTS ] = {  PAPI_L2_TCM, PAPI_L2_TCA, PAPI_L3_TCM, PAPI_L3_TCA };
	int EventFpCode[ NUM_FPEVENTS ] = { PAPI_FP_OPS };
	//Adding events to the eventset
	/*if( PAPI_add_events( EventSet, EventCode, 4 ) != PAPI_OK ){
		printf( "Problem in adding events \n" );
		exit( 1 );
	}*/
	if( PAPI_add_events( EventSet1, EventFpCode, 1 ) != PAPI_OK ){
		printf( "Problem in adding the flops event \n" );
		exit( 1 );
	}
	printf( "Success in adding events \n" );

	//Adding events to the eventset1
	/*if (( PAPI_add_event(EventSet1[0], PAPI_TOT_INS) != PAPI_OK ) &&
	   ( PAPI_add_event(EventSet1[1], PAPI_TOT_CYC) != PAPI_OK ) ){
			exit(1);
	}*/

	// Start the eventset counters
	// PAPI_start( EventSet );
	PAPI_start( EventSet1 );
	/*PAPI_start (EventSet[1]);
	PAPI_start (EventSet[2]);
	PAPI_start (EventSet[3]);*/

	startusec = PAPI_get_real_usec();

	/* initialization  */
	// read-in the input file
	int f_status;
	if( !strcmp( format, "bin" ) ){
		f_status = read_bin_formatted(file_in, &nintci, &nintcf, &nextci, &nextcf, &lcc,
					&bs, &be, &bn, &bw, &bl, &bh, &bp, &su, &nboard);
	} else if( !strcmp( format, "txt" ) ){
		f_status = read_formatted(file_in, &nintci, &nintcf, &nextci, &nextcf, &lcc,
					&bs, &be, &bn, &bw, &bl, &bh, &bp, &su, &nboard);
	}

	if (f_status != 0)
	{
		printf( "failed to initialize data! \n" );
		return EXIT_FAILURE;
	}

	// allocate arrays used in gccg
	int nomax = 3;
	/** the reference residual*/
	double resref = 0.0;
	/** the ratio between the reference and the current residual*/
	double ratio;

	/** array storing residuals */
	double* resvec = (double *) calloc(sizeof(double), (nintcf + 1));
	/** the variation vector -> keeps the result in the end */
	double* var = (double *) calloc(sizeof(double), (nextcf + 1));

	/** the computation vectors */
	double* direc1 = (double *) calloc(sizeof(double), (nextcf + 1));
	double* direc2 = (double *) calloc(sizeof(double), (nextcf + 1));

	/** additional vectors */
	double* cgup = (double *) calloc(sizeof(double), (nextcf + 1));
	double* oc = (double *) calloc(sizeof(double), (nintcf + 1));
	double* cnorm = (double *) calloc(sizeof(double), (nintcf + 1));
	double* adxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
	double* adxor2 = (double *) calloc(sizeof(double), (nintcf + 1));
	double* dxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
	double* dxor2 = (double *) calloc(sizeof(double), (nintcf + 1));

	// initialize the reference residual
	for (int nc = nintci; nc <= nintcf; nc++)
	{
		resvec[nc] = su[nc];
		resref = resref + resvec[nc] * resvec[nc];
	}
	resref = sqrt(resref);
	if (resref < 1.0e-15)
	{
		printf("i/o - error: residue sum less than 1.e-15 - %lf\n", resref);
		return EXIT_FAILURE;
	}

	// initialize the arrays
	for (int nc = 0; nc <= 10; nc++)
	{
		oc[nc] = 0.0;
		cnorm[nc] = 1.0;
	}

	for (int nc = nintci; nc <= nintcf; nc++)
	{
		cgup[nc] = 0.0;
		var[nc] = 0.0;
	}

	for (int nc = nextci; nc <= nextcf; nc++)
	{
		var[nc] = 0.0;
		cgup[nc] = 0.0;
		direc1[nc] = 0.0;
		bs[nc] = 0.0;
		be[nc] = 0.0;
		bn[nc] = 0.0;
		bw[nc] = 0.0;
		bl[nc] = 0.0;
		bh[nc] = 0.0;
	}

	for (int nc = nintci; nc <= nintcf; nc++)
		cgup[nc] = 1.0 / bp[nc];

	int if1 = 0;
	int if2 = 0;
	int iter = 1;
	int nor = 1;
	int nor1 = nor - 1;
	/* finished initalization */

	endusec = PAPI_get_real_usec();

	//Read the eventSet counters
	//PAPI_read( EventSet, eventValues );
	PAPI_read( EventSet1, eventFpValue );

	fprintf( res_fp, "Execution time in microseconds for the initialisation: %lld \n",endusec-startusec );
	fprintf( res_fp, "Initialisation.... \n" );
	/*fprintf( res_fp, "INPUT \t PAPI_L2_TCM \t %lld \n", eventValues[0] );
	fprintf( res_fp, "INPUT \t PAPI_L2_TCA \t %lld \n", eventValues[1] );
	fprintf( res_fp, "INPUT \t PAPI_L3_TCM \t %lld \n", eventValues[2] );
	fprintf( res_fp, "INPUT \t PAPI_L3_TCA \t %lld \n", eventValues[3] );*/
	fprintf( res_fp, "INPUT \t PAPI_FP_OPS \t %lld \n", eventFpValue[0] );
	/*PAPI_read ( EventSet1[0], eventFValues+0 );
	PAPI_read ( EventSet1[1], eventFValues+1 );
	printf ("%lld Floating point instructions executed \n", eventFValues[0] );
	printf ("%lld Total cycles \n", eventFValues[1] );*/

	//Resetting the event counters
	//PAPI_reset( EventSet );
	PAPI_reset( EventSet1 );
	/*PAPI_reset ( EventSet[1] );
	PAPI_reset ( EventSet[2] );
	PAPI_reset ( EventSet[3] );*/
	/*PAPI_reset ( EventSet1[0] );
	PAPI_reset ( EventSet1[1] );*/
	
	fprintf ( res_fp, "Starting with the computation part \n" );
	startusec = PAPI_get_real_usec();

	/* start computation loop */
	while (iter < 10000)
	{
		/* start phase 1 */

		// update the old values of direc
		for (int nc = nintci; nc <= nintcf; nc++)
		{
			direc1[nc] = direc1[nc] + resvec[nc] * cgup[nc];
		}

		// compute new guess (approximation) for direc
		for (int nc = nintci; nc <= nintcf; nc++)
		{
			direc2[nc] = bp[nc] * direc1[nc] - bs[nc] * direc1[lcc[0][nc]]
					- bw[nc] * direc1[lcc[3][nc]] - bl[nc] * direc1[lcc[4][nc]]
					- bn[nc] * direc1[lcc[2][nc]] - be[nc] * direc1[lcc[1][nc]]
					- bh[nc] * direc1[lcc[5][nc]];
		} /* end phase 1 */

		/*  start phase 2 */
		// execute normalization steps
		double oc1, oc2, occ;
		if (nor1 == 1)
		{
			oc1 = 0;
			occ = 0;
			for (int nc = nintci; nc <= nintcf; nc++)
			{
				occ = occ + adxor1[nc] * direc2[nc];
			}
			oc1 = occ / cnorm[1];
			for (int nc = nintci; nc <= nintcf; nc++)
			{
				direc2[nc] = direc2[nc] - oc1 * adxor1[nc];
				direc1[nc] = direc1[nc] - oc1 * dxor1[nc];
			}
			if1++;

		}
		else if (nor1 == 2)
		{
			oc1 = 0;
			occ = 0;
			for (int nc = nintci; nc <= nintcf; nc++)
				occ = occ + adxor1[nc] * direc2[nc];

			oc1 = occ / cnorm[1];
			oc2 = 0;
			occ = 0;
			for (int nc = nintci; nc <= nintcf; nc++)
				occ = occ + adxor2[nc] * direc2[nc];

			oc2 = occ / cnorm[2];
			for (int nc = nintci; nc <= nintcf; nc++)
			{
				direc2[nc] = direc2[nc] - oc1 * adxor1[nc] - oc2 * adxor2[nc];
				direc1[nc] = direc1[nc] - oc1 * dxor1[nc] - oc2 * dxor2[nc];
			}

			if2++;
		}

		cnorm[nor] = 0;
		double omega = 0;

		// compute the new residual
		for (int nc = nintci; nc <= nintcf; nc++)
		{
			cnorm[nor] = cnorm[nor] + direc2[nc] * direc2[nc];
			omega = omega + resvec[nc] * direc2[nc];
		}
		omega = omega / cnorm[nor];

		double resnew = 0.0;
		for (int nc = nintci; nc <= nintcf; nc++)
		{
			var[nc] = var[nc] + omega * direc1[nc];
			resvec[nc] = resvec[nc] - omega * direc2[nc];
			resnew = resnew + resvec[nc] * resvec[nc];
		}
		resnew = sqrt(resnew);
		ratio = resnew / resref;

		// exit on no improvements of residual
		if (ratio <= 1.0e-10)
			break;

		iter++;

		// prepare additional arrays for the next iteration step
		if (nor == nomax)
			nor = 1;
		else
		{
			if (nor == 1)
			{
				for (int nc = nintci; nc <= nintcf; nc++)
				{
					dxor1[nc] = direc1[nc];
					adxor1[nc] = direc2[nc];
				}

			} else if (nor == 2)
			{
				for (int nc = nintci; nc <= nintcf; nc++)
				{
					dxor2[nc] = direc1[nc];
					adxor2[nc] = direc2[nc];
				}
			}
			nor++;
		}
		nor1 = nor - 1;

	}/* end phase 2 */

	/* finished computation loop */
	endusec = PAPI_get_real_usec();

	//Read the eventSet counters
	//PAPI_stop( EventSet, eventValues );
	PAPI_stop( EventSet1, eventFpValue );
	/*PAPI_stop ( EventSet[1], eventValues+1 );
	PAPI_stop ( EventSet[2], eventValues+2 );
	PAPI_stop ( EventSet[3], eventValues+3 );*/
	fprintf( res_fp, "Execution time in microseconds for the computation : %lld \n",endusec-startusec);
	/*fprintf( res_fp, "CALC \t PAPI_L2_TCM \t %lld \n", eventValues[0] );
	fprintf( res_fp, "CALC \t PAPI_L2_TCA \t %lld \n", eventValues[1] );
	fprintf( res_fp, "CALC \t PAPI_L3_TCM \t %lld \n", eventValues[2] );
	fprintf( res_fp, "CALC \t PAPI_L3_TCA \t %lld \n", eventValues[3] );*/
	fprintf( res_fp, "CALC \t PAPI_FP_OPS \t %lld \n", eventFpValue[0] );

/*	PAPI_read ( EventSet1[0], eventFValues+0 );
	PAPI_read ( EventSet1[1], eventFValues+1 );
	printf ("%lld Floating point instructions executed \n", eventFValues[0] );
	printf ("%lld Total cycles \n", eventFValues[1] );*/

	/* write output file  */
	if ( write_result(file_in, file_out, nintci, nintcf, var, iter, ratio) != 0 )
		printf("error when trying to write to file %s\n", file_out);

	/* Free all the dynamically allocated memory */
	free(direc2); free(direc1); free(dxor2); free(dxor1); free(adxor2); free(adxor1);
	free(cnorm); free(oc); free(var); free(cgup); free(resvec); free(su); free(bp);
	free(bh); free(bl); free(bw); free(bn); free(be); free(bs);

	printf("Simulation completed successfully!\n");

	fclose( res_fp );
	return EXIT_SUCCESS;
}
