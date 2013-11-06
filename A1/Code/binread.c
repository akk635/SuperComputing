/*
 * binread.c
 *
 *  Created on: Nov 1, 2013
 *      Author: karthik
 */
#include "binread.h"
#include <stdio.h>
#include <stdlib.h>

int read_bin_formatted(char *fileName,
                int *NINTCI, int *NINTCF, int *NEXTCI, int *NEXTCF,
                int ***LCC,
                double **BS, double **BE, double **BN, double **BW, double **BL, double **BH, double **BP, double **SU,
                int **NBOARD
                )
{
	int i;
	FILE *fp = fopen(fileName, "rb");
	if (fp == NULL)
	{
		printf("Error opening file %s\n", fileName);
		return -1;
	}

	// Reading the volume dimensions
	fread( NINTCI, sizeof(int), 1, fp );
	fread( NINTCF, sizeof(int), 1, fp );
	fread( NEXTCI, sizeof(int), 1, fp );
	fread( NEXTCF, sizeof(int), 1, fp );

	// Allocating topological information in LCC
	if ((*LCC = (int**) malloc(6 * sizeof(int*))) == NULL)
	{
		printf("malloc failed to allocate first dimension of LCC (6)");
		return -1;
	}
	for (i = 0; i < 6; i++)
	{
		if (((*LCC)[i] = (int *) malloc((*NINTCF + 1) * sizeof(int))) == NULL)
		{
			printf("malloc(LCC) failed\n");
			return -1;
		}
	}
	//start reading LCC
	//Note that C array index starts from 0 while Fortran starts from 1!
	for (i = (*NINTCI); i <= *NINTCF; i++)
	{
		fread( &(*LCC)[0][i], sizeof(int), 1, fp );
		fread( &(*LCC)[1][i], sizeof(int), 1, fp );
		fread( &(*LCC)[2][i], sizeof(int), 1, fp );
		fread( &(*LCC)[3][i], sizeof(int), 1, fp );
		fread( &(*LCC)[4][i], sizeof(int), 1, fp );
		fread( &(*LCC)[5][i], sizeof(int), 1, fp );
	}
	// allocate other arrays
	if ((*BS = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL)
	{
		printf("malloc() failed\n");
		return -1;
	}
	if ((*BE = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL)
	{
		printf("malloc() failed\n");
		return -1;
	}
	if ((*BN = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL)
	{
		printf("malloc() failed\n");
		return -1;
	}
	if ((*BW = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL)
	{
		printf("malloc() failed\n");
		return -1;
	}
	if ((*BL = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL)
	{
		printf("malloc() failed\n");
		return -1;
	}
	if ((*BH = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL)
	{
		printf("malloc() failed\n");
		return -1;
	}
	if ((*BP = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL)
	{
		printf("malloc() failed\n");
		return -1;
	}
	if ((*SU = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL)
	{
		printf("malloc() failed\n");
		return -1;
	}

	// Reading the boundary stencils and source values
	for (i = (*NINTCI); i <= *NINTCF; i++)
	{
		fread( &((*BS)[i]), sizeof(double), 1, fp );
		fread( &((*BE)[i]), sizeof(double), 1, fp );
		fread( &((*BN)[i]), sizeof(double), 1, fp );
		fread( &((*BW)[i]), sizeof(double), 1, fp );
		fread( &((*BL)[i]), sizeof(double), 1, fp );
		fread( &((*BH)[i]), sizeof(double), 1, fp );
		fread( &((*BP)[i]), sizeof(double), 1, fp );
		fread( &((*SU)[i]), sizeof(double), 1, fp );
	}

	// Read red-black parameters
	// read board
	if ((*NBOARD = (int *) malloc((*NINTCF + 1) * sizeof(int))) == NULL)
	{
		printf("malloc() failed\n");
		return -1;
	}

	for (i = (*NINTCI); i <= *NINTCF; i++){
		fread( &((*NBOARD)[i]), sizeof(int), 1, fp );
	}

	fclose( fp );
	return 0;
}

