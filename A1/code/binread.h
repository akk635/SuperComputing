/*
 * binread.h
 *
 *  Created on: Nov 1, 2013
 *      Author: karthik
 */

#ifndef BINREAD_H_
#define BINREAD_H_

int read_bin_formatted(char *fileName,
                int *NINTCI, int *NINTCF, int *NEXTCI, int *NEXTCF,
                int ***LCC,
                double **BS, double **BE, double **BN, double **BW, double **BL, double **BH, double **BP, double **SU,
                int **NBOARD
                );

#endif /* BINREAD_H_ */
