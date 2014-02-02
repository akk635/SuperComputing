/*
 * stats.cpp
 *
 *  Created on: Feb 1, 2014
 *      Author: karthik
 */
#include <stdio.h>

int main(int argc, char *argv[]) {
	if (argc < 3) {
		printf("Usage: %s input_file output_file\n", argv[0]);
		return 1;
	}

	FILE * data_fp, *stats_fp;

	// file name pointers
	char* data_file;
	char* stats_file;

	int procsCounter = 7;
	int msr = 3;
	int results = 8;
	//  nprocs counters and measurements
	double read[procsCounter][results];
	data_file = argv[1];
	stats_file = argv[2];

	data_fp = fopen(data_file, "r");
	stats_fp = fopen(stats_file, "a");
	fprintf(stats_fp, "\n");

	for (int i = 0; i < procsCounter; i++) {
		for (int j = 0; j < msr; j++) {
			fscanf(data_fp, "%lf", &(read[i][j]));
		}
	}

	fprintf(stats_fp,
			"Initialization, Computation, Finalization, OverallTime, OverallSpeedUp, InitSpeedUp, ComputatoinSpeepUp, FinalisationSpeedUp \n");

	fclose(data_fp);

	for (int i = 0; i < procsCounter; i++) {
		read[i][msr] = 0;
		for (int j = 0; j < msr; j++ ){
			fprintf(stats_fp, "%d ," , (int) read[i][j] );
			read[i][msr]+=read[i][j];
		}
		int j = msr; //OverallTime
		fprintf(stats_fp, "%d ," , (int) read[i][j] );
		j++; //OverallSPeedup
		if(i == 0){
			read[i][j] = 1;
			fprintf(stats_fp, "%lf ," , read[i][j] );
		}else{
			read[i][j] = read[0][msr]/read[i][msr]; // T(1)/T(p)
			fprintf(stats_fp, "%lf ," , read[i][j] );
		}
		j++; //InitSpeedUp, Compute and FInal
		int k = 0;
		while(j < results){
			if(i == 0){
				read[i][j] = 1;
				if(j == (results -1)){
					fprintf(stats_fp, "%lf \n" , read[i][j] );
				} else {
					fprintf(stats_fp, "%lf," , read[i][j] );
				}

			}else{
				read[i][j] = read[0][k]/read[i][k]; // T(1)/T(p)
				if(j == (results -1)){
					fprintf(stats_fp, "%lf \n" , read[i][j] );
				} else {
					fprintf(stats_fp, "%lf," , read[i][j] );
				}
			}
			j++; k++;
		}
	}

	fclose(stats_fp);
	return 0;
}

