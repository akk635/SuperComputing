/**
 * Finalization step - write results and other computational vectors to files
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#ifndef FINALIZATION_H_
#define FINALIZATION_H_

void finalization( char* file_in, char* out_prefix, int total_iters, double residual_ratio,
                   int nintci, int nintcf, int points_count, int** points, int* elems, double* var,
                   double* cgup, double* su, int *local_global_index,int local_int_cells,
                   int elemcount, int writing_proc );

#endif /* FINALIZATION_H_ */

