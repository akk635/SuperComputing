/**
 * Helper functions for reading from input data file
 *
 * @author E. Xue, V. Petkov
 * @date 22-May-2009, 22-Oct-2012
 */
#ifndef XREAD_H_
#define XREAD_H_
#include <mpi.h>

int read_binary_geo( char *file_name, char* part_type, int *NINTCI, int *NINTCF,
                     int *NEXTCI, int *NEXTCF, int ***LCC, double **BS, double **BE,
                     double **BN, double **BW, double **BL, double **BH, double **BP,
                     double **SU, int* points_count, int*** points, int** elems,
                     int **local_global_index, int *elemcount, int *local_int_cells,
                     int ***global_local_index, int **epart, int **npart, int *objval );

#endif /* XREAD_H_ */


