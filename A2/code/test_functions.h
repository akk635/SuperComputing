/**
 * Functions to test the data distribution and communication lists creation algorithms
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#ifndef TEST_FUNCTIONS_H_
#define TEST_FUNCTIONS_H_
#include "util_write_files.h"
#include <mpi.h>

int test_distribution(char *file_in, char *file_vtk_out, int *local_global_index,
                      int **global_local_index, int nintci, int nintcf, int points_count,
                      int **points, int *elems, int local_int_cells, double *cgup, int elemcount);

int test_communication(char *file_in, char *file_vtk_out, int *local_global_index,
                       int local_num_elems, int neighbors_count, int* send_count, int** send_list,
                       int* recv_count, int** recv_list);

#endif /* TEST_FUNCTIONS_H_ */

