/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdlib.h>
#include <stdio.h>
#include "util_read_files.h"
#include "initialization.h"

int initialization(char* file_in, char* part_type, int* nintci, int* nintcf, int* nextci,
                   int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
                   double** bl, double** bh, double** bp, double** su, int* points_count,
                   int*** points, int** elems, double** var, double** cgup, double** oc,
                   double** cnorm, int** local_global_index, int*** global_local_index,
                   int* neighbors_count, int** send_count, int*** send_list, int** recv_count,
                   int*** recv_list, int** epart, int** npart, int* objval, int *elemcount, int *local_int_cells) {
    /********** START INITIALIZATION **********/
    int i = 0;

    // read-in the input file
    int f_status = read_binary_geo( file_in, part_type, &*nintci, &*nintcf, &*nextci, &*nextcf,
                                    &*lcc, &*bs, &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su,
                                    &*points_count, &*points, &*elems, &*local_global_index,
                                    elemcount, local_int_cells, global_local_index, epart, npart,
                                    objval );

    if ( f_status != 0 )
        return f_status;

    *var = (double*) calloc( sizeof(double), *elemcount );
    *cgup = (double*) calloc( sizeof(double), *elemcount );
    *cnorm = (double*) calloc( sizeof(double), *local_int_cells );

    // initialize the arrays
    for ( i = 0; i <= 10; i++ ) {
        ( *cnorm )[i] = 1.0;
    }

    for ( i = 0; i < *local_int_cells; i++ ) {
        ( *var )[i] = 0.0;
    }

    for ( i = 0; i < *local_int_cells; i++ ) {
        ( *cgup )[i] = 1.0 / ( ( *bp )[i] );
    }

    for ( i = *local_int_cells; i < *elemcount; i++ ) {
        ( *var )[i] = 0.0;
        ( *cgup )[i] = 0.0;
        ( *bs )[i] = 0.0;
        ( *be )[i] = 0.0;
        ( *bn )[i] = 0.0;
        ( *bw )[i] = 0.0;
        ( *bh )[i] = 0.0;
        ( *bl )[i] = 0.0;
    }

    return 0;
}

