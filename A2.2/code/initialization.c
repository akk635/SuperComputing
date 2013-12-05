/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "util_read_files.h"
#include "initialization.h"

int contains( int index, int *search_array, int size );

int initialization( char* file_in, char* part_type, int* nintci, int* nintcf, int* nextci,
                    int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
                    double** bl, double** bh, double** bp, double** su, int* points_count,
                    int*** points, int** elems, double** var, double** cgup, double** oc,
                    double** cnorm, int** local_global_index, int*** global_local_index,
                    int* neighbors_count, int** send_count, int*** send_list, int** recv_count,
                    int*** recv_list, int** epart, int** npart, int* objval, int *elemcount,
                    int *local_int_cells ) {
    /********** START INITIALIZATION **********/
    int i = 0;
    int my_rank, nproc;
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &nproc );

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

    int temp_rank = 0;
    // Temporary counts just for the allocation
    int *neighbors = calloc( nproc, sizeof(int) );
    // Provides with the max size of the array
    for ( int i = 0; i < ( *local_int_cells ); i++ ) {
        for ( int j = 0; j < 6; j++ ) {
            // Only appending the internal cells
            if ( ( temp_rank = ( *global_local_index )[( *lcc )[i][j]][0] ) != my_rank ) {
                neighbors[temp_rank] += 1;
            }
        }
    }

    assert( neighbors[my_rank] == 0 );

    // For the actual count
    *recv_count = calloc( nproc, sizeof(int) );
    *send_count = calloc( nproc, sizeof(int) );

    ( *send_list ) = (int **) malloc( nproc * sizeof(int *) );
    // Allocating with the max array size
    for ( int i = 0; i < nproc; i++ ) {
        ( *send_list )[i] = (int *) calloc( neighbors[i], sizeof(int) );
    }

    ( *recv_list ) = (int **) malloc( nproc * sizeof(int *) );

    // Allocating the max array size
    for ( int i = 0; i < nproc; i++ ) {
        ( *recv_list )[i] = (int *) calloc( neighbors[i], sizeof(int) );
    }

    for ( int i = 0; i < ( *local_int_cells ); i++ ) {
        for ( int j = 0; j < 6; j++ ) {
            // Only appending the internal cells
            if ( ( temp_rank = ( *global_local_index )[( *lcc )[i][j]][0] ) != my_rank ) {
                if ( !contains( ( *lcc )[i][j], ( *recv_list )[temp_rank],
                                ( *recv_count )[temp_rank] ) ) {
                    // Adding the unique global indices
                    ( *recv_list )[temp_rank][( ( *recv_count )[temp_rank] )++] = ( *lcc )[i][j];
                }

                if ( !contains( i, ( *send_list )[temp_rank], ( *send_count )[temp_rank] ) ) {
                    // Adding the global indices to be sent to each processor
                    ( *send_list )[temp_rank][( ( *send_count )[temp_rank] )++] = i;
                }
            }
        }
    }

    assert( ( *recv_count )[my_rank] == 0 );
    assert( ( *send_count )[my_rank] == 0 );
    // Now the ghost cells have to be allocated based on the recv counts

    return 0;
}

int contains( int index, int *search_array, int size ) {
    for ( int i = 0; i < size; i++ ) {
        if ( search_array[i] == index ) {
            return 1;
        }
    }
    return 0;
}

