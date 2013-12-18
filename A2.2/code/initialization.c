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

/*int contains( int index, int *search_array, int size );*/

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
    MPI_Status status[nproc];
    MPI_Request request[2 * nproc];

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
    // Counter for the internal domain cells for the recv_count calc
    int *counter_int_cells = (int *) calloc( ( ( *nintcf ) - ( *nintci ) + 1 ), sizeof(int) );
    assert( counter_int_cells != NULL );
    // Making it global access
    counter_int_cells = counter_int_cells - ( *nintci );
    // For the last send cell
    int *last_send = (int *) malloc( nproc * sizeof(int) );

    // Initialising with -1's
    for ( int i = 0; i < nproc; i++ ) {
        last_send[i] = -1;
    }

    printf( "working \n" );
    // Temporary counts just for the allocation
    /*    int *neighbors = calloc( nproc, sizeof(int) );
     // Provides with the max size of the array
     for ( int i = 0; i < ( *local_int_cells ); i++ ) {
     for ( int j = 0; j < 6; j++ ) {
     // Only appending the internal cells
     if ( ( temp_rank = ( *global_local_index )[( *lcc )[i][j]][0] ) != my_rank ) {
     neighbors[temp_rank] += 1;
     }
     }
     }*/

    // macro for counting the send_count
#define match( rank, index, list ) \
        (list[rank]) == (index) ? (int)1 : (int)0

    printf( "true : %d \n", match( my_rank, 0, last_send ) );

    // For the actual count
    *recv_count = (int *) calloc( nproc, sizeof(int) );
    *send_count = (int *) calloc( nproc, sizeof(int) );

    for ( int i = 0; i < ( *local_int_cells ); i++ ) {
        for ( int j = 0; j < 6; j++ ) {
            if ( ( temp_rank = ( *global_local_index )[( *lcc )[i][j]][0] ) != my_rank ) {
                if ( counter_int_cells[( *lcc )[i][j]] == 0 ) {
                    counter_int_cells[( *lcc )[i][j]] = 1;
                    ( ( *recv_count )[temp_rank] )++;
                }
                if ( !( match( temp_rank, i, last_send ) ) ) {
                    last_send[temp_rank] = i;
                    ( ( *send_count )[temp_rank] )++;
                }
            }
        }
    }
    if ( my_rank == 1 ) {
        printf( " send_count : %d \n", ( *send_count )[0] );
    }

    ( *send_list ) = (int **) malloc( nproc * sizeof(int *) );
    // Allocating with the max array size
    for ( int i = 0; i < nproc; i++ ) {
        ( *send_list )[i] = (int *) calloc( ( *send_count )[i], sizeof(int) );
    }

    ( *recv_list ) = (int **) malloc( nproc * sizeof(int *) );

    // Allocating the max array size
    for ( int i = 0; i < nproc; i++ ) {
        ( *recv_list )[i] = (int *) calloc( ( *recv_count )[i], sizeof(int) );
    }

    int *counter = (int *) calloc( nproc, sizeof(int) );

    for ( int i = 0; i < ( *local_int_cells ); i++ ) {
        for ( int j = 0; j < 6; j++ ) {
            // Only appending the internal cells
            if ( ( temp_rank = ( *global_local_index )[( *lcc )[i][j]][0] ) != my_rank ) {
                if ( !( match( temp_rank, i, last_send ) ) ) {
                    // Adding the global indices to be sent to each processor
                    last_send[temp_rank] = i;
                    ( *send_list )[temp_rank][( counter[temp_rank] )++] = i;
                }
            }
        }
    }

    for ( int i = 0; i < nproc; i++ ) {
        assert( counter[i] == ( *send_count )[i] );
    }

    int **index_send_list = (int **) malloc( nproc * sizeof(int *) );
    // Allocating an initialising
    for ( int i = 0; i < nproc; i++ ) {
        index_send_list[i] = (int *) malloc( ( *send_count )[i] * sizeof(int) );
        assert( index_send_list[i] != NULL );
        for ( int j = 0; j < ( *send_count )[i]; j++ ) {
            index_send_list[i][j] = ( *local_global_index )[( *send_list )[i][j]];
        }
    }
    // Correcting the list using communication
    // MPI_Sendrecv( &sendbuf, sendcount, sendtype, dest, sendtag, &recvbuf, recvcount, recvtype,
    // source, recvtag, comm, &status );
    for ( int i = 0; i < nproc; i++ ) {
        if ( ( *send_count )[i] > 0 ) {
            // MPI_Isend (&buf,count,datatype,dest,tag,comm,&request)
            MPI_Isend( index_send_list[i], ( *send_count )[i], MPI_INT, i, i, MPI_COMM_WORLD,
                       request + i );

            /*            MPI_Sendrecv( index_send_list[i], ( *send_count )[i], MPI_INT, i, i, recv_list[i],
             ( *recv_count )[i], MPI_INT, i, my_rank, MPI_COMM_WORLD, &status );*/

            // MPI_Recv (&buf,count,datatype,source,tag,comm,&status)
            /*            MPI_Recv( ( *recv_list )[i], ( *recv_count )[i], MPI_INT, i, my_rank, MPI_COMM_WORLD,
             status + i );*/
            // MPI_Irecv(buffer,count,type,source,tag,comm,request)
            MPI_Irecv( ( *recv_list )[i], ( *recv_count )[i], MPI_INT, i, my_rank, MPI_COMM_WORLD,
                       nproc + request + i );
        }
    }

    for ( int i = 0; i < nproc; i++ ) {
        if ( i != my_rank ) {
            MPI_Wait( request + i, status + i );
            MPI_Wait( nproc + request + i, status + i );
        }
    }

    counter_int_cells = counter_int_cells + ( *nintci );
    printf( "no dead lock \n" );
    // Freeing the buffers
    for ( int i = 0; i < nproc; i++ ) {
        free( index_send_list[i] );
    }

    free( counter );
    free( counter_int_cells );
    free( last_send );

    return 0;
}

/*
 int contains( int index, int *search_array, int size ) {
 for ( int i = size - 1; i >= 0; i-- ) {
 if ( search_array[i] == index ) {
 return 1;
 }
 }
 return 0;
 }
 */

