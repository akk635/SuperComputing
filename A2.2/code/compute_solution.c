/**
 * Computational loop
 *
 * @file compute_solution.c
 * @date 22-Oct-2012
 * @author V. Petkov
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <assert.h>
#include <papi.h>

double getRecv( int global_index, double **recv_buffer, int **recv_list, int temp_rank, int *count );
int compar( const void *, const void * );

int compute_solution( const int max_iters, int nintci, int nintcf, int nextcf, int** lcc,
                      double* bp, double* bs, double* bw, double* bl, double* bn, double* be,
                      double* bh, double* cnorm, double* var, double *su, double* cgup,
                      double* residual_ratio, int* local_global_index, int** global_local_index,
                      int neighbors_count, int* send_count, int** send_list, int* recv_count,
                      int** recv_list, int elemcount, int local_int_cells ) {
    /** parameters used in gccg */
    int iter = 1;
    int if1 = 0;
    int if2 = 0;
    int nor = 1;
    int nor1 = nor - 1;
    int nc = 0;
    int nomax = 3;

    int my_rank, nproc;
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &nproc );
    MPI_Request request[2 * nproc];
    MPI_Status status[nproc];

    int **blocklengths;
    blocklengths = (int **) malloc( nproc * sizeof(int *) );
    for ( int i = 0; i < nproc; i++ ) {
        ( blocklengths )[i] = (int *) calloc( sizeof(int), send_count[i] );
    }
    // Initialising with all 1's
    for ( int i = 0; i < nproc; i++ ) {
        for ( int j = 0; j < send_count[i]; j++ ) {
            blocklengths[i][j] = 1;
        }
    }

    MPI_Datatype *indextype;
    indextype = (MPI_Datatype *) malloc( nproc * sizeof(MPI_Datatype) );
    for ( int i = 0; i < nproc; i++ ) {
        MPI_Type_indexed( send_count[i], blocklengths[i], send_list[i], MPI_DOUBLE, indextype + i );
        MPI_Type_commit( indextype + i );
    }

    // For receiving the communication
    double **recv_buffer;
    recv_buffer = (double **) malloc( nproc * sizeof(double *) );
    for ( int i = 0; i < nproc; i++ ) {
        ( recv_buffer )[i] = (double*) malloc( recv_count[i] * sizeof(double) );
    }

#define solve(global_index) \
        (global_local_index[global_index][0]) == (my_rank) ? (direc1[global_local_index[global_index][1]]) \
                : getRecv( global_index, recv_buffer, recv_list, global_local_index[global_index][0], recv_count )

    /** the reference residual */
    double resref = 0.0;

    /** array storing residuals */
    double *resvec = (double *) calloc( local_int_cells, sizeof(double) );

    // initialize the reference residual
    for ( nc = 0; nc < local_int_cells; nc++ ) {
        resvec[nc] = su[nc];
        resref = resref + resvec[nc] * resvec[nc];
    }

    // MPI_Allreduce (&sendbuf,&recvbuf,count,datatype,op,comm)
    MPI_Allreduce( MPI_IN_PLACE, &resref, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    resref = sqrt( resref );
    if ( resref < 1.0e-15 ) {
        fprintf( stderr, "Residue sum less than 1.e-15 - %lf\n", resref );
        return 0;
    }

    /** the computation vectors */
    double *direc1 = (double *) calloc( elemcount, sizeof(double) );
    double *direc2 = (double *) calloc( elemcount, sizeof(double) );
    double *adxor1 = (double *) calloc( local_int_cells, sizeof(double) );
    double *adxor2 = (double *) calloc( local_int_cells, sizeof(double) );
    double *dxor1 = (double *) calloc( local_int_cells, sizeof(double) );
    double *dxor2 = (double *) calloc( local_int_cells, sizeof(double) );

    while ( iter < max_iters ) {
        /**********  START COMP PHASE 1 **********/
        // update the old values of direc
        for ( nc = 0; nc < local_int_cells; nc++ ) {
            direc1[nc] = direc1[nc] + resvec[nc] * cgup[nc];
        }

        // Communication of the neighbors values
        for ( int i = 0; i < nproc; i++ ) {
            if ( send_count[i] > 0 ) {
                // MPI_Isend (&buf,count,datatype,dest,tag,comm,&request)
                MPI_Isend( direc1, 1, indextype[i], i, i, MPI_COMM_WORLD, request + i );
            }
        }

        for ( int i = 0; i < nproc; i++ ) {
            if ( recv_count[i] > 0 ) {
                // MPI_Recv (&buf,count,datatype,source,tag,comm,&status)
                // Blocked until exact values are received at the application buffer
                MPI_Recv( recv_buffer[i], recv_count[i], MPI_DOUBLE, i, my_rank, MPI_COMM_WORLD,
                          status + i );
                // MPI_Irecv(buffer,count,type,source,tag,comm,request)
/*                MPI_Irecv( recv_buffer[i], recv_count[i], MPI_DOUBLE, i, my_rank, MPI_COMM_WORLD,
                          nproc + request + i );*/
            }
        }

/*
        for ( int i = 0; i < nproc; i++ ) {
            if ( i != my_rank ) {
                MPI_Wait( request + i, status + i );
                MPI_Wait( nproc + request + i, status + i );
            }
        }
*/

        // compute new guess (approximation) for direc
        for ( nc = 0; nc < local_int_cells; nc++ ) {
            direc2[nc] = bp[nc] * direc1[nc] - bs[nc] * ( solve( lcc[nc][0] ) )
                    - be[nc] * ( solve( lcc[nc][1] ) ) - bn[nc] * ( solve( lcc[nc][2] ) )
                    - bw[nc] * ( solve( lcc[nc][3] ) ) - bl[nc] * ( solve( lcc[nc][4] ) )
                    - bh[nc] * ( solve( lcc[nc][5] ) );
        }
        /********** END COMP PHASE 1 **********/

        /********** START COMP PHASE 2 **********/
        // execute normalization steps
        double oc1, oc2, occ;
        if ( nor1 == 1 ) {
            oc1 = 0;
            occ = 0;

            for ( nc = 0; nc < local_int_cells; nc++ ) {
                occ = occ + direc2[nc] * adxor1[nc];
            }
            // MPI_Allreduce (&sendbuf,&recvbuf,count,datatype,op,comm)
            MPI_Allreduce( MPI_IN_PLACE, &occ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

            oc1 = occ / cnorm[1];
            for ( nc = 0; nc < local_int_cells; nc++ ) {
                direc2[nc] = direc2[nc] - oc1 * adxor1[nc];
                direc1[nc] = direc1[nc] - oc1 * dxor1[nc];
            }

            if1++;
        } else {
            if ( nor1 == 2 ) {
                oc1 = 0;
                occ = 0;

                for ( nc = 0; nc < local_int_cells; nc++ ) {
                    occ = occ + direc2[nc] * adxor1[nc];
                }
                // MPI_Allreduce (&sendbuf,&recvbuf,count,datatype,op,comm)
                MPI_Allreduce( MPI_IN_PLACE, &occ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

                oc1 = occ / cnorm[1];
                oc2 = 0;
                occ = 0;
                for ( nc = 0; nc < local_int_cells; nc++ ) {
                    occ = occ + direc2[nc] * adxor2[nc];
                }
                // MPI_Allreduce (&sendbuf,&recvbuf,count,datatype,op,comm)
                MPI_Allreduce( MPI_IN_PLACE, &occ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

                oc2 = occ / cnorm[2];
                for ( nc = 0; nc < local_int_cells; nc++ ) {
                    direc1[nc] = direc1[nc] - oc1 * dxor1[nc] - oc2 * dxor2[nc];
                    direc2[nc] = direc2[nc] - oc1 * adxor1[nc] - oc2 * adxor2[nc];
                }

                if2++;
            }
        }

        // compute the new residual
        cnorm[nor] = 0;
        double omega = 0;
        for ( nc = 0; nc < local_int_cells; nc++ ) {
            cnorm[nor] = cnorm[nor] + direc2[nc] * direc2[nc];
            omega = omega + resvec[nc] * direc2[nc];
        }
        // MPI_Allreduce (&sendbuf,&recvbuf,count,datatype,op,comm)
        MPI_Allreduce( MPI_IN_PLACE, &( cnorm[nor] ), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        // MPI_Allreduce (&sendbuf,&recvbuf,count,datatype,op,comm)
        MPI_Allreduce( MPI_IN_PLACE, &omega, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

        omega = omega / cnorm[nor];
        double res_updated = 0.0;
        for ( nc = 0; nc < local_int_cells; nc++ ) {
            resvec[nc] = resvec[nc] - omega * direc2[nc];
            res_updated = res_updated + resvec[nc] * resvec[nc];
            var[nc] = var[nc] + omega * direc1[nc];
        }
        // MPI_Allreduce (&sendbuf,&recvbuf,count,datatype,op,comm)
        MPI_Allreduce( MPI_IN_PLACE, &res_updated, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

        res_updated = sqrt( res_updated );
        *residual_ratio = res_updated / resref;

        // exit on no improvements of residual
        if ( *residual_ratio <= 1.0e-10 )
            break;

        iter++;

        // prepare additional arrays for the next iteration step
        if ( nor == nomax ) {
            nor = 1;
        } else {
            if ( nor == 1 ) {
                for ( nc = 0; nc < local_int_cells; nc++ ) {
                    dxor1[nc] = direc1[nc];
                    adxor1[nc] = direc2[nc];
                }
            } else {
                if ( nor == 2 ) {
                    for ( nc = 0; nc < local_int_cells; nc++ ) {
                        dxor2[nc] = direc1[nc];
                        adxor2[nc] = direc2[nc];
                    }
                }
            }

            nor++;
        }
        nor1 = nor - 1;
        /********** END COMP PHASE 2 **********/
    }

    free( direc1 );
    free( direc2 );
    free( adxor1 );
    free( adxor2 );
    free( dxor1 );
    free( dxor2 );
    free( resvec );
    for ( int i = 0; i < nproc; i++ ) {
        free( blocklengths[i] );
        free( recv_buffer[i] );
        MPI_Type_free( indextype + i );
    }

    return iter;
}

double getRecv( int global_index, double **recv_buffer, int **recv_list, int temp_rank, int *count ) {
    const int *temp_pointer = &global_index;
    const int* base = recv_list[temp_rank];
    int arraysize = count[temp_rank];
    int *index = (int *) bsearch( (void *) temp_pointer, (void *) base, (size_t) arraysize,
                                  sizeof(int), compar );
    assert( index != NULL );
    assert( ( index - base ) == ( (char * ) index - (char* ) base ) / sizeof(int) );
    return recv_buffer[temp_rank][index - base];
}

int compar( const void * a, const void *b ) {
    int *A = (int *) a;
    int *B = (int *) b;
    return ( *A ) == ( *B ) ? 0 : ( ( *A ) - ( *B ) );
}

