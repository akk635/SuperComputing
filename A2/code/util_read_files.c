/**
 * Helper functions for reading from input data file
 *
 * @author E. Xue, V. Petkov
 * @date 22-May-2009, 22-Oct-2012
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util_read_files.h"
#include <assert.h>
#include "metis.h"

/**
 * Parse an binary input data set and initialize simulation variables
 *
 * @param file_name
 * @param NINTCI
 * @param NINTCF
 * @param NEXTCI
 * @param NEXTCF
 * @param LCC
 * @param BS
 * @param BE
 * @param BN
 * @param BW
 * @param BL
 * @param BH
 * @param BP
 * @param SU
 * @param points_count
 * @param points
 * @param elems
 * @return
 */
int read_binary_geo(char *file_name, char* part_type, int *NINTCI, int *NINTCF, int *NEXTCI, int *NEXTCF, int ***LCC,
		double **BS, double **BE, double **BN, double **BW, double **BL, double **BH, double **BP,
		double **SU, int* points_count, int*** points, int** elems, int **local_global_index,
		int *elemcount, int *local_int_cells, int ***global_local_index, int **epart, int **npart, int *objval ) {

    int i = 0;
    int my_rank, nproc;
    MPI_Status status;
    FILE *fp = fopen( file_name, "rb" );

    if ( fp == NULL ) {
        fprintf( stderr, "Error opening file %s\n", file_name );
        return -1;
    }

    // 4 variables in total!!!
    fread( NINTCI, sizeof(int), 1, fp );
    fread( NINTCF, sizeof(int), 1, fp );
    fread( NEXTCI, sizeof(int), 1, fp );
    fread( NEXTCF, sizeof(int), 1, fp );

    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &nproc );

    int tot_domain_cells = *NEXTCF - *NINTCI + 1;

    // For storing ranks according to the indices
    ( *epart ) = (int *) malloc( sizeof(int) * tot_domain_cells );
    int *distr_buffer = ( *epart );

    if ( strcmp( part_type, "classical" ) == 0 ) {

        if ( my_rank == 0 ) {
            int fpcount;
            fpcount = *NINTCI;

            // So that we can call on explicit coordinates
            distr_buffer = distr_buffer - *NINTCI;

            int local_cells_size;

            int normal_local_size = tot_domain_cells / nproc;
            int res_cells = tot_domain_cells % nproc;

            // Initializing the distribution array with -1
            for ( int i = *NINTCI; i < *NEXTCF + 1; i++ ) {
                distr_buffer[i] = -1;
            }

            // All the neighbors do not cover all the cells.
            // there are some external cells which are not neighbors of internal
            int temp_cells_size = 0;

            // To get the indices to write in the distr_buffer
            int temp_buffer;

            // equally distributing the local size for each processor
            for ( int i = nproc - 1; i >= 0; i-- ) {
                *elemcount = 0;
                local_cells_size = normal_local_size;
                // Reading the topological info and then distributing according to the locality
                if ( res_cells > 0 ) {
                    res_cells--;
                    local_cells_size++;
                }

                while ( ( *elemcount ) < local_cells_size ) {

                    if ( fpcount == *NINTCF + 1 ) {
                        break;
                    }

                    if ( distr_buffer[fpcount] == -1 ) {
                        distr_buffer[fpcount] = i;
                        ( *elemcount )++;
                        temp_cells_size++;
                    }

                    for ( int j = 0; j < 6; j++ ) {
                        fread( &temp_buffer, sizeof(int), 1, fp );
                        if ( distr_buffer[temp_buffer] == -1 ) {
                            distr_buffer[temp_buffer] = i;
                            ( *elemcount )++;
                            temp_cells_size++;
                        }
                    }
                    fpcount++;
                }

                // Sending the internal cells indices to the respective processes
                // MPI_Send (&buf,count,datatype,dest,tag,comm)
                if ( i != 0 ) {
                    MPI_Send( elemcount, 1, MPI_INT, i, i, MPI_COMM_WORLD );
                }
            }

            // Now distributing all the remaining external cells to process 0
            for ( int i = *NEXTCI; i <= *NEXTCF; i++ ) {
                if ( distr_buffer[i] == -1 ) {
                    distr_buffer[i] = 0;
                    ( *elemcount )++;
                    temp_cells_size++;
                }
            }
            assert( temp_cells_size == tot_domain_cells );
            // Now distributing the buffer to all the processors
            for ( int i = 1; i < nproc; i++ ) {
                // MPI_Send (&buf,count,datatype,dest,tag,comm)
                MPI_Send( distr_buffer + *NINTCI, tot_domain_cells, MPI_INT, i, i, MPI_COMM_WORLD );
            }
        } else {
            MPI_Recv( elemcount, 1, MPI_INT, 0, my_rank, MPI_COMM_WORLD, &status );
            // MPI_Recv (&buf,count,datatype,source,tag,comm,&status)
            MPI_Recv( distr_buffer, tot_domain_cells, MPI_INT, 0, my_rank, MPI_COMM_WORLD,
                      &status );
        }

    } else {
        if ( my_rank == 0 ) {
            int lcc_read_end = ( ( *NINTCF - *NINTCI + 1 ) * 6 + 4 ) * sizeof(int);
            int coe_read_end = lcc_read_end + 8 * ( *NINTCF - *NINTCI + 1 ) * sizeof(double);

            fseek( fp, coe_read_end, SEEK_SET );

            // read geometry
            // allocate elems
            if ( ( *elems = (int*) malloc( ( *NINTCF - *NINTCI + 1 ) * 8 * sizeof(int) ) ) == NULL ) {
                fprintf( stderr, "malloc(elems) failed" );
                return -1;
            }

            // read elems
            for ( i = ( *NINTCI ); i < ( *NINTCF + 1 ) * 8; i++ ) {
                fread( &( ( *elems )[i] ), sizeof(int), 1, fp );
            }

            fread( points_count, sizeof(int), 1, fp );

            idx_t ne = *NINTCF - *NINTCI + 1;
            idx_t nn = *points_count;
            idx_t *eptr;
            idx_t *eind;
            idx_t *vwgt = NULL;
            idx_t *vsize = NULL;
            idx_t ncommon = 4;
            idx_t nparts = nproc;
            real_t *tpwgts = NULL;
            idx_t options[METIS_NOPTIONS];
            idx_t *temp_epart;
            idx_t *temp_npart;

            METIS_SetDefaultOptions( options );

            if ( ( eptr = (idx_t *) malloc( ( ne + 1 ) * sizeof(idx_t) ) ) == NULL ) {
                fprintf( stderr, "malloc(eptr) failed\n" );
                return -1;
            }

            if ( ( eind = (idx_t *) malloc( ( ne * 8 ) * sizeof(idx_t) ) ) == NULL ) {
                fprintf( stderr, "malloc(eind) failed\n" );
                return -1;
            }

            for ( int i = 0; i < ne + 1; i++ ) {
                eptr[i] = i * 8;
            }

            for ( int i = 0; i < ne * 8; i++ ) {
                eind[i] = ( *elems )[i];
            }

            if ( ( temp_epart = (idx_t *) malloc( ( ne ) * sizeof(idx_t) ) ) == NULL ) {
                fprintf( stderr, "malloc(epart) failed\n" );
                return -1;
            }

            if ( ( temp_npart = (idx_t *) malloc( ( nn ) * sizeof(idx_t) ) ) == NULL ) {
                fprintf( stderr, "malloc(npart) failed\n" );
                return -1;
            }

            if ( ( ( *npart ) = (int *) malloc( ( nn ) * sizeof(int) ) ) == NULL ) {
                fprintf( stderr, "malloc(npart) failed\n" );
                return -1;
            }

            if ( strcmp( part_type, "dual" ) == 0 ) {
                if ( METIS_PartMeshDual( &ne, &nn, eptr, eind, vwgt, vsize, &ncommon, &nparts,
                                         tpwgts, options, (idx_t *) objval, temp_epart, temp_npart )
                        != METIS_OK ) {
                    fprintf( stderr, "Partitioning of METIS DUAL FAILED!\n" );
                    return -1;
                }

            } else if ( strcmp( part_type, "nodal" ) == 0 ) {
                if ( METIS_PartMeshNodal( &ne, &nn, eptr, eind, vwgt, vsize, &nparts, tpwgts,
                                          options, (idx_t *) objval, temp_epart, temp_npart )
                        != METIS_OK ) {
                    fprintf( stderr, "Partitioning of METIS NODAL FAILED!\n" );
                    return -1;
                }
            }

            // Initializing the distribution array with -1
            for ( int i = 0; i < *NEXTCF - *NINTCI + 1; i++ ) {
                ( *epart )[i] = -1;
            }

            for ( int i = 0; i < *NINTCF - *NINTCI + 1; i++ ) {
                ( *epart )[i] = (int) temp_epart[i];
            }
            for ( int i = 0; i < *points_count; i++ ) {
                ( *npart )[i] = (int) temp_npart[i];
            }

            // Adjusting the position of the dsitr_array to the global position
            distr_buffer = distr_buffer - *NINTCI;
            ( *elemcount ) = 0;

            free( temp_epart );
            free( temp_npart );
            free( eptr );
            free( eind );

            // Distributing the things to other processors
            for ( int i = nproc - 1; i > 0; i-- ) {
                // MPI_Send (&buf,count,datatype,dest,tag,comm)
                MPI_Send( ( *epart ), tot_domain_cells, MPI_INT, i, i, MPI_COMM_WORLD );
            }
        } else {
            // MPI_Recv (&buf,count,datatype,source,tag,comm,&status)
            MPI_Recv( ( *epart ), tot_domain_cells, MPI_INT, 0, my_rank, MPI_COMM_WORLD, &status );
            ( *elemcount ) = 0;
        }
    }

    if ( my_rank != 0 ) {
        distr_buffer = distr_buffer - *NINTCI;
    }

    ( *global_local_index ) = (int **) malloc( ( *NEXTCF - *NINTCI + 1 ) * sizeof(int *) );
    ( *global_local_index ) = ( *global_local_index ) - *NINTCI;

    for ( int i = *NINTCI; i <= *NEXTCF; i++ ) {
        ( *global_local_index )[i] = (int *) calloc( 2, sizeof(int) );
    }

    *local_int_cells = 0;

    for ( int i = *NINTCI; i <= *NINTCF; i++ ) {
        if ( distr_buffer[i] == my_rank ) {
            ( *local_int_cells )++;
        }
    }

    *local_global_index = malloc( ( *local_int_cells ) * sizeof(int) );

    // Only filling the internal cells in the global_local_index
    int j = 0;
    for ( int i = *NINTCI; i <= *NINTCF; i++ ) {
        ( *global_local_index )[i][0] = distr_buffer[i];
        if ( distr_buffer[i] == my_rank ) {
            ( *local_global_index )[j] = i;
            ( *global_local_index )[i][1] = j;
            j++;
        }
    }

    int check;
    check = distr_buffer[*NINTCF];

    assert( j == ( *local_int_cells ) );
    j = 0;

    // Allocating the LCC in individual processes
    if ( ( *LCC = (int**) malloc( ( *local_int_cells ) * sizeof(int*) ) ) == NULL ) {
        fprintf( stderr, "malloc failed to allocate first dimension of LCC" );
        return -1;
    }

    for ( i = 0; i < ( *local_int_cells ); i++ ) {
        if ( ( ( *LCC )[i] = (int *) malloc( 6 * sizeof(int) ) ) == NULL ) {
            fprintf( stderr, "malloc failed to allocate second dimension of LCC\n" );
            return -1;
        }
    }

    int index_read = 4 * sizeof(int);

    if ( strcmp( part_type, "classical" ) == 0 ) {
        // Start reading LCC
        for ( int i = 0; i < ( *local_int_cells ); i++ ) {
            fseek( fp, index_read + ( ( *local_global_index )[i] - *NINTCI ) * 6 * sizeof(int),
                   SEEK_SET );
            for ( int j = 0; j < 6; j++ ) {
                fread( &( *LCC )[i][j], sizeof(int), 1, fp );
            }
        }
    } else {
        ( *elemcount ) = ( *local_int_cells );
        // Start reading LCC
        for ( int i = 0; i < ( *local_int_cells ); i++ ) {
            fseek( fp, index_read + ( ( *local_global_index )[i] - *NINTCI ) * 6 * sizeof(int),
                   SEEK_SET );
            for ( int j = 0; j < 6; j++ ) {
                fread( &( *LCC )[i][j], sizeof(int), 1, fp );
                if ( ( ( *LCC )[i][j] > *NINTCF ) & ( distr_buffer[( ( *LCC )[i][j] )] == -1 ) ) {
                    distr_buffer[( ( *LCC )[i][j] )] = my_rank;
                    ( *elemcount )++;
                }
            }
        }
    }

    j = ( *local_int_cells );
    // Now both metis and classical are at the same position
    for ( int i = *NEXTCI; i <= *NEXTCF; i++ ) {
        ( *global_local_index )[i][0] = distr_buffer[i];
        if ( distr_buffer[i] == my_rank ) {
            ( *global_local_index )[i][1] = j;
            j++;
        }
    }

    assert( j == ( *elemcount ) );
    j = 0;

    // allocate other arrays
    if ( ( *BS = (double *) malloc( ( *elemcount ) * sizeof(double) ) ) == NULL ) {
        fprintf( stderr, "malloc(BS) failed\n" );
        return -1;
    }

    if ( ( *BE = (double *) malloc( ( *elemcount ) * sizeof(double) ) ) == NULL ) {
        fprintf( stderr, "malloc(BE) failed\n" );
        return -1;
    }

    if ( ( *BN = (double *) malloc( ( *elemcount ) * sizeof(double) ) ) == NULL ) {
        fprintf( stderr, "malloc(BN) failed\n" );
        return -1;
    }

    if ( ( *BW = (double *) malloc( ( *elemcount ) * sizeof(double) ) ) == NULL ) {
        fprintf( stderr, "malloc(BW) failed\n" );
        return -1;
    }

    if ( ( *BL = (double *) malloc( ( *elemcount ) * sizeof(double) ) ) == NULL ) {
        fprintf( stderr, "malloc(BL) failed\n" );
        return -1;
    }

    if ( ( *BH = (double *) malloc( ( *elemcount ) * sizeof(double) ) ) == NULL ) {
        fprintf( stderr, "malloc(BH) failed\n" );
        return -1;
    }

    if ( ( *BP = (double *) malloc( ( *elemcount ) * sizeof(double) ) ) == NULL ) {
        fprintf( stderr, "malloc(BP) failed\n" );
        return -1;
    }

    if ( ( *SU = (double *) malloc( ( *elemcount ) * sizeof(double) ) ) == NULL ) {
        fprintf( stderr, "malloc(SU) failed\n" );
        return -1;
    }

    // read the other arrays
    int lcc_read_end = ( ( *NINTCF - *NINTCI + 1 ) * 6 + 4 ) * sizeof(int);

    if ( my_rank == check ) {
        long int fpos = ftell( fp );
        assert( fpos == lcc_read_end );
    }

    for ( i = 0; i < ( *local_int_cells ); i++ ) {
        fseek( fp, lcc_read_end + ( ( *local_global_index )[i] - *NINTCI ) * 8 * sizeof(double),
               SEEK_SET );
        fread( &( ( *BS )[i] ), sizeof(double), 1, fp );
        fread( &( ( *BE )[i] ), sizeof(double), 1, fp );
        fread( &( ( *BN )[i] ), sizeof(double), 1, fp );
        fread( &( ( *BW )[i] ), sizeof(double), 1, fp );
        fread( &( ( *BL )[i] ), sizeof(double), 1, fp );
        fread( &( ( *BH )[i] ), sizeof(double), 1, fp );
        fread( &( ( *BP )[i] ), sizeof(double), 1, fp );
        fread( &( ( *SU )[i] ), sizeof(double), 1, fp );
    }

    int coe_read_end = lcc_read_end + 8 * ( *NINTCF - *NINTCI + 1 ) * sizeof(double);

    if ( my_rank == check ) {
        long int fpos = ftell( fp );
        assert( fpos == coe_read_end );
    }

    // read geometry i.e nodes
    // allocate elems
    if ( ( *elems = (int*) malloc( ( *local_int_cells ) * 8 * sizeof(int) ) ) == NULL ) {
        fprintf( stderr, "malloc failed to allocate elems" );
        return -1;
    }

    // read elems
    for ( i = 0; i < ( *local_int_cells ); i++ ) {
        fseek( fp, coe_read_end + ( ( *local_global_index )[i] - *NINTCI ) * 8 * sizeof(int),
               SEEK_SET );
        for ( int j = 0; j < 8; j++ ) {
            fread( &( ( *elems )[( i * 8 ) + j] ), sizeof(int), 1, fp );
        }
    }

    int nodes_read_end = coe_read_end + 8 * ( *NINTCF - *NINTCI + 1 ) * sizeof(int);

    if ( my_rank == check ) {
        long int fpos = ftell( fp );
        assert( fpos == nodes_read_end );
    }

    fseek( fp, nodes_read_end, SEEK_SET );

    fread( points_count, sizeof(int), 1, fp );

    // allocate points vec
    if ( ( *points = (int **) calloc( *points_count, sizeof(int*) ) ) == NULL ) {
        fprintf( stderr, "malloc() POINTS 1st dim. failed\n" );
        return -1;
    }

    for ( i = 0; i < *points_count; i++ ) {
        if ( ( ( *points )[i] = (int *) calloc( 3, sizeof(int) ) ) == NULL ) {
            fprintf( stderr, "malloc() POINTS 2nd dim. failed\n" );
            return -1;
        }
    }

    int coordIdx;
    int pointIdx;
    for ( pointIdx = 0; pointIdx < *points_count; pointIdx++ ) {
        for ( coordIdx = 0; coordIdx < 3; coordIdx++ ) {
            fread( &( ( *points )[pointIdx][coordIdx] ), sizeof(int), 1, fp );
        }
    }

	fclose(fp);

	return 0;
}


