/**
 * Finalization step - write results and other computational vectors to files
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdio.h>
#include "util_write_files.h"
#include <stdlib.h>
#include <mpi.h>

void finalization( char* file_in, char* out_prefix, int total_iters, double residual_ratio,
                   int nintci, int nintcf, int points_count, int** points, int* elems, double* var,
                   double* cgup, double* su, int *local_global_index,int local_int_cells,
                   int elemcount, int writing_proc ) {

    char file_out[100];
    sprintf(file_out, "%s_summary.out", out_prefix);

    MPI_Status stat;
    int my_rank, nproc;
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );    /// get current process id
    MPI_Comm_size( MPI_COMM_WORLD, &nproc );    /// get number of processes

    if ( my_rank == writing_proc ) {
        int proc_loc_int_size = 0;
        int *temp_elems = NULL;
        int *temp_loc_glo_index = NULL;
        double *temp_data_values_cgup;
        double *temp_data_values_var;

        // Gather all the elems into one processor
        int *global_elems = (int *) malloc( ( nintcf - nintci + 1 ) * 8 * sizeof(int) );

        // Gather all the data values into one processor
        double *data_values_cgup = (double *) malloc( ( nintcf - nintci + 1 ) * sizeof(double) );
        double *data_values_var = (double *) malloc( ( nintcf - nintci + 1 ) * sizeof(double) );

        for ( int i = 0; i < local_int_cells; i++ ) {
            for ( int j = 0; j < 8; j++ ) {
                global_elems[( local_global_index[i] - nintci ) * 8 + j] = elems[i * 8 + j];
            }
            data_values_cgup[local_global_index[i] - nintci] = cgup[i];
            data_values_var[local_global_index[i] - nintci] = var[i];
        }

        for ( int i = nproc - 1; i >= 0; i-- ) {
            if ( i != writing_proc ) {
                // MPI_Recv (&buf,count,datatype,source,tag,comm,&status)
                MPI_Recv( &proc_loc_int_size, 1, MPI_INT, i, i, MPI_COMM_WORLD, &stat );
                temp_elems = (int *) malloc( proc_loc_int_size * 8 * sizeof(int) );
                MPI_Recv( temp_elems, proc_loc_int_size * 8, MPI_INT, i, i, MPI_COMM_WORLD,
                          &stat );

                // Now where to place indices
                temp_loc_glo_index = (int *) malloc( proc_loc_int_size * sizeof(int) );
                MPI_Recv( temp_loc_glo_index, proc_loc_int_size, MPI_INT, i, i, MPI_COMM_WORLD,
                          &stat );

                temp_data_values_cgup = (double *) malloc( proc_loc_int_size * sizeof(double) );
                MPI_Recv( temp_data_values_cgup, proc_loc_int_size, MPI_DOUBLE, i, i, MPI_COMM_WORLD,
                          &stat );

                temp_data_values_var = (double *) malloc( proc_loc_int_size * sizeof(double) );
                MPI_Recv( temp_data_values_var, proc_loc_int_size, MPI_DOUBLE, i, i, MPI_COMM_WORLD,
                          &stat );

                // Need calculation
                for ( int i = 0; i < proc_loc_int_size; i++ ) {
                    for ( int j = 0; j < 8; j++ ) {
                        global_elems[temp_loc_glo_index[i] * 8 + j] = temp_elems[i * 8 + j];
                    }
                    data_values_cgup[temp_loc_glo_index[i]] = temp_data_values_cgup[i];
                    data_values_var[temp_loc_glo_index[i]] = temp_data_values_var[i];
                }
            }
        }

        int status = store_simulation_stats(file_in, file_out, nintci, nintcf, data_values_var,
                                            total_iters, residual_ratio );
        sprintf(file_out, "%s_data.vtk", out_prefix);
        vtk_write_unstr_grid_header( file_in, file_out, nintci, nintcf, points_count, points,
                                     global_elems );
        vtk_append_double( file_out, "CGUP", nintci, nintcf, data_values_cgup );
        if ( status != 0 ) fprintf(stderr, "Error when trying to write to file %s\n", file_out);

        free( global_elems );
        free( data_values_cgup );
        free( data_values_var );
        free( temp_elems );
        free( temp_loc_glo_index );
        free( temp_data_values_cgup );
        free( temp_data_values_var );
    } else {
        // MPI_Send (&buf,count,datatype,dest,tag,comm)
        MPI_Send( &local_int_cells, 1, MPI_INT, writing_proc, my_rank, MPI_COMM_WORLD );
        MPI_Send( elems, local_int_cells * 8, MPI_INT, writing_proc, my_rank, MPI_COMM_WORLD );
        MPI_Send( local_global_index, local_int_cells, MPI_INT, writing_proc, my_rank,
                  MPI_COMM_WORLD );
        MPI_Send( cgup, local_int_cells, MPI_DOUBLE, writing_proc, my_rank, MPI_COMM_WORLD );
        MPI_Send( var, local_int_cells, MPI_DOUBLE, writing_proc, my_rank, MPI_COMM_WORLD );
    }

}
