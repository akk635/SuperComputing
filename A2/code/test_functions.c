/**
 * Functions to test the data distribution and communication lists creation algorithms
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */
#include "test_functions.h"
#include <assert.h>
#include <stdlib.h>

int test_distribution(char *file_in, char *file_vtk_out, int *local_global_index,
                      int **global_local_index, int nintci, int nintcf, int points_count,
                      int **points, int *elems, int local_int_cells, double *cgup, int elemcount) {
    // Return an error if not implemented
    //Checking for the ranks of the elements
    
    MPI_Status status;
    int my_rank, nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);    /// get number of processes
    
    if (my_rank == 0){
        int proc_loc_int_size = 0;
        int *temp_elems = NULL;
        int *temp_loc_glo_index = NULL;
        double *temp_data_values;
        
        // Gather all the elems into one processor
        int *global_elems = (int *) malloc( ( nintcf - nintci + 1 ) * 8 * sizeof(int) );
        
        // Gather all the data values into one processor
        double *data_values = (double *) malloc((nintcf - nintci + 1) * sizeof(double));
        
        for ( int i = 0; i < local_int_cells; i++ ){
            for ( int j = 0; j < 8; j++ ){
                global_elems[ ( local_global_index[ i ] - nintci ) * 8 + j ] = elems[ i * 8 + j];
            }
            data_values[ local_global_index[ i ] - nintci ] = cgup[i];
        }
        
        for (int i = nproc - 1; i > 0; i--){
            // MPI_Recv (&buf,count,datatype,source,tag,comm,&status)
            MPI_Recv( &proc_loc_int_size, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status );
            temp_elems = (int *) malloc( proc_loc_int_size * 8 * sizeof(int));
            MPI_Recv( temp_elems, proc_loc_int_size * 8, MPI_INT, i, i, MPI_COMM_WORLD, &status );
            
            // Now where to place indices
            temp_loc_glo_index = (int *) malloc( proc_loc_int_size * sizeof(int));
            MPI_Recv( temp_loc_glo_index, proc_loc_int_size, MPI_INT, i, i, MPI_COMM_WORLD, &status );
            
            temp_data_values = (double *) malloc( proc_loc_int_size * sizeof(double));
            MPI_Recv( temp_data_values, proc_loc_int_size, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status );
            
            // Need calculation
            for ( int i = 0; i < proc_loc_int_size; i++ ){
                for ( int j = 0; j < 8; j++ ){
                    global_elems[ temp_loc_glo_index[i] * 8 + j ] = temp_elems[ i * 8 + j ];
                }
                data_values[ temp_loc_glo_index[i] ] = temp_data_values[i];
            }
        }
        
        vtk_write_unstr_grid_header(file_in, file_vtk_out, nintci, nintcf, points_count, points, global_elems);
        vtk_append_double(file_vtk_out, "CGUP", nintci, nintcf, data_values);
        free(global_elems);
        free(data_values);
        free(temp_elems);
        free(temp_loc_glo_index);
        free(temp_data_values);
    } else {
        // MPI_Send (&buf,count,datatype,dest,tag,comm)
        MPI_Send( &local_int_cells, 1, MPI_INT, 0, my_rank, MPI_COMM_WORLD );
        MPI_Send( elems, local_int_cells * 8, MPI_INT, 0, my_rank, MPI_COMM_WORLD );
        MPI_Send( local_global_index, local_int_cells, MPI_INT, 0, my_rank, MPI_COMM_WORLD );
        MPI_Send( cgup, local_int_cells, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD );
    }
    
    /*    char file_out[100];
     sprintf( file_out, "%d_%s",my_rank, file_vtk_out );*/
    /*    vtk_write_unstr_grid_header(file_in, file_vtk_out, nintci, nintcf, points_count, points, elems);
     vtk_append_double(file_vtk_out, "CGUP", nintci, nintcf, cgup);*/
    
    
    return 0;
}

int test_communication(char *file_in, char *file_vtk_out, int *local_global_index,
                       int local_num_elems, int neighbors_count, int* send_count, int** send_list,
                       int* recv_count, int** recv_list) {
    // Return an error if not implemented
    return -1;
}

