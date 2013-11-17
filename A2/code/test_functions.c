/**
 * Functions to test the data distribution and communication lists creation algorithms
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */
#include "test_functions.h"
#include <assert.h>

int test_distribution(char *file_in, char *file_vtk_out, int *local_global_index,
                      int **global_local_index, int nintci, int nintcf, int points_count,
                      int **points, int *elems, int local_int_cells, double *cgup) {
    // Return an error if not implemented
    //Checking for the ranks of the elements

    MPI_Status status;
    int my_rank, nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);    /// get number of processes

/*    if (my_rank == 0){
        int proc_local_size;
        int *temp_elems = NULL;

        // Gather all the elems into one processor
        int *global_elems = (int *) malloc((nintcf - nintci + 1) * 8 * sizeof(int));
        // For global indexing
        global_elems = global_elems - nintci * 8;

        for (int i = nproc - 1; i > 0; i--){
            // MPI_Recv (&buf,count,datatype,source,tag,comm,&status)
            MPI_Recv( &proc_local_size, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status );
            temp_elems = (int *) malloc( proc_local_size * 8 * sizeof(int));
            MPI_Recv( temp_elems, proc_local_size * 8, MPI_INT, i, i, MPI_COMM_WORLD, &status );
            // Now I even need local_global_indices
        }

    } else {
        // MPI_Send (&buf,count,datatype,dest,tag,comm)
        MPI_Send( &local_int_cells, 1, MPI_INT, 0, my_rank, MPI_COMM_WORLD );

    }*/


/*    for ( int i = nintci; i <= nintcf; i++ ){
        global_elems[i] = (int*) malloc ( 8 * sizeof(int));
    }*/

/*    for (int i = 0 ; i < local_int_cells; i++){
        for (int j = 0; j < 8; j++){
            global_elems[ local_global_index[i] ][ j ] = elems[ i * 8 + j ];
        }
    }*/
    char file_out[100];
    sprintf( file_out, "%d_%s",my_rank, file_vtk_out );

    vtk_write_unstr_grid_header(file_in, file_out, nintci, nintcf, points_count, points, elems);
    vtk_append_double(file_out, "RANK", nintci, nintcf, cgup);
    return -1;
}

int test_communication(char *file_in, char *file_vtk_out, int *local_global_index,
                       int local_num_elems, int neighbors_count, int* send_count, int** send_list,
                       int* recv_count, int** recv_list) {
    // Return an error if not implemented
    return -1;
}

