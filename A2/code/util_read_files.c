/**
 * Helper functions for reading from input data file
 *
 * @author E. Xue, V. Petkov
 * @date 22-May-2009, 22-Oct-2012
 */
#include <stdio.h>
#include <stdlib.h>
#include "util_read_files.h"
#include <assert.h>

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
int read_binary_geo(char *file_name, int *NINTCI, int *NINTCF, int *NEXTCI, int *NEXTCF, int ***LCC,
                    double **BS, double **BE, double **BN, double **BW, double **BL, double **BH, double **BP,
                    double **SU, int* points_count, int*** points, int** elems, int **local_global_index,
                    int *elemcount, int *local_int_cells, int ***global_local_index ) {
    int i = 0;
    int my_rank, nproc;
    MPI_Status status;
    FILE *fp = fopen(file_name, "rb");

    if ( fp == NULL ) {
        fprintf(stderr, "Error opening file %s\n", file_name);
        return -1;
    }

    // 4 variables in total!!!
    fread(NINTCI, sizeof(int), 1, fp);
    fread(NINTCF, sizeof(int), 1, fp);
    fread(NEXTCI, sizeof(int), 1, fp);
    fread(NEXTCF, sizeof(int), 1, fp);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    int tot_domain_cells = *NEXTCF - *NINTCI + 1;

    // For storing ranks according to the indices
    int *distr_buffer = (int *) malloc( sizeof(int) * tot_domain_cells );

    if (my_rank == 0){
        int fpcount;
        fpcount = *NINTCI;

        // So that we can call on explicit coordinates
        distr_buffer = distr_buffer - *NINTCI;

        int local_cells_size;

        int normal_local_size = tot_domain_cells / nproc;
        int res_cells =  tot_domain_cells % nproc;

        // Initializing the distribution array with -1
        for ( int i = *NINTCI; i < *NEXTCF + 1; i++){
                distr_buffer[i] = -1;
        }

        // All the neighbors do not cover all the cells.
        // there are some external cells which are not neighbors of internal
        int temp_cells_size = 0;

        // To get the indices to write in the distr_buffer
        int temp_buffer;

        // equally distributing the local size for each processor
        for( int i = nproc-1; i >= 0 ; i-- ){
            *elemcount = 0;
            local_cells_size = normal_local_size;
            // Reading the topological info and then distributing according to the locality
            if ( res_cells > 0 ){
                res_cells--;
                local_cells_size++;
            }

            while ((*elemcount) < local_cells_size) {

                if( fpcount == *NINTCF + 1){
                    break;
                }

                if( distr_buffer[ fpcount ] == -1){
                    distr_buffer[ fpcount ] = i;
                    (*elemcount)++;
                    temp_cells_size++;
                }

                for ( int j = 0;j < 6; j++ ){
                    fread( &temp_buffer, sizeof(int), 1, fp );
                    if ( distr_buffer[ temp_buffer ] == -1){
                        distr_buffer[ temp_buffer ] = i;
                        (*elemcount)++;
                        temp_cells_size++;
                    }
                }
                fpcount++;
            }

            // Sending the internal cells indices to the respective processes
            // MPI_Send (&buf,count,datatype,dest,tag,comm)
            if( i != 0 ){
                MPI_Send( elemcount, 1, MPI_INT, i, i, MPI_COMM_WORLD );
            }
        }

        // Now distributing all the remaining external cells to process 0
        for( int i = *NEXTCI; i <= *NEXTCF; i++){
            if( distr_buffer[ i ] == -1){
                distr_buffer[ i ] = 0;
                (*elemcount)++;
                temp_cells_size++;
            }
        }
        assert( temp_cells_size == tot_domain_cells );
        // Now distributing the buffer to all the processors
        for( int i = 1; i < nproc; i++ ){
            // MPI_Send (&buf,count,datatype,dest,tag,comm)
            MPI_Send( distr_buffer + *NINTCI, tot_domain_cells, MPI_INT, i, i, MPI_COMM_WORLD );
        }
    } else {
        MPI_Recv( elemcount, 1, MPI_INT, 0, my_rank, MPI_COMM_WORLD, &status );
        // MPI_Recv (&buf,count,datatype,source,tag,comm,&status)
        MPI_Recv( distr_buffer, tot_domain_cells, MPI_INT, 0, my_rank, MPI_COMM_WORLD, &status );
    }

    if( my_rank != 0){
        distr_buffer = distr_buffer - *NINTCI;
    }

    (*global_local_index) = (int **) malloc( (*NEXTCF - *NINTCI + 1) * sizeof(int *) );
    (*global_local_index) = (*global_local_index) - *NINTCI;

    for( int i = *NINTCI; i <= *NEXTCF; i++ ){
        (*global_local_index)[i] = (int *) calloc( 2,  sizeof(int));
    }

    *local_int_cells = 0;

    for (int i = *NINTCI; i <= *NINTCF ; i++){
        if( distr_buffer[i] == my_rank ){
            (*local_int_cells)++;
        }
    }

    *local_global_index =  malloc( (*elemcount) * sizeof(int));

    int j = 0;
    for (int i = *NINTCI; i <= *NEXTCF ; i++){
        (*global_local_index)[i][0] = distr_buffer[i];
        if( distr_buffer[i] == my_rank ){
            (*local_global_index)[j] = i;
            (*global_local_index)[i][1] = j;
            j++;
        }
     }
    int check;
    check = distr_buffer[*NEXTCF];

    assert( j == (*elemcount) );
    j = 0;

    // Allocating the LCC in individual processes
    if ( (*LCC = (int**) malloc( (*local_int_cells) * sizeof(int*))) == NULL ) {
            fprintf(stderr, "malloc failed to allocate first dimension of LCC");
            return -1;
    }

    for ( i = 0; i < (*local_int_cells); i++ ) {
        if ( ((*LCC)[i] = (int *) malloc(6 * sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc failed to allocate second dimension of LCC\n");
            return -1;
        }
    }

    int index_read = 4 * sizeof(int);

    // Start reading LCC
    for (int i = 0; i < (*local_int_cells); i++){
        fseek( fp, index_read + ( (*local_global_index)[i] - *NINTCI )* 6 * sizeof(int), SEEK_SET );
        for ( int j = 0; j < 6; j++){
            fread(&(*LCC)[i][j], sizeof(int), 1, fp);
        }
    }

    // allocate other arrays
    if ( (*BS = (double *) malloc( (*elemcount) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BS) failed\n");
        return -1;
    }

    if ( (*BE = (double *) malloc( (*elemcount) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BE) failed\n");
        return -1;
    }

    if ( (*BN = (double *) malloc( (*elemcount) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BN) failed\n");
        return -1;
    }

    if ( (*BW = (double *) malloc( (*elemcount) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BW) failed\n");
        return -1;
    }

    if ( (*BL = (double *) malloc( (*elemcount) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BL) failed\n");
        return -1;
    }

    if ( (*BH = (double *) malloc( (*elemcount) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BH) failed\n");
        return -1;
    }

    if ( (*BP = (double *) malloc( (*elemcount) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BP) failed\n");
        return -1;
    }

    if ( (*SU = (double *) malloc( (*elemcount) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(SU) failed\n");
        return -1;
    }

    // read the other arrays
    int lcc_read_end = ( ( *NINTCF - *NINTCI + 1 ) * 6 + 4 ) * sizeof( int );

    if (my_rank == check){
        long int fpos = ftell(fp);
        assert( fpos == lcc_read_end );
    }

    for ( i = 0; i < (*local_int_cells); i++ ) {
        fseek( fp, lcc_read_end + ( (*local_global_index)[i] - *NINTCI ) * 8 * sizeof(double), SEEK_SET );
        fread(&((*BS)[i]), sizeof(double), 1, fp);
        fread(&((*BE)[i]), sizeof(double), 1, fp);
        fread(&((*BN)[i]), sizeof(double), 1, fp);
        fread(&((*BW)[i]), sizeof(double), 1, fp);
        fread(&((*BL)[i]), sizeof(double), 1, fp);
        fread(&((*BH)[i]), sizeof(double), 1, fp);
        fread(&((*BP)[i]), sizeof(double), 1, fp);
        fread(&((*SU)[i]), sizeof(double), 1, fp);
    }

    int coe_read_end = lcc_read_end + 8 * ( *NINTCF - *NINTCI + 1 ) * sizeof(double);

    if (my_rank == check){
        long int fpos = ftell(fp);
        assert( fpos == coe_read_end );
    }

    // read geometry i.e nodes
    // allocate elems
    if ( (*elems = (int*) malloc( (*local_int_cells) * 8 * sizeof(int))) == NULL ) {
        fprintf(stderr, "malloc failed to allocate elems");
        return -1;
    }

    // read elems
    for ( i = 0; i < (*local_int_cells) ; i++ ) {
            fseek( fp, coe_read_end + ( (*local_global_index)[i] - *NINTCI ) * 8 * sizeof(int), SEEK_SET );
        for (int j = 0; j < 8; j++){
            fread(&((*elems)[( i * 8 ) + j]), sizeof(int), 1, fp);
        }
    }

    int nodes_read_end = coe_read_end + 8 * ( *NINTCF - *NINTCI + 1 ) * sizeof(int);

    if ( my_rank == check ){
        long int fpos = ftell(fp);
        assert( fpos == nodes_read_end );
    }

    fseek( fp, nodes_read_end, SEEK_SET );

    fread(points_count, sizeof(int), 1, fp);

    if( my_rank == 0 ){
        printf( "points count : %d \n", *points_count );
    }

    // allocate points vec
    if ( (*points = (int **) calloc(*points_count, sizeof(int*))) == NULL ) {
        fprintf(stderr, "malloc() POINTS 1st dim. failed\n");
        return -1;
    }

    for ( i = 0; i < *points_count; i++ ) {
        if ( ((*points)[i] = (int *) calloc(3, sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc() POINTS 2nd dim. failed\n");
            return -1;
        }
    }

    int coordIdx;
    int pointIdx;
    for ( pointIdx = 0; pointIdx < *points_count; pointIdx++ ) {
        for ( coordIdx = 0; coordIdx < 3; coordIdx++ ) {
            fread(&((*points)[pointIdx][coordIdx]), sizeof(int), 1, fp);
        }
    }


/*    // allocating LCC
    if ( (*LCC = (int**) malloc((*NINTCF + 1) * sizeof(int*))) == NULL ) {
        fprintf(stderr, "malloc failed to allocate first dimension of LCC");
        return -1;
    }

    for ( i = 0; i < *NINTCF + 1; i++ ) {
        if ( ((*LCC)[i] = (int *) malloc(6 * sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc failed to allocate second dimension of LCC\n");
            return -1;
        }
    }

    // start reading LCC
    // Note that C array index starts from 0 while Fortran starts from 1!
    for ( i = (*NINTCI); i <= *NINTCF; i++ ) {
        fread(&(*LCC)[i][0], sizeof(int), 1, fp);
        fread(&(*LCC)[i][1], sizeof(int), 1, fp);
        fread(&(*LCC)[i][2], sizeof(int), 1, fp);
        fread(&(*LCC)[i][3], sizeof(int), 1, fp);
        fread(&(*LCC)[i][4], sizeof(int), 1, fp);
        fread(&(*LCC)[i][5], sizeof(int), 1, fp);
    }

    long int temp;
    temp = ftell(fp);
    printf( "fpos : %ld \n", temp);

    // allocate other arrays
    if ( (*BS = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BS) failed\n");
        return -1;
    }

    if ( (*BE = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BE) failed\n");
        return -1;
    }

    if ( (*BN = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BN) failed\n");
        return -1;
    }

    if ( (*BW = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BW) failed\n");
        return -1;
    }

    if ( (*BL = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BL) failed\n");
        return -1;
    }

    if ( (*BH = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BH) failed\n");
        return -1;
    }

    if ( (*BP = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BP) failed\n");
        return -1;
    }

    if ( (*SU = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(SU) failed\n");
        return -1;
    }

    // read the other arrays
    for ( i = (*NINTCI); i <= *NINTCF; i++ ) {
        fread(&((*BS)[i]), sizeof(double), 1, fp);
        fread(&((*BE)[i]), sizeof(double), 1, fp);
        fread(&((*BN)[i]), sizeof(double), 1, fp);
        fread(&((*BW)[i]), sizeof(double), 1, fp);
        fread(&((*BL)[i]), sizeof(double), 1, fp);
        fread(&((*BH)[i]), sizeof(double), 1, fp);
        fread(&((*BP)[i]), sizeof(double), 1, fp);
        fread(&((*SU)[i]), sizeof(double), 1, fp);
    }

    // read geometry
    // allocate elems
    if ( (*elems = (int*) malloc((*NINTCF + 1) * 8 * sizeof(int))) == NULL ) {
        fprintf(stderr, "malloc failed to allocate elems");
        return -1;
    }

    // read elems
    for ( i = (*NINTCI); i < (*NINTCF + 1) * 8; i++ ) {
        fread(&((*elems)[i]), sizeof(int), 1, fp);
    }

    fread(points_count, sizeof(int), 1, fp);

    printf( "points_count : %d \n", points_count );

    // allocate points vec
    if ( (*points = (int **) calloc(*points_count, sizeof(int*))) == NULL ) {
        fprintf(stderr, "malloc() POINTS 1st dim. failed\n");
        return -1;
    }

    for ( i = 0; i < *points_count; i++ ) {
        if ( ((*points)[i] = (int *) calloc(3, sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc() POINTS 2nd dim. failed\n");
            return -1;
        }
    }

    int coordIdx;
    int pointIdx;
    for ( pointIdx = 0; pointIdx < *points_count; pointIdx++ ) {
        for ( coordIdx = 0; coordIdx < 3; coordIdx++ ) {
            fread(&((*points)[pointIdx][coordIdx]), sizeof(int), 1, fp);
        }
    }*/
    fclose(fp);

    return 0;
}


