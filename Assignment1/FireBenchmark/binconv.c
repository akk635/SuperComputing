/*For converting the text file into binary file*/
#include <stdio.h>

int main( int argc, char *argv[] ){

	if( argc < 3 ){
			printf("Usage: %s input_file output_file\n", argv[0]);
			return EXIT_FAILURE;
	}

	// File pointers for opening the file streams
	FILE* text_fp, *binary_fp;

	// file name pointers
	char* text_file;
	char* binary_file;

	// Parameters for reading the text files
	int nintst, nintend, nextst, nextend;
	text_file = argv[1];
	binary_file = argv[2];

	// Opening the files for reading
	text_fp = fopen( text_file, "r" );
	binary_fp = fopen( binary_file, "wb" );

	// Error measures wrt opening the files
	if( text_fp == NULL ){
		printf("Error opening the file \n");
		return -1;
	}
	if( binary_fp == NULL ){
		printf("Error opening the binary file \n");
		return -1;
	}

	// Reading the volume dimensions form the text files
	fscanf( text_fp, "%d", &nintst);
	fwrite( &nintst, sizeof(int), 1, binary_fp );
	fscanf( text_fp, "%d", &nintend);
	fwrite( &nintend, sizeof(int), 1, binary_fp );
	fscanf( text_fp, "%d", &nextst);
	fwrite( &nextst, sizeof(int), 1, binary_fp );
	fscanf( text_fp, "%d", &nextend);
	fwrite( &nextend, sizeof(int), 1, binary_fp );

	// Reading the topological info (6*sizeof(int*))
	int buffer;
	for( int i = 0; i <= 5; i++ ){
		for( int j=nintst;j<=nintend;j++){
			fscanf( text_fp, "%d", &buffer );
			fwrite( &buffer, sizeof(int), 1, binary_fp );
		}
	}

	// Reading the eight miscellaneous arrays (boundary layers & sources) --double values
	double buff;
	for( int i = 0; i <= 7; i++ ){
		for( int j = nintst; j <= nintend; j++){
			fscanf( text_fp, "%lf", &buff);
			fwrite( &buff, sizeof(double), 1, binary_fp );
		}
	}

	// The red-black integer array is the last thing
	while( 1 ) {
		if( fscanf( text_fp, "%d", &buffer ) == EOF ){
			break;
		}
		fwrite( &buffer, sizeof(int), 1, binary_fp );
	}

	printf( "Wrote to the bin file" );
	fclose( text_fp );
	fclose( binary_fp );
	return 0;
}

