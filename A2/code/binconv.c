#include <stdio.h>

int main( int argc, char **argv ) {
    // File pointers for opening the file streams
    FILE* text_fp, *binary_fp;

    // file name pointers
    char* text_file;
    char* binary_file;
    text_file = argv[1];
    binary_file = argv[2];
    // Opening the files for reading
    text_fp = fopen( text_file, "r" );
    binary_fp = fopen( binary_file, "wb" );
    int buffer;
    while ( fscanf( text_fp, "%d", &buffer ) != EOF ) {
        fwrite( &buffer, sizeof(int), 1, binary_fp );
    }

};

