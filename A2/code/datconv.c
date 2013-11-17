#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv){
        // File pointers for opening the file streams
        FILE* text_fp, *binary_fp;

        // file name pointers
        char* text_file;
        char* binary_file;
	binary_file = argv[1];
	text_file = argv[2];

	text_fp = fopen( text_file, "w");
	binary_fp = fopen( binary_file, "rb");

	int buffer[4];
	for (int i=0; i< 4; i++	){
		fread(buffer+i, sizeof(int), 1, binary_fp);
		fprintf( text_fp, "%d \t", buffer[i]);
	}
fprintf(text_fp, "\n");
int fpcount = 0;
int tot_freads = (buffer[1]-buffer[0]+1);
int temp_buffer;
while( fpcount < tot_freads ){	
	for (int i = 0; i<6;i++){
		fread(&temp_buffer, sizeof(int), 1, binary_fp);
		fprintf( text_fp, "%d \t", temp_buffer);	
	}
fpcount++;
fprintf(text_fp, "\n" );
}
	fclose(binary_fp);
	fclose(text_fp);
return 0;
}
