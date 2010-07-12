#include "salida-metricas.h"


FILE *abrir_archivo(char *path, char *modo) {
	FILE *arch;
	
	if ((arch = fopen(path, modo)) == NULL) {
		fprintf(stderr, 
				"abrir_archivo(): Error abriendo archivo '%s' en modo '%s'\n", 
				path, modo);
		perror("sys_msg");
		exit(1);
	}
	
	return arch;
}

FILE *open_csvFile(char *filename) {
	return abrir_archivo(filename, "w");
}

FILE *open_csvFile_Append(char *filename) {
	return abrir_archivo(filename, "a");
}

void close_csvFile(FILE *fd) {
	fclose(fd);
}

void print_header(FILE *csvFile,char *fileno) {
	
	struct stat info;
	
	stat(fileno, &info);
	
	if(!info.st_size) {
		fprintf(csvFile,"Hilos,Tamaño,Tiempo,Overhead\n");			
		fflush(csvFile);	
	}	
}

void print_test_result(FILE *csvFile, result_t *r) {
	fprintf(csvFile,"%s,%d,%d,%lld,%lld\n", 
			r->testId, r->q, r->size, r->T, r->O);
	
	fflush(csvFile);
}
