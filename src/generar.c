
#include "commons.h"


int main(int argc, char ** argv){
	if(argc != 4){
		printf("\nDebe especificar: FILAS COLUMNAS PATH_MATRIZ\n");
		exit(1);
	}

	int filas = atoi(argv[1]);
	int col = atoi(argv[2]);
	FILE * archivo = fopen(argv[3], "r");
	if(archivo != NULL){
		fclose(archivo);
		printf("\nEl archivo %s ya existe.\n", argv[3]);
		exit(1);		
	}
	archivo = NULL;
	archivo = fopen(argv[3], "w+");
	matrizAleatoria(archivo, filas, col);
	fclose(archivo);
	printf("\nMatriz creada exitosamente en el archivo %s\n", argv[3]);
	return 0;
}
