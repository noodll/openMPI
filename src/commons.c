#include "commons.h"

void multiplicar_secuencial(float *A , float *B, float *C) {
    int i, j, k =0;

	for (i=0; i<M; i++) {
		for (j=0; j<R; j++) {
			C[(i*R)+j] = 0;
			//printf("\n C%d = ", (i*R)+j);
			for (k=0; k<N; k++) {
				C[(i*R)+j] = C[(i*R)+j] + (A[(i*N)+k] * B[(k*R)+j]);
				//printf(" (A%d * B%d) +", (i*N)+k, (k*R)+j);
			}
			//printf(" = %f ", C[(i*R)+j]);
		}
	}
}

void matrizAleatoria(FILE *archivo, int Filas, int Columnas) {

	int i,j;
	float random_value;
	char *linea = malloc(sizeof(char) * 20);
	char *cabecera = malloc(sizeof(char) * 20);

	fprintf(archivo, "%%%MatrixMarket matrix coordinate real general\n");
	int cant = Filas * Columnas;
	sprintf(cabecera, "%d %d %d\n", Filas, Columnas, cant);
	fprintf(archivo, cabecera);
	
	for (i=0; i<Filas; i++) {
		for (j=0; j<Columnas; j++) {
			int mivalor =1+(int) (100.0*rand()/(RAND_MAX+1.0)); 
			sprintf(linea, "%d %d %d\n", i,j, mivalor);
			fprintf(archivo, linea);
		}
	}	
	free(linea);
	free(cabecera);
}

void inicializarMatriz(float *Matriz, int Filas, int Columnas) {
  int i;
  for ( i = 0; i < Filas*Columnas; ++i) {
        Matriz[i] = 0.0;
  }  
} 

struct timeval obtener_tiempo() {
	struct timeval tv;

	if (gettimeofday(&tv, NULL) == -1) {
		fprintf(stderr, "medir_tiempo(): Error al obtener el tiempo\n");
		perror("sys_msg");
		exit(1);
	}
	
	return tv;
}

double tiempo_seg() {
	double seg, mseg;
	struct timeval tv = obtener_tiempo();
	
	seg = tv.tv_sec;
	mseg += (tv.tv_usec / 1000000.0);
	
	return (seg + mseg);
}

long long tiempo_milis() {
	long long seg, mseg;
	struct timeval tv = obtener_tiempo();
	
	// Guardamos los tiempos en long long int
	seg  = tv.tv_sec;
	mseg = tv.tv_usec;
	
	// Convertimos los tiempos en milisegundos
	seg  = seg * 1000;
	mseg = mseg / 1000;
	
	// Retornamos el resultado
	return (seg + mseg);
}

void verificar_longitud_path(char *path) {
	int maxlen = strlen(path);
	
	if (maxlen > PATH_BUFF_SIZE) {
		fprintf(stderr, "Path '%s' demasiado largo (%d, max=%d)",
				path, strlen(path), 200);
		exit(1);
	}
}

char *verificar_barra_final(char *path) {
	int len = strlen(path);
	
	if (path[len - 1] == '/')
		return "";
	
	return "/";
}
