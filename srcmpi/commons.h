#ifndef COMMONS_H_
#define COMMONS_H_

#include <mpi.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#include <sys/time.h>
#include <dirent.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include "../../mmio/mmio.h"


// Conjunto de subrutinas para acceder a la matriz
// Utilizamos las matrices publicadas en math.nist.gov/MatrixMarket

// Algunas definiciones globales
#define MATRICES "../../matrices/"
#define FILEA "pA.mtx"
#define FILEB "pB.mtx"
#define FILEC "pC.mtx"


/* Macro para calcular el mínimo
 * de dos valores.
 */
#define min(a, b) (a < b ? a : b)

/* Marco para realizar el tiempo de ejecución de un proceso determinado
 * 
 */
#define MEDICION(f) {   \
   t1 = MPI_Wtime(),    \
   f,                   \
   t2 = MPI_Wtime(),    \
   t3 = t2 -t1;         \
}

// Tamaño del buffer para path
#define PATH_BUFF_SIZE 256

/* 
 * VARIABLES GLOBALES UTILIZADAS POR LOS PROCESOS
 * ===========================================================================
 */
float *MatrizA; 
float *MatrizB; 
float *MatrizC;
int genMat; 
int printResult; 
int jointResult;
int seqAlg;
int P; // cantidad de procesos
int Debug;

float *A; 
float *B; 
float *C;
int tamanho_bloque; 

double t1,t2,t3;

/*
 * Registro de identificación de cada hilo y que contiene la información
 * de que parte de la matriz procesa 
 */
typedef struct
{
  int       rank;              // identificador del proceso
                               // al que corresponde la particion
  
  int       CfilaIni;          // primera fila a procesar de C
  int       CcolIni;           // primera columna a procesar de C

  int       CfilaFin;          // última fila a procesar de C
  int       CcolFin;           // última columna a procesar de C
  
  int       AfilaIni;          // primera fila a procesar de A
  int       AfilaFin;          // última fila a procesar de A
  
  int       BcolIni;           // primera columna a procesar de B
  int       BcolFin;           // última columna a procesar de B
  
} bloque;

// OBS.: Se dejan los parámetros así para poder implementar en el futuro
//       algoritmos que soporten matrices no cuadradas. 

char *FA;       // path del archivo que contiene la matriz A
char *FB;       // path del archivo que contiene la matriz B
char *FC;       // path del archivo que contiene la matriz C
char *FCSV;     // path del archivo en que se escribirán los resultados de las métricas
int M,N,R;      // variables extras de uso global

// Funciones comunes a las dos implementaciones del algoritmo por hilos
void multiplicar_secuencial(float *A , float *B, float *C);
void matrizAleatoria(FILE *archivo, int Filas, int Columnas);
void inicializarMatriz(float *Matriz, int Filas, int Columnas);
float *asignar_mem_matriz(int Size);
/*
* Función que retorna el tiempo transcurrido
* desde Epoch (ver time(2)).
* La unidad de medida es "segundos" con
* decimales hasta la precisión de microsegundos.
*/
double tiempo_seg();

/*
* Función que retorna el tiempo transcurrido
* desde Epoch (ver time(2)).
* La unidad de medida es "milisegundos".
*/
long long tiempo_milis();

/* 
* Función que verificamos que el path del directorio no sea 
* muy largo.
* Se asume que entre el patron y el numero agregado para la
* creación archivos y directorios, no se superan los 50 
* caracteres.
*/
void verificar_longitud_path(char *path);

/*
* Función que se encarga de verificar la existencia
* de una barra final en el path de un directorio,
* Si no existe, retorna un caracter barra, en
* caso contrario retorna el caracter vacío.
*/
char *verificar_barra_final(char *path);

int procesar_parametros (int argc, char ** argv, char ** alg);
void uso(void);

void leer_matrices_aleatorio(float *A, float *B);
void leer_matrices();
void acumular_multiplicar(float A[tamanho_bloque][tamanho_bloque], 
                          float B[tamanho_bloque][tamanho_bloque], 
                          float C[tamanho_bloque][tamanho_bloque]);

void cerar_matriz(float *mat);
void cerar_bloque(float blo[tamanho_bloque][tamanho_bloque]);
void imprimirBloque(float bloq[tamanho_bloque][tamanho_bloque], int tam);

#endif /*COMMONS_H_*/
