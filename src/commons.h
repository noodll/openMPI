#ifndef COMMONS_H_
#define COMMONS_H_

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
   t1 = tiempo_milis(), \
   f,                   \
   t2 = tiempo_milis(), \
   t3 = t2 -t1;          \
}


// Tamaño del buffer para path
#define PATH_BUFF_SIZE 256

/* 
 * VARIABLES GLOBALES QUE REPRESENTAN LA MEMORIA COMPARTIDA DE TODOS LOS HILOS
 * ===========================================================================
 */
float *MatrizA; 
float *MatrizB; 
float *MatrizC;

// Parámetros de la matriz A
int M;       // cantidad de filas
int N;       // cantidad de columnas

//  Parámetros de la matriz B
int S;       // cantidad de filas --> debe ser igual a N
int R;       // canidad de columnas -->



/*
 * Registro de identificación de cada hilo y que contiene la información
 * de que parte de la matriz procesa 
 */
typedef struct
{
  int       id;              // identificador del hilo=ij= i*10+j (hiloi,j)
                             // al que corresponde la particion
  int       filaIni;         // primera fila a procesar
  int       colIni;          // primera columna a procesar

  int       filaFin;         // primera fila a procesar
  int       colFin;          // primera columna a procesar
} particion;


// Funciones comunes a las dos implementaciones del algoritmo por hilos
void multiplicar_secuencial(float *A , float *B, float *C);
void matrizAleatoria(FILE *archivo, int Filas, int Columnas);
void inicializarMatriz(float *Matriz, int Filas, int Columnas);

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




#endif /*COMMONS_H_*/
