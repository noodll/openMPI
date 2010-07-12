/***********************************************************************
* 2D.h                                                                 *
* -------------------------------------------------------------------  *
*                                                                      *
* <CABECERA> Implementación del algoritmo de multiplicación de         *
* matrices 1D con pthreads                                             *
*                                                                      *
* Esta versión, particiona por Filas                                   *
*                                                                      *
* Autores:                                                             *
* - Cristhian Parra                                                    *
* - Fernando Mancía                                                    *
*                                                                      *
* NOTA: por el momento, solo implementamos el algoritmo 1D que realiza *
* la partición de las filas en cada proceso                            *
*                                                                      *
* * Sugerencia, utilizar un archivo de configuración para colocar los  *
*   parámetros.                                                        *
*                                                                      *
* * PRUEBAS                                                            *
*   - Probar varias veces con un tamaño --> sacar promedios            *
*   - Que pasa cuando se aumenta la cantidad de hilos                  *
*                                                                      * 
* * Fundamentar el problema de la memoria, como se podría resolver     *
*                                                                      *
* CONCLUSIONES                                                         *
* 1. Sacar conclusiones fundamentadas en las pruebas                   *
*                                                                      *
************************************************************************/
#ifndef MATRIX2D_H_
#define MATRIX2D_H_

#include <pthread.h>
#include <stdio.h>
#include <sys/types.h>
#include <malloc.h>
#include <math.h>

// Conjunto de subrutinas para acceder a la matriz
// Utilizamos las matrices publicadas en math.nist.gov/MatrixMarket
#include "../../mmio/mmio.h"
#include "commons.h"


/* 
 * VARIABLES GLOBALES QUE REPRESENTAN LA MEMORIA COMPARTIDA DE TODOS LOS HILOS
 * ===========================================================================
 */
/* 
 * Código de cada hilo. De acuerdo a los parámetros de partición, realiza
 * la multiplicación de la partición que le corresponde. 
 * 
 */
void * runthread2D (void *param);

void crear_particiones2D(int M, int R, int q, particion *p);

void imprimir_particiones2D(int q, particion *p);

#endif /*MATRIX2D_H_*/
