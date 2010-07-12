/***********************************************************************
* 1D.h                                                                 *
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

#ifndef MATRIZ1D_H_
#define MATRIZ1D_H_

#include <pthread.h>
#include <stdio.h>
#include <sys/types.h>
#include <malloc.h>
#include <math.h>
#include "commons.h"

/**
 * Función que se ejecuta en cada hilo para realizar el procesamiento
 * de la salida.
 * 
 * De acuerdo a los parámetros de partición, realiza la multiplicación 
 * de la partición que le corresponde. 
 * 
 */
void* runthread1D (void *param);

/**
 * Crea las particiones que cada hilo procesará. 
 * 
 * Se realiza una partición de los datos de salida. 
 * 
 * M: Cantidad de filas de la salida (filas de matriz A)
 * R: Cantidad de columnas de la salida (columnas de B)
 * q: Cantidad de Hilos que procesarán el resultado.
 * 
 * Obs: La partición i se define así: 
 *      p->id = i
 *      p->fini = i*salto
 *      p->ffin = (i+1)*salto - 1 
 */
void crear_particiones1D(int M, int R, int q, particion *p);

void imprimir_particiones1D(particion *p, int q);


#endif /*MATRIX1D2D_H_*/
