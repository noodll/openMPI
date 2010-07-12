/* 
 * File:   2dblock.h
 * Author: cparra
 *
 * Created on 10 de julio de 2008, 05:48 PM
 */

#ifndef _2DBLOCK_H
#define	_2DBLOCK_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <mpi.h>
#include <stdio.h>
#include <sys/types.h>
#include <malloc.h>
#include <math.h>

#include "commons.h"

#define TAGMULT 100

    int SubBlockSize;  /* Define el tamaño de los bloques asignados a las tareas
                        * ---> Las tareas van de 0 ... N, y a cada uno la 
                        *      asignamos lo necesario para computar un Bloque
                        *      de tamaño BlockSize x BlockSize de C. 
                        *
                        *      Este parámetro es igual a sqrt(N)
                        */
    
    int ProcBlockSize; /* Define el tamaño de los bloques de Proceso, en los 
                        * cuales, las tareas se distribuirán de manera cíclica. 
                        * ---> Las proceso van de 0 ... P, y a cada uno la 
                        *      asignamos el cómputo de una tarea de manera cíclica
                        *      
                        *      Es igual a sqrt(P)
                        */
    
    int indiceTarea;   /* Es un número que indica que tarea estamos procesando. 
                        * ---> Va de 0 a N
                        */
        

    /*
     * Función que ejecuta la lógica del algoritmo 2D con particionamiento por
     * Bloques cíclicos.
     *
     * prank    : id del proceso que se está ejecutando. 
     * N        : cantidad de filas, columnas, de las matrices cuadradas
     * A        : Matriz A completa (solo prank = 0, en el resto es null)
     * B        : Matriz B completa (solo prank = 0, en el resto es null)
     * C        : Matriz C de salida (solo prank = 0, en el resto es null)
     */
    void ejecutar_2d(int prank);

    /*
     * Carga el Buffer que se le envía a un proceso con las filas y columnas
     * correspondientes al Bloque que tiene que procesar de C. 
     *
     * El buffer, contiene de manera continua, sqrt(P) filas de A, seguido de 
     * sqrt(P) columnas de B. 
     *
     * buffer: buffer de salida donde se guardan los datos a enviar
     * A     : matriz A
     *
     */
    void cargar_buffer_envio(float *buffer, int fini, int cantFilas, int cini, int cantCols, int BuffCut);
    
    void cargar_resultado_C(float *buffer, int fini, int cantFilas, int cini, int cantCols);
    
    void multiplicar_bloque(float *buffer, float *Salida, int BuffCut, int prank,int tarea); 
    
    void imprimir_buffer(float *buffer);
    
#ifdef	__cplusplus
}
#endif

#endif	/* _2DBLOCK_H */

