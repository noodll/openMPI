#ifndef SALIDAMETRICAS_H_
#define SALIDAMETRICAS_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>

/* 
 * Estructura donde se almacenan los resultados de una prueba
 */
typedef struct {
	char testId[20];   /* Identificador de la prueba 			*/
	int q;	     	/* Cantidad de Hilos utilizados         */
	int size;	    /* Tamaño del problema procesado (MxR)  */
	long long T;	/* Tiempo ejecución total de la multip. */ 
	long long O; 	/* Overhead de creación de hilos y parts*/
} result_t;

/* 
 * Operaciones para abrir y cerrar el archivo CSV de resultados 
 */

FILE *abrir_archivo(char *path, char *modo);
FILE *open_csvFile(char *filename);
FILE *open_csvFile_Append(char *filename);
void close_csvFile(FILE *fd);

/*
 * Imprime la cabecera del archivo de salida. 
 * 
 * Se ejecuta solo una vez. 
 */
void print_header(FILE *csvFile,char *fileno);

/* Imprime una línea de resultados en el archivo de resultados
 * 
 * Se ejecuta cada vez que se termina una operación, luego de 
 * calcular los valores de las métricas correspondientes. 
 */
void print_test_result(FILE *csvFile, result_t *p);


#endif /*SALIDAMETRICAS_H_*/
