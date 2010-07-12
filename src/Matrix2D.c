/***********************************************************************
* 2D.c                                                             *
* -------------------------------------------------------------------  *
*                                                                      *
* Implementación del algoritmo de multiplicación de matrices 1D        *
* con pthreads                                                         *
*                                                                      *
* Esta versión, particiona por Filas                                   *
*                                                                      *
* Autores:                                                             *
* - Cristhian Parra                                                    *
* - Fernando Mancía                                                    *
*                                                                      *
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


#include <pthread.h>
#include <stdio.h>
#include <sys/types.h>
#include <malloc.h>

// Conjunto de subrutinas para acceder a la matriz
// Utilizamos las matrices publicadas en math.nist.gov/MatrixMarket
#include "../../mmio/mmio.h"
#include "Matrix2D.h"

void* runthread2D (void *param)
{

	particion  *p = (particion *) param;
	//printf("Proceso %d iniciado...\n",p->id);
	int i,j,k;
	//printf("\n\nProceso %d calcula ",p->id);
	// calculo al multiplicación para las k (filaFin - filaIni) de A
	for ( i = p->filaIni; i <= p->filaFin; i++){

		// Multiplico fila i por N columnas de B
		for (j = p->colIni; j <= p->colFin; j++) {

			float Cij = 0.0;
			//printf(" C%d = ",(i * R) + j);
			for ( k = 0; k < N; k++) {
				float Aik = MatrizA[(i * N) + k];
				float Bkj = MatrizB[(k * R) + j];
				Cij = Cij + (Aik * Bkj);
				//printf("(A%d * B%d) + ", (i * N) + k , (k * R) + j );
			}
			//printf("\n");
			MatrizC[(i * R) + j] = Cij;
		}
	}

	return NULL;
}
 
void crear_particiones2D(int M, int R, int q, particion *p) {

       double step_float;
       int stepFila = 0,stepColumna = 0, i,j ,raiz_q;

       raiz_q = (int) sqrt( (double) q );

       // calculo del salto de índice para filas según la cantidad de hilos
       step_float = ( (float) M ) / ((float)raiz_q);
       stepFila = (int)rint(step_float); 

       // calculo del salto de índice para filas según la cantidad de hilos
       step_float = ( (float) R ) / ((float)raiz_q);
       stepColumna = rint(step_float);

       //printf("\n* step_columna = %d * step_fila= %d \n", stepColumna , stepFila );

       int filatmpIni = 0, filatmpFin = 0, count = 0;
       for (i=0; i< raiz_q; i++) {
           filatmpIni = stepFila * i;
           filatmpFin = (stepFila  * (i + 1)) -1;
           if(i == raiz_q -1 ){
               filatmpFin = M -1;
           }

           for (j=0; j< raiz_q; j++) {

               p[count].id = count;
               p[count].filaIni = filatmpIni ;
               p[count].filaFin = filatmpFin ;
               p[count].colIni = stepColumna * j;
               p[count].colFin = (stepColumna * (j + 1)) -1;
               if(j == raiz_q -1){
                   p[count].colFin = R -1;
               }
               count++;

          }
       }
}


void imprimir_particiones2D(int q, particion *p) {
	int i;
	for (i=0; i< q; i++) {
            printf("\n Particion %d filas (%d a %d) columnas (%d a %d)", p[i].id, p[i].filaIni, p[i].filaFin, p[i].colIni, p[i].colFin);
       }
}
