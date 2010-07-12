/***********************************************************************
* 1D.c                                                             *
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

#include "matriz1d.h"


void* runthread1D (void *param) {

  /* 
   * Recorrer A desde la filaIni hasta la filaFin y procesar el cálculo de
   * la matriz C de salida: 
   *     Desde C[filaIni,*] hasta C[filaFin,*]
   *
   */
   
   
  particion  *p = (particion *) param;
  
  //printf("Proceso %d iniciado...\n",p->id);
  int i,j,k;
  float Aik,Bkj,Cij;
  
  // calculo la multiplicación para las k (filaFin - filaIni) de C
  for ( i = p->filaIni; i <= p->filaFin; i++){
	// por las N columnas de B  
      for ( j = 0; j < R; j++) {
    	  Cij = 0.0;
    	  for (k = 0; k < N; k++) {       

    		  Aik = MatrizA[(i * N) + k];
    		  Bkj = MatrizB[(k * R) + j];
    		  
    		  Cij = Cij + Aik * Bkj;
                
    		  //printf("(%d,%d,%d)(%d,%d,%d)Multiplicando %f por %f = %f (C en %d,%d)...\n",i,j,k,M,N,R,Aik,Bkj,Cij,i,j );
          }
    	  
    	  MatrizC[(i * R) + j] = Cij;
      }
  }
  
  return NULL;
}

void crear_particiones1D(int M, int R, int q, particion *p) {

	float step_float;
	int step,i,filaFinTmp;
	
	// calculo del salto de índice según la cantidad de hilos
	step_float = ( (float) M ) / q;	
	step = rintf(step_float);         
	
	for (i=0; i< q; i++) {
		
		p[i].id = i;
		p[i].filaIni = i*step;
		filaFinTmp = (i+1)*step - 1;
		
		if (i==q-1) {
			p[i].filaFin = M-1;
		} else {
			p[i].filaFin = filaFinTmp;
		}
	}	
}

void imprimir_particiones1D(particion *p, int q) {
	int i;
	for (i=0; i< q; i++) {		
		printf("\n Particion %d filas (%d a %d) ",p[i].id,p[i].filaIni,p[i].filaFin);
	}	
}
