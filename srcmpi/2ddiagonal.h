/* 
 * File:   2ddiagonal.h
 * Author: cparra
 *
 * Created on 20 de julio de 2008, 01:56 PM
 */

#ifndef _2DDIAGONAL_H
#define	_2DDIAGONAL_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <malloc.h>
#include <mpi.h>
#include "commons.h"
#define TAGMULT 100


    /** 
     * Variables Globales del Algoritmo. 
     *
     */

    MPI_Comm Comm2d; // comunicador de la topología 2d.
    MPI_Comm commFilas;
    MPI_Comm commColumnas; 

    int Grupos; // Cantidad de Grupos de filas y columnas distribuidas.
    // en la diagonal. 

    int GrupSize; // Define la cantidad de filas o columnas del grupo.
    // Es igual que Grupos, pero se define así para diferenciar.

    float *IthGrupoA; // Buffer para el iésimo grupo de columnas de A. (A*,i)
    float *IthGrupoB; // Buffer para el iésimo grupo de filas de B. (Bi,*)
    float *BloqueB;   // Bloque de B a procesar en el producto externo (Bi,k)
    float *BloquesB;  // Bloques de B serializados para el scatter. (Concatenación de todos los Bi,k en Pii) 
    float *IProduct;  // Producto externo calculado por cada proceso. (I*,i)
    float *IProductReduced; // Producto externo recibido durante el Reduce
    
    float *BloqueC; // Bloque C*,i reducido en cada proceso diagonal. 

    int *Tags; // Array que almacena el nro. de secuencia actual para etiquetar
    // los mensajes enviados/recibidos por un proceso

    int Dims[2];
    int Periods[2];
    int mycoords[2];
    int myrank;

    /**
     * Funciones principales del Algoritmo
     *
     */
    /*
     * Realiza la distribución inicial de los grupos de filas y columnas desde 
     * el proces P0,0
     */
    void distribucion_inicial();

    /* 
     * Ejecuta la lógica principal del algoritmo, una vez distribuidas los grupos
     * iniciales. 
     *    - Broadcast de A*,j a todos los procesos p*,j desde la diagonal
     *    - Send de Bik a cada pkj                                       
     *    - Recibir A*,j y Bji                                          
     *    - Calcular producto externo I*,i = A*,j x Bji                  
     *    - Enviar I*,i en la dirección i hasta Pi,i.                    
     *    - En cada Pi,i, recibir I*,i y calcular C*,i = C*,i + I*,i     
     */
    void ejecutar2dd();

    /*
     * Carga en Buffer, en el formato serializado de filas (una fila detrás de 
     * la otra), el contenido de origen. 
     *
     * Exactamente se copian cantFilas desde la fila fini, y cantColumnas desde
     * la columna cini
     */
    void cargar_buffer_grupo(float *buffer, float *origen, int size, int fini, int cantFilas, int cini, int cantCols);

    /*
     * Carga el Iésimo bloque de B en array serializado de bloques GrupSize x GrupSize 
     * para distribuir mediante scatter con el comunicador de filas.
     */
    void cargar_bloquesB(float *buffer, float *origen, int size, int cantFilas, int cantCols);
    
    /*
     * Obtiene las coordenadas a las que debe ser enviado el iésimo grupo
     */
    void get_coords_ithgrupo(int diagonal, int *Fila, int *Columna);

    void cargar_bloque(float *Block, float *Buff, int cantFilas, int cantCols);

    void combinar_resultados();
    
    void calcular_producto_externo();
    
    void cargar_matrizC(float *buffer, int fini, int cantFilas, int cini, int cantCols);

    void imprimir_bloque(float *bloque, int size, int fini, int cantFilas, int cini, int cantCols);
  
#ifdef	__cplusplus
}
#endif

#endif	/* _2DDIAGONAL_H */

