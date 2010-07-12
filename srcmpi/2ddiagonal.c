/*****************************************************************************
 *
 * Implementación del Algoritmo Paralelo para multiplicación de matrices 
 * 2D Diagonal, propuesto en el Paper "Communication Efficient Matrix 
 * Multiplication on Hypercubes" de Himanshu Gupta y P. Sadayappan. 
 *
 *
 * AUTORES: 
 * - Cristhian Daniel Parra Trepowski
 * - Fernando Manuel Mancía Zelaya. 
 *   
 * Disponible en: http://mmulfpuna.googlecode.com/
 *
 * DESCRIPCIÓN DEL ALGORITMO IMPLEMENTADO. 
 *
 * El algoritmo realiza una distribución inicial en los procesos de la diagonal 
 * de la malla, y luego realiza dos grupos de comunicaciones más para obtener
 * el resultado final de nuevo en los procesos diagonales. 
 *
 * LIMITACIONES DE LA IMPLEMENTACIÓN: 
 * - N y P deben ser cuadrados perfecto.
 *
 * OBSERVACIÓN: 
 *
 * La implementación, utiliza la librería MPI "OpenMPI" y fue desarrollado 
 * en una plataforma Linux. La librería OpenMPI no provee actualmente el 
 * soporte adecuado para su desarrollo en Windows, por lo cual no se garantiza
 * el funcionamiento en dicha plataforma. 
 *****************************************************************************/

#include <stdio.h>


#include "commons.h"
#include "2ddiagonal.h"

void distribucion_inicial() {
    int k, fini, cini, destrank;
    int dest[2];

    /*
     * A cada Proceso Pk,k se le envía el iésimo grupo de columnas de A y 
     * el iésimo grupo de filas de B
     *
     */
    if (mycoords[0] == 0 && mycoords[1] == 0) {
        for (k = 1; k < Grupos; k++) {
            
            get_coords_ithgrupo(k, &fini, &cini);
           
            // Cargo en IthGrupoA GrupSize columnas de MatrizA desde la columna 
            // cini
            cargar_buffer_grupo(IthGrupoA, MatrizA, N, 0, N, cini, GrupSize);
            
            if (Debug > 0) {
                printf("(%d,%d) --> (%d,%d) | IthGrupoA: (%d,%d)\n",mycoords[0],mycoords[1],k,k,0,cini); 
            }
            
            // Cargo en IthGrupoB GrupSize filas de MatrizA desde la fila 
            // fini         
            cargar_buffer_grupo(IthGrupoB, MatrizB, N, fini,GrupSize, 0, N);

            if (Debug > 0) {
                printf("(%d,%d) --> (%d,%d) | IthGrupoB: (%d,%d)\n",mycoords[0],mycoords[1],k,k,fini,0);
            }
            
            dest[0] = dest[1] = k;

            MPI_Cart_rank(Comm2d, dest, &destrank);

            MPI_Send(IthGrupoA, GrupSize*N, MPI_FLOAT, destrank, k*10, Comm2d);
            MPI_Send(IthGrupoB, GrupSize*N, MPI_FLOAT, destrank, k*10+1, Comm2d);

        }
    } else {
        if (mycoords[0] == mycoords[1]) {
            MPI_Status recstatus;
            MPI_Recv(IthGrupoA, GrupSize*N, MPI_FLOAT, 0, mycoords[0]*10, Comm2d, &recstatus);
            MPI_Recv(IthGrupoB, GrupSize*N, MPI_FLOAT, 0, mycoords[0]*10+1, Comm2d, &recstatus);
            
            if (Debug > 0) {
                printf("(%d,%d) <-- (%d,%d) | Grupos A recibido!!\n", mycoords[0], mycoords[1], 0, 0, fini, 0); 
                printf("(%d,%d) <-- (%d,%d) | Grupos B recibido!!\n", mycoords[0], mycoords[1], 0, 0, fini, 0); 
            }
        }
    }

    /*
     * Cargamos los Grupos de filas y columnas para el Proces 0,0
     *
     */
    
    if (mycoords[0] == 0 && mycoords[1] == 0) {
        fini = cini = 0; 
        cargar_buffer_grupo(IthGrupoA, MatrizA, N, 0, N, cini, GrupSize);
        cargar_buffer_grupo(IthGrupoB, MatrizB, N, fini, GrupSize, 0, N);

        if (Debug > 0) {
            printf("(%d,%d) --> (%d,%d) | IthGrupoA: (%d,%d)\n", mycoords[0], mycoords[1], 0, 0, 0, cini);
            printf("(%d,%d) --> (%d,%d) | IthGrupoB: (%d,%d)\n", mycoords[0], mycoords[1], 0, 0, fini, 0);
        }
    }
    

}

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
void ejecutar2dd() {

    int remains[2];
    int i,j,k;
    
    /************************************************************************
     * Paso #0: Creamos comunicadores para realizar las comunicaciones.     *
     ************************************************************************/
    remains[0] = 0; 
    remains[1] = 1; 
    MPI_Cart_sub(Comm2d,remains,&commFilas);
    
    remains[0] = 1; 
    remains[1] = 0; 
    MPI_Cart_sub(Comm2d,remains,&commColumnas);
    
    /************************************************************************
     * Paso #1: Broadcast de A*,j a P*,j desde Pj,j                         *
     ************************************************************************/
    int rrank = mycoords[0];
    MPI_Bcast(IthGrupoA,GrupSize*N,MPI_FLOAT,rrank,commFilas);
    
    
    /************************************************************************
     * Paso #2: Send Bik a Pkj desde Pii                                    *
     ************************************************************************/
    BloquesB = malloc(sizeof(float)*GrupSize*N);
    
    for (k = 0; k < Grupos; k++) {
                
        if (mycoords[0] == mycoords[1] && mycoords[0] == k) {            
            cargar_bloquesB(BloquesB,IthGrupoB,N,GrupSize,GrupSize);  
            
            if (Debug > 0) {
                printf("(%d,%d) --> Buffer para scatter en fila %d Listo!!! ",mycoords[0],mycoords[1],k);
            }     
        }        
        
        if ( mycoords[0] == k) {
            MPI_Scatter(BloquesB, GrupSize*GrupSize, MPI_FLOAT, BloqueB,GrupSize*GrupSize, MPI_FLOAT, k , commFilas);            
        }
    }
        
    printf("(%d,%d) --> Scatter Realizado. BloqueB recibido!\n",mycoords[0],mycoords[1]);
    fflush(stdout);
    
    calcular_producto_externo();
        
    printf("(%d,%d) --> Producto Externo Calculado!\n",mycoords[0],mycoords[1]);
    fflush(stdout);
    
    int group; 
    
    for(group=0; group<Grupos; group++) {
        if (mycoords[1] == group) {
            MPI_Reduce(IProduct,IProductReduced,GrupSize*N,MPI_FLOAT,MPI_SUM,group,commColumnas);        
        }
    }
    
    printf("(%d,%d) --> Resultado Reducido!\n",mycoords[0],mycoords[1]);
    free(BloquesB);
}

void combinar_resultados() {

    int group = mycoords[0];
    
    if (mycoords[0] == mycoords[1] && mycoords[0] != 0) {
            int source[2];
            int sourcerank;
            source[0] = source[1] = group;
            MPI_Cart_rank(Comm2d, source, &sourcerank);

            int root[2];
            int rootrank;
            root[0] = root[1] = 0;
            MPI_Cart_rank(Comm2d, root, &rootrank);

            MPI_Send(IProductReduced, GrupSize*N, MPI_FLOAT, rootrank, mycoords[0]*TAGMULT, Comm2d);

    } else if (mycoords[0] == mycoords[1] && mycoords[0] == 0) {
        
        cargar_matrizC(IProductReduced,0,N,0,GrupSize);
        
        for (group = 1; group < Grupos; group++) {
            int source[2];
            int sourcerank;
            source[0] = source[1] = group;
            MPI_Cart_rank(Comm2d, source, &sourcerank);

            MPI_Status statusfinal;
            MPI_Recv(IProduct,GrupSize*N,MPI_FLOAT,sourcerank,group*TAGMULT,Comm2d,&statusfinal);
            
            int fini,cini;
            
            get_coords_ithgrupo(group,&fini,&cini);
            
            cargar_matrizC(IProduct,0,N,cini,GrupSize);
        }
    }
}

void cargar_buffer_grupo(float *buffer, float *origen, int size, int fini, int cantFilas, int cini, int cantCols) {
    int i, j, buffcount;

    if (Debug > 0) {
        printf("\n(%d,%d) Cargando Buffer (%d,%d). CantFilas %d, CantCols %d.\n", 
                mycoords[0],mycoords[1], fini, cini,cantFilas,cantCols);
        fflush(stdout);
    }

    buffcount = 0;
    for (i = fini; i<(fini+cantFilas); i++) {
        for (j = cini; j<(cini+cantCols); j++) {
            buffer[buffcount] = origen[j+(size*i)];
            buffcount++;
        }
    }

    if (Debug > 0) {
        printf("(%d,%d)BloquesB cargado Exitosamente!!!\n", mycoords[0], mycoords[1]);
        fflush(stdout);
    }
}

/*
 * Carga el Iésimo bloque de B en array serializado de bloques GrupSize x GrupSize 
 * para distribuir mediante scatter con el comunicador de filas.
 */
void cargar_bloquesB(float *buffer, float *origen, int size, int cantFilas, int cantCols) {
    int i, j, k, buffcount;

    if (Debug > 0) {
        printf("(%d,%d) --> Concatenando los BloqueB.\n", mycoords[0], mycoords[1]);
        fflush(stdout);
    }

    buffcount = 0;
    for (k = 0; k < N; k = k+GrupSize) {
        for (i = 0; i < cantFilas; i++) {
            for (j = k; j < (k+cantCols); j++) {
                int pos = j+(size*i);           
                buffer[buffcount] = origen[pos];
                buffcount++;
            }
        }
    }

    if (Debug > 0) {
        printf("(%d,%d) --> Buffer CARGADO!!!\n", mycoords[0],mycoords[0]);
        fflush(stdout);
    }
}

void get_coords_ithgrupo(int grupo, int *Fila, int *Columna) {
    *Fila = *Columna = GrupSize*grupo;
}

void calcular_producto_externo() {
    int i, j, k = 0;
    
    for (i = 0; i < N; i++) {
        for (j = 0; j < GrupSize; j++) {            
            
            IProduct[(i*GrupSize)+j] = 0;            
            for (k = 0; k < GrupSize; k++) {                                                
                float aux = IthGrupoA[(i*GrupSize)+k]*BloqueB[(k*GrupSize)+j];
                int pos = (i*GrupSize)+j;                
                IProduct[pos] += aux;
            }
        }
    }    
}

void cargar_matrizC(float *buffer, int fini, int cantFilas, int cini, int cantCols) {
    int i, j, k;
    k = 0;
    for (i = fini; i < (fini + cantFilas); i++) {
        for (j = cini; j < (cini + cantCols); j++) {
            int pos = j+(N*i);
            MatrizC[pos]=buffer[k];
            k++;
        }
    }
}

void imprimir_bloque(float *bloque, int size, int fini, int cantFilas, int cini, int cantCols) {
    int i, j, k;
    k = 0;
    printf("\n");
    for (i = fini; i < (fini + cantFilas); i++) {
        for (j = cini; j < (cini + cantCols); j++) {
            int pos = j+(size*i);
            printf(" %4f(%d,%d) * ",bloque[pos],mycoords[0],mycoords[1]);
        }
        printf("\n");
    }
}
