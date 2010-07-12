/* 
 * File:   algCannon.c
 * Author: 
 *    - Fernando Mancía
 *    - Cristhian Parra
 *
 */


/**
 * Descripcion del
 * Utiliza una malla con conexiones entre los elementos de cada lado (toro) 
 * para desplazar los elementos de A hacia la izquierda y los de B hacia arriba.

  * El algoritmo sigue los siguientes pasos:
    1. El procesador Pij tiene los elementos aij y bij.
 
    2. La fila i-ésima de A se desplaza i posiciones a la izquierda, 
     y la columna j-ésima de B se desplaza j posiciones hacia arriba, 
     y todo esto teniendo en cuenta que el elemento que sale por un extremo 
     entra por el otro. Con este paso se consigue que el procesador Pij 
     contenga los elementos aij+i y bi+jj, que son necesarios 
     para calcular cij.

    3. Cada procesador multiplica su par de elementos.
     
    4. Cada fila de A se desplaza una posición a la izquierda, y cada 
     columna de B una posición hacia arriba.
    5. Cada procesador multiplica su nuevo par de elementos, y 
      suma el resultado al anterior.
    6. Se repiten los pasos 4 y 5 hasta terminar, es decir n – 1 desplazanientos
**/


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "mmio.h"
#include "commons.h"
#include "salida-metricas.h"
#include <unistd.h>
#include <math.h>



//Funciones Adicionales
void matriz2buff(float *mat, int tamBloq, int process_i, int process_j, float * buffer);
void buff2Bloque(float bloqueX[tamanho_bloque][tamanho_bloque], int tam_bloq, float * buffer);



int main(int argc, char** argv) {
    
    //Solo se recibe N como parametro
    char * ene = argv[1];
    N = atoi(ene);
    result_t r; // estructura de resultados
    
    int p;      //cantidad de procesos. P debe ser cuadrado perfecto.
                // sqrt(p) debe ser divisible N
    int raiz_p;
    
    
    int id_proceso;  // id del proceso actual
    int malla_id;  // id del proceso actual en la topología 2d "malla"
    int dims[2];
    int periods[2];
    int micoord[2];
    
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id_proceso);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm malla;
    MPI_Status estado_msg;

    raiz_p = sqrt(p);
    tamanho_bloque = N / raiz_p;
    dims[0]=dims[1]=raiz_p;
    periods[0]=periods[1]=1;
    MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods, 0, &malla);
    MPI_Comm_rank(malla, &malla_id);
    MPI_Cart_coords(malla, malla_id, 2, micoord);
    
    //Validaciones de los procesos.
    if(raiz_p * raiz_p != p ){
        printf("\nLa cantidad de procesos no es un cuadrado perfecto\n", p, N);
        MPI_Finalize();
        return(EXIT_FAILURE);        
    }    
    
    int comp = N % raiz_p; 
    if(comp!= 0){
        
        printf("\n sqrt(p) no es divisible por N.  P=%d ; N=%d \n", p, N);
        MPI_Finalize();
        return(EXIT_FAILURE);
    }

    
    
    float bloqueA[tamanho_bloque][tamanho_bloque];
    float bloqueB[tamanho_bloque][tamanho_bloque];
    float bloqueC[tamanho_bloque][tamanho_bloque];

    
    float buffA[tamanho_bloque*tamanho_bloque]; //buffer p/ I/O de mat A
    float buffB[tamanho_bloque*tamanho_bloque]; //buffer p/ I/O de mat B
    
    

    if(micoord[0] == 0 && micoord[1] == 0){
        // La alineación inicial de A y B hacemos en el proceso 0
        // leer_matrices.
        A = malloc(sizeof(float) * N * N);
        B = malloc(sizeof(float) * N * N);
        C = malloc(sizeof(float) * N * N);
        leer_matrices();
        cerar_matriz(C);
        
        t1 = MPI_Wtime();
                
        //enviar a cada proceso su parte.
        int i,j;
        for(i=0; i< raiz_p ; i++){
            for(j=0; j< raiz_p ; j++){    
                if (!( i == 0 && j == 0)){
                    
                    if(i+j < raiz_p){
                        matriz2buff(A, tamanho_bloque, i   , i+j , buffA);    
                        matriz2buff(B, tamanho_bloque, i+j , j   , buffB);
                    }else{
                        matriz2buff(A, tamanho_bloque, i            , i+j-raiz_p, buffA);    
                        matriz2buff(B, tamanho_bloque, i+j-raiz_p   , j   , buffB);
                    }
                    
                    

                    int id_destino;
                    int destino[2];
                    destino[0]=i;
                    destino[1]=j;
                    MPI_Cart_rank(malla, destino, &id_destino);
                    
                    MPI_Send(buffA, tamanho_bloque * tamanho_bloque, MPI_FLOAT, id_destino, id_destino, malla );           
                    MPI_Send(buffB, tamanho_bloque * tamanho_bloque, MPI_FLOAT, id_destino, id_destino + (N*N), malla );
                }
            }
        }
        // Para no enviarse a si mismo
        matriz2buff(A, tamanho_bloque, 0, 0, buffA);
        matriz2buff(B, tamanho_bloque, 0, 0, buffB);        
        
        buff2Bloque(bloqueA,tamanho_bloque, buffA);
        buff2Bloque(bloqueB,tamanho_bloque, buffB);

    }else{

        t1 = MPI_Wtime();
        // Todos los procesos reciben su parte de la matriz 
        //printf("\nAntres de recibir en proceso %d. tamanho = %d \n", id_proceso, tamanho_bloque * tamanho_bloque);
        MPI_Recv(buffA, tamanho_bloque * tamanho_bloque, MPI_FLOAT, 0, malla_id , malla, &estado_msg);
        MPI_Recv(buffB, tamanho_bloque * tamanho_bloque, MPI_FLOAT, 0, (N*N)+malla_id , malla, &estado_msg);
        //printf("\nDespues de recibir en proceso %d\n", id_proceso);
    
        //Colocamos en el bloque temporal de A y B
        buff2Bloque(bloqueA,tamanho_bloque, buffA);
        buff2Bloque(bloqueB,tamanho_bloque, buffB);
    }

    
    
    
    //Multiplicamos
    cerar_bloque(bloqueC);
    acumular_multiplicar(bloqueA, bloqueB, bloqueC);
    
   
    int pasos;
    for(pasos=1; pasos< raiz_p; pasos++){
        int origen; int destino;
        // Desplazamiento Para Matriz A
        MPI_Cart_shift(malla, 1, -1, &origen, &destino);
        //printf("\nProc %d Matriz A Del proceso %d al proceso %d\n", malla_id, origen, destino);
        MPI_Sendrecv_replace(buffA, tamanho_bloque * tamanho_bloque, MPI_FLOAT, destino, 0 , origen, 0, malla, &estado_msg );           


        
        // Desplazamiento Para Matriz B
        MPI_Cart_shift(malla, 0, -1, &origen, &destino);
        MPI_Sendrecv_replace(buffB, tamanho_bloque * tamanho_bloque, MPI_FLOAT, destino, 0 , origen, 0, malla, &estado_msg );               

        buff2Bloque(bloqueA, tamanho_bloque, buffA);
        buff2Bloque(bloqueB, tamanho_bloque, buffB);


        acumular_multiplicar(bloqueA, bloqueB, bloqueC);
    }
    t2 = MPI_Wtime();
    t3 = t2 -t1;

    MPI_Barrier(MPI_COMM_WORLD);
    
    
    //Se colocan los tiempos en un archivo
    double buffTiempo[1];
    
    if(id_proceso == 0){
        int proc;
        char * FCSV = malloc(sizeof (char) * PATH_BUFF_SIZE);
        sprintf(FCSV,"%s-%s.csv","salida", "CANON");
        FILE *fcsv = open_csvFile_Append(FCSV);
        print_header(fcsv, FCSV);
        
        sprintf(r.testId, "%s", "CANNON");
        r.T = t3;
        r.q = p;
        r.size = N;
        r.prank = 0;
        print_test_result(fcsv, &r);
        
        for(proc = 1; proc < p; proc++){    
            MPI_Status stat;
            MPI_Recv(buffTiempo , 1, MPI_DOUBLE, proc, proc + (N*N*2) , MPI_COMM_WORLD, &stat);
            printf("Se envia %f seg ",buffTiempo[0] );
            
            double tmp = buffTiempo[0] ;
            r.T = tmp;
            r.q = p;
            r.size = N;
            r.prank = proc;
            
            print_test_result(fcsv, &r);
            
        }
        close_csvFile(fcsv);
    }else{
        buffTiempo[0]=t3;
        
        MPI_Send(buffTiempo, 1, MPI_DOUBLE, 0, id_proceso + (N*N*2), MPI_COMM_WORLD);
    }
    
    
    //Fin del Algoritmo Paralelo
    printf("\n#Fin proceso %d en %f segundos#\n", id_proceso, t3);
    MPI_Finalize();
    return (EXIT_SUCCESS);
}



void matriz2buff(float *mat, int tamBloq, int process_i, int process_j, float * buffer){
    
    // suponiendo que los procesos que tienen los pedazon 
    // van de 1 a p
    int fil = process_i * tamBloq;
    int col = process_j * tamBloq;
    int i,j, c=0;
    for(i = fil; i < (tamBloq +fil); i++){
        
        for(j = col; j < (tamBloq +col); j++){
            buffer[c]=mat[(i*N)+j];
            c++;
        }
    }
    
}

/**
 * Pasa los valores desde el buffer al bloque bloqueX.
 * Se supone que el buffer tiene el bloque por filas
 *
 **/
void buff2Bloque(float bloqueX[tamanho_bloque][tamanho_bloque], int tam_bloq, float * buffer){
    
    int i,j, c=0;
    for(i = 0; i < tam_bloq; i++){
        for(j = 0; j < tam_bloq; j++){
            bloqueX[i][j]=buffer[c];
            c++;
        }
    }
    
}

