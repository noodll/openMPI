/* 
 * File:   algDNS.c
 * Author: fmancia
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <malloc.h>
#include "mmio.h"
#include "commons.h"
#include "salida-metricas.h"
#include <unistd.h>
#include <math.h>

/*

 *Comunicacion:
  1 Mover las columnas de A y las filas de B a sus respectivos planos.
  2 Ejecutar un broadcast a lo largo del eje j para A y a lo largo del eje i para B.
  3 Acumular los nodos simples a lo largo del eje k.
	
 Todas estas operaciones se ejecutan en subcubos del hipercubo de n procesadores.
 */

//Funciones Adicionales
void matriz2buff(float *mat, int tamBloq, int process_i, int process_j, float * buffer);
void buff2Bloque(float bloqueX[tamanho_bloque][tamanho_bloque], int tam_bloq, float * buffer);
void bloque2Buff(float bloqueX[tamanho_bloque][tamanho_bloque], int tam_bloq, float * buffer);



int main(int argc, char** argv) {

    char * ene = argv[1];
    N = atoi(ene);
    result_t r; // estructura de resultados
    
    int p;      //cantidad de procesos. P debe ser cuadrado perfecto.
                // N debe ser divisible raiz_cubica(p) 
    int raiz_cubica_p;
    int id_proceso;  // id del proceso actual 
    
    int dims[3], micoord[3], periods[3];
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id_proceso);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    MPI_Comm comm_3d;
    int comm_3d_id;
    
    MPI_Comm comm_plano_2d;
    int comm_plano_2d_id;
    
    MPI_Status estado_msg;

    raiz_cubica_p = cbrt(p);
    tamanho_bloque = N / raiz_cubica_p;
    dims[0]=dims[1]=dims[2]=raiz_cubica_p;
    periods[0]=periods[1]=periods[2]=0;

    //Verificación del algoritmo
    if( (raiz_cubica_p * raiz_cubica_p * raiz_cubica_p) != p){
        
        printf("\n P no es cubo perfecto. No se puede ejecutar el algoritmo asi. P=%d ; N=%d \n", p, N);
        MPI_Finalize();
        return(EXIT_FAILURE);
    }
    int comp = N % raiz_cubica_p;
    if( comp!= 0){
        
        printf("\n La raiz_cubica(p) no es divisible por N. No se puede ejecutar "
                " el algoritmo asi. P=%d ; N=%d \n", p, N);
        MPI_Finalize();
        return(EXIT_FAILURE);
    }

    
    //Creamos la topología en 3D
    MPI_Cart_create(MPI_COMM_WORLD,3,dims,periods, 0, &comm_3d);
    MPI_Comm_rank(comm_3d, &comm_3d_id);
    MPI_Cart_coords(comm_3d, comm_3d_id, 3, micoord);
    
    
    float bloqueA[tamanho_bloque][tamanho_bloque];
    float bloqueB[tamanho_bloque][tamanho_bloque];
    float bloqueC[tamanho_bloque][tamanho_bloque];

    float buffA[tamanho_bloque*tamanho_bloque]; //buffer p/ I/O de mat A
    float buffB[tamanho_bloque*tamanho_bloque]; //buffer p/ I/O de mat B
    float buffC[tamanho_bloque*tamanho_bloque]; 
    
    
   
    // Creamos un comunicador comm_plano_2d
    /**int remain[3];
    remain[0] = 1;
    remain[1] = 1;
    remain[2] = 0;
    MPI_Cart_sub(comm_3d, remain, &comm_plano_2d );
    MPI_Comm_rank(comm_plano_2d, &comm_plano_2d_id);**/

        
    
    //Si soy el primer proceso
    if(micoord[0]== 0 && micoord[1]==0 && micoord[2] == 0){
        // leer_matrices.
        A = malloc(sizeof(float) * N * N);
        B = malloc(sizeof(float) * N * N);
        C = malloc(sizeof(float) * N * N);
        leer_matrices();
        cerar_matriz(C);
        
        // Cargamos las matrices A y B en los buffers 
        // para enviar deps vía Scater en el plano 0
        int i,j;
        for(i=0; i< raiz_cubica_p; i++){
            for(j=0; j< raiz_cubica_p; j++){
                matriz2buff(A,tamanho_bloque, i, j, buffA);
                matriz2buff(B,tamanho_bloque, i, j, buffB);
                if(i==0 && j==0){
                    buff2Bloque(bloqueA, tamanho_bloque, buffA);
                    buff2Bloque(bloqueB, tamanho_bloque, buffB);
                    
                }else{
                    int dest[3] = {i,j,0};
                    int id_dest;
                    MPI_Cart_rank(comm_3d, dest, &id_dest);
                    
                    MPI_Send(buffA, tamanho_bloque * tamanho_bloque, MPI_FLOAT, id_dest, id_dest, comm_3d);
                    MPI_Send(buffB, tamanho_bloque * tamanho_bloque, MPI_FLOAT, id_dest, id_dest + (N*N), comm_3d);
                }
            }
        }
    }
    else if(micoord[2] == 0)
    {
        // si no soy el proceso (0,0,0) 
        // pero estoy en el plano 0 debo recibir
        int sour[3] = {0,0,0};
        int id_sour;
        MPI_Cart_rank(comm_3d, sour, &id_sour);
        
        MPI_Recv(buffA, tamanho_bloque * tamanho_bloque, MPI_FLOAT, id_sour, comm_3d_id, comm_3d, &estado_msg);
        MPI_Recv(buffB, tamanho_bloque * tamanho_bloque, MPI_FLOAT, id_sour, comm_3d_id + (N*N), comm_3d, &estado_msg);
        
        buff2Bloque(bloqueA, tamanho_bloque, buffA);
        buff2Bloque(bloqueB, tamanho_bloque, buffB);
        
    }
    
    
    //Una vez que cada procesador del plano 0 tenga sus elementos empezamos a medir el tiempo
    t1 = MPI_Wtime();
  
    //En el plano 0, 
    if(micoord[2] == 0){
    
        //Ahora debemos enviar buffA a los planos j = k en el sgte ciclo
        if(micoord[1] != 0){
            int destino[3];
            destino[0]=micoord[0];
            destino[1]=micoord[1];       
            destino[2]=micoord[1];
            int id_dest;
            MPI_Cart_rank(comm_3d, destino, &id_dest);
            
            MPI_Send(buffA, tamanho_bloque * tamanho_bloque, MPI_FLOAT,
                     id_dest, id_dest ,comm_3d);
            
        }
        
        if(micoord[0] != 0){
            int destino[3];
            destino[0]=micoord[0];
            destino[1]=micoord[1];       
            destino[2]=micoord[0];
            int id_dest;
            MPI_Cart_rank(comm_3d, destino, &id_dest);
            
            MPI_Send(buffB, tamanho_bloque * tamanho_bloque, MPI_FLOAT,
                     id_dest, id_dest + (N*N),comm_3d);
            
        }
                

    }else{
        //En los demas planos recibimos 
        if(micoord[1] == micoord[2]){
            int fuente[3];
            fuente[0]=micoord[0]; // la misma fila;
            fuente[1]=micoord[1]; // la misma columna;   
            fuente[2]=0;          //plano cero
            int id_fuente;
            MPI_Cart_rank(comm_3d, fuente, &id_fuente);
        
            MPI_Status stA;
            MPI_Recv(buffA, tamanho_bloque * tamanho_bloque, MPI_FLOAT,
                     id_fuente, comm_3d_id  , comm_3d, &stA);        
    
            buff2Bloque(bloqueA, tamanho_bloque, buffA);
            
        }
        
        if(micoord[0] == micoord[2]){
            int fuente[3];
            fuente[0]=micoord[0]; // la misma fila;
            fuente[1]=micoord[1]; // la misma columna;   
            fuente[2]=0;          //plano cero
            int id_fuente;
            MPI_Cart_rank(comm_3d, fuente, &id_fuente);
        
            MPI_Status stB;
            MPI_Recv(buffB, tamanho_bloque * tamanho_bloque, MPI_FLOAT,
                     id_fuente, comm_3d_id + (N*N) , comm_3d, &stB);        
            
            buff2Bloque(bloqueB, tamanho_bloque, buffB);
        }
        
    }

    

    // Ahora C/ proceso que está en el plano j y tiene la columna j 
    // debe replicar a los demás en su plano (Bcast creo q se puede usar)
    int remainA[3]={0,1,0};
    MPI_Comm comm_1d_A;
    MPI_Cart_sub(comm_3d, remainA, &comm_1d_A );
    MPI_Bcast(buffA, tamanho_bloque * tamanho_bloque, MPI_FLOAT, micoord[2], comm_1d_A);
    
    int remainB[3]={1,0,0};
    MPI_Comm comm_1d_B;
    MPI_Cart_sub(comm_3d, remainB, &comm_1d_B );
    MPI_Bcast(buffB, tamanho_bloque * tamanho_bloque, MPI_FLOAT, micoord[2], comm_1d_B);
    
    // Hasta aquí cada procesador ya tiene un bloque de la matriz

    
    
    // Se multiplican los bloques de A y B q tiene el proceso
    buff2Bloque(bloqueA,tamanho_bloque, buffA );
    buff2Bloque(bloqueB,tamanho_bloque, buffB );
    
 
    cerar_bloque(bloqueC);
    acumular_multiplicar(bloqueA, bloqueB, bloqueC);
    bloque2Buff(bloqueC, tamanho_bloque, buffC);
    

    //Creamos un comunicador vertical 
    MPI_Comm comm_vertical;
    int comm_vertical_id;
    int remain[3]={0,0,1};
    MPI_Cart_sub(comm_3d, remain, &comm_vertical);
    MPI_Comm_rank(comm_vertical, &comm_vertical_id);

    
    //Para saber quien es el que recibe
    /**int elroot[3];
    elroot[0]=micoord[0];
    elroot[1]=micoord[1];
    elroot[2]=0;
    int elroot_id;
    MPI_Cart_rank(comm_3d, elroot, &elroot_id);**/
    
    
    float bloqueFinal[tamanho_bloque][tamanho_bloque];
    float bufferFinal[tamanho_bloque * tamanho_bloque];
    
    MPI_Reduce( buffC , bufferFinal, tamanho_bloque * tamanho_bloque, MPI_FLOAT,
                    MPI_SUM, 0, comm_vertical);
    
    buff2Bloque(bloqueFinal, tamanho_bloque, bufferFinal);
    
    //Tomamos el tiempo y guardamos en un archivo
    t2 = MPI_Wtime();
    t3 = t2 -t1;

    MPI_Barrier(MPI_COMM_WORLD);
    
    
    //Se colocan los tiempos en un archivo
    double buffTiempo[1];
    
    if(id_proceso == 0){
        int proc;
        char * FCSV = malloc(sizeof (char) * PATH_BUFF_SIZE);
        sprintf(FCSV,"%s-%s.csv","salida", "DNS");
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


void buff2Bloque(float bloqueX[tamanho_bloque][tamanho_bloque], int tam_bloq, float * buffer){
    
    int i,j, c=0;
    for(i = 0; i < tam_bloq; i++){
        for(j = 0; j < tam_bloq; j++){
            bloqueX[i][j]=buffer[c];
            c++;
        }
    }
    
}
void bloque2Buff(float bloqueX[tamanho_bloque][tamanho_bloque], int tam_bloq, float * buffer){
    
    int i,j, c=0;
    for(i = 0; i < tam_bloq; i++){
        for(j = 0; j < tam_bloq; j++){
            buffer[c]=bloqueX[i][j];
            c++;
        }
    }
    
}
