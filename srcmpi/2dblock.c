/*****************************************************************************
 *
 * Implementación del Algoritmo Paralelo para multiplicación de matrices 
 * 2D con distribución por Bloques Cíclico. 
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
 * Basado en la distribución por bloques cíclica presentada en el capítulo 3 
 * del Libro "Introduction to Parallel Computing, Second Edition. Ananth Grama, 
 * Anshul Gupta, George Kapyris, Vipin Kumar".
 *
 * El algoritmo particiona la matriz de salida de la multiplicación y realiza
 * una distribución cíclica de las tareas y sus datos entre los P procesos. 
 *
 * El Proceso 0 es el que realiza el recorrido en bloque (según la fuente del 
 * libro) y va realizando la asignación en un esquema maestro-esclavo, donde P0
 * es el maestro y los demás procesos son los esclavos. 
 *
 *
 * LIMITACIONES DE LA IMPLEMENTACIÓN: 
 * - N y P deben ser cuadrados perfecto. 
 * - El Proceso 0 tiene una excesiva sobrecarga porque también realiza 
 *   computación. Se podría adaptar la coputación para que exista un procesos 
 *   exclusivamente dedicado a la tarea de ser "Maestro". 
 *
 * 
 * OBSERVACIÓN: 
 *
 * La implementación, utiliza la librería MPI "OpenMPI" y fue desarrollado 
 * en una plataforma Linux. La librería OpenMPI no provee actualmente el 
 * soporte adecuado para su desarrollo en Windows, por lo cual no se garantiza
 * el funcionamiento en dicha plataforma. 
 *****************************************************************************/

#include <mpi.h>
#include <mpi.h>
#include "2dblock.h"

void ejecutar_2d(int prank) {

    int i, j; // índices de los Bloques de tareas. 
    int x, y; // índices de los bloques de procesos. 

    // Salto de Bloque para nuevo ciclo de distribución cíclica en los P procesos
    //int SaltoBloqueP = floor((double) N / ProcBlockSize);
    int SaltoBloqueP = ProcBlockSize * SubBlockSize;
    

    // Salto de Bloque para asignar una nueva Tarea a un Proceso
    int SaltoBloqueT = floor((double) N / SubBlockSize); // el fondo es lo mismo que hacer SubBlockSize :P

    int BuffSize = N * SubBlockSize * 2; // para contener a N*ProcBlockSize filas de A
    // y N*ProcBlockSize columnas de A

    int BuffCut = BuffSize / 2; // variable de corte del buffer para extraer
    // filas y columnas

    int BuffInnerCut = N;

    float *Buffer = malloc(sizeof (float) * BuffSize);
    
    float *BufferLocal = malloc(sizeof (float) * BuffSize);

    int SalidaSize = SubBlockSize * SubBlockSize;

    if (Debug > 0 && prank == 0) {
        printf("(prank:%d)Valores del algoritmo: \n", prank);
        printf("(prank:%d)=========================\n");
        printf("(prank:%d) --> N: %d\n", prank, N);
        printf("(prank:%d) --> P: %d\n", prank, P);
        printf("(prank:%d) --> SubBlockSize: %d\n", prank, SubBlockSize);
        printf("(prank:%d) --> ProcBlockSize: %d\n", prank, ProcBlockSize);
        printf("(prank:%d) --> SaltoBloqueP: %d\n", prank, SaltoBloqueP);
        printf("(prank:%d) --> SaltoBloqueT: %d\n", prank, SaltoBloqueT);
        printf("(prank:%d) --> BuffSize: %d\n", prank, BuffSize);
        printf("(prank:%d) --> BuffCut: %d\n", prank, BuffCut);
        printf("(prank:%d) --> BuffInnercut: %d\n", prank, BuffInnerCut);
        fflush(stdout);
    }

    float *Salida = malloc(sizeof (float) * SalidaSize);

    int *filas_part = malloc(sizeof (int) * N); // memoria de filas de inicio asignadas una tarea
    int *cols_part = malloc(sizeof (int) * N); // memoria de columnas de inicio asignadas a una tarea

    /* 
     * El proceso de distribución y particionamiento está basado en lo que se 
     * explica en la Sección de 3.4 del Libro "Introduction to Parallel 
     * Computing, Second Edition. Ananth Grama, Anshul Gupta, George Kapyris, 
     * Vipin Kumar
     *
     * Para ver un ejemplo, observar la figura 3.30
     *
     * El algoritmo es simple, P0, recorre la matriz C en bloques de 
     * SubBlockSize x SubBlockSize, para ir asignando las tareas en 
     * bloques de sqrt(P) x sqrt(p), de manera cíclica a cada proceso. 
     */

    int k = 0; // indice de Tarea a asignar
    int done = 0;
    int TagDone = 10000;
    int BufferDone[1];
    int Tags[P];
    int z;
    for (z = 0; z < P; z++) {
        Tags[z] = 0;
    }


    if (prank == 0) { // Si soy el proceso 0, envío datos y proceso mi parte
        int Pdestino = 0; /* Indica a quien le envío la tarea actual. 
                                               * Al final del ciclo de asignación, indica 
                                               * cuál fue el último proceso que recibió.
                                               */
        int PdestinoFinalAnterior = P - 1; /* Indica en que proceso se inició la asignación en el 
                                               * ciclo anterior. 
                                               */

        for (i = 0; i < N; i += SaltoBloqueP) { // Columnas de C
            for (j = 0; j < N; j += SaltoBloqueP) { // Filas de C
                /**
                 * Los ciclos externos nos permiten saltar de un bloque de 
                 * sqrt(P) x sqrt(P) a otro para un nuevo ciclo de asignación 
                 * de cíclica de tareas. 
                 */

                int Tag;
                for (x = i; x < (i + SaltoBloqueP); x += SaltoBloqueT) { // Cols. de C
                    for (y = j; y < (j + SaltoBloqueP); y += SaltoBloqueT) { // Filas de C
                        /**
                         * Los ciclos internos nos permiten saltar de un bloque de 
                         * sqrt(N) x sqrt(N) a otro para una nueva asignación
                         * de tarea. 
                         */

                        if (y < N && x < N) {

                            Pdestino = k % P;

                            int fini = y;
                            int cini = x;

                            filas_part[k] = fini;
                            cols_part[k] = cini;

                            int cantFilas = SubBlockSize;
                            int cantCols = SubBlockSize;


                            if (Pdestino != 0) {

                                cargar_buffer_envio(Buffer, fini, cantFilas, cini,
                                        cantCols, BuffCut);

                                Tag = Tags[Pdestino] + Pdestino*N;
                                MPI_Send(Buffer, BuffSize, MPI_FLOAT, Pdestino, Tag, MPI_COMM_WORLD);
                                Tags[Pdestino]++;

                                if (Debug > 0) {
                                    printf("(prank:%d->%d))SubBloque (%d,%d) enviado a Proceso %d. Tag %d. CantFilas = %d, Cantcols = %d\n",
                                            prank, Pdestino, fini, cini, Pdestino, Tag, cantFilas, cantCols);
                                    fflush(stdout);
                                }
                                
                            } else {

                                cargar_buffer_envio(BufferLocal, fini, cantFilas, cini,
                                        cantCols, BuffCut);
                                
                                if (Debug > 0) {
                                    printf("(prank:%d->%d))SubBloque (%d,%d) enviado a Proceso %d. Tag %d. CantFilas = %d, Cantcols = %d\n",
                                            prank, Pdestino, fini, cini, Pdestino, 0, cantFilas, cantCols);
                                    fflush(stdout);
                                }
                                
                            }                            
                            
                            k++;
                        }
                    }
                }


                if ((i + SaltoBloqueP) >= N &&
                        (j + SaltoBloqueP) >= N) { // Seguirán habiendo bloques para enviar
                    done = 1;
                } else {
                    done = 0;
                }

                // estas dos variables se usan para los casos en que N%P != 0, 
                // pudiendo haber ciclos en los que algunos procesos no tienen 
                // tarea asignada. 
                int pfrom = (PdestinoFinalAnterior + 1) % P;
                int pto = Pdestino;
                
                PdestinoFinalAnterior = pto;

                if (pfrom == 0) {
                    printf("(prank:%d)Procesando parte local\n", prank);
                    fflush(stdout);
                    multiplicar_bloque(BufferLocal, Salida, BuffCut, prank,k-P);

                    if (Debug > 0) {
                        printf("(prank:%d<-%d))SubBloque resultado (%d,%d) recibido del proceso %d. Tag: %d.\n",
                                prank, prank, filas_part[k - P], cols_part[k - P], prank, 0);
                        fflush(stdout);
                    }
                    
                    cargar_resultado_C(Salida, filas_part[k - P], SubBlockSize, cols_part[k - P], SubBlockSize);

                }

                printf("(prank:%d)Empezando a recibir resultados...\n", prank);
                fflush(stdout);
                int m;
                for (m = pfrom; m < (pto + 1); m++) {
                    if (m != 0) {

                        MPI_Status status2d;
                        Tag = Tags[m] + m*N;
                        MPI_Recv(Salida, SubBlockSize * SubBlockSize, MPI_FLOAT, m, Tag, MPI_COMM_WORLD, &status2d);
                        Tags[m]++;
                        
                        if (Debug > 0) {
                            printf("(prank:%d<-%d))SubBloque resultado (%d,%d) recibido del proceso %d. Tag: %d.\n",
                                    prank, m, filas_part[k - (P - m)], cols_part[k - (P - m)], m, Tag);
                            fflush(stdout);
                        }
                        
                        
                        cargar_resultado_C(Salida, filas_part[k - (P - m)], SubBlockSize, cols_part[k - (P - m)], SubBlockSize);

                        // Enviar confirmación de que el proceso continuará o No. 

                        BufferDone[0] = done;
                        int TagD = Tags[m] + m*N;
                        Tags[m]++;
                        MPI_Send(BufferDone, 1, MPI_INT, m, TagD, MPI_COMM_WORLD);
                        
                        if (Debug > 0) {
                            printf("(prank:%d->%d))Continuar = %d, enviado al proceso %d. Tag: %d\n",
                                    prank, m, done, m, TagD);
                            fflush(stdout);
                        }
                    }
                }
                printf("\n(prank:%d)Todo recibido..\n", prank);
            }
        }
    } else {

        int tagsr = 0;
        int tagr;
        while (done == 0) {
            tagr = tagsr + N*prank;
            tagsr++;

            MPI_Status status2d;
            MPI_Recv(Buffer, BuffSize, MPI_FLOAT, 0, tagr, MPI_COMM_WORLD, &status2d);

            if (Debug > 0) {
                printf("(prank:%d<-0)Bloques recibidos. Tag(%d). Multiplicando Matrices...\n",
                        prank, tagr);
                fflush(stdout);
            }

            multiplicar_bloque(Buffer, Salida, BuffCut, prank,N);

            tagr = tagsr + N*prank;
            tagsr++;

            if (Debug > 0) {
                printf("(prank:%d->0)Multiplicación terminada. Enviando Resultado con Tag: %d.\n",
                        prank, tagr);
                fflush(stdout);
            }
            
            MPI_Send(Salida, SubBlockSize * SubBlockSize, MPI_FLOAT, 0, tagr, MPI_COMM_WORLD);

            tagr = tagsr + N*prank;
            tagsr++;

            MPI_Recv(BufferDone, 1, MPI_INT, 0, tagr, MPI_COMM_WORLD, &status2d);

            if (Debug > 0) {
                printf("(prank:%d<-0) Confirmación de continuación recibida con Tag: %d.\n",
                        prank, tagr);
            }
            
            done = BufferDone[0];
        }

        if (Debug > 0) {
            printf("(prank:%d) Terminado. Último Tag Recibido: %d.\n",
                    prank, tagr);
            fflush(stdout);
        }
    }

    free(Buffer);
    free(Salida);
    free(cols_part);
    free(filas_part);
}

void cargar_buffer_envio(float *buffer, int fini, int cantFilas, int cini, int cantCols, int BuffCut) {
    int i, j, buffcount, ciniaux, prank;

    MPI_Comm_rank(MPI_COMM_WORLD, &prank);
    
    if (Debug > 0) {
        printf("\n(prank:%d)Cargando Buffer (%d,%d).\n", prank, fini, cini);
        fflush(stdout);
    }
    
    buffcount = 0;
    // Cargar filas
    for (i = 0; i < cantFilas; i++) {
        for (j = 0; j < N; j++) {
            buffer[j + (N * i)] = MatrizA[j + (N * fini)];
            buffcount++;
        }
        fini++;
    }

    // Cargar columnas, desde la mitad del buffer
    for (i = 0; i < cantCols; i++) {
        ciniaux = cini;
        for (j = 0; j < N; j++) {
            buffer[j + (N * i) + BuffCut] = MatrizB[i + (N * ciniaux)];
            buffcount++;
            ciniaux++;
        }
    }

    if (Debug > 0) {
        printf("(prank:%d)Buffer CARGADO!!!\n", prank);
        fflush(stdout);
    }
}

void cargar_resultado_C(float *buffer, int fini, int cantFilas, int cini, int cantCols) {
    int i, j, k;
    k = 0;
    for (i = fini; i < (fini + cantFilas); i++) {
        for (j = cini; j < (cini + cantCols); j++) {
            int pos = j + (N * i);
            MatrizC[pos] = buffer[k];
            k++;
        }
    }
}

void multiplicar_bloque(float *buffer, float *Salida, int BuffCut, int prank,int tarea) {
    int i, j, k, l, m, z;

    float current;

    m = 0;
    for (i = 0; i < BuffCut; i = i + N) { // filas de A
        for (j = BuffCut; j < (BuffCut * 2); j = j + N) { //columnas de B
            l = i;
            k = j;
            current = 0;
            Salida[m] = 0;

            for (z = 0; z < N; z++) {
                float aux = (buffer[k])*(buffer[l]);
                current = current + aux;
                l++;
                k++;
            }

            Salida[m] = current;
            m++;
        }
    }
}

void imprimir_buffer(float *buffer) {
    int i;

    for (i = 0; i < SubBlockSize * N * 2; i++) {
        printf("%f ", buffer[i]);
        fflush(stdout);
        if (i % N == 0) {
            printf("\n");
            fflush(stdout);
            if (i == SubBlockSize * N) {
                printf("\n");
                fflush(stdout);
            }
        }
    }
}
