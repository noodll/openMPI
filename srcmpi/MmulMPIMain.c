/* 
 * File:   2dblockMain.c
 * Author: cparra
 *
 * Created on 10 de julio de 2008, 09:05 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <mpi.h>

#include "commons.h"
#include "salida-metricas.h"
#include "mmio.h"
#include "2dblock.h"
#include "2ddiagonal.h"

/*
 * Código principal del Proyecto. 
 *
 * Inicializa la comunicación MPI y los parámetros iniciales. 
 *
 * En base al tipo de algoritmo seleccionado por parámetros, ejecuta el 
 * algoritmo correspondiente. 
 */
int main(int argc, char** argv) {

    int prank, i, j, k, n1, n2, namelong; // variables del proceso
    char pname[MPI_MAX_PROCESSOR_NAME]; // nombre del proceso
    result_t r; // estructura de resultados

    /* 
     * Reserva de memoria para viables de parámetros de nodo 0
     */

    FA = malloc(sizeof (char) * PATH_BUFF_SIZE);
    FB = malloc(sizeof (char) * PATH_BUFF_SIZE);
    FC = malloc(sizeof (char) * PATH_BUFF_SIZE);
    FCSV = malloc(sizeof (char) * PATH_BUFF_SIZE);
    char * algoritmo = malloc(sizeof (char) * 4);

    /* 
     * Valores por defecto para parámetros opcionales de nodo 0
     */
    FC = "matrizC.mtx";
    /* 
     * Bloque global de inicialización MPI, común a todos los algoritmos. 
     */
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &prank); // obtenemos el id del proceso
    MPI_Comm_size(MPI_COMM_WORLD, &P); // total de procesos
    MPI_Get_processor_name(pname, &namelong);
    procesar_parametros(argc, argv, &algoritmo);

    /* Funciones especiales del proceso 0 
     * 
     *
     * El proceso maestro (0) carga las matrices A y B por ejemplo. 
     */

    if (prank == 0) {

        /* Lectura de las matrices. Se mantiene la lectura de los parámetros 
         * M, n1 y n2, y R para la futura implementación de algoritmos paralelos
         * para matrices no cuadradas. 
         */
        if (genMat == 0) { // Leer matrices de Archivo
            leer_matriz(FA, &MatrizA, &M, &n1);
            leer_matriz(FB, &MatrizB, &n2, &R);

            N = M; // por el momento se asume que las matrices serán cuadradas. 

            if (n1 != n2) {
                printf("Maatrices no compatibles");
                free(MatrizA);
                free(MatrizB);
                MPI_Finalize();
                exit(EXIT_FAILURE);
            } else if (n1 != M) {
                printf("Matrices no cuadradas");
                free(MatrizA);
                free(MatrizB);
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }
        } else {
            MatrizA = malloc(sizeof (float) * N * N);
            MatrizB = malloc(sizeof (float) * N * N);
            leer_matrices_aleatorio(MatrizA, MatrizB);
        }
        // Reserva de matriz resultado
        MatrizC = malloc(sizeof (float) * N * N);
    }

    if (genMat == 0) {
        int Nbuff[1];
        if (prank == 0) {
            Nbuff[0] = N;
        }
        MPI_Bcast(Nbuff, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (prank != 0) {
            N = Nbuff[0];
        }
    }

    M = R = N; 
    /** 
     * Ejecución de algoritmos Específicos. 
     *
     */
    if (strcmp(algoritmo, "2d") == 0) {
        printf("\n(prank:%d)Se ejecutará el algoritmo 2D con particionamiento cíclico por bloques.\n", prank);
        sprintf(r.testId, "%s", "2D");

        /*
         * Condiciones iniciales a cumplir: 
         * 1. N es un cuadrado perfecto. 
         * 2. SubBlockSize será la raíz de N para tener el balalceo más preciso. 
         * 3. El cuadrado de P debe ser menor que SubBlockSize: sqrt(P) < sqrt(N). 
         * 4. P es el número de procesos, y debe ser cuadrado perfecto. 
         * 5. Se admiten matrices para las cuales N % sqrt(P) != 0, solo que en ese caso
         *    se viola en cierta medida la restricción del particionamiento por bloque. 
         */

        if (N % N != 0) { // revisión de condición 1.
            printf("N no es cuadrado perfecto\n");
            free(MatrizA);
            free(MatrizB);
            free(MatrizC);
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }

        if (P % P != 0) { // revisión de condición 2. 
            printf("P no es cuadrado perfecto\n");
            free(MatrizA);
            free(MatrizB);
            free(MatrizC);
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }

        if (sqrt(P) > sqrt(N)) { // revisión de condición 3. 
            printf("P no es menor o igual que el cuadrado de N (%d, %f)\n", P, sqrt(N));
            free(MatrizA);
            free(MatrizB);
            free(MatrizC);
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }

        SubBlockSize = sqrt(N);
        ProcBlockSize = sqrt(P);

        MEDICION(ejecutar_2d(prank));

    } else if (strcmp(algoritmo, "2dd") == 0) {
        /*
         * Condiciones iniciales a cumplir: 
         * 1. N es un cuadrado perfecto.  
         * 2. P es un cuadrado perfecto 
         * 3. P debe ser como máximo igual a NxN
         * 4. N debe ser divisible por el cuadrado de P
         */

        if (N % N != 0) { // revisión de condición 1.
            printf("(%d,%d) ---> N no es cuadrado perfecto\n", mycoords[0], mycoords[1]);
            free(MatrizA);
            free(MatrizB);
            free(MatrizC);
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }

        if (P % P != 0) { // revisión de condición 2. 
            printf("(%d,%d) ---> P no es cuadrado perfecto\n", mycoords[0], mycoords[1]);
            free(MatrizA);
            free(MatrizB);
            free(MatrizC);
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }

        if (N % ((int) sqrt(P)) != 0) { // revisión de condición 4. 
            printf("(%d,%d) ---> N no es divisible por el cuadrado de P\n", mycoords[0], mycoords[1]);
            free(MatrizA);
            free(MatrizB);
            free(MatrizC);
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }

        /*********************************************************************
         * PASOS A SEGUIR POR EL ALGORITMO:                                  *
         *                                                                   *
         * 0. Inicializar Parámetros.                                        *
         *    - Asignar memoria a Buffers                                    *
         *    - Crear topología 2D de procesos.                              *
         *                                                                   *
         * 1. Realizar la distribución inicial en la diagonla principal de   *
         *    la malla.                                                      *
         *    - Desde P0, enviar a todos los algoritmos de la malla, iésimo  *
         *      grupo de columnas y filas de A y B.                          *
         *                                                                   *
         * 2. Ejecutar la lógica principal del Algoritmo.                    *
         *    - Broadcast de A*,j a todos los procesos p*,j desde la diagonal*
         *    - Send de Bik a cada pkj                                       *
         *    - Recibir A*,j y Bji                                           *
         *    - Calcular producto externo I*,i = A*,j x Bji                  *
         *    - Enviar I*,i en la dirección i hasta Pi,i.                    *
         *    - En cada Pi,i, recibir I*,i y calcular C*,i = C*,i + I*,i     *
         *                                                                   *
         *3. Enviar todos los C*,i calculados en los procesos diagonales al  *
         *   proceso P0,0                                                    *
         *
         * Observación: Solo el paso 2 se incluye en la medición de tiempo ya*
         * que solo esta parte constituye la lógica del algoritmo según el   *
         * paper de Gupta. 
         *********************************************************************/

        /*********************************************************************
         * PRELIMINARES --> Inicialización de Variables de uso global        *
         *********************************************************************/

        Grupos = sqrt(P);
        GrupSize = (int) (floor(N / Grupos));
        IthGrupoA = malloc(sizeof (float) * GrupSize * N); // para albergar N cols de A
        IthGrupoB = malloc(sizeof (float) * GrupSize * N); // para albergar N filas de B
        BloqueB = malloc(sizeof (float) * GrupSize * GrupSize); // bloque de B a multiplicar por IthGrupoA
        IProduct = malloc(sizeof (float) * GrupSize * N); // para albergar producto externo de IthGrupoA x BloqueB
        IProductReduced = malloc(sizeof (float) * GrupSize * N); // para albergar producto externo recibido en Reduce

        //--> Creación de la topología

        Dims[0] = sqrt(P);
        Dims[1] = sqrt(P);

        Periods[0] = 0;
        Periods[1] = 0;

        MPI_Cart_create(MPI_COMM_WORLD, 2, Dims, Periods, 0, &Comm2d);
        MPI_Cart_coords(Comm2d, prank, 2, mycoords);
        MPI_Cart_rank(Comm2d, mycoords, &myrank);

        printf("(%d,%d) ---> Se ejecutara el algoritmo 2D Diagonal.\n", mycoords[0], mycoords[1]);
        fflush(stdout);
        sprintf(r.testId, "%s", "2DD");

        /*********************************************************************
         * PASO #1 --> Distribución Inicial.                                 *
         * - Se envía a cada proceso ubicado en la diagonal de la topología  *
         *   el iésimo grupo de columnas de A y el iésimo grupo de filas de  *
         *   la matriz B.                                                    *
         *********************************************************************/

        // Debug Info
        if (Debug > 0) {
            printf("(%d,%d) --> Iniciando distribución Inicial\n", mycoords[0], mycoords[1]);
            fflush(stdout);
        }

        distribucion_inicial();

        // Debug Info
        if (Debug > 0) {
            printf("(%d,%d) --> Distribución Inicial finalizada. Inicio del Algoritmo\n", mycoords[0], mycoords[1]);
            fflush(stdout);
        }

        /*********************************************************************
         * PASO #2 --> Lógica principal del algoritmo                        *
         *********************************************************************/
        MEDICION(ejecutar2dd());

        printf("(%d,%d) --> Algoritmo Finalizado. Combinando resultados. \n", mycoords[0], mycoords[1]);
        fflush(stdout);

        /*********************************************************************
         * PASO #3 --> Lógica principal del algoritmo                        *
         *********************************************************************/
        if (jointResult == 1) {
            combinar_resultados();
        }

        printf("(%d,%d) --> Combinación de resultados FINALIZADA. Guardando mediciones \n", mycoords[0], mycoords[1]);
        fflush(stdout);

        free(IthGrupoA);
        free(IthGrupoB);
        free(BloqueB);
        free(IProduct);
        free(IProductReduced);

        printf("(%d,%d) --> Memoria liberada por 2dd\n", mycoords[0], mycoords[1]);
        fflush(stdout);

        /*********************************************************************
         * FIN DE 2D DIAGONAL
         *********************************************************************/
    } else if (strcmp(algoritmo, "seq") == 0) {
        printf("(%d,%d) ---> Se ejecutara el algoritmo secuencial.\n", mycoords[0], mycoords[1]);
        fflush(stdout);
        sprintf(r.testId, "%s", "SEQ");
        MEDICION(multiplicar_secuencial(MatrizA, MatrizB, MatrizC)); 
    }


    printf("(prank:%d) --> Guardando mediciones \n", prank);
    fflush(stdout);

    double buffTiempo[1];

    if (prank == 0) {
        int proc;
        sprintf(FCSV, "%s-%s.csv", "salida", r.testId);

        FILE *fcsv = open_csvFile_Append(FCSV);
        print_header(fcsv, FCSV);

        r.T = t3;
        r.q = P;
        r.size = N;
        r.prank = 0;
        print_test_result(fcsv, &r);


        printf("(prank:%d) --> Recibiendo mediciones \n", prank);
        fflush(stdout);

        for (proc = 1; proc < P; proc++) {
            MPI_Status stat;

            MPI_Recv(buffTiempo, 1, MPI_DOUBLE, proc, proc * TAGMULT*TAGMULT, MPI_COMM_WORLD, &stat);

            printf("Se recibe %f seg de %d \n", buffTiempo[0], proc);
            fflush(stdout);

            double tmp = buffTiempo[0];
            r.T = tmp;
            r.q = P;
            r.size = N;
            r.prank = proc;

            print_test_result(fcsv, &r);
        }

        close_csvFile(fcsv);
    } else {

        printf("(prank:%d) --> Enviando mediciones \n", prank);
        fflush(stdout);

        buffTiempo[0] = t3;
        MPI_Send(buffTiempo, 1, MPI_DOUBLE, 0, prank * TAGMULT*TAGMULT, MPI_COMM_WORLD);
    }

    /* 
     * El Proceso 0 realiza las actividades de Finalización. 
     * - Imprime las matrices. 
     * - Finaliza MPI
     *
     */
    if (prank == 0 && printResult == 1 && jointResult == 1) {
        printf("(prank:%d) --> Algoritmo %s Terminado!\n\n", prank, r.testId);
        fflush(stdout);

        printf("(prank:%d) --> La multilplicacion de matrices duró %4f segundos\n", prank, t3);
        fflush(stdout);

        if (prank == 0) {
            printf("\n\n");
            imprimir_matriz(MatrizA, N, N);
            printf("\n\n\t       * \n");
            imprimir_matriz(MatrizB, N, N);
            printf("\n\n\t       = \n");
            imprimir_matriz(MatrizC, N, N);
            printf("\n\n");
        }

    }
    
    if (prank == 0) {
        free(MatrizA);
        free(MatrizB);
        free(MatrizC);

        printf("(prank:%d) --> Matrices liberadas...\n", prank);
        fflush(stdout);
    }

    MPI_Finalize();
    return (MPI_SUCCESS);
}

