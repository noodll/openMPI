#include "commons.h"

void multiplicar_secuencial(float *A, float *B, float *C) {
    int i, j, k = 0;

    for (i = 0; i < M; i++) {
        for (j = 0; j < R; j++) {
            C[(i * R) + j] = 0;
            //printf("\n C%d = ", (i*R)+j);
            for (k = 0; k < N; k++) {
                C[(i * R) + j] = C[(i * R) + j] + (A[(i * N) + k] * B[(k * R) + j]);
                //printf(" (A%d * B%d) +", (i*N)+k, (k*R)+j);
            }
            //printf(" = %f ", C[(i*R)+j]);
        }
    }
}

void matrizAleatoria(FILE *archivo, int Filas, int Columnas) {

    int i, j;
    float random_value;
    char *linea = malloc(sizeof (char) * 20);
    char *cabecera = malloc(sizeof (char) * 20);

    fprintf(archivo, "%%%MatrixMarket matrix coordinate real general\n");
    int cant = Filas * Columnas;
    sprintf(cabecera, "%d %d %d\n", Filas, Columnas, cant);
    fprintf(archivo, cabecera);

    srand((int) ( tiempo_milis()/1000 ) );
   
    for (i = 0; i < Filas; i++) {
        for (j = 0; j < Columnas; j++) {
            int mivalor = 1 + (int) (100.0 * rand() / (RAND_MAX + 1.0));
            sprintf(linea, "%d %d %d\n", i, j, mivalor);
            fprintf(archivo, linea);
        }
    }
    free(linea);
    free(cabecera);
}

void inicializarMatriz(float *Matriz, int Filas, int Columnas) {
    int i;
    for (i = 0; i < Filas * Columnas; ++i) {
        Matriz[i] = 0.0;
    }
}

struct timeval obtener_tiempo() {
    struct timeval tv;

    if (gettimeofday(&tv, NULL) == -1) {
        fprintf(stderr, "medir_tiempo(): Error al obtener el tiempo\n");
        perror("sys_msg");
        exit(1);
    }

    return tv;
}

double tiempo_seg() {
    double seg, mseg;
    struct timeval tv = obtener_tiempo();

    seg = tv.tv_sec;
    mseg += (tv.tv_usec / 1000000.0);

    return (seg + mseg);
}

long long tiempo_milis() {
    long long seg, mseg;
    struct timeval tv = obtener_tiempo();

    // Guardamos los tiempos en long long int
    seg = tv.tv_sec;
    mseg = tv.tv_usec;

    // Convertimos los tiempos en milisegundos
    seg = seg * 1000;
    mseg = mseg / 1000;

    // Retornamos el resultado
    return (seg + mseg);
}

void verificar_longitud_path(char *path) {
    int maxlen = strlen(path);

    if (maxlen > PATH_BUFF_SIZE) {
        fprintf(stderr, "Path '%s' demasiado largo (%d, max=%d)",
                path, strlen(path), 200);
        exit(1);
    }
}

char *verificar_barra_final(char *path) {
    int len = strlen(path);

    if (path[len - 1] == '/')
        return "";

    return "/";
}

int procesar_parametros(int argc, char ** argv, char ** alg) {
    int nopt = 1;
    int f1_visto = 0;
    int f2_visto = 0;
    int s_visto = 0;
    int t_visto = 0;
    int o_visto = 0;
    int v_visto = 0;
    int n_visto = 0;
    int j_visto = 0; 
    int opciones = 0;
    int p_visto = 0;

    // Valores por defecto

    FA = "matrizA.mtx";
    FB = "matrizB.mtx";
    Debug = 0; 
    genMat = 0; 
    printResult = 0;

    while (1) {
        if (argc > nopt && strcmp(argv [nopt], "-f1") == 0) {
            if (f1_visto != 0)
                fprintf(stderr, "-f1 solo puede especificarse"
                    "una vez\n");
            f1_visto = 1;
            FA = argv [nopt + 1];
            nopt += 2;

        } else if (argc > nopt && strcmp(argv [nopt], "-f2") == 0) {
            if (f2_visto != 0)
                fprintf(stderr, "-f2 solo puede especificarse"
                    "una vez\n");
            f2_visto = 1;
            FB = argv [nopt + 1];
            nopt += 2;

        } else if (argc > nopt && strcmp(argv [nopt], "-s") == 0) {
            if (s_visto != 0)
                fprintf(stderr, "-s solo puede especificarse"
                    "una vez\n");
            s_visto = 1;
            FC = argv [nopt + 1];
            nopt += 2;

        } else if (argc > nopt && strcmp(argv [nopt], "-t") == 0) {
            
            printf("Parámetro de tipo de algoritmo: ....");
            if (t_visto != 0)
                fprintf(stderr, "-t solo puede especificarse"
                    "una vez\n");
            t_visto = 1;
            *alg = malloc(sizeof (char) * 3);
            *alg = argv [nopt + 1];
            nopt += 2;
            opciones++;
        } else if (argc > nopt && strcmp(argv [nopt], "-o") == 0) {
            if (o_visto != 0)
                fprintf(stderr, "-o solo puede especificarse"
                    "una vez\n");
            o_visto = 1;
            FCSV = argv [nopt + 1];
            nopt += 2;            
        } else if (argc > nopt && strcmp(argv [nopt], "-v") == 0) {
            if (v_visto != 0)
                fprintf(stderr, "-v solo puede especificarse"
                    "una vez\n");
            v_visto = 1;
            Debug = 1;
            nopt += 1;
        } else if (argc > nopt && strcmp(argv [nopt], "-n") == 0) {
            if (n_visto != 0)
                fprintf(stderr, "-n solo puede especificarse"
                    "una vez\n");
            n_visto = 1;
            N = atoi(argv[nopt + 1]);
            nopt += 2;            
        } else if (argc > nopt && strcmp(argv [nopt], "-p") == 0) {
            if (p_visto != 0)
                fprintf(stderr, "-v solo puede especificarse"
                    "una vez\n");
            
            if (j_visto == 0)
                fprintf(stderr, "-p solo puede especificarse"
                    "si -j fue especificado antes\n");
            
            p_visto = 1;
            printResult = 1; 
            nopt += 1;
        } else if (argc > nopt && strcmp(argv [nopt], "-j") == 0) {
            if (j_visto != 0)
                fprintf(stderr, "-j solo puede especificarse"
                    "una vez\n");
            j_visto = 1;
            jointResult = 1; 
            nopt += 1;
        } else {
            break;
        }

    }
    
    if (n_visto == 1) {
        genMat = 1; 
    }
    
    if (opciones != 1) {
        uso();

        free(MatrizA);
        free(MatrizB);
        free(MatrizC);
        MPI_Finalize();
        exit(1);
    }
    return nopt;
}

void uso(void) {

    printf("El programa debe recibir: \n");
    printf("  ./mmulfpunaMain -f1 <fileA> -f2 <fileB> -s <FileC> -t [c|d|2d|2dd]\n\n -q n -o FileCSV");
    printf("  File A      : Archivo que contiene la matriz A");
    printf("  File B      : Archivo que contiene la matriz B\n");
    printf("  File C      : Archivo de salida para matriz resultado (opcional)\n");    
    printf("  -n N        : Especificar tamaño de matriz (en este caso la matriz se genera y no se leen)\n");
    printf("  -t c        : Ejecutar algoritmo de Cannon\n");
    printf("  -t d        : Ejecutar algoritmo DNS\n");
    printf("  -t 2d       : Ejecutar algoritmo 2D con particionamiento por bloque cíclico\n");
    printf("  -t 2dd      : Ejecutar algoritmo 2D Diagonal\n");
    printf("  -o FileCSV  : Archivo CSV de salida para guardar resultados de mediciones (opcional: salida.csv por defecto)\n\n");
    printf("  -v          : Imprimir información de depuración (opcional: desactivado por defecto)\n\n");
    printf("  -p          : Imprimir la matriz de resultado final. Solo con -j.\n\n");
    printf("  -j          : Combinar resultados desde todos los procesos\n\n");
    printf("  Los archivos deben estar en formato de Matrix Market Exchange\n");

}

float *asignar_mem_matriz(int Size) {
    return malloc(sizeof (float) * Size);
}


void leer_matrices_aleatorio(float *A, float *B){
    
    int itmp,jtmp = 0;
    printf("Tamaño a Generar: %d\n",N);
    srand(tiempo_milis()/1000);
    
    for(itmp = 0; itmp < N ; itmp++){
        for(jtmp = 0; jtmp < N ; jtmp++){    
            A[itmp*N+jtmp] = 1 + (int) (100.0 * rand() / (RAND_MAX + 1.0));
        }
    }
    for(itmp = 0; itmp < N ; itmp++){
        for(jtmp = 0; jtmp < N ; jtmp++){    
            B[itmp*N+jtmp] = 1 + (int) (100.0 * rand() / (RAND_MAX + 1.0));
        }
    }
}

void leer_matrices(){
    int itmp,jtmp = 0;
    for(itmp = 0; itmp < N ; itmp++){
        for(jtmp = 0; jtmp < N ; jtmp++){    
            A[(itmp * N)+jtmp] = (itmp * N) + jtmp;
        }
    }
    for(itmp = 0; itmp < N ; itmp++){
        for(jtmp = 0; jtmp < N ; jtmp++){    
            B[(itmp * N)+jtmp] = (itmp * N) + jtmp;
        }
    }
}

void acumular_multiplicar(float A[tamanho_bloque][tamanho_bloque], 
                          float B[tamanho_bloque][tamanho_bloque], 
                          float C[tamanho_bloque][tamanho_bloque]){

    int i,j,k;
    float tmp;
    for(i=0;i<tamanho_bloque;i++){
        for(j=0;j<tamanho_bloque;j++){    
            tmp = C[i][j];
            for(k=0;k<tamanho_bloque;k++){    
                tmp = tmp + (A[i][k] * B[k][j]);
            }
            C[i][j] = tmp;
        }
    }
    
}

void cerar_matriz(float *mat){

    int i,j;
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){    
            mat[(i*N) + j] = 0.0;
        }
    }
    
}

void cerar_bloque(float blo[tamanho_bloque][tamanho_bloque]){

    int i,j;
    for(i=0;i<tamanho_bloque;i++){
        for(j=0;j<tamanho_bloque;j++){    
            blo[i][j] = 0.0;
        }
    }
    
}

void imprimirBloque(float bloq[tamanho_bloque][tamanho_bloque], int tam){
    int i,j= 0;
    for(i= 0; i < tam ; i++){
        for(j= 0; j < tam ; j++){    
            printf("%f\t", bloq[i][j]);
            fflush(stdout);
        }
        printf("\n");
        fflush(stdout);
    }
}

