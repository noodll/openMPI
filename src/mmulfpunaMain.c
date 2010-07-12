/***********************************************************************
* mmul1D.c                                                             *
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

#include "commons.h"
#include "matriz1d.h"
#include "../../mmio/mmio.h"
#include "Matrix2D.h"
#include "salida-metricas.h"

char *FA;       // path del archivo que contiene la matriz A
char *FB;       // path del archivo que contiene la matriz B
char *FC;       // path del archivo que contiene la matriz C
char *FCSV;     // path del archivo en que se escribirán los resultados de las métricas
int q;          // cantidad de hilos a ejecutar

void uso(void);
int procesar_parametros (int argc, char ** argv, char ** alg);

int ejecutar_1d(int q,int *overhead);
int ejecutar_2d(int q,int *overhead);
int ejecutar_secuencial();
void uso_memoria();

int main (int argc, char **argv) {

	//Reserva de memoria
	FA   = malloc(sizeof(char) * PATH_BUFF_SIZE);
	FB   = malloc(sizeof(char) * PATH_BUFF_SIZE);
	FC   = malloc(sizeof(char) * PATH_BUFF_SIZE);
	FCSV = malloc(sizeof(char) * PATH_BUFF_SIZE);
	char * algoritmo = malloc(sizeof(char) * 4);
    	long long t1, t2, t3;  
	int n1, n2,overhead;
	result_t r;
    
	system("clear");

	// Valores por defecto para parámetros opcionales
	FCSV = "salida.csv";
	q = 9;
	FC = "matrizC.mtx";
		
	// Procesamiento de los parametros de la linea de comandos.
	procesar_parametros(argc, argv, &algoritmo);
    
	//Lectura de las matrices
	leer_matriz(FA, &MatrizA, &M, &n1);
	leer_matriz(FB, &MatrizB, &n2, &R);
	
	//TODO Verificar espacio disponible en memoria (debe ser antes de leer matrices)
	uso_memoria();
	
	if(n1 != n2){
		printf("Las matrices especificadas no pueden multiplicarse");
		exit(1);
	}
	N = n1;


	// Reserva de matriz resultado
	MatrizC = malloc (sizeof(float) * M * R );
	
	r.size = (int) M*R;
	r.O = 0;
    
	if (strcmp(algoritmo,"s")==0) {
		printf("\nSe ejecutara el algoritmo secuencial.\n");
		sprintf(r.testId, "%s", "Seq");
		q = 1; 
		MEDICION (ejecutar_secuencial());		
	} else if(strcmp(algoritmo,"1d")==0) {
		printf("\nSe ejecutara el algoritmo 1D con %d hilos y matriz C de tamaño = %d.\n",q,r.size);
		sprintf(r.testId, "%s", "1D");
		MEDICION (ejecutar_1d(q,&overhead));
		r.O = overhead;
	} else if(strcmp(algoritmo,"2d")==0) {
		printf("\nSe ejecutara el algoritmo 2D con %d hilos y matriz C de tamaño = %d.\n",q,r.size);
		sprintf(r.testId, "%s", "2D");
		MEDICION (ejecutar_2d(q,&overhead));
		r.O = overhead;
	}

	r.q = q;
	r.T = t3;
	FILE *fcsv = open_csvFile_Append(FCSV);
	
	print_header(fcsv,FCSV);
	print_test_result(fcsv,&r);
	close_csvFile(fcsv);
	
	// TODO Cambiar impresión de tiempo por guardado de métrica
	printf ("\n\nTiempo de Ejecución total: %lld\n",t3);
	fflush(stdout);		

	//printf("\nMatrix C..\n\n"); 
	//imprimir_matriz(MatrizC, M,R);

	FILE * archivo_salida = fopen(FC, "w"); 
	imprimir_matriz_archivo(archivo_salida, MatrizC, M, R);
	
	return 0; 
}

int ejecutar_secuencial() {

    multiplicar_secuencial(MatrizA, MatrizB, MatrizC);
    return 0;
}

int ejecutar_1d(int q,int *overhead) {
	
	pthread_t *hilos;           // arreglo de hilos
	particion *particiones;     // arreglo de parámetros de cada partición de datos
	int i;
	long long o1, o2, o3;       // variables para medir tiempo de overhead 

	o1 = tiempo_milis();
	// Reservar memoria para los hilos y las particiones 
	hilos = (pthread_t *) malloc(q * sizeof(pthread_t));
	particiones =(particion *) malloc(sizeof( particion )*q );

	// Configurar Particiones
	crear_particiones1D(M,R,q,particiones);

	for (i = 0; i < q; ++i) {
		pthread_create(&hilos[i], NULL, runthread1D, (void *)(particiones+i));    
	}   

	o2 = tiempo_milis();
	o3 = o2 - o1; 

	// TODO Cambiar impresión de tiempo por guardado de métrica
	*overhead = o3;
	printf ("Overhead 1D: %lld",o3);
	fflush(stdout);

	for ( i = 0; i < q; ++i){
		pthread_join(hilos[i], NULL);
	}

	free(hilos);
	free(particiones);
	return 0;
}

int ejecutar_2d(int q, int *overhead)
{

	pthread_t *hilos;           // arreglo de hilos
	particion *particiones;     // arreglo de parámetros de cada partición de datos

	int i;
	long long o1, o2, o3; 

	o1 = tiempo_milis();
	// Reservar memoria para los hilos y las particiones 
	hilos = (pthread_t *) malloc(q * sizeof(pthread_t));
	particiones =(particion *) malloc(sizeof( particion )*q );

	// Configurar Particiones
	crear_particiones2D(M, R, q, particiones);
	imprimir_particiones2D(q, particiones);

	for (i = 0; i < q; i++) {
		pthread_create(&hilos[i], NULL, runthread2D, (void *)(particiones+i)); 
	}
	
	o2 = tiempo_milis();
	o3 = o2 - o1; 
	
	// TODO Cambiar impresión de tiempo por guardado de métrica
	*overhead = o3;
	printf ("Overhead 2D: %lld",o3);
	fflush(stdout);
	
	for ( i = 0; i < q; i++){
		pthread_join(hilos[i], NULL);
	}

	free(hilos);
	free(particiones);

	return 0;
}

int procesar_parametros (int argc, char ** argv, char ** alg)
{
	int nopt     = 1;
	int f1_visto = 0;
	int f2_visto = 0;
	int s_visto  = 0;
	int t_visto  = 0;
	int q_visto  = 0;
	int o_visto  = 0;
	int opciones = 0;
	
	while ( 1 ) {
		if ( argc > nopt && strcmp (argv [nopt], "-f1") == 0 ) {
			if ( f1_visto != 0 )
				fprintf(stderr, "-f1 solo puede especificarse"
						"una vez\n");
			f1_visto = 1;
			FA = argv [nopt+1];
			nopt += 2;
			opciones++;

		} else if (argc > nopt &&  strcmp (argv [nopt], "-f2") == 0 ) {
			if ( f2_visto != 0 )
				fprintf(stderr, "-f2 solo puede especificarse"
						"una vez\n");
			f2_visto = 1;
			FB = argv [nopt+1];
			nopt += 2;
			opciones++;

		} else if (argc > nopt &&  strcmp (argv [nopt], "-s") == 0 ) {
			if ( s_visto != 0 )
				fprintf(stderr, "-s solo puede especificarse"
						"una vez\n");
			s_visto = 1;
			FC = argv [nopt+1];
			nopt += 2;

		} else if (argc > nopt &&  strcmp (argv [nopt], "-t") == 0 ) {
			if ( t_visto != 0 )
				fprintf(stderr, "-t solo puede especificarse"
						"una vez\n");
			t_visto = 1;
			*alg = malloc(sizeof(char) * 3);
			*alg = argv [nopt+1];
			nopt += 2;
			opciones++;
		} else if (argc > nopt &&  strcmp (argv [nopt], "-o") == 0 ) {
			if ( o_visto != 0 )
				fprintf(stderr, "-0 solo puede especificarse"
						"una vez\n");
			o_visto = 1;
			FCSV = argv [nopt+1];
			nopt += 2;
		} else if (argc > nopt &&  strcmp (argv [nopt], "-q") == 0 ) {
			if ( q_visto != 0 )
				fprintf(stderr, "-q solo puede especificarse"
						"una vez\n");
			q_visto = 1;
			char *snum = malloc(sizeof(char) * 3);
			snum = argv [nopt+1];
			
			long int q_aux = strtol(snum,NULL,0);
			
			q = (int) q_aux;
			nopt += 2;

		} else {
			break;
		}

	}
	if ( opciones != 3 ){
		uso();
		exit(1);
        }
	return nopt;
}

void uso(void){

        printf("El programa debe recibir: \n");
        printf("  ./mmulfpunaMain -f1 <fileA> -f2 <fileB> -s <FileC> -t [s|1d|2d]\n\n -q n -o FileCSV");
        printf("  File A      : Archivo que contiene la matriz A");  
        printf("  File B      : Archivo que contiene la matriz B\n");  
        printf("  File C      : Archivo de salida para matriz resultado (opcional)\n");
        printf("  -t s        : Ejecutar multiplicación secuencial\n");
        printf("  -t 1d       : Ejecutar multiplicación 1d\n");
        printf("  -t 2d       : Ejecutar multiplicación 2d\n");
        printf("  -q n        : Utilizar n hilos (opcional: 9 por defecto)\n");
        printf("  -o FileCSV  : Archivo CSV de salida para guardar resultados de mediciones (opcional: salida.csv por defecto)\n\n");
        printf("  Los archivos deben estar en formato de Matrix Market Exchange\n");

}

void uso_memoria() {
	long int pageDisp  = get_avphys_pages();
	int pageSize       = getpagesize();

	// Uso del programa en variables
	int ptSize         = sizeof(pthread_t)*q;
	int doSize         = sizeof(double)*0;
	int charSize       = sizeof(char)*(PATH_BUFF_SIZE+4);
	int partSize       = sizeof(particion)*q;
	int longlongSize   = sizeof(long long)*3;
	int intSize   	   = sizeof(int)*3;
	
	//Uso de memoria de las matrices
	int flSize         = sizeof(float)*(M*R+M*N+N*R);
	
	int suma = ptSize + doSize + charSize + partSize + longlongSize + intSize + flSize;
	printf("Pags. Disponibles :  %ld \n",pageDisp);
	printf("Tam. de Pagina    :  %d \n",pageSize);
	printf("Mem. Disponible   :  %ld bytes\n",pageSize*pageDisp);
	printf("Memoria Requerida :  %d bytes\n",suma);
}


