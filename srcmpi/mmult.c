/*
 * mmult.c: matrix multiplication using MPI.
 * There are some simplifications here. The main one is that matrices B and C
 * are fully allocated everywhere, even though only a portion of them is
 * used by each processor (except for processor 0)
 */

#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define SIZE 8			/* Size of matrices */

int A[SIZE][SIZE], B[SIZE][SIZE], C[SIZE][SIZE];

/* Se carga la matriz con valores desde 0 hasta SIZE*SIZE */
void fill_matrix(int m[SIZE][SIZE])
{
  static int n=0;
  int i, j;
  for (i=0; i<SIZE; i++)
    for (j=0; j<SIZE; j++)
      m[i][j] = n++;
}

/* Se imprime la matriz recibida como par�metro */
void print_matrix(int m[SIZE][SIZE])
{
  int i, j = 0;
  for (i=0; i<SIZE; i++) {
    printf("\n\t| ");
    for (j=0; j<SIZE; j++)
      printf("%2d ", m[i][j]);
    printf("|");
  }
}


int main(int argc, char *argv[])
{
  int myrank, P, from, to, i, j, k;
  int tag = 555;		/* any value will do */
  MPI_Status status;

  MPI_Init (&argc, &argv);

  /* Se estable el identificador del proceso llamador en el comunicador */
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  /* Numero de procesos que pertenecen al comunicador */
  MPI_Comm_size(MPI_COMM_WORLD, &P);

  /* Just to use the simple variants of MPI_Gather and MPI_Scatter we */
  /* impose that SIZE is divisible by P. By using the vector versions, */
  /* (MPI_Gatherv and MPI_Scatterv) it is easy to drop this restriction. */

  if (SIZE%P!=0) {
    if (myrank==0) printf("Matrix size not divisible by number of processors\n");
    MPI_Finalize();
    exit(-1);
  }
/* Se definen los limites */
  from = myrank * SIZE/P;
  to = (myrank+1) * SIZE/P;

  /* Process 0 fills the input matrices and broadcasts them to the rest */
  /* (actually, only the relevant stripe of A is sent to each process) */

/* El proceso maestro (0) carga las matrices A y B */
  if (myrank==0) {
    fill_matrix(A);
    fill_matrix(B);
  }

/* MPI_Bcast envia B del proceso maestro a todos los procesos */
  MPI_Bcast (B, SIZE*SIZE, MPI_INT, 0, MPI_COMM_WORLD);

/* MPI_Scatter envia partes diferentes de A del proceso fuente 0, a cada proceso
 * incluyendose el mismo. El dato recibido es almacenado en el buffer (A[from]) .
 * El proceso i recibe SIZE*SIZE/P elementos contiguos del tipo MPI_INT empezando
 * desde la posicion i*(SIZE*SIZE/P) del buffer enviado del proceso fuente */

 MPI_Scatter (A, SIZE*SIZE/P, MPI_INT, A[from], SIZE*SIZE/P, MPI_INT, 0, MPI_COMM_WORLD);

/* Se realiza el computo local en cada proceso */
  printf("computing slice %d (from row %d to %d)\n", myrank, from, to-1);
  for (i=from; i<to; i++)
    for (j=0; j<SIZE; j++) {
      C[i][j]=0;
      for (k=0; k<SIZE; k++)
	C[i][j] += A[i][k]*B[k][j];
    }

  /* Proceso inverso al MPI_Scatter. Cada proceso, incluyendo el objetivo (0)
   * envia el dato almacenado en el buffer de envio (C[from]), al proceso 0.
   * El dato es almacenado en el buffer de recepcion C del proceso 0 en el orden
   * en que fueron enviados a cada proceso. Es decir, el dato del proceso con
   * rank i es almacenado en C empezando en i* cantidadDeDatosEnviados.
   * El dato enviado por cada proceso debe ser del mismo tama�o y tipo */

  MPI_Gather (C[from], SIZE*SIZE/P, MPI_INT, C, SIZE*SIZE/P, MPI_INT, 0, MPI_COMM_WORLD);

/* El proceso maestro imprime las matrices */
  if (myrank==0) {
    printf("\n\n");
    print_matrix(A);
    printf("\n\n\t       * \n");
    print_matrix(B);
    printf("\n\n\t       = \n");
    print_matrix(C);
    printf("\n\n");
  }
/* Realiza tareas de limpieza para finalizar el entorno de ejecuci�n */
  MPI_Finalize();
  return 0;
}




