all: main 
	
main: 
	mpicc MmulMPIMain.c commons.c 2dblock.c 2ddiagonal.c salida-metricas.c -I ../../mmio/ ../../mmio/mmio.c -lm -o mmult 

2D:
	mpicc -c 2dblock.c ../../mmio/mmio.c -lm 
	
commons:
	mpicc -c commons.c ../../mmio/mmio.c -lpthread -lm

cannon:
	mpicc algCannon.c commons.c salida-metricas.c ../../mmio/mmio.c -I ../../mmio/ -lm -o mmultCannon

dns:
	mpicc algDNS.c commons.c  salida-metricas.c ../../mmio/mmio.c -I ../../mmio/ -lm -o mmultDNS

matriz2d:
	mpicc -c Matrix2D.c commons.c ../../mmio/mmio.c -lpthread -lm

generar:
	mpicc generar.c commons.c -o genMat
	
clean:
	rm -rf *.o mmult genMat
