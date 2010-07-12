all: main 
	
main: 
	gcc mmulfpunaMain.c commons.c matriz1d.c Matrix2D.c salida-metricas.c ../../mmio/mmio.c -lpthread -lm -o mmult 

2D:
	gcc Matrix2D.c ../../mmio/mmio.c -lpthread -lm -o m2d
	
commons:
	gcc -c commons.c ../../mmio/mmio.c -lpthread -lm

matriz1d:
	gcc -c matriz1d.c commons.c ../../mmio/mmio.c -lpthread -lm

matriz2d:
	gcc -c Matrix2D.c commons.c ../../mmio/mmio.c -lpthread -lm

generar:
	gcc generar.c commons.c -o generar
	
clean:
	rm -rf *.o
