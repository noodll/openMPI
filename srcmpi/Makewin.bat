gcc ..\..\mmio\mmio.c salida-metricas.c commons.c algCannon.c -o mmultCannon -I C:\ap -I ..\..\mmio -L c:\ap -lmpich -lm

gcc ..\..\mmio\mmio.c salida-metricas.c commons.c algDNS.c -o mmultDNS -I C:\ap -I ..\..\mmio -L c:\ap -lmpich -lm

gcc ..\..\mmio\mmio.c salida-metricas.c commons.c 2dblock.c 2ddiagonal.c MmulMPIMain.c -o mmult -I C:\ap -I ..\..\mmio -L c:\ap -lmpich -lm
