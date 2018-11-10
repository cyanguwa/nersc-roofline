##
# module load vtune
# module load sde
# module unload darshan
#
FC=ftn
CC=cc
FFLAGS= -g -dynamic -O3 -openmp
CFLAGS= -g -I $(VTUNE_AMPLIFIER_XE_2016_DIR)/include
LDFLAGS=$(VTUNE_AMPLIFIER_XE_2016_DIR)/lib64/libittnotify.a
EXE=jacobi.x
COBJ= api_itt_sde.o
F90OBJ= module_itt_sde.o
default: jacobi
clean:
	rm *.o *.mod $(EXE)
%.o:   %.c
	$(CC) $(CFLAGS) -c $<
%.o:    %.f90
	$(FC) $(FFLAGS) -c $<
jacobi: $(COBJ) $(F90OBJ) jacobi.o
	$(FC) $(FFLAGS) -o jacobi.x $(COBJ) $(F90OBJ) jacobi.o $(LDFLAGS)
