MPICC = cc
CFLAGS = -g -O3 -dynamic -qopenmp -restrict -qopt-streaming-stores always  \
         -DSTREAM_ARRAY_SIZE=400000000 -DNTIMES=50 \
         -I$(VTUNE_AMPLIFIER_XE_2018_DIR)/include
LDFLAGS = -L$(VTUNE_AMPLIFIER_XE_2018_DIR)/lib64 -littnotify

stream_mpi.exe: stream_mpi.c Makefile
	$(MPICC) $(CFLAGS) stream_mpi.c -o stream_mpi.exe $(LDFLAGS)

clean:
	rm -f stream_mpi.exe
