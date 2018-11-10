EXE = gppKer_gpuComplex.ex.div.cmplx.nofma.iw6
SRC1 = gppKer_gpuComplex.cpp 
SRC2 = GPUComplex.cu
#EXE = gppKer_double.ex
#SRC1 = gppKer_double.cpp
#SRC2 = GPUComplex_double.cu


#CXX = xlc++
#CXX = CC 
#CXX = g++
#CXX = clang++
CXX = nvcc

LINK = ${CXX}

ifeq ($(CXX),nvcc)
	CXXFLAGS= -g -O3 -std=c++11 -Wno-deprecated-gpu-targets -arch=sm_70 
#	CXXFLAGS+=-Xptxas -v --maxrregcount=150 #output usage of registers
	CXXFLAGS+=-fmad=true
	LINKFLAGS=-Wno-deprecated-gpu-targets
endif 

ifeq ($(CXX),g++)
	CXXFLAGS= -g -O3 -std=c++11 -fopenmp -foffload="-lm" -foffload=nvptx-none
	LINKFLAGS=-fopenmp
endif 

ifeq ($(CXX),xlc++)
	CXXFLAGS=-O3 -std=gnu++11 -g -qsmp=noauto:omp -qoffload #-Xptxas -v
	LINKFLAGS=-qsmp=noauto:omp -qoffload 
endif 

ifeq ($(CXX),clang++)
	CXXFLAGS=-O3 -std=gnu++11 -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda --cuda-path=${CUDA_HOME}
	LINKFLAGS=-fopenmp -fopenmp-targets=nvptx64-nvidia-cuda --cuda-path=${CUDA_HOME}
endif 

ifeq ($(CXX),icc)
	CXXFLAGS=-O3 -qopenmp -qopt-report=5
	CXXFLAGS+=xCORE_AVX2
#	CXXFLAGS+=-xMIC_AVX512
	LINKFLAGS=-qopenmp
endif 

OBJ1 = $(SRC1:.cpp=.o)
OBJ2 = $(SRC2:.cu=.o)

$(EXE): $(OBJ1) $(OBJ2) 
	$(CXX) $(OBJ1) $(OBJ2) -o $(EXE) $(LINKFLAGS)

$(OBJ1): $(SRC1) 
	$(CXX) -c $(SRC1) $(CXXFLAGS)

$(OBJ2): $(SRC2) 
	$(CXX) -c $(SRC2) $(CXXFLAGS)

clean: 
	rm -f $(OBJ1) $(OBJ2) $(EXE)

