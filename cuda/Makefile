verbose=1
dbg=0
fastmath=0

# Compilers
NVCC       = nvcc 

CUDA_GLOBAL_INSTALL = /usr/local
INCLUDES  += -I. -I$(CUDA_GLOBAL_INSTALL)/include -I./common/inc -I/share/apps/gsl-1.9/include/ -I/share/apps/cfitsio/include/
LIB	   = ${LIB} -L$(CUDA_GLOBAL_INSTALL)/lib -L./lib -L./common/lib -lcuda -lcudart -L/usr/lib

#use of cuda in the program
#CFLAGS 	+= -DUSE_CUDA
#NVCCFLAGS 	+= -maxrregcount 42 

# Device emulation configuration
ifeq ($(emu), 1)
	NVCCFLAGS   += -deviceemu -D__DEVICE_EMULATION__
	CXXFLAGS    += -D__DEVICE_EMULATION__ -ffloat-store
	CFLAGS	    += -D__DEVICE_EMULATION__ -ffloat-store
endif

ifeq ($(fastmath), 1)
	NVCCFLAGS += -use_fast_math
endif

CUDALIBFLAGS := -lcuda -lcudart -lGL -lGLU
CUDAEXTRAS   := ./common/obj/release/bank_checker.cpp_o \
	./common/obj/release/cmd_arg_reader.cpp_o \
	./common/obj/release/cutil.cpp_o \
	./common/obj/release/cutil_interop.cpp_o \
	./common/obj/release/error_checker.cpp_o \
	./common/obj/release/stopwatch.cpp_o \
	./common/obj/release/stopwatch_linux.cpp_o 

all: cmc_cuda.cu_o

%.cpp_o: ./common/src/%.cpp
	$(CXX) $(CFLAGS) -o $@ -c $<

%.cu_o : %.cu 
	$(NVCC) -arch sm_13 $(NVCCFLAGS) -o $@ -c $< $(INCLUDES)

clean:
	rm -f cmc_cuda.cu_o
