### Makefile template based on https://bitbucket.org/berkeleylab/upcxx/src/master/example/cuda_vecadd/Makefile

ifeq ($(UPCXX_INSTALL),)
$(warning UPCXX_INSTALL environment variable is not set, assuming upcxx is in the PATH)
UPCXX=upcxx
UPCXXRUN=upcxx-run
else
ifeq ($(wildcard $(UPCXX_INSTALL)/bin/upcxx),)
$(error Please set UPCXX_INSTALL=/path/to/upcxx/install)
else
UPCXX=$(UPCXX_INSTALL)/bin/upcxx
UPCXXRUN=$(UPCXX_INSTALL)/bin/upcxx-run
endif
endif

NVCC ?= nvcc

# Honor user's NVCCARCH_FLAGS (if any) warning if NVCCARCH is also set.
# Otherwise convert NVCCARCH (if any) to an -arch=... option.
ifneq ($(strip $(NVCCARCH)),)
  ifneq ($(strip $(NVCCARCH_FLAGS)),)
    $(warning Both NVCCARCH and NVCCARCH_FLAGS are set.  Ignoring NVCCARCH.)
  else
    NVCCARCH_FLAGS = -arch=$(NVCCARCH)
  endif
endif

CXXFLAGS = -Iinclude -std=c++11 -O3
DEBUGFLAGS = -Iinclude -std=c++11 -g

CPUSRCS = src/main.cpp src/options.cpp
UPCXXINC = -I/home/kirtus/Software/upcxx/install/include/

CUDALINK = -lcurand_static -lculibos -lcudart_static -lpthread -ldl -I /usr/local/cuda-12.1/include -L /usr/local/cuda-12.1/lib64

TARGETS = testSIMCoVGPU debugSIMCoVGPU timingSIMCoVGPU commDebug profileSIMCoVGPU
OBJS = kernels.o simcovcuda.o
all: $(TARGETS)
clean:
	rm -f $(TARGETS) $(OBJS)
.PHONY: all clean

testSIMCoVGPU:
	$(NVCC) $(NVCCFLAGS) $(NVCCARCH_FLAGS) $(CXXFLAGS) $(UPCXXINC) -c src/simcovcuda.cu -o kernels.o
	$(UPCXX) $(CXXFLAGS) -c src/simcovcudaMem.cpp -o simcovcuda.o
	$(UPCXX) $(CXXFLAGS) $(CPUSRCS) -o testSIMCoVGPU kernels.o simcovcuda.o $(CUDALINK)

debugSIMCoVGPU:
	$(NVCC) $(NVCCFLAGS) $(NVCCARCH_FLAGS) $(DEBUGFLAGS) $(UPCXXINC) -G -c src/simcovcuda.cu -o kernels.o
	$(UPCXX) $(DEBUGFLAGS) -c src/simcovcudaMem.cpp -o simcovcuda.o
	$(UPCXX) $(DEBUGFLAGS) $(CPUSRCS) -o debugSIMCoVGPU kernels.o simcovcuda.o $(CUDALINK)

timingSIMCoVGPU:
	$(NVCC) -DTIMING $(NVCCFLAGS) $(NVCCARCH_FLAGS) $(CXXFLAGS) $(UPCXXINC) -c src/simcovcuda.cu -o kernels.o
	$(UPCXX) -DTIMING $(CXXFLAGS) -c src/simcovcudaMem.cpp -o simcovcuda.o
	$(UPCXX) -DTIMING $(CXXFLAGS) $(CPUSRCS) -o timingSIMCoVGPU kernels.o simcovcuda.o $(CUDALINK)

commDebug:
	$(NVCC) $(NVCCFLAGS) $(NVCCARCH_FLAGS) $(DEBUGFLAGS) -DCOMMDEBUG $(UPCXXINC) -G -c src/simcovcuda.cu -o kernels.o
	$(UPCXX) $(DEBUGFLAGS) -DCOMMDEBUG -c src/simcovcudaMem.cpp -o simcovcuda.o
	$(UPCXX) $(DEBUGFLAGS) -DCOMMDEBUG $(CPUSRCS) -o commDebug kernels.o simcovcuda.o $(CUDALINK)

profileSIMCoVGPU:
	$(NVCC) -DPROFILE $(NVCCFLAGS) $(NVCCARCH_FLAGS) $(CXXFLAGS) $(UPCXXINC) -c src/simcovcuda.cu -o kernels.o
	$(UPCXX) -DPROFILE $(CXXFLAGS) -c src/simcovcudaMem.cpp -o simcovcuda.o
	$(UPCXX) -DPROFILE $(CXXFLAGS) $(CPUSRCS) -o timingSIMCoVGPU kernels.o simcovcuda.o $(CUDALINK)

rebuild: clean testSIMCoVGPU