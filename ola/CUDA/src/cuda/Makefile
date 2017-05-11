SHELL := /bin/bash # Use bash syntax

CC = sm_21

NVCC	= nvcc -arch=$(CC)
SAMPLE4 = nlmGlobal
SAMPLE5 = Shared
SAMPLE6 = SharedWithMatlab
ARGS	= -ptx

all: sample

sample:
	$(NVCC) $(ARGS) $(SAMPLE4).cu
	$(NVCC) $(ARGS) $(SAMPLE5).cu
	$(NVCC) $(ARGS) $(SAMPLE6).cu

clean:
	rm -rf *~ *.ptx
