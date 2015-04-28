

CC = gcc
LAPACK = /Users/mmccull/lib
OPTS = -O3 -ftree-vectorize

einstein: compute_diffusion_constant_einstein.c dcdlib.c dcdlib.h 
	$(CC) -c compute_diffusion_constant_einstein.c dcdlib.c $(OPTS)
	$(CC) compute_diffusion_constant_einstein.o dcdlib.o -o compute_diffusion_constant_einstein.x $(OPTS)
#	cp compute_diffusion_constant_einstein.x ../../tip3_md/nvt_run/

