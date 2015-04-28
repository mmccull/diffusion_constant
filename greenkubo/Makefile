

CC = gcc
LAPACK = /Users/mmccull/lib
OPTS = -O3 -ftree-vectorize -fopenmp

greenkubo: compute_diffusion_constant_greenkubo.c dcdlib.c dcdlib.h 
	$(CC) -c compute_diffusion_constant_greenkubo.c dcdlib.c $(OPTS)
	$(CC) compute_diffusion_constant_greenkubo.o dcdlib.o -o compute_diffusion_constant_greenkubo.x $(OPTS)
	cp compute_diffusion_constant_greenkubo.x ../../tip3_md/nvt_run/

