ifndef COMPILER
COMPILER=GPP
endif

ifeq ($(COMPILER),GPP)
	CPP=g++
	WARN=-Wall
	VERSION=-std=c++
	DEBUG=-g
	OPTIMIZE=-O3 -mavx2 -funroll-loops
	OMP=-fopenmp
	LAPACK=-llapack
	BLAS=-lblas
	ARPACK=-larpack
	GFORTRAN=-lgfortran
	LM=-lm
	endif

