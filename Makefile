COMPILER=GPP
include make.inc

all : ising.exe

ising.exe : SRC/*.cpp
	$(CPP) $(DEBUG) $(OPTIMIZE) $(VERSION)17 -o $@ $^ $(BLAS) $(LAPACK) $(ARPACK)

clean :
	rm -rf *.exe *.o *~ *.mod
