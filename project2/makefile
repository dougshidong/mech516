OUT = project2.exe
OBJ = prec.f90 glob.f90 eval.f90 riemann_solver.f90 sampling.f90 flux.f90
OFF = prec.o   glob.o   eval.o   riemann_solver.o   sampling.o   flux.f90
DBG = -Wall -Wextra -Wconversion -pedantic -fbounds-check -ffpe-trap=zero,overflow,underflow -fbacktrace -pg
#PRC = -fdefault-real-8 -fdefault-integer-8
#PRC = -fdefault-integer-8
compile:
	gfortran $(PRC) $(OBJ) -c -O3
	gfortran euler1d.f90 $(OFF) -o $(OUT)
	gfortran main.f90 $(OFF) -o exact.exe
dbg:
	gfortran $(PRC) $(DBG) $(OBJ) -c
	gfortran euler1d.f90 $(DBG) $(OFF) -o $(OUT)
	gfortran main.f90 $(DBG) $(OFF) -o exact.exe
pro:
	gfortran $(PRC) $(OBJ) -p -c -O1
	gfortran main.f90 -p $(OFF) -o $(OUT)
clean:
	rm -fv *.out *.mod *.MOD $(OUT)
