OUT = p2a.exe
OBJ =
OFF =
DBG = -Wall -Wextra -Wconversion -pedantic -fbounds-check -ffpe-trap=zero,overflow,underflow -fbacktrace -pg
#PRC = -fdefault-real-8 -fdefault-integer-8
#PRC = -fdefault-integer-8
compile:
	gfortran linear_advection.f90 -o $(OUT)
dbg:
	gfortran $(DBG) linear_advection.f90 -o $(OUT)
pro:
	gfortran $(PRC) $(OBJ) -p -c -O1
	gfortran main.f90 -p $(OFF) -o $(OUT)
clean:
	rm -fv *.out *.mod *.MOD $(OUT)
