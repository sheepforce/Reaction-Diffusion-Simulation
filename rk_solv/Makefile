FC              = gfortran
FFLAGS		= -O3

rk_solv:
	$(FC) -fopenmp -c rk_solv.f90

test_rk_solv: rk_solv
	$(FC) -c -fopenmp test_rk_solv.f90
	$(FC) -fopenmp rk_solv.o test_rk_solv.o -o test_rk_solv.bin

clean:
	rm *.o
	rm *.bin
	rm *.dat
	rm *.log
