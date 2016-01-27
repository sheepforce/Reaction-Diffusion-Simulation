FC		= mpif90
FFLAGS		= -O3

test_pipe_3devolve: routines
	$(FC) -c test_pipe_3devolve.f90
	$(FC) test_pipe_3devolve.o flow_profile.o find_pipe.o -o test_pipe_3devolve.bin

pipe_3devolve: routines
	$(FC) -c pipe_3devolve.f90 flow_profile.f90 find_pipe.f90

routines:
	$(FC) -c flow_profile.f90 find_pipe.f90

clean:
	rm *.o
	rm *.bin
