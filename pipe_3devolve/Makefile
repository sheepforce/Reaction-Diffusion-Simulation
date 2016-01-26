FC		= mpif90
FFLAGS		= -O3

test_pipe_3devolve: flow_profile
	$(FC) -c test_pipe_3devolve.f90
	$(FC) test_pipe_3devolve.o flow_profile.o -o test_pipe_3devolve.bin

pipe_3devolve: flow_profile
	$(FC) -c pipe_3devolve.f90

flow_profile:
	$(FC) -c flow_profile.f90

clean:
	rm *.o
	rm *.bin
	rm *.out
