FC		= mpif90
FFLAGS		= -O3

test_pipe_3devolve: routines pipe_3devolve
	$(FC) -c test_pipe_3devolve.f90
	$(FC) test_pipe_3devolve.o flow_profile.o find_pipe.o pipe_3devolve.o rks_int.o rks_rhs_ode.o -o test_pipe_3devolve.bin

test_rks_int:
	$(FC) -c rks_int.f90 test_rks_int.f90 rks_rhs_ode.f90
	$(FC) rks_int.o  test_rks_int.o rks_rhs_ode.o -o test_rks_int.bin

pipe_3devolve: routines
	$(FC) -c pipe_3devolve.f90

routines:
	$(FC) -c flow_profile.f90 find_pipe.f90 rks_int.f90 rks_rhs_ode.f90

clean:
	rm *.o
	rm *.bin
