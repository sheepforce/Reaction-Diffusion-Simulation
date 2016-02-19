FC		= mpif90
FFLAGS		= -O3 -march=native -mtune=native -funroll-loops -ffast-math -mfpmath=sse

test_pipe_3devolve: routines pipe_3devolve
	$(FC) $(FFLAGS) -c test_pipe_3devolve.f90
	$(FC) $(FFLAGS) test_pipe_3devolve.o flow_profile.o find_pipe.o pipe_3devolve.o rks_int.o rks_rhs_ode.o flow_int.o laplacian.o diff_int.o -o test_pipe_3devolve.bin

test_rks_int:
	$(FC) $(FFLAGS) -c rks_int.f90 test_rks_int.f90 rks_rhs_ode.f90
	$(FC) $(FFLAGS) rks_int.o test_rks_int.o rks_rhs_ode.o -o test_rks_int.bin

pipe_3devolve: routines
	$(FC) $(FFLAGS) -c pipe_3devolve.f90

routines:
	$(FC) $(FFLAGS) -c flow_profile.f90 find_pipe.f90 rks_int.f90 rks_rhs_ode.f90 flow_int.f90 laplacian.f90 diff_int.f90

clean:
	rm *.o
	rm *.bin
	rm *.dat
	rm -rf out
