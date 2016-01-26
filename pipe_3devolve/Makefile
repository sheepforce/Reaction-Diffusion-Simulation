FC		= mpif90
FFLAGS		= -O3

pipe_3devolve: flow_profile
	$(FC) -c pipe_3devolve.f90

flow_profile:
	$(FC) -c flow_profile.f90

clean:
	rm *.o
	rm *.bin
