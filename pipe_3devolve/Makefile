# Redefine this variables for custome build

# Fortran Compiler, possible options:
FC		= mpif90
#FC		= gfortran

# Compilerflags for fortran compiler
FFLAGS		= -fopenmp -cpp -O3 -march=native -mtune=native -funroll-loops -ffast-math -mfpmath=sse

# Preprocessor variables
# use MPI parallelitzation? Alternatively openMP will be used
#FPPOPTIONS	= -DopenMPI		# use MPI parallelization

all: reform_conc test_pipe_3devolve pipe_3devolve rks_solv test_rks_solv

test_pipe_3devolve: routines pipe_3devolve
	$(FC) $(FFLAGS) -c test_pipe_3devolve.f90
	$(FC) $(FFLAGS) test_pipe_3devolve.o flow_profile.o find_pipe.o pipe_3devolve.o rks_int.o rks_rhs_ode.o flow_int.o laplacian.o diff_int.o rks_rk4int.o -o test_pipe_3devolve.bin

test_rks_solv: rks_solv
	$(FC) $(FFLAGS) $(FPPOPTIONS) -c test_rks_solv.f90
	$(FC) $(FFLAGS) $(FPPOPTIONS) test_rks_solv.o rks_solv.o rks_int.o rks_rhs_ode.o rks_rk4int.o -o test_rks_solv.bin

pipe_3devolve: routines
	$(FC) $(FFLAGS) $(FPPOPTIONS) -c pipe_3devolve.f90
	$(FC) $(FFLAGS) $(FPPOPTIONS) -c init_pipe_3devolve.f90
	$(FC) $(FFLAGS) $(FPPOPTIONS) flow_profile.o find_pipe.o pipe_3devolve.o rks_int.o rks_rhs_ode.o flow_int.o laplacian.o diff_int.o rks_rk4int.o init_pipe_3devolve.o -o pipe_3devolve.bin

rks_solv: routines
	$(FC) $(FFLAGS) $(FPPOPTIONS) -c rks_solv.f90

routines:
	$(FC) $(FFLAGS) -c flow_profile.f90 find_pipe.f90 rks_int.f90 rks_rhs_ode.f90 flow_int.f90 laplacian.f90 diff_int.f90 rks_rk4int.f90

reform_conc:
	$(FC) $(FFLAGS) reform_concH.f90 -o reform_concH.bin
	$(FC) $(FFLAGS) reform_concP.f90 -o reform_concP.bin
	$(FC) $(FFLAGS) reform_concM.f90 -o reform_concM.bin

clean:
	rm *.o
	rm *.bin
	rm *.dat
	rm -rf out
