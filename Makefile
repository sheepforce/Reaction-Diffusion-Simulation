# This is the primary Makefile for the reaction diffusion flow simulator.
# It supports building the programm and installation in different configurations.
# Read the following lines and change to your needs.
#
############################
## CONFIGURATION OPTIIONS ##
############################
#
# Installation prefix
PREFIX		= ~/software

# Fortran compiler, must support Fortran90, others will propably work too but were not tested
FC		= gfortran			# gfortran from the gnu compiler collection, RECOMMENDED
#FC		= mpif90			# OpenMPI mpif90 if you want to build with MPI Support, EXPERIMENTAL

# Fortran compiler flags
FFLAGS		=  -fopenmp			# enables OpenMP SMP parallelization, RECOMMENDED
FFLAGS		+= -cpp				# use C preprocessor directives, NECESSARY
FFLAGS		+= -O3				# highest optimization level. -O2 is also fine, RECOMMENDED
FFLAGS		+= -march=native		# optimize for cpu architecture, OPTIONAL
FFLAGS		+= -mtune=native		# optimize for cpu architecture, OPTIONAL
FFLAGS		+= -ffast-math			# optimize math, OPTIONAL
FFLAGS		+= -mfpmath=sse			# optimize math, OPTIONAL

# Parallelization
#FPPOPTIONS	= -DopenMPI			# preprocessor directives, enables MPI, disbales OpenMP, EXPERIMENTAL

##########################################
## NO CHANGES SHOULD BE NECESSARY BELOW ##
##########################################

NUMERIC_SOURCE	=  flow_profile.f90
NUMERIC_SOURCE	+= find_pipe.f90
NUMERIC_SOURCE	+= rks_int.f90
NUMERIC_SOURCE	+= rks_rhs_ode.f90
NUMERIC_SOURCE	+= flow_int.f90
NUMERIC_SOURCE	+= laplacian.f90
NUMERIC_SOURCE	+= diff_int.f90
NUMERIC_SOURCE	+= rks_rk4int.f90

NUMERIC_OBJECT	=  flow_profile.o
NUMERIC_OBJECT	+= find_pipe.o
NUMERIC_OBJECT	+= rks_int.o
NUMERIC_OBJECT	+= rks_rhs_ode.o
NUMERIC_OBJECT	+= flow_int.o
NUMERIC_OBJECT	+= laplacian.o
NUMERIC_OBJECT	+= diff_int.o
NUMERIC_OBJECT	+= rks_rk4int.o

all: pipe_3devolve reform_conc

pipe_3devolve: numeric_routines
	cd pipe_3devolve && $(FC) $(FFLAGS) $(FPPOPTIONS) -c pipe_3devolve.f90
	cd pipe_3devolve && $(FC) $(FFLAGS) $(FPPOPTIONS) -c init_pipe_3devolve.f90
	cd pipe_3devolve && $(FC) $(FFLAGS) $(FPPOPTIONS) $(NUMERIC_OBJECT) pipe_3devolve.o init_pipe_3devolve.o -o pipe_3devolve.bin

reform_conc:
	cd pipe_3devolve && for i in reform_conc*.f90 ; do $(FC) $(FFLAGS) $$i -o $${i%".f90"}.bin ; done

numeric_routines:
	cd pipe_3devolve && $(FC) $(FFLAGS) -c $(NUMERIC_SOURCE)

testing: pipe_3devolve
	cd pipe_3devolve && $(FC) $(FFLAGS) -c test_pipe_3devolve.f90
	cd pipe_3devolve && $(FC) $(FFLAGS) $(NUMERIC_OBJECT) test_pipe_3devolve.o pipe_3devolve.o -o test_pipe_3devolve.bin

clean:
	cd pipe_3devolve && rm *.o
	cd pipe_3devolve && rm *.bin

install:
	mkdir -p $(PREFIX)/bin
	cd pipe_3devolve && for i in *.bin ; do cp $$i $(PREFIX)/bin/$${i%".bin"} ; done
	cp InputParser/inputparser.py $(PREFIX)/bin/inputparser
	cp ConcFileCreator/concfilecreator $(PREFIX)/bin/concfilecreator
	cp Plotter/Plotter.py $(PREFIX)/bin/plotter
	cp -r InputParser/Data $(PREFIX)/bin/.
	cp -r ConcFileCreator/Data $(PREFIX)/bin/.
	cp -r Plotter/Data $(PREFIX)/bin/.
	cd $(PREFIX)/bin && chmod +x *
