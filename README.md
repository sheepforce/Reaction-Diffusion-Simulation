# Reaction-Diffusion-Simulation
Simulating a chemical system in a 3-dimensional circular pipe. Laminar flow, reaction kinetics and diffusion are taken into account. The numerical routines are implemented in Fortran90. For preparation of the input a Python program is used, which offers a graphical interface. Output can be visualized by Gnuplot, which requires awk and sed to process the output.

Parallelization is rudimentary implemented with OpenMP and optionally and experimentally using MPI.

## Dependencies
* gfortran
* make
* Python2
* Qt 4 (optional but higly recommended)
* Gnuplot (optional but highly recommended)
	* awk
	* sed
* OpenMPI (optional)

On Debian you can install those packages with:
> \# aptitude install gfortran make python libpython-dev qt4-default libqt4-dev gnuplot gawk sed openmpi-bin libopenmpi-dev

## Build
A Makefile for use with gfortran on Unix-Systems is provided. Its head looks like this. 

    ############################
    ## CONFIGURATION OPTIIONS ##
    ############################
    #
    # Installation prefix
    PREFIX		= /usr/local
    
    # Fortran compiler, must support Fortran90, others will propably work too but were not tested
    FC			= gfortran			# gfortran from the gnu compiler collection, RECOMMENDED
    #FC			= mpif90			# OpenMPI mpif90 if you want to build with MPI Support, EXPERIMENTAL
    
    # Fortran compiler flags
    FFLAGS		=  -fopenmp			# enables OpenMP SMP parallelization, RECOMMENDED
    FFLAGS		+= -cpp				# use C preprocessor directives, NECESSARY
    FFLAGS		+= -O3				# highest optimization level. -O2 is also fine, RECOMMENDED
    FFLAGS		+= -march=native		# optimize for cpu architecture, OPTIONAL
    FFLAGS		+= -mtune=native		# optimize for cpu architecture, OPTIONAL
    FFLAGS		+= -ffast-math			# optimize math, OPTIONAL
    FFLAGS		+= -mfpmath=sse			# optimize math, OPTIONAL
    
    # Parallelization
    #FPPOPTIONS		= -DopenMPI			# preprocessor directives, enables MPI, disbales OpenMP, EXPERIMENTAL

Edit the prefix if you want a different path for the installation. Comment out

    FC		= gfortran

and uncomment 

    FC		= mpif90
    FPPOPTIONS 	= -DopenMPI

if you would like  to use OpenMPI instead of OpenMP.

You can now enter
> make all && make install

# How to prepare an input
