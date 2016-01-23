! RK_SOLV will numerical integrate the differential equations describing the kinetics
! of a rection network
!
! EdMat         --> 2-dimensional array containing stoichiometric coefficients of the
!               educts. first index are reactions, second substances
!
! ProdMat       --> 2-dimensional array containing stoichiometric coefficients of the
!               products. first index are reactions, second substances
!
! RateVec       --> 1-dimensional array containing rate coefficients for every reaction
!
! ConcVec       --> 1-dimensional array containing start concentrations for every substance
!
! Omega         --> number of substances, width of matrices EdMat, ProdMat and vector ConcVec
!
! N             --> number of reactions, height of matrices EdMat, ProdMat and vector RateVec
!
! I             --> serial number of iterators over reactions [1 ... N]
!
! lambda        --> serial number of iterators over substances [1 ... Omega]
!
! IntConcMat    --> Omega x 2 matrix. column 1 contains old concentration of substance
!               Lambda in row Lambda and column 2 new concentration of Lambda in Row Lambda
!
! SubsNum       --> number of substance to calculate differential equation from.
!               also from [1 ... Omega]
!
! starttime     --> time at which to start the integration step
!
! dtime         --> integration step length
!
! finaltime     --> time at which to stop integration
!
! NSteps        --> number of integration steps
!
! IntStep       --> serial number of integration [1 ... NSteps]
!
! method        --> integer determining which integrator to use
!
! NUMTHREAD     --> number of parallel OMP threads
!


SUBROUTINE RK_SOLV(EdMat, ProdMat, RateVec, ConcVec, dtime, finaltime, method, NUMTHREAD)

USE omp_lib                                                                                             ! OMP definitions

IMPLICIT NONE

INCLUDE 'reactionkinetics.fh'                                                                           ! include type declarations

REAL*8 ::  DiffEquation
INTERFACE
        FUNCTION DiffEquation(EdMat, ProdMat, RateVec, ConcVec, SubsNum)
                IMPLICIT NONE

                INCLUDE 'reactionkinetics.fh'

                REAL*8 :: ConcMultiEduct, ConcMultiProduct, SumEduct, SumProd                           ! interim results of ODE
        END FUNCTION DiffEquation
END INTERFACE


!!!!!!!!!!!!!!!!!!!!!!!!!!
!! INITIALIZATION PHASE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!

Omega = SIZE(ConcVec)                                                                                   ! determines number of substances
N = SIZE(RateVec)                                                                                       ! determines number of reactions
NSteps = NINT(finaltime / dtime)                                                                        ! Number of steps, calculated by time to analyse and time step length

ALLOCATE(IntConcMat(1:Omega,1:2))                                                                       ! allocate matrix, that stores input and ouput of integrator step


!!!!!!!!!!!!!!!!!
!! INTEGRATION !!
!!!!!!!!!!!!!!!!!

OPEN(UNIT=100, FILE="Concentration.dat")                                                                ! opens output file

! initialize IntConcMat for intial conditions at t = 0
!$OMP PARALLEL NUM_THREADS(NUMTHREAD)
DO lambda = 1, Omega
        IntConcMat(lambda,1) = DiffEquation(EdMat, ProdMat, RateVec, ConcVec, lambda)                   ! Calculate Value of ODE at initial conditions, t = 0
        IntConcMat(lambda,2) = 0
END DO
!$OMP END PARALLEL

DO IntStep = 1, NSteps                                                                                  ! Integrate NSteps times
        !$OMP PARALLEL NUM_THREADS(NUMTHREAD)
        DO lambda = 1, Omega                                                                            ! Iterate over every substance before doing next integration step
                IF (method == 1) THEN                                                                   ! Euler Integrator
                        IntConcMat(lambda,2) = IntConcMat(lambda,1) - dtime * IntConcMat(lambda,1)
                END IF

        END DO
        !$OMP END PARALLEL

        WRITE(UNIT=100, FMT='(D10.5, 4X)') dtime * IntStep, (IntConcMat(lambda,2), lambda = 1, Omega)   ! write result of this round to output file
	
	!$OMP PARALLEL NUM_THREADS(NUMTHREAD)
	DO lambda = 1, Omega
		IntConcMat(lambda,1) = DiffEquation(EdMat, ProdMat, RateVec, IntConcMat(:,2), lambda)	! calculate new value of differential equation with integatred values
	END DO
	!$OMP END PARALLEL

END DO

CLOSE(UNIT=100)

DEALLOCATE(IntConcMat)

END SUBROUTINE RK_SOLV



!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RIGHT HAND SIDE OF ODE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL*8 FUNCTION DiffEquation(EdMat, ProdMat, RateVec, ConcVec, SubsNum)                                 ! Calculate Initial Value of ODE by using 
IMPLICIT NONE

INCLUDE 'reactionkinetics.fh'

REAL*8 :: ConcMultiEduct, ConcMultiProduct, SumEduct, SumProd                                           ! interim results of ODE

Omega = SIZE(ConcVec)                                                                                   ! determines number of substances
N = SIZE(RateVec)                                                                                       ! determines number of reactions

DO I = 1, N
        DO lambda = 1, Omega
                ConcMultiEduct = ConcVec(lambda) ** (EdMat(I,lambda)) * ConcMultiEduct                  ! calculates product of concentrations^(stochiometric coefficient on educt side)
                ConcMultiProduct = ConcVec(lambda) ** (ProdMat(I,lambda)) * ConcMultiProduct            ! calculates product of concentrations^(stochiometric coefficient on product side)
        END DO

        SumEduct = RateVec(I) * EdMat(I,SubsNum) * ConcMultiEduct                                       ! calculates sum of rate constant of reaction I * stochiometric coefficient of 
        SumProd = RateVec(I) * ProdMat(I,SubsNum) * ConcMultiProduct                                    ! substance SubsNum * ConcMulti... --> yealding consuming and producing part of the right side DE
END DO
DiffEquation = - SumEduct + SumProd

END FUNCTION DiffEquation
