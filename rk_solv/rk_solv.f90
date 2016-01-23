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

INCLUDE 'rk_solv.fh'											! include type declarations

REAL*8 ::  DiffEquation
INTERFACE
        FUNCTION DiffEquation(EdMat, ProdMat, RateVec, ConcVec, SubsNum)
                IMPLICIT NONE

                INCLUDE 'rk_solv.fh'

                REAL*8 :: ConcMultiUse, ConcMultiProduce, SumUse, SumProduce                           ! interim results of ODE
        END FUNCTION DiffEquation
END INTERFACE


!!=============!!
!! START DEBUG !!
!!=============!!
WRITE(*, *) "======================================================="
WRITE(*, *) "|| Here is RK_SOLV. Numerical solving the ODE system ||"
WRITE(*, *) "======================================================="
WRITE(*, *)
WRITE(*, *)
WRITE(*, *) "Punching some values for debug purposes now"
WRITE(*, *) "-------------------------------------------"
WRITE(* ,*) 
WRITE(*, *) "== 0-dimensional =="
WRITE(*, *) "NUMTHREAD     ", NUMTHREAD
WRITE(*, *) "method        ", method
WRITE(*, *) "finaltime     ", finaltime
WRITE(*, *) "dtime         ", dtime
WRITE(*, *)
WRITE(*, *) "== 1-dimensional =="
WRITE(*, *) "ConcVec with size", SIZE(ConcVec)
WRITE(*, *) "        ", ConcVec(1:)
WRITE(*, *)
WRITE(*, *) "RateVec with size", SIZE(RateVec)
DO i = 1, SIZE(RateVec)
        WRITE(*, *) "        ", RateVec(i)
END DO
WRITE(*, *)
WRITE(*, *) "== 2-dimensional =="
WRITE(*, *) "ProdMat with size", SIZE(ProdMat, 1), "x", SIZE(ProdMat, 2)
DO i = 1, SIZE(ProdMat, 1)
        WRITE(*, *) "        ", ProdMat(i,:)
END DO
WRITE(*, *)
WRITE(*, *) "EdMat with size", SIZE(EdMat, 1), "x", SIZE(EdMat, 2)
DO i = 1, SIZE(EdMat, 1)
        WRITE(*, *) "        ", EdMat(i,:)
END DO

DO i = 1, 5
	WRITE(*, *)
END DO
!!=============!!
!!  END DEBUG  !!
!!=============!!

!!!!!!!!!!!!!!!!!!!!!!!!!!
!! INITIALIZATION PHASE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!

Omega = SIZE(ConcVec)                                                                                   ! determines number of substances
N = SIZE(RateVec)                                                                                       ! determines number of reactions
NSteps = NINT(finaltime / dtime)                                                                        ! Number of steps, calculated by time to analyse and time step length

!!=============!!
!! START DEBUG !!
!!=============!!
WRITE(*, *) "Omega is", Omega
WRITE(*, *) "N is", N

DO i = 1, 5
        WRITE(*, *)
END DO
!!=============!!
!!  END DEBUG  !!
!!=============!!
ALLOCATE(IntConcMat(1:Omega,1:2))                                                                       ! allocate matrix, that stores input and ouput of integrator step


!!!!!!!!!!!!!!!!!
!! INTEGRATION !!
!!!!!!!!!!!!!!!!!

OPEN(UNIT=100, FILE="Concentration.dat")                                                                ! opens output file

!!=============!!
!! START DEBUG !!
!!=============!!
WRITE(*, *) "parallel integration routine is launchin now"
WRITE(*, *) "--------------------------------------------"
WRITE(*, *)
WRITE(*, *)
!!=============!!
!!  END DEBUG  !!
!!=============!!

! initialize IntConcMat for intial conditions at t = 0
!!!$OMP PARALLEL NUM_THREADS(NUMTHREAD)

!!=============!!
!! START DEBUG !!
!!=============!!
WRITE(*, *) "intial ODE calculation, calling function DiffEquation, OpenMP thread", omp_get_thread_num()
!!=============!!
!!  END DEBUG  !!
!!=============!!
DO lambda = 1, Omega
        IntConcMat(lambda,1) = DiffEquation(EdMat, ProdMat, RateVec, ConcVec, lambda)                   ! Calculate Value of ODE at initial conditions, t = 0
        IntConcMat(lambda,2) = 0
END DO
!!!$OMP END PARALLEL

!!=============!!
!! START DEBUG !!
!!=============!!
WRITE(*, *)
WRITE(*, *) "Initial ODE values stored in IntConcMat(1,:). IntConcMat(2,:) initialized with 0. Here is IntConcMat"
DO i = 1, SIZE(IntConcMat, 1)
	WRITE(*, *) IntConcMat(i,1:)
END DO
!!=============!!
!!  END DEBUG  !!
!!=============!!

!DO IntStep = 1, NSteps                                                                                  ! Integrate NSteps times
!        !$OMP PARALLEL NUM_THREADS(NUMTHREAD)
!        DO lambda = 1, Omega                                                                            ! Iterate over every substance before doing next integration step
!                IF (method == 1) THEN                                                                   ! Euler Integrator
!                        IntConcMat(lambda,2) = IntConcMat(lambda,1) - dtime * IntConcMat(lambda,1)
!                END IF
!
!        END DO
!        !$OMP END PARALLEL
!
!        WRITE(UNIT=100, FMT='(D10.5, 4X)') dtime * IntStep, (IntConcMat(lambda,2), lambda = 1, Omega)   ! write result of this round to output file
!	
!	!$OMP PARALLEL NUM_THREADS(NUMTHREAD)
!	DO lambda = 1, Omega
!		IntConcMat(lambda,1) = DiffEquation(EdMat, ProdMat, RateVec, IntConcMat(:,2), lambda)	! calculate new value of differential equation with integatred values
!	END DO
!	!$OMP END PARALLEL
!
!END DO

CLOSE(UNIT=100)

DEALLOCATE(IntConcMat)

END SUBROUTINE RK_SOLV



!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RIGHT HAND SIDE OF ODE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL*8 FUNCTION DiffEquation(EdMat, ProdMat, RateVec, ConcVec, SubsNum)                                 ! Calculate Initial Value of ODE by using 
IMPLICIT NONE

INCLUDE 'rk_solv.fh'

REAL*8 :: ConcMultiConsume, ConcMultiProduce, SumConsume, SumProduce					! interim results of ODE

!!=============!!
!! START DEBUG !!
!!=============!!
WRITE(*, *) "================================"
WRITE(*, *) "|| FUNCTION DiffEquation here ||"
WRITE(*, *) "================================"
WRITE(*, *)
WRITE(*, *) "*** SubsNum", SubsNum, "***"
!!=============!!
!!  END DEBUG  !!
!!=============!!


Omega = SIZE(ConcVec)                                                                                   ! determines number of substances
N = SIZE(RateVec)

ConcMultiConsume = 1.0d0										! neutral values regarding their operation
ConcMultiProduce = 1.0d0
SumConsume = 0.0d0
SumProduce = 0.0d0

! consuming part of ODE
DO I = 1, N
	DO lambda = 1, Omega
		ConcMultiConsume = ( ConcVec(lambda) ** EdMat(I,lambda) ) * ConcMultiConsume		! calculating the product term of consuming part of the ODE
	END DO
	SumConsume = RateVec(I) * EdMat(I,SubsNum) * ConcMultiConsume + SumConsume			! calculating the sum term (over products) of the consuming part of the ODE
	
	ConcMultiConsume = 1.0d0									! reset product term for next loop
END DO
WRITE(*, *) "SumConsume is", SumConsume


! producing part of ODE
DO I = 1, N
	DO lambda = 1, Omega
		ConcMultiProduce = ( ConcVec(lambda) ** ProdMat(I,lambda) ) * ConcMultiProduce		! calculating the product term of producing part of the ODE
	END DO
	SumProduce = RateVec(I) * ProdMat(I,lambda) * ConcMultiProduce + SumProduce			! calculating the sum term (over products) of the producing part of the ODE

	ConcMultiProduce = 1.0d0
END DO
WRITE(*, *) "SumProduce is", SumProduce

DiffEquation = - SumConsume + SumProduce

END FUNCTION DiffEquation
