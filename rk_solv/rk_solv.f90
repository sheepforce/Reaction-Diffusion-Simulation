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
! ODEVec        --> matrix storing old and new values of the right hand side of the ODE
!               2 x Omega matrix, column 1 old ODE, column 2 new ODE
!
! ConcMat       --> matrix storing old and new values of the concentration of the substances
!               2 x Omega matrix, columnt 1 old concentration, column 2 new concentration
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
! NUMTHREAD     --> number of parallel OMP Threads
!


SUBROUTINE RK_SOLV(EdMat, ProdMat, RateVec, ConcVec, dtime, finaltime, method, NUMTHREAD)

USE omp_lib                                                                                             ! OMP definitions

IMPLICIT NONE

INCLUDE 'rk_solv.fh'											! include type declarations

REAL*8 :: RHS_ODE

INTERFACE
	FUNCTION RHS_ODE(EdMat, ProdMat, RateVec, ConcVec, SubsNum)
		IMPLICIT NONE

		INCLUDE 'rk_solv.fh'

		REAL*8, DIMENSION(1:SIZE(RateVec)) :: PIProd, PICons
		REAL*8 :: SIGMAProd, SIGMACons
	END FUNCTION RHS_ODE
END INTERFACE


!!!!!!!!!!!!!!!!!!!!!!!!
!! TESTING CONDITIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!

IF (method < 1 .OR. method > 1) THEN
	RETURN
END IF


!!!!!!!!!!!!!!!!!!!!!!!!!!
!! INITIALIZATION PHASE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!

Omega = SIZE(ConcVec)                                                                                   ! determines number of substances
N = SIZE(RateVec)                                                                                       ! determines number of reactions
NSteps = NINT(finaltime / dtime)                                                                        ! Number of steps, calculated by time to analyse and time step length
lambda = 1
I = 1

ALLOCATE(ODEVec(1:Omega))										! allocate vector, that stores values of the ODEs
ALLOCATE(ConcMat(1:Omega,1:2))										! allocate matrix, that stroes old and new values of substance concentration


!!!!!!!!!!!!!!!!!
!! INTEGRATION !!
!!!!!!!!!!!!!!!!!

OPEN(UNIT=100, FILE="Concentration.dat")                                                                ! opens output file

WRITE(UNIT=100, FMT='(A16, A4, $)') "time", "    "
DO lambda = 1, Omega
	WRITE(UNIT=100, FMT='(A10, I6, A4, $)') "Substance ", lambda, "    "
END DO
WRITE(UNIT=100, FMT='(A1)') " "

! initialize ODEVec for intial conditions at t = 0
!$OMP PARALLEL NUM_THREADS(NUMTHREAD)
DO lambda = 1, Omega
        ODEVec(lambda) = RHS_ODE(EdMat, ProdMat, RateVec, ConcVec, lambda)				! Calculate Value of ODE at initial conditions, t = 0
END DO
!$OMP END PARALLEL

! initialize ConcMat with initial concentrations at t = 0
DO lambda = 1, Omega
	ConcMat(lambda,1) = ConcVec(lambda)
	ConcMat(lambda,2) = 0
END DO

! integrate the differential equations
DO IntStep = 1, NSteps
	! calculate finite, new concentration with "cold concentration" + "rate of concentration change" * "time step"
	!$OMP PARALLEL NUM_THREADS(NUMTHREAD)
	DO lambda = 1, Omega
		IF (method == 1) THEN										! Euler Method
			ConcMat(lambda,2) = ConcMat(lambda,1) + ODEVec(lambda) * dtime
		END IF
	END DO
	!$OMP END PARALLEL
	
	! calculate new values of the ODEs with new concentrations
	!$OMP PARALLEL NUM_THREADS(NUMTHREAD)
	DO lambda = 1, Omega
		ODEVec(lambda) = RHS_ODE(EdMat, ProdMat, RateVec, ConcMat(:,2), lambda)
	END DO
	!$OMP END PARALLEL

	! Write the output for this integration step
	WRITE(UNIT=100, FMT='(D16.8, A4, $)') IntStep * dtime, "    "						! writing the time to the first column, no linebreak because $
	DO lambda = 1, Omega
		WRITE(UNIT=100, FMT='(D16.8, A4, $)') ConcMat(lambda,2), "    "
	END DO
	WRITE(UNIT=100, FMT='(A1)') " "
!	WRITE(UNIT=100, FMT='(D16.8, A4, $)') ( ConcMat(lambda,2), "    ", lambda = 1, Omega )
!	WRITE(UNIT=100, FMT='(D12.4, 4X)') IntStep * dtime, (ConcMat(lambda,2), lambda = 1, Omega)		! Writing the time in the first column
!	WRITE(UNIT=100, FMT='(D12.4, 4X)') (ConcMat(lambda,2), lambda = 1, Omega)				! Write the concentrations of substance N in column N+1


	! new concentration becomes old concentration
	!$OMP PARALLEL NUM_THREADS(NUMTHREAD)
	DO lambda = 1, Omega
		ConcMat(lambda,1) = ConcMat(lambda,2)
	END DO
	!$OMP END PARALLEL
END DO
 

CLOSE(UNIT=100)

DEALLOCATE(ODEVec)
DEALLOCATE(ConcMat)

END SUBROUTINE RK_SOLV

