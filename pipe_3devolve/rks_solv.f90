SUBROUTINE RKS_SOLV(EdMat, ProdMat, RateVec, ConcVec, dTime, FinalTime, method)
IMPLICIT NONE

INTEGER, DIMENSION(:,:), INTENT(IN) :: EdMat, ProdMat
REAL*8, DIMENSION(SIZE(EdMat, 2)) :: ConcVec
REAL*8, DIMENSION(SIZE(EdMat, 1)), INTENT(IN) :: RateVec
REAL*8, INTENT(IN) :: dtime, FinalTime
INTEGER, INTENT(IN) :: method
INTEGER ::  lambda, Omega, NSteps, N, IntStep, I
REAL*8, DIMENSION(SIZE(EdMat, 2)) :: DeltaConc

INTERFACE
        SUBROUTINE RKS_INT(EdMat, ProdMat, RateVec, ConcVec, dTime, method, DeltaConc)
                IMPLICIT NONE
                INTEGER, DIMENSION(:,:), INTENT(IN) :: EdMat, ProdMat
                REAL*8, DIMENSION(:), INTENT(IN) :: RateVec , ConcVec
                REAL*8, INTENT(IN) :: dTime
                INTEGER, INTENT(IN) :: method

                REAL*8, DIMENSION(1:SIZE(ConcVec)), INTENT(OUT) :: DeltaConc

                INTEGER :: Omega, N, I, lambda, SubsNum
                REAL*8, DIMENSION(:), ALLOCATABLE :: ODEVec
                REAL*8 :: RHS_ODE
        END SUBROUTINE
END INTERFACE


!!!!!!!!!!!!!!!!!!!!!!!!!!
!! INITIALIZE VARIABLES !!
!!!!!!!!!!!!!!!!!!!!!!!!!!

Omega = SIZE(ConcVec)
N = SIZE(RateVec)
NSteps = NINT(FinalTime / dTime)
lambda = 1
N = SIZE(RateVec)
I = 1


!!//////////////!!
!! OUTPUT START !!
!!\\\\\\\\\\\\\\!!
1000 FORMAT (A80)
1001 FORMAT (4X, A)
1002 FORMAT (1X, I2.2, 1X)
1003 FORMAT (1X, D12.6, 1X)

WRITE(*, 1000) "=================================="
WRITE(*, 1000) "||   REACTION KINETICS MODULE   ||"
WRITE(*, 1000) "=================================="
WRITE(*, FMT='(2/)')
WRITE(*, 1001) "Input Parameters"
WRITE(*, 1001) "----------------"
WRITE(*, *)
WRITE(*, 1001) "matrix of stochiometric coefficients of the educts --> EdMat"
DO I = 1, N
	WRITE(*, FMT='(8X)', ADVANCE='NO')
	DO lambda = 1, Omega
		WRITE(*, 1002, ADVANCE='NO') EdMat(I, lambda)
	END DO
	WRITE(*, *)
END DO
WRITE(*, *)
WRITE(*, 1001) "matrix of stochiometric coefficients of the products --> ProdMat"
DO I = 1, N
	WRITE(*, FMT='(8X)', ADVANCE='NO')
	DO lambda = 1, Omega
		WRITE(*, 1002, ADVANCE='NO') ProdMat(I, lambda)
	END DO
	WRITE(*, *)
END DO
WRITE(*, *)
WRITE(*, 1001) "vector of the the rate coefficients of reactions --> RateVec"
DO I = 1, N
	WRITE(*, FMT='(8X)', ADVANCE='NO')
	WRITE(*, 1003) RateVec(I)
END DO
WRITE(*, *)
WRITE(*, 1001) "vector of the initial concentrations of substances --> ConcVec"
DO lambda = 1, Omega
        WRITE(*, FMT='(8X)', ADVANCE='NO')
        WRITE(*, 1003) ConcVec(lambda)
END DO
WRITE(*, FMT='(4X, A)', ADVANCE='NO') "integrator (kin.)  : "
IF (method == 1) THEN
	WRITE(*, FMT='(4X, A)') "Euler"
ELSE IF (method == 2) THEN
	WRITE(*, FMT='(4X, A)') "Runge-Kutta 4th order"
END IF
WRITE(*, FMT='(2/)')
!!//////////////!!
!!  OUTPUT END  !!
!!\\\\\\\\\\\\\\!!


!!!!!!!!!!!!!!!!!
!! CALCULATION !!
!!!!!!!!!!!!!!!!!

!!//////////////!!
!! OUTPUT START !!
!!\\\\\\\\\\\\\\!!
OPEN(UNIT=100, FILE="Concentration.dat")
WRITE(UNIT=100, FMT='(A16, A4, $)') "time", "    "
DO lambda = 1, Omega
        WRITE(UNIT=100, FMT='(A10, I6, A4, $)') "Substance ", lambda, "    "
END DO
WRITE(UNIT=100, FMT='(A1)') " "
!!//////////////!!
!!  OUTPUT END  !!
!!\\\\\\\\\\\\\\!!

DO IntStep = 1, NSteps
        CALL RKS_INT(EdMat, ProdMat, RateVec, ConcVec, dtime, method, DeltaConc)

        DO lambda = 1, Omega
                ConcVec(lambda) = DeltaConc(lambda) + ConcVec(lambda)
        END DO

        WRITE(UNIT=100, FMT='(D16.8, A4, $)') IntStep * dtime, "    "
        DO lambda = 1, Omega
                WRITE(UNIT=100, FMT='(D16.8, A4, $)') ConcVec(lambda), "    "
        END DO
        WRITE(UNIT=100, FMT='(A1)') " "
END DO

!!//////////////!!
!! OUTPUT START !!
!!\\\\\\\\\\\\\\!!
CLOSE(UNIT=100)
WRITE(*, 1001) "Calculation finished. Check Concentration.dat for evolution of concentration"
WRITE(*, *)
WRITE(*, 1001) "Leaving REACTION KINETIC now. Shutting down parallel computing. Bye!"
WRITE(*, 1001) "--------------------------------------------------------------------"
!!//////////////!!
!!  OUTPUT END  !!
!!\\\\\\\\\\\\\\!!
END SUBROUTINE
