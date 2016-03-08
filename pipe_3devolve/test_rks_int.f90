! Calls RKS_INT to test its functionality

PROGRAM TEST_RKS_INT
IMPLICIT NONE

INTEGER, DIMENSION(2,4) :: EdMat, ProdMat
REAL*8, DIMENSION(4) :: ConcVec
REAL*8, DIMENSION(2) :: RateVec
REAL*8 :: dtime = 1.0d-2, finaltime=1.0d3
INTEGER :: method = 2, i, lambda, Omega, NSteps, N, IntStep
REAL*8, DIMENSION(4) :: DeltaConc

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

WRITE(*, *) "================================================================================"
WRITE(*, *) "|| Here is TEST_RKS_INT. Testing the routine RKS_INT with some initial values ||"
WRITE(*, *) "================================================================================"

NSteps = NINT(FinalTime / dTime)

EdMat(1,1) = 1
EdMat(1,2) = 1
EdMat(1,3) = 0
EdMat(1,4) = 0
EdMat(2,1) = 0
EdMat(2,2) = 0
EdMat(2,3) = 1
EdMat(2,4) = 1

WRITE(*, *)
WRITE(*, *) "EdMat"
DO i = 1, 2
        WRITE(*, *) "        ", EdMat(i,:)
END DO
WRITE(*, *)



ProdMat(1,1) = 0
ProdMat(1,2) = 0
ProdMat(1,3) = 1
ProdMat(1,4) = 1
ProdMat(2,1) = 1
ProdMat(2,2) = 1
ProdMat(2,3) = 0
ProdMat(2,4) = 0

WRITE(*, *)
WRITE(*, *) "ProdMat"
DO i = 1, 2
        WRITE(*, *) "        ", ProdMat(i,:)
END DO
WRITE(*, *)

ConcVec(1) = 10.0d0
ConcVec(2) = 4.0d0
ConcVec(3) = 2.0d0
ConcVec(4) = 0.2d0

WRITE(*, *)
WRITE(*, *) "ConcVec"
WRITE(*, *) "        ", ConcVec(1:)


RateVec(1) = 0.005d0
RateVec(2) = 0.01d0

WRITE(*, *)
WRITE(*, *) "RateVec"
DO i = 1, 2
        WRITE(*, *) "        ", RateVec(i)
END DO
WRITE(*, *)

WRITE(*, *) "*****************************"
WRITE(*, *) "**** now calling RKS_INT ****"
WRITE(*, *) "*****************************"

Omega = SIZE(ConcVec)                                                                                   ! determines number of subst$
N = SIZE(RateVec)
lambda = 1
I = 1

OPEN(UNIT=100, FILE="Concentration.dat")

WRITE(UNIT=100, FMT='(A16, A4, $)') "time", "    "
DO lambda = 1, Omega
        WRITE(UNIT=100, FMT='(A10, I6, A4, $)') "Substance ", lambda, "    "
END DO
WRITE(UNIT=100, FMT='(A1)') " "

DO IntStep = 1, NSteps

	CALL RKS_INT(EdMat, ProdMat, RateVec, ConcVec, dtime, method, DeltaConc)
	
	DO lambda = 1, Omega
		ConcVec(lambda) = DeltaConc(lambda) + ConcVec(lambda)
	END DO

!	WRITE(*, *) "DeltaConc(1) is", DeltaConc(1)

        WRITE(UNIT=100, FMT='(D16.8, A4, $)') IntStep * dtime, "    "                                           ! writing the time t$
        DO lambda = 1, Omega
                WRITE(UNIT=100, FMT='(D16.8, A4, $)') ConcVec(lambda), "    "
        END DO
        WRITE(UNIT=100, FMT='(A1)') " "

END DO

WRITE(*, *) "*****************************"
WRITE(*, *) "*** finished with RKS_INT ***"
WRITE(*, *) "*****************************"

CLOSE(UNIT=100)

END PROGRAM
