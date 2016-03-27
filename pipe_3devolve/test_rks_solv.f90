! Calls RKS_INT to test its functionality

PROGRAM TEST_RKS_INT
IMPLICIT NONE

INTEGER, DIMENSION(2,4) :: EdMat, ProdMat
REAL*8, DIMENSION(4) :: ConcVec
REAL*8, DIMENSION(2) :: RateVec
REAL*8 :: dtime = 1.0d-1, finaltime=1.0d5
INTEGER :: method = 2, i, lambda, Omega, NSteps, N, IntStep
REAL*8, DIMENSION(4) :: DeltaConc

INTERFACE
	SUBROUTINE RKS_SOLV(EdMat, ProdMat, RateVec, ConcVec, dTime, FinalTime, method)
		IMPLICIT NONE
		INTEGER, DIMENSION(:,:), INTENT(IN) :: EdMat, ProdMat
		REAL*8, DIMENSION(SIZE(EdMat, 2)) :: ConcVec
		REAL*8, DIMENSION(SIZE(EdMat, 1)), INTENT(IN) :: RateVec
		REAL*8, INTENT(IN) :: dtime, FinalTime
		INTEGER, INTENT(IN) :: method
		INTEGER ::  lambda, Omega, NSteps, N, IntStep, I
		REAL*8, DIMENSION(SIZE(EdMat, 2)) :: DeltaConc
	END SUBROUTINE
END INTERFACE

WRITE(*, *) "================================================================================"
WRITE(*, *) "|| Here is TEST_RKS_INT. Testing the routine RKS_SOLV with some initial values ||"
WRITE(*, *) "================================================================================"

EdMat(1,1) = 1
EdMat(1,2) = 1
EdMat(1,3) = 0
EdMat(1,4) = 0
EdMat(2,1) = 0
EdMat(2,2) = 0
EdMat(2,3) = 1
EdMat(2,4) = 1

ProdMat(1,1) = 0
ProdMat(1,2) = 0
ProdMat(1,3) = 1
ProdMat(1,4) = 1
ProdMat(2,1) = 1
ProdMat(2,2) = 1
ProdMat(2,3) = 0
ProdMat(2,4) = 0

ConcVec(1) = 10.0d0
ConcVec(2) = 4.0d0
ConcVec(3) = 2.0d0
ConcVec(4) = 0.2d0

RateVec(1) = 0.005d0
RateVec(2) = 0.01d0

WRITE(*, *) "***************************"
WRITE(*, *) "**** now calling SOLV *****"
WRITE(*, *) "***************************"

CALL RKS_SOLV(EdMat, ProdMat, RateVec, ConcVec, dTime, FinalTime, method)

WRITE(*, *) "*******************************"
WRITE(*, *) "*** finished with RKS_SOLV ***"
WRITE(*, *) "*******************************"

END PROGRAM
