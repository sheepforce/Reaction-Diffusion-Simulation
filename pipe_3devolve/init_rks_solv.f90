! Calls RKS_INT to test its functionality

PROGRAM TEST_RKS_INT
IMPLICIT NONE

INTEGER, DIMENSION(:,:), ALLOCATABLE :: EdMat, ProdMat
REAL*8, DIMENSION(:), ALLOCATABLE :: ConcVec, RateVec
REAL*8 :: dtime, FinalTime
INTEGER :: method, IntStep, EdProdDIM1, EdProdDIM2

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

OPEN(UNIT=100, FILE='Init_RKS.dat', STATUS='OLD', FORM='formatted')

READ(100, *) dTime
READ(100, *) FinalTime
READ(100, *) EdProdDIM1
READ(100, *) EdProdDIM2
READ(100, *) method

ALLOCATE(EdMat(EdProdDIM1,EdProdDIM2), ProdMat(EdProdDIM1,EdProdDIM2), RateVec(EdProdDIM1), ConcVec(EdProdDIM2))

READ(100, *)
READ(100, *) EdMat
READ(100, *)
READ(100, *) ProdMat
READ(100, *)
READ(100, *) RateVec
READ(100, *)
READ(100, *) ConcVec


CALL RKS_SOLV(EdMat, ProdMat, RateVec, ConcVec, dTime, FinalTime, method)

END PROGRAM
