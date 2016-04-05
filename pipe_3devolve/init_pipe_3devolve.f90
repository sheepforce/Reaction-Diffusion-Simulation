PROGRAM INIT_PIPE_3DEVOLVE
IMPLICIT NONE


REAL*8 :: PipeLength, PipeRadius, vmax, vadd, dTime, FinalTime
INTEGER :: NGridXY, NGridZ, method
INTEGER :: EdProdDIM1, EdProdDIM2
INTEGER, DIMENSION(:,:), ALLOCATABLE :: EdMat, ProdMat
REAL*8, DIMENSION(:), ALLOCATABLE :: RateVec, DiffVec
REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE :: PipeConc
LOGICAL :: InFlow


INTEGER :: lambda, Omega, N, I


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! READING THE INPUT FILE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

OPEN(UNIT=100, FILE='Init.dat', STATUS='OLD', FORM='formatted')

READ(100, *) PipeLength
READ(100, *) PipeRadius
READ(100, *) NGridXY
READ(100, *) NGridZ
READ(100, *) vmax
READ(100, *) vadd
READ(100, *) dTime
READ(100, *) FinalTime
READ(100, *) EdProdDIM1
READ(100, *) EdProdDIM2
READ(100, *) method
READ(100, *) Inflow


ALLOCATE(EdMat(EdProdDIM1,EdProdDIM2), ProdMat(EdProdDIM1,EdProdDIM2), RateVec(EdProdDIM1), DiffVec(EdProdDIM2))

NGridZ = NINT(PipeLength / (PipeRadius / NGridXY))
ALLOCATE(PipeConc(-NGridXY:+NGridXY,-NGridXY:+NGridXY,NGridZ,EdProdDIM2,2))


READ(100, *)
READ(100, *) EdMat
READ(100, *)
READ(100, *) ProdMat
READ(100, *)
READ(100, *) RateVec
READ(100, *)
READ(100, *) DiffVec
READ(100, *)
READ(100, *) PipeConc(:,:,:,:,1)

PipeConc(:,:,:,:,2) = 0.0d0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CALLING PIPE 3D EVOLVE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL PIPE_3DEVOLVE(PipeLength, PipeRadius, NGridXY, NGridZ, vmax, vadd, dTime, &
FinalTime, EdMat, ProdMat, RateVec, DiffVec, PipeConc, method, Inflow)




END PROGRAM
