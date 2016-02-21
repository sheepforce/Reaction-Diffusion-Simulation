PROGRAM TEST_PIPE_3DEVOLVE
IMPLICIT NONE

REAL*8 :: PipeLength, PipeRadius, vmax, vadd, dTime, FinalTime
REAL*8, DIMENSION(2) :: RateVec
REAL*8, DIMENSION(4) :: DiffVec
REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE :: PipeConc
INTEGER :: NGridXY, method, NGridZ, Omega, N, lambda, I, PointX, PointY, PointZ
INTEGER, DIMENSION(2,4) :: EdMat, ProdMat
LOGICAL :: Inflow = .TRUE.


INTERFACE
	SUBROUTINE PIPE_3DEVOLVE(PipeLength, PipeRadius, NGridXY, NGridZ, vmax, vadd, dTime, &
	FinalTime, EdMat, ProdMat, RateVec, DiffVec, PipeConc, method, Inflow)
		USE mpi

		INCLUDE 'pipe_3devolve.fh'

		INTEGER :: ierr, NProcs, ProcID, IProc, GenericTag = 1		
	END SUBROUTINE
END INTERFACE

WRITE(*, *) "================================="
WRITE(*, *) "|| Here is TEST_PIPE_3D_EVOLVE ||"
WRITE(*, *) "================================="
WRITE(*, *)
WRITE(*, *) "Initializing variables now"

NGridXY = 5
vmax = 1.0d-3
vadd = 1.0d-4
PipeLength = 5.0d-1
PipeRadius = 1.0d-2
dTime = 1.0d-3
FinalTime = 1.0d2
method = 1

NGridZ = NINT(PipeLength / (PipeRadius / NGridXY))

RateVec(1) = 0.005d0
RateVec(2) = 0.01d0

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

Omega = SIZE(EdMat,2)
N = SIZE(RateVec)

ALLOCATE(PipeConc(-NGridXY:+NGridXY,-NGridXY:+NGridXY,NGridZ,Omega,2))

PipeConc(:,:,:,:,2) = 0.0d0							! no time step has been taken yet, so this is 0
PipeConc(:,:,1,1,1) = 10.0d0							! c(A) in the whole XY plane at first grid point in z
PipeConc(:,:,1,2,1) = 4.0d0							! c(B)
PipeConc(:,:,1,3,1) = 2.0d0							! c(C)
PipeConc(:,:,1,4,1) = 0.2d0							! c(C)
PipeConc(:,:,2:,:,1) = 0.0d0							! concentration of all substances from first xy plane away are 0

DiffVec(1) = 5.0d-2
DiffVec(2) = 1.0d-5
DiffVec(3) = 3.0d-3
DiffVec(4) = 5.0d-4

!WRITE(*, *)
!WRITE(*, *) "Here are the start conditions"
!WRITE(*, *) "NGridXY = ", NGridXY
!WRITE(*, *) "vmax = ", vmax
!WRITE(*, *) "vadd = ", vadd
!WRITE(*, *) "PipeLength = ", PipeLength
!WRITE(*, *) "PipeRadius = ", PipeRadius
!WRITE(*, *) "dTime = ", dTime
!WRITE(*, *) "FinalTime = ", FinalTime
!WRITE(*, *) "method = ", method
!WRITE(*, *) "NGridXY = ", NGridXY
!WRITE(*, *) "NgridZ = ", NGridZ
!WRITE(*, *)
!WRITE(*, *)
!WRITE(*, *)
!WRITE(*, *) "EdMat"
!DO i = 1, 2
!        WRITE(*, *) "        ", EdMat(i,:)
!END DO
!WRITE(*, *)
!WRITE(*, *)
!WRITE(*, *) "ProdMat"
!DO i = 1, 2
!        WRITE(*, *) "        ", ProdMat(i,:)
!END DO
!WRITE(*, *)
!WRITE(*, *)
!WRITE(*, *) "RateVec"
!DO i = 1, 2
!        WRITE(*, *) "        ", RateVec(i)
!END DO
!WRITE(*, *)
!WRITE(*, *)
!WRITE(*, *) "PipeConc"
!WRITE(*, *) "    in xy plane at pipe start (PointZ = 0)"
!DO lambda = 1, Omega
!	WRITE(*, *) "    Substance", lambda
!	DO PointX = -NgridXY, NGridXY
!	WRITE(*, *) "    ", PipeConc(PointX,-NGridXY:,1,lambda,1)
!	END DO
!	WRITE(*, *)
!END DO
!WRITE(*, *) "    in xy plane at the pipe end (PointZ = NGridZ)"
!DO lambda = 1, Omega
!        WRITE(*, *) "    Substance", lambda
!        DO PointX = -NgridXY, NGridXY
!        WRITE(*, *) "    ", PipeConc(PointX,-NGridXY:,NGridZ,lambda,1)
!        END DO
!        WRITE(*, *)
!END DO
!WRITE(*, *)
WRITE(* ,*) "*********************************"
WRITE(*, *) "*** now calling PIPE_3DEVOLVE ***"
WRITE(*, *) "*********************************"

CALL PIPE_3DEVOLVE(PipeLength, PipeRadius, NGridXY, NGridZ, vmax, vadd, dTime, &
FinalTime, EdMat, ProdMat, RateVec, DiffVec, PipeConc, method, Inflow)

END PROGRAM
