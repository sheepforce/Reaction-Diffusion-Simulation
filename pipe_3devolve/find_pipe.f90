! Writing to a logical array of the size 2*NgridXY+1:2*NgridXY+1, that is mapping the xy plane
! of the pipe, if a grid point is inside or outside the pipe
! 
! REQUIRES AN EXPLICIT INTERFACE

SUBROUTINE FIND_PIPE(vmax, NGridXY, PipeArea)

IMPLICIT NONE
REAL*8, INTENT(IN) ::  vmax
INTEGER, INTENT(IN) :: NGridXY
REAL*8, PARAMETER :: vadd = 0.0d0
LOGICAL, DIMENSION(-NGridXY:+NGridXY,-NGridXY:NGridXY), INTENT(OUT) :: PipeArea
INTEGER :: PointX, PointY
REAL*8 :: PipeTest, FLOW_PROFILE


DO PointX = -NGridXY, NGridXY
	DO PointY = -NGridXY, NGridXY
		PipeTest = FLOW_PROFILE(PointX, PointY, vmax, vadd, NGridXY)
		IF (PipeTest < 0.0d0) THEN
			PipeArea(PointX,PointY) = .FALSE.
		ELSE IF (PipeTest >= 0.0d0) THEN
			PipeArea(PointX,PointY) = .TRUE.
		ELSE
			WRITE(*, *) "ERROR: Numerical instability in FIND_PIPE encountered. Exiting now"
			EXIT
		END IF
	END DO
END DO

END SUBROUTINE
