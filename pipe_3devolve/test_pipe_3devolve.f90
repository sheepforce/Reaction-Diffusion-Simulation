PROGRAM TEST_PIPE_3DEVOLVE
IMPLICIT NONE
INTEGER :: NGridXY, PointX, PointY
REAL*8 :: vmax, vadd, PipeRadius, FLOW_PROFILE
REAL*8, DIMENSION(:,:), ALLOCATABLE :: FlowMatrix
LOGICAL, DIMENSION(:,:), ALLOCATABLE :: PipeArea


INTERFACE
	SUBROUTINE FIND_PIPE(vmax, NGridXY, PipeArea)
		IMPLICIT NONE
		REAL*8, INTENT(IN) ::  vmax
		INTEGER, INTENT(IN) :: NGridXY
		REAL*8, PARAMETER :: vadd = 0.0d0
		LOGICAL, DIMENSION(-NGridXY:+NGridXY,-NGridXY:NGridXY), INTENT(OUT) :: PipeArea
		INTEGER :: PointX, PointY
		REAL*8 :: PipeTest, FLOW_PROFILE
	END SUBROUTINE
END INTERFACE

NGridXY = 5
vmax = 1.0d0
vadd = 0.0d0

ALLOCATE(FlowMatrix(-NGridXY:+NgridXY,-NGridXY:+NGridXY))


!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TESTING FLOW_PROFILE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(*, *) "Testing function FLOW_PROFILE now"
DO PointX = -NgridXY, +NGridXY
	DO PointY = -NgridXY, +NGridXY
		WRITE(*, *) "Now on point", PointX, PointY
		FlowMatrix(PointX,PointY) = FLOW_PROFILE(PointX, PointY, vmax, vadd, NGridXY)
	END DO
END DO

WRITE(*, *) "Here is the FlowMatrix"
WRITE(*, *)
DO PointX = -NgridXY, NGridXY
	WRITE(*, *) FlowMatrix(PointX,-NgridXY:)
END DO

!!!!!!!!!!!!!!!!!!!!!!
!! TESING FIND_PIPE !!
!!!!!!!!!!!!!!!!!!!!!!

WRITE(*, *)
WRITE(*, *)
WRITE(*, *) "Testing subroutine FIND_PIPE now"

ALLOCATE(PipeArea(-NGridXY:+NGridXY,-NGridXY:+NGridXY))

WRITE(*, *) "look at PipeArea, should be random now"

DO PointX = -NGridXY, +NGridXY
	WRITE(*, *) "Row", PointX, ":", PipeArea(PointX,-NgridXY:)
END DO

WRITE(*, *) "Calling FIND_PIPE now, PipeArea should become meaningfull"
CALL FIND_PIPE(vmax, NGridXY, PipeArea)

WRITE(*, *) "Here is PipeArea as retourned by FIND_PIPE"
DO PointX = -NgridXY, NGridXY
	WRITE(*, *) PipeArea(PointX,-NGridXY:) 
END DO

DEALLOCATE(PipeArea)
END PROGRAM
