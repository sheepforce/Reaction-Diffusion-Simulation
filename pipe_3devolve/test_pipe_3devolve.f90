PROGRAM TEST_PIPE_3DEVOLVE
IMPLICIT NONE
INTEGER :: NGridXY, PointX, PointY, a, b
REAL*8 :: vmax, vadd, PipeRadius, FLOW_PROFILE
REAL*8, DIMENSION(:,:), ALLOCATABLE :: FlowMatrix


NGridXY = 5
vmax = 1.0d0
vadd = 0.0d0

ALLOCATE(FlowMatrix(-NGridXY:+NgridXY,-NGridXY:+NGridXY))

DO a = -NgridXY, NGridXY
	DO b = -NgridXY, NGridXY
		WRITE(*, *) "Now on point", a, b
		FlowMatrix(a,b) = FLOW_PROFILE(a, b, vmax, vadd, NGridXY)
	END DO
END DO


WRITE(*, *) "Here is the FlowMatrix"
WRITE(*, *)
DO a = -NgridXY, NGridXY
	WRITE(*, *) FlowMatrix(a,-NgridXY:)
END DO

END PROGRAM
