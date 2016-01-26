! The FLOW_PROFILE function gives the speed of flow on a given point on the XY Grid
!

REAL*8 FUNCTION FLOW_PROFILE(PointX, PointY, vmax, vadd, PipeRadius, NGridXY)

IMPLICIT NONE
INTEGER, INTENT(IN) :: PointX, PointY, NGridXY
REAL*8, INTENT(IN) :: vmax, vadd, PipeRadius


FLOW_PROFILE = - vmax * (PointX ** 2 + PointY ** 2) + vmax * NGridXY ** 2 + vadd

END FUNCTION
