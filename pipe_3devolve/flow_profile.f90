! The FLOW_PROFILE function gives the speed of flow on a given point on the XY Grid
!

REAL*8 FUNCTION FLOW_PROFILE(PointX, PointY, vmax, vadd, NGridXY)

IMPLICIT NONE
INTEGER, INTENT(IN) :: PointX, PointY, NGridXY
REAL*8, INTENT(IN) :: vmax, vadd


FLOW_PROFILE = vmax * (1 - ((SQRT(DBLE(PointX ** 2 + PointY ** 2))) / DBLE(NGridXY)) ** 2) + vadd

END FUNCTION
