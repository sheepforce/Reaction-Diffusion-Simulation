! Calculates numerical value for the 3D spatial Laplacian with a 3x3x3 Matrix as input
! Only direct neighbours are considered an therefor a 3D cross is used on the space for
! calculating the laplacian

REAL*8 FUNCTION LAPLACIAN(SpatConcDist)

IMPLICIT NONE
REAL*8, DIMENSION(-1:+1,-1:+1,-1:+1), INTENT(IN) :: SpatConcDist

LAPLACIAN =	(-6) * SpatConcDist(0,0,0) &
		+ SpatConcDist(-1,0,0) + SpatConcDist(+1,0,0) &
		+ SpatConcDist(0,-1,0) + SpatConcDist(0,+1,0) &
		+ SpatConcDist(0,0,-1) + SpatConcDist(0,0,+1)

END FUNCTION
