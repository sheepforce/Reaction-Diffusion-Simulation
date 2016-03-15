
! Subroutine FLOW_INT calculates changes in concentration due to laminar flow.
! Takes local flow speed, vLoc, time step size, dTime, the Length of the pipe,
! PipeLenght, the number of Grid points in z direction, NGridZ, the number of 
! sunstances in the system, Omega.
! LocalConc stores prevoius (1) and current point (2) in first index, substances
! in the second index
! returns the CHANGES in concentration for every substances at the current point
! in DeltaConc

SUBROUTINE FLOW_INT(LocalConc, vLoc, dTime, PipeLength, NGridZ, Omega, DeltaConc)
IMPLICIT NONE
REAL*8, DIMENSION(1:2,1:Omega), INTENT(IN) :: LocalConc
REAL*8, INTENT(IN) :: vLoc, dTime, PipeLength
INTEGER, INTENT(IN) :: NGridZ, Omega
REAL*8, DIMENSION(1:Omega), INTENT(OUT) :: DeltaConc

INTEGER :: lambda
REAL*8 :: AbsMove, GridMove

!! DEBUG START
!WRITE(*, *) "vLoc is", vLoc
!! DEBUG END

!!!!!!!!!!!!!!!!!!!!
!! INITIALIZATION !!
!!!!!!!!!!!!!!!!!!!!

AbsMove = vLoc * dTime											! movement of liquid in real units
GridMove = AbsMove / (PipeLength / DBLE(NGridZ))							! movement of liquid in grid points

IF (GridMove > 1.0d0) THEN
	WRITE(*, *) "ERROR : too large movement due to laminar flow, would be", GridMove, "points"
	WRITE(*, *) "ERROR : would move concentrations more than 1 grid point. Not reliable"
	WRITE(*, *) "ERROR : redefine flow speed or integration step size in your input"
	STOP 221
END IF

!!!!!!!!!!!!!!!!!
!! CALCULATION !!
!!!!!!!!!!!!!!!!!

DO lambda = 1, Omega
	DeltaConc(lambda) = +LocalConc(1,lambda) * GridMove - LocalConc(2,lambda) * GridMove		! inflow from previous point, outflow from current point
END DO

END SUBROUTINE
