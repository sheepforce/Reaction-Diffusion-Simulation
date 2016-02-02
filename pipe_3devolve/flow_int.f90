! FLOW_INT will calculate the changes in conentrations of each substance due to flow of the liquid.
! The concentration of a cell changes because of outflow and inflow. Outflow needs the concentration
! of the current cell, inflow the concentration of current cell - 1 and therefore FlowConcMat stores in
! the first index the current coentrations (2) and the conentrations of the precursor cell (1) and in the
! second index the concentrations of each Substance. 
! 
! PipeLength	--> length of the pipe as in the calling programm
!
! vFlowXY	--> flow speed at current point
!
! FlowConcMat	--> Matrix of the Size (2,Omega), containinig concentrations of substances at the current (2)
! 		and previous point (1)
!
! dTime		--> Integrations step size
!
! DeltaConc	--> Output of the CHANGES in concentration in the current cell to the calling programm
!
! DeltaN	--> displacement of the concentration. Is the number of cells the concentrations should be
! 		displaced and should allway be smaller than 1
!
! DeltaL	--> displacement of the concentration. Is the real distance concentration would be displaced
!

SUBROUTINE FLOW_INT(PipeLength, NGridZ, vFlowXY, FlowConcMat, dTime, DeltaConc)

IMPLICIT NONE
REAL*8, INTENT(IN) :: PipeLength, vFlowXY, dtime
REAL*8, DIMENSION(:,:), INTENT(IN) :: FlowConcMat
INTEGER, INTENT(IN) :: NGridZ

REAL*8, DIMENSION(SIZE(FlowConcMat,2)), INTENT(OUT) :: DeltaConc

REAL*8 :: DeltaL, DeltaN
INTEGER :: Omega, I, lambda


!!!!!!!!!!!!!!!!!!!!!!!!!!
!! INITIALIZATION PHASE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!

Omega = SIZE(FlowConcMat, 2)
lambda = 1


!!!!!!!!!!!!!!!!!!!!!!!
!! CALCULATION PHASE !!
!!!!!!!!!!!!!!!!!!!!!!!

DeltaL = vFlowXY * dTime										! calculate the displacement in absolut space
DeltaN = DeltaL / (PipeLength / DBLE(NGridZ))								! calculate the displacement in relative grid space

! if the movement would be larger than 1 cell the values are not reliable, therefore kill the 
! programm, raice an exception and 
IF (DeltaN > 1.0d0) THEN
	WRITE(*, *) "ERROR : Integration step too large. You must choose smaller integration steps."
	WRITE(*, *) "ERROR : Movement im Grid points would be", DeltaN
	WRITE(*, *) "ERROR : Will terminate now, look over your input again!"
	STOP 221
END IF

! decrease in concentration due to outflow from current cell and increase due to inflow from previous cell
DO lambda = 1, Omega
	DeltaConc(lambda) = - DeltaN * FlowConcMat(2,lambda) + DeltaN * FlowConcMat(1,lambda)
END DO


END SUBROUTINE
