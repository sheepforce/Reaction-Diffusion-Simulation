! DIFF_INT calculates the changes in concentrations on a given grid point due to diffusion.
! It uses a discret laplacian
!
! ConcEnv	--> Array of size (-1:+1,-1:+1,-1:+1,Omega) that stores the conentration enviroment
!		the current point for all substances
! DiffVec	--> Vector containing the diffussion coefficient of all substances
!
! DiffConc	--> Output that stores the CHANGES in concentrations of every substance
!
! dTime		--> Integration Step Width
!

SUBROUTINE DIFF_INT(ConcEnv, DiffVec, DeltaConc, dTime)

IMPLICIT NONE
REAL*8, DIMENSION(-1:,-1:,-1:,:), INTENT(IN) :: ConcEnv
REAL*8, DIMENSION(:), INTENT(IN) :: DiffVec

REAL*8, DIMENSION(1:SIZE(DiffVec)), INTENT(OUT) :: DeltaConc

REAL*8 :: LAPLACIAN, dTime
INTEGER :: Omega, lambda

INTERFACE
	FUNCTION LAPLACIAN(SpatConcDist)
		IMPLICIT NONE
		REAL*8, DIMENSION(-1:+1,-1:+1,-1:+1), INTENT(IN) :: SpatConcDist
	END FUNCTION
END INTERFACE


!!!!!!!!!!!!!!!!!!!!!!!!!!
!! INITIALIZATION PHASE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!

Omega = SIZE(DiffVec)
lambda = 1
DeltaConc = 0.0d0

!!!!!!!!!!!!!!!!!
!! CALCULATION !!
!!!!!!!!!!!!!!!!!

DO lambda = 1, Omega
	DeltaConc(lambda) = LAPLACIAN(ConcEnv(-1:+1,-1:+1,-1:+1,lambda)) * DiffVec(lambda) * dTime
END DO

END SUBROUTINE
