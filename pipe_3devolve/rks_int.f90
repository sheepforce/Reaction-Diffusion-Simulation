! RKS_INT will calculate the conecentration change in a finite time interval with the
! reaction kinetics for 1 time step at 1 grid point
!
! EdMat         --> 2-dimensional array containing stoichiometric coefficients of the
!               educts. first index are reactions, second substances
!
! ProdMat       --> 2-dimensional array containing stoichiometric coefficients of the
!               products. first index are reactions, second substances
!
! RateVec       --> 1-dimensional array containing rate coefficients for every reaction
!
! ConcVec       --> 1-dimensional array containing start concentrations for every substance
!		subset of PipeConc -> PipeConc(PointX,PointY,PointZ,:,1)
!
! Omega         --> number of substances, width of matrices EdMat, ProdMat and vector ConcVec
!
! N             --> number of reactions, height of matrices EdMat, ProdMat and vector RateVec
!
! I             --> serial number of iterators over reactions [1 ... N]
!
! lambda        --> serial number of iterators over substances [1 ... Omega]
!
! ODEVec        --> matrix storing old and new values of the right hand side of the ODE
!               2 x Omega matrix, column 1 old ODE, column 2 new ODE
!
! ConcMat       --> matrix storing old and new values of the concentration of the substances
!               2 x Omega matrix, columnt 1 old concentration, column 2 new concentration
!
! SubsNum       --> number of substance to calculate differential equation from.
!               also from [1 ... Omega]
!
! starttime     --> time at which to start the integration step
!
! dTime         --> integration step length
!
!
! method        --> integer determining which integrator to use
!
! 

SUBROUTINE RKS_INT(EdMat, ProdMat, RateVec, ConcVec, dTime, method, DeltaConc)

IMPLICIT NONE
INTEGER, DIMENSION(:,:), INTENT(IN) :: EdMat, ProdMat
REAL*8, DIMENSION(:), INTENT(IN) :: RateVec , ConcVec
REAL*8, INTENT(IN) :: dTime
INTEGER, INTENT(IN) :: method

REAL*8, DIMENSION(:), INTENT(OUT) :: DeltaConc

INTEGER :: Omega, N, I, lambda, SubsNum
REAL*8, DIMENSION(:), ALLOCATABLE :: ODEVec
REAL*8 :: RHS_ODE

INTERFACE
        FUNCTION RHS_ODE(EdMat, ProdMat, RateVec, ConcVec, SubsNum)
		IMPLICIT NONE

		INTEGER, DIMENSION(:,:), INTENT(IN) :: EdMat, ProdMat
		INTEGER :: Omega, N, I, lambda, SubsNum, NSteps, IntStep
		REAL*8, DIMENSION(:), INTENT(IN) :: RateVec , ConcVec
		REAL*8, DIMENSION(1:SIZE(RateVec)) :: ConcProd
		REAL*8 :: SIGMAProd, SIGMACons
        END FUNCTION RHS_ODE
END INTERFACE


!!!!!!!!!!!!!!!!!!!!!!!!!!
!! INITIALIZATION PHASE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!

Omega = SIZE(ConcVec)                                                                                   ! determines number of substances
N = SIZE(RateVec)                                                                                       ! determines number of reactions
lambda = 1
I = 1

ALLOCATE(ODEVec(1:Omega))                                                                               ! allocate vector, that stores values of the ODEs


!!!!!!!!!!!!!!!!!
!! INTEGRATION !!
!!!!!!!!!!!!!!!!!

! calculate the values of the ODE right hand side with start concentrations
DO lambda = 1, Omega
        ODEVec(lambda) = RHS_ODE(EdMat, ProdMat, RateVec, ConcVec, lambda) 
END DO

! calculate the __CHANGE__ in concentration of every substance
DO lambda = 1, Omega
	DeltaConc(lambda) = ODEVec(lambda) * dtime
END DO

DEALLOCATE(ODEVec)

END SUBROUTINE
