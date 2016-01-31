! RHS_ODE calculates numerical value of the right hand side of the ODE

REAL*8 FUNCTION RHS_ODE(EdMat, ProdMat, RateVec, ConcVec, SubsNum)
IMPLICIT NONE

INTEGER, DIMENSION(:,:), INTENT(IN) :: EdMat, ProdMat
INTEGER :: Omega, N, I, lambda, SubsNum, NSteps, IntStep
REAL*8, DIMENSION(:), INTENT(IN) :: RateVec , ConcVec

REAL*8, DIMENSION(1:SIZE(RateVec)) :: ConcProd
REAL*8 :: SIGMAProd, SIGMACons	


N = SIZE(RateVec)
Omega = SIZE(ConcVec)

ConcProd = 1.0d0

! calculate products of concentrations^stochiometric factor for every reaction
DO I = 1, N
	DO lambda = 1, Omega
		ConcProd(I) = ConcVec(lambda) ** EdMat(I,lambda) * ConcProd(I)
	END DO 
END DO 


SIGMAProd = 0.0d0
SIGMACons = 0.0d0

! calculate the sums in the ODE

DO I = 1, N
	SIGMAProd = SIGMAProd + ( RateVec(I) * ProdMat(I,SubsNum) * ConcProd(I) )
	SIGMACons = SIGMACons + ( RateVec(I) * EdMat(I,SubsNum) * ConcProd(I) )
END DO

RHS_ODE = SIGMAProd - SIGMACons

END FUNCTION RHS_ODE
