! RK4_INT implements 4th order Runge-Kutta integrator

REAL*8 FUNCTION RK4_INT(RHSInit, dTime)
IMPLICIT NONE
REAL*8, INTENT(IN) :: RHSInit, dTime
REAL*8 :: k1, k2, k3, k4

k1 = RHSInit
k2 = RHSInit * (dTime / 2.0d0) + (dTime / 2.0d0) * k1
k3 = RHSInit * (dTime / 2.0d0) + (dTime / 2.0d0) * k2
k4 = RHSInit * dTime + dTime * k3

RK4_INT = (dTime / 6.0d0) * (k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4)


END FUNCTION
