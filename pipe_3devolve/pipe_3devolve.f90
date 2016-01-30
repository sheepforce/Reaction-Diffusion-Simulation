! Subroutine PIPE_3DEVOLVE will numerical solve the time evolution of 
! the concentrations of all substances taking reaction kinetics, 
! diffusion and flow into account
!
!*EdMat         --> 2-dimensional array containing stoichiometric coefficients of the
!               educts. first index are reactions, second substances
!
!*ProdMat       --> 2-dimensional array containing stoichiometric coefficients of the
!               products. first index are reactions, second substances
!
!*RateVec       --> 1-dimensional array containing rate coefficients for every reaction
!
! Omega         --> number of substances, width of matrices EdMat, ProdMat and vector ConcVec
!
! N             --> number of reactions, height of matrices EdMat, ProdMat and vector RateVec
!
! I             --> serial number of iterators over reactions [1 ... N]
!
! lambda        --> serial number of iterators over substances [1 ... Omega]
!
! SubsNum       --> number of substance to calculate differential equation from.
!               also from [1 ... Omega]
!
! StartTime     --> time at which to start the integration step
!
!*dTime         --> integration step length
!
!*Finaltime     --> time at which to stop integration
!
! NSteps        --> number of integration steps
!
! IntStep       --> serial number of integration [1 ... NSteps]
!
!*PipeLength    --> length of the pipe z = [0 ... l]
!
!*PipeRadius    --> radius of the pipe x = [-r, ... , 0 , ... , +r] and y = [-r, ... , 0 , ... , +r]
!
!*NGridXY       --> Number of points in x and y direction or in how many segments to split the radius
!
!*PipeConc      --> Array of type (2*NGridXY+1:2*NGridXY+1:INT(Pipelength/Piperadius/NGridXY):Omega:2)
!
!*vmax          --> maximum speed of flow in the middle of the pipe (laminar flow)
!
!*vadd          --> constant to be added to the flow at every point
!
! PointX/Y/Z    --> integer describing position on grid
!
! PipeArea      --> 2D Array in xy plane, containing logical values, telling you if you are inside or outside the pipe
!
!*method        --> numerical integration method
!
!*Inflow        --> is the concentration in the first xy plane constant == is there inflow of educts
!


SUBROUTINE PIPE_3DEVOLVE(PipeLength, PipeRadius, NGridXY, vmax, vadd, dTime, &
FinalTime, EdMat, ProdMat, RateVec, PipeConc, method, Inflow)

USE mpi

INCLUDE 'pipe_3devolve.fh'


INTEGER :: ierr, NProcs, ProcID, IProc, GenericTag = 1

INTERFACE
        SUBROUTINE FIND_PIPE(vmax, NGridXY, PipeArea)
                IMPLICIT NONE
                REAL*8, INTENT(IN) ::  vmax
                INTEGER, INTENT(IN) :: NGridXY
                REAL*8, PARAMETER :: vadd = 0.0d0
                LOGICAL, DIMENSION(-NGridXY:+NGridXY,-NGridXY:NGridXY), INTENT(OUT) :: PipeArea
                INTEGER :: PointX, PointY
                REAL*8 :: PipeTest, FLOW_PROFILE
        END SUBROUTINE
END INTERFACE

CALL MPI_INIT(ierr)

!!=============!!
!! DEBUG START !!
!!=============!!

WRITE(*, *) "============================"
WRITE(*, *) "|| Here is PIPE_3DEVOLVE  ||"
WRITE(*, *) "============================"

!!=============!!
!!  DEBUG END  !!
!!=============!!

!!!!!!!!!!!!!!!!!!!!!!!!!!
!! INITIALIZE VARIABLES !!
!!!!!!!!!!!!!!!!!!!!!!!!!!

NGridZ = SIZE(PipeConc,3)
Omega = SIZE(PipeConc,4)
N = SIZE(RateVec)
NSteps = NINT(FinalTime / dTime)
lambda = 1
I = 1
SubsNum = 1
IF (method < 1 .OR. method > 1) THEN
        RETURN
END IF

!!=============!!
!! DEBUG START !!
!!=============!!

WRITE(*, *) 
WRITE(*, *) "Now showing you which parameters i got as input"
WRITE(*, *)
WRITE(*, *)
WRITE(*, *) "Here are the start conditions"
WRITE(*, *) "NGridXY = ", NGridXY
WRITE(*, *) "vmax = ", vmax
WRITE(*, *) "vadd = ", vadd
WRITE(*, *) "PipeLength = ", PipeLength
WRITE(*, *) "PipeRadius = ", PipeRadius
WRITE(*, *) "dTime = ", dTime
WRITE(*, *) "FinalTime = ", FinalTime
WRITE(*, *) "method = ", method
WRITE(*, *) "NGridXY = ", NGridXY
WRITE(*, *) "NgridZ = ", NGridZ
WRITE(*, *) "Omega = ", Omega
WRITE(*, *) "N = ", N
WRITE(*, *)
WRITE(*, *)
WRITE(*, *) "EdMat"
DO i = 1, 2
        WRITE(*, *) "        ", EdMat(i,:)
END DO
WRITE(*, *)
WRITE(*, *)
WRITE(*, *) "ProdMat"
DO i = 1, 2
        WRITE(*, *) "        ", ProdMat(i,:)
END DO
WRITE(*, *)
WRITE(*, *)
WRITE(*, *) "RateVec"
DO i = 1, 2
        WRITE(*, *) "        ", RateVec(i)
END DO
WRITE(*, *)
WRITE(*, *)
WRITE(*, *) "PipeConc"
WRITE(*, *) "    Size in x", SIZE(PipeConc,1)
WRITE(*, *) "    Size in y", SIZE(PipeConc,2)
WRITE(*, *) "    Size in z", SIZE(PipeConc,3)
WRITE(*, *) "    Size in Substances", SIZE(PipeConc,4)
WRITE(*, *) "    Size in time", SIZE(PipeConc,5)

WRITE(*, *) "    in xy plane at pipe start (PointZ = 0)"
DO lambda = 1, Omega
        WRITE(*, *) "    Substance", lambda
        DO PointX = -NgridXY, NGridXY
        WRITE(*, *) "    ", PipeConc(PointX,-NGridXY:,1,lambda,1)
        END DO
        WRITE(*, *)
END DO
WRITE(*, *) "    in xy plane at the pipe end (PointZ = NGridZ)"
DO lambda = 1, Omega
        WRITE(*, *) "    Substance", lambda
        DO PointX = -NgridXY, NGridXY
        WRITE(*, *) "    ", PipeConc(PointX,-NGridXY:,NGridZ,lambda,1)
        END DO
        WRITE(*, *)
END DO

!!=============!!
!!  DEBUG END  !!
!!=============!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! INITIALIZE GRID AND PIPE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




CALL MPI_FINALIZE(ierr)

END SUBROUTINE
