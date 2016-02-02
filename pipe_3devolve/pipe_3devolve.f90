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


SUBROUTINE PIPE_3DEVOLVE(PipeLength, PipeRadius, NGridXY, NGridZ, vmax, vadd, dTime, &
FinalTime, EdMat, ProdMat, RateVec, PipeConc, method, Inflow)

USE mpi

IMPLICIT NONE

INCLUDE 'pipe_3devolve.fh'

INTEGER :: ierr, NProcs, ProcID, IProc, GenericTag = 1, BCastSendCount

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

INTERFACE
        SUBROUTINE RKS_INT(EdMat, ProdMat, RateVec, ConcVec, dTime, method, DeltaConc)
                IMPLICIT NONE
                INTEGER, DIMENSION(:,:), INTENT(IN) :: EdMat, ProdMat
                REAL*8, DIMENSION(:), INTENT(IN) :: RateVec , ConcVec
                REAL*8, INTENT(IN) :: dTime
                INTEGER, INTENT(IN) :: method

		REAL*8, DIMENSION(1:SIZE(ConcVec)), INTENT(OUT) :: DeltaConc

                INTEGER :: Omega, N, I, lambda, SubsNum
                REAL*8, DIMENSION(:), ALLOCATABLE :: ODEVec
                REAL*8 :: RHS_ODE
        END SUBROUTINE
END INTERFACE

CALL MPI_INIT(ierr)											! fork the processes with MPI
CALL MPI_COMM_RANK(MPI_COMM_WORLD, ProcID, ierr)							! get the ID of the current process
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NProcs, ierr)							! how many MPI processes are there over all

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

Omega = SIZE(PipeConc,4)
N = SIZE(RateVec)
NSteps = NINT(FinalTime / dTime)
lambda = 1
I = 1
SubsNum = 1
IF (method < 1 .OR. method > 1) THEN
        RETURN
END IF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! INITIALIZE GRID AND PIPE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL FIND_PIPE(vmax,NGridXY, PipeArea)									! find the grid points inside the pipe and outside the pipe


FlowMat = 0.0d0

DO PointX = -NgridXY, NGridXY
	DO PointY = -NgridXY, NGridXY
		IF (PipeArea(PointX,PointY) .EQV. .TRUE.) THEN
			FlowMat(PointX,PointY) = FLOW_PROFILE(PointX, PointY, vmax, vadd, NGridXY)	! calculate flow profile inside the pipe, stored in FlowMat
		END IF
	END DO
END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! INTEGRATION OVER TIME !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

! distribute PipeConc(:,:,:,:,1) to all MPI threads, number of elements to send
BCastSendCount = SIZE(PipeConc) / 2

! MPI threads > 0 calculate the end of the pipe in z direction 
ZGridElementsMPI = NgridZ / NProcs

DO IntStep = 1, NSteps											! Integrate over time
	CALL MPI_BCAST(PipeConc, BCastSendCount, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)		! distribute the new start concentrations to all threads
	
	DeltaConc = 0.0d0										! reset DeltaConc after each integration step
	PipeConc(:,:,:,:,2) = 0.0d0									! reset PipeConc(:,:,:,:,2) after each integration step

	DO PointX = -NgridXY, +NGridXY									! iterate over points in X
	DO PointY = -NGridXY, +NgridXY									! iterate over points in Y
	DO PointZ = (ZGridElementsMPI * ProcID) + 1, ZGridElementsMPI * (ProcID + 1)			! iterate over points in Z

		! if we are not in the pipe area at the current point we do not need to solve nothing...
		IF (PipeArea(PointX,PointY) .EQV. .FALSE.) THEN
			EXIT
		END IF		


		!!!!!!!!!!!!!!!!!!!!!!!
		!! REACTION KINETICS !!
		!!!!!!!!!!!!!!!!!!!!!!!

		CALL RKS_INT(EdMat, ProdMat, RateVec, PipeConc(PointX,PointY,PointZ,:,1), &		! calculates the change in concentrations by reaction 
		dTime, method, DeltaConc)								! kinetics and stores them in DeltaConc

		DO lambda = 1, Omega
			PipeConc(PointX,PointY,PointZ,lambda,2) = DeltaConc(lambda) &
			+ PipeConc(PointX,PointY,PointZ,lambda,2)					! store concentration changes due to reaction kinetics in PipeConc, Layer 2
		END DO

		DeltaConc = 0.0d0									! reset DeltaConc for next integration step


		!!!!!!!!!!!!!!!!!!!!!
		!! FLOW OF LIQUIDS !!
		!!!!!!!!!!!!!!!!!!!!!

		IF (PointZ == 1) THEN
			DummyFlowInMat(1,:) = 0.0d0							! the concentrations before the pipe starts are 0, so that there can be no inflow
			DummyFlowInMat(1,:) = PipeConc(PointX,PointY,PointZ,:,1)			! the concentrations at the start of the pipe are the choes concentrations of the pipe...
			
			CALL FLOW_INT(PipeLength, NGridZ, FlowMat(PointX,PointY), &			! calculates the change in concentration by laminar flow
			DummyFlowInMat, dTime, DeltaConc)						! using the current point and the previous point in z direction
		ELSE
			CALL FLOW_INT(PipeLength, NGridZ, FlowMat(PointX,PointY), &			! calculates the change in concentration by laminar flow
			PipeConc(PointX,PointY,PointZ-1:PointZ,:,1), dTime, DeltaConc)			! using the current point and the previous point in z direction
		END IF

		DO lambda = 1, Omega
			PipeConc(PointX,PointY,PointZ,lambda,2) = DeltaConc(lambda) &
			+ PipeConc(PointX,PointY,PointZ,lambda,2)					! add changes in concetration due to laminar flow to the changes in concentrations
		END DO

		DeltaConc = 0.0d0									! reset DeltaConc for next integration step


		!!!!!!!!!!!!!!!!
		!! DIFFUSSION !!
		!!!!!!!!!!!!!!!!

	END DO												! stop iterate over points in Z	
	END DO												! stop iterate over points in Y
	END DO												! stop iterate over points in X


END DO													! end integrate over time

CALL MPI_FINALIZE(ierr)

END SUBROUTINE
