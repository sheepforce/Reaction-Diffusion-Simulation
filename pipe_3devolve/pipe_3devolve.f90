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
FinalTime, EdMat, ProdMat, RateVec, DiffVec, PipeConc, method, Inflow)

USE mpi

IMPLICIT NONE

INCLUDE 'pipe_3devolve.fh'

INTEGER :: ierr, NProcs, ProcID, IProc, GenericTag = 1, BCastSendCount
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status

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

INTERFACE
	SUBROUTINE FLOW_INT(PipeLength, NGridZ, vFlowXY, FlowConcMat, dTime, DeltaConc)

		IMPLICIT NONE
		REAL*8, INTENT(IN) :: PipeLength, vFlowXY, dtime
		REAL*8, DIMENSION(:,:), INTENT(IN) :: FlowConcMat
		INTEGER, INTENT(IN) :: NGridZ

		REAL*8, DIMENSION(SIZE(FlowConcMat,2)), INTENT(OUT) :: DeltaConc

		REAL*8 :: DeltaL, DeltaN
		INTEGER :: Omega, I, lambda
	END SUBROUTINE
END INTERFACE

INTERFACE
	SUBROUTINE DIFF_INT(ConcEnv, DiffVec, DeltaConc)

		IMPLICIT NONE
		REAL*8, DIMENSION(-1:,-1:,-1:,:), INTENT(IN) :: ConcEnv
		REAL*8, DIMENSION(:), INTENT(IN) :: DiffVec

		REAL*8, DIMENSION(1:SIZE(DiffVec)), INTENT(OUT) :: DeltaConc

		REAL*8 :: LAPLACIAN
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

DO PointX = -NGridXY, +NGridXY
	DO PointY = -NgridXY, +NGridXY
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

DO IntStep = 1, NSteps											! Integrate over time
	!!//////////////!!
	!! OUTPUT START !!
	!!\\\\\\\\\\\\\\!!
	IF (ProcID == 0 .AND. MOD(IntStep, NSteps / 100) == 0) THEN
		WRITE(*, *) "Now at time step", IntStep, " of ", NSteps					! write the integration step to output
	END IF
	!!//////////////!!
	!!  OUTPUT END  !!
	!!\\\\\\\\\\\\\\!!

	CALL MPI_BCAST(PipeConc, BCastSendCount, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)		! distribute the new start concentrations to all threads
	
	DeltaConc = 0.0d0										! reset DeltaConc after each integration step
	PipeConc(:,:,:,:,2) = 0.0d0									! reset PipeConc(:,:,:,:,2) after each integration step

	DO PointX = -NgridXY, +NGridXY									! iterate over points in X
	DO PointY = -NGridXY, +NgridXY									! iterate over points in Y

	! if we are not in the pipe area at the current point we do not need to solve nothing...
	IF (PipeArea(PointX,PointY) .EQV. .FALSE.) THEN
		CYCLE
	END IF		

	DO PointZ = ProcID + 1, NGridZ, NProcs								! iterate over points in Z, distributed with MPI


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
			DummyFlowInMat(2,:) = PipeConc(PointX,PointY,PointZ,:,1)			! the concentrations at the start of the pipe are the chosen concentrations of the pipe...
			
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

		DummyDiffInMat = 0.0d0									! intialize dummy input for DIFF_INT, that will be adapted to constraints due to boundaries of the pipe
		DummyDiffInMat = PipeConc(PointX-1:PointX+1,PointY-1:PointY+1,PointZ-1:PointZ+1,:,1)	! if neighbouring points are not at the boundaries this is simply the corresponding segment of PipeConc

		! now testing if any of the relevant neighbouring points is outside an replace its values with the value of the center
		IF (PipeArea(PointX-1,PointY) .EQV. .FALSE.) THEN					! is PointX-1 outside
			DummyDiffInMat(-1,:,:,:) = PipeConc(PointX,PointY-1:PointY+1,PointZ-1:PointZ+1,:,1)
		END IF
		IF (PipeArea(PointX+1,PointY) .EQV. .FALSE.) THEN					! is PointX+1 outside
			DummyDiffInMat(+1,:,:,:) = PipeConc(PointX,PointY-1:PointY+1,PointZ-1:PointZ+1,:,1)
		END IF
		IF (PipeArea(PointX,PointY-1) .EQV. .FALSE.) THEN					! is PointY-1 outside
			DummyDiffInMat(:,-1,:,:) = PipeConc(PointX-1:PointX+1,PointY,PointZ-1:PointZ+1,:,1)
		END IF
		IF (PipeArea(PointX,PointY+1) .EQV. .FALSE.) THEN					! is PointY+1 outside
			DummyDiffInMat(:,+1,:,:) = PipeConc(PointX-1:PointX+1,PointY,PointZ-1:PointZ+1,:,1)
		END IF
		IF (PointZ - 1 < 1) THEN								! are we at the start of the pipe
			DummyDiffInMat(:,:,PointZ-1,:) = PipeConc(PointX-1:PointX+1,PointY-1:PointY+1,PointZ,:,1)
		END IF
		IF (PointZ + 1 > NGridZ) THEN								! are we at the end of the pipe
			DummyDiffInMat(:,:,PointZ+1,:) = PipeConc(PointX-1:PointX+1,PointY-1:PointY+1,PointZ,:,1)
		END IF

		CALL DIFF_INT(DummyDiffInMat, DiffVec, DeltaConc)

		DO lambda = 1, Omega
                        PipeConc(PointX,PointY,PointZ,lambda,2) = DeltaConc(lambda) &
                        + PipeConc(PointX,PointY,PointZ,lambda,2)                                       ! add changes in concetration due to laminar flow to the changes in concentrations
                END DO

		DeltaConc = 0.0d0


	END DO												! stop iterate over points in Z	
	END DO												! stop iterate over points in Y
	END DO												! stop iterate over points in X


	! MPI_SEND for all processes that are not root to send values in PipeConc(:,:,:,:,2) back to root
	! Root thread sums them up and calculates new concentrations
	IF (ProcID == 0) THEN
		DO IProc = 1, NProcs - 1
			CALL MPI_RECV(LocConcChange, BCastSendCount, MPI_DOUBLE_PRECISION, &		! receive array with concentration changes of every thread
			MPI_ANY_SOURCE, GenericTag, MPI_COMM_WORLD, status, ierr)

			PipeConc(:,:,:,:,2) = PipeConc(:,:,:,:,2) + LocConcChange			! sum up the concentration changes of root thread and the other threads

			LocConcChange = 0.0d0
		END DO
	ELSE IF (ProcID /= 0) THEN
			LocConcChange = PipeConc(:,:,:,:,2)
			CALL MPI_SEND(LocConcChange, BCastSendCount, MPI_DOUBLE_PRECISION, &		! send array with changes in concentrations to the root thread
			0, GenericTag, MPI_COMM_WORLD, status, ierr)
			
	END IF

	! calculating the concentrations for the new round by using the concentrations of this round and the concentration changes
	IF (ProcID == 0) THEN
		PipeConc(:,:,:,:,1) = PipeConc(:,:,:,:,1) + PipeConc(:,:,:,:,2)
	END IF


END DO													! end integrate over time

CALL MPI_FINALIZE(ierr)

END SUBROUTINE
