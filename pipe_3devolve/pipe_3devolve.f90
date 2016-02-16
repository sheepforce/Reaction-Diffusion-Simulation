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
CHARACTER*12 :: IntName, SubsName

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
	WRITE(*, *) "ERROR :  Requested not implemented integration method"
        STOP 201
END IF


!!//////////////!!
!! OUTPUT START !!
!!\\\\\\\\\\\\\\!!
IF (ProcID == 0) THEN
	1000 FORMAT (A80)
	1001 FORMAT (4X, A)
	1002 FORMAT (1X, I2.2, 1X)
	1003 FORMAT (1X, D12.6, 1X)

	WRITE(*, 1000) "==============================="
	WRITE(*, 1000) "||   PIPE 3D EVOLVE MODULE   ||"
	WRITE(*, 1000) "==============================="
	WRITE(*, FMT='(2/)')
	WRITE(*, 1001) "Input Parameters"
	WRITE(*, 1001) "----------------"
	WRITE(* ,*)
	WRITE(*, 1001) "matrix of stochiometric coefficients of the educts --> EdMat"
	DO I = 1, N
		WRITE(*, FMT='(8X)', ADVANCE='NO')
		DO lambda = 1, Omega
			WRITE(*, 1002, ADVANCE='NO') EdMat(I, lambda)
		END DO
		WRITE(*, *)
	END DO
	WRITE(*, *)
	WRITE(*, 1001) "matrix of stochiometric coefficients of the products --> ProdMat"
        DO I = 1, N
                WRITE(*, FMT='(8X)', ADVANCE='NO')
        	DO lambda = 1, Omega
                	WRITE(*, 1002, ADVANCE='NO') ProdMat(I, lambda)
	        END DO
                WRITE(*, *)
        END DO
	WRITE(*, *)
	WRITE(*, 1001) "vector of the the rate coefficients of reactions --> RateVec"
	DO I = 1, N
		WRITE(*, FMT='(8X)', ADVANCE='NO')
		WRITE(*, 1003) RateVec(I)
	END DO
	WRITE(*, *)
	WRITE(*, 1001) "vector of the diffusion coefficients of substances --> DiffVec"
	DO lambda = 1, Omega
		WRITE(*, FMT='(8X)', ADVANCE='NO')
		WRITE(*, 1003) DiffVec(lambda)

	END DO
	WRITE(*, *)
	WRITE(*, *)
	WRITE(*, FMT='(4X, A, D12.6)') "length of the pipe : ", PipeLength
	WRITE(*, FMT='(4X, A, D12.6)') "radius of the pipe : ", PipeRadius
	WRITE(*, FMT='(4X, A, I12)')   "radial grid points : ", NGridXY
	WRITE(*, FMT='(4X, A, I12)')   "grid points in z   : ", NGridZ
	WRITE(*, FMT='(4X, A, D12.6)') "size of time steps : ", dTime
	WRITE(*, FMT='(4X, A, I12)')   "integration steps  : ", NSteps
	WRITE(*, FMT='(4X, A, D12.6)') "flow speed vmax    : ", vmax
	WRITE(*, FMT='(4X, A, D12.6)') "flow speed vadd    : ", vadd
	WRITE(*, FMT='(4X, A, L12, /,/,/)') "pipe has inflow	   : ", InFlow
END IF
!!//////////////!!
!!  OUTPUT END  !!
!!\\\\\\\\\\\\\\!!

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

!!//////////////!!
!! OUTPUT START !!
!!\\\\\\\\\\\\\\!!
IF (ProcID == 0) THEN
	WRITE(*, 1001) "Calculation"
	WRITE(*, 1001) "-----------"
	WRITE(*, *)
	WRITE(*, 1001) "Found pipe area"
	DO PointX = -NGridXY, +NGridXY
		WRITE(*, FMT='(8X)', ADVANCE='NO')
		DO PointY = -NGridXY, +NGridXY
			WRITE(*, FMT='(L1, 1X)', ADVANCE='NO') PipeArea(PointX,PointY)
		END DO
		WRITE(*, *)
	END DO
	WRITE(*, *)
	WRITE(*, 1001) "Writing Speed Profile in the pipe to SpeedProfile.dat"
	OPEN(UNIT=201, FILE="SpeedProfile.dat", ACTION='write')
        DO PointX = -NGridXY, +NGridXY
	        DO PointY = -NGridXY, +NGridXY
			IF (PipeArea(PointX,PointY) .EQV. .TRUE.) THEN
				WRITE(UNIT=201, FMT='(D12.6, 1X)', ADVANCE='NO') FlowMat(PointX,PointY)
			ELSE
				WRITE(UNIT=201, FMT='(6X, A1, 6X)', ADVANCE='NO') "F"
			END IF
        	END DO
                WRITE(201, *)
        END DO	
	CLOSE(UNIT=201)
	WRITE(*, FMT='(/,/)')
END IF

!!//////////////!!
!!  OUTPUT END  !!
!!\\\\\\\\\\\\\\!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! INTEGRATION OVER TIME !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

! distribute PipeConc(:,:,:,:,1) to all MPI threads, number of elements to send
BCastSendCount = SIZE(PipeConc) / 2


!! DEBUG START
OPEN(UNIT=999, FILE="DeltaConc.debug")
!! DEBUG END

DO IntStep = 1, NSteps											! Integrate over time
	!!//////////////!!
	!! OUTPUT START !!
	!!\\\\\\\\\\\\\\!!
	IF (ProcID == 0 .AND. MOD(IntStep, NSteps / 100) == 0) THEN
		WRITE(*, FMT='(4X, A17, I12, A4, I12)') "Integration step ", IntStep, " of ", NSteps	! write integration step to output
	END IF

	! create the concentration outputs
	IF (ProcID == 0 .AND. MOD(IntStep, 5) == 0) THEN						! write all 5 time steps output
		WRITE(IntName, FMT='(I0)') IntStep							! write current integration step to String IntName
		CALL SYSTEM('mkdir -p out/t_'//TRIM(IntName))						! create output directory and in this directory directories for time steps
		
		DO lambda = 1, Omega									! iterate over substances
			WRITE(SubsName, FMT='(I0)') lambda						
			OPEN(UNIT=202, FILE="out/t_"//TRIM(IntName)//"/Subs"//TRIM(SubsName)//".dat")	! write different outputs for different susbstances 
			DO PointZ = 1, NGridZ
				DO PointX = -NGridXY, +NGridXY
					WRITE(202, *) PipeConc(PointX,-NGridXY:,PointZ,lambda,1)
				END DO
			WRITE(202, *)
			END DO
			CLOSE(UNIT=202)
		END DO

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
		
		DummyFLowInMat = 0.0d0									! make sure DummyFLowInMat is reseted

!		IF (PointZ == 1) THEN
!			DummyFlowInMat(1,:) = 0.0d0
!		ELSE
!			DummyFlowInMat(1,:) = PipeConc(PointX,PointY,PointZ-1,:,1)
!		END IF
!		DummyFlowInMat(2,:) = PipeConc(PointX,PointY,PointZ,:,1)
!
!		CALL FLOW_INT(PipeLength, NGridZ, FlowMat(PointX,PointY), &				! calculates the change in concentration by laminar flow
!		DummyFlowInMat, dTime, DeltaConc)							! using the current point and the previous point in z direction
!
!		DO lambda = 1, Omega
!			PipeConc(PointX,PointY,PointZ,lambda,2) = DeltaConc(lambda) &
!			+ PipeConc(PointX,PointY,PointZ,lambda,2)					! add changes in concetration due to laminar flow to the changes in concentrations
!		END DO
!
!		!! DEBUG START
!		IF (ProcID == 0 .AND. IntStep == 1 .AND. PointZ == 1) THEN
!			WRITE(999, FMT='(D12.6, 2X)', ADVANCE='NO') DeltaConc(1)
!			
!		END IF
!		
!		IF (ProcID == 0 .AND. IntStep == 1 .AND. PointX == 0 .AND. PointY == 0) THEN
!			WRITE(*, *) "DummyFlowInMat(:,1)", DummyFlowInMat(:,1)
!			WRITE(*, *) FlowMat(PointX,PointY)
!		END IF
		!! DEBUG END 

		DummyFlowInMat = 0.0d0									! reset DummyFlowInMat for next steps
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
	!! DEBUG START
	WRITE(999, *) 
	!! DEBUG END
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
		IF (InFlow .EQV. .TRUE.) THEN
			PipeConc(:,:,2:,:,1) = PipeConc(:,:,2:,:,1) + PipeConc(:,:,2:,:,2)		! all except the inflow layer changes
			PipeConc(:,:,1,:,1) = PipeConc(:,:,1,:,1)					! the inflow layer stays constant
		ELSE
			PipeConc(:,:,:,:,1) = PipeConc(:,:,:,:,1) + PipeConc(:,:,:,:,2)
		END IF



	END IF

END DO													! end integrate over time

!! DEBUG START
CLOSE(UNIT=999)
!! DEBUG END

CALL MPI_FINALIZE(ierr)

END SUBROUTINE
