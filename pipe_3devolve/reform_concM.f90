PROGRAM REFORM_CONC
IMPLICIT NONE
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: PipeWrite
INTEGER :: NGridXY, NGridZ, lambda, Omega, NSteps, PointX, PointY, PointZ, IntStep, WriteInt, UOpen
CHARACTER*12 :: IntName, SubsName


WRITE(*, *) "Converting the binary output to plottable ASCII, python version"
WRITE(*, *) "---------------------------------------------------------------"

! open the output from PIPE_3DEVOLVE and read the binary array
OPEN(UNIT=100, FILE='Pipe.dat', FORM='UNFORMATTED', ACCESS='STREAM', ACTION='READ', STATUS='OLD')

! reading the header an determine parameters
READ(100) NGridXY
READ(100) NGridZ
READ(100) Omega
READ(100) NSteps
READ(100) WriteInt

! write the read header out
WRITE(*, *) "NGridXY  :", NGridXY
WRITE(*, *) "NGridZ   :", NGridZ
WRITE(*, *) "Omega    :", Omega
WRITE(*, *) "NSteps   :", NSteps
WRITE(*, *) "WriteInt :", WriteInt

ALLOCATE(PipeWrite(-NGridXY:NGridXY,-NGridXY:+NGridXY,NgridZ,Omega))
PipeWrite = 0.0d0


CALL SYSTEM('mkdir -p out')


DO lambda = 1, Omega
	UOpen = 200 + lambda
	WRITE(SubsName, FMT='(I0)') lambda
	OPEN(UNIT=Uopen, FILE="out/Subs"//TRIM(SubsName)//".dat")
END DO


DO IntStep = WriteInt, NSteps, WriteInt
	READ(100) PipeWrite
	WRITE(IntName, FMT='(I0)') IntStep

	DO lambda = 1, Omega
		WRITE(SubsName, FMT='(I0)') lambda
		UOpen = lambda + 200

		IF (IntStep == WriteInt) THEN
			WRITE(UOpen, *) "Steps output :", NSteps / WriteInt
			WRITE(UOpen, *) "X            :", NGridXY * 2 + 1
			WRITE(UOpen, *) "Y            :", NGridXY * 2 + 1
			WRITE(UOpen, *) "Z            :", NGridZ
			WRITE(UOpen, *)
		END IF

		DO PointX = -NgridXY, +NGridXY
			DO PointY = -NgridXY, +NgridXY
				DO PointZ = 1, NGridZ
					WRITE(UOpen, *) PipeWrite(PointX,PointY,PointZ,lambda)
				END DO
			END DO
		END DO
		CLOSE(UNIT=202)
	END DO
END DO

DO lambda = 1, Omega
	UOpen = 200 + lambda
	CLOSE(UNIT=UOpen)
END DO


CLOSE(UNIT=100)

END PROGRAM
