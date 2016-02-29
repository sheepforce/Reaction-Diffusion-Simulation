PROGRAM REFORM_CONC
IMPLICIT NONE
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: PipeWrite
INTEGER :: NGridXY, NGridZ, lambda, Omega, NSteps, PointX, PointY, PointZ, IntStep, WriteInt
CHARACTER*12 :: IntName, SubsName


WRITE(*, *) "Converting the binary output to plottable ASCII"
WRITE(*, *) "-----------------------------------------------"

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

CALL SYSTEM('mkdir -p out/plot')

DO IntStep = WriteInt, NSteps, WriteInt
	
	READ(100) PipeWrite

	WRITE(IntName, FMT='(I0)') IntStep

	OPEN(UNIT=202, FILE="out/plot/t_"//TRIM(IntName)//".dat")
	DO PointX = -NgridXY, +NGridXY
		DO PointY = -NgridXY, +NgridXY
			DO PointZ = 1, NGridZ
				WRITE(202, *) PointX, PointY, PointZ, PipeWrite(PointX,PointY,PointZ,1:)
			END DO

		END DO
	END DO
	CLOSE(UNIT=202)

	PipeWrite = 0.0d0
END DO

CLOSE(UNIT=100)

END PROGRAM
