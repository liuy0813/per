PROGRAM TWMODEL
!LIUY TRY TO MAKE IT OPENMP IN 2016
USE MOD_GLOBAL
USE MOD_KINDS
USE MOD_READDATA
IMPLICIT NONE

CHARACTER(LEN=100)  :: PARAMFILE
CHARACTER(LEN=100)  :: AA
PARAMFILE='run_data.dat'
AA='ALL'
CALL READINPUT(PARAMFILE,AA)
PAUSE
END PROGRAM TWMODEL

