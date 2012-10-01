MODULE modfd_openclose_ops

PRIVATE
PUBLIC :: fd_open,fd_open_sres

CONTAINS
!===========================================================
SUBROUTINE fd_open

USE parameters,  ONLY : max_char_len,set_unit,out_unit,grd_unit,eres_unit
USE shared_data, ONLY : problem_name,problem_len
IMPLICIT NONE
INTEGER :: ierror
CHARACTER(LEN = 1) :: dummy
PRINT *, 'ENTER PROBLEM NAME:  '
READ(*,'(A20)') problem_name

problem_len = LEN_TRIM(ADJUSTL(problem_name))
OPEN (UNIT=set_unit,FILE=problem_name(1:problem_len)//'.set',IOSTAT=ierror)
IF(ierror /= 0)GOTO 100
OPEN (UNIT=out_unit,FILE=problem_name(1:problem_len)//'.out',IOSTAT=ierror)
IF(ierror /= 0)GOTO 100
OPEN (UNIT=grd_unit,FILE=problem_name(1:problem_len)//'.grd',IOSTAT=ierror)
IF(ierror /= 0)GOTO 100
OPEN (UNIT=eres_unit,FILE=problem_name(1:problem_len)//'.ere',FORM='UNFORMATTED',IOSTAT=ierror)
IF(ierror /= 0)GOTO 100

100 CONTINUE
IF(ierror /= 0)THEN
  WRITE(*,*)'fd_open: cannot open a required file'
  WRITE(*,*)'Press any key...'
  READ (*,*)dummy
  STOP
ENDIF
END SUBROUTINE fd_open

SUBROUTINE fd_open_sres

USE parameters,  ONLY : sres_unit
USE shared_data, ONLY : problem_name,problem_len

IMPLICIT NONE

OPEN(UNIT=sres_unit,FILE=problem_name(1:problem_len)//'.ers',FORM='UNFORMATTED')

END SUBROUTINE fd_open_sres

END MODULE modfd_openclose_ops