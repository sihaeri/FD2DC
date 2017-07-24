!    This file is part of FD2DC.
!
!    FD2DC is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    FD2DC is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with FD2DC.  If not, see <http://www.gnu.org/licenses/>.

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
