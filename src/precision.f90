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

MODULE precision
IMPLICIT NONE
!====================never change these values
INTEGER, PARAMETER :: r_double       = SELECTED_REAL_KIND(12,100) ! don't touch
INTEGER, PARAMETER :: r_single_fixed = SELECTED_REAL_KIND(6,30)   ! don't touch
INTEGER, PARAMETER :: i_single       = SELECTED_INT_KIND(8)       ! don't touch
INTEGER, PARAMETER :: i_double       = SELECTED_INT_KIND(18)      ! don't touch
! INTEGER(8) values -9223372036854775808 to 9223372036854775807

!====THE ONLY COMPLILE TIME SETTING===================
INTEGER, PARAMETER :: r_single = r_double !--don't touch for this program
INTEGER, PARAMETER :: r_single_p = 12 ! set to the value chosen for r_single
INTEGER, PARAMETER :: r_single_r = 100 ! set to the value chosen for r_single
!=====================================================
END MODULE precision

