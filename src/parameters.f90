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

MODULE parameters

!--Number of variables to solve 
INTEGER,PARAMETER       :: nphi = 4
INTEGER,PARAMETER       :: max_char_len = 20

!--File units
INTEGER,PARAMETER       :: set_unit = 500 !--Problem setup file
INTEGER,PARAMETER       :: out_unit = set_unit + 1 !--Output information at each iteration
INTEGER,PARAMETER       :: grd_unit = out_unit + 1 !--grid file
INTEGER,PARAMETER       :: sres_unit = grd_unit + 3240 !--read restart file
INTEGER,PARAMETER       :: eres_unit = sres_unit + 1 !--write restart file
INTEGER,PARAMETER       :: plt_unit = eres_unit + 1 !-- tect plot results
INTEGER,PARAMETER       :: ubar_unit = plt_unit + 1 
INTEGER,PARAMETER       :: vbar_unit = ubar_unit + 1 
INTEGER,PARAMETER       :: u2bar_unit = vbar_unit + 1 
INTEGER,PARAMETER       :: v2bar_unit = u2bar_unit + 1 
INTEGER,PARAMETER       :: uvbar_unit = v2bar_unit + 1 
INTEGER,PARAMETER       :: init_field_unit = uvbar_unit + 1 
INTEGER,PARAMETER       :: alloc_create = 1
INTEGER,PARAMETER       :: alloc_destroy = alloc_create+1
INTEGER,PARAMETER       :: max_len_tecline  = 32000

INTEGER,PARAMETER       :: no_force = -1
INTEGER,PARAMETER       :: force_predict = 1
INTEGER,PARAMETER       :: force_correct = force_predict + 1

INTEGER,PARAMETER       :: solver_sparsekit = 1
INTEGER,PARAMETER       :: solver_sip       = solver_sparsekit + 1
INTEGER,PARAMETER       :: solver_cg        = solver_sip + 1
INTEGER,PARAMETER       :: solver_hypre     = solver_cg + 1

INTEGER,PARAMETER       :: maxmove = 10

INTEGER,PARAMETER       :: do_collect_stat = 1
INTEGER,PARAMETER       :: end_collect_stat = do_collect_stat + 1

INTEGER,PARAMETER       :: use_GPU_no  = 0
INTEGER,PARAMETER       :: use_GPU_yes = use_GPU_no + 1
INTEGER,PARAMETER       :: OPERROR = -1
INTEGER,PARAMETER       :: SOLVER_FAILED = OPERROR-1
INTEGER,PARAMETER       :: BCG_PTR_NOT_INITIALIZED=SOLVER_FAILED-1
INTEGER,PARAMETER       :: OPSUCCESS = 0
INTEGER,PARAMETER       :: SOLVER_DONE = 1

INTEGER,PARAMETER       :: OUTER_ITR_DONE = -1

END MODULE parameters
