MODULE modcusp_library_intrf

  USE, INTRINSIC :: iso_c_binding
  PRIVATE
  PUBLIC :: getInstance_cusp_biCGSTAB_solver, cusp_biCGSTAB_initDevice, &
       cusp_biCGSTAB_allocDevice, cusp_biCGSTAB_copyH2D_AInds, &
       cusp_biCGSTAB_copyH2D_system, cusp_biCGSTAB_copyD2H_x, &
       cusp_biCGSTAB_solveDev_system, cusp_biCGSTAB_getMonitor, &
       cusp_biCGSTAB_shutdown


  TYPE(c_ptr) :: cusp_biCGSTAB_solver_ptr = c_NULL_ptr
  INTEGER, PARAMETER :: CUSP_PTR_NOT_INITIALIZED = -4

  INTERFACE

     FUNCTION getInstance_cusp_biCGSTAB_solver_intrf() RESULT(ptr) BIND (C,name='getInstance_cusp_biCGSTAB_solver')
       USE, INTRINSIC :: iso_c_binding, ONLY : c_ptr
       IMPLICIT NONE
       TYPE(c_ptr) :: ptr
     END FUNCTION getInstance_cusp_biCGSTAB_solver_intrf

     FUNCTION cusp_biCGSTAB_initDevice_intrf(ptr, devID) BIND (C,name='cusp_biCGSTAB_initDevice_intrf')
       USE, INTRINSIC :: iso_c_binding, ONLY : c_ptr, c_int
       IMPLICIT NONE
       TYPE(c_ptr),VALUE    :: ptr
       INTEGER(c_int),VALUE :: devID
       INTEGER(c_int)       :: cusp_biCGSTAB_initDevice_intrf
     END FUNCTION cusp_biCGSTAB_initDevice_intrf

     FUNCTION cusp_biCGSTAB_allocDevice_intrf(ptr, n, m) BIND (C,name='cusp_biCGSTAB_allocDevice_intrf')
       USE, INTRINSIC :: iso_c_binding, ONLY : c_ptr, c_int
       IMPLICIT NONE
       TYPE(c_ptr),VALUE    :: ptr
       INTEGER(c_int),VALUE :: m,n
       INTEGER(c_int)       :: cusp_biCGSTAB_allocDevice_intrf
     END FUNCTION cusp_biCGSTAB_allocDevice_intrf

     FUNCTION cusp_biCGSTAB_copyH2D_AInds_intrf(ptr, rows, cols) &
          BIND (C,name='cusp_biCGSTAB_copyH2D_AInds_intrf')
       USE, INTRINSIC :: iso_c_binding, ONLY : c_ptr, c_int, c_double
       IMPLICIT NONE
       TYPE(c_ptr),VALUE    :: ptr
       INTEGER(c_int)       :: rows(*), cols(*)
       INTEGER(c_int)       :: cusp_biCGSTAB_copyH2D_AInds_intrf
     END FUNCTION cusp_biCGSTAB_copyH2D_AInds_intrf

     FUNCTION cusp_biCGSTAB_copyH2D_system_intrf(ptr, Avals, xHost, bHost) &
          BIND (C,name='cusp_biCGSTAB_copyH2D_system_intrf')
       USE, INTRINSIC :: iso_c_binding, ONLY : c_ptr, c_int, c_double
       IMPLICIT NONE
       TYPE(c_ptr),VALUE    :: ptr
       REAL(c_double)       :: xHost(*), bHost(*), Avals(*)
       INTEGER(c_int)       :: cusp_biCGSTAB_copyH2D_system_intrf
     END FUNCTION cusp_biCGSTAB_copyH2D_system_intrf

     FUNCTION cusp_biCGSTAB_copyD2H_x_intrf(ptr, xHost) &
          BIND (C,name='cusp_biCGSTAB_copyD2H_x_intrf')
       USE, INTRINSIC :: iso_c_binding, ONLY : c_ptr, c_int, c_double
       IMPLICIT NONE
       TYPE(c_ptr),VALUE    :: ptr
       REAL(c_double)       :: xHost(*)
       INTEGER(c_int)       :: cusp_biCGSTAB_copyD2H_x_intrf
     END FUNCTION cusp_biCGSTAB_copyD2H_x_intrf

     FUNCTION cusp_biCGSTAB_solveDev_system_intrf(ptr, relTol, absTol, maxItr) &
          BIND (C,name='cusp_biCGSTAB_solveDev_system_intrf')
       USE, INTRINSIC :: iso_c_binding, ONLY : c_ptr, c_int, c_double
       IMPLICIT NONE
       TYPE(c_ptr),VALUE    :: ptr
       REAL(c_double),VALUE :: relTol, absTol
       INTEGER(c_int),VALUE :: maxItr
       INTEGER(c_int)       :: cusp_biCGSTAB_solveDev_system_intrf
     END FUNCTION cusp_biCGSTAB_solveDev_system_intrf

     FUNCTION cusp_biCGSTAB_getMonitor_intrf(ptr, residual, nItr) &
          BIND (C,name='cusp_biCGSTAB_getMonitor_intrf')
       USE, INTRINSIC :: iso_c_binding, ONLY : c_ptr, c_int, c_double
       IMPLICIT NONE
       TYPE(c_ptr),VALUE    :: ptr
       REAL(c_double)       :: residual
       INTEGER(c_int)       :: nItr
       INTEGER(c_int)       :: cusp_biCGSTAB_getMonitor_intrf
     END FUNCTION cusp_biCGSTAB_getMonitor_intrf

     FUNCTION cusp_biCGSTAB_shutdown_intrf(ptr) &
          BIND (C,name='cusp_biCGSTAB_shutdown_intrf')
       USE, INTRINSIC :: iso_c_binding, ONLY : c_ptr, c_int, c_double
       IMPLICIT NONE
       TYPE(c_ptr),VALUE    :: ptr
       INTEGER(c_int)       :: cusp_biCGSTAB_shutdown_intrf
     END FUNCTION cusp_biCGSTAB_shutdown_intrf

  END INTERFACE

CONTAINS

  SUBROUTINE getInstance_cusp_biCGSTAB_solver()
    IMPLICIT NONE
    cusp_biCGSTAB_solver_ptr = getInstance_cusp_biCGSTAB_solver_intrf();
  END SUBROUTINE getInstance_cusp_biCGSTAB_solver

  SUBROUTINE cusp_biCGSTAB_initDevice(devID, err)
    IMPLICIT NONE
    INTEGER :: devID
    INTEGER, INTENT(out) :: err
    IF(C_ASSOCIATED(cusp_biCGSTAB_solver_ptr,c_NULL_ptr))THEN
       err = CUSP_PTR_NOT_INITIALIZED
       RETURN
    END IF

    err = cusp_biCGSTAB_initDevice_intrf(cusp_biCGSTAB_solver_ptr, devID)

  END SUBROUTINE cusp_biCGSTAB_initDevice

  SUBROUTINE cusp_biCGSTAB_allocDevice(n, m, err)
    IMPLICIT NONE
    INTEGER :: m, n
    INTEGER, INTENT(out) :: err
    IF(C_ASSOCIATED(cusp_biCGSTAB_solver_ptr,c_NULL_ptr))THEN
       err = CUSP_PTR_NOT_INITIALIZED
       RETURN
    END IF
    err = cusp_biCGSTAB_allocDevice_intrf(cusp_biCGSTAB_solver_ptr, n, m)

  END SUBROUTINE cusp_biCGSTAB_allocDevice

  SUBROUTINE cusp_biCGSTAB_copyH2D_AInds(rows, cols, err)
    IMPLICIT NONE
    INTEGER :: rows(*), cols(*)
    INTEGER, INTENT(out) :: err
    IF(C_ASSOCIATED(cusp_biCGSTAB_solver_ptr,c_NULL_ptr))THEN
       err = CUSP_PTR_NOT_INITIALIZED
       RETURN
    END IF

    err = cusp_biCGSTAB_copyH2D_AInds_intrf(cusp_biCGSTAB_solver_ptr, rows, cols)

  END SUBROUTINE cusp_biCGSTAB_copyH2D_AInds

  SUBROUTINE cusp_biCGSTAB_copyH2D_system(Avals, xHost, bHost, err)
    IMPLICIT NONE
    DOUBLE PRECISION :: xHost(*), bHost(*), Avals(*)
    INTEGER, INTENT(out) :: err
    IF(C_ASSOCIATED(cusp_biCGSTAB_solver_ptr,c_NULL_ptr))THEN
       err = CUSP_PTR_NOT_INITIALIZED
       RETURN
    END IF

    err = cusp_biCGSTAB_copyH2D_system_intrf(cusp_biCGSTAB_solver_ptr, Avals, xHost, bHost)

  END SUBROUTINE cusp_biCGSTAB_copyH2D_system

  SUBROUTINE cusp_biCGSTAB_copyD2H_x(xHost, err)
    IMPLICIT NONE
    DOUBLE PRECISION :: xHost(*)
    INTEGER, INTENT(out) :: err
    IF(C_ASSOCIATED(cusp_biCGSTAB_solver_ptr,c_NULL_ptr))THEN
       err = CUSP_PTR_NOT_INITIALIZED
       RETURN
    END IF

    err = cusp_biCGSTAB_copyD2H_x_intrf(cusp_biCGSTAB_solver_ptr, xHost)

  END SUBROUTINE cusp_biCGSTAB_copyD2H_x

  SUBROUTINE cusp_biCGSTAB_solveDev_system(relTol, absTol, maxItr, err)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: relTol, absTol
    INTEGER, INTENT(IN)          :: maxItr
    INTEGER, INTENT(out)         :: err
    IF(C_ASSOCIATED(cusp_biCGSTAB_solver_ptr,c_NULL_ptr))THEN
       err = CUSP_PTR_NOT_INITIALIZED
       RETURN
    END IF

    err = cusp_biCGSTAB_solveDev_system_intrf(cusp_biCGSTAB_solver_ptr, &
         relTol, absTol, maxItr)

  END SUBROUTINE cusp_biCGSTAB_solveDev_system

  SUBROUTINE cusp_biCGSTAB_getMonitor(residual, nItr, err)
    IMPLICIT NONE
    DOUBLE PRECISION    :: residual
    INTEGER             :: nItr
    INTEGER, INTENT(out):: err
    IF(C_ASSOCIATED(cusp_biCGSTAB_solver_ptr,c_NULL_ptr))THEN
       err = CUSP_PTR_NOT_INITIALIZED
       RETURN
    END IF

    err = cusp_biCGSTAB_getMonitor_intrf(cusp_biCGSTAB_solver_ptr, residual, nItr)

  END SUBROUTINE cusp_biCGSTAB_getMonitor

  SUBROUTINE cusp_biCGSTAB_shutdown(err)
    IMPLICIT NONE
    INTEGER, INTENT(out):: err
    IF(C_ASSOCIATED(cusp_biCGSTAB_solver_ptr,c_NULL_ptr))THEN
       err = CUSP_PTR_NOT_INITIALIZED
       RETURN
    END IF

    err = cusp_biCGSTAB_shutdown_intrf(cusp_biCGSTAB_solver_ptr)

  END SUBROUTINE cusp_biCGSTAB_shutdown

END MODULE modcusp_library_intrf

!---------------------------UNIT TESTS---------------------------------------------------------
!PROGRAM TESTER
!  USE modcusp_library_intrf, ONLY : getInstance_cusp_biCGSTAB_solver, cusp_biCGSTAB_initDevice, &
!       cusp_biCGSTAB_allocDevice, cusp_biCGSTAB_copyH2D_AInds, &
!       cusp_biCGSTAB_copyH2D_system, cusp_biCGSTAB_copyD2H_x, &
!       cusp_biCGSTAB_solveDev_system, cusp_biCGSTAB_getMonitor, &
!       cusp_biCGSTAB_shutdown
!
!  IMPLICIT NONE
!  INTEGER :: devID = 0
!  INTEGER :: err, n = 4, m = 10
!
!  INTEGER ,DIMENSION(10) :: rows = (/0, 0, 1, 1, 1, 2, 2, 2, 3, 3/)
!  INTEGER ,DIMENSION(10):: cols = (/0, 1, 0, 1, 2, 1, 2, 3, 2, 3/)
!  DOUBLE PRECISION ,DIMENSION(10):: vals = (/2.D0, -1.D0, -1.D0, 2.D0, -1.D0 , &
!       -1.D0, 2.D0, -1.D0, -1.D0, 2.D0/)
!
!  DOUBLE PRECISION ,DIMENSION(4) :: x = (/0.D0, 0.D0, 0.D0, 0.D0/)
!  DOUBLE PRECISION ,DIMENSION(4):: b = (/1.D0, 2.D0, 2.D0, 1.D0/)
!
!  DOUBLE PRECISION ,DIMENSION(4):: x2= (/0.D0, 0.D0, 0.D0, 0.D0/)
!
!  DOUBLE PRECISION :: relTol = 1.d-12, absTol = 0.D0, residual
!  INTEGER :: maxItr = 10000, nItr
!
!
!  CALL getInstance_cusp_biCGSTAB_solver()
!
!  CALL cusp_biCGSTAB_initDevice(devID, err)
!  IF(err == 0)THEN
!     WRITE(*,*) "initDevice completed successfully!"
!  END IF
!
!  CALL cusp_biCGSTAB_allocDevice(n, m, err)
!  IF(err == 0)THEN
!     WRITE(*,*) "initDevice completed successfully!"
!  END IF
!
!  CALL cusp_biCGSTAB_copyH2D_AInds(rows, cols, err)
!  IF(err == 0)THEN
!     WRITE(*,*) "copyH2D_AInds completed successfully!"
!  END IF
!
!  CALL cusp_biCGSTAB_copyH2D_system(vals,x,b, err)
!  IF(err == 0)THEN
!     WRITE(*,*) "copyH2D_system completed successfully!"
!  END IF
!
!  CALL cusp_biCGSTAB_solveDev_system(relTol, absTol, maxItr, err)
!
!  CALL cusp_biCGSTAB_copyD2H_x(x2, err)
!  IF(err == 0)THEN
!     WRITE(*,*) "copyD2H_x completed successfully!"
!  END IF
!
!  CALL cusp_biCGSTAB_getMonitor(residual, nItr, err)
!
!  WRITE(*,*) x2
!  WRITE(*,*) x
!
!  CALL cusp_biCGSTAB_shutdown(err)
!  IF(err == 0)THEN
!     WRITE(*,*) "cusp shutdown completed successfully!"
!  END IF
!
!
!END PROGRAM TESTER
