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

MODULE modfd_set_bc

PRIVATE
PUBLIC :: fd_bctime,fd_bcuv,fd_bcpressure,fd_bctemp,fd_bcout,fd_calc_bcuv_grad

CONTAINS
!====================================================
SUBROUTINE fd_bctime(nt)

USE precision,    ONLY : r_single
USE shared_data,  ONLY : nim,li,nj,om,time,ulid,u,flomas,flomom,densit,f1,y,r,duct,njm,ulid,ndt,ft1,cpf,ft1,&
                         ni
USE real_parameters, ONLY : zero,one,half

IMPLICIT NONE
INTEGER,INTENT(IN)   :: nt
INTEGER             :: i,ij,j
REAL(KIND = r_single):: fact,dummy

IF(ndt>0)THEN
  fact = tanh(REAL(nt,r_single)/REAL(ndt,r_single))
ELSE
  fact=one
ENDIF

IF(duct)THEN
  flomas=zero
  flomom=zero
  DO j=2,njm
    ij=li(1)+j 
    u(ij)=ulid
    f1(ij)=half*densit*(y(j)-y(j-1))*(r(j)+r(j-1))*u(ij)
    ft1(ij)=half*cpf*densit*(y(j)-y(j-1))*(r(j)+r(j-1))*u(ij)
    flomas=flomas+f1(ij)
    flomom=flomom+f1(ij)*u(ij)
  END DO
ELSE
  DO i=2,nim
    ij=li(i)+nj
    u(ij)=ulid !*sin(om*time)
  END DO
  DO i=2,nim
    ij=li(i)+1
    u(ij)=-ulid !*sin(om*time)
  END DO
ENDIF

END SUBROUTINE fd_bctime

SUBROUTINE fd_bcpressure(phi)

!--This routine calculates boundary values of pressure or
!--pressure-correction by extrapolating (linearly) from inside.

USE shared_data,  ONLY : li,fy,fx,nim,nj,njm,ni,yPeriodic,xPeriodic
USE precision,    ONLY : r_single

IMPLICIT NONE 
REAL(KIND = r_single),DIMENSION(:),INTENT(INOUT) :: phi
INTEGER :: i,j,ij,nj2

!--SOUTH AND NORTH BOUNDARIES
IF(yPeriodic == 0)THEN
  DO i=2,nim
    ij=li(i)+1
    phi(ij)=phi(ij+1)+(phi(ij+1)-phi(ij+2))*fy(2)
    ij=li(i)+nj
    phi(ij)=phi(ij-1)+(phi(ij-1)-phi(ij-2))*(1.-fy(njm-1)) 
  END DO 
ENDIF
!--WEST AND EAST BOUNDARIES
IF(xPeriodic == 0)THEN
  nj2=2*nj
  DO j=2,njm
    ij=li(1)+j
    phi(ij)=phi(ij+nj)+(phi(ij+nj)-phi(ij+nj2))*fx(2) 
    ij=li(ni)+j
    phi(ij)=phi(ij-nj)+(phi(ij-nj)-phi(ij-nj2))*(1.-fx(nim-1))
  END DO
ENDIF

END SUBROUTINE fd_bcpressure

SUBROUTINE fd_bcuv

!--In this routine, boundary conditions for U and V equations
!--are implemented, i.e. fluxes through boundary cell faces
!--are approximated. Here, the boundaries encountered in 
!--cavity flows are considered

USE shared_data,   ONLY : r,x,y,xc,yc,su,apu,apv,sv,visc,li,nim,njm,u,v,ni,nj,&
                          f1,duct,ae,movingmesh,lamvisc,xPeriodic,yPeriodic
USE precision,     ONLY : r_single
USE real_parameters,ONLY : zero,half 
IMPLICIT NONE

REAL(KIND = r_single) :: d,awc
INTEGER               :: i,j,ij


IF(yPeriodic == 0)THEN
!  IF(duct)THEN
!  !--Use slip walls
!    DO i=2,nim
!      ij=li(i)+2
!      d=lamvisc(ij-1)*(x(i)-x(i-1))*r(1)/(yc(2)-yc(1))
!      apv(ij)=apv(ij)+d
!    END DO
!
!    DO i=2,nim
!      ij=li(i)+njm
!      d=lamvisc(ij+1)*(x(i)-x(i-1))*r(njm)/(yc(nj)-yc(njm))
!      apv(ij)=apv(ij)+d
!    END DO
!
!  ELSE
    DO i=2,nim
      ij=li(i)+2
      d=lamvisc(ij-1)*(x(i)-x(i-1))*r(1)/(yc(2)-yc(1))
      apu(ij)=apu(ij)+d
      su(ij) =su(ij) +d*u(ij-1)
    END DO

    DO i=2,nim
      ij=li(i)+njm
      d=lamvisc(ij+1)*(x(i)-x(i-1))*r(njm)/(yc(nj)-yc(njm))
      apu(ij)=apu(ij)+d
      su(ij) =su(ij) +d*u(ij+1)
    END DO
 ! ENDIF
ENDIF
!--west
IF(xPeriodic == 0)THEN
  IF(duct)THEN
    DO j=2,njm
      ij=li(2)+j
      d=half*lamvisc(ij-nj)*(y(j)-y(j-1))*(r(j)+r(j-1))/(xc(2)-xc(1))
      awc=d+f1(ij-nj)
      apu(ij)=apu(ij)+awc
      apv(ij)=apv(ij)+awc
      su(ij) =su(ij) +awc*u(ij-nj)
      sv(ij) =sv(ij) +awc*v(ij-nj)
    ENDDO
  ELSE
    DO j=2,njm
      ij=li(2)+j
      d=half*lamvisc(ij-nj)*(y(j)-y(j-1))*(r(j)+r(j-1))/(xc(2)-xc(1))
      apv(ij)=apv(ij)+d
      sv(ij) =sv(ij) +d*v(ij-nj)
    END DO 
  ENDIF

  !--east
  IF(duct)THEN !--outlet
    DO j=2,njm
      ij=li(nim)+j
      ae(ij)=zero
    END DO
  ELSE
    IF(movingmesh)THEN !--Assume outlet
      DO j=2,njm
        ij=li(nim)+j
        ae(ij)=zero
      END DO
    ELSE
      DO j=2,njm
        ij=li(nim)+j
        d=half*lamvisc(ij+nj)*(y(j)-y(j-1))*(r(j)+r(j-1))/(xc(ni)-xc(nim))
        apv(ij)=apv(ij)+d
        sv(ij) =sv(ij) +d*v(ij+nj)
      END DO
    ENDIF
  ENDIF
ENDIF
END SUBROUTINE fd_bcuv

SUBROUTINE fd_bctemp

!--In this routine, boundary conditions for the temperature 
!--equation are implemented, i.e. heat fluxes through the
!--boundary cell faces are calculated. Here, specified wall 
!--temperature and adiabatic wall (zero heat flux) are considered;
!--treatment at symmetry planes is the same as for an adiabatic
!--wall
USE real_parameters,ONLY: half,zero
USE shared_data,   ONLY : ni,nim,li,t,nj,njm,visc,celkappa,y,yc,r,ap,su,xc,x,duct,ae,movingmesh,&
                          yPeriodic,xPeriodic,f1,celcp,ft1
USE precision,     ONLY : r_single

IMPLICIT NONE

INTEGER :: i,ij,j
REAL(KIND = r_single) :: d,awc

IF(yPeriodic == 0)THEN
  IF(movingmesh)THEN !--ISOTHERMAL WALL
    !--NORTH BOUNDARY (ISOTHERMAL WALL, NON-ZERO DIFFUSIVE FLUX)
    DO i=2,nim
      ij=li(i)+njm
      d=celkappa(ij)*(x(i)-x(i-1))*r(njm)/(yc(nj)-yc(njm))
      ap(ij)=ap(ij)+d
      su(ij)=su(ij)+d*t(ij+1)
    ENDDO

    !--SOUTH BOUNDARY (ISOTHERMAL WALL, NON-ZERO DIFFUSIVE FLUX)
    DO i=2,nim
      ij=li(i)+2
      d=celkappa(ij)*(x(i)-x(i-1))*r(1)/(yc(2)-yc(1))
      ap(ij)=ap(ij)+d
      su(ij)=su(ij)+d*t(ij-1)
    ENDDO

  ELSE

    !--SOUTH BOUNDARY (ADIABATIC WALL, DT/DY=0, ZERO FLUX)
    DO i=2,nim
      ij=li(i)+1
      t(ij)=t(ij+1)
    END DO

    !--NORTH BOUNDARY (ADIABATIC WALL, DT/DY=0, ZERO FLUX)
    DO i=2,nim
      ij=li(i)+nj
      t(ij)=t(ij-1)
    END DO

  ENDIF
ENDIF

IF(xPeriodic == 0)THEN
  IF(duct)THEN !--Inlet
    DO j=2,njm
      ij=li(2)+j
      d=half*celkappa(ij-nj)*(y(j)-y(j-1))*(r(j)+r(j-1))/(xc(2)-xc(1))
      awc=d+ft1(ij-nj)
      ap(ij)=ap(ij)+awc
      su(ij) =su(ij) +awc*t(ij-nj)
    ENDDO
  ELSE
    !--WEST BOUNDARY (ISOTHERMAL WALL, NON-ZERO DIFFUSIVE FLUX)
    DO j=2,njm
      ij=li(2)+j
      d=half*celkappa(ij)*(y(j)-y(j-1))*(r(j)+r(j-1))/(xc(2)-xc(1))
      ap(ij)=ap(ij)+d
      su(ij)=su(ij)+d*t(ij-nj)
    ENDDO
  !    !--West BOUNDARY (ADIABATIC WALL, DT/DX=0, ZERO FLUX)
  !  DO j=2,njm
  !    t(j)=t(j+nj)
  !  END DO
  ENDIF

  IF(duct)THEN !--Outlet
    DO j=2,njm
      ij=li(nim)+j
      t(ij+nj) = t(ij) !-Set to calculate the grads
      ae(ij)=zero
    END DO
  ELSE
    IF(movingmesh)THEN !--Outlet
      DO j=2,njm
        ij=li(nim)+j
        ae(ij)=zero
      END DO
    ELSE
      !--EAST BOUNDARY (ISOTHERMAL WALL)
      DO j=2,njm
        ij=li(nim)+j
        d=half*celkappa(ij)*(y(j)-y(j-1))*(r(j)+r(j-1))/(xc(ni)-xc(nim))
        ap(ij)=ap(ij)+d
        su(ij)=su(ij)+d*t(ij+nj)
      ENDDO
  !    !--EAST BOUNDARY (ADIABATIC WALL, DT/DX=0, ZERO FLUX)
  !    DO j=2,njm
  !      ij = li(ni) + j
  !      t(ij)=t(ij-nj)
  !    END DO
    ENDIF
  ENDIF

ENDIF

END SUBROUTINE fd_bctemp 

SUBROUTINE fd_bcout

USE precision,          ONLY : r_single
USE real_parameters,    ONLY : zero,tiny,half,one
USE shared_data,        ONLY : njm,nim,li,densit,y,u,v,r,f1,flomas,nj,movingmesh

IMPLICIT NONE 

REAL(KIND = r_single) :: fmout,fac
INTEGER               :: j,ij

fmout=zero

DO j=2,njm
  ij=li(nim)+j
  f1(ij)=half*densit*(y(j)-y(j-1))*(r(j)+r(j-1))*u(ij)
  fmout=fmout+f1(ij)
END DO
!CALCULATE SCALING FACTOR TO MAKE MASS FLUX EQUAL INCOMING ONE
IF(movingmesh)THEN
  fac = one
ELSE
  fac=flomas/(fmout+tiny)
ENDIF

!CORRECT VELOCITY AND MASS FLUXES AT OUTLET

DO j=2,njm
  ij=li(nim)+j
  f1(ij)=f1(ij)*fac
  u(ij+nj)=u(ij)*fac
  v(ij+nj)=v(ij) !--to get the gradients right
END DO

END SUBROUTINE fd_bcout

SUBROUTINE fd_calc_bcuv_grad

USE shared_data,        ONLY : ni,nim,nj,njm,yc,y,xc,x,u,v,dux,duy,dvx,dvy,fy,fx,li,duct,movingmesh
USE real_parameters,    ONLY : zero,one

IMPLICIT NONE

INTEGER     :: i,j,ij

IF(duct)THEN
!--Use slip walls dU/dY = 0, dV/dX = 0
  DO i=2,nim
    ij=li(i)+1
    duy(ij) = zero
    dvx(ij) = zero
    dux(ij) = dux(ij + 1)
    dvy(ij) = (v(ij + 1) - v(ij))/(yc(2)-yc(1))
  END DO
!
  DO i=2,nim
    ij=li(i)+nj
    duy(ij) = zero
    dvx(ij) = zero
    dux(ij) = dux(ij - 1)
    dvy(ij) = (v(ij) - v(ij - 1))/(yc(nj)-yc(njm))
  END DO
!
ELSE
  DO i=2,nim
    ij=li(i)+1
    dux(ij) = zero
    dvx(ij) = zero
    dvy(ij) = zero
    duy(ij) = (u(ij + 1) - u(ij))/(yc(2)-yc(1))
  END DO
!
  DO i=2,nim
    ij=li(i)+nj
    dux(ij) = zero
    dvx(ij) = zero
    dvy(ij) = zero
    duy(ij) = (u(ij) - u(ij - 1 ))/(yc(nj)-yc(njm))
  END DO
ENDIF

!--west
IF(duct)THEN
  DO j=2,njm
    dux(j) = (u(j + nj) - u(j)) / (xc(2)-xc(1))
    dvx(j) = (v(j + nj) - v(j)) / (xc(2)-xc(1))
    !--Assume uniform inlet !--Shoule be calculated otherwise
    duy(j) = zero
    dvy(j) = zero
  ENDDO
ELSE
  DO j=2,njm
    dux(j) = zero
    duy(j) = zero
    dvy(j) = zero
    dvx(j) = (v(j+nj) - v(j))/(xc(2)-xc(1))
  END DO 
ENDIF

!
!!--east
IF(duct)THEN !--outlet
  DO j=2,njm
    ij=li(ni)+j
    dux(ij) = zero
    dvx(ij) = zero
    duy(ij) = duy(ij - nj)
    dvy(ij) = dvy(ij - nj)  
  END DO
ELSE
  IF(movingmesh)THEN !--Assume outlet
    DO j=2,njm
      ij=li(ni)+j
      dux(ij) = zero
      dvx(ij) = zero
      duy(ij) = duy(ij - nj)
      dvy(ij) = dvy(ij - nj)  
    END DO
  ELSE
    DO j=2,njm
      ij=li(ni)+j
      dux(j) = zero
      duy(j) = zero
      dvy(j) = zero
      dvx(j) = (v(j) - v(j-nj))/(xc(ni)-xc(nim))
    END DO
  ENDIF
ENDIF

END SUBROUTINE fd_calc_bcuv_grad

END MODULE modfd_set_bc
