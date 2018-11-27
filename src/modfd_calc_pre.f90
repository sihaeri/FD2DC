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

MODULE modfd_calc_pre

PRIVATE
PUBLIC :: fd_calc_pre
CONTAINS
!============================================
SUBROUTINE fd_calc_pre

USE parameters,         ONLY : out_unit,solver_sparsekit,solver_sip,solver_cg,solver_hypre,OPSUCCESS,SOLVER_DONE,&
                               use_GPU_yes,SOLVER_FAILED
USE real_parameters,    ONLY : one,zero,half,two,four
USE precision,          ONLY : r_single
USE modfd_set_bc,       ONLY : fd_bcpressure
USE shared_data,        ONLY : urf,ip,p,su,apu,apv,nim,njm,&
                               xc,ni,nj,li,fx,y,r,visc,su,&
                               u,v,ae,aw,an,as,fy,yc,x,lcal,ien,&
                               den,deno,denoo,laxis,p,nij,xPeriodic,yPeriodic,&
                               sor,resor,nsw,f1,f2,ft1,ft2,celcp,dpx,dpy,&
                               ltime,ap,ipr,jpr,ltest,pp,fdsu,fdsv,&
                               dux,duy,dvx,dvy,fdsuc,fdsvc,fdsub,fdsvb,putobj,&
                               Ncel,NNZ,Acoo,Arow,Acol,Acsr,Aclc,Arwc, &
                               solver_type,rhs,sol,work,alu,jlu,ju,jw,ct1,ct2,ct3,&
                               Hypre_A,Hypre_b,Hypre_x,mpi_comm,lli,celcp,use_GPU
!handle, config, platformOpts, formatOpts, solverOpts, precondOpts, culaStat,res

!use cula_sparse_type
!use cula_sparse

USE  modfd_solve_linearsys,  ONLY : fd_solve_sip2d,fd_spkit_interface,copy_solution,calc_residual,&
                                    fd_solve_cgs2d,fd_cooVal_create
USE modcusp_library_intrf,   ONLY : cusp_biCGSTAB_copyH2D_system,cusp_BiCGSTAB_solveDev_system , &
     cusp_biCGSTAB_getMonitor, cusp_biCGSTAB_copyD2H_x
IMPLICIT NONE
REAL(KIND = r_single) :: fxe,fxp,dxpe,s,d,fyn,fyp,dypn,&
                         dx,dy,rp,ppe,ppw,ppn,pps,dpxel,&
                         uel,apue,ue,vole,voln,dpynl,vnl,apvn,&
                         vn,sum,ppo,dpxe,dpyn,pcor,aet(nij),ant(nij),difft !,&
                         !duxel,duxe,dvynl,dvyn,fdsuce,fdsvcn,dfyn,dfxe,c,K

INTEGER               :: ij,i,j,ije,ijn,ijw,ijs,ijpref,ipar(16),debug,error,xend,yend
REAL                  :: fpar(16)
DOUBLE PRECISION      :: absTol=0.D0

debug = 0
!--EAST CV FACES (S - AREA, VOLE - VOLUME BETWEEN P AND E)
xend = nim+xPeriodic-1 !--Apply periodic conditions
DO i=2,xend

  dxpe=xc(i+1)-xc(i)
  fxe=fx(i)
  fxp=one-fxe

  DO j=2,njm
    ij=li(i)+j
    ije=ij+nj-i/nim*((i-1)*nj)

    s=(y(j)-y(j-1))*(r(j)+r(j-1))*half
    vole=dxpe*s
    d=(den(ije)*fxe+den(ij)*fxp)*s
    difft=(celcp(ije)*fxe+celcp(ij)*fxp)*s ! celcp contains both cp and den
    !--INTERPOLATED CELL FACE QUANTITIES (PRESSURE GRAD., U AND 1/AP)
    !--Note: pressure gradient is interpolated midway between P and E,
    !--since the gradient calculated at cell face is second order
    !--accurate at that location; the velocity is interpolated linearly,
    !--to achieve second order accuracy at cell face center.

    dpxel=half*(dpx(ije)+dpx(ij))
    uel=u(ije)*fxe+u(ij)*fxp
    apue=apu(ije)*fxe+apu(ij)*fxp

    !--CELL FACE GRADIENT, VELOCITY AND MASS FLUX
    dpxe=(p(ije)-p(ij))/dxpe

    pcor=apue*vole*(dpxe-dpxel)

    ue=uel - pcor

    f1(ij) = d*ue
    ft1(ij)= difft*ue
    !--COEFFICIENTS OF P' EQUATION, AE(P) AND AW(E)

    ae(ij)=-d*apue*s
    aet(ij)=-difft*apue*s
    aw(ije)=ae(ij)

  END DO
END DO

!--NORTH CV FACES (S - AREA, VOLN - VOLUME BETWEEN P AND N)
yend = njm+yPeriodic-1 !--Apply periodic conditions
DO j=2,yend
  dypn=yc(j+1)-yc(j)
  fyn=fy(j)
  fyp=one-fyn

  DO i=2,nim
    ij=li(i)+j
    ijn=ij+1-j/njm*(j-1)
    !ijs=ij-1

    s=(x(i)-x(i-1))*r(j)
    voln=s*dypn
    d=(den(ijn)*fyn+den(ij)*fyp)*s
    difft=(celcp(ijn)*fyn+celcp(ij)*fyp)*s
    !--INTERPOLATED CELL-FACE QUANTITIES (PRESSURE GRAD., U AND 1/AP)
    dpynl=half*(dpy(ijn)+dpy(ij))
    vnl=v(ijn)*fyn+v(ij)*fyp
    apvn=apv(ijn)*fyn+apv(ij)*fyp

    !dvynl=half*(dvy(ijn)+dvy(ij))
    !--CELL-FACE GRADIENT, VELOCITY AND MASS FLUX

    dpyn=(p(ijn)-p(ij))/dypn
    !dvyn=(v(ijn)-v(ij))/dypn
    !dfyn=(fdsv(ijn)-fdsv(ij))/dypn
    pcor=apvn*voln*(dpyn-dpynl)

    vn=vnl - pcor

    f2(ij)=d*vn
    ft2(ij)=difft*vn
    !--COEFFICIENTS OF P' EQUATION, AN(P) AND AS(N)

    an(ij)=-d*apvn*s
    ant(ij)=-difft*apvn*s
    as(ijn)=an(ij)

  END DO
END DO


!--BOUNDARY CONDITIONS: PRESCRIBED MASS FLUXES, ZERO CORRECTION
!--(EQUIVALENT TO ZERO NORMAL GRADIENT FOR P'; COEFFICIENT FOR
!--THE BOUNDARY NODE IS ZERO, NO SPECIAL TREATMENT REQUIRED)


!--SORCE TERM AND COEFFICIENT OF NODE P

sum=zero
DO i=2,nim
  DO j=2,njm
    ij=li(i)+j
    ijs=ij-1+yPeriodic*Max(3-j,0)*(njm-1)
    ijw=ij-nj+xPeriodic*Max(3-i,0)*(nim-1)*nj
    su(ij)=f1(ijw)-f1(ij)+f2(ijs)-f2(ij)
    IF(putobj)THEN
      su(ij) = su(ij) + (ct3*denoo(ij) + ct2*deno(ij) - ct1*den(ij)) *  ((x(i)-x(i-1))*(y(j)-y(j-1))*half*(r(j)+r(j-1)))
    ENDIF
    ap(ij)=-(ae(ij)+aw(ij)+an(ij)+as(ij))
    sum=sum+su(ij)
    pp(ij)=zero
  END DO
END DO

!--SUM MUST BE ZERO IF GLOBAL MASS CONSERVATION IS ASSURED!
IF(ltest) WRITE(2,*) '       SUM = ',SUM

!--SOLVE EQUATIONS SYSTEM FOR P' AND APPLY CORRECTIONS
IF(solver_type == solver_sparsekit)THEN

   SOL = 0.0
  IF(use_GPU == use_GPU_yes)THEN

    CALL fd_cooVal_create(ap,as,an,aw,ae,su,pp)

!    config%relativeTolerance = sor(ip)
!    config%maxIterations = nsw(ip)
!    culaStat = culaSparseCudaDcooBicgstabJacobi(handle, config, platformOpts, formatOpts, &
!                                solverOpts, precondOpts, NCel, NNZ, Acoo, Arow, Acol, sol, rhs, res)

    CALL cusp_biCGSTAB_copyH2D_system(Acoo, SOL, RHS, error)
    IF(error /= OPSUCCESS)GOTO 100

    CALL cusp_BiCGSTAB_solveDev_system(sor(ip),absTol,nsw(ip),error)
    IF(error /= OPSUCCESS)GOTO 100

    CALL cusp_BiCGSTAB_getMonitor(resor(ip),ipar(1),error)
    IF(error /= OPSUCCESS)GOTO 100

    CALL cusp_biCGSTAB_copyD2H_x(sol,error)
    IF(error /= OPSUCCESS)GOTO 100

!    IF(culaStat /= culaSparseNoError)THEN
!      WRITE(*,*)'CULA encoutered an error.'
!      culaStat=culaSparseDestroy(handle)
!      STOP
!    ENDIF

    CALL copy_solution(pp,sol)
    100 CONTINUE
    IF(error /= OPSUCCESS)THEN
      WRITE(*,*)'GPU--Pressure solver problem in CUDA operations.'
      STOP
      STOP
    ENDIF
 ELSE

      CALL fd_cooVal_create(ap,as,an,aw,ae,su,pp)
      CALL coocsr(Ncel,NNZ,Acoo,Arow,Acol,Acsr,Aclc,Arwc)

      ipar  =  0
      fpar  = 0.0

      !--preconditioner
      !ipar(2) = 1     ! 1=left, 2=right
      ipar(2) = 0           ! 0 == no preconditioning
      ipar(3) = 0
      ipar(4) = NNZ*8       ! workspace for BGStab
      ipar(6) = nsw(ip)     ! maximum number of matrix-vector multiplies

      fpar(1) = sor(ip)  ! relative tolerance, must be between (0, 1)
      fpar(2) = 0.0      ! absolute tolerance, must be positive

          !call solve_spars(debug,IOdbg,Ncel,RHS,SOL,ipar,fpar,&
          !                       WORK,Acsr,Aclc,Arwc)
      CALL solve_spars(debug,out_unit,Ncel,RHS,SOL,ipar,fpar,&
                            WORK, Acsr,Aclc,Arwc,NNZ,Alu,Jlu,Ju,Jw)

     CALL copy_solution(pp,sol)

     CALL calc_residual(ap,as,an,aw,ae,su,pp,resor(ip))

   ENDIF


ELSEIF(solver_type == solver_sip)THEN

  CALL fd_solve_sip2d(ae,an,aw,as,ap,su,pp,li,ni,nj,nsw(ip),resor(ip),sor(ip))

ELSEIF(solver_type == solver_cg)THEN

  CALL fd_solve_cgs2d(ae,an,aw,as,ap,su,pp,li,ni,nj,nsw(ip),resor(ip),sor(ip))

ELSEIF(solver_type == solver_hypre)THEN
  !CALL solve_hypre(mpi_comm,li,lli,ncel,ni,nj,ap,as,an,aw,ae,su,pp,sor(ip),resor(ip),&
  !                  Hypre_A,Hypre_b,Hypre_x,1,nsw(ip))

ENDIF

!--CALCULATE PRESSURE CORRECTION AT BOUNDARIES

CALL fd_bcpressure(pp)

!--VALUE OF P' AT REFERENCE LOCATION TO BE SUBTRACTED FROM ALL P'

ijpref=li(ipr)+jpr
ppo=pp(ijpref)

!--CORRECT EAST MASS FLUXES
DO i=2,xend
  DO j=2,njm
    ij = li(i) + j
    ije=ij+nj-i/nim*((i-1)*nj)
    f1(ij)=f1(ij)+ae(ij)*(pp(ije)-pp(ij))
    ft1(ij)=ft1(ij)+aet(ij)*(pp(ije)-pp(ij))
  END DO
END DO

!--CORRECT NORTH MASS FLUXES

DO i=2,nim
  DO j=2,yend
    ij = li(i) + j
    ijn=ij+1-j/njm*(j-1)
    f2(ij)=f2(ij)+an(ij)*(pp(ijn)-pp(ij))
    ft2(ij)=ft2(ij)+ant(ij)*(pp(ijn)-pp(ij))
  END DO
END DO

sum=zero
DO i=2,nim
  DO j=2,njm
    ij=li(i)+j
    ijs=ij-1+yPeriodic*Max(3-j,0)*(njm-1)
    ijw=ij-nj+xPeriodic*Max(3-i,0)*(nim-1)*nj
    su(ij)=f1(ijw)-f1(ij)+f2(ijs)-f2(ij)
    IF(putobj)THEN
      su(ij) = su(ij) +  (ct3*denoo(ij) + ct2*deno(ij) - ct1*den(ij)) *  ((x(i)-x(i-1))*(y(j)-y(j-1))*half*(r(j)+r(j-1)))
    ENDIF
    sum=sum+ABS(su(ij))
  END DO
END DO
resor(ip) = sum
!--CORRECT PRESSURE AND VELOCITIES AT CELL CENTER
DO i=2,nim
  dx=x(i)-x(i-1)

  DO j=2,njm
    ij=li(i)+j
    rp=half*(r(j)+r(j-1))
    dy=y(j)-y(j-1)

    ije=ij+nj-xPeriodic*i/nim*((i-1)*nj)
    ijw=ij-nj+xPeriodic*Max(3-i,0)*(nim-1)*nj
    ijn=ij+1-yPeriodic*j/njm*(j-1)
    ijs=ij-1+yPeriodic*Max(3-j,0)*(njm-1)

    ppe=pp(ije)*fx(i)+pp(ij)*(one-fx(i))
    ppw=pp(ij)*fx(i-1)+pp(ijw)*(one-fx(i-1))
    ppn=pp(ijn)*fy(j)+pp(ij)*(one-fy(j))
    pps=pp(ij)*fy(j-1)+pp(ijs)*(one-fy(j-1))

    u(ij)=u(ij)-(ppe-ppw)*dy*rp*apu(ij)
    v(ij)=v(ij)-(ppn-pps)*dx*rp*apv(ij)
    p(ij)=p(ij)+urf(ip)*(pp(ij)-ppo)

  END DO
END DO

END SUBROUTINE fd_calc_pre

END MODULE modfd_calc_pre
