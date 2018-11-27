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

MODULE modfd_calc_mom

PRIVATE
PUBLIC :: fd_calc_mom
CONTAINS
!==================================================
SUBROUTINE fd_calc_mom

!--This routine sets the coefficient matrix for the U and
!--V equations, and calls the linear equation solver to
!--update the velocity components. Constant fluid
!--properties are assumed (parts of diffusive fluxes
!--cancel out)
USE parameters,         ONLY : out_unit,solver_sparsekit,solver_sip,solver_cg,solver_hypre,use_GPU_yes,OPSUCCESS,SOLVER_DONE,&
                               SOLVER_FAILED
USE real_parameters,    ONLY : one,zero,half,small,two,betam
USE shared_data,        ONLY : urf,iu,iv,p,su,sv,apu,apv,nim,njm,&
                               xc,ni,nj,li,fx,y,r,su,sv,gds,&
                               u,v,ae,aw,an,as,fy,yc,x,lcal,ien,&
                               gravx,gravy,beta,den,deno,denoo,laxis,p,ct1,ct2,ct3,&
                               vo,uo,voo,uoo,sor,resor,nsw,f1,f2,&
                               dpx,dpy,t,tref,ltime,ap,fdsu,fdsv,putobj,&
                               dvx,dvy,dux,duy,densit,temp_visc,lamvisc,&
                               Ncel,NNZ,Acoo,Arow,Acol,Acsr,Aclc,Arwc, &
                               solver_type,rhs,sol,work,alu,jlu,ju,jw,&
                               Hypre_A,Hypre_b,Hypre_x,mpi_comm,lli,celbeta,&
                               fdfcu,fdfcv,use_GPU,xPeriodic,yPeriodic
!handle, config, platformOpts, formatOpts, solverOpts, precondOpts, res, culaStat

!use cula_sparse_type
!use cula_sparse

USE modfd_set_bc,       ONLY : fd_bcpressure,fd_bcuv,fd_calc_bcuv_grad
USE precision,          ONLY : r_single
USE modfd_solve_linearsys,     ONLY : fd_solve_sip2d,fd_spkit_interface,copy_solution,calc_residual,fd_cooVal_create
USE modcusp_library_intrf,     ONLY : cusp_biCGSTAB_copyH2D_system,cusp_BiCGSTAB_solveDev_system , &
     cusp_biCGSTAB_getMonitor, cusp_biCGSTAB_copyD2H_x

IMPLICIT NONE

REAL(KIND = r_single) :: urfu,urfv,fxe,fxp,dxpe,s,d,cp,ce,&
                         fuuds,fvuds,fucds,fvcds,fyn,fyp,dypn,&
                         dx,dy,rp,vol,sb,pe,pw,pn,ps,apt,apto,aptoo,cn,&
                         vele,velw,veln,vels,due,dve,dun,dvn,visi,&
                         duv(2),dduv(2),uvct,uc,ud,vc,vd,ufcd,vfcd,uf,vf,fg,graduv,gradcf
INTEGER               :: ij,i,j,ije,ijn,ijs,ijw,ipar(16),debug,error,xend,yend,ic,ii,jj
REAL                  :: fpar(16)
DOUBLE PRECISION      :: absTol=0.D0

debug = 0

!--RECIPROCAL VALUES OF UNDER-RELAXATION FACTORS FOR U AND V

urfu=one/urf(iu)
urfv=one/urf(iv)

!--SET BOUNDARY PRESSURE (LINEAR EXTRAPOLATION FROM INSIDE)

CALL fd_bcpressure(p)

!--INITIALIZE TEMPORARILY STORED VARIABLES
su=zero
sv=zero
apu=zero
apv=zero

!--FLUXES THROUGH INTERNAL EAST CV FACES

!--F1(IJ) is the mass flux through the east face (outward
!--normal directed to E); FX(I) is the ratio of distance
!--from P to cell face, to distance from P to E; IJ
!--denotes node P and IJE node E.
!--Contribution of convective and diffusive fluxes from
!--east face to AE(P), AW(E), and source terms at both
!--P and E are calculated; contributions to AP(P) and
!--AP(E) come through the sum of neighbor coefficients
!--and are not explicitly calculated.
xend = nim+xPeriodic-1 !--Apply periodic conditions
DO i=2,xend

!--INTERPOLATION FACTORS, DISTANCE FROM P TO E (SAME FOR ALL J)
  fxe =fx(i)
  fxp =one-fxe
  dxpe=xc(i+1)-xc(i)

  DO j=2,njm
    ij=li(i)+j
    ije=ij+nj-i/nim*((i-1)*nj)

    !--CELL FACE AREA S = DY*RE*1

    s=(y(j)-y(j-1))*(r(j)+r(j-1))*half

    !--COEFFICIENT RESULTING FROM DIFFUSIVE FLUX (SAME FOR U AND V)
    visi = lamvisc(ije)*fxe+lamvisc(ij)*fxp
    d=visi*s/dxpe

    !--EXPLICIT CONVECTIVE FLUXES FOR UDS AND CDS
    ce=MIN(f1(ij),zero)
    cp=MAX(f1(ij),zero)

    fuuds=cp*u(ij)+ce*u(ije)
    fvuds=cp*v(ij)+ce*v(ije)

    IF( f1(ij) >= zero )THEN

      ic    = ij
      ii    = i
      uc    = u(ij)
      vc    = v(ij)

      ud    = u(ije)
      vd    = v(ije)

      duv(1) = ud-uc
      duv(2) = vd-vc

      graduv = DOT_PRODUCT(duv,duv)

      dduv(1) = two * dux(ic) * dxpe
      dduv(2) = two * dvx(ic) * dxpe

      gradCf   = DOT_PRODUCT(duv,dduv)

      uvct    = one - graduv/(gradCf+small)

    ELSE

      ic    = ije
      ii    = i+1
      uc    = u(ije)
      vc    = v(ije)

      ud    = u(ij)
      vd    = v(ij)

      duv(1) = ud-uc
      duv(2) = vd-vc

      graduv = DOT_PRODUCT(duv,duv)

      dduv(1) = -two * dux(ic) * dxpe
      dduv(2) = -two * dvx(ic) * dxpe

      gradCf   = DOT_PRODUCT(duv,dduv)

      uvCt    = one - graduv/(gradCf+small)

    ENDIF

    ufcd = u(ije)*fxe+u(ij)*fxp
    vfcd = v(ije)*fxe+v(ij)*fxp
    !--Gamma
    IF( uvct <= zero .OR. uvct >= one )THEN        ! 0.1 <= Beta_m <= 0.5 upwind
       uf = uc
       vf = vc
    ELSEIF( betam <= uvct .AND. uvct < one  )THEN  ! Central
       uf = ufcd
       vf = vfcd
    ELSEIF(  zero < uvct .AND. uvct < betam )THEN  ! Gamma
       fg = uvct / betam
       uf = fg*ufcd + (one-fg)*uc
       vf = fg*vfcd + (one-fg)*vc
    ENDIF
    !--Minmod
!    IF( uvct <= zero .OR. uvct >= one )THEN        ! 0.1 <= Beta_m <= 0.5 upwind
!      uf = uc
!      vf = vc
!    ELSEIF(  uvct < half )THEN  ! Second order Upwind
!      uf = uc + dux(ic)*(x(i) - xc(ii))
!      vf = vc + dvx(ic)*(x(i) - xc(ii))
!    ELSE  ! Central
!      uf = ufcd
!      vf = vfcd
!    ENDIF

    fucds=f1(ij)*uf
    fvcds=f1(ij)*vf

    !--COEFFICIENTS AE(P) AND AW(E) DUE TO UDS

    ae(ij) = ce-d
    aw(ije)=-cp-d

    !--SOURCE TERM CONTRIBUTIONS AT P AND E DUE TO DEFERRED CORRECTION

    su(ij) =su(ij) +gds(iu)*(fuuds-fucds)
    su(ije)=su(ije)-gds(iu)*(fuuds-fucds)
    sv(ij) =sv(ij) +gds(iu)*(fvuds-fvcds)
    sv(ije)=sv(ije)-gds(iu)*(fvuds-fvcds)

    IF(temp_visc)THEN
      due = visi*s*(dux(ije)*fxe+dux(ij)*fxp)
      su(ij)  = su(ij) + due
      su(ije) = su(ije)- due
      dve = visi*s*(duy(ije)*fxe+duy(ij)*fxp)
      sv(ij)  = sv(ij) + dve
      sv(ije) = sv(ije)- dve
    ENDIF

  END DO
END DO

!--FLUXES THROUGH INTERNAL NORTH CV FACES
!--F2(IJ) is the mass flux through the north face (outward
!--normal directed to N); FY(J) is the ratio of distance
!--from P to cell face, to distance from P to N; IJ
!--denotes node P and IJN node N.
!--Contribution of convective and diffusive fluxes from
!--north face to AN(P), AS(N), and source terms at both
!--P and N are calculated; contributions to AP(P) and
!--AP(N) come through the sum of neighbor coefficients
!--and are not explicitly calculated.
yend = njm+yPeriodic-1 !--Apply periodic conditions
DO j=2,yend

!--INTERPOLATION FACTORS, DISTANCE FROM P TO N (SAME FOR ALL J)
  fyn =fy(j)
  fyp =one-fyn
  dypn=yc(j+1)-yc(j)

  DO i=2,nim
    ij =li(i)+j
    ijn=ij+1-j/njm*(j-1)

    !--CELL FACE AREA S = DX*RN*1
    s=(x(i)-x(i-1))*r(j)
    !--COEFFICIENT RESULTING FROM DIFFUSIVE FLUX (SAME FOR U AND V)
    visi = lamvisc(ijn)*fyn+lamvisc(ij)*fyp
    d=visi*s/dypn

    !--EXPLICIT CONVECTIVE FLUXES FOR UDS AND CDS
    cn=MIN(f2(ij),zero)
    cp=MAX(f2(ij),zero)

    fuuds=cp*u(ij)+cn*u(ijn)
    fvuds=cp*v(ij)+cn*v(ijn)

    IF( f2(ij) >= zero )THEN

      ic    = ij
      jj    = j
      uc    = u(ij)
      vc    = v(ij)

      ud    = u(ijn)
      vd    = v(ijn)

      duv(1) = ud-uc
      duv(2) = vd-vc

      graduv = dot_product(duv,duv)

      dduv(1) = two * duy(ic)*dypn
      dduv(2) = two * dvy(ic)*dypn

      gradCf   = dot_product(duv,dduv)

      uvct    = one - graduv/(gradCf+small)

    ELSE

      ic    = ijn
      jj    = j+1
      uc    = u(ijn)
      vc    = v(ijn)

      ud    = u(ij)
      vd    = v(ij)

      duv(1) = ud-uc
      duv(2) = vd-vc

      graduv = dot_product(duv,duv)

      dduv(1) = -two * duy(ic) * dypn
      dduv(2) = -two * dvy(ic) * dypn

      gradCf   = dot_product(duv,dduv)

      uvCt    = one - graduv/(gradCf+small)

    ENDIF

    ufcd = u(ijn)*fyn+u(ij)*fyp
    vfcd = v(ijn)*fyn+v(ij)*fyp

    IF( uvct <= zero .OR. uvct >= one )THEN        ! 0.1 <= Beta_m <= 0.5 upwind
       uf = uc
       vf = vc
    ELSEIF( betam <= uvct .AND. uvct < one  )THEN  ! Central
       uf = ufcd
       vf = vfcd
    ELSEIF(  zero < uvct .AND. uvct < betam )THEN  ! Gamma
       fg = uvct / betam
       uf = fg*ufcd + (one-fg)*uc
       vf = fg*vfcd + (one-fg)*vc
    ENDIF

!    IF( uvct <= zero .OR. uvct >= one )THEN        ! 0.1 <= Beta_m <= 0.5 upwind
!      uf = uc
!      vf = vc
!    ELSEIF(  uvct <= half )THEN  ! Second order Upwind
!      uf = uc + duy(ic)*(y(j) - yc(jj))
!      vf = vc + dvy(ic)*(y(j) - yc(jj))
!    ELSE  ! Central
!      uf = ufcd
!      vf = vfcd
!    ENDIF

    fucds=f2(ij)*uf
    fvcds=f2(ij)*vf

    !--COEFFICIENTS AN(P) AND AS(N) DUE TO UDS
    an(ij) = cn-d
    as(ijn)=-cp-d

    !--SOURCE TERM CONTRIBUTIONS AT P AND N DUE TO DEFERRED CORRECTION

    su(ij) =su(ij) +gds(iu)*(fuuds-fucds)
    su(ijn)=su(ijn)-gds(iu)*(fuuds-fucds)
    sv(ij) =sv(ij) +gds(iu)*(fvuds-fvcds)
    sv(ijn)=sv(ijn)-gds(iu)*(fvuds-fvcds)

    IF(temp_visc)THEN
      dun = visi*s*(dvx(ijn)*fyn+dvx(ij)*fyp)
      su(ij)  = su(ij) + dun
      su(ijn) = su(ijn)- dun
      dvn = visi*s*(dvy(ijn)*fyn+dvy(ij)*fyp)
      sv(ij)  = sv(ij) + dvn
      sv(ijn) = sv(ijn)- dvn
    ENDIF

  END DO
END DO

!--VOLUME INTEGRALS (SOURCE TERMS)
!--Cell-face pressure calculated using linear interpolation;
!--cell volume is VOL, RP is the radius at node P; DX and DY
!--are the width and height of the cell. Contribution to AP
!--coefficient from volume integrals is stored temporarily
!--in arrays APU and APV for U and V, respectively; these
!--arrays are later used to store 1/AP, which is needed in
!--the pressure-correction equation.
DO i=2,nim

  dx=x(i)-x(i-1)

  DO j=2,njm
    dy=y(j)-y(j-1)
    rp=half*(r(j)+r(j-1))
    vol=dx*dy*rp
    ij=li(i)+j
    ije=ij+nj-xPeriodic*i/nim*((i-1)*nj)
    ijw=ij-nj+xPeriodic*Max(3-i,0)*(nim-1)*nj
    ijn=ij+1-yPeriodic*j/njm*(j-1)
    ijs=ij-1+yPeriodic*Max(3-j,0)*(njm-1)

    !--CELL-FACE PRESSURE, CELL-CENTER GRADIENT, SOURCE

    pe=p(ije)*fx(i)+p(ij)*(one-fx(i))
    pw=p(ij)*fx(i-1)+p(ijw)*(one-fx(i-1))
    pn=p(ijn)*fy(j)+p(ij)*(one-fy(j))
    ps=p(ij)*fy(j-1)+p(ijs)*(one-fy(j-1))
    dpx(ij)=(pe-pw)/dx
    dpy(ij)=(pn-ps)/dy
    su(ij)=su(ij)-dpx(ij)*vol
    sv(ij)=sv(ij)-dpy(ij)*vol

    !--BUOYANCY SOURCE CONTRIBUTION

    IF(lcal(ien)) THEN
      sb=-celbeta(ij)*vol*(t(ij)-tref)
      su(ij)=su(ij)+gravx*sb
      sv(ij)=sv(ij)+gravy*sb
    ENDIF

    !--AXISYMMETRIC CONTRIBUTION

    IF(laxis) THEN
      apv(ij)=apv(ij)+lamvisc(ij)*vol/rp**2
    ENDIF

    !--UNSTEADY TERM CONTRIBUTION TO AP AND SU

    IF(LTIME) THEN
      aptoo=denoo(ij)*vol*ct3
      apto =deno(ij)*vol*ct2
      apt  =den(ij)*vol*ct1
      su(ij)=su(ij)+apto*uo(ij)+aptoo*uoo(ij)
      sv(ij)=sv(ij)+apto*vo(ij)+aptoo*voo(ij)
      apv(ij)=apv(ij)+apt
      apu(ij)=apu(ij)+apt
    ENDIF

    IF(putobj)THEN
      su(ij) = su(ij) + fdsu(ij) + (den(ij)-densit) * gravx * vol !+ fdfcu(ij)
      sv(ij) = sv(ij) + fdsv(ij) + (den(ij)-densit) * gravy * vol !+ fdfcv(ij)
    ENDIF
  END DO
END DO

!--PROBLEM MODIFICATIONS - BOUNDARY CONDITIONS

CALL fd_bcuv

!--UNDER-RELAXATION, SOLVING EQUATION SYSTEM FOR U-VELOCITY
DO i=2,nim
  DO ij=li(i)+2,li(i)+njm
    ap(ij)=(-ae(ij)-aw(ij)-an(ij)-as(ij)+apu(ij))*urfu
    su(ij)=su(ij)+(one-urf(iu))*ap(ij)*u(ij)
    apu(ij)=one/ap(ij)
  END DO
END DO

IF(solver_type == solver_sparsekit)THEN

  IF(use_GPU == use_GPU_yes)THEN

    CALL fd_cooVal_create(ap,as,an,aw,ae,su,u)

!    config%relativeTolerance = sor(iu)
!    config%maxIterations = nsw(iu)
!    culaStat = culaSparseCudaDcooBicgstabJacobi(handle, config, platformOpts, formatOpts, &
!                                solverOpts, precondOpts, NCel, NNZ, Acoo, Arow, Acol, sol, rhs, res)

!    resor(iu) = res%residual%relative
     CALL cusp_biCGSTAB_copyH2D_system(Acoo, SOL, RHS, error)
     IF(error /= OPSUCCESS)GOTO 100

     CALL cusp_BiCGSTAB_solveDev_system(sor(iu),absTol,nsw(iu),error)
     IF(error /= OPSUCCESS)GOTO 100

     CALL cusp_BiCGSTAB_getMonitor(resor(iu),ipar(1),error)
     IF(error /= OPSUCCESS)GOTO 100

     CALL cusp_biCGSTAB_copyD2H_x(sol,error)
     IF(error /= OPSUCCESS)GOTO 100

!    IF(culaStat /= culaSparseNoError)THEN
!      WRITE(*,*)'CULA encoutered an error.'
!      culaStat=culaSparseDestroy(handle)
!      STOP
!    ENDIF

    CALL copy_solution(u,sol)
    
    100 CONTINUE
    IF(error /= OPSUCCESS)THEN
      WRITE(*,*)'GPU-- x-Momentum solver problem in CUDA operations.'
      STOP
      STOP
    ENDIF
 ELSE

    CALL fd_cooVal_create(ap,as,an,aw,ae,su,u)
    CALL coocsr(Ncel,NNZ,Acoo,Arow,Acol,Acsr,Aclc,Arwc)

    ipar  =  0
    fpar  = 0.0

    !--preconditioner
    !ipar(2) = 1     ! 1=left, 2=right
    ipar(2) = 0           ! 0 == no preconditioning
    ipar(3) = 0
    ipar(4) = NNZ*8       ! workspace for BGStab
    ipar(6) = nsw(iu)         ! maximum number of matrix-vector multiplies

    fpar(1) = sor(iu)  ! relative tolerance, must be between (0, 1)
    fpar(2) = 0.0  ! absolute tolerance, must be positive

        !call solve_spars(debug,IOdbg,Ncel,RHS,SOL,ipar,fpar,&
        !                       WORK,Acsr,Aclc,Arwc)
    CALL solve_spars(debug,out_unit,Ncel,RHS,SOL,ipar,fpar,&
                          WORK, Acsr,Aclc,Arwc,NNZ,Alu,Jlu,Ju,Jw)


    CALL copy_solution(u,sol)

    CALL calc_residual(ap,as,an,aw,ae,su,u,resor(iu))

  ENDIF

ELSEIF(solver_type == solver_sip)THEN

  CALL fd_solve_sip2d(ae,an,aw,as,ap,su,u,li,ni,nj,nsw(iu),resor(iu),sor(iu))

ELSEIF(solver_type == solver_cg)THEN !--cann't be use

  CALL fd_solve_sip2d(ae,an,aw,as,ap,su,u,li,ni,nj,nsw(iu),resor(iu),sor(iu))

ELSEIF(solver_type == solver_hypre)THEN

!  CALL solve_hypre(mpi_comm,li,lli,ncel,ni,nj,ap,as,an,aw,ae,su,u,sor(iu),resor(iu),&
!                    Hypre_A,Hypre_b,Hypre_x,2,nsw(iu))
ENDIF
!--UNDER-RELAXATION, SOLVING EQUATION SYSTEM FOR V-VELOCITY

DO i=2,nim
  DO ij=li(i)+2,li(i)+njm
    ap(ij)=(-ae(ij)-aw(ij)-an(ij)-as(ij)+apv(ij))*urfv
    sv(ij)=sv(ij)+(one-urf(iv))*ap(ij)*v(ij)
    apv(ij)=one/ap(ij)
  END DO
END DO

IF(solver_type == solver_sparsekit)THEN

  IF(use_GPU == use_GPU_yes)THEN

    CALL fd_cooVal_create(ap,as,an,aw,ae,sv,v)

!    config%relativeTolerance = sor(iv)
!    config%maxIterations = nsw(iv)
!    culaStat = culaSparseCudaDcooBicgstabJacobi(handle, config, platformOpts, formatOpts, &
!                                solverOpts, precondOpts, NCel, NNZ, Acoo, Arow, Acol, sol, rhs, res)

!    resor(iv) = res%residual%relative
     CALL cusp_biCGSTAB_copyH2D_system(Acoo, SOL, RHS, error)
     IF(error /= OPSUCCESS)GOTO 200

     CALL cusp_BiCGSTAB_solveDev_system(sor(iv),absTol,nsw(iv),error)
     IF(error /= OPSUCCESS)GOTO 200

     CALL cusp_BiCGSTAB_getMonitor(resor(iv),ipar(1),error)
     IF(error /= OPSUCCESS)GOTO 200

     CALL cusp_biCGSTAB_copyD2H_x(sol,error)
     IF(error /= OPSUCCESS)GOTO 200


!    IF(culaStat /= culaSparseNoError)THEN
!    WRITE(*,*)'CULA encoutered an error.'
!     culaStat=culaSparseDestroy(handle)
!      STOP
!    ENDIF

    CALL copy_solution(v,sol)

    200 CONTINUE
    IF(error /= OPSUCCESS)THEN
      WRITE(*,*)'GPU-- y-Momentum solver problem in CUDA operations.'
      STOP
    ENDIF
  ELSE

    CALL fd_cooVal_create(ap,as,an,aw,ae,sv,v)
    CALL coocsr(Ncel,NNZ,Acoo,Arow,Acol,Acsr,Aclc,Arwc)

    ipar  =  0
    fpar  = 0.0

    !--preconditioner
    !ipar(2) = 1     ! 1=left, 2=right
    ipar(2) = 0           ! 0 == no preconditioning
    ipar(3) = 0
    ipar(4) = NNZ*8       ! workspace for BGStab
    ipar(6) = nsw(iv)         ! maximum number of matrix-vector multiplies

    fpar(1) = sor(iv)  ! relative tolerance, must be between (0, 1)
    fpar(2) = 0.0  ! absolute tolerance, must be positive

        !call solve_spars(debug,IOdbg,Ncel,RHS,SOL,ipar,fpar,&
        !                       WORK,Acsr,Aclc,Arwc)
    CALL solve_spars(debug,out_unit,Ncel,RHS,SOL,ipar,fpar,&
                          WORK, Acsr,Aclc,Arwc,NNZ,Alu,Jlu,Ju,Jw)

    CALL copy_solution(v,sol)

    CALL calc_residual(ap,as,an,aw,ae,sv,v,resor(iv))

  ENDIF

ELSEIF(solver_type == solver_sip)THEN

  CALL fd_solve_sip2d(ae,an,aw,as,ap,sv,v,li,ni,nj,nsw(iv),resor(iv),sor(iv))

ELSEIF(solver_type == solver_cg)THEN !--cann't be used

  CALL fd_solve_sip2d(ae,an,aw,as,ap,sv,v,li,ni,nj,nsw(iv),resor(iv),sor(iv))

ELSEIF(solver_type == solver_hypre)THEN
  !CALL solve_hypre(mpi_comm,li,lli,ncel,ni,nj,ap,as,an,aw,ae,sv,v,sor(iv),resor(iv),&
  !                  Hypre_A,Hypre_b,Hypre_x,2,nsw(iv))
ENDIF

!--Calculate velocity gradients !--This is fine as long as no symmetry BC present
!--otherwise boundary velocities should be set
DO i=2,nim
  dx=x(i)-x(i-1)
  DO j=2,njm
    dy=y(j)-y(j-1)
    ij=li(i)+j
    ije=ij+nj-xPeriodic*i/nim*((i-1)*nj)
    ijw=ij-nj+xPeriodic*Max(3-i,0)*(nim-1)*nj
    ijn=ij+1-yPeriodic*j/njm*(j-1)
    ijs=ij-1+yPeriodic*Max(3-j,0)*(njm-1)
    !--CELL-CENTER GRADIENT

    vele=u(ije)*fx(i)+u(ij)*(one-fx(i))
    velw=u(ij)*fx(i-1)+u(ijw)*(one-fx(i-1))
    veln=u(ijn)*fy(j)+u(ij)*(one-fy(j))
    vels=u(ij)*fy(j-1)+u(ijs)*(one-fy(j-1))
    dux(ij)=(vele-velw)/dx
    duy(ij)=(veln-vels)/dy

    vele=v(ije)*fx(i)+v(ij)*(one-fx(i))
    velw=v(ij)*fx(i-1)+v(ijw)*(one-fx(i-1))
    veln=v(ijn)*fy(j)+v(ij)*(one-fy(j))
    vels=v(ij)*fy(j-1)+v(ijs)*(one-fy(j-1))
    dvx(ij)=(vele-velw)/dx
    dvy(ij)=(veln-vels)/dy

  ENDDO
ENDDO

!IF(temp_visc)CALL fd_calc_bcuv_grad

END SUBROUTINE fd_calc_mom

END MODULE modfd_calc_mom
