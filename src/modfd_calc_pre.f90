MODULE modfd_calc_pre

PRIVATE
PUBLIC :: fd_calc_pre
CONTAINS
!============================================
SUBROUTINE fd_calc_pre

USE parameters,         ONLY : out_unit,solver_sparsekit,solver_sip,solver_cg,solver_hypre
USE real_parameters,    ONLY : one,zero,half,two,four
USE precision,          ONLY : r_single
USE modfd_set_bc,       ONLY : fd_bcpressure
USE shared_data,        ONLY : urf,ip,p,su,apu,apv,nim,njm,&
                               xc,ni,nj,li,fx,y,r,visc,su,&
                               u,v,ae,aw,an,as,fy,yc,x,lcal,ien,&
                               den,deno,laxis,p,dtr,&
                               gamt,sor,resor,nsw,f1,f2,dpx,dpy,&
                               ltime,ap,ipr,jpr,ltest,pp,fdsu,fdsv,&
                               dux,duy,dvx,dvy,fdsuc,fdsvc,fdsub,fdsvb,putobj,&
                               Ncel,NNZ,Acoo,Arow,Acol,Acsr,Aclc,Arwc,&
                               solver_type,rhs,sol,work,alu,jlu,ju,jw,dt,&
                               Hypre_A,Hypre_b,Hypre_x,mpi_comm,lli
USE  modfd_solve_linearsys,  ONLY : fd_solve_sip2d,fd_spkit_interface,copy_solution,calc_residual,&
                                    fd_solve_cgs2d
IMPLICIT NONE
REAL(KIND = r_single) :: fxe,fxp,dxpe,s,d,fyn,fyp,dypn,&
                         dx,dy,rp,ppe,ppw,ppn,pps,dpxel,&
                         uel,apue,ue,vole,voln,dpynl,vnl,apvn,&
                         vn,sum,ppo,dpxe,dpyn,pcor !,&
                         !duxel,duxe,dvynl,dvyn,fdsuce,fdsvcn,dfyn,dfxe,c,K

INTEGER               :: ij,i,j,ije,ijn,ijpref,ipar(16),debug
REAL                  :: fpar(16)

debug = 0
!--EAST CV FACES (S - AREA, VOLE - VOLUME BETWEEN P AND E)
DO i=2,nim-1
  
  dxpe=xc(i+1)-xc(i)
  fxe=fx(i)
  fxp=one-fxe

  DO j=2,njm
    ij=li(i)+j
    ije=ij+nj
    !ijw=ij-nj
    s=(y(j)-y(j-1))*(r(j)+r(j-1))*half
    vole=dxpe*s
    d=(den(ije)*fxe+den(ij)*fxp)*s
     
    !--INTERPOLATED CELL FACE QUANTITIES (PRESSURE GRAD., U AND 1/AP)
    !--Note: pressure gradient is interpolated midway between P and E,
    !--since the gradient calculated at cell face is second order
    !--accurate at that location; the velocity is interpolated linearly,
    !--to achieve second order accuracy at cell face center.
 
    dpxel=half*(dpx(ije)+dpx(ij))
    uel=u(ije)*fxe+u(ij)*fxp
    apue=apu(ije)*fxe+apu(ij)*fxp
    !duxel=half*(dux(ije)+dux(ij))
    !--CELL FACE GRADIENT, VELOCITY AND MASS FLUX
    dpxe=(p(ije)-p(ij))/dxpe
    !duxe=(u(ije)-u(ij))/dxpe
    !dfxe=(fdsu(ije)-fdsu(ij))/dxpe
    pcor=apue*vole*(dpxe-dpxel)
!    IF(uel >= pcor)THEN
!      c = one
!    ELSE
!      c = zero
!    ENDIF
!    IF(.NOT. (ABS(u(ije)) + two*ABS(u(ij)) + ABS(u(ijw))) > zero)THEN
!      K = zero
!    ELSE
!      K = ABS(dux(ij))/(ABS(u(ije)) + two*ABS(u(ij)) + ABS(u(ijw)))
!    ENDIF
    !Rahman correction : + dxpe/four*(apu(ij) - apu(ije))*dfxe
    !Gu correction: - K*dxpe**2*(duxe - duxel)
    !test correction: + apue*vole*((fdsu(ij)+fdsu(ije))/two - (fdsub(ij)+fdsub(ije))/two)
    ue=uel - pcor 
    
    f1(ij) = d*ue
!    IF(putobj)THEN
!      fdsuce = fdsuc(ije)*fxe+fdsuc(ij)*fxp
!      fdsub(ij) = fdsuce*d*apue
!    ENDIF
    !--COEFFICIENTS OF P' EQUATION, AE(P) AND AW(E)
    
    ae(ij)=-d*apue*s
    aw(ije)=ae(ij)

  END DO
END DO

!--NORTH CV FACES (S - AREA, VOLN - VOLUME BETWEEN P AND N)
DO j=2,njm-1
  dypn=yc(j+1)-yc(j)
  fyn=fy(j)
  fyp=one-fyn

  DO i=2,nim
    ij=li(i)+j
    ijn=ij+1
    !ijs=ij-1

    s=(x(i)-x(i-1))*r(j)
    voln=s*dypn
    d=(den(ijn)*fyn+den(ij)*fyp)*s
    
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
!    IF(vnl >= pcor)THEN
!      c = one
!    ELSE
!      c = zero
!    ENDIF
!    IF(.NOT. (ABS(v(ijn)) + two*ABS(v(ij)) + ABS(v(ijs))) > zero)THEN
!      K = zero
!    ELSE
!      K = ABS(dvy(ij))/(ABS(v(ijn)) + two*ABS(v(ij)) + ABS(v(ijs)))
!    ENDIF
    !Rahman correction: +dypn/four*(apv(ij) - apv(ijn))*dfyn 
    !Gu correction: - K*dxpe**2*(dvyn - dvynl)
    !test correction: + apue*vole*((fdsu(ij)+fdsu(ije))/two - (fdsub(ij)+fdsub(ije))/two) 
    vn=vnl - pcor 
    
    f2(ij)=d*vn
!    IF(putobj)THEN
!      fdsvcn = fdsvc(ijn)*fyn+fdsvc(ij)*fyp    
!      fdsvb(ij) = fdsvcn*d*apvn
!    ENDIF
    !--COEFFICIENTS OF P' EQUATION, AN(P) AND AS(N)

    an(ij)=-d*apvn*s
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
    su(ij)=f1(ij-nj)-f1(ij)+f2(ij-1)-f2(ij) 
    IF(putobj)THEN
      su(ij) = su(ij) +  (deno(ij) - den(ij)) *  ((x(i)-x(i-1))*(y(j)-y(j-1))*half*(r(j)+r(j-1))) /dt
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
  
  CALL fd_spkit_interface(ap,as,an,aw,ae,su,pp)
  CALL coocsr(Ncel,NNZ,Acoo,Arow,Acol,Acsr,Aclc,Arwc)
  
  SOL = 0.0
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
                        WORK, Acsr,Aclc,Arwc,        &
	                         NNZ,  Alu,Jlu,Ju,Jw          )
  
  CALL copy_solution(pp,sol)

  !CALL calc_residual(ap,as,an,aw,ae,su,pp,resor(ip))

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
DO i=2,nim-1
  DO ij=li(i)+2,li(i)+njm
    f1(ij)=f1(ij)+ae(ij)*(pp(ij+nj)-pp(ij))
  END DO
END DO

!--CORRECT NORTH MASS FLUXES 

DO i=2,nim
  DO ij=li(i)+2,li(i)+njm-1
    f2(ij)=f2(ij)+an(ij)*(pp(ij+1)-pp(ij))
  END DO
END DO

sum=zero
DO i=2,nim
  DO j=2,njm
    ij=li(i)+j
    su(ij)=f1(ij-nj)-f1(ij)+f2(ij-1)-f2(ij) 
    IF(putobj)THEN
      su(ij) = su(ij) +  (deno(ij) - den(ij)) *  ((x(i)-x(i-1))*(y(j)-y(j-1))*half*(r(j)+r(j-1))) /dt
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

    ppe=pp(ij+nj)*fx(i)+pp(ij)*(one-fx(i))
    ppw=pp(ij)*fx(i-1)+pp(ij-nj)*(one-fx(i-1))
    ppn=pp(ij+1)*fy(j)+pp(ij)*(one-fy(j))
    pps=pp(ij)*fy(j-1)+pp(ij-1)*(one-fy(j-1))

    u(ij)=u(ij)-(ppe-ppw)*dy*rp*apu(ij)
    v(ij)=v(ij)-(ppn-pps)*dx*rp*apv(ij)
    p(ij)=p(ij)+urf(ip)*(pp(ij)-ppo)

  END DO
END DO

END SUBROUTINE fd_calc_pre

END MODULE modfd_calc_pre
