MODULE modfd_calc_temp

PRIVATE
PUBLIC :: fd_calc_temp
CONTAINS
!=====================================
SUBROUTINE fd_calc_temp
!--Best to solve enthalpy, this is not very stable, but for now ....
USE parameters,         ONLY : out_unit,solver_sparsekit,solver_sip,solver_cg,solver_hypre,use_GPU_yes,OPSUCCESS,SOLVER_DONE,&
                               SOLVER_FAILED
USE real_parameters,    ONLY : one,zero,half
USE precision,          ONLY : r_single
USE modfd_set_bc,       ONLY : fd_bctemp
USE shared_data,        ONLY : urf,ien,t,to,too,su,nim,njm,&
                               xc,ni,nj,li,fx,y,r,visc,su,&
                               u,v,ae,aw,an,as,fy,yc,x,lcal,ien,&
                               deno,den,laxis,dtr,gamt,sor,resor,&
                               nsw,f1,f2,dpx,dpy,ltime,ap,ltest,celkappa,celcp,celcpo,gds,&
                               putobj,fdst,dtx,dty,apt,xperiodic,yperiodic,&
                               Ncel,NNZ,Acoo,Arow,Acol,Acsr,Aclc,Arwc,&
                               solver_type,rhs,sol,work,alu,jlu,ju,jw,&
                               Hypre_A,Hypre_b,Hypre_x,mpi_comm,lli,dux,dvy,p,use_GPU
USE  modfd_solve_linearsys,  ONLY : fd_solve_sip2d,fd_spkit_interface,copy_solution,calc_residual,fd_cooVal_create
USE modcu_BiCGSTAB,          ONLY : cu_cpH2D_sysDP,cu_BiCGSTAB_setStop,cu_BiCGSTAB_itr,&
                                  cu_cpD2H_solDP

IMPLICIT NONE
REAL(KIND = r_single) :: urfi,fxe,fxp,dxpe,s,d,fyn,fyp,dypn,cn,&
                         ce,cp,dx,dy,rp,fuds,fcds,vol,aptt,te,tw,ts,tn,cpe,cpn
INTEGER               :: ij,i,j,ije,ijn,ijs,ijw,ipar(16),debug,error,xend,yend
REAL                  :: fpar(16)

debug = 0

su = zero
ap = zero
    
urfi=one/urf(ien)

!--FLUXES THROUGH INTERNAL EAST CV-FACES
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

    !--COEFFICIENT RESULTING FROM DIFFUSIVE FLUX
    d=(celkappa(ije)*fxe+celkappa(ij)*fxp)*s/dxpe
    cpe = (celcp(ije)*fxe+celcp(ij)*fxp)

    !--EXPLICIT CONVECTIVE FLUX FOR UDS AND CDS
    ce=MIN(f1(ij),zero)
    cp=MAX(f1(ij),zero)

    fuds=cp*t(ij)+ce*t(ije)
    fcds=cpe*f1(ij)*(t(ije)*fxe+t(ij)*fxp)

    !--COEFFICIENTS AE(P) AND AW(E) DUE TO UDS

    ae(ij) = ce-d
    aw(ije)=-cp-d

    !--SOURCE TERM CONTRIBUTIONS AT P AND E DUE TO DEFERRED CORRECTION

    su(ij) =su(ij) +gds(ien)*(fuds-fcds)
    su(ije)=su(ije)-gds(ien)*(fuds-fcds)
  END DO
END DO

!--FLUXES THROUGH INTERNAL NORTH CV FACES 
yend = njm+yPeriodic-1 !--Apply periodic conditions
DO j=2,yend

  !--INTERPOLATION FACTORS, DISTANCE FROM P TO N (SAME FOR ALL J)

  fyn =fy(j)
  fyp =one-fyn
  dypn=yc(j+1)-yc(j)

  DO i=2,nim
    ij =li(i)+j
    ijn=ij+1
    ijn=ij+1-j/njm*(j-1)
    !--CELL FACE AREA S = DX*RN*1

    s=(x(i)-x(i-1))*r(j)
    
    !--COEFFICIENT RESULTING FROM DIFFUSIVE FLUX (SAME FOR U AND V)
    d=(celkappa(ijn)*fyn+celkappa(ij)*fyp)*s/dypn
    cpn=(celcp(ijn)*fyn+celcp(ij)*fyp)

    !--EXPLICIT CONVECTIVE FLUXES FOR UDS AND CDS
    cn=MIN(f2(ij),zero)
    cp=MAX(f2(ij),zero)

    fuds=cp*t(ij)+cn*t(ijn)
    fcds=cpn*f2(ij)*(t(ijn)*fyn+t(ij)*fyp)

    !--COEFFICIENTS AE(P) AND AW(E) DUE TO UDS

    an(ij) = cn-d
    as(ijn)=-cp-d
    
    !--SOURCE TERM CONTRIBUTIONS AT P AND E DUE TO DEFERRED CORRECTION

    su(ij) =su(ij) +gds(ien)*(fuds-fcds)
    su(ijn)=su(ijn)-gds(ien)*(fuds-fcds)

  END DO
END DO

!--VOLUME INTEGRALS (SOURCE TERMS)
DO i=2,nim
  dx=x(i)-x(i-1)

  DO j=2,njm
    ij=li(i)+j
    dy=y(j)-y(j-1)
    rp=half*(r(j)+r(j-1))
    vol=dx*dy*rp

    !--UNSTEADY TERM CONTRIBUTION TO AP AND SU
    IF(ltime) THEN
      aptt=celcpo(ij)*deno(ij)*vol*dtr
      su(ij)=su(ij)+(one+gamt)*aptt*to(ij)-half*gamt*aptt*too(ij)
      ap(ij)=ap(ij)+(one+half*gamt)*aptt
    ENDIF

    IF(putobj)THEN
      su(ij) = su(ij) + fdst(ij)  !- p(ij)*(dux(ij) + dvy(ij))*vol
    ENDIF
  END DO
END DO

!--PROBLEM MODIFICATIONS - BOUNDARY CONDITIONS
CALL fd_bctemp

!--UNDER-RELAXATION, SOLVING EQUATION SYSTEM FOR TEMPERATURE

DO i=2,nim
  DO ij=li(i)+2,li(i)+njm
    ap(ij)=(ap(ij)-aw(ij)-ae(ij)-an(ij)-as(ij))*urfi
    su(ij)=su(ij)+(one-urf(ien))*ap(ij)*t(ij)
  END DO
END DO

IF(solver_type == solver_sparsekit)THEN

  
  IF(use_GPU == use_GPU_yes)THEN
    
    CALL fd_cooVal_create(ap,as,an,aw,ae,su,t)

    CALL cu_cpH2D_sysDP(Acoo, RHS, SOL,error)
    IF(error /= OPSUCCESS)GOTO 100
    
    CALL cu_BiCGSTAB_setStop(nsw(ien),sor(ien),error)
    IF(error /= OPSUCCESS)GOTO 100
    
    CALL  cu_BiCGSTAB_itr(resor(ien),ipar(1),error)
    IF(error /= SOLVER_DONE)GOTO 100
    
    CALL  cu_cpD2H_solDP(sol,error)
    IF(error /= OPSUCCESS)GOTO 100
    
    CALL copy_solution(t,sol)

    100 CONTINUE
    IF(error /= OPSUCCESS .OR. error == SOLVER_FAILED)THEN
      WRITE(*,*)'GPU--solver problem.'
      PAUSE
    ENDIF
  ELSE
    
    CALL fd_cooVal_create(ap,as,an,aw,ae,su,t)
    CALL coocsr(Ncel,NNZ,Acoo,Arow,Acol,Acsr,Aclc,Arwc)
  
    ipar  =  0
    fpar  = 0.0
  
    !--preconditioner
    !ipar(2) = 1     ! 1=left, 2=right
    ipar(2) = 0           ! 0 == no preconditioning
    ipar(3) = 0
    ipar(4) = NNZ*8       ! workspace for BGStab
    ipar(6) = nsw(ien)         ! maximum number of matrix-vector multiplies

    fpar(1) = sor(ien)  ! relative tolerance, must be between (0, 1)
    fpar(2) = 0.0  ! absolute tolerance, must be positive

        !call solve_spars(debug,IOdbg,Ncel,RHS,SOL,ipar,fpar,&
        !                       WORK,Acsr,Aclc,Arwc)
    CALL solve_spars(debug,out_unit,Ncel,RHS,SOL,ipar,fpar,&
                          WORK, Acsr,Aclc,Arwc,        &
	                           NNZ,  Alu,Jlu,Ju,Jw          )
  
    CALL copy_solution(t,sol)
 
    CALL calc_residual(ap,as,an,aw,ae,su,t,resor(ien))
  
  ENDIF

ELSEIF(solver_type == solver_sip)THEN
  
  CALL fd_solve_sip2d(ae,an,aw,as,ap,su,t,li,ni,nj,nsw(ien),resor(ien),sor(ien))

ELSEIF(solver_type == solver_cg)THEN !--cann't be use

  CALL fd_solve_sip2d(ae,an,aw,as,ap,su,t,li,ni,nj,nsw(ien),resor(ien),sor(ien))

ELSEIF(solver_type == solver_hypre)THEN

  !CALL solve_hypre(mpi_comm,li,lli,ncel,ni,nj,ap,as,an,aw,ae,su,t,sor(ien),resor(ien),&
  !                  Hypre_A,Hypre_b,Hypre_x,2,nsw(ien))
ENDIF

!--Calculate temperature gradients
DO i=2,nim 
  dx=x(i)-x(i-1)
  DO j=2,njm
    dy=y(j)-y(j-1)
    ij=li(i)+j
    ije=ij+nj-i/nim*((i-1)*nj)
    ijw=ij-nj+xPeriodic*Max(3-i,0)*(nim-1)*nj
    ijn=ij+1-j/njm*(j-1)
    ijs=ij-1+yPeriodic*Max(3-j,0)*(njm-1)
    !--CELL-CENTER GRADIENT 
    
    te=t(ije)*fx(i)+t(ij)*(one-fx(i))
    tw=t(ij)*fx(i-1)+t(ijw)*(one-fx(i-1))
    tn=t(ijn)*fy(j)+t(ij)*(one-fy(j))
    ts=t(ij)*fy(j-1)+t(ijs)*(one-fy(j-1))
    dtx(ij)=(te-tw)/dx
    dty(ij)=(tn-ts)/dy

  ENDDO
ENDDO

END SUBROUTINE fd_calc_temp

END MODULE modfd_calc_temp