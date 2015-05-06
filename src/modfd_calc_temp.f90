MODULE modfd_calc_temp

PRIVATE
PUBLIC :: fd_calc_temp
CONTAINS
!=====================================
SUBROUTINE fd_calc_temp
!--Best to solve enthalpy, this is not very stable, but for now ....
USE parameters,         ONLY : out_unit,solver_sparsekit,solver_sip,solver_cg,solver_hypre,use_GPU_yes,OPSUCCESS,SOLVER_DONE,&
                               SOLVER_FAILED
USE real_parameters,    ONLY : one,zero,half,two,small,betam
USE precision,          ONLY : r_single
USE modfd_set_bc,       ONLY : fd_bctemp
USE shared_data,        ONLY : urf,ien,t,to,too,su,nim,njm,&
                               xc,ni,nj,li,fx,y,r,visc,su,&
                               u,v,ae,aw,an,as,fy,yc,x,lcal,ien,&
                               denoo,deno,den,laxis,sor,resor,ct1,ct2,ct3,&
                               nsw,ft1,ft2,dpx,dpy,ltime,ap,ltest,celkappa,celcp,celcpo,celcpoo,gds,&
                               putobj,fdst,dtx,dty,apt,xperiodic,yperiodic,&
                               Ncel,NNZ,Acoo,Arow,Acol,Acsr,Aclc,Arwc,&
                               solver_type,rhs,sol,work,alu,jlu,ju,jw,&
                               Hypre_A,Hypre_b,Hypre_x,mpi_comm,lli,dux,dvy,p,use_GPU,&
                               handle, config, platformOpts, formatOpts, solverOpts, precondOpts, res, culaStat

use cula_sparse_type
use cula_sparse

USE  modfd_solve_linearsys,  ONLY : fd_solve_sip2d,fd_spkit_interface,copy_solution,calc_residual,fd_cooVal_create
!USE modcu_BiCGSTAB,          ONLY : cu_cpH2D_sysDP,cu_BiCGSTAB_setStop,cu_BiCGSTAB_itr,&
!                                  cu_cpD2H_solDP

IMPLICIT NONE
REAL(KIND = r_single) :: urfi,fxe,fxp,dxpe,s,d,fyn,fyp,dypn,cn,&
                         ce,cp,dx,dy,rp,fuds,fcds,vol,aptt,aptto,apttoo,te,tw,ts,tn,&
                         phic,phid,phict,phif,phifcd,fg,dphi,ddphi
INTEGER               :: ij,i,j,ije,ijn,ijs,ijw,ipar(16),debug,error,xend,yend,ic,ii,jj
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

    !--EXPLICIT CONVECTIVE FLUX FOR UDS AND CDS
    ce=MIN(ft1(ij),zero)
    cp=MAX(ft1(ij),zero)

    fuds=cp*t(ij)+ce*t(ije)

    IF( ft1(ij) >= zero )THEN
      
      ic    = ij 
      ii    = i
      phic  = t(ij)
      phid  = t(ije)
            
      dphi  = phid - phic
      
      ddphi = two * dtx(ic) * dxpe 
      phiCt = one - dphi/(ddphi+small)

    ELSE
      
      ic    = ije
      ii    = i+1
      phic  = t(ije)
      phid  = t(ij)
      
      dPhi  = phid - phic
     
      ddPhi = - two * dtx(ic) * dxpe
      phict = one - dphi/(ddphi+Small)

    ENDIF

    phifcd = t(ije)*fxe+t(ij)*fxp

    IF( phict <= zero .OR. phict >= one )THEN       ! Upwind
                                                      
      phif = phic                              
                                                       
    ELSEIF( betam <= phict .AND. phict < one )THEN  ! select CDS
                                                       
      phif = phifcd                              
                                                       
    ELSEIF( zero < phict .AND. phict < betam )then  ! select gamma
                                                       
      fg = phict / betam                           
      phif = fg*phifcd + (one-fg)*phic           
   
    ENDIF
!    IF( phict <= zero .OR. phict >= one )THEN       ! Upwind
!      phif = phic
!    ELSEIF(  phict <= half )THEN  ! Second order Upwind 
!      phif = phic + dtx(ic)*(x(i) - xc(ii)) 
!    ELSE  ! Central  
!      phif = phifcd                                          
!    ENDIF

    fcds=ft1(ij)*phif

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

    !--EXPLICIT CONVECTIVE FLUXES FOR UDS AND CDS
    cn=MIN(ft2(ij),zero)
    cp=MAX(ft2(ij),zero)

    fuds=cp*t(ij)+cn*t(ijn)

    IF( ft2(ij) >= zero )THEN
      
      ic    = ij
      jj    = j
      phic  = t(ij)
      phid  = t(ijn)

      dphi  = phid - phic
      
      ddphi = two * dty(ij) * dypn 
      phiCt = one - dphi/(ddphi+small)

    ELSE
      
      ic    = ijn
      jj    = j+1
      phic  = t(ijn)
      phid  = t(ij)

      dPhi  = phid - phic
     
      ddPhi = - two * dty(ijn) * dypn
      phict = one - dphi/(ddphi+small)

    ENDIF

    phifcd = t(ijn)*fyn+t(ij)*fyp

    IF( phict <= zero .or. phict >= one )THEN       ! Upwind
                                                      
      phif = phic                              
                                                       
    ELSEIF( betam <= phict .AND. phict < one )THEN  ! select CDS
                                                       
      phif = phifcd                              
                                                       
    ELSEIF( zero < phict .AND. phict < betam )then  ! select gamma
                                                       
      fg = phict / betam                           
      phif = fg*phifcd + (one-fg)*phic           
   
    ENDIF
!    IF( phict <= zero .OR. phict >= one )THEN       ! Upwind
!      phif = phic
!    ELSEIF(  phict <= half )THEN  ! Second order Upwind 
!      phif = phic + dty(ic)*(y(j) - yc(jj)) 
!    ELSE  ! Central  
!      phif = phifcd                                          
!    ENDIF

    fcds=ft2(ij)*phif

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
      apttoo=celcp(ij)*vol*ct3
      aptto =celcp(ij)*vol*ct2
      aptt  =celcp(ij)*vol*ct1
      
      su(ij)=su(ij)+aptto*to(ij)+apttoo*too(ij)
      ap(ij)=ap(ij)+aptt
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
    
    config%relativeTolerance = sor(ien)
    config%maxIterations = nsw(ien)
    culaStat = culaSparseCudaDcooBicgstabJacobi(handle, config, platformOpts, formatOpts, &
                                solverOpts, precondOpts, NCel, NNZ, Acoo, Arow, Acol, sol, rhs, res)
  
    resor(ien) = res%residual%relative

!    CALL cu_cpH2D_sysDP(Acoo, RHS, SOL,error)
!    IF(error /= OPSUCCESS)GOTO 100
    
!    CALL cu_BiCGSTAB_setStop(nsw(ien),sor(ien),error)
!    IF(error /= OPSUCCESS)GOTO 100
    
!    CALL  cu_BiCGSTAB_itr(resor(ien),ipar(1),error)
!    IF(error /= SOLVER_DONE)GOTO 100
    
!    CALL  cu_cpD2H_solDP(sol,error)
!    IF(error /= OPSUCCESS)GOTO 100
    
    IF(culaStat /= culaSparseNoError)THEN
      WRITE(*,*)'CULA encoutered an error.'
      culaStat=culaSparseDestroy(handle)
      STOP
    ENDIF

    CALL copy_solution(t,sol)

!    100 CONTINUE
!    IF(error == SOLVER_FAILED)THEN
!     WRITE(*,*)'GPU--Temperature solver failed to find solution.'
!      STOP
!      STOP
!    ELSEIF(error /= OPSUCCESS)THEN
!      WRITE(*,*)'GPU--Temperature solver problem in CUDA operations.'
!      STOP
!    ENDIF
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
                          WORK, Acsr,Aclc,Arwc,NNZ,Alu,Jlu,Ju,Jw)
  
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
