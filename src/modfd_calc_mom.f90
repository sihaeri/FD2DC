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
USE real_parameters,    ONLY : one,zero,half
USE shared_data,        ONLY : urf,iu,iv,p,su,sv,apu,apv,nim,njm,&
                               xc,ni,nj,li,fx,y,r,su,sv,gds,&
                               u,v,ae,aw,an,as,fy,yc,x,lcal,ien,&
                               gravx,gravy,beta,den,deno,laxis,p,dtr,&
                               gamt,vo,uo,voo,uoo,sor,resor,nsw,f1,f2,&
                               dpx,dpy,t,tref,ltime,ap,fdsu,fdsv,putobj,ibsu,ibsv,nfil,&
                               dvx,dvy,dux,duy,densit,temp_visc,lamvisc,&
                               Ncel,NNZ,Acoo,Arow,Acol,Acsr,Aclc,Arwc,&
                               solver_type,rhs,sol,work,alu,jlu,ju,jw,&
                               Hypre_A,Hypre_b,Hypre_x,mpi_comm,lli,celbeta,&
                               fdfcu,fdfcv,use_GPU
USE modfd_set_bc,       ONLY : fd_bcpressure,fd_bcuv,fd_calc_bcuv_grad
USE precision,          ONLY : r_single
USE  modfd_solve_linearsys,  ONLY : fd_solve_sip2d,fd_spkit_interface,copy_solution,calc_residual,fd_cooVal_create
USE modcu_BiCGSTAB,          ONLY : cu_cpH2D_sysDP,cu_BiCGSTAB_setStop,cu_BiCGSTAB_itr,&
                                  cu_cpD2H_solDP
IMPLICIT NONE

REAL(KIND = r_single) :: urfu,urfv,fxe,fxp,dxpe,s,d,cp,ce,&
                         fuuds,fvuds,fucds,fvcds,fyn,fyp,dypn,&
                         dx,dy,rp,vol,sb,pe,pw,pn,ps,apt,cn,&
                         vele,velw,veln,vels,due,dve,dun,dvn,visi
INTEGER               :: ij,i,j,ije,ijn,ipar(16),debug,error
REAL                  :: fpar(16)

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
DO i=2,nim-1

!--INTERPOLATION FACTORS, DISTANCE FROM P TO E (SAME FOR ALL J)
  fxe =fx(i)
  fxp =one-fxe
  dxpe=xc(i+1)-xc(i)

  DO j=2,njm
    ij=li(i)+j
    ije=ij+nj

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
    fucds=f1(ij)*(u(ije)*fxe+u(ij)*fxp)
    fvcds=f1(ij)*(v(ije)*fxe+v(ij)*fxp)

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

DO j=2,njm-1

!--INTERPOLATION FACTORS, DISTANCE FROM P TO N (SAME FOR ALL J)
  fyn =fy(j)
  fyp =one-fyn
  dypn=yc(j+1)-yc(j)

  DO i=2,nim
    ij =li(i)+j
    ijn=ij+1

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
    fucds=f2(ij)*(u(ijn)*fyn+u(ij)*fyp)
    fvcds=f2(ij)*(v(ijn)*fyn+v(ij)*fyp)

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
    !--CELL-FACE PRESSURE, CELL-CENTER GRADIENT, SOURCE 
    
    pe=p(ij+nj)*fx(i)+p(ij)*(one-fx(i))
    pw=p(ij)*fx(i-1)+p(ij-nj)*(one-fx(i-1))
    pn=p(ij+1)*fy(j)+p(ij)*(one-fy(j))
    ps=p(ij)*fy(j-1)+p(ij-1)*(one-fy(j-1))
    dpx(ij)=(pe-pw)/dx
    dpy(ij)=(pn-ps)/dy
    su(ij)=su(ij)-dpx(ij)*vol
    sv(ij)=sv(ij)-dpy(ij)*vol

    !--BUOYANCY SOURCE CONTRIBUTION

    IF(lcal(ien)) THEN
      sb=-celbeta(ij)*den(ij)*vol*(t(ij)-tref)
      su(ij)=su(ij)+gravx*sb
      sv(ij)=sv(ij)+gravy*sb
    ENDIF

    !--AXISYMMETRIC CONTRIBUTION

    IF(laxis) THEN
      apv(ij)=apv(ij)+lamvisc(ij)*vol/rp**2
    ENDIF

    !--UNSTEADY TERM CONTRIBUTION TO AP AND SU

    IF(LTIME) THEN
      apt=deno(ij)*vol*dtr
      su(ij)=su(ij)+(one+gamt)*apt*uo(ij)-half*gamt*apt*uoo(ij)
      sv(ij)=sv(ij)+(one+gamt)*apt*vo(ij)-half*gamt*apt*voo(ij)
      apv(ij)=apv(ij)+(one+half*gamt)*apt
      apu(ij)=apu(ij)+(one+half*gamt)*apt
    ENDIF

    IF(putobj)THEN
      IF(nfil > 0)THEN
        su(ij) = su(ij) + ibsu(ij) + fdsu(ij) + (den(ij)-densit) * gravx * vol !+ fdfcu(ij) 
        sv(ij) = sv(ij) + ibsv(ij) + fdsv(ij) + (den(ij)-densit) * gravy * vol !+ fdfcv(ij)
      ELSE
        su(ij) = su(ij) + fdsu(ij) + (den(ij)-densit) * gravx * vol !+ fdfcu(ij)
        sv(ij) = sv(ij) + fdsv(ij) + (den(ij)-densit) * gravy * vol !+ fdfcv(ij)
      ENDIF
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

    CALL cu_cpH2D_sysDP(Acoo, RHS, SOL,error)
    IF(error /= OPSUCCESS)GOTO 100
    
    CALL cu_BiCGSTAB_setStop(nsw(iu),sor(iu),error)
    IF(error /= OPSUCCESS)GOTO 100
    
    CALL  cu_BiCGSTAB_itr(resor(iu),ipar(1),error)
    IF(error /= SOLVER_DONE)GOTO 100
    
    CALL  cu_cpD2H_solDP(sol,error)
    IF(error /= OPSUCCESS)GOTO 100
    
    CALL copy_solution(u,sol)

    100 CONTINUE
    IF(error /= OPSUCCESS .OR. error == SOLVER_FAILED)THEN
      WRITE(*,*)'GPU--solver problem.'
      PAUSE
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
                          WORK, Acsr,Aclc,Arwc,        &
	                           NNZ,  Alu,Jlu,Ju,Jw          )
  

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

    CALL cu_cpH2D_sysDP(Acoo, RHS, SOL,error)
    IF(error /= OPSUCCESS)GOTO 200
    
    CALL cu_BiCGSTAB_setStop(nsw(iv),sor(iv),error)
    IF(error /= OPSUCCESS)GOTO 200
    
    CALL  cu_BiCGSTAB_itr(resor(iv),ipar(1),error)
    IF(error /= SOLVER_DONE)GOTO 200
    
    CALL  cu_cpD2H_solDP(sol,error)
    IF(error /= OPSUCCESS)GOTO 200
    
    CALL copy_solution(v,sol)

    200 CONTINUE
    IF(error /= OPSUCCESS .OR. error == SOLVER_FAILED)THEN
      WRITE(*,*)'GPU--solver problem.'
      PAUSE
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
                          WORK, Acsr,Aclc,Arwc,        &
	                           NNZ,  Alu,Jlu,Ju,Jw          )

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
    !--CELL-CENTER GRADIENT 
    
    vele=u(ij+nj)*fx(i)+u(ij)*(one-fx(i))
    velw=u(ij)*fx(i-1)+u(ij-nj)*(one-fx(i-1))
    veln=u(ij+1)*fy(j)+u(ij)*(one-fy(j))
    vels=u(ij)*fy(j-1)+u(ij-1)*(one-fy(j-1))
    dux(ij)=(vele-velw)/dx
    duy(ij)=(veln-vels)/dy

    vele=v(ij+nj)*fx(i)+v(ij)*(one-fx(i))
    velw=v(ij)*fx(i-1)+v(ij-nj)*(one-fx(i-1))
    veln=v(ij+1)*fy(j)+v(ij)*(one-fy(j))
    vels=v(ij)*fy(j-1)+v(ij-1)*(one-fy(j-1))
    dvx(ij)=(vele-velw)/dx
    dvy(ij)=(veln-vels)/dy

  ENDDO
ENDDO

!IF(temp_visc)CALL fd_calc_bcuv_grad

END SUBROUTINE fd_calc_mom

END MODULE modfd_calc_mom
