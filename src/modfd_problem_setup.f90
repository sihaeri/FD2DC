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

MODULE modfd_problem_Setup

PRIVATE
PUBLIC :: fd_problem_setup,fd_problem_restart,fd_copy_time_arrays,fd_print_var,fd_alloc_surfforce_arrays, fd_write_restart,&
          fd_update_visc,fd_calc_timeCoefs
CONTAINS
!============================================
SUBROUTINE fd_problem_setup

USE real_parameters,        ONLY : zero,one,two,three,four,pi,half,quarter,fifth,vlarge,vsmall,oneandhalf,fourthird,small
USE parameters,             ONLY : out_unit,set_unit,alloc_create,nphi,grd_unit,max_char_len,solver_sparsekit,solver_sip,&
                                   solver_cg,solver_hypre,plt_unit,init_field_unit,use_GPU_yes,use_GPU_no,&
                                   OPSUCCESS
USE precision,              ONLY : r_single
USE shared_data,            ONLY : devID,solver_type,NNZ,NCel,lli,itim,time,lread,lwrite,ltest,louts,loute,ltime,&
                                   maxit,imon,jmon,ipr,jpr,sormax,slarge,alfa,minit,initFromFile,&
                                   densit,visc,gravx,gravy,beta,tref,problem_name,problem_len,lamvisc,&
                                   ulid,tper,itst,nprt,dt,dto,gamt,iu,iv,ip,ien,tper,f1,ft1,f2,voo,uoo,too,&
                                   nsw,lcal,sor,resor,urf,gds,om,cpf,kappaf,nim,njm,ni,nj,calcwalnusselt,&
                                   calcwalnusselt_ave,naverage_wsteps,xPeriodic,yPeriodic,Ldomainx,Ldomainy,&
                                   li,x,y,xc,yc,fx,fy,laxis,r,nij,u,v,nsw,sor,lcal,p,t,th,tc,den,deno,denoo,dxmean,&
                                   vo,uo,to,title,duct,flomas,flomom,f1,stationary,celbeta,celkappa,celcp,celcpo,celcpoo,&  !--FD aspect
                                   objcentx,objcenty,objcentxo,objcentyo,objradius,putobj,nsurfpoints,read_fd_geom,&
                                   betap,movingmesh,forcedmotion,up,vp,omp,isotherm,objqp,objbradius,&
                                   densitp,fd_urf,mcellpercv,objtp,nsphere,calcsurfforce,rigidForce_contrib,& !--heat
                                   nnusseltpoints,calclocalnusselt,deltalen,sphnprt,betap,kappap,cpp,&
                                   mpi_comm,Hypre_A,Hypre_b,Hypre_x,nprocs_mpi,myrank_mpi,& !--viscosity temperature dependance
                                   temp_visc,viscgamma,viscN,calclocalnusselt_ave,naverage_steps,&
                                   use_GPU,arow,acol,acoo,acsr,aclc,arwc,ndt,&!--DEM Shared data
                                   dem_kn,dem_kt,dem_gamman,dem_gammat,dem_mut,dem_dt,subTimeStep
!--CULA Variables handle, platformOpts, formatOpts,  precondOpts, solverOpts, config

USE modfd_set_bc,           ONLY : fd_bctime
USE modfd_create_geom,      ONLY : fd_create_geom,fd_calc_mi,fd_calc_physprops,fd_calc_physprops,&
                                   fd_calc_quality,fd_init_temp,fd_find_rigidforce_contrib
USE modfd_tecwrite,         ONLY : fd_tecwrite_sph_s,fd_tecwrite_sph_v,fd_tecwrite_eul

USE modcusp_library_intrf,  ONLY : getInstance_cusp_biCGSTAB_solver,&
     cusp_biCGSTAB_initDevice,cusp_biCGSTAB_allocDevice,cusp_biCGSTAB_copyH2D_AInds, &
     cusp_biCGSTAB_shutdown

USE modfd_solve_linearsys,  ONLY : fd_cooInd_create
IMPLICIT NONE

REAL(KIND = r_single) :: uin,vin,pin,tin,dx,dy,rp,domain_vol,rdummy,massi,massj,mstar
INTEGER               :: i,j,ij,ncell,idummy,error
!REAL(KIND = r_single),ALLOCATABLE :: volp(:,:),q(:)

error = OPSUCCESS
solver_type = solver_sparsekit
itim = 0
time = zero
naverage_steps = 0

CALL fd_alloc_solctrl_arrays(alloc_create)

READ(set_unit,6)title
READ(set_unit,*)xPeriodic,yPeriodic,use_GPU,devID

IF(use_GPU /= use_GPU_no .AND. use_GPU /= use_GPU_yes)THEN
  WRITE(*,'(A, I5, A, I5)')'use_gpu should be:', use_GPU_yes,', or:', use_GPU_no
  STOP
ENDIF

IF((xPeriodic /= 0 .AND. xPeriodic /= 1) .OR. (yPeriodic /= 0 .AND. yPeriodic /= 1))THEN
  WRITE(*,*)'Periodic boundary controller variables should be 0 or 1.'
  STOP
ENDIF

IF((xPeriodic == 1 .OR. yPeriodic == 1) .AND. solver_type /= solver_sparsekit)THEN
  WRITE(*,*)'Periodic boundary controller only works with BiCGSTAB Solver set solver_type = solver_sparsekit.'
  STOP
ENDIF

READ(set_unit,*)lread,initFromFile,lwrite,ltest,laxis,louts,loute,ltime,duct
READ(set_unit,*)maxit,minit,imon,jmon,ipr,jpr,sormax,slarge,alfa
READ(set_unit,*)densit,visc,cpf,kappaf,gravx,gravy,beta,th,tc,tref
READ(set_unit,*)uin,vin,pin,tin,ulid,ndt,tper
READ(set_unit,*)itst,nprt,dt,gamt

dto = dt
IF(gamt /= 1 .AND. gamt /= 2)THEN
  WRITE(out_unit,*)'GAMT:: Is the time integration order = 1 or 2.'
  STOP
ENDIF

READ(set_unit,*)(lcal(i),i=1,nphi)
READ(set_unit,*)(urf(i),i=1,nphi)
READ(set_unit,*)(sor(i),i=1,nphi)
READ(set_unit,*)(nsw(i),i=1,nphi)
READ(set_unit,*)(gds(i),i=1,nphi)

READ(set_unit,*)temp_visc
IF(temp_visc)READ(set_unit,*)viscN
READ(set_unit,*)calcwalnusselt
IF(calcwalnusselt)READ(set_unit,*)calcwalnusselt_ave,naverage_wsteps
READ(set_unit,*)putobj,read_fd_geom,stationary,forcedmotion,movingmesh,calcsurfforce,calclocalnusselt,isotherm
IF(calclocalnusselt)READ(set_unit,*)calclocalnusselt_ave,naverage_steps
IF(putobj)THEN
  READ(set_unit,*)deltalen
  READ(set_unit,*)fd_urf
  READ(set_unit,*)nsphere
  IF(nsphere > 0)THEN
    READ(set_unit,*)dem_kn,dem_kt,dem_gamman,dem_gammat,dem_mut

    READ(set_unit,*)sphnprt
    CALL fd_alloc_objprop_arrays(alloc_create)
    DO i = 1,nsphere
      IF(isotherm)THEN
        READ(set_unit,*)objcentx(i),objcenty(i),objradius(i),objbradius(i),densitp(i),objtp(i),betap(i),cpp(i),kappap(i)

      ELSE
        READ(set_unit,*)objcentx(i),objcenty(i),objradius(i),objbradius(i),densitp(i),objtp(i),betap(i),cpp(i),kappap(i),&
                        objqp(i)

      ENDIF
      READ(set_unit,*)nsurfpoints(i),nnusseltpoints(i),mcellpercv(i)
      objcentxo = objcentx
      objcentyo = objcenty
    ENDDO
    !--Estimate the collision time for the spring-dashpot model
    dem_dt = four/three*pi*(SUM(densitp)/nsphere)*(SUM(objradius)/nsphere)**3
    dem_dt = fifth*pi/SQRT(two*dem_kn/dem_dt-(dem_gamman/two)**2)
    subTimeStep = MAX(1,FLOOR(dt/dem_dt))
    WRITE(*,*) 'DEM dt= ', dem_dt, ' , flow dt= ', dt, ', sub-steps= ',subTimeStep
    IF(.NOT. lread)THEN
      CALL fd_alloc_objgeom_arrays(alloc_create,MAXVAL(nsurfpoints(:),1))
      IF(calclocalnusselt)CALL fd_alloc_nusselt_arrays(alloc_create,MAXVAL(nnusseltpoints(:),1))
      IF(calcsurfforce)CALL fd_alloc_surfforce_arrays(alloc_create,MAXVAL(nsurfpoints(:),1))
    ENDIF
  ENDIF
ENDIF

!-- Init solved variable index
iu = 1
iv = 2
ip = 3
ien = 4


om = two*pi/tper !--omega (for oscilatory lid)
!prr = one/prm !--inverse prandtl !not needed anymore

READ(grd_unit,*) i
READ(grd_unit,*) i
READ(grd_unit,*) ni
READ(grd_unit,*) nj
READ(grd_unit,*) ij

nim=ni-1
njm=nj-1
nij=ni*nj

CALL fd_alloc_geom_arrays(alloc_create)

READ(grd_unit,*) (x(i),i=1,ni)
READ(grd_unit,*) (y(j),j=1,nj)
DO i=1,ni
  li(i)=(i-1)*nj
END DO

!-----number of s-n faces + number of e-w faces + central coefficients
NCel = 0
DO i=2,nim
  DO j=2,njm
    NCel = NCel + 1
    lli(li(i) + j) = NCel
  ENDDO
ENDDO

NNZ = 2*(njm - 2 + yPeriodic)*(nim - 1) + 2*(nim - 2 + xPeriodic)*(njm - 1) + NCel
IF(solver_type == solver_sparsekit .OR. solver_type == solver_hypre)THEN

  CALL fd_alloc_spkit_arrays(alloc_create)
  CALL fd_cooInd_create() !--Create coordinate arrays Arow, Acol

  IF(use_GPU == use_GPU_yes)THEN
    WRITE(*,*)'Initializing GPU using ...'
    CALL getInstance_cusp_biCGSTAB_solver()

    CALL cusp_biCGSTAB_initDevice(devID,error)
    IF(error /= OPSUCCESS)GOTO 100

    CALL cusp_biCGSTAB_allocDevice(NCel,NNZ,error)
    IF(error /= OPSUCCESS)GOTO 100
    Arow = Arow - 1
    Acol = Acol - 1
    CALL cusp_biCGSTAB_copyH2D_AInds(Arow,Acol,error)
    IF(error /= OPSUCCESS)GOTO 100

!     IF(culaSparseCreate(handle) /= culaSparseNoError)GOTO 100
!     IF(culaSparseCudaOptionsInit(handle, platformOpts) /= culaSparseNoError)GOTO 100
!     platformOpts%deviceId=devId
!     platformOpts%useHybridFormat=1
!    IF(culaSparseCooOptionsInit(handle, formatOpts) /= culaSparseNoError)GOTO 100
!    formatOpts%indexing = 1
!    IF(culaSparseJacobiOptionsInit(handle, precondOpts) /= culaSparseNoError)GOTO 100
    !IF(culaSparseEmptyOptionsInit(handle, precondOpts) /= culaSparseNoError)GOTO 100
!    IF(culaSparseConfigInit(handle, config) /= culaSparseNoError)GOTO 100
!    IF(culaSparseBicgstabOptionsInit(handle, solverOpts) /= culaSparseNoError)GOTO 100
    WRITE(*,*)'GPU Initialized.'
  ENDIF

  100 CONTINUE
  IF(error /= OPSUCCESS)THEN
    WRITE(*,*)'Unable to start the GPU. Running the solve on the CPU instead...'
    use_GPU = use_GPU_no
    CALL cusp_biCGSTAB_shutdown(error)
  ENDIF
ENDIF

!--cv centres
DO i=2,nim
  xc(i)=half*(x(i)+x(i-1))
END DO

xc(1)=x(1) - xPeriodic*(x(nim) - xc(nim))
xc(ni)=x(nim) + xPeriodic*(xc(2) - x(1))

DO j=2,njm
  yc(j)=half*(y(j)+y(j-1))
END DO

yc(1)=y(1) - yPeriodic*(y(njm) - yc(njm))
yc(nj)=y(njm) + yPeriodic*(yc(2) - y(1))

Ldomainx = x(nim) - x(1)
Ldomainy = y(njm) - y(1)

!--interpolation factors
DO i=1,nim
  fx(i)=(x(i)-xc(i))/(xc(i+1)-xc(i))
END DO

DO j=1,njm
  fy(j)=(y(j)-yc(j))/(yc(j+1)-yc(j))
END DO

!--radius (significant only if axi-symmetric)
IF(laxis) THEN
  DO j=1,nj
    r(j)=y(j)
  END DO
ELSE
  DO j=1,nj
    r(j)=one
  END DO
ENDIF

!--Setup the object
IF(putobj)THEN
  IF(.NOT.lread)THEN
    IF(nsphere > 0)CALL fd_create_geom
    IF(nsphere > 0)CALL fd_alloc_objcell_vars(alloc_create)
    CALL fd_alloc_sources(alloc_create)
  ENDIF
  CALL fd_calc_mi

  ncell = 0
  domain_vol = zero
  DO i=2,nim
    dx=x(i)-x(i-1)
    DO j=2,njm
      ncell = ncell + 1
      dy=y(j)-y(j-1)
      rp=half*(r(j)+r(j-1))
      domain_vol = domain_vol + (dx*dy*rp)**half
    ENDDO
  ENDDO

  dxmean = domain_vol/REAL(ncell,r_single)
ENDIF

!---------------------------------------------------
!--BOUNDARY AND INITIAL CONDITIONS
!---------------------------------------------------
!--Allocate working arrays (matrix diagnals, sources etc...)
CALl fd_alloc_work_arrays(alloc_create)

!--WEST AND EAST ISOTHERMAL BOUNDARIES
!
IF(.NOT.movingmesh)THEN
  DO j=1,nj
      t(j)=th
  END DO

  DO j=1,nj
    t(li(ni)+j)=tc
  END DO
ELSE
  DO j=1,nj
    t(j)=th
  END DO

  DO i=1,ni
    t(li(i)+1)=th
  END DO

  DO i=1,ni
    t(li(i)+nj)=th
  END DO

ENDIF

!--NORTH WALL VELOCITY (FOR LID-DRIVEN CAVITY)
IF(ltime) THEN
  CALL fd_bctime(0)
ELSE
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
    u(li(i)+nj)=ulid
  END DO
  DO i=2,nim
    u(li(i)+1)=zero
  END DO
ENDIF
ENDIF


!--INITIAL VARIBLE VALUES (INITIAL CONDITIONS)
IF(.NOT. initFromFile)THEN
  DO i=2,nim
    DO ij=li(i)+2,li(i)+njm
      u(ij)=uin
      v(ij)=vin
      t(ij)=tin
      p(ij)=pin
      uo(ij)=uin
      vo(ij)=vin
      to(ij)=tin
    END DO
  END DO
  den  = densit
  deno = densit
  denoo= densit
  celbeta = beta*densit
  celcp = cpf*densit
  celcpo = cpf*densit
  celcpoo= cpf*densit
  celkappa = kappaf
  lamvisc = visc
ELSE
  OPEN(UNIT = init_field_unit,FILE=problem_name(1:problem_len)//'.ini',STATUS='OLD',FORM='UNFORMATTED')
    READ(init_field_unit) idummy,rdummy,idummy,idummy,idummy,idummy,idummy,dto,&
         ((x(i),j=1,nj),i=1,ni),((y(j),j=1,nj),i=1,ni),&
         ((xc(i),j=1,nj),i=1,ni),((yc(j),j=1,nj),i=1,ni),&
         (f1(ij),ij=1,nij),(f2(ij),ij=1,nij),&
         (u(ij),ij=1,nij),(v(ij),ij=1,nij),(p(ij),ij=1,nij),(t(ij),ij=1,nij),&
         (uo(ij),ij=1,nij),(vo(ij),ij=1,nij),(to(ij),ij=1,nij),&
         (uoo(ij),ij=1,nij),(voo(ij),ij=1,nij),(too(ij),ij=1,nij),&
         (den(ij),ij=1,nij),(deno(ij),ij=1,nij),(denoo(ij),ij=1,nij),(celbeta(ij),ij=1,nij),&
         (celcp(ij),ij=1,nij),(celcpo(ij),ij=1,nij),(celcpoo(ij),ij=1,nij),(celkappa(ij),ij=1,nij)
  CLOSE(init_field_unit)
ENDIF

IF(putobj)THEN
  IF(nsphere > 0)THEN
    !--Mark affected cells
    CALL fd_find_rigidforce_contrib
    CALL fd_calc_physprops(1)
    CALL fd_init_temp(tin)
  ENDIF
  deno = den
  denoo = deno

  celcpo  = celcp
  celcpoo = celcpo

  to = t
  too = to

ENDIF

CALL fd_print_problem_setup

OPEN(UNIT = plt_unit,FILE=problem_name(1:problem_len)//'_init.plt',STATUS='NEW')
IF(temp_visc)THEN
  CALL fd_tecwrite_eul(plt_unit,lamvisc)
ELSE
  CALL fd_tecwrite_eul(plt_unit)
ENDIF

CLOSE(plt_unit)
!ALLOCATE(volp(nij,nsphere),q(nsphere))
!
!volp = zero
!CALL fd_calc_geom_volf(volp)
!OPEN(UNIT = plt_unit,FILE=problem_name(1:problem_len)//'_init'//'_geomvolp.plt',STATUS='NEW')
!CALL fd_tecwrite_eul(plt_unit,volp)
!CLOSE(plt_unit)
!
!volp = zero
!CALL fd_calc_physprops(volp)
!OPEN(UNIT = plt_unit,FILE=problem_name(1:problem_len)//'_init'//'_volp.plt',STATUS='NEW')
!CALL fd_tecwrite_eul(plt_unit,volp)
!CLOSE(plt_unit)
!CALL fd_calc_quality(volp,q)
!WRITE(*,*)q
!DEALLOCATE(volp,q)

6 FORMAT(A80)

END SUBROUTINE fd_problem_setup

SUBROUTINE fd_problem_restart

USE real_parameters,        ONLY : zero
USE parameters,             ONLY : sres_unit,alloc_create,out_unit
USE shared_data

IMPLICIT NONE

INTEGER :: lni,lnj,lnim,lnjm,lnij,ij,i,j,lnsphere,nn,ik,maxnobjcell,lxPeriodic,&
           lyPeriodic
LOGICAL :: lputobj,lread_fd_geom,lstationary,lforcedmotion,lmovingmesh,&
                  lcalcsurfforce,lcalclocalnusselt,lisotherm,ltemp_visc
IF(putobj)THEN
  READ(sres_unit) lxPeriodic,lyPeriodic
  READ(sres_unit) itim,time,lni,lnj,lnim,lnjm,lnij,dto
  IF(lni /= ni .OR. lnj /= nj .OR. lnij /= nij .OR. lnim /= nim .OR. lnjm /= njm .OR. &
     lxPeriodic /= xPeriodic .OR. lyPeriodic /= yPeriodic)THEN
      WRITE(out_unit,*)'fd_problem_restart: restart file inconsistency.'
      STOP
  ENDIF

  CALL  fd_alloc_sources(alloc_create)

  READ(sres_unit)((x(i),j=1,nj),i=1,ni),((y(j),j=1,nj),i=1,ni),&
         ((xc(i),j=1,nj),i=1,ni),((yc(j),j=1,nj),i=1,ni),&
         (f1(ij),ij=1,nij),(f2(ij),ij=1,nij),&
         (ft1(ij),ij=1,nij),(ft2(ij),ij=1,nij),&
         (u(ij),ij=1,nij),&
         (v(ij),ij=1,nij),(p(ij),ij=1,nij),(t(ij),ij=1,nij),&
         (uo(ij),ij=1,nij),(vo(ij),ij=1,nij),(to(ij),ij=1,nij),&
         (uoo(ij),ij=1,nij),(voo(ij),ij=1,nij),(too(ij),ij=1,nij),&
         (dpx(ij),ij=1,nij),(dpy(ij),ij=1,nij),(dux(ij),ij=1,nij),&
         (duy(ij),ij=1,nij),(dvx(ij),ij=1,nij),(dvy(ij),ij=1,nij),&
         (dtx(ij),ij=1,nij),(dty(ij),ij=1,nij),(den(ij),ij=1,nij),(deno(ij),ij=1,nij),&
         (denoo(ij),ij=1,nij),(celbeta(ij),ij=1,nij),(celcp(ij),ij=1,nij),&
         (celcpo(ij),ij=1,nij),(celcpoo(ij),ij=1,nij),&
         (celkappa(ij),ij=1,nij),(fdsu(ij),ij=1,nij),(fdsv(ij),ij=1,nij),(fdsub(ij),ij=1,nij),&
         (fdsvb(ij),ij=1,nij),(fdsuc(ij),ij=1,nij),(fdsvc(ij),ij=1,nij),&
         (fdst(ij),ij=1,nij),(fdstc(ij),ij=1,nij),(lamvisc(ij),ij=1,nij),&
         (fdfcu(ij),ij=1,nij),(fdfcv(ij),ij=1,nij)

  READ(sres_unit)ltemp_visc
  IF(ltemp_visc .NEQV. temp_visc)THEN
    WRITE(out_unit,*)'fd_problem_restart: restart file inconsistency in viscosity dependant switches.'
    STOP
  ENDIF

  READ(sres_unit)lputobj,lread_fd_geom,lstationary,lforcedmotion,lmovingmesh,&
                  lcalcsurfforce,lcalclocalnusselt,lisotherm

  IF(lputobj .NEQV. putobj .OR. lread_fd_geom .NEQV. read_fd_geom .OR. lstationary .NEQV. stationary .OR. &
     lforcedmotion .NEQV. forcedmotion .OR. lmovingmesh .NEQV. movingmesh .OR. lcalcsurfforce .NEQV. calcsurfforce .OR. &
     lcalclocalnusselt .NEQV. calclocalnusselt .OR. lisotherm .NEQV. isotherm)THEN
      WRITE(out_unit,*)'fd_problem_restart: restart file inconsistency in lagrangian switches.'
      STOP
  ENDIF

  READ(sres_unit)lnsphere

  IF(lnsphere /= nsphere)THEN
      WRITE(out_unit,*)'fd_problem_restart: restart file inconsistency in the number of objects.'
      STOP
  ENDIF

  READ(sres_unit)(objtp(ij),ij=1,nsphere),(densitp(ij),ij=1,nsphere),(objcentx(ij),ij=1,nsphere),&
                  (objcenty(ij),ij=1,nsphere),(objradius(ij),ij=1,nsphere),(objbradius(ij),ij=1,nsphere),&
                  (objcentu(ij),ij=1,nsphere),&
                  (objcentv(ij),ij=1,nsphere),(objcentom(ij),ij=1,nsphere),(nsurfpoints(ij),ij=1,nsphere),&
                  (nobjcells(ij),ij=1,nsphere),(mcellpercv(ij),ij=1,nsphere),(nnusseltpoints(ij),ij=1,nsphere),&
                  (objcentmi(ij),ij=1,nsphere),(objvol(ij),ij=1,nsphere),&
                  (objcentxo(ij),ij=1,nsphere),(objcentyo(ij),ij=1,nsphere),(objcentvo(ij),ij=1,nsphere),&
                  (objcentuo(ij),ij=1,nsphere),(objcentomo(ij),ij=1,nsphere),(betap(ij),ij=1,nsphere),&
                  (cpp(ij),ij=1,nsphere),(kappap(ij),ij=1,nsphere),(Fpq(1,ij),ij=1,nsphere),&
                  (Fpq(2,ij),ij=1,nsphere),(zobjcentx(ij),ij=1,nsphere),(zobjcenty(ij),ij=1,nsphere),&
                  (zobjcentxo(ij),ij=1,nsphere),(zobjcentyo(ij),ij=1,nsphere)
  DO nn=1,nsphere
    READ(sres_unit)(objcento(ij,nn),ij=1,4)
  ENDDO

  CALL fd_alloc_objgeom_arrays(alloc_create,MAXVAL(nsurfpoints(:),1))

  IF(forcedmotion)THEN
    DO nn = 1,nsphere
      READ(sres_unit)(surfpointx(ij,nn),ij=1,nsurfpoints(nn)),&
                      (surfpointy(ij,nn),ij=1,nsurfpoints(nn)),&
                      (surfpointxinit(ij,nn),ij=1,nsurfpoints(nn)),&
                      (surfpointyinit(ij,nn),ij=1,nsurfpoints(nn)),&
                      (zsurfpointx(ij,nn),ij=1,nsurfpoints(nn)),&
                      (zsurfpointy(ij,nn),ij=1,nsurfpoints(nn))
    ENDDO
  ELSE
    DO nn = 1,nsphere
      READ(sres_unit)(surfpointx(ij,nn),ij=1,nsurfpoints(nn)),&
                      (surfpointy(ij,nn),ij=1,nsurfpoints(nn)),&
                      (zsurfpointx(ij,nn),ij=1,nsurfpoints(nn)),&
                      (zsurfpointy(ij,nn),ij=1,nsurfpoints(nn))
    ENDDO
  ENDIF

  IF(calclocalnusselt)THEN
    CALL fd_alloc_nusselt_arrays(alloc_create,MAXVAL(nnusseltpoints(:),1))
    DO nn=1,nsphere
      READ(sres_unit)(nusseltpointx(ij,nn),ij=1,nnusseltpoints(nn)),&
                      (nusseltpointy(ij,nn),ij=1,nnusseltpoints(nn)),&
                      (nusseltnx(ij,nn),ij=1,nnusseltpoints(nn)),&
                      (nusseltny(ij,nn),ij=1,nnusseltpoints(nn)),&
                      (nusseltds(ij,nn),ij=1,nnusseltpoints(nn)),&
                      (nusseltcentx(ij,nn),ij=1,nnusseltpoints(nn)),&
                      (nusseltcenty(ij,nn),ij=1,nnusseltpoints(nn)),&
                      (localnusselt(ij,nn),ij=1,nnusseltpoints(nn))
    ENDDO
    DO nn=1,nsphere
      DO ij=1,nnusseltpoints(nn)
        READ(sres_unit)(nusseltpoint_cvy(i,ij,nn),i=1,3),&
                        (nusseltpoint_cvx(i,ij,nn),i=1,3),&
                        (nusseltinterpx(i,ij,nn),i=1,2),&
                        (nusseltinterpy(i,ij,nn),i=1,2)
        DO ik=1,2
           READ(sres_unit)(nusseltpoint_interp(i,ik,ij,nn),i=1,4)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF(calcsurfforce)THEN
    CALL fd_alloc_surfforce_arrays(alloc_create,MAXVAL(nsurfpoints(:),1))
    DO nn = 1,nsphere
      READ(sres_unit)(surfds(ij,nn),ij=1,nsurfpoints(nn)),&
                      (surfnx(ij,nn),ij=1,nsurfpoints(nn)),&
                      (surfny(ij,nn),ij=1,nsurfpoints(nn)),&
                      (surfcentrex(ij,nn),ij=1,nsurfpoints(nn)),&
                      (surfcentrey(ij,nn),ij=1,nsurfpoints(nn))
      DO ij=1,nsurfpoints(nn)
        READ(sres_unit)(surfpoint_cvx(i,ij,nn),i=1,3),&
                        (surfpoint_cvy(i,ij,nn),i=1,3)
        DO ik = 1,2
          READ(sres_unit)(surfpoint_interp(i,ik,ij,nn),i=1,4),surfinterpx(ik,ij,nn),surfinterpy(ik,ij,nn)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  maxnobjcell = MAXVAL(nobjcells(:),1)
  CALL fd_alloc_objcell_vars(alloc_create)
  ALLOCATE(objcellx(maxnobjcell,nsphere),objcelly(maxnobjcell,nsphere),zobjcellx(maxnobjcell,nsphere),&
           zobjcelly(maxnobjcell,nsphere),objcellvol(maxnobjcell,nsphere),objpoint_cvx(maxnobjcell,nsphere),&
           objpoint_cvy(maxnobjcell,nsphere),objcellvertx(4,maxnobjcell,nsphere),objcellverty(4,maxnobjcell,nsphere),&
           zobjcellvertx(4,maxnobjcell,nsphere),zobjcellverty(4,maxnobjcell,nsphere),&
           objpoint_interpx(2,maxnobjcell,nsphere),objpoint_interpy(2,maxnobjcell,nsphere))
  objcellx = zero;objcelly = zero
  objcellvol = zero;objcellvertx = zero;objcellverty = zero
  objpoint_cvx = 0;objpoint_cvy = 0;objpoint_interpx = 0;objpoint_interpy = 0
  zobjcellx = 0; zobjcelly = 0
  zobjcellvertx = 0; zobjcellverty = 0

  DO nn=1,nsphere
    READ(sres_unit)(objcellx(ij,nn),ij=1,nobjcells(nn)),&
                    (objcelly(ij,nn),ij=1,nobjcells(nn)),&
                    (objcellvol(ij,nn),ij=1,nobjcells(nn)),&
                    (objpoint_cvx(ij,nn),ij=1,nobjcells(nn)),&
                    (objpoint_cvy(ij,nn),ij=1,nobjcells(nn)),&
                    (zobjcellx(ij,nn),ij=1,nobjcells(nn)),&
                    (zobjcelly(ij,nn),ij=1,nobjcells(nn))

    DO ij=1,nobjcells(nn)
      READ(sres_unit)(objcellvertx(i,ij,nn),i=1,4),&
                      (objcellverty(i,ij,nn),i=1,4),&
                      (objpoint_interpx(i,ij,nn),i=1,2),&
                      (objpoint_interpy(i,ij,nn),i=1,2),&
                      (zobjcellvertx(i,ij,nn),i=1,4),&
                      (zobjcellverty(i,ij,nn),i=1,4)
    ENDDO

  ENDDO

  IF(forcedmotion)THEN
    ALLOCATE(objcentxinit(nsphere),objcentyinit(nsphere),objcellxinit(maxnobjcell,nsphere),&
           objcellyinit(maxnobjcell,nsphere),objcellvertxinit(4,maxnobjcell,nsphere),&
           objcellvertyinit(4,maxnobjcell,nsphere))
    objcentxinit = zero;objcentyinit = zero;objcellxinit = zero
    objcellyinit = zero;objcellvertxinit = zero
    objcellvertyinit = zero

    READ(sres_unit)(objcentxinit(nn),nn=1,nsphere),(objcentyinit(nn),nn=1,nsphere)
    DO nn=1,nsphere
      READ(sres_unit)(objcellxinit(ij,nn),ij=1,nobjcells(nn)),&
                      (objcellyinit(ij,nn),ij=1,nobjcells(nn))
      DO ij=1,nobjcells(nn)
        READ(sres_unit)(objcellvertxinit(i,ij,nn),i=1,4),&
                        (objcellvertyinit(i,ij,nn),i=1,4)
      ENDDO

    ENDDO

  ENDIF

ELSE
  READ(sres_unit) itim,time,lni,lnj,lnim,lnjm,lnij
  IF(lni /= ni .OR. lnj /= nj .OR. lnij /= nij .OR. lnim /= nim .OR. lnjm /= njm)THEN
      WRITE(out_unit,*)'fd_problem_restart: restart file inconsistency.'
      STOP
  ENDIF
  READ(sres_unit)((x(i),j=1,nj),i=1,ni),((y(j),j=1,nj),i=1,ni),&
      ((xc(i),j=1,nj),i=1,ni),((yc(j),j=1,nj),i=1,ni),&
      (f1(ij),ij=1,nij),(f2(ij),ij=1,nij),&
      (ft1(ij),ij=1,nij),(ft2(ij),ij=1,nij),&
      (u(ij),ij=1,nij),&
      (v(ij),ij=1,nij),(p(ij),ij=1,nij),(t(ij),ij=1,nij),&
      (uo(ij),ij=1,nij),(vo(ij),ij=1,nij),(to(ij),ij=1,nij)
  REWIND sres_unit
ENDIF
END SUBROUTINE fd_problem_restart

SUBROUTINE fd_alloc_solctrl_arrays(create_or_destroy)

USE real_parameters,        ONLY : zero
USE parameters,             ONLY : alloc_create,alloc_destroy,nphi
USE shared_data,            ONLY : nsw,lcal,sor,resor,urf,gds

IMPLICIT NONE
INTEGER,INTENT(IN)      :: create_or_destroy

IF(create_or_destroy == alloc_create)THEN
  ALLOCATE(nsw(nphi),lcal(nphi),sor(nphi),resor(nphi),urf(nphi),gds(nphi))
  urf = zero
  gds = zero
  resor = zero
  sor = zero
  nsw = 0
  lcal = .FALSE.
ELSEIF(create_or_destroy == alloc_destroy)THEN
  DEALLOCATE(nsw,lcal,sor,resor,urf,gds)
ENDIF

END SUBROUTINE fd_alloc_solctrl_arrays

SUBROUTINE fd_alloc_geom_arrays(create_or_destroy)

USE parameters,             ONLY : alloc_create,alloc_destroy,nphi
USE shared_data,            ONLY : r,x,y,xc,yc,fx,fy,li,nj,ni,lli,nij
USE real_parameters,        ONLY : zero

IMPLICIT NONE
INTEGER,INTENT(IN)      :: create_or_destroy

IF(create_or_destroy == alloc_create)THEN
  ALLOCATE(x(ni),y(nj),r(nj),xc(ni),yc(nj),fx(ni-1),fy(nj-1),li(ni),lli(nij))
  x = zero
  y = zero
  yc = zero
  xc = zero
  fx = zero
  fy = zero
  r = zero
  lli = -1
ELSEIF(create_or_destroy == alloc_destroy)THEN
  DEALLOCATE(x,y,r,xc,yc,fx,fy,li,lli)
ENDIF

END SUBROUTINE fd_alloc_geom_arrays

SUBROUTINE fd_alloc_spkit_arrays(create_or_destroy)

USE parameters,             ONLY : alloc_create,alloc_destroy,nphi
USE shared_data,            ONLY : Work,Acoo,Arow,Acol,Acsr,Arwc,Aclc,RHS,SOL,NCel,NNZ,&
                                   alu,jlu,ju,jw
USE real_parameters,        ONLY : zero

IMPLICIT NONE
INTEGER,INTENT(IN)      :: create_or_destroy

IF(create_or_destroy == alloc_create)THEN
  ALLOCATE(Work(NNZ,8),Acoo(NNZ),Arow(NNZ),Acol(NNZ),Acsr(NNZ),Arwc(NCel+1),Aclc(NNZ),RHS(NCel),SOL(NCel),&
           Alu(NNZ),jlu(NNZ),Ju(NCel),Jw(2*NCel))
  Work = 0.D0
  Acoo = 0.D0
  Acsr = 0.D0
  RHS  = 0.D0
  SOL  = 0.D0
  Alu  = 0.D0
  Arow = 0
  Acol = 0
  Arwc = 0
  Aclc = 0
  jlu  = 0
  Ju   = 0
  Jw   = 0
ELSEIF(create_or_destroy == alloc_destroy)THEN
  DEALLOCATE(Work,Acoo,Arow,Acol,Acsr,Arwc,Aclc,RHS,SOL,Alu,jlu,Ju,Jw)
ENDIF

END SUBROUTINE fd_alloc_spkit_arrays


SUBROUTINE fd_print_problem_setup

USE real_parameters, ONLY : zero
USE shared_data,  ONLY : ulid,gravx,gravy,ni,li,t,alfa,urf,gds,dt,gamt,tper,densit,visc,ien,&
                         ltime,iu,iv,ip,lcal,cpf,kappaf
USE parameters,   ONLY : out_unit
IMPLICIT NONE

WRITE(out_unit,601) densit,visc

IF(ULID /= zero) THEN
  WRITE(out_unit,*) '          MAX. LID VELOCITY:  ',ULID
ENDIF
IF(LCAL(IEN)) THEN
  WRITE(out_unit,*) '          GRAVITY IN X-DIR.:  ',GRAVX
  WRITE(out_unit,*) '          GRAVITY IN Y-DIR.:  ',GRAVY
  WRITE(out_unit,*) '          HOT  WALL TEMPER.:  ',t(2)
  WRITE(out_unit,*) '          COLD WALL TEMPER.:  ',t(li(ni)+1)
  WRITE(out_unit,*) '          Cpf, Kappaf      :  ',cpf, kappaf
ENDIF
WRITE(out_unit,*) '  '
WRITE(out_unit,*) '          ALFA  PARAMETER  :  ',ALFA
WRITE(out_unit,*) '  '
WRITE(out_unit,*) '          UNDERRELAXATION  FACTORS'
WRITE(out_unit,*) '          ========================'
WRITE(out_unit,*) '          U-VELOCITY  :  ',URF(IU)
WRITE(out_unit,*) '          V-VELOCITY  :  ',URF(IV)
WRITE(out_unit,*) '          PRESSURE    :  ',URF(IP)
WRITE(out_unit,*) '          TEMPERATURE :  ',URF(IEN)
WRITE(out_unit,*) '  '
WRITE(out_unit,*) '          SPATIAL BLENDING FACTORS (CDS-UDS)'
WRITE(out_unit,*) '          =================================='
WRITE(out_unit,*) '          U-VELOCITY  :  ',GDS(IU)
WRITE(out_unit,*) '          V-VELOCITY  :  ',GDS(IV)
WRITE(out_unit,*) '          TEMPERATURE :  ',GDS(IEN)
WRITE(out_unit,*) '  '
IF(LTIME) THEN
  WRITE(out_unit,*) '          UNSTEADY FLOW SIMULATION'
  WRITE(out_unit,*) '          ================================='
  WRITE(out_unit,*) '          TIME STEP SIZE       : ',DT
  WRITE(out_unit,*) '          TIME STEP ORDER      : ',GAMT
  WRITE(out_unit,*) '          OSCILLATION PERIOD   : ',TPER
ENDIF
WRITE(out_unit,*) '  '
WRITE(out_unit,*) '  '
RETURN

601 FORMAT(10X,' FLUID DENSITY    :  ',1P1E12.4,/,10X,' DYNAMIC VISCOSITY:  ',1P1E12.4)
END SUBROUTINE fd_print_problem_setup

SUBROUTINE fd_alloc_work_arrays(create_or_destroy)

USE parameters,             ONLY : alloc_create,alloc_destroy,out_unit
USE shared_data,            ONLY : nij,u,v,p,t,pp,uo,vo,to,uoo,voo,too,&
                                   ap,an,as,ae,aw,su,sv,apu,apv,apt,f1,f2,ft1,ft2,dpx,dpy,&
                                   dux,duy,dvx,dvy,dtx,dty,apt,den,deno,denoo,celbeta,celcp,celcpo,celcpoo,celkappa,lamvisc
USE real_parameters,        ONLY : zero

IMPLICIT NONE
INTEGER,INTENT(IN)      :: create_or_destroy
INTEGER                 :: ierror

IF(create_or_destroy == alloc_create)THEN
  ALLOCATE(u(nij),v(nij),t(nij),p(nij),pp(nij),to(nij),uo(nij),vo(nij),voo(nij),too(nij),uoo(nij),&
           ap(nij),an(nij),as(nij),ae(nij),aw(nij),su(nij),sv(nij),apu(nij),apv(nij),f1(nij),f2(nij),&
           ft1(nij),ft2(nij),dpx(nij),dpy(nij),dux(nij),duy(nij),dvx(nij),dvy(nij),dtx(nij),dty(nij),den(nij),deno(nij),&
           denoo(nij),celbeta(nij),celcp(nij),celcpo(nij),celcpoo(nij),celkappa(nij),lamvisc(nij),STAT=ierror)
  IF(ierror /= 0)WRITE(out_unit,*)'Not enough memory to allocate working arrays.'
  u = zero
  v = zero
  t = zero
  p = zero
  pp = zero
  uo = zero;vo = zero;to = zero
  uoo = zero;voo = zero;too = zero
  ap = zero;an = zero;as = zero;ae = zero;aw = zero
  su = zero;sv = zero
  apu = zero;apv = zero
  f1 = zero;f2 = zero
  ft1 = zero;ft2 = zero
  dpx = zero;dpy = zero
  dux = zero;duy = zero
  dvx = zero;dvy = zero
  dtx = zero;dty = zero
  den = zero;deno = zero;denoo = zero
  celbeta = zero;celcp = zero;celcpo = zero;celcpoo = zero;celkappa = zero;lamvisc = zero
ELSEIF(create_or_destroy == alloc_destroy)THEN
  DEALLOCATE(u,v,t,p,pp,to,uo,vo,voo,too,uoo,&
           ap,an,as,ae,aw,su,sv,apu,apv,f1,f2,ft1,ft2,dpx,dpy,dux,duy,dvx,dvy,dtx,dty,apt,den,deno,denoo,celbeta,lamvisc,&
           celcp,celcpo,celcpoo)
ENDIF

END SUBROUTINE fd_alloc_work_arrays

SUBROUTINE fd_copy_time_arrays(do_extrap)

USE precision,              ONLY : r_single
USE shared_data,            ONLY : v,t,u,vo,uo,to,voo,uoo,too,nij,iu,iv,ip,ien,lcal,&
                                   nim,njm,li
USE real_parameters,        ONLY : two

IMPLICIT NONE
LOGICAL,INTENT(IN) :: do_extrap
REAL(KIND = r_single)   :: uvt
INTEGER                 :: i,ij

IF(do_extrap)THEN
  DO i=2,nim
    DO ij=li(i)+2,li(i)+njm
      uvt = two*u(ij) - uo(ij)
      uoo(ij) = uo(ij)
      uo(ij) = u(ij)
      u(ij) = uvt

      uvt = two*v(ij) - vo(ij)
      voo(ij) = vo(ij)
      vo(ij) = v(ij)
      v(ij) = uvt

      uvt = two*t(ij) - to(ij)
      too(ij) = to(ij)
      to(ij) = t(ij)
      t(ij) = uvt
    ENDDO
  ENDDO
ELSE
  !--Array copy
  TOO=TO
  UOO=UO
  VOO=VO
  TO=T
  UO=U
  VO=V
ENDIF

END SUBROUTINE fd_copy_time_arrays

SUBROUTINE fd_print_var(phi,str)

USE precision,      ONLY : r_single
USE parameters,     ONLY : out_unit
USE shared_data,    ONLY : ni,nj,li
IMPLICIT NONE

CHARACTER(LEN = 6),INTENT(IN) :: str
REAL(KIND = r_single),DIMENSION(:),INTENT(IN) :: phi
INTEGER :: is,ie,nl,i,j,l

WRITE(2,20) str

nl=(ni-1)/12+1

DO l=1,nl
  is=(l-1)*12+1
  ie=MIN(ni,l*12)
  WRITE(out_unit,21) (i,i=is,ie)
  WRITE(out_unit,22)

  DO j=nj,1,-1
    WRITE(out_unit,23) j,(phi(li(i)+j),i=is,ie)
  END DO
END DO

20 FORMAT(2X,26('*-'),5X,A6,5X,26('-*'))
21 FORMAT(3X,'I = ',I3,11I10)
22 FORMAT(2X,'J')
23 FORMAT(1X,I3,1P12E10.2)

END SUBROUTINE fd_print_var


SUBROUTINE fd_alloc_objprop_arrays(create_or_destroy)

USE parameters,       ONLY : alloc_create,alloc_destroy,out_unit
USE shared_data,      ONLY : nsphere,objtp,objqp,fd_urf,densitp,objcentx,objcenty,&
                             objradius,objcentu,objcentv,objcentom,nsurfpoints,&
                             nobjcells,mcellpercv,nnusseltpoints,objcentmi,objvol,&
                             objcento,objcentxo,objcentyo,objcentvo,objcentuo,objcentomo,&
                             betap,up,vp,omp,forcedmotion,stationary,isotherm,objbradius,cpp,kappap,&
                             Fpq,Fpw,zobjcentx,zobjcenty,zobjcentxo,zobjcentyo,Fpqt,Fpwt
USE real_parameters,  ONLY : zero

IMPLICIT NONE

INTEGER,INTENT(IN)      :: create_or_destroy
INTEGER                 :: ierror

IF(create_or_destroy == alloc_create)THEN
  ALLOCATE(objtp(nsphere),densitp(nsphere),objcentx(nsphere),objcenty(nsphere),&
          objradius(nsphere),objcentu(nsphere),objcentv(nsphere),objcentom(nsphere),&
          nsurfpoints(nsphere),nobjcells(nsphere),mcellpercv(nsphere),nnusseltpoints(nsphere),&
          objcentmi(nsphere),objvol(nsphere),objcento(4,nsphere),&
          objcentxo(nsphere),objcentyo(nsphere),objcentvo(nsphere),objcentuo(nsphere),&
          objcentomo(nsphere),betap(nsphere),objbradius(nsphere),cpp(nsphere),kappap(nsphere),&
          Fpq(2,nsphere),Fpw(2,nsphere),zobjcentx(nsphere),zobjcenty(nsphere),&
          zobjcentxo(nsphere),zobjcentyo(nsphere),Fpqt(2,nsphere),Fpwt(2,nsphere),STAT=ierror)
  IF(ierror /= 0)WRITE(out_unit,*)'Not enough memory to allocate onject arrays.'
  IF(forcedmotion .OR. stationary)THEN
    ALLOCATE(up(nsphere),vp(nsphere),omp(nsphere))
    up = zero;vp = zero;omp = zero
  ENDIF
  IF(.NOT. isotherm)THEN
    ALLOCATE(objqp(nsphere))
    objqp = zero
  ENDIF
  objtp = zero;objcentx = zero;objcenty = zero
  objradius = zero; objcentu = zero
  densitp = zero;objcentv = zero; objcentom = zero; objcentmi = zero
  nsurfpoints = 0;nobjcells = 0;mcellpercv = zero;nnusseltpoints = 0;objvol = zero
  objcento = zero;objcentxo = zero ;objcentyo = zero
  objcentvo = zero;objcentuo = zero;objcentomo = zero;betap = zero
  objbradius = zero;cpp = zero;kappap = zero
  Fpq = zero; Fpw = zero; Fpqt = zero; Fpwt = zero 
  zobjcentx = 0; zobjcenty = 0; zobjcentxo = 0; zobjcentyo = 0
ELSEIF(create_or_destroy == alloc_destroy)THEN
  DEALLOCATE(objtp,densitp,objcentx,objcenty,objradius,objcentu,objcentv,objcentom,nsurfpoints,nobjcells,&
             mcellpercv,objcentmi,objvol,objcento,objcentxo,objcentyo,objcentvo,objcentuo,objcentomo,betap,&
             objbradius,cpp,kappap,Fpq,Fpw,zobjcentx,zobjcenty,zobjcentxo,zobjcentyo)
  IF(forcedmotion .OR. stationary)DEALLOCATE(up,vp,omp)
  IF(.NOT. isotherm)DEALLOCATE(objqp)
ENDIF

END SUBROUTINE fd_alloc_objprop_arrays

SUBROUTINE fd_alloc_nusselt_arrays(create_or_destroy,max_point)

USE parameters,       ONLY : alloc_create,alloc_destroy,out_unit
USE shared_data,      ONLY : nusseltpointx,nusseltpointy,nsphere,nusseltpoint_cvx,nusseltpoint_cvy,&
                             nusseltds,nusseltnx,nusseltny,nusseltcentx,nusseltcenty,&
                             nusseltinterpx,nusseltinterpy,nusseltpoint_interp,localnusselt
USE real_parameters,  ONLY : zero

IMPLICIT NONE

INTEGER,INTENT(IN)      :: create_or_destroy,max_point
INTEGER                 :: ierror

IF(create_or_destroy == alloc_create)THEN
  ALLOCATE(nusseltpointx(max_point,nsphere),nusseltpointy(max_point,nsphere),&
           nusseltpoint_cvy(3,max_point,nsphere),nusseltpoint_cvx(3,max_point,nsphere),&
           nusseltnx(max_point,nsphere),nusseltny(max_point,nsphere),nusseltds(max_point,nsphere),&
           nusseltcentx(max_point,nsphere),nusseltcenty(max_point,nsphere),&
           nusseltinterpx(2,max_point,nsphere),nusseltinterpy(2,max_point,nsphere),&
           nusseltpoint_interp(4,2,max_point,nsphere),localnusselt(max_point,nsphere),STAT=ierror)
  IF(ierror /= 0)WRITE(out_unit,*)'Not enough memory to allocate onject arrays.'
  nusseltpointx = zero
  nusseltpointy = zero
  nusseltpoint_cvy = 0
  nusseltpoint_cvx = 0
  nusseltds = zero
  nusseltnx = zero
  nusseltny = zero
  nusseltcentx = zero
  nusseltcenty = zero
  nusseltpoint_interp = zero
  nusseltinterpy = zero
  nusseltinterpx = zero
  localnusselt = zero
ELSEIF(create_or_destroy == alloc_destroy)THEN
  DEALLOCATE(nusseltpointx,nusseltpointy,nusseltpoint_cvy,nusseltpoint_cvx,&
           nusseltnx,nusseltny,nusseltds,nusseltcentx,nusseltcenty,&
           nusseltinterpx,nusseltinterpy,nusseltpoint_interp,localnusselt)
ENDIF

END SUBROUTINE fd_alloc_nusselt_arrays

SUBROUTINE fd_alloc_objgeom_arrays(create_or_destroy,max_point)

USE parameters,       ONLY : alloc_create,alloc_destroy,out_unit
USE shared_data,      ONLY : surfpointx,surfpointy,&
                             surfpointxo,surfpointyo,&
                             nsurfpoints,nsphere,&
                             surfpointxinit,surfpointyinit,forcedmotion,&
                             zsurfpointx,zsurfpointy,&
                             zsurfpointxo,zsurfpointyo

USE real_parameters,  ONLY : zero

IMPLICIT NONE

INTEGER,INTENT(IN)      :: create_or_destroy,max_point
INTEGER                 :: ierror

IF(create_or_destroy == alloc_create)THEN
    ALLOCATE(surfpointx(max_point,nsphere),surfpointy(max_point,nsphere),&
             zsurfpointx(max_point,nsphere),zsurfpointy(max_point,nsphere),&
             surfpointxo(max_point,nsphere),surfpointyo(max_point,nsphere),&
             zsurfpointxo(max_point,nsphere),zsurfpointyo(max_point,nsphere),&
             STAT=ierror)
    IF(ierror /= 0)WRITE(out_unit,*)'Not enough memory to allocate onject arrays.'
    surfpointx   = zero
    surfpointy   = zero
    zsurfpointx  = 0
    zsurfpointy  = 0
    surfpointxo  = zero
    surfpointyo  = zero
    zsurfpointxo = 0
    zsurfpointyo = 0
ELSEIF(create_or_destroy == alloc_destroy)THEN
    DEALLOCATE(surfpointx,surfpointy,zsurfpointx,zsurfpointy,&
               surfpointxo,surfpointyo,zsurfpointxo,zsurfpointyo)
ENDIF

END SUBROUTINE fd_alloc_objgeom_arrays

SUBROUTINE fd_alloc_surfforce_arrays(create_or_destroy,max_point)

USE parameters,       ONLY : alloc_create,alloc_destroy,out_unit
USE shared_data,      ONLY : surfds,surfnx,surfny,nsphere,surfpoint_cvx,surfpoint_cvy,&
                             surfcentrex,surfcentrey,surfpoint_interp,surfinterpx,surfinterpy
USE real_parameters,  ONLY : zero

IMPLICIT NONE

INTEGER,INTENT(IN)      :: create_or_destroy,max_point
INTEGER                 :: ierror

IF(create_or_destroy == alloc_create)THEN
  ALLOCATE(surfds(max_point,nsphere),surfnx(max_point,nsphere),surfny(max_point,nsphere),surfcentrex(max_point,nsphere),&
           surfcentrey(max_point,nsphere),surfpoint_cvx(3,max_point,nsphere),surfpoint_cvy(3,max_point,nsphere),&
           surfpoint_interp(4,2,max_point,nsphere),surfinterpx(2,max_point,nsphere),surfinterpy(2,max_point,nsphere),&
           STAT=ierror)
  IF(ierror /= 0)WRITE(out_unit,*)'Not enough memory to allocate surface arrays.'
  surfds = zero
  surfnx = zero
  surfny = zero
  surfcentrex = zero
  surfcentrey = zero
  surfinterpx = zero
  surfinterpy = zero
  surfpoint_cvx = 0
  surfpoint_cvy = 0
  surfpoint_interp = 0
ELSEIF(create_or_destroy == alloc_destroy)THEN
  DEALLOCATE(surfds,surfnx,surfny,surfpoint_cvx,surfpoint_cvy)
ENDIF

END SUBROUTINE fd_alloc_surfforce_arrays

SUBROUTINE fd_alloc_sources(create_or_destroy)

USE parameters,             ONLY : alloc_create,alloc_destroy,out_unit
USE shared_data,            ONLY : fdsu,fdsv,nij,fdsub,fdsvb,fdst,&
                                   fdsuc,fdsvc,fdstc,&
                                   fdfcu,fdfcv,rigidForce_contrib
USE real_parameters,        ONLY : zero

IMPLICIT NONE
INTEGER,INTENT(IN)      :: create_or_destroy
INTEGER                 :: ierror

IF(create_or_destroy == alloc_create)THEN
  ALLOCATE(fdsu(nij),fdsv(nij),fdsub(nij),fdsvb(nij),fdsuc(nij),fdsvc(nij),&
            fdst(nij),fdstc(nij),fdfcu(nij),fdfcv(nij),rigidForce_contrib(nij),STAT=ierror)
  IF(ierror /= 0)WRITE(out_unit,*)'Not enough memory to allocate working arrays.'
  fdsu = zero
  fdsv = zero
  fdsub= zero
  fdsvb= zero
  fdsuc=zero
  fdsvc=zero
  fdst =zero
  fdstc =zero
  fdfcu = zero
  fdfcv = zero
  rigidForce_contrib = 0
ELSEIF(create_or_destroy == alloc_destroy)THEN
  DEALLOCATE(fdsu,fdsv,fdsub,fdsvb,fdst,fdsuc,fdsvc,fdfcu,fdfcv,rigidForce_contrib)
ENDIF

END SUBROUTINE fd_alloc_sources

SUBROUTINE fd_alloc_objcell_vars(create_or_destroy)

USE parameters,             ONLY : alloc_create,alloc_destroy,out_unit
USE shared_data,            ONLY : nobjcells,obju,objv,objfx,objfy,objru,objrv,&
                                   objq,objrt,objt,nsphere,objapu,objapv,objapt
USE real_parameters,        ONLY : zero

IMPLICIT NONE
INTEGER,INTENT(IN)      :: create_or_destroy
INTEGER                 :: ierror,maxpoint

maxpoint = MAXVAL(nobjcells(:),1)
IF(create_or_destroy == alloc_create)THEN
  ALLOCATE(obju(maxpoint,nsphere),objv(maxpoint,nsphere),objfx(maxpoint,nsphere),objfy(maxpoint,nsphere),&
           objru(maxpoint,nsphere),objrv(maxpoint,nsphere),&
           objrt(maxpoint,nsphere),objt(maxpoint,nsphere),objq(maxpoint,nsphere),STAT=ierror)
  IF(ierror /= 0)WRITE(out_unit,*)'Not enough memory to allocate working arrays.'
  obju = zero
  objv = zero
  objru = zero
  objrv = zero
  objfx = zero
  objfy = zero
  objrt = zero
  objt  = zero
  objq  = zero
ELSEIF(create_or_destroy == alloc_destroy)THEN
  DEALLOCATE(obju,objv,objru,objrv,objfx,objfy,objrt,objt,objq)
ENDIF

END SUBROUTINE fd_alloc_objcell_vars

SUBROUTINE fd_write_restart

USE shared_data
USE parameters,     ONLY : eres_unit

IMPLICIT NONE

INTEGER       :: ij,nn,ik,i,j

  WRITE(eres_unit) xPeriodic,yPeriodic

  WRITE(eres_unit) itim,time,ni,nj,nim,njm,nij,dt

  WRITE(eres_unit)((x(i),j=1,nj),i=1,ni),((y(j),j=1,nj),i=1,ni),&
         ((xc(i),j=1,nj),i=1,ni),((yc(j),j=1,nj),i=1,ni),&
         (f1(ij),ij=1,nij),(f2(ij),ij=1,nij),&
         (ft1(ij),ij=1,nij),(ft2(ij),ij=1,nij),&
         (u(ij),ij=1,nij),&
         (v(ij),ij=1,nij),(p(ij),ij=1,nij),(t(ij),ij=1,nij),&
         (uo(ij),ij=1,nij),(vo(ij),ij=1,nij),(to(ij),ij=1,nij),&
         (uoo(ij),ij=1,nij),(voo(ij),ij=1,nij),(too(ij),ij=1,nij),&
         (dpx(ij),ij=1,nij),(dpy(ij),ij=1,nij),(dux(ij),ij=1,nij),&
         (duy(ij),ij=1,nij),(dvx(ij),ij=1,nij),(dvy(ij),ij=1,nij),&
         (dtx(ij),ij=1,nij),(dty(ij),ij=1,nij),(den(ij),ij=1,nij),(deno(ij),ij=1,nij),&
         (denoo(ij),ij=1,nij),(celbeta(ij),ij=1,nij),(celcp(ij),ij=1,nij),&
         (celcpo(ij),ij=1,nij),(celcpoo(ij),ij=1,nij),&
         (celkappa(ij),ij=1,nij),(fdsu(ij),ij=1,nij),(fdsv(ij),ij=1,nij),(fdsub(ij),ij=1,nij),&
         (fdsvb(ij),ij=1,nij),(fdsuc(ij),ij=1,nij),(fdsvc(ij),ij=1,nij),&
         (fdst(ij),ij=1,nij),(fdstc(ij),ij=1,nij),(lamvisc(ij),ij=1,nij),&
         (fdfcu(ij),ij=1,nij),(fdfcv(ij),ij=1,nij)

  WRITE(eres_unit)temp_visc

  WRITE(eres_unit)putobj,read_fd_geom,stationary,forcedmotion,movingmesh,&
                  calcsurfforce,calclocalnusselt,isotherm

  WRITE(eres_unit)nsphere

  WRITE(eres_unit)(objtp(ij),ij=1,nsphere),(densitp(ij),ij=1,nsphere),(objcentx(ij),ij=1,nsphere),&
                  (objcenty(ij),ij=1,nsphere),(objradius(ij),ij=1,nsphere),(objbradius(ij),ij=1,nsphere),&
                  (objcentu(ij),ij=1,nsphere),&
                  (objcentv(ij),ij=1,nsphere),(objcentom(ij),ij=1,nsphere),(nsurfpoints(ij),ij=1,nsphere),&
                  (nobjcells(ij),ij=1,nsphere),(mcellpercv(ij),ij=1,nsphere),(nnusseltpoints(ij),ij=1,nsphere),&
                  (objcentmi(ij),ij=1,nsphere),(objvol(ij),ij=1,nsphere),&
                  (objcentxo(ij),ij=1,nsphere),(objcentyo(ij),ij=1,nsphere),(objcentvo(ij),ij=1,nsphere),&
                  (objcentuo(ij),ij=1,nsphere),(objcentomo(ij),ij=1,nsphere),(betap(ij),ij=1,nsphere),&
                  (cpp(ij),ij=1,nsphere),(kappap(ij),ij=1,nsphere),(Fpq(1,ij),ij=1,nsphere),&
                  (Fpq(2,ij),ij=1,nsphere),(zobjcentx(ij),ij=1,nsphere),(zobjcenty(ij),ij=1,nsphere),&
                  (zobjcentxo(ij),ij=1,nsphere),(zobjcentyo(ij),ij=1,nsphere)
  DO nn=1,nsphere
    WRITE(eres_unit)(objcento(ij,nn),ij=1,4)
  ENDDO

  IF(forcedmotion)THEN
    DO nn = 1,nsphere
      WRITE(eres_unit)(surfpointx(ij,nn),ij=1,nsurfpoints(nn)),&
                      (surfpointy(ij,nn),ij=1,nsurfpoints(nn)),&
                      (surfpointxinit(ij,nn),ij=1,nsurfpoints(nn)),&
                      (surfpointyinit(ij,nn),ij=1,nsurfpoints(nn)),&
                      (zsurfpointx(ij,nn),ij=1,nsurfpoints(nn)),&
                      (zsurfpointy(ij,nn),ij=1,nsurfpoints(nn))
    ENDDO
  ELSE
    DO nn = 1,nsphere
      WRITE(eres_unit)(surfpointx(ij,nn),ij=1,nsurfpoints(nn)),&
                      (surfpointy(ij,nn),ij=1,nsurfpoints(nn)),&
                      (zsurfpointx(ij,nn),ij=1,nsurfpoints(nn)),&
                      (zsurfpointy(ij,nn),ij=1,nsurfpoints(nn))
    ENDDO
  ENDIF

  IF(calclocalnusselt)THEN
    DO nn=1,nsphere
      WRITE(eres_unit)(nusseltpointx(ij,nn),ij=1,nnusseltpoints(nn)),&
                      (nusseltpointy(ij,nn),ij=1,nnusseltpoints(nn)),&
                      (nusseltnx(ij,nn),ij=1,nnusseltpoints(nn)),&
                      (nusseltny(ij,nn),ij=1,nnusseltpoints(nn)),&
                      (nusseltds(ij,nn),ij=1,nnusseltpoints(nn)),&
                      (nusseltcentx(ij,nn),ij=1,nnusseltpoints(nn)),&
                      (nusseltcenty(ij,nn),ij=1,nnusseltpoints(nn)),&
                      (localnusselt(ij,nn),ij=1,nnusseltpoints(nn))
    ENDDO
    DO nn=1,nsphere
      DO ij=1,nnusseltpoints(nn)
        WRITE(eres_unit)(nusseltpoint_cvy(i,ij,nn),i=1,3),&
                        (nusseltpoint_cvx(i,ij,nn),i=1,3),&
                        (nusseltinterpx(i,ij,nn),i=1,2),&
                        (nusseltinterpy(i,ij,nn),i=1,2)
        DO ik=1,2
           WRITE(eres_unit)(nusseltpoint_interp(i,ik,ij,nn),i=1,4)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF(calcsurfforce)THEN
    DO nn = 1,nsphere
      WRITE(eres_unit)(surfds(ij,nn),ij=1,nsurfpoints(nn)),&
                      (surfnx(ij,nn),ij=1,nsurfpoints(nn)),&
                      (surfny(ij,nn),ij=1,nsurfpoints(nn)),&
                      (surfcentrex(ij,nn),ij=1,nsurfpoints(nn)),&
                      (surfcentrey(ij,nn),ij=1,nsurfpoints(nn))
      DO ij=1,nsurfpoints(nn)
        WRITE(eres_unit)(surfpoint_cvx(i,ij,nn),i=1,3),&
                        (surfpoint_cvy(i,ij,nn),i=1,3)
        DO ik = 1,2
          WRITE(eres_unit)(surfpoint_interp(i,ik,ij,nn),i=1,4),surfinterpx(ik,ij,nn),surfinterpy(ik,ij,nn)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  DO nn=1,nsphere
    WRITE(eres_unit)(objcellx(ij,nn),ij=1,nobjcells(nn)),&
                    (objcelly(ij,nn),ij=1,nobjcells(nn)),&
                    (objcellvol(ij,nn),ij=1,nobjcells(nn)),&
                    (objpoint_cvx(ij,nn),ij=1,nobjcells(nn)),&
                    (objpoint_cvy(ij,nn),ij=1,nobjcells(nn)),&
                    (zobjcellx(ij,nn),ij=1,nobjcells(nn)),&
                    (zobjcelly(ij,nn),ij=1,nobjcells(nn))

    DO ij=1,nobjcells(nn)
      WRITE(eres_unit)(objcellvertx(i,ij,nn),i=1,4),&
                      (objcellverty(i,ij,nn),i=1,4),&
                      (objpoint_interpx(i,ij,nn),i=1,2),&
                      (objpoint_interpy(i,ij,nn),i=1,2),&
                      (zobjcellvertx(i,ij,nn),i=1,4),&
                      (zobjcellverty(i,ij,nn),i=1,4)
    ENDDO

  ENDDO

  IF(forcedmotion)THEN
    WRITE(eres_unit)(objcentxinit(nn),nn=1,nsphere),(objcentyinit(nn),nn=1,nsphere)
    DO nn=1,nsphere
      WRITE(eres_unit)(objcellxinit(ij,nn),ij=1,nobjcells(nn)),&
                      (objcellyinit(ij,nn),ij=1,nobjcells(nn))
      DO ij=1,nobjcells(nn)
        WRITE(eres_unit)(objcellvertxinit(i,ij,nn),i=1,4),&
                        (objcellvertyinit(i,ij,nn),i=1,4)
      ENDDO

    ENDDO

  ENDIF

END SUBROUTINE fd_write_restart


SUBROUTINE fd_update_visc
! This is adopted to model a power-law fluid - themperature dependance is
! commented out
USE precision,        ONLY : r_single
USE shared_data,      ONLY : lamvisc,nim,njm,visc,tref,li,viscgamma,t,duct,movingmesh,ni,nj,th,tc,dux,duy,dvx,dvy,viscN
USE real_parameters,  ONLY : one,two,real_1e2,real_1em6,zero

IMPLICIT NONE
REAL(KIND = r_single) :: nm1,gamdot
INTEGER       :: i,j,ij

nm1 = viscN - one
IF(nm1 == zero)THEN
  lamvisc = visc
ELSE

  DO i = 2,nim
    DO j = 2,njm
      ij = li(i) + j
      gamdot = sqrt(two*dux(ij)**2+two*dvy(ij)**2+(duy(ij)+dvx(ij))**2)
      IF(gamdot > zero)THEN
        gamdot = visc*SQRT(two*dux(ij)**2+two*dvy(ij)**2+(duy(ij)+dvx(ij))**2)**nm1
        lamvisc(ij) = MIN(real_1e2,MAX(gamdot,real_1em6))
      ELSE
        lamvisc(ij) = visc
      ENDIF
    ENDDO
  ENDDO

  !--NORTH BOUNDARY
  DO i=2,nim
    ij=li(i)+nj
    gamdot = sqrt(two*dux(ij)**2+two*dvy(ij)**2+(duy(ij)+dvx(ij))**2)
    IF(gamdot > zero)THEN
      gamdot = visc*SQRT(two*dux(ij)**2+two*dvy(ij)**2+(duy(ij)+dvx(ij))**2)**nm1
      lamvisc(ij) = MIN(real_1e2,MAX(gamdot,real_1em6))
    ELSE
      lamvisc(ij) = visc
    ENDIF
  ENDDO

  !--SOUTH BOUNDARY 
  DO i=2,nim
    ij=li(i)+1
    gamdot = sqrt(two*dux(ij)**2+two*dvy(ij)**2+(duy(ij)+dvx(ij))**2)
    IF(gamdot > zero)THEN
      gamdot = visc*SQRT(two*dux(ij)**2+two*dvy(ij)**2+(duy(ij)+dvx(ij))**2)**nm1
      lamvisc(ij) = MIN(real_1e2,MAX(gamdot,real_1em6))
    ELSE
      lamvisc(ij) = visc
    ENDIF
  ENDDO

  !--WEST BOUNDARY
  DO j=2,njm
    ij=li(1)+j
    gamdot = sqrt(two*dux(ij)**2+two*dvy(ij)**2+(duy(ij)+dvx(ij))**2)
    IF(gamdot > zero)THEN
      gamdot = visc*SQRT(two*dux(ij)**2+two*dvy(ij)**2+(duy(ij)+dvx(ij))**2)**nm1
      lamvisc(ij) = MIN(real_1e2,MAX(gamdot,real_1em6))
    ELSE
      lamvisc(ij) = visc
    ENDIF  
  ENDDO

  !--EAST BOUDARY
  DO j=2,njm
    ij=li(ni)+j
    gamdot = sqrt(two*dux(ij)**2+two*dvy(ij)**2+(duy(ij)+dvx(ij))**2)
    IF(gamdot > zero)THEN
      gamdot = visc*SQRT(two*dux(ij)**2+two*dvy(ij)**2+(duy(ij)+dvx(ij))**2)**nm1
      lamvisc(ij) = MIN(real_1e2,MAX(gamdot,real_1em6))
    ELSE
      lamvisc(ij) = visc
    ENDIF   
  ENDDO

ENDIF
!DO i = 2,nim
!  DO j = 2,njm
!    ij = li(i) + j
!    lamvisc(ij) = -26.99 + 0.09*t(ij) !/(one + viscgamma*(t(ij) - tref))
!  ENDDO
!ENDDO
!
!!--NORTH BOUNDARY
!DO i=2,nim
!  ij=li(i)+nj
!  lamvisc(ij) = -26.99 + 0.09*t(ij) !visc !/(one + viscgamma*(t(ij) - tref))
!ENDDO
!
!!--SOUTH BOUNDARY (ISOTHERMAL WALL, NON-ZERO DIFFUSIVE FLUX)
!DO i=2,nim
!  ij=li(i)+1
!  lamvisc(ij) = -26.99 + 0.09*t(ij) !/(one + viscgamma*(t(ij) - tref))
!ENDDO
!
!!--WEST BOUNDARY
!DO j=2,njm
!  ij=li(1)+j
!  lamvisc(ij) = -26.99 + 0.09*t(ij) !/(one + viscgamma*(t(ij) - tref))
!ENDDO
!
!
!!--EAST BOUNDARY
!IF(duct)THEN !--Outlet (zero grdient but t is not available use previous node)
!  DO j=2,njm
!    ij=li(ni)+j
!    lamvisc(ij) = -26.99 + 0.09*t(ij-nj) !/(one + viscgamma*(t(ij-nj) - tref))
!  END DO
!ELSE
!  IF(movingmesh)THEN !--Outlet
!    DO j=2,njm
!      ij=li(ni)+j
!      lamvisc(ij) = -26.99 + 0.09*t(ij-nj) !/(one + viscgamma*(t(ij-nj) - tref))
!    END DO
!  ELSE
!    !--EAST BOUNDARY
!    DO j=2,njm
!      ij=li(ni)+j
!      lamvisc(ij) = -26.99 + 0.09*t(ij) !/(one + viscgamma*(t(ij) - tref))
!    ENDDO
!  ENDIF
!ENDIF


END SUBROUTINE fd_update_visc

SUBROUTINE fd_calc_timeCoefs

USE shared_data,      ONLY : dt,dto,ct1,ct2,ct3,time,gamt
USE real_parameters,  ONLY : one,two,zero

IMPLICIT NONE

!--If dt should be changed change it here {Begin}
 !...
!--{End}

time = time + dt

IF(gamt == 1)THEN
  ct1 = one/dt
  ct2 = one/dt !--This is actually -ct2 (since it goes to sources)
  ct3 = zero
ELSEIF(gamt == 2)THEN
  ct1 = (two*dt + dto)/(dt*(dt+dto))
  ct2 = (dt + dto)/(dt*dto) !--This is -ct2 (see above)
  ct3 = -dt/(dto*(dt+dto))  !--This is -ct3
ENDIF

dto = dt

END SUBROUTINE fd_calc_timeCoefs

END MODULE modfd_problem_setup
