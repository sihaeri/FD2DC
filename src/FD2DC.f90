PROGRAM FD2DC

USE real_parameters,     ONLY : zero,half
USE precision,           ONLY : r_single
USE modfd_openclose_ops, ONLY : fd_open,fd_open_sres
USE modfd_problem_setup, ONLY : fd_problem_setup,fd_problem_restart,fd_copy_time_arrays,&
                                fd_print_var,fd_write_restart,fd_update_visc
USE shared_data,         ONLY : lread,itim,itst,louts,lcal,iu,iv,ip,ien,u,p,v,t,ltime,time,minit,&
                                maxit,u,v,p,t,resor,nprt, filenum,problem_name,problem_len,&
                                imon,jmon,slarge,sormax,ni,nj,loute,lwrite,nim,njm,nij,&
                                x,y,xc,yc,li,f1,f2,vo,uo,to,putobj,dt,duct,filnprt,spherenum,isotherm,objtp,&
                                nfil,nfilpoints,filpointx,filpointy,filru,filrv,calcsurfforce,&
                                localnusselt,calclocalnusselt,nusseltcentx,nusseltcenty,nnusseltpoints,&
                                nsphere,nfil,filfirstposy,sphnprt,surfpointx,surfpointy,objcellx,objcelly,&
                                nobjcells,nsurfpoints,objcellvertx,objcellverty,objcentu,objcentv,objcentx,&
                                objcenty,stationary,th,movingmesh,dxmeanmoved,dxmean,dxmeanmovedtot,&
                                lamvisc,temp_visc,tc,ulid,objradius,naverage_steps,calclocalnusselt_ave,objcentom,&
                                voo,uoo,ft1,ft2,den,deno,celbeta,celcp,celkappa,too,calcwalnusselt
USE parameters,          ONLY : out_unit,eres_unit,plt_unit,force_predict,force_correct,do_collect_stat,end_collect_stat
USE modfd_set_bc,        ONLY : fd_bctime,fd_bcout
USE modfd_calc_pre,      ONLY : fd_calc_pre
USE modfd_calc_temp,     ONLY : fd_calc_temp
USE modfd_calc_mom,      ONLY : fd_calc_mom
USE mod_create_filenum,  ONLY : create_filenum
USE modfd_calc_integrals,ONLY : fd_calc_integrals,fd_calc_surf_force,fd_calc_surf_nusselt,fd_calc_wall_nusselt,&
                                fd_calc_lwall_nusselt,fd_calc_surf_nusselt_ave
USE modfd_tecwrite,      ONLY : fd_tecwrite_eul,fd_tecwrite_fil,fd_tecwrite_sph_v,fd_tecwrite_sph_s
USE modfd_create_geom,   ONLY : fd_calc_sources,fd_calc_mi,fd_calc_ori,fd_calc_pos,fd_calc_physprops,fd_copy_oldvel,&
                                fd_move_mesh,fd_calculate_stats,fd_calc_part_collision
USE modfil_create_geom,  ONLY : fil_sol_tension,fil_sol_positions,fil_calc_sources
use omp_lib
IMPLICIT NONE
LOGICAL         :: ismoved
INTEGER         :: itims,itime,ijmon,iter,n,i,j,ij,m,maxfl,l,nn,npoint
REAL(KIND = r_single) :: source,fd_resor,cd,cl,avenusselt,vb_resor,dummy
REAL(KIND = r_single),ALLOCATABLE :: wallocnusselt(:)
DOUBLE PRECISION      :: delt1,delt2,T1,T2

fd_resor = zero 
vb_resor = zero

CALL fd_open
CALL fd_problem_setup

IF(lread)THEN
  CALL fd_open_sres
  CALL fd_problem_restart
ENDIF

IF(putobj .AND. .NOT. ltime)THEN
  WRITE(out_unit,*) 'Cannot perform steady FD calculation.'
  STOP
ENDIF

itims=itim+1
itime=itim+itst
n = 0
m = 0
l = 0

IF(ltime .AND. lwrite)THEN
  IF(putobj)THEN
    IF(nfil > 0 .AND. nsphere > 0)THEN
      maxfl = MAX(itime/nprt,itime/filnprt,itime/sphnprt)
    ELSEIF(nfil > 0 .AND. nsphere == 0)THEN
      maxfl = MAX(itime/nprt,itime/filnprt)
    ELSEIF(nfil == 0 .AND. nsphere > 0)THEN
      maxfl = MAX(itime/nprt,itime/sphnprt)
    ENDIF
    ALLOCATE(filenum(maxfl),spherenum(nsphere))
    CALL create_filenum(maxfl,filenum)
    CALL create_filenum(nsphere,spherenum)
  ELSE
    maxfl = itime/nprt
    ALLOCATE(filenum(maxfl))
    CALL create_filenum(maxfl,filenum)
  ENDIF
ENDIF
npoint = 0 !--for time steps used in time averaging 
DO itim = itims,itime

  time = time + dt
  IF(putobj .AND. nsphere > 0)CALL fd_copy_oldvel
  IF(ltime)CALL fd_copy_time_arrays(.FALSE.)

  WRITE(out_unit,*) '     TIME = ',time
  WRITE(out_unit,*) '     *****************************'
  WRITE(out_unit,*) '  '

  IF(louts.AND.(itim == itims)) THEN
    IF(lcal(iu)) CALL fd_print_var(u,'U VEL.')
    IF(lcal(iv)) CALL fd_print_var(v,'V VEL.')
    IF(lcal(ip)) CALL fd_print_var(p,'PRESS.')
    IF(lcal(ien)) CALL fd_print_var(t,'TEMPER')
  ENDIF

  !--DEFINE MONITORING LOCATION (NODE WITH I=IMON, J=JMON)

  ijmon=li(imon)+jmon
  WRITE(out_unit,600) imon,jmon

  !--SET BOUNDARY CONDITIONS FOR THE NEW TIME STEP
  IF(ltime) CALL fd_bctime
  IF(putobj)THEN
    !IF(nsphere>0)CALL fd_calc_sources(force_correct,fd_resor,0)
    IF(nfil>0 )CALL fil_calc_sources(delt1,vb_resor,0)
  ENDIF
  !--OUTER ITERATIONS (SIMPLE RELAXATIONS)

  DO iter=1,maxit
!    T1 = omp_get_wtime()
    IF(lcal(iu)) CALL fd_calc_mom
    IF(lcal(iu).AND.duct) CALL fd_bcout
    IF(lcal(iu).AND.movingmesh) CALL fd_bcout
    IF(lcal(ip)) CALL fd_calc_pre
    IF(lcal(ien)) CALL fd_calc_temp
    IF(temp_visc)CALL fd_update_visc
    IF(putobj)THEN
      IF(nsphere>0)CALL fd_calc_sources(force_correct,fd_resor,iter)
      IF(nfil>0)CALL fil_calc_sources(delt1,vb_resor,iter)
    ENDIF
!    T2 = omp_get_wtime()
!    OPEN(UNIT = 12345,FILE='TIME.dat',STATUS = 'UNKNOWN',ACCESS = 'APPEND')
!    WRITE(12345,*) T2-T1, delt1
!    CLOSE(12345)
    !--CHECK CONVERGENCE OF OUTER ITERATIONS
    WRITE(out_unit,606) iter,resor(iu),resor(iv),resor(ip),&
                resor(ien),vb_resor,u(ijmon),v(ijmon),p(ijmon),t(ijmon)
    source=MAX(resor(iu),resor(iv),resor(ip),resor(ien))
    IF(source>slarge)THEN
        PRINT *,'  *** TERMINATED - OUTER ITERATIONS DIVERGING ***'
        STOP
    ENDIF
    IF(source < sormax .AND. iter > minit) EXIT
  END DO
  
  IF(nfil > 0)THEN
    !CALL fil_sol_tension
    !CALL fil_sol_positions
  ENDIF

  IF(nsphere > 0)THEN
   !CALL fd_calc_surf_force(densref,0.1,1.0)
   IF(.NOT. stationary)THEN
     DO nn = 1,nsphere
       OPEN(1234,FILE='vel_'//spherenum(nn)//'.dat',STATUS = 'UNKNOWN',ACCESS = 'APPEND')
       OPEN(1235,FILE='pos_'//spherenum(nn)//'.dat',STATUS = 'UNKNOWN',ACCESS = 'APPEND')
       IF(.NOT. isotherm)THEN
         OPEN(1236,FILE='temp_'//spherenum(nn)//'.dat',STATUS = 'UNKNOWN',ACCESS = 'APPEND')
         WRITE(1236,*)objtp(nn)
         CLOSE(1236)
       ENDIF
       WRITE(1234,*)objcentom(nn),objcentu(nn),objcentv(nn) !
       WRITE(1235,*)objcentx(nn)-dxmeanmovedtot,objcenty(nn)
       CLOSE(1234);CLOSE(1235)
     ENDDO
     CALL fd_calc_mi
     CALL fd_calc_ori
     CALL fd_calc_pos(itim-itims+1)
   ENDIF
   CALL fd_calc_physprops
   IF(nsphere > 0 .AND. movingmesh)CALL fd_move_mesh(ismoved)
   IF(stationary)THEN
     CALL fd_calc_sources(force_correct,fd_resor,iter)
   ELSE
     CALL fd_calc_sources(force_predict,fd_resor,iter)
   ENDIF
  ENDIF
  !--CONVERGED: IF UNSTEADY FLOW, PRINT AND SAVE NPRTth SOLUTION
  IF((.NOT.ltime).OR.(ltime.AND.MOD(itim,nprt)==0)) THEN
    IF(loute) THEN
      IF(lcal(iu)) CALL fd_print_var(U,'U VEL.')
      IF(lcal(iv)) CALL fd_print_var(V,'V VEL.')
      IF(lcal(ip)) CALL fd_print_var(P,'PRESS.')
      IF(lcal(ien)) CALL fd_print_var(T,'TEMPER')
    ENDIF
    
    CALL fd_calc_integrals

    IF(lwrite .AND. ltime) THEN
      n = n +1
      OPEN(UNIT = plt_unit,FILE=problem_name(1:problem_len)//filenum(n)//'.plt',STATUS='NEW')
      
      IF(temp_visc)THEN
        CALL fd_tecwrite_eul(plt_unit,lamvisc)
      ELSE
        CALL fd_tecwrite_eul(plt_unit)
      ENDIF
      CLOSE(plt_unit)
      IF(calcwalnusselt .AND. .NOT. duct)THEN
        OPEN(UNIT = plt_unit,FILE='wallnusselt'//filenum(n)//'.dat',STATUS='NEW')
        ALLOCATE(wallocnusselt(nj))
        CALL fd_calc_lwall_nusselt(wallocnusselt,th,tc)
        DO nn=2,njm
          WRITE(plt_unit,'(2(E12.5,1X))') yc(nn),wallocnusselt(nn)
        ENDDO
        CLOSE(plt_unit)

        DEALLOCATE(wallocnusselt)
      ENDIF

    ELSEIF(lwrite .AND. .NOT. ltime)THEN
     
      OPEN(UNIT = plt_unit,FILE=problem_name(1:problem_len)//'.plt',STATUS='NEW')
      
      IF(temp_visc)THEN
        CALL fd_tecwrite_eul(plt_unit,lamvisc)
      ELSE
        CALL fd_tecwrite_eul(plt_unit)
      ENDIF
       
      CLOSE(plt_unit)
    ENDIF
  ENDIF
  IF(calcsurfforce .AND. stationary)THEN
    CALL fd_calc_surf_force(1,cd,cl)
    OPEN(UNIT = 1234,FILE='drag.dat',STATUS = 'UNKNOWN',ACCESS = 'APPEND')
    !######Only for vertical oscillating sphere########
    !WRITE(1234,'(4(E12.5,1X))') filfirstposy(1),cd,cl,vb_resor
    !######Only for vertical oscillating sphere########
    WRITE(1234,'(2(E12.5,1X))') cd,cl
    CLOSE(1234)
  ENDIF 
  IF(nsphere > 0)THEN
    IF((ltime.AND.MOD(itim,sphnprt)==0)) THEN
      l = l +1
      OPEN(UNIT = plt_unit,FILE=problem_name(1:problem_len)//filenum(l)//'_sph_v.plt',STATUS='NEW')
      CALL fd_tecwrite_sph_v(plt_unit,nsphere,nobjcells,objcellvertx,objcellverty)
      CLOSE(plt_unit)
      OPEN(UNIT = plt_unit,FILE=problem_name(1:problem_len)//filenum(l)//'_sph_s.plt',STATUS='NEW')
      CALL fd_tecwrite_sph_s(plt_unit,nsphere,nsurfpoints,surfpointx,surfpointy)
      CLOSE(plt_unit)
    ENDIF
  ENDIF
  IF(nfil > 0)THEN
    IF((ltime.AND.MOD(itim,filnprt)==0)) THEN
      m = m +1
      OPEN(UNIT = plt_unit,FILE=problem_name(1:problem_len)//filenum(m)//'_fil.plt',STATUS='NEW')
      CALL fd_tecwrite_fil(plt_unit,nfil,nfilpoints,filpointx,filpointy,filru,filrv)
      CLOSE(plt_unit)
    ENDIF
  ENDIF
  !IF( time*ulid/objradius(1) > 90.0)CALL fd_calculate_stats(1,do_collect_stat) 
  !CALL fd_calculate_stats(1,do_collect_stat)
  IF(calclocalnusselt_ave)THEN
    IF( (itime - itims + 1) <= naverage_steps .OR. (itime - itim + 1) <= naverage_steps)THEN
      DO nn = 1,nsphere
        CALL fd_calc_surf_nusselt_ave(nn,npoint)
      ENDDO
      !WRITE(*,*)npoint, (itime - itims + 1), (itime - itim + 1), naverage_steps
    ENDIF
  ENDIF
ENDDO
!CALL fd_calculate_stats(1,end_collect_stat)
IF(calclocalnusselt)THEN

  OPEN(UNIT = 1234,FILE='cylnusselt.dat',STATUS = 'UNKNOWN')

  DO n = 1,nsphere
    IF(calclocalnusselt_ave)THEN
      CALL fd_calc_surf_nusselt(n,localnusselt(:,n),avenusselt,npoint)
    ELSE
      CALL fd_calc_surf_nusselt(n,localnusselt(:,n),avenusselt)
    ENDIF
    DO i = 1,nnusseltpoints(n)-1
      WRITE(1234,'(3(E12.5,1X))') nusseltcentx(i,n),nusseltcenty(i,n),localnusselt(i,1)
    ENDDO
    WRITE(1234,'(E12.5,1X)')avenusselt
  ENDDO
  CLOSE(1234)
  !OPEN(UNIT = 1234,FILE='walnusselt.dat',STATUS = 'UNKNOWN')
  !ALLOCATE(wallocnusselt(nj))
  !CALL fd_calc_wall_nusselt(wallocnusselt,th,1)
  !DO n=2,njm
  !  WRITE(1234,'(2(E12.5,1X))') yc(n),wallocnusselt(n)
  !ENDDO
  !CLOSE(1234)
ENDIF
IF(putobj)THEN
  CALL fd_write_restart
ELSE
  WRITE(eres_unit) itim,time,ni,nj,nim,njm,nij,&
         ((x(i),j=1,nj),i=1,ni),((y(j),j=1,nj),i=1,ni),&
         ((xc(i),j=1,nj),i=1,ni),((yc(j),j=1,nj),i=1,ni),&
         (f1(ij),ij=1,nij),(f2(ij),ij=1,nij),(ft1(ij),ij=1,nij),(ft2(ij),ij=1,nij),&
         (u(ij),ij=1,nij),(v(ij),ij=1,nij),(p(ij),ij=1,nij),(t(ij),ij=1,nij),&
         (uo(ij),ij=1,nij),(vo(ij),ij=1,nij),(to(ij),ij=1,nij),&
         (uoo(ij),ij=1,nij),(voo(ij),ij=1,nij),(too(ij),ij=1,nij),&
         (den(ij),ij=1,nij),(deno(ij),ij=1,nij),(celbeta(ij),ij=1,nij),&
         (celcp(ij),ij=1,nij),(celkappa(ij),ij=1,nij)
ENDIF



600 FORMAT(1X,'ITER.',3X,&
           'I---------ABSOLUTE RESIDUAL SOURCE SUMS--------I',3X,&
           'I----FIELD VALUES AT MONITORING LOCATION (',I3,',',I3,&
           ')----I',/,2X,'NO.',9X,'U',11X,'V',9X,'MASS',10X,'T',&
           16X,'U',11X,'V',11X,'P',11X,'T',/)

606 FORMAT(1X,I4,2X,1P5E12.4,5X,1P4E12.4)

END PROGRAM FD2DC