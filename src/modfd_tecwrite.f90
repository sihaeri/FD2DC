MODULE modfd_tecwrite

PRIVATE
PUBLIC :: fd_tecwrite_eul,fd_tecwrite_fil,fd_tecwrite_sph_v,fd_tecwrite_sph_s
CONTAINS
!=====================================================
SUBROUTINE fd_tecwrite_eul(tec_unit,extra_Var)

USE precision,      ONLY : r_single
USE shared_data,    ONLY : title,nim,njm,xc,u,v,p,yc,li,t,x,y,den,celbeta,celcp,celkappa,nij
USE parameters,     ONLY : max_len_tecline

IMPLICIT NONE

INTEGER,INTENT(IN)  :: tec_unit
REAL(KIND = r_single),OPTIONAL      :: extra_var(nij)
CHARACTER(LEN = max_len_tecline)    :: variableline
INTEGER             :: i,j,ij

WRITE(tec_unit,'(A)')'TITLE="'//TRIM(title)//'"'

IF(PRESENT(extra_var))THEN

  variableline = 'VARIABLES ="x", "y", "u", "v", "p" ,"t", "den", "beta", "cp", "kappa", "extra"'
  WRITE(tec_unit,*)TRIM(variableline)
  WRITE(tec_unit,'(A,I5,A,I5)')'ZONE DATAPACKING=BLOCK, VARLOCATION = ([3,4,5,6,7,8,9,10,11]=CELLCENTERED), I=',nim, ',J=',njm

ELSE

  variableline = 'VARIABLES ="x", "y", "u", "v", "p" ,"t", "den", "beta", "cp" , "kappa" '
  WRITE(tec_unit,*)TRIM(variableline)
  WRITE(tec_unit,'(A,I5,A,I5)')'ZONE DATAPACKING=BLOCK, VARLOCATION = ([3,4,5,6,7,8,9,10]=CELLCENTERED), I=',nim, ',J=',njm

ENDIF
DO j = 1,njm
  DO i = 1,nim
    WRITE(tec_unit,*)x(i)
  ENDDO
ENDDO

DO j = 1,njm
  DO i = 1,nim
    WRITE(tec_unit,*)y(j)
  ENDDO
ENDDO

DO j = 2,njm
  DO i = 2,nim
    ij = li(i)+j
    WRITE(tec_unit,*)u(ij)
  ENDDO
ENDDO

DO j = 2,njm
  DO i = 2,nim
    ij = li(i)+j
    WRITE(tec_unit,*)v(ij)
  ENDDO
ENDDO

DO j = 2,njm
  DO i = 2,nim
    ij = li(i)+j
    WRITE(tec_unit,*)p(ij)
  ENDDO
ENDDO

DO j = 2,njm
  DO i = 2,nim
    ij = li(i)+j
    WRITE(tec_unit,*)t(ij)
  ENDDO
ENDDO

DO j = 2,njm
  DO i = 2,nim
    ij = li(i)+j
    WRITE(tec_unit,*)den(ij)
  ENDDO
ENDDO

DO j = 2,njm
  DO i = 2,nim
    ij = li(i)+j
    WRITE(tec_unit,*)celbeta(ij)
  ENDDO
ENDDO

DO j = 2,njm
  DO i = 2,nim
    ij = li(i)+j
    WRITE(tec_unit,*)celcp(ij)
  ENDDO
ENDDO

DO j = 2,njm
  DO i = 2,nim
    ij = li(i)+j
    WRITE(tec_unit,*)celkappa(ij)
  ENDDO
ENDDO

IF(PRESENT(extra_var))THEN
  DO j = 2,njm
    DO i = 2,nim
      ij = li(i)+j
      WRITE(tec_unit,*)extra_var(ij)
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE fd_tecwrite_eul

SUBROUTINE fd_tecwrite_fil(tec_unit,nfil,npoints,posx,posy,u,v)

USE precision,      ONLY : r_single
USE shared_data,    ONLY : title

IMPLICIT NONE

INTEGER,INTENT(IN)  :: tec_unit,nfil
INTEGER,INTENT(IN),DIMENSION(:)   :: npoints
REAL(KIND = r_single),INTENT(IN),DIMENSION(:,:) :: posx,posy,u,v
INTEGER :: j,n

WRITE(tec_unit,'(A)')'TITLE="'//TRIM(title)//'"'
WRITE(tec_unit,9012)

DO n=1,nfil
  WRITE(tec_unit,9013)npoints(n)
  DO j=1,npoints(n)
    WRITE(tec_unit,9015) posx(j,n),posy(j,n),u(j,n),v(j,n)
  ENDDO
ENDDO

CLOSE(tec_unit)

 9012   FORMAT('variables="x","y","vx","vy"')
 9013   FORMAT('zone ,i=',I5, ',DATAPACKING=POINT')
 9015   FORMAT(4E14.6)

END SUBROUTINE fd_tecwrite_fil

SUBROUTINE fd_tecwrite_sph_v(tec_unit,nsphere,nobjcells,vertx,verty,zcellx,zcelly,zvertx,zverty)

USE precision,      ONLY : r_single
USE shared_data,    ONLY : title,LDomainx,LDomainy

IMPLICIT NONE

INTEGER,INTENT(IN)  :: tec_unit,nsphere
INTEGER,INTENT(IN),DIMENSION(:)   :: nobjcells
REAL(KIND = r_single),INTENT(IN),DIMENSION(:,:,:) :: vertx,verty
INTEGER,INTENT(IN),DIMENSION(:,:)   :: zcellx,zcelly
INTEGER,INTENT(IN),DIMENSION(:,:,:) :: zvertx,zverty

!--Locals
REAL(KIND = r_single) :: xx,yy
INTEGER :: j,n,nspzone(4,nsphere),czonex,czoney,nz,jj
INTEGER,ALLOCATABLE,SAVE :: spzone(:,:,:)
LOGICAL,SAVE             :: entered = .FALSE.

IF(.NOT.entered)THEN
  ALLOCATE(spzone(4,MAXVAL(nobjcells(1:nsphere)),nsphere))
  entered = .TRUE.
ENDIF

spzone = 0
nspzone = 0

DO n = 1,nsphere
  czonex = zcellx(1,n)
  czoney = zcelly(1,n)
  DO j = 1,nobjcells(n)
    IF(zcellx(j,n) /= czonex .AND. zcelly(j,n) /= czoney)THEN !--Crossed edge
      nspzone(4,n) = nspzone(4,n) + 1
      spzone(4,nspzone(4,n),n) = j
    ELSEIF(zcellx(j,n) /= czonex .AND. zcelly(j,n) == czoney)THEN !--Crossed e-w boundary but no edge
      nspzone(3,n) = nspzone(3,n) + 1
      spzone(3,nspzone(3,n),n) = j
    ELSEIF(zcellx(j,n) == czonex .AND. zcelly(j,n) /= czoney)THEN !--Crossed n-s boundary but no edge 
      nspzone(2,n) = nspzone(2,n) + 1
      spzone(2,nspzone(2,n),n) = j
    ELSE
      nspzone(1,n) = nspzone(1,n) + 1
      spzone(1,nspzone(1,n),n) = j
    ENDIF         
  ENDDO
ENDDO
 
WRITE(tec_unit,'(A)')'TITLE="'//TRIM(title)//'"'
WRITE(tec_unit,9012)

DO n = 1,nsphere
  DO nz = 1,4
    IF(nspzone(nz,n) /= 0)THEN  
      WRITE(tec_unit,*)'ZONE T=Test,DATAPACKING=POINT,NODES=',4*nspzone(nz,n),&
                   ',ELEMENTS=',nspzone(nz,n),',ZONETYPE=FEQUADRILATERAL,'
      DO j=1,nspzone(nz,n)
        jj = spzone(nz,j,n)
        xx = vertx(1,jj,n)+(zvertx(1,jj,n) - zcellx(jj,n))*LDomainx
        yy = verty(1,jj,n)+(zverty(1,jj,n) - zcelly(jj,n))*LDomainy
        WRITE(tec_unit,9015)xx,yy
      ENDDO

      DO j=1,nspzone(nz,n)
        jj = spzone(nz,j,n)
        xx = vertx(2,jj,n)+(zvertx(2,jj,n) - zcellx(jj,n))*LDomainx
        yy = verty(2,jj,n)+(zverty(2,jj,n) - zcelly(jj,n))*LDomainy
        WRITE(tec_unit,9015)xx,yy
      ENDDO

      DO j=1,nspzone(nz,n) 
        jj = spzone(nz,j,n)
        xx = vertx(3,jj,n)+(zvertx(3,jj,n) - zcellx(jj,n))*LDomainx
        yy = verty(3,jj,n)+(zverty(3,jj,n) - zcelly(jj,n))*LDomainy
        WRITE(tec_unit,9015)xx,yy
      ENDDO               

      DO j=1,nspzone(nz,n) 
        jj = spzone(nz,j,n)
        xx = vertx(4,jj,n)+(zvertx(4,jj,n) - zcellx(jj,n))*LDomainx
        yy = verty(4,jj,n)+(zverty(4,jj,n) - zcelly(jj,n))*LDomainy
        WRITE(tec_unit,9015)xx,yy
      ENDDO
   
      DO j=1,nspzone(nz,n)
        WRITE(tec_unit,9013)j,j+nspzone(nz,n),j+2*nspzone(nz,n),j+3*nspzone(nz,n)
      ENDDO

    ENDIF
  ENDDO
ENDDO


 9012   FORMAT('variables="x","y"')
 9013   FORMAT(4I8)
 9015   FORMAT(2E14.6)

END SUBROUTINE fd_tecwrite_sph_v

SUBROUTINE fd_tecwrite_sph_s(tec_unit,nsphere,nsurfpoints,posx,posy,zonex,zoney)

USE precision,      ONLY : r_single
USE shared_data,    ONLY : title

IMPLICIT NONE

INTEGER,INTENT(IN)  :: tec_unit,nsphere
INTEGER,INTENT(IN),DIMENSION(:)   :: nsurfpoints
INTEGER,INTENT(IN),DIMENSION(:,:) :: zonex,zoney
REAL(KIND = r_single),INTENT(IN),DIMENSION(:,:) :: posx,posy
INTEGER :: j,n,czonex,czoney,nspzone(4,nsphere),nz
INTEGER,ALLOCATABLE,SAVE :: spzone(:,:,:)
LOGICAL,SAVE             :: entered = .FALSE.

IF(.NOT.entered)THEN
  ALLOCATE(spzone(4,MAXVAL(nsurfpoints(1:nsphere)),nsphere))
  entered = .TRUE.
ENDIF

spzone = 0
nspzone = 0

DO n = 1,nsphere
  czonex = zonex(1,n)
  czoney = zoney(1,n)
  DO j = 1,nsurfpoints(n)
    IF(zonex(j,n) /= czonex .AND. zoney(j,n) /= czoney)THEN !--Crossed edge
      nspzone(4,n) = nspzone(4,n) + 1
      spzone(4,nspzone(4,n),n) = j
    ELSEIF(zonex(j,n) /= czonex .AND. zoney(j,n) == czoney)THEN !--Crossed e-w boundary but no edge
      nspzone(3,n) = nspzone(3,n) + 1
      spzone(3,nspzone(3,n),n) = j
    ELSEIF(zonex(j,n) == czonex .AND. zoney(j,n) /= czoney)THEN !--Crossed n-s boundary but no edge 
      nspzone(2,n) = nspzone(2,n) + 1
      spzone(2,nspzone(2,n),n) = j
    ELSE
      nspzone(1,n) = nspzone(1,n) + 1
      spzone(1,nspzone(1,n),n) = j
    ENDIF         
  ENDDO
ENDDO

WRITE(tec_unit,'(A)')'TITLE="'//TRIM(title)//'"'
WRITE(tec_unit,9012)

DO n=1,nsphere
  DO nz = 1,4
    IF(nspzone(nz,n) > 0)THEN
      WRITE(tec_unit,9013)nspzone(nz,n)
      DO j=1,nspzone(nz,n)
        WRITE(tec_unit,9015) posx(spzone(nz,j,n),n),posy(spzone(nz,j,n),n)
      ENDDO
    ENDIF
  ENDDO
ENDDO

 9012   FORMAT('variables="x","y"')
 9013   FORMAT('zone ,i=',I8, ',DATAPACKING=POINT')
 9015   FORMAT(2E14.6)

END SUBROUTINE fd_tecwrite_sph_s

END MODULE modfd_tecwrite