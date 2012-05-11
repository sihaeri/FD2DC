MODULE modfd_tecwrite

PRIVATE
PUBLIC :: fd_tecwrite_eul,fd_tecwrite_fil,fd_tecwrite_sph_v,fd_tecwrite_sph_s
CONTAINS
!=====================================================
SUBROUTINE fd_tecwrite_eul(tec_unit,extra_Var)

USE precision,      ONLY : r_single
USE shared_data,    ONLY : title,nim,njm,xc,u,v,p,yc,li,t,x,y,den,celbeta,celprr,nij
USE parameters,     ONLY : max_len_tecline

IMPLICIT NONE

INTEGER,INTENT(IN)  :: tec_unit
REAL(KIND = r_single),OPTIONAL      :: extra_var(nij)
CHARACTER(LEN = max_len_tecline)    :: variableline
INTEGER             :: i,j,ij

WRITE(tec_unit,'(A)')'TITLE="'//TRIM(title)//'"'

IF(PRESENT(extra_var))THEN

  variableline = 'VARIABLES ="x", "y", "u", "v", "p" ,"t", "den", "beta", "prandtle", "extra"'
  WRITE(tec_unit,*)TRIM(variableline)
  WRITE(tec_unit,'(A,I5,A,I5)')'ZONE DATAPACKING=BLOCK, VARLOCATION = ([3,4,5,6,7,8,9,10]=CELLCENTERED), I=',nim, ',J=',njm

ELSE

  variableline = 'VARIABLES ="x", "y", "u", "v", "p" ,"t", "den", "beta", "prandtle"'
  WRITE(tec_unit,*)TRIM(variableline)
  WRITE(tec_unit,'(A,I5,A,I5)')'ZONE DATAPACKING=BLOCK, VARLOCATION = ([3,4,5,6,7,8,9]=CELLCENTERED), I=',nim, ',J=',njm

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
    WRITE(tec_unit,*)celprr(ij)
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

SUBROUTINE fd_tecwrite_sph_v(tec_unit,nsphere,nobjcells,vertx,verty)

USE precision,      ONLY : r_single
USE shared_data,    ONLY : title

IMPLICIT NONE

INTEGER,INTENT(IN)  :: tec_unit,nsphere
INTEGER,INTENT(IN),DIMENSION(:)   :: nobjcells
REAL(KIND = r_single),INTENT(IN),DIMENSION(:,:,:) :: vertx,verty
INTEGER :: j,n,totp,snp
    
WRITE(tec_unit,'(A)')'TITLE="'//TRIM(title)//'"'
WRITE(tec_unit,9012)
WRITE(tec_unit,*)'ZONE T=Test,DATAPACKING=POINT,NODES=',4*SUM(nobjcells(1:nsphere)),&
             ',ELEMENTS=',SUM(nobjcells(1:nsphere)),',ZONETYPE=FEQUADRILATERAL,'

DO n = 1,nsphere
      
  DO j=1,nobjcells(n)
    WRITE(tec_unit,9015)vertx(1,j,n),verty(1,j,n)
  ENDDO

  DO j=1,nobjcells(n)
    WRITE(tec_unit,9015)vertx(2,j,n),verty(2,j,n)
  ENDDO
 
  DO j=1,nobjcells(n)
    WRITE(tec_unit,9015)vertx(3,j,n),verty(3,j,n)
  ENDDO               

  DO j=1,nobjcells(n)
    WRITE(tec_unit,9015)vertx(4,j,n),verty(4,j,n)
  ENDDO
   
ENDDO

totp = 0
DO n = 1,nsphere
  snp = nobjcells(n)
  DO j=1,snp
    WRITE(tec_unit,9013)j+totp,j+snp+totp,j+2*snp+totp,j+3*snp+totp 
  ENDDO
  totp = totp + 4*snp
ENDDO


 9012   FORMAT('variables="x","y"')
 9013   FORMAT(4I8)
 9015   FORMAT(2E14.6)

END SUBROUTINE fd_tecwrite_sph_v

SUBROUTINE fd_tecwrite_sph_s(tec_unit,nsphere,nsurfpoints,posx,posy)

USE precision,      ONLY : r_single
USE shared_data,    ONLY : title

IMPLICIT NONE

INTEGER,INTENT(IN)  :: tec_unit,nsphere
INTEGER,INTENT(IN),DIMENSION(:)   :: nsurfpoints
REAL(KIND = r_single),INTENT(IN),DIMENSION(:,:) :: posx,posy
INTEGER :: j,n

WRITE(tec_unit,'(A)')'TITLE="'//TRIM(title)//'"'
WRITE(tec_unit,9012)

DO n=1,nsphere
  WRITE(tec_unit,9013)nsurfpoints(n)
  DO j=1,nsurfpoints(n)
    WRITE(tec_unit,9015) posx(j,n),posy(j,n)
  ENDDO
ENDDO

 9012   FORMAT('variables="x","y"')
 9013   FORMAT('zone ,i=',I8, ',DATAPACKING=POINT')
 9015   FORMAT(2E14.6)

END SUBROUTINE fd_tecwrite_sph_s

END MODULE modfd_tecwrite