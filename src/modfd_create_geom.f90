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

MODULE modfd_create_geom

PRIVATE
PUBLIC :: fd_create_geom,fd_calc_sources,fd_calc_delta,find_single_point,fd_calc_mi,fd_calc_ori,fd_copy_oldvel,fd_calc_pos,&
          fd_calc_physprops,fd_move_mesh,fd_calc_quality,fd_calculate_stats,fd_calc_part_collision,fd_init_temp,&
          fd_track_single_point, fd_calc_forcedmotion, fd_update_forcedvel, fd_create_interpmol,fd_find_rigidforce_contrib,fd_update_fieldvel

CONTAINS
!===========================================
SUBROUTINE fd_create_geom

USE shared_data,        ONLY : objcentx,objcenty,objradius,objbradius,nsurfpoints,&
                               surfpointx,surfpointy,njm,nim,nobjcells,&
                               objcellx,objcelly,objcelldx,objcelldy,objcellvol,&
                               x,y,r,objpoint_cvx,objpoint_cvy,densitp,mcellpercv,&
                               nsphere,surfds,surfnx,surfny,surfcentrex,surfcentrey,&
                               calcsurfforce,xc,yc,surfpoint_cvx,surfpoint_cvy,&
                               surfpoint_interp,surfinterpx,surfinterpy,li,&
                               nusseltpointx,nusseltpointy,calclocalnusselt,nnusseltpoints,&
                               nusseltcentx,nusseltcenty,nusseltds,nusseltnx,nusseltny,nusseltpoint_interp,&
                               nusseltinterpx,nusseltinterpy,nusseltpoint_cvx,nusseltpoint_cvy,objvol,&
                               surfpointxinit,surfpointyinit,zsurfpointx,zsurfpointy,zobjcellx,zobjcelly,&
                               zobjcellvertx,zobjcellverty,xPeriodic,yPeriodic,LDomainx,LDomainy,&
                               objcellvertx,objcellverty,objpoint_interpx,objpoint_interpy,&
                               problem_name,problem_len,read_fd_geom,dxmean,hypre_A,Hypre_b,Hypre_x,mpi_comm ,nij,&
                               dxmeanmovedtot,objcentxinit,objcentyinit,objcellxinit,dxmeanmoved,dxmin,&
                               objcellyinit,objcellvertxinit,objcellvertyinit,forcedmotion,objcell_bndflag
USE precision,          ONLY : r_single
USE parameters,         ONLY : plt_unit
USE real_parameters,    ONLY : two,pi,half,zero,four,large,three,small,one
USE modfd_tecwrite,     ONLY : fd_tecwrite_sph_s,fd_tecwrite_sph_v

IMPLICIT NONE
INTEGER                 :: i,j,ncell,ii,jj,mini,maxnn,maxnnb,n,maxnobjcell,maxnsurfpoint,ip1,error
INTEGER,ALLOCATABLE     :: nn(:),nnb(:)
REAL(KIND = r_single)   :: dtheta,theta,domain_vol,dx,dy,mindist,unitvectx,unitvecty,&
                           dist,posx,posy,rp,vecx,vecy,costheta,sintheta,ds
REAL(KIND = r_single),ALLOCATABLE :: dummyposarray(:,:,:),dummyobjcellx(:,:),dummyobjcelly(:,:),&
                                     dummydx(:,:),dummydy(:,:)                       


ALLOCATE(nn(nsphere),nnb(nsphere))
dxmeanmoved = zero
dxmeanmovedtot = zero
domain_vol = zero
ncell = 0
dxmin = large 
DO i=2,nim
  dx=x(i)-x(i-1)
  dxmin = MIN(dx,dxmin)
  DO j=2,njm
    ncell = ncell + 1
    dy=y(j)-y(j-1)
    dxmin = min(dy,dxmin)
    rp=half*(r(j)+r(j-1))
    domain_vol = domain_vol + (dx*dy*rp)**half
  ENDDO
ENDDO
            
dxmean = domain_vol/REAL(ncell)
WRITE(*,*) "dxmean = ",dxmean," ,dxmin = ", dxmin
DO n = 1,nsphere
  nn(n) = CEILING(two*objradius(n) / (dxmean/REAL(mcellpercv(n),r_single)))
  nnb(n) = CEILING(two*objbradius(n) / (dxmean/REAL(mcellpercv(n),r_single)))
ENDDO

maxnn = MAXVAL(nn(:),1)
maxnnb = MAXVAL(nnb(:),1)
maxnn = MAX(maxnn,maxnnb)
ALLOCATE(dummyposarray(maxnn+1,6,nsphere))
ALLOCATE(dummyobjcellx(maxnn**2,nsphere),dummyobjcelly(maxnn**2,nsphere),dummydx(maxnn**2,nsphere),&
         dummydy(maxnn**2,nsphere))
dummyposarray = zero

DO n = 1,nsphere
  dtheta = two*pi / REAL(nsurfpoints(n))
  DO i=1,nsurfpoints(n)
  
    theta = REAL(i-1)*dtheta
    costheta = COS(theta)
    sintheta = SIN(theta)
    surfpointx(i,n) = objcentx(n)+(objradius(n)*objbradius(n)/ &
                      SQRT((objradius(n)*costheta)**2+(objbradius(n)*sintheta)**2))*costheta
    surfpointy(i,n) = objcenty(n)+(objradius(n)*objbradius(n)/ &
                      SQRT((objradius(n)*costheta)**2+(objbradius(n)*sintheta)**2))*sintheta

  ENDDO
!  CALL fd_create_ellipse_surf(surfpointx(:,n),surfpointy(:,n),nsurfpoints(n),objbradius(n),objradius(n))
!  surfpointx(:,n) = surfpointx(:,n) + objcentx(n)
!  surfpointy(:,n) = surfpointy(:,n) + objcenty(n)
  IF(calclocalnusselt)CALL fd_create_interpmol_nus(n)
    
  IF(calcsurfforce)CALL fd_create_interpmol(n)
  
  IF(.NOT. read_fd_geom)THEN  
    DO i=1,nnb(n)+1
      dummyPosArray(i,1,n)= (objcentx(n) - objbradius(n)) + two*objbradius(n)/REAL(nnb(n),r_single)*(i-1)
    ENDDO

    DO i=1,nn(n)+1
      dummyPosArray(i,2,n)= (objcenty(n) - objradius(n)) + two*objradius(n)/REAL(nn(n),r_single)*(i-1)
    ENDDO

    DO i=1,nnb(n)
      dummyPosArray(i,3,n) = half*(dummyPosArray(i,1,n)+ dummyPosArray(i+1,1,n))
    ENDDO
      
    DO i=1,nn(n)
      dummyPosArray(i,4,n) = half*(dummyPosArray(i,2,n)+ dummyPosArray(i+1,2,n))
    ENDDO
    
    DO i=1,nnb(n)-1
      dummyPosArray(i,5,n)=dummyPosArray(i+1,3,n)-dummyPosArray(i,3,n)
    ENDDO
    dummyPosArray(i,5,n) = two*objbradius(n) + dummyPosArray(1,3,n) - dummyPosArray(i,3,n)

    DO i=1,nn(n)-1
      dummyPosArray(i,6,n)=dummyPosArray(i+1,4,n)-dummyPosArray(i,4,n)
    ENDDO
    dummyPosArray(i,6,n) = two*objradius(n) + dummyPosArray(1,4,n) - dummyPosArray(i,4,n)
            
    nobjcells(n) = 0

    DO ii = 1,nnb(n)
      DO jj = 1,nn(n)
        mindist = large
        DO i = 1,nsurfpoints(n)
          posx = dummyPosArray(ii,3,n)
          posy = dummyPosArray(jj,4,n)
          dist = SQRT((posx-surfpointx(i,n))**2 + (posy - surfpointy(i,n))**2)
          IF(dist < minDist)THEN
            minDist = dist
            mini   = i
          ENDIF
        ENDDO

        IF(mini < nsurfpoints(n))THEN
          ip1 = mini + 1
        ELSE
          ip1 = 1
        ENDIF
      

        vecx = surfpointx(mini,n) - surfpointx(ip1,n)
        vecy = surfpointy(mini,n) - surfpointy(ip1,n)
        
            
        ds = SQRT(vecx**2 + vecy**2)

        !--Surface element normals (outward)
        unitvectx = -vecy/ds
        unitvecty = vecx/ds
    
        IF(((surfpointx(mini,n)-dummyposarray(ii,3,n))*unitvectx+&
            (surfpointy(mini,n)-dummyPosArray(jj,4,n))*unitvecty)>zero) THEN
          nobjcells(n) = nobjcells(n) + 1
          dummyobjcellx(nobjcells(n),n) = dummyposarray(ii,3,n)
          dummyobjcelly(nobjcells(n),n) = dummyposarray(jj,4,n)
          dummydx(nobjcells(n),n)  = dummyPosArray(ii,5,n)
          dummydy(nobjcells(n),n)  = dummyPosArray(jj,6,n)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
ENDDO

IF(.NOT. read_fd_geom)THEN
  maxnobjcell = MAXVAL(nobjcells(:),1)
  ALLOCATE(objcellx(maxnobjcell,nsphere),objcelly(maxnobjcell,nsphere),zobjcellx(maxnobjcell,nsphere),&
           zobjcelly(maxnobjcell,nsphere),objcellvol(maxnobjcell,nsphere),objpoint_cvx(maxnobjcell,nsphere),&
           objpoint_cvy(maxnobjcell,nsphere),objcellvertx(4,maxnobjcell,nsphere),objcellverty(4,maxnobjcell,nsphere),&
           zobjcellvertx(4,maxnobjcell,nsphere),zobjcellverty(4,maxnobjcell,nsphere),&
           objpoint_interpx(2,maxnobjcell,nsphere),objpoint_interpy(2,maxnobjcell,nsphere),&
           objcell_bndFlag(maxnobjcell,nsphere))
  zobjcellx = 0; zobjcelly = 0
  zobjcellvertx = 0; zobjcellverty = 0
  objcell_bndFlag = 0
  DO n = 1,nsphere
   
    objcellx(1:nobjcells(n),n) = dummyobjcellx(1:nobjcells(n),n)
    objcelly(1:nobjcells(n),n) = dummyobjcelly(1:nobjcells(n),n) 
    
    objcellvertx(1,1:nobjcells(n),n) = objcellx(1:nobjcells(n),n) - half*dummydx(1:nobjcells(n),n)
    objcellverty(1,1:nobjcells(n),n) = objcelly(1:nobjcells(n),n) - half*dummydy(1:nobjcells(n),n)

    objcellvertx(2,1:nobjcells(n),n) = objcellx(1:nobjcells(n),n) + half*dummydx(1:nobjcells(n),n)
    objcellverty(2,1:nobjcells(n),n) = objcelly(1:nobjcells(n),n) - half*dummydy(1:nobjcells(n),n)

    objcellvertx(3,1:nobjcells(n),n) = objcellx(1:nobjcells(n),n) + half*dummydx(1:nobjcells(n),n)
    objcellverty(3,1:nobjcells(n),n) = objcelly(1:nobjcells(n),n) + half*dummydy(1:nobjcells(n),n)

    objcellvertx(4,1:nobjcells(n),n) = objcellx(1:nobjcells(n),n) - half*dummydx(1:nobjcells(n),n)
    objcellverty(4,1:nobjcells(n),n) = objcelly(1:nobjcells(n),n) + half*dummydy(1:nobjcells(n),n)

  ENDDO 
  
  IF(xPeriodic == 1)THEN
    DO n = 1,nsphere
      DO i = 1,nsurfpoints(n)
        IF(surfpointx(i,n) > x(nim))THEN
          surfpointx(i,n) = surfpointx(i,n) - LDomainx
          zsurfpointx(i,n) = zsurfpointx(i,n) + 1
        ELSEIF(surfpointx(i,n) < x(1))THEN
          surfpointx(i,n) = surfpointx(i,n) + LDomainx
          zsurfpointx(i,n) = zsurfpointx(i,n) - 1
        ENDIF
      ENDDO
    ENDDO
    DO n = 1, nsphere
      DO i = 1,nobjcells(n)
        
        IF(objcellx(i,n) > x(nim))THEN
          objcellx(i,n) = objcellx(i,n) - LDomainx !--Went through the east boundary
          zobjcellx(i,n) = zobjcellx(i,n) + 1
        ELSEIF(objcellx(i,n) < x(1))THEN
          objcellx(i,n) = objcellx(i,n) + LDomainx !--Went through the west boundary
          zobjcellx(i,n) = zobjcellx(i,n) - 1
        ENDIF

        DO j = 1,4
          IF(objcellvertx(j,i,n) > x(nim))THEN
            objcellvertx(j,i,n) = objcellvertx(j,i,n) - LDomainx !--Went through the east boundary
            zobjcellvertx(j,i,n) = zobjcellvertx(j,i,n) + 1
          ELSEIF(objcellvertx(j,i,n) < x(1))THEN
            objcellvertx(j,i,n) = objcellvertx(j,i,n) + LDomainx !--Went through the west boundary
            zobjcellvertx(j,i,n) = zobjcellvertx(j,i,n) - 1
          ENDIF          
        ENDDO

      ENDDO
    ENDDO
  ENDIF

  IF(yPeriodic == 1)THEN
    
    DO n = 1,nsphere
      DO i = 1,nsurfpoints(n)
        IF(surfpointy(i,n) > y(njm))THEN
          surfpointy(i,n) = surfpointy(i,n) - LDomainy
          zsurfpointy(i,n) = zsurfpointy(i,n) + 1
        ELSEIF(surfpointy(i,n) < y(1))THEN
          surfpointy(i,n) = surfpointy(i,n) + LDomainy
          zsurfpointy(i,n) = zsurfpointy(i,n) - 1
        ENDIF
      ENDDO
    ENDDO

    DO n = 1, nsphere
      DO i = 1,nobjcells(n)
        
        IF(objcelly(i,n) > y(njm))THEN
          objcelly(i,n) = objcelly(i,n) - LDomainy !--Went through the east boundary
          zobjcelly(i,n) = zobjcelly(i,n) + 1
        ELSEIF(objcelly(i,n) < y(1))THEN
          objcelly(i,n) = objcelly(i,n) + LDomainy !--Went through the west boundary
          zobjcelly(i,n) = zobjcelly(i,n) - 1
        ENDIF

        DO j = 1,4
          IF(objcellverty(j,i,n) > y(njm))THEN
            objcellverty(j,i,n) = objcellverty(j,i,n) - LDomainy !--Went through the north boundary
            zobjcellverty(j,i,n) = zobjcellverty(j,i,n) + 1
          ELSEIF(objcellverty(j,i,n) < y(1))THEN
            objcellverty(j,i,n) = objcellverty(j,i,n) + LDomainy !--Went through the south boundary
            zobjcellverty(j,i,n) = zobjcellverty(j,i,n) - 1
          ENDIF          
        ENDDO

      ENDDO
    ENDDO
  ENDIF

  DO n = 1,nsphere
    DO i=1,nobjcells(n)
      objcellvol(i,n) = dummydx(i,n)*dummydy(i,n)
      objvol(n) = objvol(n) + objcellvol(i,n)
    ENDDO
     WRITE(*,*)'Real Area for Obj: ',n,' is: ',pi*objradius(n)*objbradius(n),', estimated: ', objvol(n),'.'
  ENDDO

ELSE
  CALL fd_geom_reader(error)
  IF(error /= 0)RETURN
ENDIF

IF(forcedmotion)THEN
  maxnobjcell = MAXVAL(nobjcells(:),1)
  maxnsurfpoint = MAXVAL(nsurfpoints(:),1)
  ALLOCATE(objcentxinit(nsphere),objcentyinit(nsphere),objcellxinit(maxnobjcell,nsphere),&
           objcellyinit(maxnobjcell,nsphere),objcellvertxinit(4,maxnobjcell,nsphere),&
           objcellvertyinit(4,maxnobjcell,nsphere),&
           surfpointxinit(maxnsurfpoint,nsphere),surfpointyinit(maxnsurfpoint,nsphere))
  objcentxinit = objcentx
  objcentyinit = objcenty
  objcellxinit = objcellx
  objcellyinit = objcelly
  objcellvertxinit = objcellvertx
  objcellvertyinit = objcellverty
  surfpointxinit = surfpointx
  surfpointyinit = surfpointy
ENDIF

DEALLOCATE(dummyobjcellx,dummyobjcelly,dummydx,dummydy,dummyposarray)

!--Find the cells
CALL find_objcells

!--draw mesh using tecplot
OPEN(UNIT = plt_unit,FILE=problem_name(1:problem_len)//'_init'//'_sph_v.plt',STATUS='NEW')
CALL fd_tecwrite_sph_v(plt_unit,nsphere,nobjcells,objcellvertx,objcellverty,zobjcellx,zobjcelly,zobjcellvertx,zobjcellverty)
CLOSE(plt_unit)
OPEN(UNIT = plt_unit,FILE=problem_name(1:problem_len)//'_init'//'_sph_s.plt',STATUS='NEW')
CALL fd_tecwrite_sph_s(plt_unit,nsphere,nsurfpoints,surfpointx,surfpointy,zsurfpointx,zsurfpointy)
CLOSE(plt_unit)


END SUBROUTINE fd_create_geom

SUBROUTINE fd_calc_unitvector(x1, y1, x2, y2, unitvectx, unitvecty, veclen)

USE precision,      ONLY : r_single
USE real_parameters,ONLY : small
IMPLICIT NONE

REAL(KIND = r_single),INTENT(IN) :: x1,y1,x2,y2
REAL(KIND = r_single),INTENT(OUT):: unitvectx, unitvecty
REAL(KIND = r_single),INTENT(OUT):: veclen

unitvectx = x2 - x1
unitvecty = y2 - y1

veclen = SQRT(unitvectx**2 + unitvecty**2)

IF(ABS(vecLen) < small) THEN
  PRINT *,'Error - Coordinates co-exist - undefined normal!'
  STOP
ENDIF

unitvectx = unitvectx/veclen
unitvecty = unitvecty/veclen

END SUBROUTINE fd_calc_unitvector

SUBROUTINE find_objcells

USE shared_data,    ONLY : objpoint_cvx,objpoint_cvy,x,y,nim,njm,&
                           objcellx,objcelly,nobjcells,nsphere,deltalen,&
                           objpoint_interpx,objpoint_interpy,nobjcells,xc,yc

IMPLICIT NONE

INTEGER            :: i,j,n

DO n = 1,nsphere
  DO j = 1,nobjcells(n)
    DO i = 2,nim
      IF(objcellx(j,n) <= x(i) .AND. objcellx(j,n) > x(i-1))THEN
        objpoint_cvx(j,n) = i
        EXIT
      ENDIF
    ENDDO

    IF(i>nim)THEN
      WRITE(*,*)'Could not find cell,', j,', for particle,', n,', in x-dirextion.'
    ENDIF
     
    DO i = 2,njm
      IF(objcelly(j,n) <= y(i) .AND. objcelly(j,n) > y(i-1))THEN
        objpoint_cvy(j,n) = i
        EXIT
      ENDIF
    ENDDO
    IF(i>njm)THEN
      WRITE(*,*)'Could not find cell,', j,', for particle,', n,', in y-dirextion.'
    ENDIF

  ENDDO


  IF(deltalen == 3)THEN
    DO j = 1,nobjcells(n)
      objpoint_interpx(1,j,n) = objpoint_cvx(j,n) - 1
      objpoint_interpx(2,j,n) = objpoint_cvx(j,n) + 1
      objpoint_interpy(1,j,n) = objpoint_cvy(j,n) - 1
      objpoint_interpy(2,j,n) = objpoint_cvy(j,n) + 1
    ENDDO
  ELSEIF(deltalen == 2 .OR. deltalen == 4 .OR. deltalen == 6)THEN
    DO j = 1,nobjcells(n)
      IF(objcellx(j,n) <= xc(objpoint_cvx(j,n)))THEN
        objpoint_interpx(1,j,n) = objpoint_cvx(j,n) - deltalen/2
        objpoint_interpx(2,j,n) = objpoint_cvx(j,n) + deltalen/2 - 1
      ELSE
        objpoint_interpx(1,j,n) = objpoint_cvx(j,n) - deltalen/2 + 1
        objpoint_interpx(2,j,n) = objpoint_cvx(j,n) + deltalen/2
      ENDIF
      IF(objcelly(j,n) <= yc(objpoint_cvy(j,n)))THEN
        objpoint_interpy(1,j,n) = objpoint_cvy(j,n) - deltalen/2
        objpoint_interpy(2,j,n) = objpoint_cvy(j,n) + deltalen/2 - 1
      ELSE
        objpoint_interpy(1,j,n) = objpoint_cvy(j,n) - deltalen/2 + 1
        objpoint_interpy(2,j,n) = objpoint_cvy(j,n) + deltalen/2
      ENDIF
    ENDDO
  ENDIF
ENDDO

END SUBROUTINE find_objcells

SUBROUTINE fd_track_single_point(bndFlag,xp,yp,ipo,jpo,ip,jp,interpx,interpy,error)

USE precision,      ONLY : r_single
USE shared_data,    ONLY : x,y,deltalen,xc,yc,nim,njm,xPeriodic,yPeriodic,objcell_bndFlag
USE parameters,     ONLY : maxmove

IMPLICIT NONE
INTEGER,INTENT(INOUT)   :: ip,jp,ipo,jpo,error,interpx(2),interpy(2),bndFlag
REAL(KIND = r_single),INTENT(IN)   :: xp,yp
INTEGER                            :: i,nim2,njm2,ii

nim2 = nim - 1
njm2 = njm - 1

error = 0
bndFlag = 0
DO i = ipo - maxmove,ipo + maxmove
  IF(xPeriodic == 1)THEN
    ii = 2+MOD(nim2+MOD(i-2,nim2),nim2)
  ELSE
    ii = i 
    IF(xp < x(2) .OR. xp > x(nim-1))THEN
        ip = ipo
        bndFlag = 1
        EXIT
    ENDIF
    IF(ii <= 1 .OR. ii > nim)CYCLE
  ENDIF
  IF(xp <= x(ii) .AND. xp > x(ii-1))THEN
    ip = ii
    EXIT
  ENDIF
  
ENDDO
      
IF( i == ipo + maxmove + 1)THEN
  error = 1
  RETURN
ENDIF
      
DO i = jpo - maxmove,jpo + maxmove
  IF(yPeriodic == 1)THEN
    ii = 2+MOD(njm2+MOD(i-2,njm2),njm2)
  ELSE
    ii = i
    IF(yp < y(2) .OR. yp > y(njm-1))THEN
        jp = jpo
        bndFlag = bndFlag + 2
        EXIT
    ENDIF
    IF(ii <= 1 .OR. ii > njm)CYCLE
  ENDIF
  IF(yp <= y(ii) .AND. yp > y(ii-1))THEN
    jp = ii
    EXIT
  ENDIF
ENDDO

IF( i == jpo + maxmove + 1)THEN
  error = 2
  RETURN
ENDIF

IF(jp /= jpo .OR. ip /= ipo)THEN
  IF(deltalen == 3)THEN
    interpx(1) = ip - 1
    interpx(2) = ip + 1
    interpy(1) = jp - 1
    interpy(2) = jp + 1
  ELSEIF(deltalen == 2 .OR. deltalen == 4 .OR. deltalen == 6)THEN
    IF(xp <= xc(ip))THEN
      interpx(1) = ip - deltalen/2
      interpx(2) = ip + deltalen/2 - 1
    ELSE
      interpx(1) = ip - deltalen/2 + 1
      interpx(2) = ip + deltalen/2
    ENDIF
    IF(yp <= yc(jp))THEN
      interpy(1) = jp - deltalen/2
      interpy(2) = jp + deltalen/2 - 1
    ELSE
      interpy(1) = jp - deltalen/2 + 1
      interpy(2) = jp + deltalen/2
    ENDIF
  ENDIF
ENDIF

END SUBROUTINE fd_track_single_point

SUBROUTINE find_single_point(xp,yp,ip,jp)

USE precision,      ONLY : r_single
USE shared_data,    ONLY : x,y,nim,njm

IMPLICIT NONE
INTEGER,INTENT(INOUT)   :: ip,jp
REAL(KIND = r_single),INTENT(IN)   :: xp,yp

INTEGER               :: i

!bndFlag = 0

!IF(xp < x(1) .OR. xp > x(nim))THEN
!  ip = -1
!  bndFlag = 1
!  RETURN
!ENDIF

DO i = 2,nim
  IF(xp <= x(i) .AND. xp > x(i-1))THEN
    ip = i
    EXIT
  ENDIF
ENDDO

!IF(yp < y(1) .OR. yp > y(njm))THEN
!  jp = -1
!  bndFlag = bndFlag + 2
!  RETURN
!ENDIF

DO i = 2,njm
  IF(yp <= y(i) .AND. yp > y(i-1))THEN
    jp = i
    EXIT
  ENDIF
ENDDO

END SUBROUTINE find_single_point

FUNCTION fd_calc_delta(x, y , objx, objy, hx, hy)

!--A three point delta function
USE precision,      ONLY : r_single
USE real_parameters,ONLY : one,zero,eight,three,four,seven,twelve,five,two,half,six
USE shared_data,    ONLY : deltalen
IMPLICIT NONE

REAL(KIND = r_single) :: x,y,objx,objy,hx,hy
REAL(KIND = r_single) :: rx,ry
REAL(KIND = r_single) :: fd_calc_delta
REAL(KIND = r_single) :: deltaX, deltaY

rx = (x - objx)/hx
ry = (y - objy)/hy

IF(deltalen == 3)THEN  

  IF (ABS(rx) <= half) THEN
    deltax = one/three*(one + SQRT(-three*rx**2 + one))
    !deltax = three/four-rx**2
  ELSE IF(ABS(rx) > half .AND. ABS(rx) <= (one+half)) THEN
    deltax = one/six*(five - three*ABS(rx) - SQRT(-three*(one-ABS(rx))**2 + one))
    !deltax = nine/eight-three/two*ABS(rx)+rx**2/two
  ELSE IF(ABS(rx) > (one+half)) THEN
    deltax = zero
  END IF
        
  IF (ABS(ry) <= half) THEN
    deltay = one/three*(one + SQRT(-three*ry**2 + one))
    !deltay = three/four-ry**2
  ELSE IF(ABS(ry) > half .AND. ABS(ry) <= (one+half)) THEN
    deltay = one/six*(five - three*ABS(ry) - SQRT(-three*(one-ABS(ry))**2 + one))
    !deltay = nine/eight-three/two*ABS(ry)+ry**2/two
  ELSE IF(ABS(ry) > one+half) THEN
    deltay = zero
  END IF

ELSEIF(deltalen == 2)THEN
  
  IF (ABS(rx) <= one) THEN
    deltax = one - ABS(rx)
  ELSE IF(ABS(rx) > one) THEN
    deltax = zero
  END IF
  
  IF (ABS(ry) <= one) THEN
    deltay = one - ABS(ry)
  ELSE IF(ABS(ry) > one) THEN
    deltay = zero
  END IF

ELSEIF(deltalen == 4)THEN

  IF (ABS(rx) <= one) THEN
    deltax = one/eight*(three - two*ABS(rx) + SQRT(one+four*ABS(rx)-four*rx**2))
  ELSE IF(ABS(rx) > one .AND. ABS(rx) <= two) THEN
    deltax = one/eight*(five - two*ABS(rx) - SQRT(-seven+twelve*ABS(rx)-four*rx**2))
  ELSE IF(ABS(rx) > two) THEN
    deltax = zero
  END IF
        
  IF (ABS(ry) <= one) THEN
    deltay = one/eight*(three - two*ABS(ry) + SQRT(one+four*ABS(ry)-four*ry**2))
  ELSE IF(ABS(ry) > one .AND. ABS(ry) <= two) THEN
    deltay = one/eight*(five - two*ABS(ry) - SQRT(-seven+twelve*ABS(ry)-four*ry**2))
  ELSE IF(ABS(ry) > two) THEN
    deltay = zero
  END IF 

ELSEIF(deltalen == 6)THEN

  deltax = phir6(rx)
  deltay = phir6(ry)

ENDIF

fd_calc_delta = deltaX/hx*deltaY/hy

END FUNCTION fd_calc_delta

RECURSIVE FUNCTION phir6(r) RESULT(resPhir6)

USE precision,        ONLY : r_single
IMPLICIT NONE

REAL(KIND = r_single) :: r,absr,resPhir6

absr = ABS(r)

If (absr <= 1.) THEN
  resPhir6 = 61. / 112. - 11. / 42. * absr - 11. / 56. * absr**2 + 1. / 12. * absr**3 + SQRT(3.) / 336. * &
        (243. + 1584. * absr - 748. * absr**2 - 1560. * absr**3 + 500. * absr**4 + 336. * absr**5 - 112. * absr**6)**(0.5)
ELSEIF (absr > 1. .AND. absr <= 2.) THEN
  resPhir6 = 21. / 16. + 7. / 12. * absr - 7. / 8. * absr**2 + 1. / 6. * absr**3 - 3. / 2. * phir6(absr - 1.)
ELSEIf (absr >= 2. .AND. absr < 3.) THEN
  resPhir6 = 9. / 8. - 23. / 12. * absr + 3. / 4. * absr**2 - 1. / 12. * absr**3 + 1. / 2. * phir6(absr - 2.)
ELSE
  resPhir6 = 0.
ENDIF

END FUNCTION phir6

SUBROUTINE fd_calc_sources(predict_or_correct,resor,iter)

USE parameters,     ONLY : force_predict,force_correct,no_force
USE shared_data,    ONLY : fdsu,fdsv,objfx,objfy,objru,objrv,obju,objv,xPeriodic,yPeriodic,&
                           objpoint_cvx, objpoint_cvy,x,y,objcellx,objcelly,&
                           xc,yc,objcentu,objcentv,objcentom,objcentx,objcenty,&
                           dt,li,objcellvol,nobjcells,u,v,densitp,nij,fd_urf,objvol,&
                           nim,njm,fx,fy,fdsub,fdsvb,nj,objtp,objt,objrt,objq,fdst,t,&
                           nsphere,fdsvc,fdsuc,fdstc,apu,apv,objapu,objapv,objapt,apt,&
                           objpoint_interpx,objpoint_interpy,stationary,isotherm,objqp,cpp,&
                           zobjcenty,zobjcentx,zobjcelly,zobjcellx,LDomainx,LDomainy,rigidforce_contrib,&
                           objcell_bndFlag
USE precision,      ONLY : r_single
USE real_parameters,ONLY : zero,one
IMPLICIT NONE

!--Inputs
INTEGER              :: predict_or_correct,iter
REAL(KIND = r_single):: resor 
!--locals
INTEGER             :: i,j,n,nn,ii,jj,nim2,njm2,ij
REAL(KIND = r_single) :: dx,dy,del,curresor,maxresor,xp,yp

nim2 = nim - 1
njm2 = njm - 1
fdsuc = zero
fdsvc = zero
fdstc = zero

obju = zero
objv = zero
objt = zero

DO nn = 1,nsphere
  DO n = 1,nobjcells(nn)
    IF(objcell_bndFlag(n,nn) == 0)THEN
      DO i = objpoint_interpx(1,n,nn),objpoint_interpx(2,n,nn)
        ii = i + xPeriodic*(MAX(2-i,0)*nim2 + MIN(nim-i,0)*nim2)
        xp = objcellx(n,nn) + xPeriodic*(MAX(2-i,0)*LDomainx + MIN(nim-i,0)*LDomainx)
        dx = x(ii) - x(ii-1)
        DO j = objpoint_interpy(1,n,nn),objpoint_interpy(2,n,nn)
          jj = j + yPeriodic*(MAX(2-j,0)*njm2 + MIN(njm-j,0)*njm2) 
          ij = li(ii)+jj 
          !IF(rigidForce_contrib(ij) == nn)THEN
            yp = objcelly(n,nn) + yPeriodic*(MAX(2-j,0)*LDomainy + MIN(njm-j,0)*LDomainy)
            dy = y(jj) - y(jj-1)
            del = fd_calc_delta(xc(ii),yc(jj),xp,yp,dx, dy)*dx*dy
            obju(n,nn) = obju(n,nn) + u(ij)*del
            objv(n,nn) = objv(n,nn) + v(ij)*del
            objt(n,nn) = objt(n,nn) + t(ij)*del
          !ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDDO
ENDDO

CALL fd_calc_rigidvel

IF(predict_or_correct /= no_force)THEN
IF(isotherm)THEN
  DO nn = 1,nsphere
    DO n = 1,nobjcells(nn)
      objru(n,nn) = objcentu(nn) - objcentom(nn)*(objcelly(n,nn) - objcenty(nn)+(zobjcelly(n,nn) - zobjcenty(nn))*LDomainy)
      objrv(n,nn) = objcentv(nn) + objcentom(nn)*(objcellx(n,nn) - objcentx(nn)+(zobjcellx(n,nn) - zobjcentx(nn))*LDomainx)
      objrt(n,nn) = objtp(nn)
      objfx(n,nn) = densitp(nn)/dt*(objru(n,nn) - obju(n,nn))
      objfy(n,nn) = densitp(nn)/dt*(objrv(n,nn) - objv(n,nn))
      objq(n,nn) = cpp(nn)*densitp(nn)/dt*(objrt(n,nn) - objt(n,nn))
    ENDDO
  ENDDO
ELSE
  objtp = zero
  DO nn = 1,nsphere
    DO n = 1,nobjcells(nn)
      objru(n,nn) = objcentu(nn) - objcentom(nn)*(objcelly(n,nn) - objcenty(nn)+(zobjcelly(n,nn) - zobjcenty(nn))*LDomainy)
      objrv(n,nn) = objcentv(nn) + objcentom(nn)*(objcellx(n,nn) - objcentx(nn)+(zobjcellx(n,nn) - zobjcentx(nn))*LDomainx)
      objfx(n,nn) = densitp(nn)/dt*(objru(n,nn) - obju(n,nn))
      objfy(n,nn) = densitp(nn)/dt*(objrv(n,nn) - objv(n,nn))
      objq(n,nn) = objqp(nn)
      objtp(nn) = objtp(nn) + objt(n,nn)*objcellvol(n,nn)
    ENDDO
    objtp(nn) = objtp(nn)/objvol(nn)
  ENDDO
ENDIF  
DO nn = 1,nsphere
  DO n = 1,nobjcells(nn)
    IF(objcell_bndFlag(n,nn) == 0)THEN
      DO i = objpoint_interpx(1,n,nn),objpoint_interpx(2,n,nn)
        ii = i + xPeriodic*(MAX(2-i,0)*nim2 + MIN(nim-i,0)*nim2)
        xp = objcellx(n,nn) + xPeriodic*(MAX(2-i,0)*LDomainx + MIN(nim-i,0)*LDomainx)
        dx = x(ii) - x(ii-1)
        DO j =objpoint_interpy(1,n,nn),objpoint_interpy(2,n,nn)
          jj = j + yPeriodic*(MAX(2-j,0)*njm2 + MIN(njm-j,0)*njm2)
          yp = objcelly(n,nn) + yPeriodic*(MAX(2-j,0)*LDomainy + MIN(njm-j,0)*LDomainy)
          dy = y(jj) - y(jj-1)
          del = fd_calc_delta(xc(ii),yc(jj),xp,yp,dx,dy)*objcellvol(n,nn)
          ij = li(ii)+jj
          !IF(rigidForce_contrib(ij) == nn)THEN
            fdsuc(ij) = fdsuc(ij) + objfx(n,nn)*del*dx*dy
            fdsvc(ij) = fdsvc(ij) + objfy(n,nn)*del*dx*dy
            fdstc(ij) = fdstc(ij) + objq(n,nn)*del*dx*dy
            !ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDDO
ENDDO 

IF(predict_or_correct==force_predict)THEN

  fdsu = fdsuc
  fdsv = fdsvc
  fdst = fdstc
  resor = zero

ELSEIF(predict_or_correct==force_correct)THEN

  maxresor = zero
  IF(iter > 1)THEN
    DO nn = 1,nsphere
      DO n = 1,nobjcells(nn)
        IF(objcell_bndFlag(n,nn) == 0)THEN
          DO i = objpoint_interpx(1,n,nn),objpoint_interpx(2,n,nn)
            ii = i + xPeriodic*(MAX(2-i,0)*nim2 + MIN(nim-i,0)*nim2)
            DO j = objpoint_interpy(1,n,nn),objpoint_interpy(2,n,nn)
              jj = j + yPeriodic*(MAX(2-j,0)*njm2 + MIN(njm-j,0)*njm2)
              curresor = ABS(fdsuc(li(ii)+jj)/fdsu(li(ii)+jj)) 
              IF(curresor > maxresor)THEN
                maxresor = curresor
!                imax=i
!                jmax=j
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  resor = maxresor    

  fdsu = fdsu + fd_urf*fdsuc
  fdsv = fdsv + fd_urf*fdsvc
  IF(isotherm)THEN
    fdst = fdst + fd_urf*fdstc
  ELSE
    fdst = fdstc
  ENDIF
ENDIF
ENDIF


END SUBROUTINE fd_calc_sources

SUBROUTINE fd_calc_rigidvel

USE shared_data,        ONLY : objcentu,objcentv,objcentom,nobjcells,objcellx,objcelly,&
                               objcentx,objcenty,objcellvol,obju,objv,nsphere,densitp,objvol,&
                               objcentmi,up,vp,omp,stationary,forcedmotion,dt,Fpq,LDomainx,LDomainy,&
                               zobjcellx,zobjcelly,zobjcentx,zobjcenty,objcell_bndFlag
USE precision,          ONLY : r_single
USe real_parameters,    ONLY : zero
IMPLICIT NONE
!-------------------------------Locals
INTEGER                        :: nn,n
REAL(KIND = r_single)          :: rxmc,rymc

IF(stationary .OR. forcedmotion)THEN
  objcentu = up
  objcentv = vp
  objcentom = omp
ELSE
  DO nn = 1,nsphere
    objcentu(nn) = zero
    objcentv(nn) = zero
    objcentom(nn) = zero
    DO n= 1,nobjcells(nn)
      IF(objcell_bndFlag(n,nn) == 0)THEN
        rxmc = objcellx(n,nn) - objcentx(nn) + (zobjcellx(n,nn) - zobjcentx(nn))*LDomainx
        rymc = objcelly(n,nn) - objcenty(nn) + (zobjcelly(n,nn) - zobjcenty(nn))*LDomainy
        objcentu(nn) = objcentu(nn) + objcellvol(n,nn) * (densitp(nn) * obju(n,nn) )
        objcentv(nn) = objcentv(nn) + objcellvol(n,nn) * (densitp(nn) * objv(n,nn) )
        objcentom(nn) = objcentom(nn) + objcellvol(n,nn) * densitp(nn) * (rxmc*objv(n,nn) - rymc*obju(n,nn)) 
      ENDIF
    ENDDO
    !--Convert Momentum to velocity
    objcentu(nn) = objcentu(nn) / (objvol(nn) * densitp(nn)) !+ Fpq(1,nn)*dt
    objcentv(nn) = objcentv(nn) / (objvol(nn) * densitp(nn)) !+ Fpq(2,nn)*dt
    objcentom(nn) = objcentom(nn) / objcentmi(nn)
  ENDDO
ENDIF

END SUBROUTINE fd_calc_rigidvel

SUBROUTINE fd_calc_physprops(itr,outvolp)

USE precision,          ONLY : r_single
USE parameters,         ONLY : OUTER_ITR_DONE
USE shared_data,        ONLY : nij,x,y,objpoint_cvx,objpoint_cvy,nsphere,nobjcells,objcellvol,den,deno,denoo,densit,&
                               xc,yc,objcellx,objcelly,li,nim,njm,densitp,objpoint_interpx,objpoint_interpy,&
                               celbeta,beta,betap,celkappa,celcp,celcpo,celcpoo,cpf,cpp,kappaf,kappap,lcal,ien,&
                               xPeriodic,yPeriodic,LDomainx,LDomainy,objcell_bndFlag
USE real_parameters,    ONLY : zero,one
IMPLICIT NONE

INTEGER,INTENT(IN)     :: itr

INTEGER                :: j,i,n,nn,ij,ii,jj,nim2,njm2
REAL(KIND = r_single)  :: dx,dy,del,volp(nij),xp,yp
REAL(KIND = r_single),OPTIONAL  :: outvolp(nij,nsphere)

nim2 = nim - 1
njm2 = njm - 1

IF(itr ==  OUTER_ITR_DONE)THEN

  denoo= deno
  deno = den

  celcpoo= celcpo
  celcpo = celcp

ENDIF

den = densit

celcp = cpf*densit

celkappa = kappaf
celbeta = beta*densit

DO nn = 1,nsphere
  volp = zero
  DO n = 1,nobjcells(nn)
    IF(objcell_bndFlag(n,nn) == 0)THEN
      DO i = objpoint_interpx(1,n,nn),objpoint_interpx(2,n,nn)
        ii = i + xPeriodic*(MAX(2-i,0)*nim2 + MIN(nim-i,0)*nim2)
        xp = objcellx(n,nn) + xPeriodic*(MAX(2-i,0)*LDomainx + MIN(nim-i,0)*LDomainx)
        dx = x(ii) - x(ii-1)
        DO j =objpoint_interpy(1,n,nn),objpoint_interpy(2,n,nn)
          jj = j + yPeriodic*(MAX(2-j,0)*njm2 + MIN(njm-j,0)*njm2) 
          ij = li(ii) + jj
          !IF(rigidForce_contrib(ij) == nn)THEN
            yp = objcelly(n,nn) + yPeriodic*(MAX(2-j,0)*LDomainy + MIN(njm-j,0)*LDomainy)
            dy = y(jj) - y(jj-1)
            del = fd_calc_delta(xc(ii),yc(jj),xp,yp,dx, dy)*objcellvol(n,nn)
            volp(ij) = volp(ij) + del
            !ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  IF(PRESENT(outvolp))THEN
    outvolp(:,nn) = volp(:)
  ENDIF
  IF(.NOT.lcal(ien))THEN
    DO i=2,nim
      DO j=2,njm
        ij=li(i)+j
        IF(volp(ij) /= zero)THEN
          den(ij) = (one - volp(ij))*densit + volp(ij)*densitp(nn)
        ENDIF
      ENDDO
    ENDDO
  ELSE
    DO i=2,nim
      DO j=2,njm
        ij=li(i)+j
        IF(volp(ij) /= zero)THEN
          den(ij) = (one - volp(ij))*densit + volp(ij)*densitp(nn)
          celbeta(ij) = (one - volp(ij))*beta*densit + volp(ij)*betap(nn)*densitp(nn)
          celcp(ij) = (one - volp(ij))*cpf*densit + volp(ij)*cpp(nn)*densitp(nn)
          celkappa(ij) = (one - volp(ij))*kappaf + volp(ij)*kappap(nn)
        ENDIF
      ENDDO
    ENDDO
  ENDIF    
ENDDO

END SUBROUTINE fd_calc_physprops

SUBROUTINE fd_init_temp(tin)

USE precision,          ONLY : r_single
USE shared_data,        ONLY : nij,x,y,objpoint_cvx,objpoint_cvy,nsphere,nobjcells,objcellvol,objtp,&
                               xc,yc,objcellx,objcelly,li,nim,njm,objpoint_interpx,objpoint_interpy,&
                               xPeriodic,yPeriodic,LDomainx,LDomainy,t
USE real_parameters,    ONLY : zero,one
IMPLICIT NONE

INTEGER                :: j,i,n,nn,ij,ii,jj,nim2,njm2
REAL(KIND = r_single)  :: dx,dy,del,volp(nij),xp,yp
REAL(KIND = r_single),INTENT(IN) :: tin

nim2 = nim - 1
njm2 = njm - 1

DO nn = 1,nsphere
  volp = zero
  DO n = 1,nobjcells(nn)
      DO i = objpoint_interpx(1,n,nn),objpoint_interpx(2,n,nn)
        ii = i + xPeriodic*(MAX(2-i,0)*nim2 + MIN(nim-i,0)*nim2)
        xp = objcellx(n,nn) + xPeriodic*(MAX(2-i,0)*LDomainx + MIN(nim-i,0)*LDomainx)
        dx = x(ii) - x(ii-1)
        DO j =objpoint_interpy(1,n,nn),objpoint_interpy(2,n,nn)
          jj = j + yPeriodic*(MAX(2-j,0)*njm2 + MIN(njm-j,0)*njm2)
          ij = li(ii) + jj
          !IF(rigidForce_contrib(ij) == nn)THEN
            yp = objcelly(n,nn) + yPeriodic*(MAX(2-j,0)*LDomainy + MIN(njm-j,0)*LDomainy)
            dy = y(jj) - y(jj-1)
            del = fd_calc_delta(xc(ii),yc(jj),xp,yp,dx, dy)*objcellvol(n,nn)
            volp(ij) = volp(ij) + del
            !ENDIF
        ENDDO
      ENDDO
  ENDDO
  DO i=2,nim
    DO j=2,njm
      ij=li(i)+j
      IF(volp(ij) /= zero)THEN
        t(ij) = (one - volp(ij))*tin + volp(ij)*objtp(nn)
      ENDIF
    ENDDO
  ENDDO   
ENDDO

END SUBROUTINE fd_init_temp

SUBROUTINE fd_calc_mi

!--Calculates the moment of inertia tensor for a general object defined by a centre
!--point, and discretised to set of control volumes. I_p = sum_{n=1}^N (rho_n*V_n*[(r_kr_k)sigma_{ij} - r_ir_j])
!--for a planar rotation it enough to calculate sum{n=1}^N (rho_n*V_n *(r_x^2+r_y^2) )
USE precision,        ONLY : r_single
USE shared_data,      ONLY : nsphere,nobjcells,densitp,objcellvol,objcentmi,objcentx,objcenty,objcellx,objcelly,&
                             LDomainx,LDomainy,zobjcellx,zobjcelly,zobjcentx,zobjcenty 
USE real_parameters,  ONLY : zero

IMPLICIT NONE
!-------------------------Locals
INTEGER :: nn,n
REAL(KIND = r_single)    :: rxcp,rycp           !--CV centre - centroid vector

DO nn = 1,nsphere 
  objcentmi(nn) = zero
  DO n = 1,nobjcells(nn)
    rxcp = objcellx(n,nn) - objcentx(nn) + (zobjcellx(n,nn) - zobjcentx(nn))*LDomainx
    rycp = objcelly(n,nn) - objcenty(nn) + (zobjcelly(n,nn) - zobjcenty(nn))*LDomainy
    objcentmi(nn) = objcentmi(nn) + ( densitp(nn) * objcellvol(n,nn) * (rxcp*rxcp + rycp*rycp) )
  ENDDO
ENDDO

END SUBROUTINE fd_calc_mi

SUBROUTINE fd_calc_ori

!--Calculate orientation tensor, in 2D its a 2x2 matrix [cos(a) -sin(a); sin(a) cos(a)]
USE precision,        ONLY : r_single
USE shared_data,      ONLY : nsphere,objcentom,objcento,dt
USE real_parameters,  ONLY : zero

IMPLICIT NONE
!-------------------------Locals
INTEGER :: nn
REAL(KIND = r_single) :: alpha,c,s,sigma

DO nn = 1,nsphere
  !--angular velocity vector is (0,0,omega)
  alpha = ABS(objcentom(nn))*dt
  IF(ABS(objcentom(nn)) > zero)THEN
    sigma = objcentom(nn)/ABS(objcentom(nn)) 
  ELSE
    sigma = zero
  ENDIF
  c     = COS(alpha)
  s     = SIN(alpha)
  objcento(1,nn) = c
  objcento(2,nn) = -sigma*s
  objcento(3,nn) = sigma*s
  objcento(4,nn) = c
ENDDO

END SUBROUTINE fd_calc_ori

SUBROUTINE fd_move_mesh(ismoved)

USE precision,        ONLY : r_single
USE shared_data,      ONLY : dxmeanmoved,dxmeanmovedtot,dxmean,u,v,p,t,&
                             to,uo,vo,too,uoo,voo, &
                              surfpointx,objcentx,objcentxo,objcellvertx,&
                              objpoint_cvx,objpoint_cvy,ni,nj,lcal,iu,ip,ien,&
                              li,f1,f2,nsphere,objcellx,objcelly,nobjcells,&
                              objpoint_interpx,objpoint_interpy,nij,celbeta,celcp,celkappa,&
                              den,deno,objcell_bndFlag
USE real_parameters,  ONLY : zero,three,two,half,one
USE modfd_set_bc,     ONLY : fd_bcout

IMPLICIT NONE
LOGICAL       :: ismoved
INTEGER       :: i,j,ij,n,nn,iip,jp,error
REAL(KIND = r_single) :: var(nij),varo(nij),varoo(nij)
error = 0
ismoved = .FALSE.
IF(ABS(dxmeanmoved) > dxmean)THEN
  ismoved = .TRUE.
  dxmeanmovedtot = dxmeanmovedtot + dxmean
  dxmeanmoved = zero
  IF(lcal(iu))THEN
    var = u
    varo = uo
    varoo = uoo
    DO i = 3,ni-1 !--Assume second line remains the same?
      DO j = 1,nj
        ij = li(i) + j
        u(ij) = var(ij - nj)
        uo(ij) = varo(ij - nj)
        uoo(ij) = varoo(ij - nj)
      ENDDO
    ENDDO
    var = v
    varo = vo
    varoo = voo
    DO i = 3,ni !--Assume second line remains the same?
      DO j = 1,nj
        ij = li(i) + j
        v(ij) = var(ij - nj)
        vo(ij) = varo(ij - nj)
        voo(ij) = varoo(ij - nj)
      ENDDO
    ENDDO
    var = f1
    DO i = 3,ni-2 !--Assume second line remains the same?
      DO j = 1,nj
        ij = li(i) + j
        f1(ij) = var(ij - nj)
      ENDDO
    ENDDO
    CALL fd_bcout
    var = f2
    DO i = 3,ni-1 !--Assume second line remains the same?
      DO j = 1,nj
        ij = li(i) + j
        f2(ij) = var(ij - nj)
      ENDDO
    ENDDO
  ENDIF
  IF(lcal(ip))THEN
    var = p
    DO i = 3,ni !--Assume second line remains the same?
      DO j = 1,nj
        ij = li(i) + j
        p(ij) = var(ij - nj)
      ENDDO
    ENDDO
  ENDIF
  IF(lcal(ien))THEN
    var = t
    varo = to
    varoo = too
    DO i = 3,ni !--Assume second line remains the same?
      DO j = 1,nj
        ij = li(i) + j
        t(ij) = var(ij - nj)
        to(ij) = varo(ij - nj)
        too(ij) = varoo(ij - nj)
      ENDDO
    ENDDO
  ENDIF
  var = den
  varo = deno
  DO i = 3,ni !--Assume second line remains the same?
    DO j = 1,nj
      ij = li(i) + j
      den(ij) = var(ij - nj)
      deno(ij) = varo(ij - nj) 
    ENDDO
  ENDDO
  var = celbeta
  varo = celcp
  DO i = 3,ni !--Assume second line remains the same?
    DO j = 1,nj
      ij = li(i) + j
      celbeta(ij) = var(ij - nj)
      celcp(ij) = varo(ij - nj) 
    ENDDO
  ENDDO
  varo = celkappa
  DO i = 3,ni !--Assume second line remains the same?
    DO j = 1,nj
      ij = li(i) + j
      celkappa(ij) = varo(ij - nj) 
    ENDDO
  ENDDO
  DO nn = 1,nsphere
    objcentx(nn) = objcentx(nn) + dxmean
    objcentxo(nn) = objcentxo(nn) + dxmean
  
    DO n = 1,nobjcells(nn)
      objcellx(n,nn) = objcellx(n,nn) + dxmean

      CALL fd_track_single_point(objcell_bndFlag(n,nn),objcellx(n,nn),objcelly(n,nn),objpoint_cvx(n,nn),objpoint_cvy(n,nn),iip,jp,&
                               objpoint_interpx(:,n,nn),objpoint_interpy(:,n,nn),error)
      IF(error /= 0)GOTO 100

      IF(objpoint_cvx(n,nn) /= iip - 1 .OR. objpoint_cvy(n,nn) /= jp)THEN
        WRITE(*,*)'Seems somthing is wrong.'
      ENDIF

      IF(objpoint_cvx(n,nn) /= iip .OR. objpoint_cvy(n,nn) /= jp)THEN
        objpoint_cvx(n,nn) = iip
        objpoint_cvy(n,nn) = jp      
      ENDIF

      objcellvertx(:,n,nn) =  objcellvertx(:,n,nn) + dxmean

    ENDDO
    !--Only rotate it (normals, centers, interpolation stencil etc... &
    !--not updated if want to explicitly calculate hyd forces add them)
    surfpointx(:,nn) = surfpointx(:,nn) + dxmean
  ENDDO   
ENDIF

100 CONTINUE
IF(error /= 0)WRITE(*,*)'fd_move_mesh: particle lost'

END SUBROUTINE fd_move_mesh

SUBROUTINE fd_calc_pos(itr,titr)

USE parameters,       ONLY : OUTER_ITR_DONE
USE precision,        ONLY : r_single
USE shared_data,      ONLY : nsphere,nobjcells,objcentmi,objcentxo,objcentyo,objcellx,objcelly,&
                             objcentxo,objcentyo,dt,objcentu,objcentv,objcento,objpoint_interpx,objpoint_interpy,&
                             objcellvertx,objcellverty,objcentx,objcenty,objcentuo,objcentvo,calcsurfforce,&
                             objpoint_cvx,objpoint_cvy,nsurfpoints,surfpointx,surfpointy,dxmeanmoved,up,vp,omp,&
                             forcedmotion,lread,Fpq,Fpw,objvol,densitp,x,y,xc,yc,objcellvol,li,fdfcu,fdfcv,&
                             zobjcentx,zobjcentxo,zobjcenty,zobjcentyo,zobjcellx,zobjcelly,LDomainx,LDomainy,&
                             zobjcellvertx,zobjcellverty,nim,njm,zsurfpointx,zsurfpointy,objcell_bndFlag,subTimeStep
USE real_parameters,  ONLY : zero,three,two,half,one,four

IMPLICIT NONE

INTEGER,INTENT(IN)    :: itr,titr

INTEGER               :: n,nn,ip,jp,error,moved,k,m
REAL(KIND = r_single) :: lindispx,lindispy,rmcx,rmcy,rmcyr,rmcxr,dxmoved(nsphere),&
                          dtk,xShift,yShift,objcentuk(nsphere),objcentvk(nsphere),&
                          diffu,diffv,dumfx,dumfy
REAL(KIND = r_single),ALLOCATABLE,SAVE      :: tempFpq(:,:),tempFpw(:,:) 
REAL(KIND = r_single),ALLOCATABLE,SAVE :: temp2Fpq(:,:) ,temp2Fpw(:,:)
REAL(KIND = r_single),ALLOCATABLE,SAVE      :: objcellxo(:,:),objcellyo(:,:)
INTEGER,ALLOCATABLE,SAVE                    :: zobjcellxo(:,:),zobjcellyo(:,:)


IF(titr == 1 )THEN !.AND. itr ==1 
ALLOCATE( objcellxo(MAXVAL(nobjcells,1),nsphere),objcellyo(MAXVAL(nobjcells,1),nsphere),&
          zobjcellxo(MAXVAL(nobjcells,1),nsphere),zobjcellyo(MAXVAL(nobjcells,1),nsphere),&
          tempFpq(2,nsphere),tempFpw(2,nsphere),temp2Fpq(2,nsphere),temp2Fpw(2,nsphere)  )

tempFpq = zero 
tempFpq = zero
tempFpq = zero
tempFpq = zero

ENDIF

IF(forcedmotion)THEN
  IF(itr == OUTER_ITR_DONE)THEN
  DO nn = 1,nsphere
    !--Move the particle 
    objcentxo(nn) = objcentx(nn)
    objcentyo(nn) = objcenty(nn) 
    CALL fd_calc_forcedmotion(nn)
    CALL fd_update_forcedvel(nn)
    IF(calcsurfforce)CALL fd_create_interpmol(nn)
    DO n = 1,nobjcells(nn)
   
      CALL fd_track_single_point(objcell_bndFlag(n,nn),objcellx(n,nn),objcelly(n,nn),objpoint_cvx(n,nn),objpoint_cvy(n,nn),ip,jp,&
                                 objpoint_interpx(:,n,nn),objpoint_interpy(:,n,nn),error)
      IF(error /= 0)GOTO 100

      IF(objpoint_cvx(n,nn) /= ip .OR. objpoint_cvy(n,nn) /= jp)THEN
        objpoint_cvx(n,nn) = ip
        objpoint_cvy(n,nn) = jp      
      ENDIF
    ENDDO
  ENDDO
  ENDIF
ELSE
      
  objcentxo(1:nsphere) = objcentx(1:nsphere)
  objcentyo(1:nsphere) = objcenty(1:nsphere)
  zobjcentxo(1:nsphere)= zobjcentx(1:nsphere)
  zobjcentyo(1:nsphere)= zobjcenty(1:nsphere)

  dtk = dt/REAL(subTimeStep,r_single)
  
  DO k = 1,subTimeStep

    tempFpq(1:2,1:nsphere) = Fpq(1:2,1:nsphere)
    tempFpw(1:2,1:nsphere) = Fpw(1:2,1:nsphere)

    DO nn = 1,nsphere    
         
        lindispx = objcentu(nn) * dtk
        lindispy = objcentv(nn) * dtk

      objcentx(nn) = objcentx(nn) + lindispx
      objcenty(nn) = objcenty(nn) + lindispy

      IF(objcentx(nn) > x(nim))THEN
        objcentx(nn) = objcentx(nn) - LDomainx
        zobjcentx(nn) = zobjcentx(nn) + 1
      ELSEIF(objcentx(nn) < x(1))THEN
        objcentx(nn) = objcentx(nn) + LDomainx 
        zobjcentx(nn) = zobjcentx(nn) - 1
      ENDIF

      IF(objcenty(nn) > y(njm))THEN
        objcenty(nn) = objcenty(nn) - LDomainy
        zobjcenty(nn) = zobjcenty(nn) + 1
      ELSEIF(objcenty(nn) < y(1))THEN
        objcenty(nn) = objcenty(nn) + LDomainy
        zobjcenty(nn) = zobjcenty(nn) - 1
      ENDIF

    ENDDO

      CALL fd_calc_part_collision(Fpq,Fpw)

    DO nn = 1,nsphere

        dumfx = tempFpq(1,nn)+Fpq(1,nn)+tempFpw(1,nn)+Fpw(1,nn)
        dumfy = tempFpq(2,nn)+Fpq(2,nn)+tempFpw(2,nn)+Fpw(2,nn)

        lindispx = one/objvol(nn)/densitp(nn)/four*dumfx*dtk**2
        lindispy = one/objvol(nn)/densitp(nn)/four*dumfy*dtk**2

        diffu = one/objvol(nn)/densitp(nn)/two*dumfx*dtk
        diffv = one/objvol(nn)/densitp(nn)/two*dumfy*dtk

        objcentu(nn) = objcentu(nn) + diffu
        objcentv(nn) = objcentv(nn) + diffv

      objcentx(nn) = objcentx(nn) + lindispx
      objcenty(nn) = objcenty(nn) + lindispy

      IF(objcentx(nn) > x(nim))THEN
        objcentx(nn) = objcentx(nn) - LDomainx
        zobjcentx(nn) = zobjcentx(nn) + 1
      ELSEIF(objcentx(nn) < x(1))THEN
        objcentx(nn) = objcentx(nn) + LDomainx 
        zobjcentx(nn) = zobjcentx(nn) - 1
      ENDIF

      IF(objcenty(nn) > y(njm))THEN
        objcenty(nn) = objcenty(nn) - LDomainy
        zobjcenty(nn) = zobjcenty(nn) + 1
      ELSEIF(objcenty(nn) < y(1))THEN
        objcenty(nn) = objcenty(nn) + LDomainy
        zobjcenty(nn) = zobjcenty(nn) - 1
      ENDIF

    ENDDO  
  
  ENDDO

    !CALL fd_calc_cellTemp

  DO nn = 1,nsphere
    
    moved = 0
    
    lindispx = objcentx(nn) - objcentxo(nn) + (zobjcentx(nn) - zobjcentxo(nn))*LDomainx
    lindispy = objcenty(nn) - objcentyo(nn) + (zobjcenty(nn) - zobjcentyo(nn))*LDomainy
    
      objcentu(nn) = lindispx/dt
      objcentv(nn) = lindispy/dt
                             
    !--Update the material points
    !--orientation updated in a seperate subroutine
    DO n = 1,nobjcells(nn)
      
      xShift = (zobjcellx(n,nn) - zobjcentxo(nn))*LDomainx
      yShift = (zobjcelly(n,nn) - zobjcentyo(nn))*LDomainy

      rmcx = objcellx(n,nn) - objcentxo(nn) + xShift
      rmcy = objcelly(n,nn) - objcentyo(nn) + yShift
    
      objcellx(n,nn) = objcentxo(nn) + (objcento(1,nn)*rmcx + objcento(2,nn)*rmcy) + lindispx - xShift
      objcelly(n,nn) = objcentyo(nn) + (objcento(3,nn)*rmcx + objcento(4,nn)*rmcy) + lindispy - yShift
      
      IF(objcellx(n,nn) > x(nim))THEN
        objcellx(n,nn) = objcellx(n,nn) - LDomainx
        zobjcellx(n,nn) = zobjcellx(n,nn) + 1
      ELSEIF(objcellx(n,nn) < x(1))THEN
        objcellx(n,nn) = objcellx(n,nn) + LDomainx 
        zobjcellx(n,nn) = zobjcellx(n,nn) - 1
      ENDIF

      IF(objcelly(n,nn) > y(njm))THEN
        objcelly(n,nn) = objcelly(n,nn) - LDomainy
        zobjcelly(n,nn) = zobjcelly(n,nn) + 1
      ELSEIF(objcelly(n,nn) < y(1))THEN
        objcelly(n,nn) = objcelly(n,nn) + LDomainy
        zobjcelly(n,nn) = zobjcelly(n,nn) - 1
      ENDIF

      CALL fd_track_single_point(objcell_bndFlag(n,nn),objcellx(n,nn),objcelly(n,nn),objpoint_cvx(n,nn),objpoint_cvy(n,nn),ip,jp,&
                                 objpoint_interpx(:,n,nn),objpoint_interpy(:,n,nn),error)
      IF(error /= 0)GOTO 100

      IF(objpoint_cvx(n,nn) /= ip .OR. objpoint_cvy(n,nn) /= jp)THEN
        !--Create the box
        moved = moved + 1
        objpoint_cvx(n,nn) = ip
        objpoint_cvy(n,nn) = jp      
      ENDIF

      DO m = 1,4
        
        xShift = (zobjcellvertx(m,n,nn) - zobjcentxo(nn))*LDomainx
        yShift = (zobjcellverty(m,n,nn) - zobjcentyo(nn))*LDomainy

        rmcxr = objcellvertx(m,n,nn) - objcentxo(nn) + xShift
        rmcyr = objcellverty(m,n,nn) - objcentyo(nn) + yShift

        objcellvertx(m,n,nn) =  objcentxo(nn) + (objcento(1,nn)*rmcxr + objcento(2,nn)*rmcyr) + lindispx - xshift
        objcellverty(m,n,nn) =  objcentyo(nn) + (objcento(3,nn)*rmcxr + objcento(4,nn)*rmcyr) + lindispy - yShift
      
        IF(objcellvertx(m,n,nn) > x(nim))THEN
          objcellvertx(m,n,nn) = objcellvertx(m,n,nn) - LDomainx
          zobjcellvertx(m,n,nn) = zobjcellvertx(m,n,nn) + 1
        ELSEIF(objcellvertx(m,n,nn) < x(1))THEN
          objcellvertx(m,n,nn) = objcellvertx(m,n,nn) + LDomainx 
          zobjcellvertx(m,n,nn) = zobjcellvertx(m,n,nn) - 1
        ENDIF

        IF(objcellverty(m,n,nn) > y(njm))THEN
          objcellverty(m,n,nn) = objcellverty(m,n,nn) - LDomainy
          zobjcellverty(m,n,nn) = zobjcellverty(m,n,nn) + 1
        ELSEIF(objcellverty(m,n,nn) < y(1))THEN
          objcellverty(m,n,nn) = objcellverty(m,n,nn) + LDomainy
          zobjcellverty(m,n,nn) = zobjcellverty(m,n,nn) - 1
        ENDIF

      ENDDO
 
    ENDDO
   
    DO n=1,nsurfpoints(nn)
      
      xShift = (zsurfpointx(n,nn) - zobjcentxo(nn))*LDomainx
      yShift = (zsurfpointy(n,nn) - zobjcentyo(nn))*LDomainy

      rmcx = surfpointx(n,nn) - objcentxo(nn) + xShift  
      rmcy = surfpointy(n,nn) - objcentyo(nn) + yShift

      surfpointx(n,nn) = objcentxo(nn) + (objcento(1,nn)*rmcx + objcento(2,nn)*rmcy) + lindispx - xShift
      surfpointy(n,nn) = objcentyo(nn) + (objcento(3,nn)*rmcx + objcento(4,nn)*rmcy) + lindispy - yShift

      IF(surfpointx(n,nn) > x(nim))THEN
        surfpointx(n,nn) = surfpointx(n,nn) - LDomainx
        zsurfpointx(n,nn) = zsurfpointx(n,nn) + 1
      ELSEIF(surfpointx(n,nn) < x(1))THEN
        surfpointx(n,nn) = surfpointx(n,nn) + LDomainx 
        zsurfpointx(n,nn) = zsurfpointx(n,nn) - 1
      ENDIF

      IF(surfpointy(n,nn) > y(njm))THEN
        surfpointy(n,nn) = surfpointy(n,nn) - LDomainy
        zsurfpointy(n,nn) = zsurfpointy(n,nn) + 1
      ELSEIF(surfpointy(n,nn) < y(1))THEN
        surfpointy(n,nn) = surfpointy(n,nn) + LDomainy
        zsurfpointy(n,nn) = zsurfpointy(n,nn) - 1
      ENDIF

    ENDDO

    IF(calcsurfforce)CALL fd_create_interpmol(nn)

    WRITE(*,*)moved, ' Particles moved in this time step for sphere, ', nn

  ENDDO
ENDIF

dxmoved = objcentx - objcentxo + (zobjcentx - zobjcentxo)*LDomainx

dxmeanmoved = dxmeanmoved + SUM(dxmoved) /REAL(nsphere)
100 CONTINUE
IF(error == 1)THEN
  WRITE(*,*)'fd_calc_pos: Node Lost in x-dir...'
ELSEIF(error == 2)THEN 
  WRITE(*,*)'fd_calc_pos: Node Lost in y-dir...' 
ENDIF

END SUBROUTINE fd_calc_pos

SUBROUTINE fd_copy_oldvel

USE shared_data,      ONLY : nsphere,objcentu,objcentv,objcentom,objcentuo,objcentvo,objcentomo

IMPLICIT NONE

objcentuo(1:nsphere) = objcentu(1:nsphere)
objcentvo(1:nsphere) = objcentv(1:nsphere)
objcentomo(1:nsphere) = objcentom(1:nsphere)
 
END SUBROUTINE fd_copy_oldvel

SUBROUTINE fd_geom_reader(ierror)

USE shared_data,        ONLY : objcellx,objcelly,objcellvol,objpoint_cvx,objpoint_cvy,objcellvertx,objcellverty,&
                               objpoint_interpx,objpoint_interpy,problem_name,problem_len,nsphere,objcentx,objcenty,&
                               nobjcells,objvol,objradius
USE mod_create_filenum, ONLY : create_filenum
USE real_parameters,    ONLY : third,quarter,half,pi
USE precision,          ONLY : r_single
IMPLICIT NONE
    INTEGER,INTENT(OUT) :: ierror

    CHARACTER(LEN=5) :: filenum(nsphere+1)
    CHARACTER(LEN=80)::dummy,geomName
    INTEGER :: NUMNP,NELEM,NGRPS,NBSETS,NDFCD,NDFVL,MaxNELEM,j,n,i,eType,np,n1,n2,n3,n4
    REAL(KIND = r_single),ALLOCATABLE:: nodePos(:,:)
    REAL(KIND = r_single)            :: x,y,vec1(2),vec2(2)

    ierror = 0

    CALL create_filenum(nsphere+1,filenum)

    MaxNelem = 0
    WRITE(*,*)'Note: Object will be translated to (objCentx,objcenty).' 
    DO n = 1,nsphere
      OPEN(UNIT = 1, FILE = problem_name(1:problem_len)//'_'//filenum(n+1)//'.neu', STATUS = 'OLD',IOSTAT=ierror)
      IF(ierror /= 0)GOTO 1001
      READ(1,101)GeomName
      READ(1,101)dummy
      READ(1,101)dummy
      READ(1,101)dummy
      READ(1,101)dummy
      READ(1,102)NUMNP,NELEM,NGRPS,NBSETS,NDFCD,NDFVL
      IF(NELEM > MaxNelem)MaxNelem = NELEM
      CLOSE(1)
    ENDDO

    ALLOCATE(objcellx(MaxNelem,nsphere),objcelly(MaxNelem,nsphere),&
           objcellvol(MaxNelem,nsphere),objpoint_cvx(MaxNelem,nsphere),&
           objpoint_cvy(MaxNelem,nsphere),objcellvertx(4,MaxNelem,nsphere),objcellverty(4,MaxNelem,nsphere),&
           objpoint_interpx(2,MaxNelem,nsphere),objpoint_interpy(2,MaxNelem,nsphere)) 
    DO n = 1,nsphere
      OPEN(UNIT = 1, FILE = problem_name(1:problem_len)//'_'//filenum(n+1)//'.neu', STATUS = 'OLD',IOSTAT=ierror)
      READ(1,101)GeomName
      READ(1,101)dummy
      READ(1,101)dummy
      READ(1,101)dummy
      READ(1,101)dummy
      READ(1,102)NUMNP,nobjcells(n),NGRPS,NBSETS,NDFCD,NDFVL
      READ(1,101)dummy
      IF(dummy .ne. 'ENDOFSECTION')THEN
        WRITE(*,*)'ERROR : Geometry Format Problem!'
        ierror = 1
        RETURN
      END IF
      
      dummy = ''

      READ(1,101)dummy

      IF(ADJUSTL(dummy) .ne. 'NODAL COORDINATES 2.3.16')THEN
        WRITE(*,*)'ERROR : Geometry Format Problem!'
        ierror = 1
        RETURN
      END IF

      dummy = ''

      ALLOCATE(nodePos(NDFVL,NUMNP))
      DO i=1,NUMNP 
        READ(1,103)j,x,y
        nodePos(1,j) = x
        nodePos(2,j) = y
      ENDDO

      READ(1,101)dummy
      IF(dummy .ne. 'ENDOFSECTION')THEN
        WRITE(*,*)'ERROR : Geometry Format Problem!'
        ierror = 1
        RETURN
      END IF

      dummy = ''

      READ(1,101)dummy

      IF(ADJUSTL(dummy) .NE. 'ELEMENTS/CELLS 2.3.16')THEN
        WRITE(*,*)'ERROR : Geometry Format Problem!'
        ierror = 1
        RETURN
      END IF

      dummy = ''

      DO i=1,NELEM  
        READ(1,1041)j,eType,nP
        IF(eType == 3)THEN
          IF(np /= 3)THEN
            WRITE(*,*)'Error in element type!'
            ierror = 1
            RETURN
          ENDIF
          BACKSPACE(1)
          READ(1,104)j,eType,nP,n1,n2,n3
          objcellvertx(1,i,n) = nodePos(1,n1) + objcentx(n)
          objcellvertx(2,i,n) = nodePos(1,n2) + objcentx(n)
          objcellvertx(3,i,n) = nodePos(1,n3) + objcentx(n)
          objcellvertx(4,i,n) = nodePos(1,n3) + objcentx(n)

          objcellverty(1,i,n) = nodePos(2,n1) + objcenty(n)
          objcellverty(2,i,n) = nodePos(2,n2) + objcenty(n)
          objcellverty(3,i,n) = nodePos(2,n3) + objcenty(n)
          objcellverty(4,i,n) = nodePos(2,n3) + objcenty(n)

          objcellx(i,n) = third * (objcellvertx(1,i,n) + objcellvertx(2,i,n) + objcellvertx(3,i,n))
          objcelly(i,n) = third * (objcellverty(1,i,n) + objcellverty(2,i,n) + objcellverty(3,i,n))

          vec1 = (/objcellvertx(2,i,n) - objcellvertx(1,i,n), objcellverty(2,i,n) - objcellverty(1,i,n)/)
          vec2 = (/objcellvertx(3,i,n) - objcellvertx(1,i,n), objcellverty(3,i,n) - objcellverty(1,i,n)/)

          objcellvol(i,n) = half * ABS(vec1(1)*vec2(2) - vec1(2)*vec2(1))

          objvol(n) = objvol(n) + objcellvol(i,n)

        ELSEIF(eType == 2)THEN
          IF(np /= 4)THEN
            WRITE(*,*)'Error in element type!'
            ierror = 1
            RETURN
          ENDIF
          BACKSPACE(1)
          READ(1,104)j,eType,nP,n1,n2,n3,n4
          objcellvertx(1,i,n) = nodePos(1,n1) + objcentx(n)
          objcellvertx(2,i,n) = nodePos(1,n2) + objcentx(n)
          objcellvertx(3,i,n) = nodePos(1,n3) + objcentx(n)
          objcellvertx(4,i,n) = nodePos(1,n4) + objcentx(n)

          objcellverty(1,i,n) = nodePos(2,n1) + objcenty(n)
          objcellverty(2,i,n) = nodePos(2,n2) + objcenty(n)
          objcellverty(3,i,n) = nodePos(2,n3) + objcenty(n)
          objcellverty(4,i,n) = nodePos(2,n4) + objcenty(n)

          objcellx(i,n) = quarter * (objcellvertx(1,i,n) + objcellvertx(2,i,n) + &
                                     objcellvertx(3,i,n) + objcellvertx(4,i,n))
          objcelly(i,n) = quarter * (objcellverty(1,i,n) + objcellverty(2,i,n) + &
                                     objcellverty(3,i,n) + objcellverty(4,i,n))

          vec1 = (/objcellvertx(2,i,n) - objcellvertx(1,i,n), objcellverty(2,i,n) - objcellverty(1,i,n)/)
          vec2 = (/objcellvertx(4,i,n) - objcellvertx(1,i,n), objcellverty(4,i,n) - objcellverty(1,i,n)/)

          objcellvol(i,n) = ABS(vec1(1)*vec2(2) - vec1(2)*vec2(1))

          objvol(n) = objvol(n) + objcellvol(i,n)

        ELSE
          WRITE(*,*)'Element type not supported'
          ierror = 1
          RETURN
        ENDIF
      ENDDO 

      READ(1,101)dummy
      IF(dummy .ne. 'ENDOFSECTION')THEN
        WRITE(*,*)'ERROR : Geometry Format Problem!'
        ierror = 1
        RETURN
      END IF
      CLOSE(1)
      
      DEALLOCATE(nodePos)

      WRITE(*,*)'Real Area for Obj: ',n,' is: ',pi*objradius(n)**2,', estimated: ', objvol(n),'.' 
    ENDDO    
    
1001 CONTINUE
    IF(ierror /= 0)THEN
      WRITE(*,*)'Object file not found.'
      RETURN
    ENDIF    
    100 FORMAT(A)
    101 FORMAT(A80)
    102 FORMAT(/6(1X,I9))
    103 FORMAT(I10,3E20.11)
    104 FORMAT(I8,1X,I2,1X,I2,1X,7I8:/(15X,7I8:))
    1041 FORMAT(I8,1X,I2,1X,I2,1X)
    105 FORMAT('GROUP:',I11,' ELEMENTS: ',I11,'MATERIAL: ',I11,' NFLAGS:',I11)
    106 FORMAT(A36)
    107 FORMAT(10I8)
    108 FORMAT(A120)
    109 FORMAT(A32, 8I10)
    
END SUBROUTINE fd_geom_reader 

!--Some Auxiliary Subroutines for timing etc...
SUBROUTINE fd_calc_quality(volp,q)

USE precision,        ONLY : r_single
USE shared_data,      ONLY : celltype,nim,njm,nsphere,li,nij
USE real_parameters,  ONLY : one,zero
IMPLICIT NONE

INTEGER               :: ij,i,j,n,nn
REAL(KIND = r_single) :: volp(nij,nsphere),q(nsphere)

q = zero
DO n = 1,nsphere
  nn = 0
  DO i = 2,nim
    DO j = 2,njm
      ij = li(i) + j
      IF(celltype(ij) == n)THEN
        nn = nn + 1
        q(n) = q(n) + (one - volp(ij,n))**2
      ENDIF
    ENDDO
  ENDDO
  q(n) = one - SQRT(q(n) / REAL(nn,r_single))
ENDDO

END SUBROUTINE fd_calc_quality

SUBROUTINE fd_calc_forcedmotion(n)

!--used for forced motion

USE precision,              ONLY : r_single
USE real_parameters,        ONLY : pi,fifth,five,two,zero,six,ten,one
USE shared_data,            ONLY : up,vp,omp,nsurfpoints,nobjcells,&
                                   objcentx,surfpointx,objcenty,surfpointy,objradius,time,&
                                   objcentxinit,surfpointxinit,objcellx,objcellxinit,&
                                   objcellvertxinit,objcellvertx,&
                                   objcentyinit,surfpointyinit,objcelly,objcellyinit,&
                                   objcellvertyinit,objcellverty

IMPLICIT NONE

INTEGER,INTENT(IN)    :: n
INTEGER               :: j

REAL(KIND = r_single) :: omega,amp,freq,KC,dx,starttime,t,freq0,freqfact,dy

!--horizontal motion in stationary fluid
starttime = fifth
freq0 = 0.95_r_single; freqfact = 1.1_r_single
freq = freq0 * freqfact
amp = 0.4_r_single * objradius(n)
omega = two*pi*freq
t = time - starttime
IF(t < zero) t = zero 
dy = amp*SIN(omega * t)
objcenty(n) = objcentyinit(n) - dy
DO j=1,nobjcells(n) 
  objcelly(j,n) = objcellyinit(j,n) - dy
  objcellverty(:,j,n) = objcellvertyinit(:,j,n) - dy
ENDDO

DO j=1,nsurfpoints(n) 
  surfpointy(j,n) = surfpointyinit(j,n) - dy
ENDDO 

!freq = fifth
!KC   = five
!omega = two*pi*freq
!amp = objradius(n) * KC / pi
!
!dx = amp*SIN(omega * time)
! 
!objcentx(n) = objcentxinit(n) - dx
!
!DO j=1,nobjcells(n) 
!  objcellx(j,n) = objcellxinit(j,n) - dx
!  objcellvertx(:,j,n) = objcellvertxinit(:,j,n) - dx
!ENDDO
!
!DO j=1,nsurfpoints(n) 
!  surfpointx(j,n) = surfpointxinit(j,n) - dx
!ENDDO 

END SUBROUTINE fd_calc_forcedmotion

SUBROUTINE fd_update_forcedvel(n)
!--used for forced motion
USE shared_data,            ONLY : objcentx,objcentxo,up,vp,dt,objcenty,objcentyo
IMPLICIT NONE

INTEGER,INTENT(IN)   :: n

up(n) = (objcentx(n)-objcentxo(n))/dt
vp(n) = (objcenty(n)-objcentyo(n))/dt
  
END SUBROUTINE fd_update_forcedvel

SUBROUTINE fd_create_ellipse_surf(x,y,n,majr,minr)

!--Evenly distributes n points on the surface of an ellipse

USE precision,      ONLY : r_single
USE real_parameters,ONLY : small,pi,one,four,three,ten,two,degtorad,zero,large,vsmall,radtodeg
IMPLICIT NONE

INTEGER                               :: n         !--number of required points
REAL(KIND = r_single),INTENT(INOUT)   :: x(n),y(n) !--x, y coordinated of the equidistant points
REAL(KIND = r_single),INTENT(IN)      :: majr,minr !--minor and major radii
REAL(KIND = r_single)                 :: ecc,Stot,s(n),aspratio,theta(n),thetao(n),tol,costheta,sintheta
DOUBLE PRECISION                      :: eccD, StotD,FE,EE,mtheta,EED
INTEGER                               :: i,j

!--Eccentricity
ecc = SQRT(one-minr**2/majr**2)
!--Calculate total area
aspratio = ((majr - minr)/(majr+minr))**2
Stot = pi*(majr + minr)*(one + (three*aspratio)/(ten+SQRT(four - three*aspratio)))

!--Test ELIT 
eccD = DBLE(ecc)
mtheta = 90.D0
CALL ELIT(eccD, mtheta,FE,EE)
StotD = four * majr * EE
IF( ABS(StotD - DBLE(Stot)) > 1.D-6)THEN
  WRITE(*,*)'Tot area from Incomplete elliptic:', StotD, 'Estimate Eq:' ,Stot
ENDIF

!--Required Arc Lengths
DO i = 1,n
  s(i) = REAL(i-1,r_single)*Stot/REAL(n,r_single)
ENDDO

!--Find an angle corresponding to an arc lengths (Cumulative values)
DO i= 1,n/4
  tol = large
  EE = large
  theta(i) = zero
  DO WHILE(tol > vsmall .AND. ABS(EE) > vsmall)
    CALL ELIT(eccD,DBLE(theta(i)),FE,EE)
    EE = majr * EE - s(i)
    CALL DELIT(eccD,DBLE(theta(i)),EED)
    thetao(i) = theta(i)
    theta(i) = theta(i) - EE/(majr * EED)
    tol = ABS(thetao(i) - theta(i))
  ENDDO

ENDDO

theta(n/4+1) = 90.0

j = n/4
DO i = n/4+2,n/2+1
  theta(i) = 90.0 + (90.0 - theta(j) )
  j = j - 1
ENDDO

j = 2
DO i = n/2 + 2, 3*n/4
  theta(i) = 180.0 + theta(j)
  j = j + 1
ENDDO

theta(3*n/4+1) = 270.0

j = n/4
DO i = 3*n/4+2, n
  theta(i) = 270.0 + (90.0 - theta(j) )
  j = j - 1 
ENDDO

theta = theta * degtorad

DO i=1,n
  
  costheta = COS(theta(i))
  sintheta = SIN(theta(i))
  x(i) = (majr*minr/SQRT((minr*costheta)**2+(majr*sintheta)**2))*costheta
  y(i) = (majr*minr/SQRT((minr*costheta)**2+(majr*sintheta)**2))*sintheta

ENDDO

END SUBROUTINE fd_create_ellipse_surf

SUBROUTINE ELIT(HK,PHI,FE,EE)

!       ==================================================
!       Purpose: Compute complete and incomplete elliptic
!                integrals F(k,phi) and E(k,phi)
!       Input  : HK  --- Modulus k ( 0 < k < 1 )
!                Phi --- Argument ( in degrees )
!       Output : FE  --- F(k,phi)
!                EE  --- E(k,phi)
!       ==================================================

IMPLICIT DOUBLE PRECISION (A-H,O-Z)
G=0.0D0
PI=3.14159265358979D0
A0=1.0D0
B0=DSQRT(1.0D0-HK*HK)
D0=(PI/180.0D0)*PHI
R=HK*HK
IF (HK.EQ.1.0D0.AND.PHI.EQ.90.0D0) THEN
  FE=1.0D+300
  EE=1.0D0
ELSE IF (HK.EQ.1.0D0) THEN
 FE=DLOG((1.0D0+DSIN(D0))/DCOS(D0))
 EE=DSIN(D0)
ELSE
  FAC=1.0D0
  DO N=1,40
    A=(A0+B0)/2.0D0
    B=DSQRT(A0*B0)
    C=(A0-B0)/2.0D0
    FAC=2.0D0*FAC
    R=R+FAC*C*C
    IF (PHI.NE.90.0D0) THEN
      D=D0+DATAN((B0/A0)*DTAN(D0))
      G=G+C*DSIN(D)
      D0=D+PI*INT(D/PI+.5D0)
    ENDIF
    A0=A
    B0=B
    IF (C.LT.1.0D-7) EXIT
  ENDDO
  CK=PI/(2.0D0*A)
  CE=PI*(2.0D0-R)/(4.0D0*A)
  IF (PHI.EQ.90.0D0) THEN
    FE=CK
    EE=CE
  ELSE
    FE=D/(FAC*A)
    EE=FE*CE/CK+G
  ENDIF
ENDIF

END SUBROUTINE ELIT

SUBROUTINE DELIT(HK,PHI,EED)

!--First differential of the incomplete elliptic

IMPLICIT DOUBLE PRECISION (A-H,O-Z)
PI=3.14159265358979D0
DEG = 180.D0

EED = SQRT(1.D0 - HK**2*DSIN(PHI*PI/DEG)**2 )

END SUBROUTINE DELIT

SUBROUTINE fd_calculate_stats(nsph,collect_stat)

!--Calculates wake statistics
USE parameters,   ONLY : do_collect_stat,end_collect_stat,ubar_unit,vbar_unit,u2bar_unit,v2bar_unit,uvbar_unit
USE precision,    ONLY : r_single
USE shared_data,  ONLY : objradius,objcentx,objcenty,u,v,itst,li,ulid,x,y,xc,yc
USE real_parameters,  ONLY : four,two,zero

IMPLICIT NONE
INTEGER,INTENT(IN) :: nsph,collect_stat 
INTEGER,SAVE    :: n_entered = 0,ncollect=0
INTEGER,SAVE    :: ijloc(2,-20:20,5),ijinterpx(2,-20:20,5),ijinterpy(2,-20:20,5)
REAL(KIND = r_single),SAVE  :: yloc_array(-20:20),xloc_array(5)
REAL(KIND = r_single),SAVE,ALLOCATABLE :: uvtime(:,:,:,:)
REAL(KIND = r_single),SAVE  :: uvbar(-20:20,5),ubar(-20:20,5),vbar(-20:20,5),u2bar(-20:20,5),v2bar(-20:20,5)
INTEGER         :: i,j,k,n,nn
REAL(KIND = r_single)       :: yloc,del,dx,dy

n_entered = n_entered + 1
IF(n_entered == 1)THEN
  ALLOCATE(uvtime(2,-20:20,5,itst))
  uvtime = zero
  DO i = -20,20
    yloc_array(i) = objcenty(nsph) + REAL(i,r_single)*four*objradius(nsph)/REAL(20,r_single) 
  ENDDO
  xloc_array(1) = objcentx(nsph) + 1.2*two*objradius(nsph)
  xloc_array(2) = objcentx(nsph) + 1.5*two*objradius(nsph)
  xloc_array(3) = objcentx(nsph) + 2.0*two*objradius(nsph)
  xloc_array(4) = objcentx(nsph) + 2.5*two*objradius(nsph)
  xloc_array(5) = objcentx(nsph) + 3.0*two*objradius(nsph)
  DO i = 1,5
    DO j = -20,20
      CALL find_single_point(xloc_array(i),yloc_array(j),ijloc(1,j,i),ijloc(2,j,i))
    ENDDO
  ENDDO

  DO i=1,5
    DO j = -20,20
     IF(xloc_array(i) <= xc(ijloc(1,j,i)))THEN
       ijinterpx(1,j,i) = ijloc(1,j,i) - 1
       ijinterpx(2,j,i) = ijloc(1,j,i)
      ELSE
        ijinterpx(1,j,i) = ijloc(1,j,i)
        ijinterpx(2,j,i) = ijloc(1,j,i) + 1
      ENDIF
      IF(yloc_array(j) <= yc(ijloc(2,j,i)))THEN
        ijinterpy(1,j,i) = ijloc(2,j,i) - 1
        ijinterpy(2,j,i) = ijloc(2,j,i) 
      ELSE
        ijinterpy(1,j,i) = ijloc(2,j,i) 
        ijinterpy(2,j,i) = ijloc(2,j,i) + 1
      ENDIF
    ENDDO
  ENDDO

ENDIF

IF(collect_stat == do_collect_stat)THEN
  ncollect = ncollect + 1
  DO nn = 1,5
    DO n = -20,20
      DO i = ijinterpx(1,n,nn),ijinterpx(2,n,nn)
        dx = x(i) - x(i-1)
        DO j = ijinterpy(1,n,nn),ijinterpy(2,n,nn)
          dy = y(j) - y(j-1)
          del = fd_calc_delta2(xc(i),yc(j),xloc_array(nn),yloc_array(n),dx, dy)*dx*dy
          uvtime(1,n,nn,ncollect) = uvtime(1,n,nn,ncollect) + u(li(i)+j)*del
          uvtime(2,n,nn,ncollect) = uvtime(2,n,nn,ncollect) + v(li(i)+j)*del
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ELSEIF(collect_stat == end_collect_stat)THEN
  OPEN(UNIT=ubar_unit,FILE='ubar.dat',STATUS = 'UNKNOWN',ACCESS = 'APPEND')
  OPEN(UNIT=vbar_unit,FILE='vbar.dat',STATUS = 'UNKNOWN',ACCESS = 'APPEND')
  OPEN(UNIT=uvbar_unit,FILE='uvbar.dat',STATUS = 'UNKNOWN',ACCESS = 'APPEND')
  OPEN(UNIT=u2bar_unit,FILE='u2bar.dat',STATUS = 'UNKNOWN',ACCESS = 'APPEND')
  OPEN(UNIT=v2bar_unit,FILE='v2bar.dat',STATUS = 'UNKNOWN',ACCESS = 'APPEND')
  uvbar = zero
  ubar = zero;v2bar=zero
  vbar = zero;u2bar=zero
  DO i = 1,5
    DO j = -20,20 
      DO k = 1,ncollect
        ubar(j,i) = ubar(j,i) + uvtime(1,j,i,k)
        vbar(j,i) = vbar(j,i) + uvtime(2,j,i,k)
      ENDDO
    ENDDO
  ENDDO
  ubar = ubar/ncollect
  vbar = vbar/ncollect
  DO i = 1,5
    DO j = -20,20 
      DO k = 1,ncollect
        uvbar(j,i) = uvbar(j,i) + (uvtime(1,j,i,k)-ubar(j,i))*(uvtime(2,j,i,k)-vbar(j,i))
        u2bar(j,i) = u2bar(j,i) + uvtime(1,j,i,k)**2
        v2bar(j,i) = v2bar(j,i) + uvtime(2,j,i,k)**2
      ENDDO
    ENDDO
  ENDDO
  uvbar = uvbar/ncollect
  DO i = 1,5
    DO j = -20,20 
      u2bar(j,i) = u2bar(j,i)/ncollect - ubar(j,i)**2
      v2bar(j,i) = v2bar(j,i)/ncollect - vbar(j,i)**2 
    ENDDO
  ENDDO
  DO j = -20,20
    yloc = (yloc_array(j)-objcenty(nsph))/two/objradius(nsph)
    WRITE(ubar_unit,'(6(E12.5,1X))') yloc,ubar(j,1)/ulid,ubar(j,2)/ulid,ubar(j,3)/ulid,ubar(j,4)/ulid,ubar(j,5)/ulid
    WRITE(vbar_unit,'(6(E12.5,1X))') yloc,vbar(j,1)/ulid,vbar(j,2)/ulid,vbar(j,3)/ulid,vbar(j,4)/ulid,vbar(j,5)/ulid
    WRITE(uvbar_unit,'(6(E12.5,1X))') yloc,uvbar(j,1)/ulid**2,uvbar(j,2)/ulid**2,uvbar(j,3)/ulid**2,uvbar(j,4)/ &
                                                                                                ulid**2,uvbar(j,5)/ulid**2
    WRITE(u2bar_unit,'(6(E12.5,1X))') yloc,u2bar(j,1)/ulid**2,u2bar(j,2)/ulid**2,u2bar(j,3)/ulid**2,u2bar(j,4)/ &
                                                                                                ulid**2,u2bar(j,5)/ulid**2
    WRITE(v2bar_unit,'(6(E12.5,1X))') yloc,v2bar(j,1)/ulid**2,v2bar(j,2)/ulid**2,v2bar(j,3)/ulid**2,v2bar(j,4)/ &
                                                                                                ulid**2,v2bar(j,5)/ulid**2 
  ENDDO

  CLOSE(ubar_unit)
  CLOSE(vbar_unit)
  CLOSE(uvbar_unit)
  CLOSE(u2bar_unit)
  CLOSE(v2bar_unit)

ENDIF

END SUBROUTINE fd_calculate_stats

FUNCTION fd_calc_delta2(x, y , objx, objy, hx, hy)

!--A three point delta function---Only to be used for simple linear interpolation
USE precision,      ONLY : r_single
USE real_parameters,ONLY : one,zero,eight,three,four,seven,twelve,five,two,half,six

IMPLICIT NONE

REAL(KIND = r_single) :: x,y,objx,objy,hx,hy
REAL(KIND = r_single) :: rx,ry
REAL(KIND = r_single) :: fd_calc_delta2
REAL(KIND = r_single) :: deltaX, deltaY

rx = (x - objx)/hx
ry = (y - objy)/hy

  
  IF (ABS(rx) <= one) THEN
    deltax = one - ABS(rx)
  ELSE IF(ABS(rx) > one) THEN
    deltax = zero
  END IF
  
  IF (ABS(ry) <= one) THEN
    deltay = one - ABS(ry)
  ELSE IF(ABS(ry) > one) THEN
    deltay = zero
  END IF

fd_calc_delta2 = deltaX/hx*deltaY/hy

END FUNCTION fd_calc_delta2

SUBROUTINE fd_calc_part_collision(Fpq,Fpw)

!--A naive collision strategy, Only for O(10^3) spherical particles  

USE shared_data,      ONLY : nsphere,objcentx,objcenty,gravx,gravy,x,y,objradius,&
                             xPeriodic,yPeriodic,LDomainx,LDomainy

USE real_parameters,  ONLY : zero,two,five,ten,three,six,one
USE precision,        ONLY : r_single

IMPLICIT NONE
REAL(KIND = r_single),INTENT(INOUT) :: Fpq(2,nsphere),Fpw(2,nsphere)
INTEGER :: q,p
LOGICAL,SAVE :: Entered = .FALSE.
REAL(KIND = r_single),SAVE :: MaxRadius,miny,maxy,minx,maxx,g

IF(.NOT.Entered)THEN
  MaxRadius = MAXVAL(objradius(1:nsphere),1)
  minx = MINVAL(x,1)
  maxx = MAXVAL(x,1)
  miny = MINVAL(y,1)
  maxy = MAXVAL(y,1)
  g = SQRT(gravx**2+gravy**2)
  Entered = .TRUE.
ENDIF

Fpq = zero
Fpw = zero

!--Particle-particle collision
DO p = 1,nsphere
  DO q = 1,nsphere
    IF(p /= q)THEN
      Fpq(:,p) = Fpq(:,p) + fd_calc_binary_colforce(objcentx(p),objcenty(p),objcentx(q),objcenty(q),p,q)
    ENDIF
  ENDDO
ENDDO

!--Periodic Colisions
IF(xPeriodic == 1)THEN
  DO p = 1,nsphere
    IF(objcentx(p) > maxx - six * MaxRadius)THEN
      DO q = 1,nsphere
        IF(objcentx(q) < minx + six * MaxRadius)THEN
          Fpq(:,p) = Fpq(:,p) + fd_calc_binary_colforce(objcentx(p),objcenty(p),objcentx(q)+ &
                                                        LDomainx,objcenty(q),p,q)
        ENDIF
      ENDDO
    ELSEIF(objcentx(p) < minx + six * MaxRadius)THEN
      DO q = 1,nsphere
        IF(objcentx(q) > maxx - six * MaxRadius)THEN
          Fpq(:,p) = Fpq(:,p) + fd_calc_binary_colforce(objcentx(p),objcenty(p),objcentx(q)- &
                                                        LDomainx,objcenty(q),p,q)
        ENDIF
      ENDDO
    ENDIF
  ENDDO
ENDIF

!--Periodic Colisions
IF(yPeriodic == 1)THEN
  DO p = 1,nsphere
    IF(objcenty(p) > maxy - six * MaxRadius)THEN
      DO q = 1,nsphere
        IF(objcenty(q) < miny + six * MaxRadius)THEN
          Fpq(:,p) = Fpq(:,p) + fd_calc_binary_colforce(objcentx(p),objcenty(p),objcentx(q),objcenty(q)+ &
                                                        LDomainy,p,q)
        ENDIF
      ENDDO
    ELSEIF(objcenty(p) < miny + six * MaxRadius)THEN
      DO q = 1,nsphere
        IF(objcenty(q) > maxy - six * MaxRadius)THEN
          Fpq(:,p) = Fpq(:,p) + fd_calc_binary_colforce(objcentx(p),objcenty(p),objcentx(q),objcenty(q)- &
                                                        LDomainy,p,q)
        ENDIF
      ENDDO
    ENDIF
  ENDDO
ENDIF

!--Particle-wall collision, (Mirror image)
DO p = 1,nsphere 

  IF(xPeriodic == 0)THEN
    !--e wall
    Fpw(:,p) = fd_calc_binary_colforce(objcentx(p),objcenty(p),maxx+objradius(p),objcenty(p),p,p)
    !--w wall
    Fpw(:,p) = Fpw(:,p) + fd_calc_binary_colforce(objcentx(p),objcenty(p),minx-objradius(p),objcenty(p),p,p)
  ENDIF
  IF(yPeriodic == 0)THEN
    !--n wall
    Fpw(:,p) = Fpw(:,p) + fd_calc_binary_colforce(objcentx(p),objcenty(p),objcentx(p),maxy+objradius(p),p,p)
    !--s wall
    Fpw(:,p) = Fpw(:,p) + fd_calc_binary_colforce(objcentx(p),objcenty(p),objcentx(p),miny-objradius(p),p,p)
  ENDIF

ENDDO

END SUBROUTINE fd_calc_part_collision

FUNCTION fd_calc_binary_colforce(xp,yp,xq,yq,p,q)

!--Suggested by Glowinsky 1995, Only forces in normal direction  

USE shared_data,    ONLY : dxmin,densitp,objvol,objradius,dem_kn,dem_kt,dem_gamman,dem_gammat,dem_mut,dem_dt
USE precision,      ONLY : r_single
USE real_parameters,ONLY : distfactor,zero,half,one

IMPLICIT NONE 

REAL(KIND = r_single) :: fd_calc_binary_colforce(2)
REAL(KIND = r_single),INTENT(IN)  :: xp,yp,xq,yq
INTEGER              ,INTENT(IN)  :: p,q
REAL(KIND = r_single)             :: dummy,dpq,delta,npq(2) !-unit vect and distance between two centres

CALL fd_calc_unitvector(xq, yq, xp, yp, npq(1), npq(2), dpq)

delta = (objradius(p) + objradius(q) + distfactor*dxmin) - dpq

dummy = MAX(zero, delta)
 
fd_calc_binary_colforce = dem_kn*npq*dummy

END FUNCTION fd_calc_binary_colforce

SUBROUTINE fd_find_rigidforce_contrib

USE precision,      ONLY : r_single
USE shared_data,    ONLY : nobjcells,nsphere,xPeriodic,yPeriodic,&
                           nim,njm,li,objcentx,objcenty,rigidforce_contrib,xc,yc,&
                           objpoint_interpx,objpoint_interpy,LDomainy,LDomainx

IMPLICIT NONE

INTEGER              :: nn,n,mm,i,ii,j,jj,ij,nim2,njm2
REAL(KIND = r_single):: distn,distm,ypn,xpn,ypm,xpm

nim2 = nim - 1
njm2 = njm - 1

rigidforce_contrib = 0

DO nn = 1,nsphere
  DO n = 1,nobjcells(nn)
    DO i = objpoint_interpx(1,n,nn),objpoint_interpx(2,n,nn)
      ii = i + xPeriodic*(MAX(2-i,0)*nim2 + MIN(nim-i,0)*nim2)
      DO j = objpoint_interpy(1,n,nn),objpoint_interpy(2,n,nn)
        jj = j + yPeriodic*(MAX(2-j,0)*njm2 + MIN(njm-j,0)*njm2)
        ij = li(ii)+jj
        IF(rigidforce_contrib(ij) == 0)THEN
            rigidforce_contrib(ij) = nn
        ELSEIF(rigidforce_contrib(ij) /= nn)THEN
            ypn= objcenty(nn) + yPeriodic*(MAX(2-j,0)*LDomainy + MIN(njm-j,0)*LDomainy)
            xpn= objcentx(nn) + xPeriodic*(MAX(2-i,0)*LDomainx + MIN(nim-i,0)*LDomainx) 
            distn = SQRT( (xpn - xc(ii))**2 + (ypn - yc(jj))**2)
            mm = rigidforce_contrib(ij) 
            xpm= objcentx(mm) + xPeriodic*(MAX(2-i,0)*LDomainx + MIN(nim-i,0)*LDomainx)
            ypm= objcenty(mm) + yPeriodic*(MAX(2-j,0)*LDomainy + MIN(njm-j,0)*LDomainy)   
            distm = SQRT( (xpm - xc(ii))**2 + (ypm - yc(jj))**2)
            IF(distn < distm) rigidforce_contrib(ij) = nn
        ENDIF
       ENDDO
    ENDDO
   ENDDO
ENDDO

END SUBROUTINE fd_find_rigidforce_contrib

SUBROUTINE fd_create_interpmol(n)
  
USE real_parameters,ONLY : one,two,small,half,zero
USE precision,      ONLY : r_single
USE shared_data,    ONLY : surfpointx,surfpointy,objcentx,objcenty,&
                               surfds,surfnx,surfny,surfcentrex,surfcentrey,&
                               xc,yc,surfpoint_cvx,surfpoint_cvy,&
                               surfpoint_interp,surfinterpx,surfinterpy,xc,yc,x,y,nsphere,&
                               nsurfpoints,objradius,nsphere

IMPLICIT NONE
INTEGER,INTENT(IN) :: n
INTEGER            :: ip1
REAL(KIND = r_single) :: vecx,vecy,coeft,dxmeanl,xp,yp
LOGICAL               :: failed
INTEGER               :: i,ip,jp,is,ie,js,je,ii,jj


DO i=1,nsurfpoints(n)
  IF(i < nsurfpoints(n))THEN
    ip1 = i + 1
  ELSE
    ip1 = 1
  ENDIF

  !--Calc centre points
  surfcentrex(i,n) = (surfpointx(i,n) + surfpointx(ip1,n))/two
  surfcentrey(i,n) = (surfpointy(i,n) + surfpointy(ip1,n))/two

  !--Locate the centre point on the Eul grid
  CALL find_single_point(surfcentrex(i,n), surfcentrey(i,n), surfpoint_cvx(1,i,n),surfpoint_cvy(1,i,n))
  !--Ignor if outside
  !IF(bndFlag > 0)CYCLE

  vecx = surfpointx(i,n) - surfpointx(ip1,n)
  vecy = surfpointy(i,n) - surfpointy(ip1,n)

  surfds(i,n) = SQRT(vecx**2 + vecy**2)
  IF(surfds(i,n) < small) THEN
    PRINT *,'Error - Coordinates co-exist - undefined normal!'
    STOP
  ENDIF

  !--Surface element normals (outward)
  surfnx(i,n) = -vecy/surfds(i,n)
  surfny(i,n) = vecx/surfds(i,n)

  !--Ensure the calculation is corrent
  IF((surfcentrex(i,n) - objcentx(n))*surfnx(i,n) + &
    (surfcentrey(i,n) - objcenty(n))*surfny(i,n) < zero)THEN

    PRINT *,'Error - outward normal not calculated!'
    STOP
  ENDIF

ENDDO

DO i=1,nsurfpoints(n)

  coeft = one
  dxmeanl = ( (x(surfpoint_cvx(1,i,n)) - x(surfpoint_cvx(1,i,n)-1)) + &
              (y(surfpoint_cvy(1,i,n)) - y(surfpoint_cvy(1,i,n)-1)) )/two
  DO
    failed = .FALSE.
    !--first interpolation point
    xp = surfcentrex(i,n) + coeft * surfnx(i,n) * dxmeanl
    yp = surfcentrey(i,n) + coeft * surfny(i,n) * dxmeanl

    !--Locate the point
    CALL find_single_point(xp,yp,ip,jp)
    !IF(bndFlag > 0)THEN
           
    !--indentify the neighbours
    IF(xp > xc(ip))THEN
      is = ip
      ie = ip + 1
    ELSE
      is = ip - 1
      ie = ip
    ENDIF

    IF(yp > yc(jp))THEN
      js = jp
      je = jp + 1
    ELSE
      js = jp - 1
      je = jp
    ENDIF

    !--Ensure all neighbours are outside of the object
    DO ii = is,ie
      DO jj = js,je
        IF( ( (xc(ii) - surfcentrex(i,n))*surfnx(i,n) + &
          (yc(jj) - surfcentrey(i,n))*surfny(i,n) ) < zero )THEN

          failed = .TRUE.
          coeft = coeft + half * dxmeanl
        ENDIF

      ENDDO
    ENDDO
    IF(.NOT. failed)EXIT
  ENDDO

  surfpoint_interp(1:4,1,i,n) = (/is,ie,js,je/)
  surfpoint_cvx(2,i,n) = ip
  surfpoint_cvy(2,i,n) = jp
  surfinterpx(1,i,n) = xp
  surfinterpy(1,i,n) = yp

  !--second point has 2 times the distance of the first point to the surface centre
  !--Assuming a convex surface no need to check for the neighbours
  surfinterpx(2,i,n) = two * surfinterpx(1,i,n) - surfcentrex(i,n)
  surfinterpy(2,i,n) = two * surfinterpy(1,i,n) - surfcentrey(i,n)

  CALL find_single_point(surfinterpx(2,i,n),surfinterpy(2,i,n),ip,jp)

  IF(xp > xc(ip))THEN
    is = ip
    ie = ip + 1
  ELSE
    is = ip - 1
    ie = ip
  ENDIF

  IF(yp > yc(jp))THEN
    js = jp
    je = jp + 1
  ELSE
    js = jp - 1
    je = jp
  ENDIF

  surfpoint_interp(1:4,2,i,n) = (/is,ie,js,je/)
  surfpoint_cvx(3,i,n) = ip
  surfpoint_cvy(3,i,n) = jp

ENDDO

END SUBROUTINE fd_create_interpmol

SUBROUTINE fd_create_interpmol_nus(n)

USE precision,          ONLY : r_single
USE real_parameters,    ONLY : pi,one,two,small,half,zero
USE shared_data,        ONLY : nnusseltpoints,nusseltpointx,nusseltpointy,objcentx,objradius,objcenty,&
                               nusseltcentx,nusseltcenty,nusseltpoint_cvx,nusseltpoint_cvy,nusseltds,&
                               nusseltnx,nusseltny,xc,yc,nusseltpoint_interp,nusseltinterpx,nusseltinterpy,&
                               objbradius,x,y


IMPLICIT NONE

INTEGER,INTENT(IN)          :: n
INTEGER                     :: i,ip1,ip,jp,is,ie,js,je,ii,jj
REAL(KIND = r_single)       :: dtheta,theta,costheta,sintheta,vecx,vecy,dxmeanl,coeft,xp,yp
LOGICAL                     :: failed

!--Used for single tube
  dtheta = pi / REAL(nnusseltpoints(n)-1)
  !--Ued for tube bank
  !dtheta = two*pi / REAL(nnusseltpoints(n)-1)
  DO i=1,nnusseltpoints(n)
  
    theta = REAL(i-1)*dtheta
    costheta = COS(pi-theta)
    sintheta = SIN(pi-theta)
    nusseltpointx(i,n) = objcentx(n)+(objradius(n)*objbradius(n)/ &
                    SQRT((objradius(n)*costheta)**2+(objbradius(n)*sintheta)**2))*costheta
    nusseltpointy(i,n) = objcenty(n)+(objradius(n)*objbradius(n)/ &
                    SQRT((objradius(n)*costheta)**2+(objbradius(n)*sintheta)**2))*sintheta
!    nusseltpointx(i,n) = objcentx(n)+objradius(n)*COS(theta)
!    nusseltpointy(i,n) = objcenty(n)+objradius(n)*SIN(theta)
  ENDDO

  DO i=1,nnusseltpoints(n)-1
      
    ip1 = i + 1

    !--Calc centre points
    nusseltcentx(i,n) = (nusseltpointx(i,n) + nusseltpointx(ip1,n))/two
    nusseltcenty(i,n) = (nusseltpointy(i,n) + nusseltpointy(ip1,n))/two

    !--Locate the centre point on the Eul grid
    CALL find_single_point(nusseltcentx(i,n), nusseltcenty(i,n), nusseltpoint_cvx(1,i,n),nusseltpoint_cvy(1,i,n))
        
    !--For single tube     
    vecx = nusseltpointx(ip1,n) - nusseltpointx(i,n) 
    vecy = nusseltpointy(ip1,n) - nusseltpointy(i,n) 
    !--For bank
    !vecx = nusseltpointx(i,n) - nusseltpointx(ip1,n) 
    !vecy = nusseltpointy(i,n) - nusseltpointy(ip1,n) 
             
    nusseltds(i,n) = SQRT(vecx**2 + vecy**2)
    IF(nusseltds(i,n) < small) THEN
      PRINT *,'Error - Coordinates co-exists - undefined normal!'
      STOP
    ENDIF
      
    !--Surface element normals (outward)
    nusseltnx(i,n) = -vecy/nusseltds(i,n)
    nusseltny(i,n) = vecx/nusseltds(i,n)

    !--Ensure the calculation is corrent
    IF((nusseltcentx(i,n) - objcentx(n))*nusseltnx(i,n) + &
       (nusseltcenty(i,n) - objcenty(n))*nusseltny(i,n) < zero)THEN

       PRINT *,'Error - outward normal not calculated!'
       STOP
    ENDIF

  ENDDO

  DO i=1,nnusseltpoints(n)-1
      
    coeft = one
    dxmeanl = ( (x(nusseltpoint_cvx(1,i,n)) - x(nusseltpoint_cvx(1,i,n)-1)) + & 
                (y(nusseltpoint_cvy(1,i,n)) - y(nusseltpoint_cvy(1,i,n)-1))   )/two
    DO          
      failed = .FALSE. 
      !--first interpolation point
      xp = nusseltcentx(i,n) + coeft * nusseltnx(i,n) * dxmeanl 
      yp = nusseltcenty(i,n) + coeft * nusseltny(i,n) * dxmeanl 
      
      !--Locate the point
      CALL find_single_point(xp,yp,ip,jp)
        
      !--indentify the neighbours      
      IF(xp > xc(ip))THEN
        is = ip
        ie = ip + 1
      ELSE
        is = ip - 1
        ie = ip
      ENDIF

      IF(yp > yc(jp))THEN
        js = jp
        je = jp + 1
      ELSE
        js = jp - 1
        je = jp
      ENDIF

      !--Ensure all neighbours are outside of the object
      DO ii = is,ie
        DO jj = js,je
          IF( ( (xc(ii) - nusseltcentx(i,n))*nusseltnx(i,n) + &
                (yc(jj) - nusseltcenty(i,n))*nusseltny(i,n)    ) < zero )THEN

             failed = .TRUE.
             coeft = coeft + half * dxmeanl
          ENDIF
             
        ENDDO
      ENDDO
      IF(.NOT. failed)EXIT
    ENDDO 

    nusseltpoint_interp(1:4,1,i,n) = (/is,ie,js,je/)
    nusseltpoint_cvx(2,i,n) = ip
    nusseltpoint_cvy(2,i,n) = jp
    nusseltinterpx(1,i,n) = xp
    nusseltinterpy(1,i,n) = yp

    !--second point has 2 times the distance of the first point to the surface centre
    !--Assuming a convex surface no need to check for the neighbours
    nusseltinterpx(2,i,n) = two * nusseltinterpx(1,i,n) - nusseltcentx(i,n) 
    nusseltinterpy(2,i,n) = two * nusseltinterpy(1,i,n) - nusseltcenty(i,n)
    
    CALL find_single_point(nusseltinterpx(2,i,n),nusseltinterpy(2,i,n),ip,jp)
    
    IF(xp > xc(ip))THEN
      is = ip
      ie = ip + 1
    ELSE
      is = ip - 1
      ie = ip
    ENDIF

    IF(yp > yc(jp))THEN
      js = jp
      je = jp + 1
    ELSE
      js = jp - 1
      je = jp
    ENDIF

      nusseltpoint_interp(1:4,2,i,n) = (/is,ie,js,je/)
      nusseltpoint_cvx(3,i,n) = ip
      nusseltpoint_cvy(3,i,n) = jp 
       
  ENDDO

END SUBROUTINE fd_create_interpmol_nus


SUBROUTINE fd_update_fieldvel

USE shared_data,    ONLY : objru,objrv,obju,objv,xPeriodic,yPeriodic,&
                           objpoint_cvx, objpoint_cvy,x,y,objcellx,objcelly,&
                           xc,yc,objcentu,objcentv,objcentom,objcentx,objcenty,&
                           dt,li,objcellvol,nobjcells,u,v,densitp,nij,fd_urf,objvol,&
                           nim,njm,fx,fy,fdsub,fdsvb,nj,objtp,objt,objrt,objq,fdst,t,&
                           nsphere,fdsvc,objpoint_interpx,objpoint_interpy,&
                           zobjcenty,zobjcentx,zobjcelly,zobjcellx,LDomainx,LDomainy,rigidforce_contrib,&
                           objcell_bndFlag,u,v,f1,f2,fx,fy,ft1,ft2,celcp,den,r
USE precision,      ONLY : r_single
USE real_parameters,ONLY : zero,one,half
IMPLICIT NONE

!--locals
INTEGER               :: i,j,n,nn,ii,jj,nim2,njm2,ij,ije,ijn,xend,yend
REAL(KIND = r_single) :: dx,dy,del,xp,yp,fxe,fxp,uel,fyn,fyp,vnl,s

nim2 = nim - 1
njm2 = njm - 1

obju = zero
objv = zero

DO nn = 1,nsphere
  DO n = 1,nobjcells(nn)
    IF(objcell_bndFlag(n,nn) == 0)THEN
      DO i = objpoint_interpx(1,n,nn),objpoint_interpx(2,n,nn)
        ii = i + xPeriodic*(MAX(2-i,0)*nim2 + MIN(nim-i,0)*nim2)
        xp = objcellx(n,nn) + xPeriodic*(MAX(2-i,0)*LDomainx + MIN(nim-i,0)*LDomainx)
        dx = x(ii) - x(ii-1)
        DO j = objpoint_interpy(1,n,nn),objpoint_interpy(2,n,nn)
          jj = j + yPeriodic*(MAX(2-j,0)*njm2 + MIN(njm-j,0)*njm2) 
          ij = li(ii)+jj 
          !IF(rigidForce_contrib(ij) == nn)THEN
            yp = objcelly(n,nn) + yPeriodic*(MAX(2-j,0)*LDomainy + MIN(njm-j,0)*LDomainy)
            dy = y(jj) - y(jj-1)
            del = fd_calc_delta(xc(ii),yc(jj),xp,yp,dx, dy)*dx*dy
            obju(n,nn) = obju(n,nn) + u(ij)*del
            objv(n,nn) = objv(n,nn) + v(ij)*del
          !ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDDO
ENDDO

DO nn = 1,nsphere
  DO n = 1,nobjcells(nn)
      objru(n,nn) = objcentu(nn) - objcentom(nn)*(objcelly(n,nn) - objcenty(nn)+(zobjcelly(n,nn) - zobjcenty(nn))*LDomainy)
      objrv(n,nn) = objcentv(nn) + objcentom(nn)*(objcellx(n,nn) - objcentx(nn)+(zobjcellx(n,nn) - zobjcentx(nn))*LDomainx)
  ENDDO
ENDDO

DO nn = 1,nsphere
  DO n = 1,nobjcells(nn)
    IF(objcell_bndFlag(n,nn) == 0)THEN
      DO i = objpoint_interpx(1,n,nn),objpoint_interpx(2,n,nn)
        ii = i + xPeriodic*(MAX(2-i,0)*nim2 + MIN(nim-i,0)*nim2)
        xp = objcellx(n,nn) + xPeriodic*(MAX(2-i,0)*LDomainx + MIN(nim-i,0)*LDomainx)
        dx = x(ii) - x(ii-1)
        DO j =objpoint_interpy(1,n,nn),objpoint_interpy(2,n,nn)
          jj = j + yPeriodic*(MAX(2-j,0)*njm2 + MIN(njm-j,0)*njm2)
          yp = objcelly(n,nn) + yPeriodic*(MAX(2-j,0)*LDomainy + MIN(njm-j,0)*LDomainy)
          dy = y(jj) - y(jj-1)
          del = fd_calc_delta(xc(ii),yc(jj),xp,yp,dx,dy)*objcellvol(n,nn)
          ij = li(ii)+jj
          !IF(rigidForce_contrib(ij) == nn)THEN
            u(ij) = u(ij) + (objru(n,nn)-obju(n,nn))*del
            v(ij) = v(ij) + (objrv(n,nn)-objv(n,nn))*del
          !ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDDO
ENDDO 

xend = nim+xPeriodic-1 !--Apply periodic conditions
DO i=2,xend

  fxe=fx(i)
  fxp=one-fxe

  DO j=2,njm
    ij=li(i)+j
    ije=ij+nj-i/nim*((i-1)*nj)
    
    uel=u(ije)*fxe+u(ij)*fxp
    
    s=(y(j)-y(j-1))*(r(j)+r(j-1))*half
    f1(ij)=(den(ije)*fxe+den(ij)*fxp)*s*uel
    ft1(ij)=(celcp(ije)*fxe+celcp(ij)*fxp)*s*uel
    
  ENDDO
ENDDO

yend = njm+yPeriodic-1 !--Apply periodic conditions
DO j=2,yend

  fyn=fy(j)
  fyp=one-fyn

  DO i=2,nim
    ij=li(i)+j
    ijn=ij+1-j/njm*(j-1)
    
    vnl=v(ijn)*fyn+v(ij)*fyp
    
    s=(x(i)-x(i-1))*r(j)
    f2(ij)=(den(ijn)*fyn+den(ij)*fyp)*s*vnl
    ft2(ij)=(celcp(ijn)*fyn+celcp(ij)*fyp)*s*vnl

  ENDDO
ENDDO

END SUBROUTINE fd_update_fieldvel
ENDMODULE modfd_create_geom
