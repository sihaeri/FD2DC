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

MODULE modlmp_particle

PRIVATE
PUBLIC :: lmp_calc_pos,lmp_setup

CONTAINS
!====================================================================================
SUBROUTINE lmp_setup


USE lammps
USE cula_sparse
USE real_parameters,    ONLY : zero,two,pi,half,quarter,vlarge,threequarter,one,fourthird,small,vsmall
USE parameters,         ONLY : subTimeStep,use_GPU_yes,use_GPU_no,USE_LAMMPS_YES,USE_LAMMPS_NO,&
                               LMP_COMMANDS,RANK_ONE,RANK_THREE,USE_GPU_YES
USE precision,          ONLY : r_single
USE shared_data,        ONLY : handle,use_lammps, lmp_ptr, lmp_kn, lmp_kt, lmp_gamman, lmp_gammat, lmp_mut, lmp_flag,&
                           lmp_command, lmp_dt, lmp_nstep, lmp_x, lmp_v, lmp_o, lmp_objRadius, lmp_objMass, lmp_fname,&
                           lmp_dppDummy, lmp_iDummy, nsphere, x, y, nim, njm, xperiodic, yperiodic, objradius, densitp, &
                           dxmean, use_gpu, dt, objvol
IMPLICIT NONE
!--Local Data
INTEGER  :: i,error,j 
REAL(KIND=r_single)  :: rdummy,massi,massj,mstar

ALLOCATE(lmp_objRadius(1:nsphere))
CALL lammps_open_no_mpi('lmp -log lammps.log',lmp_ptr )

!--Some fixed commands
DO i = 1,SIZE(LMP_COMMANDS)
  WRITE(lmp_command,'(A)')LMP_COMMANDS(i)
  CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))
ENDDO

WRITE(lmp_command,'(A)')'atom_modify map array'
CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))

!--Setup box and atoms
WRITE(lmp_command,'(A,1X,A)')'read_data',lmp_fname
CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))

!--Check everythin is fine
error = 0
CALL lammps_extract_global(lmp_dppDummy,lmp_ptr,'boxxlo')
IF(abs(lmp_dppDummy) -  abs(x(1)) > small)THEN
    WRITE(*,*)'Lammps data file not consistent: xlo'
    error = error + 1
ENDIF
WRITE(*,*)'xlo: ',lmp_dppDummy

CALL lammps_extract_global(lmp_dppDummy,lmp_ptr,'boxxhi')
IF(abs(lmp_dppDummy) -  abs(x(nim)) > small)THEN
    WRITE(*,*)'Lammps data file not consistent: xhi'
    error = error + 1
ENDIF
WRITE(*,*)'xhi: ',lmp_dppDummy

CALL lammps_extract_global(lmp_dppDummy,lmp_ptr,'boxylo')
IF(abs(lmp_dppDummy) -  abs(y(1)) > small)THEN
    WRITE(*,*)'Lammps data file not consistent: ylo'
    error = error + 1
ENDIF
WRITE(*,*)'ylo: ',lmp_dppDummy

CALL lammps_extract_global(lmp_dppDummy,lmp_ptr,'boxyhi')
IF(abs(lmp_dppDummy) -  abs(y(njm)) > small)THEN
    WRITE(*,*)'Lammps data file not consistent: yhi'
    error = error + 1
ENDIF
WRITE(*,*)'yhi: ',lmp_dppDummy

CALL lammps_extract_global(lmp_dppDummy,lmp_ptr,'boxzlo')
IF(abs(lmp_dppDummy) - abs(quarter*minval(objradius)) > small)THEN
    WRITE(*,*)'Lammps data file not consistent: zlo'
    error = error + 1
ENDIF
WRITE(*,*)'zlo: ',lmp_dppDummy

CALL lammps_extract_global(lmp_dppDummy,lmp_ptr,'boxzhi')
IF(abs(lmp_dppDummy) -  abs(quarter*maxval(objradius)) > small)THEN
    WRITE(*,*)'Lammps data file not consistent: zhi'
    error = error + 1
ENDIF
WRITE(*,*)'zhi: ',lmp_dppDummy

lmp_iDummy = lammps_get_natoms(lmp_ptr)
IF(lmp_iDummy /= nsphere)THEN
    WRITE(*,*)'Lammps data file not consistent: nsphere'
    error = error + 1
ENDIF
WRITE(*,*)'ntypes=',lammps_get_ntypes(lmp_ptr),', nspheres=',lmp_iDummy

IF(error /=0)THEN
  CALL lammps_close(lmp_ptr)
  IF(use_gpu == USE_GPU_YES)THEN
    IF(culaSparseDestroy(handle) == culaSparseNoError)THEN
      STOP
    ELSE
      WRITE(*,*)'Could not stop cula. Shutting down anyways ...'
      STOP
    ENDIF
  ELSE
    STOP
  ENDIF
ENDIF

!--Set neighbour list props
WRITE(lmp_command,'(A,1X,F12.4,1X,A)')'neighbor',half*minval(objradius),'bin'
CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))
   
WRITE(lmp_command,'(A)')'neigh_modify check no every 1 delay 0'
CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))

IF(xperiodic==1 .AND. yperiodic==1)THEN
  WRITE(lmp_command,'(A)')'change_box all boundary p p p'
  CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))

ELSEIF(xperiodic==1 .AND. yperiodic==0)THEN
  
  WRITE(lmp_command,'(A)')'change_box all boundary p f p'
  CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))
  
  WRITE(lmp_command,'(A,1X,5(F12.4,1X),I5,A,2(F12.4,1X))')'fix ygranWall all wall/gran', &
  lmp_kn,lmp_kt,lmp_gamman,lmp_gammat,lmp_mut,lmp_flag,' yplane ',y(1),y(njm)
  CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))

ELSEIF(xperiodic==0 .AND. yperiodic==1)THEN
  
  WRITE(lmp_command,'(A)')'change_box all boundary f p p'
  CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))
  
  WRITE(lmp_command,'(A,1X,5(F12.4,1X),I5,A,2(F12.4,1X))')'fix xgranWall all wall/gran', &
  lmp_kn,lmp_kt,lmp_gamman,lmp_gammat,lmp_mut,lmp_flag,' xplane ',x(1),x(nim)
  CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))

ELSE

  WRITE(lmp_command,'(A)')'change_box all boundary f f p'
  CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))
  
  WRITE(lmp_command,'(A,1X,5(F12.4,1X),I5,A,2(F12.4,1X))')'fix xgranWall all wall/gran', &
  lmp_kn,lmp_kt,lmp_gamman,lmp_gammat,lmp_mut,lmp_flag,' xplane ',x(1),x(nim)
  CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))
  
  WRITE(lmp_command,'(A,1X,5(F12.4,1X),I5,A,2(F12.4,1X))')'fix ygranWall all wall/gran', &
  lmp_kn,lmp_kt,lmp_gamman,lmp_gammat,lmp_mut,lmp_flag,' yplane ',y(1),y(njm)
  CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))

ENDIF
     
DO i=1,nsphere
  !--Calculate a radius
  lmp_ObjRadius(i) = objRadius(i)+(one+threequarter)*dxmean
  WRITE(*,*)'dxmean: ',dxmean,', LAMMPS Radius:', lmp_ObjRadius(i)
  WRITE(lmp_command,'(A,I5,1X,A,F24.16)')'set atom ', i, 'diameter ', two*lmp_ObjRadius(i)
  CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))
  
  WRITE(lmp_command,'(A,I5,1X,A,F24.16)')'set atom ', i, 'mass ', densitp(i)*objvol(i)
  CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))
ENDDO

CALL lammps_gather_atoms(lmp_ptr, 'rmass', RANK_ONE, lmp_objMass )
DO i=1,nsphere
  WRITE(*,'(A,I5,A,F12.4)')'Mass of particle ',i,' ,is', lmp_objMass(i)
ENDDO

!--Minimum contact time
rdummy = vsmall
DO i=1,nsphere
  massi = lmp_objMass(i)
  WRITE(*,'(A,I5,A,F12.4)')'Code calculated mass of particle ',i,' ,is', massi
  DO j=i,nsphere
    massj = lmp_objMass(j)
    mstar = (massi*massj)/(massi+massj)
    rdummy = MAXVAL((/rdummy,SQRT(lmp_kn/mstar - (lmp_gamman/two)**2)/))
  ENDDO
ENDDO

lmp_dt = pi/rdummy
WRITE(*,*)'Lammps time step before adjustment: ',lmp_dt,', fluid dt: ', dt

!--Fix dt to get exactly from t to t+dt
IF(lmp_dt > dt)THEN
  lmp_dt = dt / subTimeStep
ELSE
  lmp_dt = lmp_dt / subTimeStep
ENDIF
lmp_nstep = FLOOR(dt/lmp_dt);
lmp_dt = dt/REAL(lmp_nstep,r_single)
WRITE(*,*)'Lammps time step: ',lmp_dt,', fluid dt: ', dt, ' lammps steps per fluid step:', lmp_nstep

WRITE(lmp_command,'(A,1X,F24.16)')'timestep',lmp_dt 
CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))

!--Setup DEM 
WRITE(lmp_command,'(A,1X,5(F12.4,1X),I5)')'pair_style gran/hooke/history',lmp_kn,lmp_kt,&
                                          lmp_gamman,lmp_gammat,lmp_mut,lmp_flag 
CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))

WRITE(lmp_command,'(A)')'pair_coeff * *'
CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))

!--Setup integrators and enforce 2d simulation
WRITE(lmp_command, '(A)')'fix fix2D all enforce2d'
CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))

WRITE(lmp_command, '(A)')'fix integrateFix all nve/sphere'
CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))

END SUBROUTINE lmp_setup

SUBROUTINE lmp_calc_pos(itr,titr)

USE lammps
USE parameters,       ONLY : OUTER_ITR_DONE,RANK_THREE,RANK_ONE
USE precision,        ONLY : r_single
USE shared_data,      ONLY : nsphere,nobjcells,objcentmi,objcentxo,objcentyo,objcellx,objcelly,objradius,objcentom,&
                             objcentxo,objcentyo,dt,objcentu,objcentv,objcento,objpoint_interpx,objpoint_interpy,&
                             objcellvertx,objcellverty,objcentx,objcenty,objcentuo,objcentvo,calcsurfforce,&
                             objpoint_cvx,objpoint_cvy,nsurfpoints,surfpointx,surfpointy,dxmeanmoved,&
                             forcedmotion,x,y,xc,zobjcentx,zobjcentxo,zobjcenty,zobjcentyo,zobjcellx,zobjcelly,&
                             LDomainx,LDomainy,zobjcellvertx,zobjcellverty,nim,njm,zsurfpointx,zsurfpointy,&
                             lmp_ptr,lmp_x,lmp_v,lmp_o,lmp_img,lmp_command,lmp_nstep,lmp_objRadius,objcell_bndFlag

USE real_parameters,  ONLY : zero,three,two,half,one,four,five
USE modfd_create_geom,ONLY : fd_calc_ori, fd_track_single_point, fd_calc_forcedmotion, fd_update_forcedvel,&
                             fd_create_interpmol,fd_find_rigidforce_contrib

IMPLICIT NONE

INTEGER,INTENT(IN)    :: itr,titr

INTEGER               :: n,nn,ip,jp,error,moved,m
REAL(KIND = r_single) :: lindispx,lindispy,rmcx,rmcy,rmcyr,rmcxr,dxmoved(nsphere),&
                         xShift,yShift,objcentuk(nsphere),objcentvk(nsphere)

REAL(KIND = r_single),ALLOCATABLE,SAVE      :: objcellxo(:,:),objcellyo(:,:),fd2lamOmConv(:)
INTEGER,ALLOCATABLE,SAVE                    :: zobjcellxo(:,:),zobjcellyo(:,:)


IF(titr == 1 )THEN !.AND. itr ==1 
ALLOCATE( objcellxo(MAXVAL(nobjcells,1),nsphere),objcellyo(MAXVAL(nobjcells,1),nsphere),&
          zobjcellxo(MAXVAL(nobjcells,1),nsphere),zobjcellyo(MAXVAL(nobjcells,1),nsphere),&
          fd2lamOmConv(nsphere))

fd2lamOmConv = five/four * objRadius / lmp_objRadius
!fd2lamVConv  = three/four * objRadius**2 / lmp_objRadius**3 

WRITE(*,*)'Omega Conversion factor: ',fd2lamOmConv
!WRITE(*,*)'Vel Conversion factor: ',fd2lamVConv
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
  
  IF(.NOT. ALLOCATED(lmp_x))THEN;ALLOCATE(lmp_x(RANK_THREE*nsphere));ENDIF
  IF(.NOT. ALLOCATED(lmp_v))THEN;ALLOCATE(lmp_v(RANK_THREE*nsphere));ENDIF
  IF(.NOT. ALLOCATED(lmp_o))THEN;ALLOCATE(lmp_o(RANK_THREE*nsphere));ENDIF
  IF(.NOT. ALLOCATED(lmp_img))THEN;ALLOCATE(lmp_img(RANK_THREE*nsphere));ENDIF

  !--Initiallize the particle 
  DO nn = 1,nsphere
    lmp_x((nn-1)*RANK_THREE+1) = objcentx(nn)
    lmp_x((nn-1)*RANK_THREE+2) = objcenty(nn)
    lmp_x((nn-1)*RANK_THREE+3) = zero

    lmp_v((nn-1)*RANK_THREE+1) = objcentu(nn)
    lmp_v((nn-1)*RANK_THREE+2) = objcentv(nn)
    lmp_v((nn-1)*RANK_THREE+3) = zero
    
    lmp_o((nn-1)*RANK_THREE+1) = zero
    lmp_o((nn-1)*RANK_THREE+2) = zero
    lmp_o((nn-1)*RANK_THREE+3) = fd2lamOmConv(nn)*objcentom(nn)

  ENDDO
  
  lmp_img = 0
 
  WRITE(*,*)'Before lammps Linear Velcoities: ', lmp_v
  WRITE(*,*)'Before lammps Positions: ', lmp_x
  WRITE(*,*)'Before lammps Angular Velocity: ', lmp_o
  WRITE(*,*)'Before lammps Image Flags: ',lmp_img

  CALL lammps_scatter_atoms(lmp_ptr, 'x', lmp_x)
  CALL lammps_scatter_atoms(lmp_ptr, 'v', lmp_v)
  CALL lammps_scatter_atoms(lmp_ptr, 'omega', lmp_o)
  
  WRITE(lmp_command,'(A)')'set type 1 image 0 0 0'
  CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))

  WRITE(lmp_command,'(A)')'reset_timestep 0'
  CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))

  WRITE(lmp_command,'(A,1X,I5)')'run',lmp_nstep
  CALL lammps_command(lmp_ptr,TRIM(ADJUSTL(lmp_command)))
  
  CALL lammps_gather_atoms(lmp_ptr, 'x', RANK_THREE, lmp_x)
  CALL lammps_gather_atoms(lmp_ptr, 'v', RANK_THREE, lmp_v)
  CALL lammps_gather_atoms(lmp_ptr, 'omega', RANK_THREE, lmp_o)
  CALL lammps_gather_atoms_image(lmp_ptr, lmp_img)
  
  WRITE(*,*)'Linear Velcoities: ', lmp_v
  WRITE(*,*)'Positions: ', lmp_x
  WRITE(*,*)'Angular Velocity: ', lmp_o
  WRITE(*,*)'Image Flags: ',lmp_img
  
  DO nn = 1,nsphere
    objcentx(nn) = lmp_x((nn-1)*RANK_THREE+1)
    objcenty(nn) = lmp_x((nn-1)*RANK_THREE+2)

    objcentu(nn) = lmp_v((nn-1)*RANK_THREE+1)
    objcentv(nn) = lmp_v((nn-1)*RANK_THREE+2)
    
    objcentom(nn) = lmp_o((nn-1)*RANK_THREE+3)/fd2lamOmConv(nn)
    
    zobjcentx(nn) = zobjcentx(nn) + lmp_img((nn-1)*RANK_THREE+1)
    zobjcenty(nn) = zobjcenty(nn) + lmp_img((nn-1)*RANK_THREE+2)
  ENDDO
 
  CALL fd_calc_ori

  DO nn = 1,nsphere
    
    moved = 0
    
    lindispx = objcentx(nn) - objcentxo(nn) + (zobjcentx(nn) - zobjcentxo(nn))*LDomainx
    lindispy = objcenty(nn) - objcentyo(nn) + (zobjcenty(nn) - zobjcentyo(nn))*LDomainy
    
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

CALL fd_find_rigidforce_contrib

dxmoved = objcentx - objcentxo + (zobjcentx - zobjcentxo)*LDomainx

dxmeanmoved = dxmeanmoved + SUM(dxmoved) /REAL(nsphere)
100 CONTINUE
IF(error == 1)THEN
  WRITE(*,*)'fd_calc_pos: Node Lost in x-dir...'
ELSEIF(error == 2)THEN 
  WRITE(*,*)'fd_calc_pos: Node Lost in y-dir...' 
ENDIF

END SUBROUTINE lmp_calc_pos

END MODULE modlmp_particle
