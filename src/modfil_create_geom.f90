Module modfil_create_geom

PRIVATE
PUBLIC :: fil_create_geom,fil_sol_tension,fil_sol_positions,fil_calc_sources

CONTAINS
!================================================================
SUBROUTINE fil_create_geom

USE parameters,         ONLY : alloc_create
USE precision,          ONLY : r_single
USE real_parameters,    ONLY : two,pi,small,one,half,zero
USE shared_data,        ONLY : filds,nfil,nfilpoints,filpointx,filpointy,fillen,&
                               nfilpoints,filfirstposx,filfirstposy,fillasttheta,&
                               filpointxo,filpointyo,filpointxpen,filpointypen,&
                               !###############For oscillating sphere test##############
                               nsurfpoints,objradius,calcsurfforce,nsphere
                               !########################################################
IMPLICIT NONE
INTEGER :: n,j
REAL(KIND = r_single) :: ds,s,dtheta,theta

!################For oscillating sphere test###############
!nsphere = nfil
!##########################################################
DO n = 1,nfil
  ds = fillen(n) / REAL(nfilpoints(n)-1)
  DO j = 1,nfilpoints(n)
    s = REAL(j-1)*ds
    filpointx(j,n) = filfirstposx(n) + (fillen(n) - s)*COS(fillasttheta(n))
    filpointy(j,n) = filfirstposy(n) + (fillen(n) - s)*SIN(fillasttheta(n))
  ENDDO
 !################For oscillating sphere test###############
!  dtheta = two*pi / REAL(nfilpoints(n))
!  DO j=1,nfilpoints(n)
!  
!    theta = REAL(j-1)*dtheta
!    filpointx(j,n) = filfirstposx(n)+fillen(n)*COS(theta)
!    filpointy(j,n) = filfirstposy(n)+fillen(n)*SIN(theta)
!
!  ENDDO
!##########################################################
ENDDO

filpointxpen = filpointx
filpointypen = filpointy  

filpointxo = filpointx
filpointyo = filpointy

DO n = 1,nfil
  DO j = 1,nfilpoints(n)-1
	filds(j,n)=SQRT((filpointx(j,n) - filpointx(j+1,n))**2 + (filpointy(j,n) - filpointy(j+1,n))**2)
  ENDDO
  filds(j,n) = filds(j-1,n)
ENDDO

!#########################For oscillating sphere test ONLY###########################
!allocate(nsurfpoints(nfil),objradius(nfil))
!objradius = fillen
!nsurfpoints = nfilpoints
!if(calcsurfforce)then
!  call fd_alloc_surfforce_arrays(alloc_create,maxval(nfilpoints(:),1))
!  call fd_create_interpmol
!endif
!nsphere = 0
!####################################################################################

END SUBROUTINE fil_create_geom

SUBROUTINE fd_create_interpmol
  
USE real_parameters,ONLY : one,two,small,half,zero
USE precision,      ONLY : r_single
USE shared_data,    ONLY :   nfil,nfilpoints,filpointx,filpointy,filfirstposx,filfirstposy,&
                               surfds,surfnx,surfny,surfcentrex,surfcentrey,&
                               xc,yc,surfpoint_cvx,surfpoint_cvy,&
                               surfpoint_interp,surfinterpx,surfinterpy,xc,yc,x,y,nsphere,&
                               nsurfpoints,objradius,nsphere
USE modfd_create_geom,  ONLY : find_single_point

IMPLICIT NONE
INTEGER :: n,ip1
REAL(KIND = r_single) :: vecx,vecy,coeft,dxmeanl,xp,yp
LOGICAL               :: failed
INTEGER               :: i,ip,jp,is,ie,js,je,ii,jj

DO n = 1,nfil
  DO i=1,nfilpoints(n)
      
      IF(i < nfilpoints(n))THEN
        ip1 = i + 1
      ELSE
        ip1 = 1
      ENDIF

      !--Calc centre points
      surfcentrex(i,n) = (filpointx(i,n) + filpointx(ip1,n))/two
      surfcentrey(i,n) = (filpointy(i,n) + filpointy(ip1,n))/two

      !--Locate the centre point on the Eul grid
       CALL find_single_point(surfcentrex(i,n), surfcentrey(i,n), surfpoint_cvx(1,i,n),surfpoint_cvy(1,i,n))
      

      vecx = filpointx(i,n) - filpointx(ip1,n)
      vecy = filpointy(i,n) - filpointy(ip1,n)
       
      surfds(i,n) = SQRT(vecx**2 + vecy**2)
      IF(surfds(i,n) < small) THEN
        PRINT *,'Error - Coordinates co-exits - undefined normal!'
        PAUSE
      ENDIF
      
      !--Surface element normals (outward)
      surfnx(i,n) = -vecy/surfds(i,n)
      surfny(i,n) = vecx/surfds(i,n)

      !--Ensure the calculation is corrent
      IF((surfcentrex(i,n) - filfirstposx(n))*surfnx(i,n) + &
         (surfcentrey(i,n) - filfirstposy(n))*surfny(i,n) < zero)THEN
         PRINT *,'Error - outward normal not calculated!'
         PAUSE
      ENDIF

    ENDDO

    DO i=1,nfilpoints(n)
      
      coeft = one
      dxmeanl = ( (x(surfpoint_cvx(1,i,n)) - x(surfpoint_cvx(1,i,n)-1)) + & 
                  (y(surfpoint_cvy(1,i,n)) - y(surfpoint_cvy(1,i,n)-1))   )/two
      DO          
        failed = .FALSE. 
        !--first interpolation point
        xp = surfcentrex(i,n) + coeft * surfnx(i,n) * dxmeanl 
        yp = surfcentrey(i,n) + coeft * surfny(i,n) * dxmeanl 

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
            IF( ( (xc(ii) - surfcentrex(i,n))*surfnx(i,n) + &
                  (yc(jj) - surfcentrey(i,n))*surfny(i,n)    ) < zero )THEN

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
  ENDDO

END SUBROUTINE fd_create_interpmol

SUBROUTINE fil_calc_tension_rhs(n)

!--prepare right hand side of the filament equation 
!--W.-X. Huang et al. / Journal of Computational Physics 226 (2007) 2206–2228

USE precision,              ONLY : r_single
USE real_parameters,        ONLY : zero,five,two,ten,three,nine,one,four
USE shared_data,            ONLY : filst,filpointxpen,filpointx,filpointxo,&
                                   filpointypen,filpointy,filpointyo,densitfil,&
                                   nfilpoints,filds,dt,filru,filrv,filfx,filfy,&
                                   filgravx,filgravy,filbenrig,filfr

IMPLICIT NONE
INTEGER,INTENT(IN)  :: n !--current filament

INTEGER :: j

REAL(KIND = r_single) :: vec1x,vec2x,vec1y,vec2y,bndvecx,bndvecy

  filst = zero 

  filpointxpen(:,n) = two*filpointx(:,n) - filpointxo(:,n)
  filpointypen(:,n) = two*filpointy(:,n) - filpointyo(:,n)
  
  !first term contribution
  DO j=1,nfilpoints(n)-1
    vec1x = filpointx(j+1,n)-filpointx(j,n)
    vec1y = filpointy(j+1,n)-filpointy(j,n)
    vec2x = filpointxo(j+1,n)-filpointxo(j,n)
    vec2y = filpointyo(j+1,n)-filpointyo(j,n)
    filst(j) = densitfil(n)*(one - two* fil_calc_dotprod(vec1x,vec1x,vec1y,vec1y)/filds(j,n)**2 + &
                  fil_calc_dotprod(vec2x,vec2x,vec2y,vec2y)/filds(j,n)**2)/two/dt**2
  ENDDO
  
  !second term contribution
  DO j = 1,nfilpoints(n)-1
    vec1x = filru(j+1,n) - filru(j,n)
    vec1y = filrv(j+1,n) - filrv(j,n)
    filst(j) = filst(j) -  densitfil(n)*fil_calc_dotprod(vec1x,vec1x,vec1y,vec1y)/filds(j,n)**2
  ENDDO

  !Third term
  DO j=1,nfilpoints(n)-1
    vec1x = filpointxpen(j+1,n) - filpointxpen(j,n)
    vec1y = filpointypen(j+1,n) - filpointypen(j,n)
    
    IF(j > 2 .AND. j < nfilpoints(n)-2)THEN
            
      vec2x = filpointxpen(j+3,n) - five*filpointxpen(j+2,n) + &
              ten*filpointxpen(j+1,n) - ten*filpointxpen(j,n) + &
              five*filpointxpen(j-1,n) - filpointxpen(j-2,n)
      vec2y = filpointypen(j+3,n) - five*filpointypen(j+2,n) + &
              ten*filpointypen(j+1,n) - ten*filpointypen(j,n) + &
              five*filpointypen(j-1,n) - filpointypen(j-2,n)
        
    ELSEIF(j == 1 )THEN
      
      vec2x = - filpointxpen(j+2,n) + two*filpointxpen(j+1,n) - filpointxpen(j,n)
      vec2y = - filpointypen(j+2,n) + two*filpointypen(j+1,n) - filpointypen(j,n)
        
    ELSEIF(j == 2 )THEN
        
      vec2x = filpointxpen(j+3,n) - five*filpointxpen(j+2,n) + &
             ten*filpointxpen(j+1,n) - nine*filpointxpen(j,n) + &
             three*filpointxpen(j-1,n)
      vec2y = filpointypen(j+3,n) - five*filpointypen(j+2,n) + &
             ten*filpointypen(j+1,n) - nine*filpointypen(j,n) + &
             three*filpointypen(j-1,n)

    ELSEIF(j == nfilpoints(n)-2 )THEN
        
      vec2x = - three*filpointxpen(j+2,n) + nine*filpointxpen(j+1,n) - &
             ten*filpointxpen(j,n) + five*filpointxpen(j-1,n) - filpointxpen(j-2,n) 
      vec2y = - three*filpointypen(j+2,n) + nine*filpointypen(j+1,n) - &
             ten*filpointypen(j,n) + five*filpointypen(j-1,n) - filpointypen(j-2,n) 
        
     ELSEIF(j == nfilpoints(n)-1 )THEN
            
       bndvecx = filfx(j+1,n) - densitfil(n)*filgravx !filfr/SQRT(fil_calc_dotprod(filgravx,filgravx,filgravy,filgravy))*filgravx
       bndvecy = filfy(j+1,n) - densitfil(n)*filgravy !filfr/SQRT(fil_calc_dotprod(filgravx,filgravx,filgravy,filgravy))*filgravy
       
       vec2x =  - two*filpointxpen(j+1,n) + five*filpointxpen(j,n) - &
               four*filpointxpen(j-1,n) + filpointxpen(j-2,n) 
       vec2y =  - two*filpointypen(j+1,n) + five*filpointypen(j,n) - &
               four*filpointypen(j-1,n) + filpointypen(j-2,n) 
     END IF
        
     IF( j /= nfilpoints(n)-1)THEN
       filst(j) = filst(j) + filbenrig(n) * fil_calc_dotprod(vec1x,vec2x,vec1y,vec2y)/filds(j,n)**6
     ELSE
       filst(j) = filst(j) - (filbenrig(n) * fil_calc_dotprod(vec1x,vec2x,vec1y,vec2y)/filds(j,n)**6 + &
                              fil_calc_dotprod(vec1x,bndvecx,vec1y,bndvecy)/filds(j,n)**2)
     ENDIF
   ENDDO

   !Fourth term 
   DO j=1,nfilpoints(n)-1
     bndvecx = filfx(j+1,n) - filfx(j,n)
     bndvecy = filfy(j+1,n) - filfy(j,n)
     vec1x = filpointxpen(j+1,n) - filpointxpen(j,n)
     vec1y = filpointypen(j+1,n) - filpointypen(j,n)
     filst(j) = filst(j) + fil_calc_dotprod(vec1x,bndvecx,vec1y,bndvecy)/filds(j,n)**2
   ENDDO

END SUBROUTINE fil_calc_tension_rhs

SUBROUTINE fil_calc_tension_coeftmat(n)

USE real_parameters,        ONLY : two,three,zero
USE precision,              ONLY : r_single
USE shared_data,            ONLY : filae,filaw,filap,nfilpoints,filpointxpen,filpointypen,filds

IMPLICIT NONE 

INTEGER,INTENT(IN)      :: n !--current filamnet
REAL(KIND = r_single)   :: vec1x,vec1y,vec2x,vec2y,vec3x,vec3y
INTEGER :: j

    
 filae = zero
 filaw = zero
 filap = zero
   
 DO j=1,nfilpoints(n) - 1
   IF(j > 1 .AND. j < nfilpoints(n) - 1)THEN
            
     vec1x = filpointxpen(j+2,n) - filpointxpen(j+1,n)
     vec1y = filpointypen(j+2,n) - filpointypen(j+1,n)
     vec2x = filpointxpen(j+1,n) - filpointxpen(j,n)
     vec2y = filpointypen(j+1,n) - filpointypen(j,n)
     vec3x = filpointxpen(j,n) - filpointxpen(j-1,n)
     vec3y = filpointypen(j,n) - filpointypen(j-1,n)

     filaw(j) = fil_calc_dotprod(vec3x,vec2x,vec3y,vec2y)/filds(j,n)**4
     filap(j) = -two*fil_calc_dotprod(vec2x,vec2x,vec2y,vec2y)/filds(j,n)**4
     filae(j) = fil_calc_dotprod(vec1x,vec2x,vec1y,vec2y)/filds(j,n)**4
        
   ELSE IF(j == 1 )THEN
        
     vec1x = filpointxpen(j+2,n) - filpointxpen(j+1,n)
     vec1y = filpointypen(j+2,n) - filpointypen(j+1,n)
     vec2x = filpointxpen(j+1,n) - filpointxpen(j,n)
     vec2y = filpointypen(j+1,n) - filpointypen(j,n)
     filaw(j) = zero
     filap(j) = -three*fil_calc_dotprod(vec2x,vec2x,vec2y,vec2y)/filds(j,n)**4
     filae(j) = fil_calc_dotprod(vec1x,vec2x,vec1y,vec2y)/filds(j,n)**4
   ELSE IF(j == nfilpoints(n) -1 )THEN
        
     vec2x = filpointxpen(j+1,n) - filpointxpen(j,n)
     vec2y = filpointypen(j+1,n) - filpointypen(j,n)
     vec3x = filpointxpen(j,n) - filpointxpen(j-1,n)
     vec3y = filpointypen(j,n) - filpointypen(j-1,n) 
              
     filaw(j) = fil_calc_dotprod(vec3x,vec2x,vec3y,vec2y)/filds(j,n)**4
     filap(j) = - fil_calc_dotprod(vec2x,vec2x,vec2y,vec2y)/filds(j,n)**4
     filae(j) = zero
   END IF
        
 ENDDO

END SUBROUTINE fil_calc_tension_coeftmat

SUBROUTINE fil_sol_tension

USE shared_data,            ONLY : nfil,nfilpoints,filst,filae,filap,filaw,filten
USE modfd_solve_linearsys,  ONLY : fil_sol_tdma

IMPLICIT NONE

INTEGER         :: n

DO n = 1,nfil
  CALL fil_calc_tension_rhs(n)
  CALL fil_calc_tension_coeftmat(n)
  CALL fil_sol_tdma(filten(:,n),filaw,filap,filae,filst,nfilpoints(n)-1)

ENDDO

END SUBROUTINE fil_sol_tension

SUBROUTINE fil_calc_pos_rhs(n)

USE precision,                  ONLY : r_single
USE real_parameters,            ONLY : zero,six,three,four,five,two
USE shared_data,                ONLY : densitfil,filap,filaw,filae,filsx,filbenrig,&
                                       filpointx,filpointy,filpointyo,filpointxo,&
                                       filds,filpointxpen,filpointypen,filsy,dt,&
                                       filgravx,filgravy,filfx,filfy,dt,nfilpoints,filfr
IMPLICIT NONE 
INTEGER,INTENT(IN)      :: n !--current filament

REAL(KIND = r_single)   :: gx,gy
INTEGER :: j
    
filsx = zero
filsy = zero 

gx = filgravx/SQRT(fil_calc_dotprod(filgravx,filgravx,filgravy,filgravy))
gy = filgravy/SQRT(fil_calc_dotprod(filgravx,filgravx,filgravy,filgravy))
     
DO j=1,nfilpoints(n)

  IF(j > 2 .AND. j < nfilpoints(n)-1)THEN
            
    filsx(j) = -filbenrig(n)*(filpointxpen(j+2,n)-four*filpointxpen(j+1,n)+&
               six*filpointxpen(j,n)-four*filpointxpen(j-1,n)+filpointxpen(j-2,n))*dt**2/filds(j,n)**4+&
               (densitfil(n)*filgravx-filfx(j,n))*dt**2+ densitfil(n)*filpointxpen(j,n)
            
    filsy(j) = -filbenrig(n)*(filpointypen(j+2,n)-four*filpointypen(j+1,n)+&
               six*filpointypen(j,n)-four*filpointypen(j-1,n)+filpointypen(j-2,n))*dt**2/filds(j,n)**4+&
               (densitfil(n)*filgravy-filfy(j,n))*dt**2+ densitfil(n)*filpointypen(j,n)
            
  ELSEIF(j == 1 )THEN
        
    filsx(j) = -filbenrig(n)*(filpointxpen(j+3,n)-three*filpointxpen(j+2,n)+&
                three*filpointxpen(j+1,n)-filpointxpen(j,n))*dt**2/filds(j,n)**4+&
                (densitfil(n)*filgravx-filfx(j,n))*dt**2+ densitfil(n)*filpointxpen(j,n)
                !(filfr*gx-filfx(j,n))*dt**2+ filpointxpen(j,n) !--nondimmed form
            
    filsy(j) = -filbenrig(n)*(filpointypen(j+3,n)-three*filpointypen(j+2,n)+&
                three*filpointypen(j+1,n)-filpointypen(j,n))*dt**2/filds(j,n)**4+&
                (densitfil(n)*filgravy-filfy(j,n))*dt**2+ densitfil(n)*filpointypen(j,n)
                !(filfr*gy-filfy(j,n))*dt**2+ filpointypen(j,n)
            
  ELSEIF(j == 2 )THEN
        
    filsx(j) = -filbenrig(n)*(filpointxpen(j+2,n)-four*filpointxpen(j+1,n)+&
                five*filpointxpen(j,n)-two*filpointxpen(j-1,n))*dt**2/filds(j,n)**4+&
                (densitfil(n)*filgravx-filfx(j,n))*dt**2+ densitfil(n)*filpointxpen(j,n)
                !(filfr*gx-filfx(j,n))*dt**2+ filpointxpen(j,n)
            
    filsy(j) = -filbenrig(n)*(filpointypen(j+2,n)-four*filpointypen(j+1,n)+&
                five*filpointypen(j,n)-two*filpointypen(j-1,n))*dt**2/filds(j,n)**4+&
                (densitfil(n)*filgravy-filfy(j,n))*dt**2+ densitfil(n)*filpointypen(j,n)
                !(filfr*gy-filfy(j,n))*dt**2+ filpointypen(j,n)
            
  ELSEIF(j == nfilpoints(n)-1 )THEN
        
    filsx(j) = -filbenrig(n)*(-two*filpointxpen(j+1,n)+&
                five*filpointxpen(j+1,n)-four*filpointxpen(j-1,n)+filpointxpen(j-2,n))*dt**2/filds(j,n)**4+&
                (densitfil(n)*filgravx-filfx(j,n))*dt**2+ densitfil(n)*filpointxpen(j,n)
                !(filfr*gx-filfx(j,n))*dt**2+ filpointxpen(j,n)
            
    filsy(j) = -filbenrig(n)*(-two*filpointypen(j+1,n)+&
                five*filpointypen(j,n)-four*filpointypen(j-1,n)+filpointypen(j-2,n))*dt**2/filds(j,n)**4+&
                (densitfil(n)*filgravy-filfy(j,n))*dt**2+ densitfil(n)*filpointypen(j,n)
                !(filfr*gy-filfy(j,n))*dt**2+ filpointypen(j,n)
        
  ELSEIF(j == nfilpoints(n))THEN
        
    filsx(j) = densitfil(n)*filpointxpen(j,n)
            
    filsy(j) = densitfil(n)*filpointypen(j,n)
        
  ENDIF
ENDDO
         
END SUBROUTINE fil_calc_pos_rhs

SUBROUTINE fil_calc_pos_coeftmat(n)

USE real_parameters,        ONLY : zero,two,one
USE shared_data,            ONLY : filaw,filap,filae,nfilpoints,filten,dt,filds,&
                                   densitfil,filten

IMPLICIT NONE

INTEGER,INTENT(IN)  :: n
INTEGER             :: j

filae = zero
filaw = zero
filap = zero

DO j=1,nfilpoints(n)

  IF(j > 1 .AND. j < nfilpoints(n))THEN
            
    filaw(j) = -filten(j-1,n)*(dt/filds(j,n))**2
    filap(j) = densitfil(n)+(filten(j-1,n)+filten(j,n))*(dt/filds(j,n))**2
    filae(j) = -filten(j,n)*(dt/filds(j,n))**2
                    
  ELSEIF(j == 1 )THEN
        
    filaw(j) = zero
    filap(j) = densitfil(n)+two*filten(j,n)*(dt/filds(j,n))**2
    filae(j) = -two*filten(j,n)*(dt/filds(j,n))**2
            
  ELSE IF(j == nfilpoints(n) )THEN
            
    filaw(j) = zero
    filap(j) = densitfil(n)
    filae(j) = zero
            
  END IF
        
ENDDO

END SUBROUTINE fil_calc_pos_coeftmat

SUBROUTINE fil_sol_positions

USE shared_data,            ONLY : nfil,filae,filaw,filap,filpointx,filsx,nfilpoints,&
                                   filpointxo,filpointy,filpointyo,filsy
USE modfd_solve_linearsys,  ONLY : fil_sol_tdma2

IMPLICIT NONE

INTEGER     :: n

DO n = 1,nfil
    
    CALL fil_calc_pos_rhs(n)
    CALL fil_calc_pos_coeftmat(n)

    filpointxo(:,n) = filpointx(:,n)
    filpointyo(:,n) = filpointy(:,n)

    CALL fil_sol_tdma2(filpointx(:,n),filpointy(:,n),filaw,filap,filae,filsx,filsy,nfilpoints(n))
    
ENDDO

END SUBROUTINE fil_sol_positions

FUNCTION fil_calc_dotprod(vec1x,vec2x,vec1y,vec2y)

USE precision,          ONLY : r_single

IMPLICIT NONE

REAL(KIND = r_single) :: fil_calc_dotprod
REAL(KIND = r_single),INTENT(IN) :: vec1x,vec2x,vec1y,vec2y

fil_calc_dotprod = vec1x*vec2x + vec1y*vec2y

END FUNCTION fil_calc_dotprod

SUBROUTINE  fil_find_nodes(n)

USE shared_data,    ONLY : filpoint_cvx,filpoint_cvy,x,y,nim,njm,&
                           filpointx,filpointy,nfilpoints,nfil,&
                           filpoint_interpx,filpoint_interpy,xc,yc,deltalen

IMPLICIT NONE
INTEGER         :: i,j,n

DO j = 1,nfilpoints(n)
  DO i = 2,nim
    IF(filpointx(j,n) <= x(i) .AND. filpointx(j,n) > x(i-1))THEN
      filpoint_cvx(j,n) = i
      EXIT
    ENDIF
  ENDDO

  DO i = 2,njm
    IF(filpointy(j,n) <= y(i) .AND. filpointy(j,n) > y(i-1))THEN
      filpoint_cvy(j,n) = i
      EXIT
    ENDIF
  ENDDO

ENDDO

IF(deltalen == 3)THEN
  DO j = 1,nfilpoints(n)
    filpoint_interpx(1,j,n) = filpoint_cvx(j,n) - 1
    filpoint_interpx(2,j,n) = filpoint_cvx(j,n) + 1
    filpoint_interpy(1,j,n) = filpoint_cvy(j,n) - 1
    filpoint_interpy(2,j,n) = filpoint_cvy(j,n) + 1
  ENDDO
ELSEIF(deltalen == 2 .OR. deltalen == 4 .OR. deltalen == 6)THEN
  DO j = 1,nfilpoints(n)
    IF(filpointx(j,n) <= xc(filpoint_cvx(j,n)))THEN
      filpoint_interpx(1,j,n) = filpoint_cvx(j,n) - deltalen/2
      filpoint_interpx(2,j,n) = filpoint_cvx(j,n) + deltalen/2 - 1
    ELSE
      filpoint_interpx(1,j,n) = filpoint_cvx(j,n) - deltalen/2 + 1
      filpoint_interpx(2,j,n) = filpoint_cvx(j,n) + deltalen/2
    ENDIF
    IF(filpointy(j,n) <= yc(filpoint_cvy(j,n)))THEN
      filpoint_interpy(1,j,n) = filpoint_cvy(j,n) - deltalen/2
      filpoint_interpy(2,j,n) = filpoint_cvy(j,n) + deltalen/2 - 1
    ELSE
      filpoint_interpy(1,j,n) = filpoint_cvy(j,n) - deltalen/2 + 1
      filpoint_interpy(2,j,n) = filpoint_cvy(j,n) + deltalen/2
    ENDIF
  ENDDO
ENDIF  
END SUBROUTINE fil_find_nodes

SUBROUTINE fil_update_realvel(n)

USE shared_data,            ONLY : filru,filrv,filpointx,filpointy,nfilpoints,dt,&
                                   filpointxo,filpointyo
IMPLICIT NONE

INTEGER,INTENT(IN)   :: n
INTEGER :: j

DO j=1,nfilpoints(n)
    filru(j,n) = (filpointx(j,n)-filpointxo(j,n))/dt
    filrv(j,n) = (filpointy(j,n)-filpointyo(j,n))/dt
END DO
   
END SUBROUTINE fil_update_realvel

SUBROUTINE fil_calc_len(n)

USE real_parameters,            ONLY : zero
USE shared_data,                ONLY : fillen,filpointx,filpointy,nfilpoints

IMPLICIT NONE

INTEGER,INTENT(IN)  :: n
INTEGER :: j

fillen(n) = zero

DO j=1,nfilpoints(n)-1

fillen(n) = fillen(n) + SQRT((filpointx(j,n) - filpointx(j+1,n))**2 + (filpointy(j,n) - filpointy(j+1,n))**2)

END DO

END SUBROUTINE fil_calc_len

SUBROUTINE fil_calc_sources(delt,resor,iter)
use omp_lib
USE shared_data,    ONLY : fdsu,fdsv,xc,yc,dt,li,u,v,nij,nim,njm,nj,nfil,filpointx,filpointy,&
                           filds,filu,filv,filru,filrv,filfx,filfy,filintegfx,filintegfy,&
                           filcintegfx,filcintegfy,filalpha,filbeta,filpoint_cvx,filpoint_cvy,&
                           fdsuc,fdsvc,filalpha,filbeta,nfilpoints,x,y,filpoint_interpx,filpoint_interpy,&
                           ibsu,ibsv
USE precision,      ONLY : r_single
USE real_parameters,ONLY : zero,one,two
USE modfd_create_geom,ONLY : fd_calc_delta

IMPLICIT NONE
INTEGER,INTENT(IN) :: iter
INTEGER         :: nn,n,j,i
REAL(KIND = r_single),INTENT(OUT) :: resor
DOUBLE PRECISION,INTENT(OUT)      :: delt
REAL(KIND = r_single)             :: del,dx,dy,temp
DOUBLE PRECISION                  :: T1,T2

resor = zero
filu = zero
filv = zero
fdsuc = zero
fdsvc = zero
delt = zero

IF(iter == 0)THEN
  DO nn = 1,nfil
    !CALL fil_calc_forcedmotion(nn) !--this is just for testing the oscillating circular cylinder to test the VB aspect of the code
    CALL fd_create_interpmol
    CALL fil_update_realvel(nn)
    CALL fil_find_nodes(nn)
    DO n = 1,nfilpoints(nn)
      DO i = filpoint_interpx(1,n,nn),filpoint_interpx(2,n,nn)
        dx = x(i) - x(i-1)
        DO j = filpoint_interpy(1,n,nn),filpoint_interpy(2,n,nn)
          dy = y(j) - y(j-1)
          del = fd_calc_delta(xc(i),yc(j),filpointx(n,nn),filpointy(n,nn),dx, dy)*dx*dy
          filu(n,nn) = filu(n,nn) + u(li(i)+j)*del
          filv(n,nn) = filv(n,nn) + v(li(i)+j)*del
        ENDDO
      ENDDO 
    ENDDO
    DO n = 1,nfilpoints(nn)
      temp = ABS(filu(n,nn)) ! - filru(n,nn))/ABS(filru(n,nn))
      IF(temp > resor) resor = temp 
      filcintegfx(n,nn) = filalpha * (filu(n,nn) - filru(n,nn)) * dt
      filcintegfy(n,nn) = filalpha * (filv(n,nn) - filrv(n,nn)) * dt
      filintegfx(n,nn) = filintegfx(n,nn) + filcintegfx(n,nn) 
      filfx(n,nn) = filintegfx(n,nn) + filbeta*(filu(n,nn) - filru(n,nn))
      filintegfy(n,nn) = filintegfy(n,nn) + filcintegfy(n,nn) 
      filfy(n,nn) = filintegfy(n,nn) + filbeta*(filv(n,nn) - filrv(n,nn))
    ENDDO

    DO n = 1,nfilpoints(nn)
      DO i = filpoint_interpx(1,n,nn),filpoint_interpx(2,n,nn)
        dx = x(i) - x(i-1)
        DO j = filpoint_interpy(1,n,nn),filpoint_interpy(2,n,nn)
          dy = y(j) - y(j-1)
          del = fd_calc_delta(xc(i),yc(j),filpointx(n,nn),filpointy(n,nn),dx, dy)*filds(n,nn)
          fdsuc(li(i)+j) = fdsuc(li(i)+j) + filfx(n,nn)*del*dx*dy
          fdsvc(li(i)+j) = fdsvc(li(i)+j) + filfy(n,nn)*del*dx*dy
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  ibsu = fdsuc
  ibsv = fdsvc

ELSEIF(iter > 0)THEN
  !T1 = omp_get_wtime()
  DO nn = 1,nfil
    DO n = 1,nfilpoints(nn)
      DO i = filpoint_interpx(1,n,nn),filpoint_interpx(2,n,nn)
        dx = x(i) - x(i-1)
        DO j = filpoint_interpy(1,n,nn),filpoint_interpy(2,n,nn)
          dy = y(j) - y(j-1)
          del = fd_calc_delta(xc(i),yc(j),filpointx(n,nn),filpointy(n,nn),dx, dy)*dx*dy
          filu(n,nn) = filu(n,nn) + u(li(i)+j)*del
          filv(n,nn) = filv(n,nn) + v(li(i)+j)*del
        ENDDO
      ENDDO 
    ENDDO
    !T2 = omp_get_wtime()
    !delt = T2-T1
    DO n = 1,nfilpoints(nn)
      temp = ABS(filu(n,nn)) !- filru(n,nn))/ABS(filru(n,nn))
      IF(temp > resor) resor = temp
      !filintegfx(n,nn) = filintegfx(n,nn) - filcintegfx(n,nn)
      !filintegfy(n,nn) = filintegfy(n,nn) - filcintegfy(n,nn)
!      filcintegfx(n,nn) = filalpha * (filu(n,nn) - filru(n,nn)) * dt
!      filcintegfy(n,nn) = filalpha * (filv(n,nn) - filrv(n,nn)) * dt
!      filintegfx(n,nn) = filintegfx(n,nn) + filcintegfx(n,nn) 
!      filfx(n,nn) = filintegfx(n,nn) + filbeta*(filu(n,nn) - filru(n,nn))
!      filintegfy(n,nn) = filintegfy(n,nn) + filcintegfy(n,nn) 
!      filfy(n,nn) = filintegfy(n,nn) + filbeta*(filv(n,nn) - filrv(n,nn))
    ENDDO
!
!    DO n = 1,nfilpoints(nn)
!      DO i = filpoint_cvx(n,nn) - 1,filpoint_cvx(n,nn) + 1
!        dx = x(i) - x(i-1)
!        DO j = filpoint_cvy(n,nn) - 1,filpoint_cvy(n,nn) + 1
!          dy = y(j) - y(j-1)
!          del = fd_calc_delta(xc(i),yc(j),filpointx(n,nn),filpointy(n,nn),dx, dy)*filds(n,nn)
!          fdsuc(li(i)+j) = fdsuc(li(i)+j) + filfx(n,nn)*del*dx*dy
!          fdsvc(li(i)+j) = fdsvc(li(i)+j) + filfy(n,nn)*del*dx*dy
!        ENDDO
!      ENDDO
!    ENDDO
  ENDDO
!
!  fdsu = fdsuc
!  fdsv = fdsvc
!    OPEN(UNIT = 1234,FILE='TIME.dat',STATUS = 'UNKNOWN',ACCESS = 'APPEND')
!    WRITE(1234,*) T2-T1
!    CLOSE(1234)
ENDIF
END SUBROUTINE fil_calc_sources

SUBROUTINE fil_calc_forcedmotion(n)

!--Move the cylinder and calculate the velocity, here filament geometry is actually a cylinder

USE precision,              ONLY : r_single
USE real_parameters,        ONLY : pi,fifth,five,two,zero,six,ten,one
USE shared_data,            ONLY : nfilpoints,filpointx,filpointy,filpointxo,filpointyo,&
                                   filfirstposx,filfirstposy,time,fillen,filpointxpen,filpointypen,&
                                   nfil

IMPLICIT NONE
INTEGER,INTENT(IN)    :: n
INTEGER               :: j

REAL(KIND = r_single) :: omega,amp,freq,KC,dx,starttime,t,freq0,freqfact

filpointxo(:,n) = filpointx(:,n)
filpointyo(:,n) = filpointy(:,n)

!--horizontal motion in stationary fluid
!freq = fifth
!KC   = five
!omega = two*pi*freq
!amp = fillen(n) * KC / pi
!
!dx = amp*SIN(omega * time) 
!DO j=1,nfilpoints(n) 
!  filpointx(j,n) = filpointxpen(j,n) - dx
!ENDDO 

!--tranverse motion 
!--natural shedding freq (depends on Re)
starttime = fifth
freq0 = 0.95; freqfact = 1.1
freq = freq0 * freqfact
amp = 0.4 * fillen(n)
omega = two*pi*freq
t = time - starttime
IF(t < zero) t = zero 
!dx = amp - amp*COS(omega * t)
!filfirstposy(n) = 0.52 - dx
!DO j=1,nfilpoints(n) 
!  filpointy(j,n) = filpointypen(j,n) - dx
!ENDDO
dx = amp*SIN(omega * t)
filfirstposy(n) = 1.0 - dx
DO j=1,nfilpoints(n) 
  filpointy(j,n) = filpointypen(j,n) - dx
ENDDO

END SUBROUTINE fil_calc_forcedmotion

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

END MODULE modfil_create_geom