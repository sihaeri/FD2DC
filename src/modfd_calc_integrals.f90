MODULE modfd_calc_integrals

PRIVATE
PUBLIC :: fd_calc_integrals,fd_calc_surf_force,fd_calc_surf_nusselt,fd_calc_wall_nusselt,&
          fd_calc_lwall_nusselt,fd_calc_surf_nusselt_ave,fd_calc_lwall_nusselt_ave
CONTAINS
!===============================================
SUBROUTINE fd_calc_integrals

!--Output of some integral quantities for checking
!--convergence towards grid-independent solution and
!--estimation of discretization errors
USE precision,        ONLY : r_single
USE real_parameters,  ONLY : zero,real_1e3,pi,one,three,ten,four
USE shared_data,      ONLY : ien,lcal,li,nim,njm,visc,r,y,xc,pp,f1,f2,ni,nj,t,imon,p,&
                             nsphere,lcal,ien,putobj,objtp,cpf,kappaf,flomas,dt,objradius,cpf,u,objbradius,celkappa
USE parameters,       ONLY : out_unit

IMPLICIT NONE 
REAL(KIND = r_single) :: qwall,s,d,psimin,psimax,pmeani,area,pmeano,tmeani,tmeano,tw,&
                         lmtd,dpmean,qm,num,hm,aspratio
INTEGER,Target               :: i,j,ij,n
!--HEAT FLUXES AT WEST AND EAST ISOTHERMAL WALLS

WRITE(out_unit,*) '  '

IF(lcal(ien)) THEN
  qwall=zero
  DO j=2,njm
    ij=li(1)+j
    s=0.5*(r(j)+r(j-1))*(y(j)-y(j-1))
    d=celkappa(ij)*s/(xc(2)-xc(1))
    qwall=qwall+d*(t(ij+nj)-t(ij))
  END DO
  WRITE(out_unit,*) '  HEAT FLUX THROUGH WEST WALL: ',qwall

  qwall=zero
  DO j=2,njm
    ij=li(ni)+j
    s=0.5*(r(j)+r(j-1))*(y(j)-y(j-1))
    d=celkappa(ij)*s/(xc(ni)-xc(nim))
    qwall=qwall+d*(t(ij)-t(ij-nj))
  END DO
  WRITE(out_unit,*) '  HEAT FLUX THROUGH EAST WALL: ',qwall
ENDIF

!--STREAMFUNCTION VALUES AT CV-VERTICES (ZERO AT SOUTH-WEST CORNER)

  pp(li(1)+1)=zero
  
  !--WEST BOUNDARY (APPLICABLE FOR INLET OR OUTLET)
  DO j=2,njm
    ij=li(1)+j
    pp(ij)=pp(ij-1)+f1(ij)
  END DO

  !--SOUTH BOUNDARY (APPLICABLE FOR INLET OR OUTLET)

  DO i=2,nim
    ij=li(i)+1
    pp(ij)=pp(ij)-f2(ij)
    !--INNER REGION
    DO j=2,njm
      ij=li(i)+j
      pp(ij)=pp(ij-1)+f1(ij)
    END DO
  END DO
  
  !--STRENGTH OF PRIMARY AND SECONDARY EDDY (MIN and MAX values)

  psimin= 1.e20
  psimax=-1.e20
  
  DO i=1,nim
    DO j=1,njm
      ij=li(i)+j
      psimin=min(psimin,pp(ij))
      psimax=max(psimax,pp(ij))
    END DO
  END DO
  
  IF(putobj .AND. nsphere > 0 .AND. lcal(ien))THEN
      !--Mean pressure at imon 
      pmeani = zero
      pmeano = zero
      area = zero
      DO j=2,njm
        ij=li(imon)+j
        s=0.5*(r(j)+r(j-1))*(y(j)-y(j-1))
        pmeani = pmeani + p(ij)*s
        area = area + s
      ENDDO
      pmeani = pmeani / area

      !--Mean pressure at outlet
      area = zero
      DO j=2,njm
        ij=li(nim)+j   
        s=0.5*(r(j)+r(j-1))*(y(j)-y(j-1))
        pmeano = pmeano + p(ij)*s
        area = area + s
      ENDDO
      pmeano = pmeano / area

      !--Mean temp at imon
      tmeani = zero
      tmeano = zero
      area = zero
      DO j=2,njm
        ij=li(imon)+j
        s=0.5*(r(j)+r(j-1))*(y(j)-y(j-1))*u(ij)
        tmeani = tmeani + t(ij)*s
        area = area + s
      ENDDO
      tmeani = tmeani / area

      area = zero
      DO j=2,njm
        ij=li(nim)+j
        s=0.5*(r(j)+r(j-1))*(y(j)-y(j-1))*u(ij)
        tmeano = tmeano + t(ij)*s
        area = area + s
      ENDDO
      tmeano = tmeano / area  
  
      dpmean = pmeani - pmeano

      !--Assuming all cylinders have same temperature
      tw = objtp(1)
      lmtd = ((tw - tmeani) - (tw - tmeano))/LOG((tw - tmeani)/(tw - tmeano)) 
      
      qm = cpf*flomas*(tmeano - tmeani)
      
      !--calc total heat transfer area
      area = zero
      DO n = 1,nsphere
        aspratio = ((objradius(n)-objbradius(n))/(objradius(n)+objbradius(n)))**2
        area = area + pi*(objradius(n)+objbradius(n))*(one + (three*aspratio)/(ten+SQRT(four - three*aspratio)))
      ENDDO

      !--calculate mean heat transfer coefficient and mean nusselt number
      hm = qm/area/lmtd

      num = hm * (y(nj) - y(1))/kappaf

      WRITE(out_unit,*) '  '
      WRITE(out_unit,*) '  Mean temp in, out, lmtd:  ',tmeani,tmeano,lmtd 
      WRITE(out_unit,*) '  Total pressure drop:  ',dpmean 
      WRITE(out_unit,*) '  Mean nusselt number:  ',num
  ENDIF

  WRITE(out_unit,*) '  '
  WRITE(out_unit,*) '  MAXIMUM STREAMFUNCTION VALUE:  ',psimax 
  WRITE(out_unit,*) '  MINIMUM STREAMFUNCTION VALUE:  ',psimin

END SUBROUTINE fd_calc_integrals

SUBROUTINE fd_calc_surf_force(n,cd,cl)

USE precision,              ONLY : r_single
USE real_parameters,        ONLY : one,zero,two
USE modfd_solve_linearsys,  ONLY : linsolve_gelem
USE shared_data,            ONLY : dtx,dty,dux,duy,dvx,dvy,surfpoint_cvx,surfpoint_cvy,surfpoint_interp,ulid,&
                                   surfinterpx,surfinterpy,xc,yc,surfds,surfnx,surfny,surfcentrex,surfcentrey,&
                                   visc,densit,objradius,nsurfpoints,p,li,lamvisc,temp_visc

IMPLICIT NONE
INTEGER,INTENT(IN)                :: n
REAL(KIND = r_single),INTENT(OUT) :: cd,cl

INTEGER                     :: i,j,k,myind(4),myind2(4,2)
REAL(KIND = r_single)       :: vandermonde(4,4),phiarr(4,6),coefts(4,6),detmt
REAL(KIND = r_single)       :: pressurei(2),duxi(2),duyi(2),dvxi(2),dvyi(2),pressures,duxs,duys,dvxs,dvys,&
                               fpx,fsx,fpy,fsy,extr_fact,hh,visi(2),viss

fpx = zero
fsx = zero
fpy = zero
fsy = zero

DO i = 1,nsurfpoints(n) 
  !--Set the coeft matrix for this facet for each interpolation point in turn
  DO k = 1,2
   myind(1) = li(surfpoint_interp(1,k,i,n)) + surfpoint_interp(3,k,i,n)
   myind(2) = li(surfpoint_interp(2,k,i,n)) + surfpoint_interp(3,k,i,n)
   myind(3) = li(surfpoint_interp(1,k,i,n)) + surfpoint_interp(4,k,i,n)
   myind(4) = li(surfpoint_interp(2,k,i,n)) + surfpoint_interp(4,k,i,n)
   myind2(1,1) = surfpoint_interp(1,k,i,n)
   myind2(1,2) = surfpoint_interp(3,k,i,n)
   myind2(2,1) = surfpoint_interp(2,k,i,n)
   myind2(2,2) = surfpoint_interp(3,k,i,n)
   myind2(3,1) = surfpoint_interp(1,k,i,n)
   myind2(3,2) = surfpoint_interp(4,k,i,n)
   myind2(4,1) = surfpoint_interp(2,k,i,n)
   myind2(4,2) = surfpoint_interp(4,k,i,n)
   DO j = 1,4
    vandermonde(j,1) = xc(myind2(j,1))*yc(myind2(j,2))
    vandermonde(j,2) = xc(myind2(j,1))
    vandermonde(j,3) = yc(myind2(j,2))
    vandermonde(j,4) = one
   ENDDO
  
  !--set the rhs (first one pressure)
  DO j = 1,4 
    phiarr(j,1) = p(myind(j))
  ENDDO
  !--2:5 are the 4 components of vel grad 
  DO j = 1,4
     phiarr(j,2) = dux(myind(j))
     phiarr(j,3) = duy(myind(j))
     phiarr(j,4) = dvx(myind(j))
     phiarr(j,5) = dvy(myind(j))
  ENDDO
  IF(temp_visc)THEN
    DO j = 1,4 
      phiarr(j,6) = lamvisc(myind(j))
    ENDDO
  ENDIF
                                                 
  IF(temp_visc)THEN
    CALL linsolve_gelem(vandermonde,phiarr,4,6,coefts,detmt)
  ELSE
    CALL linsolve_gelem(vandermonde,phiarr,4,5,coefts,detmt)
  ENDIF 
   
   pressurei(k) = coefts(1,1)*surfinterpx(k,i,n)*surfinterpy(k,i,n) + &
               coefts(2,1)*surfinterpx(k,i,n) + &
               coefts(3,1)*surfinterpy(k,i,n) + coefts(4,1)
  
   duxi(k) = coefts(1,2)*surfinterpx(k,i,n)*surfinterpy(k,i,n) + &
               coefts(2,2)*surfinterpx(k,i,n) + &
               coefts(3,2)*surfinterpy(k,i,n) + coefts(4,2)
   duyi(k) = coefts(1,3)*surfinterpx(k,i,n)*surfinterpy(k,i,n) + &
               coefts(2,3)*surfinterpx(k,i,n) + &
               coefts(3,3)*surfinterpy(k,i,n) + coefts(4,3)

   dvxi(k) = coefts(1,4)*surfinterpx(k,i,n)*surfinterpy(k,i,n) + &
               coefts(2,4)*surfinterpx(k,i,n) + &
               coefts(3,4)*surfinterpy(k,i,n) + coefts(4,4)
   dvyi(k) = coefts(1,5)*surfinterpx(k,i,n)*surfinterpy(k,i,n) + &
               coefts(2,5)*surfinterpx(k,i,n) + &
               coefts(3,5)*surfinterpy(k,i,n) + coefts(4,5)
  
    IF(temp_visc)THEN
      visi(k) = coefts(1,6)*surfinterpx(k,i,n)*surfinterpy(k,i,n) + &
                 coefts(2,6)*surfinterpx(k,i,n) + &
                 coefts(3,6)*surfinterpy(k,i,n) + coefts(4,6)
    ENDIF
  ENDDO
  
  !--Now extrapolate to the facet centre
  hh        = SQRT((surfinterpx(1,i,n) - surfinterpx(2,i,n))**2 + (surfinterpy(1,i,n) - surfinterpy(2,i,n))**2)
  extr_fact = SQRT((surfcentrex(i,n) - surfinterpx(2,i,n))**2 + (surfcentrey(i,n) - surfinterpy(2,i,n))**2) / hh
             

  pressures = pressurei(2) + extr_fact*(pressurei(1) - pressurei(2))
  duxs      = duxi(2)      + extr_fact*(duxi(1)      - duxi(2)     )
  duys      = duyi(2)      + extr_fact*(duyi(1)      - duyi(2)     )
  dvxs      = dvxi(2)      + extr_fact*(dvxi(1)      - dvxi(2)     )
  dvys      = dvyi(2)      + extr_fact*(dvyi(1)      - dvyi(2)     )

  IF(temp_visc)THEN
    viss = visi(2) + extr_fact*(visi(1) - visi(2))
  ELSE
    viss = visc
  ENDIF

  fpx = fpx + (- pressures * surfds(i,n) * surfnx(i,n))
  fpy = fpy + (- pressures * surfds(i,n) * surfny(i,n))
  fsx = fsx + viss * (two * duxs * surfnx(i,n) + (duys + dvxs) * surfny(i,n)) * surfds(i,n)
  fsy = fsy + viss * (two * dvys * surfny(i,n) + (duys + dvxs) * surfnx(i,n)) * surfds(i,n)
ENDDO

cd = (fpx + fsx)/(densit * ulid **2 * objradius(n))
cl = (fpy + fsy)/(densit * ulid **2 * objradius(n))

END SUBROUTINE fd_calc_surf_force

SUBROUTINE fd_calc_surf_nusselt_ave(n,n_called)

USE precision,              ONLY : r_single
USE real_parameters,        ONLY : one,zero,two,three,four,pi,ten
USE modfd_solve_linearsys,  ONLY : linsolve_gelem
USE shared_data,            ONLY : dtx,dty,dux,duy,dvx,dvy,surfpoint_cvx,nusseltpoint_cvy,nusseltpoint_interp,ulid,&
                                  nusseltinterpx,nusseltinterpy,xc,yc,nusseltds,nusseltnx,nusseltny,nusseltcentx,&
                                  nusseltcenty,objtp,th,objradius,objbradius,nnusseltpoints,li,t,t_locnusselt,naverage_steps,&
                                  nsphere

IMPLICIT NONE
INTEGER,INTENT(IN)                :: n !-- current cyl
INTEGER,INTENT(OUT)               :: n_called !-- current cyl

INTEGER                     :: i,j,k,myind(4),myind2(4,2)
REAL(KIND = r_single)       :: vandermonde(4,4),phiarr(4,1),coefts(4,1),detmt,dtdn
REAL(KIND = r_single)       :: tempi(2),temps,extr_fact,hh,aspratio,area
INTEGER,SAVE                :: n_entered = 0

n_entered = n_entered + 1
IF(n_entered == 1)THEN
  ALLOCATE(t_locnusselt(MAXVAL(nnusseltpoints(:),1),nsphere,naverage_steps))
  t_locnusselt = zero
ENDIF

DO i = 1,nnusseltpoints(n)-1
  !--Set the coeft matrix for this facet for each interpolation point in turn
  DO k = 1,2
   myind(1) = li(nusseltpoint_interp(1,k,i,n)) + nusseltpoint_interp(3,k,i,n)
   myind(2) = li(nusseltpoint_interp(2,k,i,n)) + nusseltpoint_interp(3,k,i,n)
   myind(3) = li(nusseltpoint_interp(1,k,i,n)) + nusseltpoint_interp(4,k,i,n)
   myind(4) = li(nusseltpoint_interp(2,k,i,n)) + nusseltpoint_interp(4,k,i,n)
   myind2(1,1) = nusseltpoint_interp(1,k,i,n)
   myind2(1,2) = nusseltpoint_interp(3,k,i,n)
   myind2(2,1) = nusseltpoint_interp(2,k,i,n)
   myind2(2,2) = nusseltpoint_interp(3,k,i,n)
   myind2(3,1) = nusseltpoint_interp(1,k,i,n)
   myind2(3,2) = nusseltpoint_interp(4,k,i,n)
   myind2(4,1) = nusseltpoint_interp(2,k,i,n)
   myind2(4,2) = nusseltpoint_interp(4,k,i,n)
   DO j = 1,4
    vandermonde(j,1) = xc(myind2(j,1))*yc(myind2(j,2))
    vandermonde(j,2) = xc(myind2(j,1))
    vandermonde(j,3) = yc(myind2(j,2))
    vandermonde(j,4) = one
   ENDDO
  
  !--set the rhs (temperature)
  DO j = 1,4 
    phiarr(j,1) = t(myind(j))
  ENDDO
                                                 
   CALL linsolve_gelem(vandermonde,phiarr,4,1,coefts,detmt)
   
   tempi(k) = coefts(1,1)*nusseltinterpx(k,i,n)*nusseltinterpy(k,i,n) + &
               coefts(2,1)*nusseltinterpx(k,i,n) + &
               coefts(3,1)*nusseltinterpy(k,i,n) + coefts(4,1)
  
  ENDDO
  
  !--Now extrapolate to the facet centre
  hh        = SQRT((nusseltinterpx(1,i,n) - nusseltinterpx(2,i,n))**2 + (nusseltinterpy(1,i,n) - nusseltinterpy(2,i,n))**2)
  extr_fact = SQRT((nusseltcentx(i,n) - nusseltinterpx(2,i,n))**2 + (nusseltcenty(i,n) - nusseltinterpy(2,i,n))**2) / hh
             

  temps = tempi(2) + extr_fact*(tempi(1) - tempi(2))

  !--second order forward approximation to the normal gradient at the corrent point
  dtdn  =(  -three*temps + four*tempi(1) - tempi(2) )/(two*hh)
  !dtdn = (tempi(2) - objtp(n))/(two*hh)
  t_locnusselt(i,n,n_entered) = - two*objradius(n)/(objtp(n) - th) * dtdn

ENDDO

 n_called = n_entered

END SUBROUTINE fd_calc_surf_nusselt_ave

SUBROUTINE fd_calc_surf_nusselt(n,locnusselt,avenusselt,npoints)

USE precision,              ONLY : r_single
USE real_parameters,        ONLY : one,zero,two,three,four,pi,ten
USE modfd_solve_linearsys,  ONLY : linsolve_gelem
USE shared_data,            ONLY : dtx,dty,dux,duy,dvx,dvy,surfpoint_cvx,nusseltpoint_cvy,nusseltpoint_interp,ulid,&
                                  nusseltinterpx,nusseltinterpy,xc,yc,nusseltds,nusseltnx,nusseltny,nusseltcentx,&
                                  nusseltcenty,objtp,th,objradius,objbradius,nnusseltpoints,li,t,t_locnusselt

IMPLICIT NONE
INTEGER,INTENT(IN)                :: n !-- current cyl
INTEGER,INTENT(IN),OPTIONAL       :: npoints !-- current cyl
REAL(KIND = r_single),INTENT(OUT) :: locnusselt(:),avenusselt !--local and average nusselt numbers for the current cyl

INTEGER                     :: i,j,k,myind(4),myind2(4,2),nn
REAL(KIND = r_single)       :: vandermonde(4,4),phiarr(4,1),coefts(4,1),detmt,dtdn
REAL(KIND = r_single)       :: tempi(2),temps,extr_fact,hh,aspratio,area

avenusselt = zero
IF(.NOT. PRESENT(npoints))THEN
  DO i = 1,nnusseltpoints(n)-1
    !--Set the coeft matrix for this facet for each interpolation point in turn
    DO k = 1,2
     myind(1) = li(nusseltpoint_interp(1,k,i,n)) + nusseltpoint_interp(3,k,i,n)
     myind(2) = li(nusseltpoint_interp(2,k,i,n)) + nusseltpoint_interp(3,k,i,n)
     myind(3) = li(nusseltpoint_interp(1,k,i,n)) + nusseltpoint_interp(4,k,i,n)
     myind(4) = li(nusseltpoint_interp(2,k,i,n)) + nusseltpoint_interp(4,k,i,n)
     myind2(1,1) = nusseltpoint_interp(1,k,i,n)
     myind2(1,2) = nusseltpoint_interp(3,k,i,n)
     myind2(2,1) = nusseltpoint_interp(2,k,i,n)
     myind2(2,2) = nusseltpoint_interp(3,k,i,n)
     myind2(3,1) = nusseltpoint_interp(1,k,i,n)
     myind2(3,2) = nusseltpoint_interp(4,k,i,n)
     myind2(4,1) = nusseltpoint_interp(2,k,i,n)
     myind2(4,2) = nusseltpoint_interp(4,k,i,n)
     DO j = 1,4
      vandermonde(j,1) = xc(myind2(j,1))*yc(myind2(j,2))
      vandermonde(j,2) = xc(myind2(j,1))
      vandermonde(j,3) = yc(myind2(j,2))
      vandermonde(j,4) = one
     ENDDO
  
    !--set the rhs (temperature)
    DO j = 1,4 
      phiarr(j,1) = t(myind(j))
    ENDDO
                                                 
     CALL linsolve_gelem(vandermonde,phiarr,4,1,coefts,detmt)
   
     tempi(k) = coefts(1,1)*nusseltinterpx(k,i,n)*nusseltinterpy(k,i,n) + &
                 coefts(2,1)*nusseltinterpx(k,i,n) + &
                 coefts(3,1)*nusseltinterpy(k,i,n) + coefts(4,1)
  
    ENDDO
  
    !--Now extrapolate to the facet centre
    hh        = SQRT((nusseltinterpx(1,i,n) - nusseltinterpx(2,i,n))**2 + (nusseltinterpy(1,i,n) - nusseltinterpy(2,i,n))**2)
    extr_fact = SQRT((nusseltcentx(i,n) - nusseltinterpx(2,i,n))**2 + (nusseltcenty(i,n) - nusseltinterpy(2,i,n))**2) / hh
             

    temps = tempi(2) + extr_fact*(tempi(1) - tempi(2))

    !--second order forward approximation to the normal gradient at the corrent point
    dtdn  =(  -three*temps + four*tempi(1) - tempi(2) )/(two*hh)
    !dtdn = (tempi(2) - objtp(n))/(two*hh)
    locnusselt(i) = - two*objradius(n)/(objtp(n) - th) * dtdn

    avenusselt = avenusselt + locnusselt(i)*nusseltds(i,n)

  ENDDO
  
ELSE !-- Do averaging on t_locnusselt for npoints
  DO i = 1,nnusseltpoints(n)-1
    locnusselt(i) = zero
    DO nn = 1,npoints
     locnusselt(i) = locnusselt(i) + t_locnusselt(i,n,nn)
    ENDDO
    locnusselt(i) = locnusselt(i)/npoints
    avenusselt = avenusselt + locnusselt(i)*nusseltds(i,n)
  ENDDO
  
ENDIF

!--For both Circular and elliptic cyl
aspratio = ((objradius(n)-objbradius(n))/(objradius(n)+objbradius(n)))**2
area = pi*(objradius(n)+objbradius(n))*(one + (three*aspratio)/(ten+SQRT(four - three*aspratio)))
!--tube bank
!avenusselt = avenusselt / area
!--single tube 

avenusselt = two * avenusselt / area

END SUBROUTINE fd_calc_surf_nusselt

SUBROUTINE fd_calc_wall_nusselt(walnusselt,tref,n)

USE real_parameters,  ONLY : zero,three,four,one
USE shared_data,      ONLY : njm,li,dtx,xc,objtp,x,y,t,fx,nj,nim
USE precision,        ONLY : r_single

IMPLICIT NONE

REAL(KIND = r_single),INTENT(INOUT) :: walnusselt(nj)
REAL(KIND = r_single),INTENT(IN)    :: tref
INTEGER                             :: n !--object to use
INTEGER               :: j,ij
REAL(KIND = r_single) :: dx,dtdx,tw

walnusselt(1) = zero
walnusselt(nj) = zero
DO j = 2,njm
  ij = li(nim) + j
  tw=t(ij)*fx(nim-1)+t(ij-nj)*(one-fx(nim-1))
  dx=x(nim)-x(nim-1) 
  dtdx = (three * tref - four * t(nim) + tw)/dx
  walnusselt(j) = (y(njm) - y(1))/(objtp(n) - tref) * dtdx
ENDDO

END SUBROUTINE fd_calc_wall_nusselt

SUBROUTINE fd_calc_lwall_nusselt(walnusselt,tw,tref)

USE real_parameters,  ONLY : zero,three,four,one
USE shared_data,      ONLY : njm,li,dtx,xc,x,y,t,fx,nj
USE precision,        ONLY : r_single

IMPLICIT NONE

REAL(KIND = r_single),INTENT(INOUT) :: walnusselt(nj)
REAL(KIND = r_single),INTENT(IN)    :: tw,tref

INTEGER               :: j,ij
REAL(KIND = r_single) :: dx,dtdx,tf

walnusselt(1) = zero
walnusselt(nj) = zero

DO j = 2,njm
  ij = li(2) + j
  tf=t(ij+nj)*fx(2)+t(ij)*(one-fx(2))
  dx=x(2)-x(1) 
  dtdx = -(three * tw - four * t(2) + tf)/dx
  walnusselt(j) = (y(njm) - y(1))/(tw - tref) * dtdx
ENDDO

END SUBROUTINE fd_calc_lwall_nusselt

SUBROUTINE fd_calc_lwall_nusselt_ave(walnusselt,ncalled,tw,tref)

USE real_parameters,  ONLY : zero,three,four,one
USE shared_data,      ONLY : njm,li,dtx,xc,x,y,t,fx,nj
USE precision,        ONLY : r_single

IMPLICIT NONE

REAL(KIND = r_single),INTENT(INOUT) :: walnusselt(nj)
REAL(KIND = r_single),INTENT(IN)    :: tw,tref
INTEGER,INTENT(INOUT)               :: ncalled

INTEGER               :: j,ij
REAL(KIND = r_single) :: dx,dtdx,tf
INTEGER,SAVE          :: n_entered = 0

n_entered = n_entered + 1

DO j = 2,njm
  ij = li(2) + j
  tf=t(ij+nj)*fx(2)+t(ij)*(one-fx(2))
  dx=x(2)-x(1) 
  dtdx = -(three * tw - four * t(2) + tf)/dx
  walnusselt(j) = walnusselt(j) + (y(njm) - y(1))/(tw - tref) * dtdx
ENDDO

ncalled = n_entered

END SUBROUTINE fd_calc_lwall_nusselt_ave

END MODULE modfd_calc_integrals