MODULE shared_data

USE parameters,     ONLY : nphi, max_char_len
USE precision,      ONLY : r_single

INTEGER     :: solver_type
INTEGER     :: use_GPU
LOGICAL     :: temp_visc !--Set to true in the problem setup for temperature dependant viscosity
!-------------------------Part-1: INTEGERS-----------------------------------------!
!-- ni, nj: number of grid points in x, y direction respectively (boundary included)
!-- nim, njm: last live cell in x, y direction
!-- nij: ni*nj
!-- li: conversion array from 2D to 1D
!-- imon, jmon monitor location for printing residuals
!-- ijmon: location of the monitor point in 1D format
!-- maxit,minit: maximum and minimum number of outer iteration
!-- iu,iv,ip,ien: a unique value for each solved variable ( 1 <  and <= nphi)
!-- ipr,jpr,ijpref: grid index of the pressure correction reference location, ijpref 1D index
!-- nsw: is the maximum allowed number of inner iterations for the ith variable
!-- itst: numebr of time steps to perform for unsteady calculations
!-- nprt: results are written to the file every nprt time steps
!-- nwrite: number of files to be written 
!-- itim: current tim step
INTEGER :: ni,nj,nim,njm,nij,imon,jmon,ijmon,maxit,minit,iu,iv,ip,ien,ipr,jpr,ijpref,itst,nprt,nwrite,itim
INTEGER,DIMENSION(:),ALLOCATABLE :: li(:),nsw(:)
!-------------------------Part-2: LOGICALS-----------------------------------------!
!-- lread set to true: results from previous run are read,before starting computation
!-- lwrite set to true: the results of calculation are written onto a file for post-
!--      processing or continuation in a later run
!-- ltest set to true, additional output is printed (global mass conservation
!--      check, convergence of inner iterations) 
!-- louts set to true, the initial field values of all variables are printed out
!-- loute set to true, the final values of all variables are printed out
!-- ltime set true, unsteady calculation is performed.
!-- laxis set to true, axi-symmetric geometry is assumed
!-- lcal(i) set to true, ith variable is solved for
LOGICAL :: initFromFile !--This just initializes the field from the output of a simple Eulerian run, 
                        !--but it is only an initialization not a restart.
                        !--having both lread = true and initFromFile = true, lread = true i.e it will 
                        !--overwrite the initialized values
LOGICAL :: lwrite,lread,ltest,laxis,louts,loute,ltime
LOGICAL,ALLOCATABLE,DIMENSION(:) :: lcal

!-------------------------Part-3: REALS-----------------------------------------!
!-- 3.1- solution control

!-- sor(i) is the required ratio of reduction of the residual norm
!--       during inner iterations for Ith variable before they are stoped 
!-- resor(i) is the residual of the ith variable
!-- urf(i) is the under-relaxation factor for the ith variable
!-- gamt is the blending factor for time differencing schemes (gamt=1 -> three time levels scheme,
!--       GAMT=0 -> Euler implicit scheme).
!-- beta is the volumetric expansion factor for the fluid
!-- gravx and gravy are the X and Y component of the gravity vector
!-- gds(I) is the blending factor for UDS and CDS in the equation for
!--       ith variable (convective terms; value 1.0 means CDS (second order), 
!--       0.0 means UDS (first order), any value between 0.0 and 1.0 can
!--       be used). The value 1.0 is recomended, except for coarse grids,
!--       in case convergence problems are encountered.
!-- sormax is the level of the residual norm at which
!--       outer iterations are stoped (converged solution)
!-- slarge is the level of the residual norm at which iterations are stoped because of divergence
!-- alfa is the parameter in the SIP solver ~ 0.9
!-- ulid is the lid velocity for the lid-driven flow, inlet velocity for duct
!-- TPER is the oscillation period in the case of unsteady flow with oscillating lid.
!-- om is the frequency of the oscillation for the oscillating lid case
REAL(KIND = r_single) :: gamt,beta,gravx,gravy,sormax,slarge,alfa,om,tper,ulid
REAL(KIND = r_single),DIMENSION(:),ALLOCATABLE :: sor,resor,urf,gds

!-- 3.2 flow variables

!-- densit: fluid density
!-- visc: fluid dynamic viscosity
!-- cpf: fluid mean heat capacity, kappaf: fluid thermal conductivity
!-- tref: reference temperature
!-- u: u-velocity field
!-- v: v-velocity field
!-- p: pressure field, pp: pressure correction
!-- t: temperature field. f1: mass flux e face. f2: mass flux n face = rho_e u_e A_e
!-- dpx,dpy: cell centre gradient 
!-- dux,duy,dvx,dvy: cell centre gradient 
!-- th,tc: hot and cold wall temperature (th input temp for flow problem)
!-- celbeta, celcp, celkappa: cell values for beta cp and kappa
!-- viscgamma: for temperature dependant viscosity, to be used in conjunction with temp_visc
REAL(KIND = r_single) :: densit,visc,tref,th,tc,cpf,kappaf,viscgamma
REAL(KIND = r_single),DIMENSION(:),ALLOCATABLE :: celbeta,celcp,celkappa,den,deno,u,v,p,pp,t,f1,f2,dpx,&
                                                  dpy,dux,duy,dvx,dvy,dtx,dty,lamvisc,ft1,ft2


!-- 3.3 geometry

!-- x,y: (x,y) of each grid point. 
!-- xc, yc: (xc,yc) of each grid center point. (2..nim, 2..njm) xc(1) = x(1), xc(ni)=x(nim) same for yc
!-- fx,fy: inerpolation factors (2..nim, 2..njm) (1,ni or nj) = 0.0
!-- r: radius = y(i) for axi-symmetric, 1.0 otherwise
REAL(KIND = r_single),DIMENSION(:),ALLOCATABLE :: x,y,xc,yc,fx,fy,r

!-- 3.4 time variables

!-- time: current time for unsteady calculation
!-- dt: time step size
!-- dtr = 1.0/dt
REAL(KIND = r_single) :: time, dt, dtr

!-- 3.5 old variables

!-- o -> previous time step
!-- oo -> two time sptes ago
REAL(KIND = r_single),DIMENSION(:),ALLOCATABLE ::uo,vo,to,uoo,voo,too

!-- 3.6 strings

!-- filenum: keeps the time step index on the end of the ascii files. 00000, 00001, 00002 etc
CHARACTER(LEN=5),SAVE, ALLOCATABLE     :: filenum(:)

!-- spherenum: keeps the sphere index on the end of the ascii files. 00000, 00001, 00002 etc
CHARACTER(LEN=5),SAVE, ALLOCATABLE     :: spherenum(:)

!--working arrays
REAL(KIND = r_single),DIMENSION(:),ALLOCATABLE :: ae,aw,an,as,ap,su,sv,apu,apv,apt

CHARACTER(LEN = max_char_len) :: problem_name
CHARACTER(LEN = 80)           :: title
INTEGER                       :: problem_len

REAL(KIND = r_single)         :: flomas,flomom
LOGICAL                       :: duct !--if true inlet outlet boundary
!-------------------------Part-4: FD ASPECT-----------------------------------------!

!--putobj: put a circle inside the domain or not 
!--objcentx,objcenty: centre of the circles
!--objradius: radius of the circles
!--objbradius:  other radius if ellipse used, should be equal to the objradius for circle
!--npoints: number of points to define the surfaces
!--surfpointx,surfpointy: x,y coordinates of the surface points
!--objcellx,objcelly: x,y coordinates of the material volumes
!--objcelldx,objcelldy: dx,dy of each material cv (mcv)
!--objcellvol: volume of MCVs
!--fdsu,fdsv: source terms for the momentum equation. objfx,objfy are the lag forces and fdsu,fdsv
!  are the corresponding eulerian forces.  
!--objpoint_cvx (y) : eulerian cv on which mcv resides
!--obju,objv: interpolated components of u,v to MCV centre
!--objfx,objfy: lagrangian force
!--objru,objrv: rigid velocity of the MCVs
!--objcentu,objcentv,objcentom: object center velocity (u,v) and angular vel
!--densitp: particle density (for now = densit)
!--fd_urf: fd force under relaxation
!--mcellpercv: number of material cv per cv
!--objtp: particle temperature
!--nsphere: number of spheres
!--filnprt: similar to nprt (for filament geommetry)
!--sphnprt: similar to nprt (for moving spheres)
!--nnusseltpoints: number ofnusseltpoints to calculate the local nusselt number 
!                  this is defined seperatly from nsurfacepoints since we what 
!                  half sphere in direction of the flow
!--calclocalnusselt: calculate local nusselt or not
!--objcellvertx,objcellverty: x and y coordinates of 4 vertices of each material cell (4,MAX(nobjcells),nsphere)
!--objpoint_interpx,objpoint_interpy: (1:2,MAX(nobjcells),nsphere) start-end of the interpolation (depends on delta width)
!--read_fd_geom: if true geometry is read other wise a circlular disc is created (stair-step grid) 
!--stationary:  set to TRUE for forced motion or stationary objects
!--forcedmotion:  set to TRUE if it is forced motion and change the forced_motion function for the type of motion
!--dxmean:      mean dx of the domain calculated in fd_create_geom
!--dxmeanmoved: mean motion of the particles in x-dir since the last mesh motion
!--cpp:  mean heat capacity particle, kappap: particle conductivity 
!--movingmesh:  does not work properly, for simulation of infinitly large duct
!--isotherm:   set to false if particle has source term (variable temp), If true objqp should be set and objtp is the initial temp
!--ibsu,ibsv: momentum sources similar to fdsu,fdsv
!--naverage_steps: number of time steps to be included in time averaging (it is done on the last naverage_steps steps)
!--calcwalnusselt: calculate the local nusselt number on the left wall.
REAL(KIND = r_single),ALLOCATABLE,DIMENSION(:)     :: ibsu,ibsv
INTEGER                       :: nsphere,nfil,filnprt,sphnprt,naverage_steps
LOGICAL                       :: putobj,calcsurfforce,calclocalnusselt,read_fd_geom,stationary,movingmesh,&
                                 forcedmotion,isotherm,calclocalnusselt_ave,calcwalnusselt
REAL(KIND = r_single)         :: fd_urf,filgravx,filgravy,filfr,filalpha,filbeta,filgamma
REAL(KIND = r_single)         :: dxmean,dxmeanmoved,dxmeanmovedtot

!--Current particle velocity for the forced motion case
REAL(KIND = r_single),ALLOCATABLE,DIMENSION(:)     :: up,vp,omp,objcentxinit,objcentyinit
REAL(KIND = r_single),ALLOCATABLE,DIMENSION(:,:)   :: objcellxinit,objcellyinit,surfpointxinit,surfpointyinit
REAL(KIND = r_single),ALLOCATABLE,DIMENSION(:,:,:) :: objcellvertxinit,objcellvertyinit
!--allocatable for 1..nsphere
INTEGER,DIMENSION(:),ALLOCATABLE :: nsurfpoints,nobjcells,nnusseltpoints
REAL(KIND = r_single),DIMENSION(:),ALLOCATABLE   :: objqp,objtp,densitp,betap,objcentx,objcenty,objradius,objbradius,&
                                                    objcentu,objcentv,objcentom,objcentmi,objvol,&
                                                    objcentxo,objcentyo,objcentvo,objcentuo,objcentomo,&
                                                    mcellpercv,cpp,kappap
REAL(KIND = r_single),DIMENSION(:,:),ALLOCATABLE :: objcento
!--allocatable for 1..Max(nobjcells),1..nsphere
REAL(KIND = r_single),DIMENSION(:,:),ALLOCATABLE :: surfpointx,surfpointy,objcellx,objcelly,&
                                                    objcelldx,objcelldy,objcellvol,obju,objv,&
                                                    objfx,objfy,objru,objrv,objt,objrt,objq,&
                                                    objapu,objapv,objapt,nusseltpointx,nusseltpointy,&
                                                    localnusselt
INTEGER,DIMENSION(:,:),ALLOCATABLE               :: objcell_overlap
!--Particle Collision forces, 
!--Fpq, sum of forces on particle p due to all other particles q /= p, 
!--Fpw, sum of forces on particle p due to all walls  
REAL(KIND = r_single),DIMENSION(:,:),ALLOCATABLE :: Fpq,Fpw

REAL(KIND = r_single),DIMENSION(:,:,:),ALLOCATABLE :: objcellvertx,objcellverty,t_locnusselt
INTEGER,DIMENSION(:,:),ALLOCATABLE               :: objpoint_cvx,objpoint_cvy
INTEGER,DIMENSION(:,:,:),ALLOCATABLE             :: objpoint_interpx,objpoint_interpy

!--Allocatable for 1..nfil; These are similar to sphere counterparts
!--densitfil is the difference between the densities of filament and fluid
!--filbenrig: is the filament bending rigidity 
!--filsb,filsy,filsx: used to keep the rhs for the tension eq and position eq (x,y equations).
!--filinitposx,*y: are the initial positions of the first node (fixed)
!--filtheta: angle of the last point meassured counterclockwise with respect to x-axis (in radians)
!--filintegfx,filintegfy: integral of the velocity error
!--filcintegfx,filcintegfy: current time step contribution to the integral of the error, used to correct the integral
INTEGER                          :: deltalen 
INTEGER,DIMENSION(:),ALLOCATABLE :: nfilpoints
REAL(KIND = r_single),DIMENSION(:),ALLOCATABLE   :: densitfil,filap,filaw,filae,filst,filbenrig,filsx,filsy,fillen,fillenc,&
                                                    filfirstposx,filfirstposy,fillasttheta
!--Allocatable for 1..Max(nfilpoints),1..nfil
REAL(KIND = r_single),DIMENSION(:,:),ALLOCATABLE :: filpointx,filpointy,filpointyo,filpointxo,&
                                                    filten,filds,filpointxpen,filpointypen,filu,filv,&
                                                    filru,filrv,filfx,filfy,filintegfx,filintegfy,&
                                                    filcintegfx,filcintegfy

INTEGER,DIMENSION(:,:),ALLOCATABLE               :: filpoint_cvx,filpoint_cvy
INTEGER,DIMENSION(:,:,:),ALLOCATABLE             :: filpoint_interpx,filpoint_interpy

!==allocatable for 1..ncv (Eulerian) (no everlaps between objects), used both for VB and FD
!--these are body forces, s means source, then variable, b is averaged c is correction
!-fdFc collision force
REAL(KIND = r_single),DIMENSION(:),ALLOCATABLE :: fdsu,fdsv,fdsuc,fdsvc,fdsub,fdsvb,fdst,fdstc,FdFcu,FdFcv

!--Surface property calculation arrays (1..Max(nsurfpoint),1..nsphere
REAL(KIND = r_single),DIMENSION(:,:),ALLOCATABLE   :: surfds,surfnx,surfny,surfcentrex,surfcentrey,nusseltds,nusseltnx,nusseltny,&
                                                      nusseltcentx,nusseltcenty
!--location of the intepolation points 1..2,1..Max(nsurfpoint),1..nsphere
REAL(KIND = r_single),DIMENSION(:,:,:),ALLOCATABLE :: surfinterpx,surfinterpy,nusseltinterpx,nusseltinterpy 
!--1..3, 1..Max(nsurfpoint),1..nsphere (1,:,:) is the location of the surface point 2 is the first interpolation point, 
!--   is the third interpolation point. 
INTEGER,DIMENSION(:,:,:),ALLOCATABLE               :: surfpoint_cvx,surfpoint_cvy

!--1..3, 1..Max(nnusseltpoint),1..nsphere (1,:,:) is the location of the surface point 2 is the first interpolation point, 
!--   is the third interpolation point. These are defined for easier calculation otherwise surface points could be used
!--   for nusselt calculations
INTEGER,DIMENSION(:,:,:),ALLOCATABLE               :: nusseltpoint_cvx,nusseltpoint_cvy

!--1..4, 1..2, 1..Max(nsurfpoint),1..nsphere (1,:,:) are the index of four cells around each intepolation points. 
INTEGER,DIMENSION(:,:,:,:),ALLOCATABLE               :: surfpoint_interp,nusseltpoint_interp

!--sparsekit and hypre variables
 DOUBLE PRECISION, ALLOCATABLE :: Alu(:)
 INTEGER, ALLOCATABLE :: Jlu(:), Ju(:), Jw(:)

 DOUBLE PRECISION, ALLOCATABLE :: Work(:,:)

 DOUBLE PRECISION, ALLOCATABLE :: Acoo(:)   ! A-matrix in COO-format
 INTEGER, ALLOCATABLE          :: Arow(:)   ! A-matrix row entries
 INTEGER, ALLOCATABLE          :: Acol(:)   ! A-matrix column entries

 DOUBLE PRECISION, ALLOCATABLE :: Acsr(:)   ! A-matrix in CSR-format
 INTEGER, ALLOCATABLE          :: Arwc(:)   ! A-matrix row entries
 INTEGER, ALLOCATABLE          :: Aclc(:)   ! A-matrix column entries

 DOUBLE PRECISION, ALLOCATABLE :: RHS(:)    ! righthand side
 DOUBLE PRECISION, ALLOCATABLE :: SOL(:)    ! solution

 INTEGER                       :: NNZ,Ncel
 INTEGER,ALLOCATABLE           :: lli(:)

 INTEGER*8                     :: hypre_A,Hypre_b,Hypre_x,mpi_comm !--hypre pointers to the linear system
 INTEGER                       :: myrank_mpi,nprocs_mpi 

!--Some Auxiliary arrays 
INTEGER,ALLOCATABLE            :: celltype(:)

END MODULE shared_data