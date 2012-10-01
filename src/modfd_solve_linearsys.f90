MODULE modfd_solve_linearsys

PRIVATE
PUBLIC :: fd_solve_sip2d,linsolve_gelem,fil_sol_tdma2,fil_sol_tdma,fd_spkit_interface,&
          calc_residual,copy_solution,fd_solve_cgs2d,fd_cooInd_create,fd_cooVal_create
CONTAINS
!=====================================
SUBROUTINE fd_solve_sip2d(ae,an,aw,as,ap,su,phi,li,ni,nj,nsweep,resor,sor)

USE precision,      ONLY : r_single
USE parameters,     ONLY : out_unit
USE shared_data,    ONLY : alfa,ltest
USE real_parameters,ONLY : vsmall,one,zero

IMPLICIT NONE

!--Inputs
INTEGER,INTENT(IN)                              :: ni,nj,nsweep
INTEGER,DIMENSION(:),INTENT(IN)                 :: li
REAL(KIND = r_single),DIMENSION(:),INTENT(IN)   :: ae,an,aw,as,ap,su
REAL(KIND = r_single),DIMENSION(:),INTENT(INOUT):: phi
REAL(KIND = r_single)             ,INTENT(IN)   :: sor
REAL(KIND = r_single)             ,INTENT(INOUT):: resor
!--Internals
REAL(KIND = r_single),DIMENSION(:),ALLOCATABLE  :: lw,ls,un,ue,lpr,res
REAL(KIND = r_single)                           :: p1,p2,resl,rsm
INTEGER                                         :: ierror,nim,njm,i,ij,L

ALLOCATE(lw(ni*nj),ls(ni*nj),un(ni*nj),ue(ni*nj),lpr(ni*nj),res(ni*nj),STAT = ierror)

IF(ierror /= 0)THEN
  WRITE(out_unit,*)'fd_solve_sip2d: not enough memory'
  STOP
ENDIF

lw = zero
ls = zero
un = zero
ue = zero
lpr= zero
res= zero
nim = ni-1
njm = nj-1

DO i=2,nim
  DO ij=li(i)+2,li(i)+njm
    lw(ij)=aw(ij)/(one+alfa*un(ij-nj))
    ls(ij)=as(ij)/(one+alfa*ue(ij-1))
    p1=alfa*lw(ij)*un(ij-nj)
    p2=alfa*ls(ij)*ue(ij-1)
    lpr(ij)=one/(ap(ij)+p1+p2-lw(ij)*ue(ij-nj)-ls(ij)*un(ij-1))
    un(ij)=(an(ij)-p1)*lpr(ij)
    ue(ij)=(ae(ij)-p2)*lpr(ij)
  ENDDO
ENDDO

DO L=1,nsweep
  resl=zero
      
  !--CALCULATE RESIDUAL AND OVERWRITE IT BY INTERMEDIATE VECTOR
  DO i=2,nim
    DO ij=li(i)+2,li(i)+njm
      res(ij)=su(ij)-an(ij)*phi(ij+1)-as(ij)*phi(ij-1)-&
              ae(ij)*phi(ij+nj)-aw(ij)*phi(ij-nj)-ap(ij)*phi(ij)
      resl=resl+abs(res(ij))
      res(ij)=(res(ij)-ls(ij)*res(ij-1)-lw(ij)*res(ij-nj))*lpr(ij)
    END DO
  END DO
  !--STORE INITIAL RESIDUAL SUM FOR CHECKING CONV. OF OUTER ITER.
  
  IF(L == 1) resor=resl
  rsm=resl/(resor+vsmall)

  !--BACK SUBSTITUTION AND CORRECTION

  DO i=nim,2,-1
    DO ij=li(i)+njm,li(i)+2,-1
      res(ij)=res(ij)-un(ij)*res(ij+1)-ue(ij)*res(ij+nj)
      phi(ij)=phi(ij)+res(ij)
    END DO
  END DO

  !--CHECK CONVERGENCE OF INNER ITERATIONS

  IF(ltest) WRITE(out_unit,*) '  ',L,'INNER ITER, RESL = ',RESL
  IF(rsm < sor)THEN
    DEALLOCATE(lw,ls,un,ue,lpr,res)
    RETURN
  ENDIF
END DO

DEALLOCATE(lw,ls,un,ue,lpr,res)

END SUBROUTINE fd_solve_sip2d

SUBROUTINE fd_solve_cgs2d(ae,an,aw,as,ap,su,phi,li,ni,nj,nsweep,resor,sor)

USE precision,      ONLY : r_single
USE parameters,     ONLY : out_unit
USE shared_data,    ONLY : alfa,ltest
USE real_parameters,ONLY : vsmall,one,zero,vlarge

IMPLICIT NONE

!--Inputs
INTEGER,INTENT(IN)                              :: ni,nj,nsweep
INTEGER,DIMENSION(:),INTENT(IN)                 :: li
REAL(KIND = r_single),DIMENSION(:),INTENT(IN)   :: ae,an,aw,as,ap,su
REAL(KIND = r_single),DIMENSION(:),INTENT(INOUT):: phi
REAL(KIND = r_single)             ,INTENT(IN)   :: sor
REAL(KIND = r_single)             ,INTENT(INOUT):: resor
!--Locals
REAL(KIND = r_single) :: res0,s0,sk,bet,s2,resn,alf,rsm
REAL(KIND = r_single) :: pk(ni*nj),zk(ni*nj),d(ni*nj),res(ni*nj)
INTEGER               :: i,ij,n,nim,njm

pk = zero
zk = zero
d = zero
res = zero
nim = ni-1
njm = nj-1

 res0=zero
 DO i=2,nim
   DO ij=li(i)+2,li(i)+njm
     res(ij)=su(ij)-ap(ij)*phi(ij)-ae(ij)*phi(ij+nj)-&
          aw(ij)*phi(ij-nj)-an(ij)*phi(ij+1)-as(ij)*phi(ij-1)
     res0=res0+abs(res(ij))
   END DO
 END DO

 resor = res0
 !--Preconditioned diagonal
 DO i=2,nim
   DO ij=li(i)+2,li(i)+njm
     d(ij)=one/(ap(ij)-d(ij-nj)*aw(ij)**2-d(ij-1)*as(ij)**2)
   END DO
 END DO

 s0=vlarge
 DO n=1,nsweep

   DO i=2,nim
     DO ij=li(i)+2,li(i)+njm
       zk(ij)=(res(ij)-aw(ij)*zk(ij-nj)-as(ij)*zk(ij-1))*d(ij)
     END DO
   END DO

   DO i=2,nim
     DO ij=li(i)+2,li(i)+njm
       zk(ij)=zk(ij)/(d(ij)+vsmall)
     END DO
   END DO

   sk=zero
   DO i=nim,2,-1
     DO ij=li(i)+njm,li(i)+2,-1
       zk(ij)=(zk(ij)-ae(ij)*zk(ij+nj)-an(ij)*zk(ij+1))*d(ij)
       sk=sk+res(ij)*zk(ij) 
     END DO
   END DO

   bet=sk/s0
   DO i=2,nim
     DO ij=li(i)+2,li(i)+njm
       pk(ij)=zk(ij)+bet*pk(ij)
     END DO
   END DO

   s2=zero
   DO i=2,nim
     DO ij=li(i)+2,li(i)+njm
       zk(ij)=ap(ij)*pk(ij)+ae(ij)*pk(ij+nj)+aw(ij)*pk(ij-nj)+ &
              an(ij)*pk(ij+1)+as(ij)*pk(ij-1) 
       s2=s2+pk(ij)*zk(ij)
     END DO
   END DO

   alf=sk/s2
   resn=zero
   DO i=2,nim
     DO ij=li(i)+2,li(i)+njm
       phi(ij)=phi(ij)+alf*pk(ij) 
       res(ij)=res(ij)-alf*zk(ij)
       resn=resn+ABS(res(ij))
     END DO
   END DO
   
   s0=sk

   rsm=resn/(res0+vsmall)
   IF(ltest) WRITE(out_unit,*) '  ',N,'INNER ITER, RESL = ',RSM

   IF(rsm < sor) RETURN 

 END DO

END SUBROUTINE fd_solve_cgs2d

SUBROUTINE linsolve_gelem(a,b,row,soln,ans,detmt)

!--A naive direct gaussian elemination for very small systems

! INPUT: a(row,row) -> Matrix of coefficients of system.
! INPUT: b(row,soln) -> RHS Vectors(matrix)
!   row = # of equations,      soln = # of solution vectors (multiple rhs)
! OUTPUT:ans(row,soln) -> #soln unknown(answer) vectors
! OUTPUT:detmt is determinanat

USE precision,      ONLY : r_single
USE real_parameters,ONLY : one,zero,small,minusone
IMPLICIT NONE
INTEGER, INTENT(IN) :: row,soln
REAL(KIND = r_single) , INTENT(IN), TARGET :: a(row,row), b(row,soln)
REAL(KIND = r_single) , INTENT(OUT) :: ans(row,soln)
REAL(KIND = r_single) , INTENT(OUT) :: detmt

INTEGER i,j, k(1),p  ! need to specify p with dimension one to
                     ! use MAXLOC intrinsic function
INTEGER :: imm, memoryvec(row), iw,ix,iy,iz, itmp
REAL(KIND = r_single) , POINTER    :: swap_ij(:), c(:,:)
REAL(KIND = r_single)              :: tmp
REAL(KIND = r_single), ALLOCATABLE :: ctmp(:,:)
REAL(KIND = r_single),PARAMETER    :: EPS=3.0e-15

ALLOCATE (c(row,row+soln), swap_ij(row+soln))
ALLOCATE (ctmp(row,row+soln))
detmt = 1.0d0

! ERROR Check #1
! Checking for rows of zeros (Theorem 1)
DO i = 1, row
  tmp = zero
  DO j = 1, row 
      tmp = tmp + ABS( a(i,j) )
  END DO
  IF (tmp < eps*eps) THEN
     PRINT *, 'Error the sum of row ',i,' is less than ',eps*eps
     PRINT *, 'High possibility that matrix is singular!!!'
     RETURN
  END IF
END DO

! Constructing ONE single augmented matrix (actually c is a pointer)
DO i = 1, row
  c(i, 1:row) = a(i,:)
  c(i, row+1 : row+soln) = b(i,:)
END DO

! Initializing memory vector for the index/position for solution 
! Must be used if there is column swapping
DO i = 1, row
  memoryvec(i) = i
END DO


! Begin Row-Reduction Operations
DO j = 1, row                 ! for the j-th column

!  PIVOTING - of the submatrix for the j-th operation
!  i.e. for every row, normalize the entire row by the element with the largest magnitude
!  ONLY do so provided that largest value in a row is < 1.0
  DO i = j, row
    tmp = MAXVAL( ABS(c(i,j:row)) )
    IF ( (tmp<one) .AND. (tmp>eps) ) THEN
      DO p = j, row+soln  ! the preceeding columns in this row are ZEROS anyway!
        c(i,p) = c(i,p) / tmp
      END DO
      detmt = detmt * tmp
    END IF
  END DO         ! i loop


  ! Finding location of max value along row j
  k = j - 1 + MAXLOC( ABS(c(j , j:row)) )  ! assigning to the k-array
  imm = k(1)  ! Swap row/col  imm, j

  IF (imm /= j) THEN
    DO iw = 1, row
      IF (iw==imm) THEN
        ix=j
      ELSEIF (iw==j) THEN
        ix=imm
      ELSE
        ix=iw
      ENDIF
      DO iy = 1, row
        IF (iy==imm) THEN
          iz=j
        ELSEIF (iy==j) THEN
          iz=imm
        ELSE
          iz=iy
        ENDIF
        ctmp(ix,iz) = c(iw,iy)     ! Building up temporary matrix to represent
      END DO     ! END loop iy      ! the matrix where the row/col is changed imm<->i
      DO iy = row+1 , row+soln
        ctmp(ix,iy) = c(iw,iy)
      END DO
    END DO  ! END for iw
    itmp = memoryvec(j)  ! Keep track of swapped variables/solutions
    memoryvec(j) = memoryvec(imm)
    memoryvec(imm) = itmp
! Actual swapping of row/col imm<->j
    c(:,:) = ctmp(:,:)
  ENDIF


! Find row p that has the largest element in column j
! among the elements below the diagonal element.

!  Pivotal strategy - SWAP rows to make pivotal element c[j][j] have the
!  greatest magnitude in its column. This prevents unnecessary division by
!  a small number.
  k = j - 1 + MAXLOC( ABS(c(j:row , j)) )        ! assigning to the k-array
  p = k(1)

 IF (p /= j) THEN
    swap_ij(:) = c(j,:)
    c(j,:) = c(p,:)
    c(p,:) = swap_ij(:)
    detmt = minusone * detmt
 ENDIF
! Determinant change signs if rows are swapped



! ERROR Check #2
! If after swapping rows the diagonal element(now having the largest value)
! is still very small then Matrix is singular
 IF ( (ABS(c(j,j)) < small*eps).AND.(j/=row)) THEN
  PRINT *, 'ERROR: Matrix is Singular. Found at ',j,'-th row operation'
  PRINT *, 'diagonal element is ', c(j,j)
  RETURN
 ENDIF

! ERROR Check #3
! If at the j-th row, all other elements are zero and element(j,j) is very small
! THEN might be singular or inconsistent matrix
  tmp = zero
  DO i = 1, j-1
    tmp = tmp + ABS(c(j,i))
  END DO
  DO i = j+1, row
    tmp = tmp + ABS(c(j,i))
  END DO                 
 IF ( (tmp<eps) .AND. (ABS(c(j,j))<100_r_single*eps) ) THEN
    PRINT *, 'ERROR: Matrix is Singular/Inconsistent. Found at ',j,'-th row operation'
    PRINT *, 'Diagonal element is too small ', c(j,j)
    PRINT *, 'And all other elements in this row are almost zero'
    RETURN
  ENDIF

! Divide the j-th row by leading diagonal
  tmp = c(j,j)
  c(j,j) = one
  DO i = j+1, row+soln
    c(j,i) = c(j,i) / tmp
  END DO  

! Finding Determinant as the factors of the diagoanals
  detmt = detmt * tmp

! Subtract multiple of j-th row from all rows (except j-th row itself)
! This leaves the j-th column with only one "1" in the diagonal position.
  DO i = 1, row                                      ! 1 0 ~ ~
   IF (i /= j) THEN                                  ! 0 1 ~ ~
      tmp = c(i,j)                                   ! 0 0 1 ~
      c(i,j) = zero                                  ! 0 0 ~ ~
      DO p = j+1, row+soln
        c(i,p) = c(i,p) - tmp * c(j,p)
      END DO
    ENDIF
  END DO


END DO      ! j loop

ans(memoryvec(:),:) = c(:,row+1:row+soln)

DEALLOCATE(c, ctmp, swap_ij)

END subroutine linsolve_gelem

SUBROUTINE fil_sol_tdma(x,aw,ap,ae,b,nx)

USE precision,          ONLY : r_single

IMPLICIT NONE

REAL(KIND = r_single), DIMENSION(:), INTENT(inout) :: b,aw,ap,ae,x
INTEGER,                             INTENT(IN)    :: nx
INTEGER :: i

DO i = 2,nx
  ap(i) = ap(i) - aw(i)/ap(i-1)*ae(i-1)
  b(i) = b(i) - aw(i)/ap(i-1)*b(i-1)
EnDDO

x(nx) = b(nx)/ap(nx)

DO i = nx-1,1,-1
  x(i) = (b(i) - ae(i)*x(i+1))/ap(i)
ENDDO

END SUBROUTINE fil_sol_tdma

SUBROUTINE fil_sol_tdma2(x,y,aw,ap,ae,bx,by,nx)

USE precision,          ONLY : r_single

IMPLICIT NONE

REAL(KIND = r_single), DIMENSION(:), INTENT(inout) :: bx,by,aw,ap,ae,x,y
INTEGER,                             INTENT(IN)    :: nx
INTEGER :: i

DO i = 2,nx
  ap(i) = ap(i) - aw(i)/ap(i-1)*ae(i-1)
  bx(i) = bx(i) - aw(i)/ap(i-1)*bx(i-1)
  by(i) = by(i) - aw(i)/ap(i-1)*by(i-1)
EnDDO

x(nx) = bx(nx)/ap(nx)
y(nx) = by(nx)/ap(nx)

DO i = nx-1,1,-1
  x(i) = (bx(i) - ae(i)*x(i+1))/ap(i)
  y(i) = (by(i) - ae(i)*y(i+1))/ap(i)
ENDDO

END SUBROUTINE fil_sol_tdma2

SUBROUTINE fd_spkit_interface(ap,as,an,aw,ae,su,phi)
  
USE shared_data,      ONLY : Acoo,Arow,Acol,nj,nim,njm,lli,li,nij,Ncel,SOL,RHS,NNZ
USE precision,        ONLY : r_single
IMPLICIT NONE

REAL(KIND = r_single) :: ap(nij),as(nij),an(nij),aw(nij),ae(nij),su(nij),phi(nij)
INTEGER :: ia,i,j,ijn,ije,ijw,ijs,ij

  ia = 0
  DO i = 2,nim
    IF( i == 2)THEN
      DO j = 2,njm
        ij = li(i) + j; ije=ij+nj; ijn=ij+1; ijw=ij-nj; ijs=ij-1
        rhs(lli(ij)) = DBLE( su(ij) )
        sol(lli(ij)) = DBLE( phi(ij))
        IF(j == 2)THEN
          ia = ia + 1
          Acoo(ia) = DBLE( ae(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ije)
          ia = ia + 1
          Acoo(ia) = DBLE( an(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ijn)
          ia = ia + 1
          Acoo(ia) = DBLE( ap(ij) ); Arow(ia) = lli(ij); Acol(ia) = Arow(ia)
        ELSEIF(j == njm)THEN
          ia = ia + 1
          Acoo(ia) = DBLE( ae(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ije)
          ia = ia + 1
          Acoo(ia) = DBLE( as(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ijs)
          ia = ia + 1
          Acoo(ia) = DBLE( ap(ij) ); Arow(ia) = lli(ij); Acol(ia) = Arow(ia)
        ELSE
          ia = ia + 1
          Acoo(ia) = DBLE( ae(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ije)
          ia = ia + 1
          Acoo(ia) = DBLE( an(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ijn)
          ia = ia + 1
          Acoo(ia) = DBLE( as(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ijs)
          ia = ia + 1
          Acoo(ia) = DBLE( ap(ij) ); Arow(ia) = lli(ij); Acol(ia) = Arow(ia)
        ENDIF
      ENDDO
    ELSEIF(i == nim)THEN
       DO j = 2,njm
        ij = li(i) + j; ije=ij+nj; ijn=ij+1; ijw=ij-nj; ijs=ij-1
        rhs(lli(ij)) = DBLE( su(ij))
        sol(lli(ij)) = DBLE(phi(ij))
        IF(j == 2)THEN
          ia = ia + 1
          Acoo(ia) = DBLE( an(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ijn)
          ia = ia + 1
          Acoo(ia) = DBLE( aw(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ijw)
          ia = ia + 1
          Acoo(ia) = DBLE( ap(ij) ); Arow(ia) = lli(ij); Acol(ia) = Arow(ia)
        ELSEIF(j == njm)THEN
          ia = ia + 1
          Acoo(ia) = DBLE( aw(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ijw)
          ia = ia + 1
          Acoo(ia) = DBLE( as(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ijs)
          ia = ia + 1
          Acoo(ia) = DBLE( ap(ij) ); Arow(ia) = lli(ij); Acol(ia) = Arow(ia)
        ELSE
          ia = ia + 1
          Acoo(ia) = DBLE( an(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ijn)
          ia = ia + 1
          Acoo(ia) = DBLE( aw(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ijw)
          ia = ia + 1
          Acoo(ia) = DBLE( as(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ijs)
          ia = ia + 1
          Acoo(ia) = DBLE( ap(ij) ); Arow(ia) = lli(ij); Acol(ia) = Arow(ia)
        ENDIF
      ENDDO                       
    ELSE
       DO j = 2,njm
        ij = li(i) + j; ije=ij+nj; ijn=ij+1; ijw=ij-nj; ijs=ij-1
        rhs(lli(ij)) = DBLE( su(ij) )
        sol(lli(ij)) = DBLE( phi(ij))
        IF(j == 2)THEN
          ia = ia + 1
          Acoo(ia) = DBLE( ae(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ije)
          ia = ia + 1
          Acoo(ia) = DBLE( an(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ijn)
          ia = ia + 1
          Acoo(ia) = DBLE( aw(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ijw)
          ia = ia + 1
          Acoo(ia) = DBLE( ap(ij) ); Arow(ia) = lli(ij); Acol(ia) = Arow(ia)
        ELSEIF(j == njm)THEN
          ia = ia + 1
          Acoo(ia) = DBLE( ae(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ije)
          ia = ia + 1
          Acoo(ia) = DBLE( aw(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ijw)
          ia = ia + 1
          Acoo(ia) = DBLE( as(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ijs)
          ia = ia + 1
          Acoo(ia) = DBLE( ap(ij) ); Arow(ia) = lli(ij); Acol(ia) = Arow(ia)
        ELSE
          ia = ia + 1
          Acoo(ia) = DBLE( ae(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ije)
          ia = ia + 1
          Acoo(ia) = DBLE( an(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ijn)
          ia = ia + 1
          Acoo(ia) = DBLE( aw(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ijw)
          ia = ia + 1
          Acoo(ia) = DBLE( as(ij) ); Arow(ia) = lli(ij); Acol(ia) = lli(ijs)
          ia = ia + 1
          Acoo(ia) = DBLE( ap(ij) ); Arow(ia) = lli(ij); Acol(ia) = Arow(ia)
        ENDIF
      ENDDO 
    ENDIF
  ENDDO

  IF(ia /= nnz)WRITE(*,*)'fd_spkit_interface:: ia = ',ia,' ,NNZ = ',NNZ
END SUBROUTINE fd_spkit_interface

SUBROUTINE calc_residual(ap,as,an,aw,ae,su,phi,resor)

USE shared_data,      ONLY : ni,nj,nim,njm,li,nij
USE precision,        ONLY : r_single
USE real_parameters,  ONLY : zero
IMPLICIT NONE

REAL(KIND = r_single),INTENT(IN)    :: ap(nij),as(nij),an(nij),aw(nij),ae(nij),su(nij),phi(nij)
REAL(KIND = r_single),INTENT(INOUT) :: resor

REAL(KIND = r_single) :: res(nij)
INTEGER               :: ij,i

  res = zero
  DO i=2,nim
    DO ij=li(i)+2,li(i)+njm
      res(ij)=su(ij)-an(ij)*phi(ij+1)-as(ij)*phi(ij-1)-&
              ae(ij)*phi(ij+nj)-aw(ij)*phi(ij-nj)-ap(ij)*phi(ij)
    END DO
  END DO

  resor = SUM(ABS(Res)) 

ENDSUBROUTINE calc_residual

SUBROUTINE copy_solution(a,sol)

USE shared_data,    ONLY : NCel,nij,nim,njm,li
USE precision,      ONLY : r_single

IMPLICIT NONE

REAL(KIND = r_single) :: a(nij)
DOUBLE PRECISION      :: sol(NCel)

INTEGER :: ic,i,j
 
  ic = 0
  DO i= 2,nim
    DO j= 2,njm
      ic = ic + 1
      a(li(i) + j) = SOL(ic)
    ENDDO
  ENDDO

END SUBROUTINE copy_solution

SUBROUTINE fd_cooInd_create()
  
USE shared_data,      ONLY : Arow,Acol,nj,nim,njm,lli,li,nij,Ncel,SOL,RHS,NNZ
USE precision,        ONLY : r_single
IMPLICIT NONE

INTEGER :: ia,i,j,ijn,ije,ijw,ijs,ij

  ia = 0
  DO i = 2,nim
    IF( i == 2)THEN
      DO j = 2,njm
        ij = li(i) + j; ije=ij+nj; ijn=ij+1; ijw=ij-nj; ijs=ij-1
        IF(j == 2)THEN
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ije)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ijn)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = Arow(ia)
        ELSEIF(j == njm)THEN
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ije)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ijs)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = Arow(ia)
        ELSE
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ije)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ijn)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ijs)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = Arow(ia)
        ENDIF
      ENDDO
    ELSEIF(i == nim)THEN
       DO j = 2,njm
        ij = li(i) + j; ije=ij+nj; ijn=ij+1; ijw=ij-nj; ijs=ij-1
        IF(j == 2)THEN
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ijn)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ijw)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = Arow(ia)
        ELSEIF(j == njm)THEN
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ijw)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ijs)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = Arow(ia)
        ELSE
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ijn)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ijw)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ijs)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = Arow(ia)
        ENDIF
      ENDDO                       
    ELSE
       DO j = 2,njm
        ij = li(i) + j; ije=ij+nj; ijn=ij+1; ijw=ij-nj; ijs=ij-1
        IF(j == 2)THEN
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ije)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ijn)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ijw)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = Arow(ia)
        ELSEIF(j == njm)THEN
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ije)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ijw)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ijs)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = Arow(ia)
        ELSE
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ije)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ijn)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ijw)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = lli(ijs)
          ia = ia + 1
          Arow(ia) = lli(ij); Acol(ia) = Arow(ia)
        ENDIF
      ENDDO 
    ENDIF
  ENDDO

  IF(ia /= nnz)WRITE(*,*)'fd_spkit_interface:: ia = ',ia,' ,NNZ = ',NNZ

END SUBROUTINE fd_cooInd_create

SUBROUTINE fd_cooVal_create(ap,as,an,aw,ae,su,phi)
  
USE shared_data,      ONLY : Acoo,nj,nim,njm,lli,li,nij,Ncel,SOL,RHS,NNZ
USE precision,        ONLY : r_single
IMPLICIT NONE

REAL(KIND = r_single) :: ap(nij),as(nij),an(nij),aw(nij),ae(nij),su(nij),phi(nij)
INTEGER :: ia,i,j,ijn,ije,ijw,ijs,ij

  ia = 0
  DO i = 2,nim
    IF( i == 2)THEN
      DO j = 2,njm
        ij = li(i) + j; ije=ij+nj; ijn=ij+1; ijw=ij-nj; ijs=ij-1
        rhs(lli(ij)) = DBLE( su(ij) )
        sol(lli(ij)) = DBLE( phi(ij))
        IF(j == 2)THEN
          ia = ia + 1
          Acoo(ia) = DBLE( ae(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( an(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( ap(ij) )
        ELSEIF(j == njm)THEN
          ia = ia + 1
          Acoo(ia) = DBLE( ae(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( as(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( ap(ij) )
        ELSE
          ia = ia + 1
          Acoo(ia) = DBLE( ae(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( an(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( as(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( ap(ij) )
        ENDIF
      ENDDO
    ELSEIF(i == nim)THEN
       DO j = 2,njm
        ij = li(i) + j; ije=ij+nj; ijn=ij+1; ijw=ij-nj; ijs=ij-1
        rhs(lli(ij)) = DBLE( su(ij))
        sol(lli(ij)) = DBLE(phi(ij))
        IF(j == 2)THEN
          ia = ia + 1
          Acoo(ia) = DBLE( an(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( aw(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( ap(ij) )
        ELSEIF(j == njm)THEN
          ia = ia + 1
          Acoo(ia) = DBLE( aw(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( as(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( ap(ij) )
        ELSE
          ia = ia + 1
          Acoo(ia) = DBLE( an(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( aw(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( as(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( ap(ij) )
        ENDIF
      ENDDO                       
    ELSE
       DO j = 2,njm
        ij = li(i) + j; ije=ij+nj; ijn=ij+1; ijw=ij-nj; ijs=ij-1
        rhs(lli(ij)) = DBLE( su(ij) )
        sol(lli(ij)) = DBLE( phi(ij))
        IF(j == 2)THEN
          ia = ia + 1
          Acoo(ia) = DBLE( ae(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( an(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( aw(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( ap(ij) )
        ELSEIF(j == njm)THEN
          ia = ia + 1
          Acoo(ia) = DBLE( ae(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( aw(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( as(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( ap(ij) )
        ELSE
          ia = ia + 1
          Acoo(ia) = DBLE( ae(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( an(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( aw(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( as(ij) )
          ia = ia + 1
          Acoo(ia) = DBLE( ap(ij) )
        ENDIF
      ENDDO 
    ENDIF
  ENDDO

  IF(ia /= nnz)WRITE(*,*)'fd_spkit_interface:: ia = ',ia,' ,NNZ = ',NNZ
END SUBROUTINE fd_cooVal_create


END MODULE  modfd_solve_linearsys