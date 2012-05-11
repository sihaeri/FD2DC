      SUBROUTINE solve_hypre_sys_init(mpi_comm,myrank_mpi,nprocs_mpi,
     x ncel,A,b,x)

      IMPLICIT NONE
      
      INCLUDE 'mpif.h'

      INTEGER*8 A,b,x,mpi_comm
      INTEGER   myrank_mpi,nprocs_mpi
      INTEGER   ncel,ilower,iupper,ierr,HYPRE_PARCSR
      !--From HYPRE.C
      PARAMETER (HYPRE_PARCSR = 5555)
      ilower = 1
      iupper = ncel

      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank_mpi,ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs_mpi,ierr)

      mpi_comm = MPI_COMM_WORLD

      !--Create The matrix
      CALL HYPRE_IJMatrixCreate( mpi_comm, ilower,
     x iupper, ilower, iupper, A, ierr)

      !--Choose The format
      CALL HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR, ierr)
      
      !--Initialize the Matrix 
      CALL HYPRE_IJMatrixInitialize(A, ierr)

      !--Create rhs and sol
      CALL HYPRE_IJVectorCreate(mpi_comm, ilower, iupper,b, ierr)
      CALL HYPRE_IJVectorCreate(mpi_comm, ilower, iupper,x, ierr)

      CALL HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR, ierr)
      CALL HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR, ierr)

      CALL HYPRE_IJVectorInitialize(b, ierr)
      CALL HYPRE_IJVectorInitialize(x, ierr)
  
      END 
      
      SUBROUTINE solve_hypre(mpi_comm,li,lli,ncel,ni,nj,ap,as,an,aw,ae,
     x su,phi,tol,res,A,b,x,solver_type,maxit)

      IMPLICIT NONE

      INTEGER ncel,ni,nj,lli((ni*nj)),li((ni*nj)),solver_type,maxit
      INTEGER*8  mpi_comm,A,b,x
      DOUBLE PRECISION  ap((ni*nj)),as((ni*nj)),an((ni*nj)),aw((ni*nj)),
     x ae((ni*nj)),su((ni*nj)),phi((ni*nj)),res,tol

      INTEGER    ierr
      INTEGER    solver_id
      INTEGER    nnz,ilower,iupper,i,nim,njm,ije,ijw,ijs,ijn,ij,nc,j
      INTEGER    precond_id
      DOUBLE PRECISION rhs_values(ncel)
      DOUBLE PRECISION x_values(ncel)
      INTEGER    rows(ncel)
      INTEGER    cols(5)
      DOUBLE PRECISION values(5)
      INTEGER    num_iterations,nlog
      INTEGER*8  precond,ParCSR_A,Par_b,Par_x,solver

      ilower = 1
      iupper = ncel
      nim = ni - 1
      njm = nj - 1
      nc = 0
      DO i = 2,nim
        IF( i == 2)THEN
          DO j = 2,njm
            nnz = 1;nc = nc + 1
            ij = li(i) + j; ije=ij+nj; ijn=ij+1; ijw=ij-nj; ijs=ij-1
            rhs_values(lli(ij)) = su(ij)
            x_values(lli(ij)) = phi(ij)
            rows(lli(ij))     = lli(ij)
            IF(j == 2)THEN
              cols(nnz) = lli(ije); values(nnz) = ae(ij); nnz = nnz + 1
              cols(nnz) = lli(ijn); values(nnz) = an(ij); nnz = nnz + 1
              cols(nnz) = lli(ij);  values(nnz) = ap(ij)
              CALL HYPRE_IJMatrixSetValues(A, 1, nnz , lli(ij), cols, 
     x                                     values, ierr)
            ELSEIF(j == njm)THEN
              cols(nnz) = lli(ije); values(nnz) = ae(ij); nnz = nnz + 1
              cols(nnz) = lli(ijs); values(nnz) = as(ij); nnz = nnz + 1
              cols(nnz) = lli(ij);  values(nnz) = ap(ij)
              CALL HYPRE_IJMatrixSetValues(A, 1, nnz , lli(ij), cols, 
     x                                     values, ierr)

            ELSE
              cols(nnz) = lli(ije); values(nnz) = ae(ij); nnz = nnz + 1
              cols(nnz) = lli(ijn); values(nnz) = an(ij); nnz = nnz + 1
              cols(nnz) = lli(ijs); values(nnz) = as(ij); nnz = nnz + 1
              cols(nnz) = lli(ij);  values(nnz) = ap(ij)
              CALL HYPRE_IJMatrixSetValues(A, 1, nnz , lli(ij), cols, 
     x                                     values, ierr)

            ENDIF
          ENDDO
        ELSEIF(i == nim)THEN
          DO j = 2,njm
            nnz = 1;nc = nc + 1
            ij = li(i) + j; ije=ij+nj; ijn=ij+1; ijw=ij-nj; ijs=ij-1
            rhs_values(lli(ij)) = su(ij)
            x_values(lli(ij)) = phi(ij)
            rows(lli(ij))     = lli(ij)
            IF(j == 2)THEN
              cols(nnz) = lli(ijw); values(nnz) = aw(ij); nnz = nnz + 1
              cols(nnz) = lli(ijn); values(nnz) = an(ij); nnz = nnz + 1
              cols(nnz) = lli(ij);  values(nnz) = ap(ij)
              CALL HYPRE_IJMatrixSetValues(A, 1, nnz , lli(ij), cols, 
     x                                     values, ierr)
            ELSEIF(j == njm)THEN
              cols(nnz) = lli(ijw); values(nnz) = aw(ij); nnz = nnz + 1
              cols(nnz) = lli(ijs); values(nnz) = as(ij); nnz = nnz + 1
              cols(nnz) = lli(ij);  values(nnz) = ap(ij)
              CALL HYPRE_IJMatrixSetValues(A, 1, nnz , lli(ij), cols, 
     x                                     values, ierr)
         
            ELSE
              cols(nnz) = lli(ijw); values(nnz) = aw(ij); nnz = nnz + 1
              cols(nnz) = lli(ijn); values(nnz) = an(ij); nnz = nnz + 1
              cols(nnz) = lli(ijs); values(nnz) = as(ij); nnz = nnz + 1
              cols(nnz) = lli(ij);  values(nnz) = ap(ij)
              CALL HYPRE_IJMatrixSetValues(A, 1, nnz , lli(ij), cols, 
     x                                     values, ierr)
            ENDIF
          ENDDO                       
        ELSE
          DO j = 2,njm
            nnz = 1;nc = nc + 1
            ij = li(i) + j; ije=ij+nj; ijn=ij+1; ijw=ij-nj; ijs=ij-1
            rhs_values(lli(ij)) = su(ij)
            x_values(lli(ij)) = phi(ij)
            rows(lli(ij))     = lli(ij)
            IF(j == 2)THEN
              cols(nnz) = lli(ije); values(nnz) = ae(ij); nnz = nnz + 1
              cols(nnz) = lli(ijw); values(nnz) = aw(ij); nnz = nnz + 1
              cols(nnz) = lli(ijn); values(nnz) = an(ij); nnz = nnz + 1
              cols(nnz) = lli(ij);  values(nnz) = ap(ij)
              CALL HYPRE_IJMatrixSetValues(A, 1, nnz , lli(ij), cols, 
     x                                     values, ierr)
            ELSEIF(j == njm)THEN
              cols(nnz) = lli(ije); values(nnz) = ae(ij); nnz = nnz + 1
              cols(nnz) = lli(ijw); values(nnz) = aw(ij); nnz = nnz + 1
              cols(nnz) = lli(ijs); values(nnz) = as(ij); nnz = nnz + 1
              cols(nnz) = lli(ij);  values(nnz) = ap(ij)
              CALL HYPRE_IJMatrixSetValues(A, 1, nnz , lli(ij), cols, 
     x                                     values, ierr)
            ELSE
              cols(nnz) = lli(ije); values(nnz) = ae(ij); nnz = nnz + 1
              cols(nnz) = lli(ijw); values(nnz) = aw(ij); nnz = nnz + 1
              cols(nnz) = lli(ijn); values(nnz) = an(ij); nnz = nnz + 1
              cols(nnz) = lli(ijs); values(nnz) = as(ij); nnz = nnz + 1
              cols(nnz) = lli(ij);  values(nnz) = ap(ij)
              CALL HYPRE_IJMatrixSetValues(A, 1, nnz , lli(ij), cols, 
     x                                     values, ierr)
            ENDIF
          ENDDO 
        ENDIF
      ENDDO

      IF(nc /= ncel)THEN
        WRITE(*,*)'Error in solver_hypre.'
        RETURN
      ENDIF

      CALL HYPRE_IJVectorSetValues(b, nc, rows, rhs_values, ierr)
      CALL HYPRE_IJVectorSetValues(x, nc, rows, x_values, ierr)
      
      CALL HYPRE_IJMatrixAssemble(A, ierr)
      CALL HYPRE_IJVectorAssemble(b, ierr)
      CALL HYPRE_IJVectorAssemble(x, ierr)
      
      CALL HYPRE_IJMatrixGetObject(A, ParCSR_A,ierr)
      CALL HYPRE_IJVectorGetObject(b, Par_b, ierr)
      CALL HYPRE_IJVectorGetObject(x, Par_x, ierr)
      
      IF(solver_type == 1)THEN
        
        !--Create solver
        CALL HYPRE_BoomerAMGCreate(solver, ierr)

        !--print solve info + parameters 
        CALL HYPRE_BoomerAMGSetPrintLevel(solver, 0, ierr)  
        !--Falgout coarsening
        CALL HYPRE_Boome rAMGSetCoarsenType(solver, 6, ierr) 
        !--G-S/Jacobi hybrid relaxation 
        CALL HYPRE_BoomerAMGSetRelaxType(solver, 0, ierr)     
        !--Sweeeps on each level
        CALL HYPRE_BoomerAMGSetNumSweeps(solver, 1, ierr)  
        !--conv. tolerance
        CALL HYPRE_BoomerAMGSetTol(solver, tol, ierr)
        !--max. itr
        CALL HYPRE_BoomerAMGSetMaxIter(solver, maxit, ierr)
        CALL HYPRE_BoomerAMGSetMaxLevels(solver, 10, ierr)

        !CALL HYPRE_ParCSRBiCGSTABSetPrecond(solver, precond_id, precond, 
        !x ierr)

        !--Now setup and solve!
        CALL HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x, 
     x  ierr )
        CALL HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x,
     x  ierr )

        !--Run info - needed logging turned on 
        !CALL HYPRE_BoomerAMGGetNumIterations(solver, num_iterations,ierr)
        CALL HYPRE_BoomerAMGGetFinalReltvRes(solver, res,ierr)
      
      ELSEIF(solver_type == 2)THEN
        
        CALL HYPRE_ParCSRBiCGSTABCreate(mpi_comm, solver, ierr)
        CALL HYPRE_ParCSRBiCGSTABSetMaxIter(solver, maxit, ierr)
        CALL HYPRE_ParCSRBiCGSTABSetTol(solver, tol, ierr)
        CALL HYPRE_ParCSRBiCGSTABSetLogging(solver, nlog, ierr)
        
        precond_id = 1
        precond = 0

        CALL CHYPRE_ParCSRBiCGSTABSetPrecond(solver,precond_id,
     x       precond,ierr)

        CALL HYPRE_ParCSRBiCGSTABSetup(solver, parcsr_A, par_b, par_x,
     x       ierr)

        CALL HYPRE_ParCSRBiCGSTABSolve(solver, parcsr_A, par_b, par_x,
     x       ierr)

        CALL HYPRE_ParCSRBiCGSTABGetFinalRel(solver, res, ierr)
      ENDIF

      !--Copyback the results
      DO i = 1,ncel
        x_values(i) = 0.
      ENDDO

      CALL HYPRE_IJVectorGetValues(x, nc, rows, x_values, ierr) 
      DO i = 2,nim
        DO j = 2,njm
          ij = li(i) + j
          phi(ij) = x_values(lli(ij))
        ENDDO
      ENDDO
      !--Destroy solver
      IF(solver_type == 1)CALL HYPRE_BoomerAMGDestroy(solver, ierr )
      IF(solver_type == 2)CALL HYPRE_ParCSRBiCGSTABDestroy(solver, ierr)

      RETURN

      END 
