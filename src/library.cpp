#include "library.h"

int cuBlas_biCGSTAB::cuBlas_initDevice()
{

  if (cublasInit() != CUBLAS_STATUS_SUCCESS) return OPERROR;

  if (cusparseCreate(&handle) != CUSPARSE_STATUS_SUCCESS) return OPERROR;

  // matrix descriptor 
  if (cusparseCreateMatDescr(&descra) != CUSPARSE_STATUS_SUCCESS) return OPERROR;

  cusparseSetMatType(descra,CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(descra,CUSPARSE_INDEX_BASE_ONE);

  return OPSUCCESS;
};
//
int cuBlas_biCGSTAB::cuBlas_allocDevice(int nz, int m)
{
	nnz = nz;
	N   = m;
	if(cudaMalloc((void**)&cooRowIndAdev,nnz*sizeof(cooRowIndAdev[0])) !=cudaSuccess)return OPERROR;
	if(cudaMalloc((void**)&cooColIndAdev,nnz*sizeof(cooColIndAdev[0])) !=cudaSuccess)return OPERROR;
	if(cudaMalloc((void**)&cooValAdev, nnz*sizeof(cooValAdev[0])) !=cudaSuccess)return OPERROR;
	if(cudaMalloc((void**)&r, N*sizeof(r[0])) !=cudaSuccess)return OPERROR;
	if(cudaMalloc((void**)&s, N*sizeof(s[0])) !=cudaSuccess)return OPERROR;
	if(cudaMalloc((void**)&t, N*sizeof(t[0])) !=cudaSuccess)return OPERROR;
	if(cudaMalloc((void**)&r_tld, N*sizeof(r_tld[0]))!=cudaSuccess)return OPERROR;
	if(cudaMalloc((void**)&p, N*sizeof(p[0]))!=cudaSuccess)return OPERROR;
	if(cudaMalloc((void**)&p_hat, N*sizeof(p_hat[0]))!=cudaSuccess)return OPERROR;
	if(cudaMalloc((void**)&s_hat, N*sizeof(s_hat[0]))!=cudaSuccess)return OPERROR;
	if(cudaMalloc((void**)&v, N*sizeof(v[0]))!=cudaSuccess)return OPERROR;
	if(cudaMalloc((void**)&x, N*sizeof(x[0]))!=cudaSuccess)return OPERROR;
	if(cudaMalloc((void**)&csrRowPtrAdev,(N+1)*sizeof(csrRowPtrAdev[0]))!=cudaSuccess)return OPERROR;
	return OPSUCCESS;
};
//
int cuBlas_biCGSTAB::cuBlas_copyH2D_coords(int *cooRowIndAhost,
								int *cooColIndAhost)
{

	if(cudaMemcpy(cooRowIndAdev, cooRowIndAhost,
		 (size_t)(nnz*sizeof(cooRowIndAdev[0])),
		 cudaMemcpyHostToDevice) != cudaSuccess)return OPERROR;
	if(cudaMemcpy(cooColIndAdev, cooColIndAhost,
		 (size_t)(nnz*sizeof(cooColIndAdev[0])),
		 cudaMemcpyHostToDevice) != cudaSuccess)return OPERROR;
	return OPSUCCESS;
};
int cuBlas_biCGSTAB::cuBlas_copyH2D_systemDP(double *cooValAhost, double *bhost, double *xhost)
{
	if(cudaMemcpy(cooValAdev, cooValAhost,
			 (size_t)(nnz*sizeof(cooValAdev[0])),
			 cudaMemcpyHostToDevice) != cudaSuccess)return OPERROR;	
	if(cudaMemcpy(r, bhost, (size_t)(N*sizeof(r[0])),cudaMemcpyHostToDevice) != cudaSuccess)return OPERROR;
	if(cudaMemcpy(x, xhost, (size_t)(N*sizeof(x[0])),cudaMemcpyHostToDevice) != cudaSuccess)return OPERROR;
	return OPSUCCESS;
};

int cuBlas_biCGSTAB::cuBlas_copyH2D_jacobPreconDP(double *P)
{
	if(cudaMemcpy(jacPreconDev, P,
			 (size_t)(N*sizeof(jacPreconDev[0])),
			 cudaMemcpyHostToDevice) != cudaSuccess)return OPERROR;	
	return OPSUCCESS;
};

int cuBlas_biCGSTAB::cuBlas_copyD2H_solDP(double *xhost)
{
	if(cudaMemcpy(xhost, x, (size_t)(N*sizeof(cooValAdev[0])),cudaMemcpyDeviceToHost) != cudaSuccess)return OPERROR;	
	return OPSUCCESS;
};

int cuBlas_biCGSTAB::cuBlas_coo2csr()
{
	if(cusparseXcoo2csr(handle, cooRowIndAdev, nnz,N,csrRowPtrAdev,CUSPARSE_INDEX_BASE_ONE) != CUSPARSE_STATUS_SUCCESS)return OPERROR;
	return OPSUCCESS;
};

int cuBlas_biCGSTAB::cuBlas_biCGSTAB_shutdown()
{

	if(cudaFree(cooRowIndAdev)!= cudaSuccess)return OPERROR;
	if(cudaFree(cooColIndAdev)!= cudaSuccess)return OPERROR;
	if(cudaFree(cooValAdev)!= cudaSuccess)return OPERROR;
	if(cudaFree(r)!= cudaSuccess)return OPERROR;
	if(cudaFree(s)!= cudaSuccess)return OPERROR;
	if(cudaFree(t)!= cudaSuccess)return OPERROR;
	if(cudaFree(r_tld)!= cudaSuccess)return OPERROR;
	if(cudaFree(p)!= cudaSuccess)return OPERROR;
	if(cudaFree(p_hat)!= cudaSuccess)return OPERROR;
	if(cudaFree(s_hat)!= cudaSuccess)return OPERROR;
	if(cudaFree(v)!= cudaSuccess)return OPERROR;
	if(cudaFree(x)!= cudaSuccess)return OPERROR;
	if(cudaFree(csrRowPtrAdev)!= cudaSuccess)return OPERROR;
	
	if(cusparseDestroy(handle)!= CUSPARSE_STATUS_SUCCESS)return OPERROR;
	cublasShutdown();

	return OPSUCCESS;
};

int cuBlas_biCGSTAB::cuBlas_biCGSTAB_setStop(int MaxIt, double TOL )
{
	solverTol = TOL;
	solverMaxIt = MaxIt;
	return OPSUCCESS;
};

int cuBlas_biCGSTAB::cuBlas_biCGSTAB_iterate(double *resid,int *iter)
{
	double bnrm2,omega,rho,rho_1,alpha,beta,snrm2;

	bnrm2 = cublasDnrm2(N, r, 1);

	if(cusparseDcsrmv(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, 
		N, N, -1.0,descra, cooValAdev, csrRowPtrAdev, cooColIndAdev, x,1.0, r) != CUSPARSE_STATUS_SUCCESS)return OPERROR;

	if(cublasDnrm2 (N, r, 1) < solverTol*bnrm2) return SOLVER_DONE; //x is already close enough to r
	
	omega = 1.0;

	if (cudaMemcpy(r_tld, r, (size_t)(N*sizeof(r[0])),cudaMemcpyDeviceToDevice) != cudaSuccess)return OPERROR; 

	for( *iter = 0; *iter < solverMaxIt; (*iter)++)
	{
		rho = cublasDdot(N,r_tld,1,r,1);

		if(rho == 0) break;

		if(*iter > 0)
		{
		  beta = (rho/rho_1) * (alpha/omega);
		  cublasDaxpy (N, -omega, v, 1, p, 1);
		  cublasDaxpy (N, 1.0/beta, r, 1, p, 1);
		  cublasDscal (N, beta, p, 1);         /* p = r + beta*( p - omega*v ) */
		}
		else
		{
			//set p = r
			if(cudaMemcpy(p, r,	 (size_t)(N*sizeof(r[0])),cudaMemcpyDeviceToDevice) != cudaSuccess)return OPERROR;    
		}

		//Assume no preconditioner i.e M = I => p_hat = p
		if(cudaMemcpy(p_hat, p, (size_t)(N*sizeof(x[0])),cudaMemcpyDeviceToDevice) != cudaSuccess)return OPERROR;    

		// v = A*p_hat
		if(cusparseDcsrmv(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, 1.0, descra, 
						  cooValAdev, csrRowPtrAdev, cooColIndAdev, p_hat, 0.0, v) != CUSPARSE_STATUS_SUCCESS)return OPERROR;  
		
		// alph = rho / ( r_tld'*v ) 
		alpha = rho / cublasDdot (N, r_tld, 1, v, 1); 
		
		if(cudaMemcpy(s, r, (size_t)(N*sizeof(r[0])),cudaMemcpyDeviceToDevice) != cudaSuccess)return OPERROR;

		cublasDaxpy (N, -alpha, v, 1, s, 1);
		snrm2 = cublasDnrm2( N, s, 1);


		if ( snrm2 < solverTol*bnrm2 )
		{
			//cublasDaxpy (N, alpha, p_hat, 1, s, 1);
			*resid = snrm2;
			break;
		}

		if(cudaMemcpy(s_hat, s, (size_t)(N*sizeof(s[0])),cudaMemcpyDeviceToDevice) != cudaSuccess)return OPERROR;

		/* t = A*s_hat */
		if(cusparseDcsrmv(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, 1.0,
					     descra, cooValAdev, csrRowPtrAdev, cooColIndAdev, s_hat,0.0, t) != CUSPARSE_STATUS_SUCCESS)return OPERROR;

		omega = cublasDdot(N, t, 1, s, 1)/cublasDdot(N, t, 1, t, 1);
		// x = x + alph*p_hat + omega*s_ha
		cublasDaxpy (N, alpha, p_hat, 1, x, 1);
		cublasDaxpy (N, omega, s_hat, 1, x, 1);

		cublasDaxpy (N, -omega, t, 1, s, 1);

		if(cudaMemcpy(r, s, (size_t)(N*sizeof(r[0])),cudaMemcpyDeviceToDevice) != cudaSuccess)return OPERROR;

		*resid =  cublasDnrm2( N, r, 1);
		
		if ( *resid <= solverTol*bnrm2 )break;
		if (omega == 0.0) break;
		rho_1 = rho;
	}


	if(*resid <= solverTol*bnrm2)return SOLVER_DONE;
	if(omega == 0.0)return SOLVER_FAILED;
	if(rho == 0.0)return SOLVER_FAILED;

	return SOLVER_DONE;
}

extern "C" void* instantiate_cuBlasBiCGSTAB()
{
	cuBlas_biCGSTAB *cuBlas_biCGSTAB_=new cuBlas_biCGSTAB(0,0,0.0,0);
	return static_cast<void *>(cuBlas_biCGSTAB_);
};

extern "C" int wrap_cuBlas_initDevice(void *cuBlas_biCGSTAB_ptr)
{
	cuBlas_biCGSTAB *ptr = static_cast<cuBlas_biCGSTAB*>(cuBlas_biCGSTAB_ptr);
	return ptr->cuBlas_initDevice();
};

extern "C" int wrap_cuBlas_allocDevice(void *cuBlas_biCGSTAB_ptr,int nz, int m)
{
	cuBlas_biCGSTAB *ptr = static_cast<cuBlas_biCGSTAB*>(cuBlas_biCGSTAB_ptr);
	return ptr->cuBlas_allocDevice(nz,m);
};

extern "C" int wrap_BiCGSTAB_setStop(void *cuBlas_biCGSTAB_ptr, int MaxIt, double Tol)
{
	cuBlas_biCGSTAB *ptr = static_cast<cuBlas_biCGSTAB*>(cuBlas_biCGSTAB_ptr);
	return ptr->cuBlas_biCGSTAB_setStop(MaxIt,Tol);
};

extern "C" int wrap_cuBlas_coo2csr(void *cuBlas_biCGSTAB_ptr)
{
	cuBlas_biCGSTAB *ptr = static_cast<cuBlas_biCGSTAB*>(cuBlas_biCGSTAB_ptr);
	return ptr->cuBlas_coo2csr();
};

extern "C" int wrap_cuBlas_cpH2D_sysDP(void *cuBlas_biCGSTAB_ptr,double *cooValAhost, double *bhost, double *xhost)
{
	cuBlas_biCGSTAB *ptr = static_cast<cuBlas_biCGSTAB*>(cuBlas_biCGSTAB_ptr);
	return ptr->cuBlas_copyH2D_systemDP(cooValAhost, bhost, xhost);
};

extern "C" int wrap_cuBlas_cpH2D_coords(void *cuBlas_biCGSTAB_ptr,int *cooRowPtrAhost,int *cooColPtrAhost)
{
	cuBlas_biCGSTAB *ptr = static_cast<cuBlas_biCGSTAB*>(cuBlas_biCGSTAB_ptr);
	return ptr->cuBlas_copyH2D_coords(cooRowPtrAhost, cooColPtrAhost);
};

extern "C" int wrap_cuBlas_biCGSTAB_itr(void *cuBlas_biCGSTAB_ptr,double *resid,int *iter)
{
	cuBlas_biCGSTAB *ptr = static_cast<cuBlas_biCGSTAB*>(cuBlas_biCGSTAB_ptr);
	return ptr->cuBlas_biCGSTAB_iterate(resid,iter);
};

extern "C" int wrap_cuBlas_cpD2H_solDP(void *cuBlas_biCGSTAB_ptr,double *xhost)
{
	cuBlas_biCGSTAB *ptr = static_cast<cuBlas_biCGSTAB*>(cuBlas_biCGSTAB_ptr);
	return ptr->cuBlas_copyD2H_solDP(xhost);
};

extern "C" int wrap_cuBlas_shutdown(void *cuBlas_biCGSTAB_ptr)
{
	cuBlas_biCGSTAB *ptr = static_cast<cuBlas_biCGSTAB*>(cuBlas_biCGSTAB_ptr);
	if(ptr->cuBlas_biCGSTAB_shutdown() != OPSUCCESS)return OPERROR;
	delete cuBlas_biCGSTAB_ptr;
	return OPSUCCESS;
};