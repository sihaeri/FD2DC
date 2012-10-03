#include <stdlib.h>
#include "cublas.h"
#include "cusparse.h"

#define OPERROR -1
#define SOLVER_FAILED -2
#define OPSUCCESS 0
#define SOLVER_DONE 1

class cuBlas_biCGSTAB
{
	int N, nnz;
	int *cooRowIndAdev;
	int *cooColIndAdev;
	int *csrRowPtrAdev;
	double *cooValAdev;
	double *x, *r, *r_tld, *p, *p_hat, *s, *s_hat, *t, *v;
	double *jacPreconDev;
	cusparseHandle_t handle;
    cusparseMatDescr_t descra;
	double solverTol;
	int solverMaxIt;
public:
	cuBlas_biCGSTAB(int n_,int nnz_,double solTol_,int solMaxIt_):N(n_),nnz(nnz_),solverTol(solTol_),solverMaxIt(solMaxIt_),
	cooRowIndAdev(),cooColIndAdev(),csrRowPtrAdev(),cooValAdev(),x(),r(),r_tld(),p(),p_hat(),s(),s_hat(),t(),v(),
	jacPreconDev(){}
	int cuBlas_allocDevice(int, int);
	int cuBlas_initDevice();
	int cuBlas_copyH2D_coords(int *,int *);
	int cuBlas_copyH2D_systemDP(double *,double *,double *);
	int cuBlas_copyH2D_jacobPreconDP(double *);
	int cuBlas_coo2csr();
	int cuBlas_biCGSTAB_iterate(double *,int *);
	int cuBlas_biCGSTAB_setStop(int, double);
	int cuBlas_copyD2H_solDP(double *);
	int cuBlas_biCGSTAB_shutdown();
};

