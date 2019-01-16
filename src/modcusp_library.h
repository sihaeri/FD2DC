#include <stdlib.h>
#include <memory>
#include <iostream>

#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <cusp/coo_matrix.h>
#include <cusp/monitor.h>
#include <cusp/krylov/bicgstab.h>
#include <cusp/array1d.h>
#include <cusp/precond/diagonal.h>
#include <cusp/precond/aggregation/smoothed_aggregation.h>
#include <cusp/precond/ainv.h>

#define SOLVER_FAILED -2
#define OPSUCCESS 0
#define OPERROR -1
#define NODEVICE -2
#define SOLVER_DONE 1

typedef cusp::device_memory devMemorySpace;
typedef double valueType;
typedef int indexType;
typedef typename cusp::array1d_view< thrust::device_ptr<indexType> > deviceIndexArrayView;
typedef typename cusp::array1d_view< thrust::device_ptr<valueType> > deviceValueArrayView;
typedef cusp::coo_matrix_view<deviceIndexArrayView, deviceIndexArrayView, deviceValueArrayView> deviceView;

class cusp_biCGSTAB_solver
{
  int N, nnz, devID;
  //cusp::coo_matrix<indexType, valueType, devMemorySpace> A; // system to be solved A.x=b
  //cusp::array1d<valueType, devMemorySpace> x;
  //cusp::array1d<valueType, devMemorySpace> b;
  indexType *cooRowIndADev;
  indexType *cooColIndADev;
  valueType *cooValADev;
  valueType *xDev;
  valueType *bDev;
  cusp::monitor<valueType> monitor();
  valueType residuals;
  indexType solverItr;
public:
  cusp_biCGSTAB_solver(indexType n_,indexType nnz_,valueType solTol_,indexType solMaxIt_):
    N(n_),nnz(nnz_){}
  int cusp_biCGSTAB_initDevice(indexType);
  int cusp_biCGSTAB_allocDevice(indexType, indexType );
  int cusp_biCGSTAB_copyH2D_AInds(indexType *, indexType *);
  int cusp_biCGSTAB_copyH2D_system(valueType*, valueType*, valueType*);
  int cusp_biCGSTAB_copyD2H_x(valueType*);
  int cusp_biCGSTAB_solveDev_system(valueType, valueType, indexType);
  //int cusp_biCGSTAB_setMonitor(valueType, indexType);
  int cusp_biCGSTAB_getMonitor(valueType &, indexType &);
  int cusp_biCGSTAB_shutdown();
};

/* N: Number of unknowns
   nnz: Number of non-zero elements
   devID: device ID to run the simulations on
   cusp_biCGSTAB_copyH2D_A copy matrix A from host to device using coo format (rows, cols, vals)
   cusp_biCGSTAB_copyH2D_x/b copy x and rhs from host to device (A.x=b)
   cusp_biCGSTAB_copyD2H_x copy the solution back to the host
   cusp_biCGSTAB_solveDev_sys wraps the linear system arrays and solves the system
*/
