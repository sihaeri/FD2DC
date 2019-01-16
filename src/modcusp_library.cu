#include "modcusp_library.h"

int cusp_biCGSTAB_solver::cusp_biCGSTAB_initDevice(indexType devID)
{
  int deviceCount = 0;
  cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
  if (error_id != cudaSuccess)
  {
    printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)error_id, cudaGetErrorString(error_id));
    return OPERROR;
  }
  if (deviceCount == 0)
    {
      printf("There are no available CUDA device(s). Reverting to a CPU Solver\n");
      return NODEVICE;
    }
  else
    {
      printf("Detected %d CUDA Capable device(s)\n", deviceCount);
      if ( devID >= deviceCount)
        {
          printf("Device id=$d not found. Maximum id=%d. Reverting to a CPU solver \n", devID, deviceCount);
          return NODEVICE;
        }
      else
        {
          cudaSetDevice(devID);
          cudaDeviceProp deviceProp;
          cudaGetDeviceProperties(&deviceProp, devID);
          printf("\nRunning on device %d: \"%s\"\n", devID, deviceProp.name);
          char msg[256];
          SPRINTF(msg, "Total amount of global memory:                 %.0f MBytes (%llu bytes)\n", (float)deviceProp.totalGlobalMem/1048576.0f, (unsigned long long) deviceProp.totalGlobalMem);
          printf("%s", msg);
          printf("  (%2d) Multiprocessors, (%3d) CUDA Cores/MP:     %d CUDA Cores\n",
                 deviceProp.multiProcessorCount,
                 _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor),
                 _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) * deviceProp.multiProcessorCount);
        }
    }

  return OPSUCCESS;
}

int cusp_biCGSTAB_solver::cusp_biCGSTAB_allocDevice(indexType n, indexType m)
{
  N = n;
  nnz = m;
  if(N <= 0 or nnz <= N)
    {
      printf("The size of the coeffcient matrix is not set correctly, N=%d, NNZ=%d\n.", n, m);
      return OPERROR;
    }
  if(cudaMalloc(&cooRowIndADev, nnz*sizeof(indexType)) != cudaSuccess)return OPERROR;
  if(cudaMalloc(&cooColIndADev, nnz*sizeof(indexType)) != cudaSuccess)return OPERROR;
  if(cudaMalloc(&cooValADev, nnz*sizeof(valueType)) != cudaSuccess)return OPERROR;
  if(cudaMalloc(&xDev, N*sizeof(valueType)) != cudaSuccess)return OPERROR;
  if(cudaMalloc(&bDev, N*sizeof(valueType)) != cudaSuccess)return OPERROR;
  return OPSUCCESS;
}

int cusp_biCGSTAB_solver::cusp_biCGSTAB_copyH2D_AInds(indexType *rows, indexType *cols)
{
  if(cudaMemcpy(cooRowIndADev, rows, nnz*sizeof(cooRowIndADev[0]),
                cudaMemcpyHostToDevice) != cudaSuccess) return OPERROR;
  if(cudaMemcpy(cooColIndADev, cols, nnz*sizeof(cooColIndADev[0]),
                cudaMemcpyHostToDevice) != cudaSuccess) return OPERROR;
  return OPSUCCESS;
}

int cusp_biCGSTAB_solver::cusp_biCGSTAB_copyH2D_system(valueType *Avals, valueType *xHost, valueType *bHost)
{
  if(cudaMemcpy(bDev, bHost, N*sizeof(bDev[0]),cudaMemcpyHostToDevice) != cudaSuccess) return OPERROR;
  if(cudaMemcpy(xDev, xHost, N*sizeof(xDev[0]),cudaMemcpyHostToDevice) != cudaSuccess) return OPERROR;
  if(cudaMemcpy(cooValADev, Avals, nnz*sizeof(cooValADev[0]),
                cudaMemcpyHostToDevice) != cudaSuccess) return OPERROR;

  return OPSUCCESS;
}

int cusp_biCGSTAB_solver::cusp_biCGSTAB_copyD2H_x(valueType *xHost)
{
  if(cudaMemcpy(xHost, xDev, N*sizeof(xHost[0]),cudaMemcpyDeviceToHost) != cudaSuccess) return OPERROR;
  return OPSUCCESS;
}

int cusp_biCGSTAB_solver::cusp_biCGSTAB_solveDev_system(valueType relTol,
                                                      valueType absTol,
                                                      indexType maxItr)
{
  // Wrap device pointers
  thrust::device_ptr<indexType> wrapped_cooRowIndADev(cooRowIndADev);
  thrust::device_ptr<indexType> wrapped_cooColIndADev(cooColIndADev);
  thrust::device_ptr<valueType> wrapped_cooValADev(cooValADev);
  thrust::device_ptr<valueType> wrapped_xDev(xDev);
  thrust::device_ptr<valueType> wrapped_bDev(bDev);

  // Wrap in cusp array1d
  deviceIndexArrayView rowInds (wrapped_cooRowIndADev, wrapped_cooRowIndADev + nnz);
  deviceIndexArrayView colInds (wrapped_cooColIndADev, wrapped_cooColIndADev + nnz);
  deviceValueArrayView values  (wrapped_cooValADev, wrapped_cooValADev + nnz);
  deviceValueArrayView x       (wrapped_xDev, wrapped_xDev + N);
  deviceValueArrayView b       (wrapped_bDev, wrapped_bDev + N);

  // Create coo_matrix_view from the 3 array1d views
  deviceView A(N, N, nnz, rowInds, colInds, values);

  // Setup a monitor and solve
  cusp::monitor<valueType> monitor(b, maxItr, relTol, absTol, false);
  cusp::precond::diagonal<valueType, devMemorySpace> M(A);

  cusp::krylov::bicgstab(A, x, b, monitor, M);

  residuals = monitor.residual_norm();
  solverItr = monitor.iteration_count();

  return OPSUCCESS;

}

int cusp_biCGSTAB_solver::cusp_biCGSTAB_getMonitor(valueType &res, indexType &nItr)
{
  nItr = solverItr;
  res = residuals;
  return OPSUCCESS;
}

int cusp_biCGSTAB_solver::cusp_biCGSTAB_shutdown()
{

  if(cudaFree(cooRowIndADev) != cudaSuccess)return OPERROR;
  if(cudaFree(cooColIndADev) != cudaSuccess)return OPERROR;
  if(cudaFree(cooValADev) != cudaSuccess)return OPERROR;
  if(cudaFree(xDev) != cudaSuccess)return OPERROR;
  if(cudaFree(bDev) != cudaSuccess)return OPERROR;

  return OPSUCCESS;
}

/******************************************************/
//External Interfaces                                  /
/******************************************************/
extern "C" void* getInstance_cusp_biCGSTAB_solver()
{
  cusp_biCGSTAB_solver *cusp_biCGSTAB_solver_ = new cusp_biCGSTAB_solver(0,0,0.,0.);
  return static_cast<void *>(cusp_biCGSTAB_solver_);
}

extern "C" int cusp_biCGSTAB_initDevice_intrf(void *cusp_biCGSTAB_solver_ptr, indexType devID)
{
  cusp_biCGSTAB_solver *ptr = static_cast<cusp_biCGSTAB_solver*>(cusp_biCGSTAB_solver_ptr);
  return ptr->cusp_biCGSTAB_initDevice(devID);
}

extern "C" int cusp_biCGSTAB_allocDevice_intrf(void *cusp_biCGSTAB_solver_ptr, indexType n, indexType m)
{
  cusp_biCGSTAB_solver *ptr = static_cast<cusp_biCGSTAB_solver*>(cusp_biCGSTAB_solver_ptr);
  return ptr->cusp_biCGSTAB_allocDevice(n,m);
}

extern "C" int cusp_biCGSTAB_copyH2D_AInds_intrf(void *cusp_biCGSTAB_solver_ptr,
                                             indexType *rows, indexType *cols)
{
  cusp_biCGSTAB_solver *ptr = static_cast<cusp_biCGSTAB_solver*>(cusp_biCGSTAB_solver_ptr);
  return ptr->cusp_biCGSTAB_copyH2D_AInds(rows,cols);
}

extern "C" int cusp_biCGSTAB_copyH2D_system_intrf(void *cusp_biCGSTAB_solver_ptr,valueType *Avals, valueType *xHost, valueType *bHost)
{
  cusp_biCGSTAB_solver *ptr = static_cast<cusp_biCGSTAB_solver*>(cusp_biCGSTAB_solver_ptr);
  return ptr->cusp_biCGSTAB_copyH2D_system(Avals, xHost, bHost);

}

extern "C" int cusp_biCGSTAB_copyD2H_x_intrf(void *cusp_biCGSTAB_solver_ptr,valueType *xHost)
{
  cusp_biCGSTAB_solver *ptr = static_cast<cusp_biCGSTAB_solver*>(cusp_biCGSTAB_solver_ptr);
  return ptr->cusp_biCGSTAB_copyD2H_x(xHost);

}

extern "C" int cusp_biCGSTAB_solveDev_system_intrf(void *cusp_biCGSTAB_solver_ptr,
                                        valueType relTol, valueType absTol,
                                        indexType maxItr)
{
  cusp_biCGSTAB_solver *ptr = static_cast<cusp_biCGSTAB_solver*>(cusp_biCGSTAB_solver_ptr);
  return ptr->cusp_biCGSTAB_solveDev_system(relTol, absTol, maxItr);
}

extern "C" int cusp_biCGSTAB_getMonitor_intrf(void *cusp_biCGSTAB_solver_ptr,
                                        valueType &residual, indexType &nItr)
{
  cusp_biCGSTAB_solver *ptr = static_cast<cusp_biCGSTAB_solver*>(cusp_biCGSTAB_solver_ptr);
  return ptr->cusp_biCGSTAB_getMonitor(residual, nItr);
}

extern "C" int cusp_biCGSTAB_shutdown_intrf(void *cusp_biCGSTAB_solver_ptr)
{
  cusp_biCGSTAB_solver *ptr = static_cast<cusp_biCGSTAB_solver*>(cusp_biCGSTAB_solver_ptr);
  return ptr->cusp_biCGSTAB_shutdown();
}
