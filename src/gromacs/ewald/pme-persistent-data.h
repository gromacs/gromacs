#ifndef GMX_EWALD_PME_PERSISTENT_DATA_H
#define GMX_EWALD_PME_PERSISTENT_DATA_H

#include "config.h"

#include "gromacs/gpu_utils/gputraits_ocl.h"

#include <memory>


struct PmeGpu; //FIXME remove
/*! \brief \internal
 * PME persistent host data, which is initialized once-per-context for the whole execution.
 * Will help to not recompilekernels for each OpenCL unit test.
 * //TODO - move the stream here as well?
 */
struct PmeGpuPersistentData
{
    Context context;

  //! Conveniently all the PME kernels use the same single argument type
#if GMX_GPU == GMX_GPU_OPENCL
    using PmeKernelHandle = cl_kernel;

      cl_program program;
#elif GMX_GPU == GMX_GPU_CUDA
    using PmeProgramHandle = void *;
#endif
    // TODO: All these kernels are compiled during pme_gpu_init() only for the given PME order!
    // (and only order of 4 is supported now, anyway).
    // spreading kernels also have hardcoded X/Y indices wrapping parameters as a placeholder for implementing
    // 1/2D decomposition.
    PmeKernelHandle splineKernel;
    PmeKernelHandle spreadKernel;
    PmeKernelHandle splineAndSpreadKernel;
    // Same for gather: hardcoded X/Y unwrap parameters, order of 4,
    // + it can reduce with previous forces in the host buffer, or ignore it.
    PmeKernelHandle gatherReduceWithInputKernel;
    PmeKernelHandle gatherKernel;
    // solve kernel doesn't care about spline order, but can optionally compute energy and virial,
    // and supports XYZ and YZX grid orderings.
    PmeKernelHandle solveYZXKernel;
    PmeKernelHandle solveXYZKernel;
    PmeKernelHandle solveYZXEnergyKernel;
    PmeKernelHandle solveXYZEnergyKernel;
    // There are also FFT kernels which are managed entirely by cu/clFFT

    PmeGpuPersistentData() = delete;
    explicit PmeGpuPersistentData(PmeGpu *pmeGpu);
    ~PmeGpuPersistentData();

private:
    void pme_gpu_compile_kernels(PmeGpu *pmeGpu);
};

using PmePersistentDataHandle = std::shared_ptr<PmeGpuPersistentData>;

#endif // PMEPERSISTENTDATA_H
