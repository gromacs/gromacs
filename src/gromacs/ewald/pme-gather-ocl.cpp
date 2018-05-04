#include "gmxpre.h"

#include <cassert>

#include "gromacs/gpu_utils/gputraits_ocl.h"

#include "gromacs/ewald/pme.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "pme-types-ocl.h"
//#include "pme-timings.cuh"

//#include "pme-ocl-types-kernel.clh"
//CAN_USE_BUFFERS_IN_STRUCTS
//FIXME copy
constexpr bool c_canUseBuffersInStructs = (GMX_GPU != GMX_GPU_OPENCL); //&& (__OPENCL_C_VERSION__ <= 200))


void pme_gpu_gather(PmeGpu                *pmeGpu,
                    PmeForceOutputHandling forceTreatment,
                    const float           *h_grid
                    )
{
    /* Copying the input CPU forces for reduction */
    if (forceTreatment != PmeForceOutputHandling::Set)
    {
        pme_gpu_copy_input_forces(pmeGpu);
    }

    const int    order           = pmeGpu->common->pme_order;
    const auto  *kernelParamsPtr = pmeGpu->kernelParams.get();

    if (!pme_gpu_performs_FFT(pmeGpu) || pme_gpu_is_testing(pmeGpu))
    {
        pme_gpu_copy_input_gather_grid(pmeGpu, const_cast<float *>(h_grid));
    }

    if (pme_gpu_is_testing(pmeGpu))
    {
        pme_gpu_copy_input_gather_atom_data(pmeGpu);
    }

    const int atomsPerBlock  =  (c_gatherMaxThreadsPerBlock / PME_SPREADGATHER_THREADS_PER_ATOM);
    GMX_ASSERT(!c_usePadding || !(PME_ATOM_DATA_ALIGNMENT % atomsPerBlock), "inconsistent atom data padding vs. gathering block size");

    KernelLaunchConfig config;
    config.stream      = pmeGpu->archSpecific->pmeStream;

    const int blockCount = pmeGpu->nAtomsPadded / atomsPerBlock;
    config.blockSize.x = order;
    config.blockSize.y = order;
    config.blockSize.z = atomsPerBlock;
    config.gridSize.x = blockCount;

    //    //FIXME pickup Fermi fix pmeGpuCreateGrid(pmeGpu, blockCount);

    const bool wrapX = true;
    const bool wrapY = true;
    GMX_UNUSED_VALUE(wrapX);
    GMX_UNUSED_VALUE(wrapY);

    const auto &kernel = (forceTreatment == PmeForceOutputHandling::Set) ? pmeGpu->archSpecific->persistent->gatherKernel : pmeGpu->archSpecific->persistent->gatherReduceWithInputKernel;

    if (order != 4)
    {
        // just an assert? should be caught during init
        GMX_THROW(gmx::NotImplementedError("The code for pme_order != 4 is not implemented"));
    }
    pme_gpu_start_timing(pmeGpu, gtPME_GATHER);
    if (c_canUseBuffersInStructs)
    {
        launchGpuKernel(config, kernel, kernelParamsPtr);
    }
    else
    {
        launchGpuKernel(config, kernel, kernelParamsPtr,
                        &kernelParamsPtr->atoms.d_coefficients,
                        &kernelParamsPtr->grid.d_realGrid,
                        &kernelParamsPtr->atoms.d_theta,
                        &kernelParamsPtr->atoms.d_dtheta,
                        &kernelParamsPtr->atoms.d_gridlineIndices,
                        &kernelParamsPtr->atoms.d_forces
                        );
    }
    //CU_LAUNCH_ERR("pme_gather_kernel");
    pme_gpu_stop_timing(pmeGpu, gtPME_GATHER);

    pme_gpu_copy_output_forces(pmeGpu);
}
