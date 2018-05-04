#include "pme-persistent-data.h"
#include "pme-types-ocl.h"

#include "gromacs/gpu_utils/ocl_compiler.h"

PmeGpuPersistentData::PmeGpuPersistentData(PmeGpu *pmeGpu)
{
  printf("hello consstructor\n");

    cl_context_properties     context_properties[3];
    cl_platform_id            platform_id;
    cl_device_id              device_id;
    cl_int                    clError;

    platform_id      = pmeGpu->deviceInfo->ocl_gpu_id.ocl_platform_id;
    device_id        = pmeGpu->deviceInfo->ocl_gpu_id.ocl_device_id;

    context_properties[0] = CL_CONTEXT_PLATFORM;
    context_properties[1] = (cl_context_properties) platform_id;
    context_properties[2] = 0; /* Terminates the list of properties */

    context = clCreateContext(context_properties, 1, &device_id, nullptr, nullptr, &clError);
    throwUponFailure(clError);

     //GMX_RELEASE_ASSERT(CL_SUCCESS == clError, "whatever");
                     /*
                       gmx::formatString("Failed to create context for PME on GPU #%s:\n OpenCL error %d: %s",
                  pmeGpu->deviceInfo->device_name,
                  clError, ocl_get_error_string(clError).c_str())*/


    pme_gpu_compile_kernels(pmeGpu); //FIXME make a constructor
}

#if GMX_GPU == GMX_GPU_OPENCL

#include <string>
#include <vector>

// based on nbnxn_gpu_compile_kernels
void PmeGpuPersistentData::pme_gpu_compile_kernels(PmeGpu *pmeGpu)
{
    program  = nullptr;
    /* Need to catch std::bad_alloc here and during compilation string
       handling. */
    try
    {
        /* Here we pass macros and static const int variables defined in include
         * files outside the nbnxn_ocl as macros, to avoid including those files
         * in the JIT compilation that happens at runtime.
         */


        //FIXME thsi should be deduced with pme_gpu_copy_common_data_from()
        const int order = 4;

        const std::string spreadGatherDefines = gmx::formatString(
                    //All those are not needed for solve, but whatever
                    "-Dwarp_size=%d "
                    "-Dorder=%d "
                    "-DPME_SPREADGATHER_ATOMS_PER_WARP=%d "
                    "-DPME_SPREADGATHER_THREADS_PER_ATOM=%d "
                    "-DPME_SPLINE_THETA_STRIDE=%d "
                    "-Dc_usePadding=%d "
                    "-Dc_skipNeutralAtoms=%d "  //TODO stringify
                    ,
            warp_size,
                    order,
                    PME_SPREADGATHER_ATOMS_PER_WARP,
                    PME_SPREADGATHER_THREADS_PER_ATOM,
                    PME_SPLINE_THETA_STRIDE,
                    c_usePadding,
                    c_skipNeutralAtoms);
        const std::string spreadOnlyDefines = gmx::formatString(
                    //FIXME "-DatomsPerBlock=%d "
                    "-Dc_pmeMaxUnitcellShift=%f "
                    // unused template params for decomposition
                    "-DwrapX=true -DwrapY=true",
                    //c_spreadMaxThreadsPerBlock / PME_SPREADGATHER_THREADS_PER_ATOM,
                    static_cast<float>(c_pmeMaxUnitcellShift)
                    );
    /*
        const std::string gatherOnlyDefines = gmx::formatString(
                    "-DatomsPerBlock=%d "
                    // unused template params for decomposition
                    "-DwrapX=true -DwrapY=true",
                    c_gatherMaxThreadsPerBlock / PME_SPREADGATHER_THREADS_PER_ATOM
                    );
                    */



        const std::string spreadDefines = spreadGatherDefines + " " + spreadOnlyDefines;
        //const std::string gatherDefines = spreadGatherDefines + " " + gatherOnlyDefines;


#if 0
                gmx::formatString(
                    " -DCENTRAL=%d "
                    "-DNBNXN_GPU_NCLUSTER_PER_SUPERCLUSTER=%d -DNBNXN_GPU_CLUSTER_SIZE=%d -DNBNXN_GPU_JGROUP_SIZE=%d "
                    "-DGMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY=%d "
                    "-DNBNXN_MIN_RSQ=%s %s",
                    CENTRAL,                                                /* Defined in ishift.h */
                    c_nbnxnGpuNumClusterPerSupercluster,                    /* Defined in nbnxn_pairlist.h */
                    c_nbnxnGpuClusterSize,                                  /* Defined in nbnxn_pairlist.h */
                    c_nbnxnGpuJgroupSize,                                   /* Defined in nbnxn_pairlist.h */
                    getOclPruneKernelJ4Concurrency(nb->dev_info->vendor_e), /* In nbnxn_ocl_types.h  */
                    STRINGIFY_MACRO(NBNXN_MIN_RSQ)                          /* Defined in nbnxn_consts.h */
                                                                            /* NBNXN_MIN_RSQ passed as string to avoid
                                                                                floating point representation problems with sprintf */
                    , (nb->bPrefetchLjParam) ? "-DIATYPE_SHMEM" : ""
                    );
#endif

        try
        {
            /* TODO when we have a proper MPI-aware logging module,
               the log output here should be written there */
            program = gmx::ocl::compileProgram(stderr,
                                               "../../ewald/pme-program.cl", //FIXME
                                               spreadDefines,
                                               context,
                                               pmeGpu->deviceInfo->ocl_gpu_id.ocl_device_id,
                                               pmeGpu->deviceInfo->vendor_e);
        }
        catch (gmx::GromacsException &e)
        {
            e.prependContext(gmx::formatString("Failed to compile PME kernels for GPU #%s\n",
                                               pmeGpu->deviceInfo->device_name));
            throw;
        }
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    std::vector<cl_kernel> kernels;
    cl_uint justEnough = 9; //?
    kernels.resize(justEnough);
    cl_uint actualKernelCount = 0;
    cl_int status = clCreateKernelsInProgram(program, justEnough,
                          kernels.data(), &actualKernelCount);
    throwUponFailure(status);
    //fprintf(stderr, "got me soem kernels %u\n", actualKernelCount);
    kernels.resize(actualKernelCount);

    std::array<char, 100> kernelNamesBuffer;
    for (const auto &kernel: kernels)
    {
        status = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME,
                 kernelNamesBuffer.size(), kernelNamesBuffer.data(), nullptr);
    throwUponFailure(status);
    //FIXME fprintf(stderr, "got a nice kernel: %s\n", kernelNamesBuffer.data());
    if (!strcmp(kernelNamesBuffer.data(), "pmeSplineKernel"))
      splineKernel = kernel;
    if (!strcmp(kernelNamesBuffer.data(), "pmeSplineAndSpreadKernel"))
      splineAndSpreadKernel = kernel;
    if (!strcmp(kernelNamesBuffer.data(), "pmeSpreadKernel"))
      spreadKernel = kernel;
    if (!strcmp(kernelNamesBuffer.data(), "pmeGatherKernel"))
      gatherKernel = kernel;
    if (!strcmp(kernelNamesBuffer.data(), "pmeGatherReduceWithInputKernel"))
      gatherReduceWithInputKernel = kernel;
    if (!strcmp(kernelNamesBuffer.data(), "pmeSolveYZXKernel"))
      solveYZXKernel = kernel;
    if (!strcmp(kernelNamesBuffer.data(), "pmeSolveYZXEnergyKernel"))
      solveYZXEnergyKernel = kernel;
    if (!strcmp(kernelNamesBuffer.data(), "pmeSolveXYZKernel"))
      solveXYZKernel = kernel;
    if (!strcmp(kernelNamesBuffer.data(), "pmeSolveXYZEnergyKernel"))
      solveXYZEnergyKernel = kernel;
    throwUponFailure(status);
    }

    //TODO put those guys in a map with string name as a key
}
#endif

PmeGpuPersistentData::~PmeGpuPersistentData()
{
    printf("die die die\n");
    clReleaseKernel(splineAndSpreadKernel);
    clReleaseKernel(splineKernel);
    clReleaseKernel(spreadKernel);
    clReleaseKernel(gatherKernel);
    clReleaseKernel(gatherReduceWithInputKernel);
    clReleaseKernel(solveXYZKernel);
    clReleaseKernel(solveXYZEnergyKernel);
    clReleaseKernel(solveYZXKernel);
    clReleaseKernel(solveYZXEnergyKernel);
    clReleaseProgram(program);
    clReleaseContext(context);
    //throwUponFailure(clError);
    //GMX_RELEASE_ASSERT(clError == CL_SUCCESS, "PME context destruction error");
}
