#ifndef PMEOCLDEFINITELYCOMMON_H
#define PMEOCLDEFINITELYCOMMON_H

/*
inline int actualExecutionWidth()
{
    return 32;
}
*/

//#define warp_size
//actualExecutionWidth()

#if defined __OPENCL_C_VERSION__
#define OPENCL_COMPILATION 1
//FIXME duplicate, rename to OPENCL_DEVICE
#else
#define OPENCL_COMPILATION 0
#endif

#define USE_C99_ONLY (OPENCL_COMPILATION && (__OPENCL_C_VERSION__ <= 200))

//FIXME check if macros are duplicate, and use them correctly
#define CAN_USE_BUFFERS_IN_STRUCTS !USE_C99_ONLY

#define c_solveMaxWarpsPerBlock 8
//! Solving kernel max block size in threads
#define c_solveMaxThreadsPerBlock (c_solveMaxWarpsPerBlock * warp_size)

#define c_virialAndEnergyCount 7

// copy of all kernel param structures (non-OCL ones)
// All deviceBuffers are shielded. Again: check if 2.x is viable with those at all.

//FIXME
#define POINTER_SIZE 8

#ifndef DIM
#define DIM 3
#endif

#define PACKED  __attribute__ ((aligned(POINTER_SIZE)))

/*! \internal \brief
 * A GPU data structure for storing the constant PME data.
 * This only has to be initialized once.
 */
struct PACKED PmeGpuConstParams
{
    /*! \brief Electrostatics coefficient = ONE_4PI_EPS0 / pme->epsilon_r */
    float elFactor;
#if CAN_USE_BUFFERS_IN_STRUCTS
    /*! \brief Virial and energy GPU array. Size is PME_GPU_ENERGY_AND_VIRIAL_COUNT (7) floats.
     * The element order is virxx, viryy, virzz, virxy, virxz, viryz, energy. */
    DeviceBuffer<float> d_virialAndEnergy;
#else
    char dummy[1 * POINTER_SIZE];    
#endif
};

/*! \internal \brief
 * A GPU data structure for storing the PME data related to the grid sizes and cut-off.
 * This only has to be updated at every DD step.
 */
struct PACKED PmeGpuGridParams
{
    /* Grid sizes */
    /*! \brief Real-space grid data dimensions. */
    int realGridSize[DIM];
    /*! \brief Real-space grid dimensions, only converted to floating point. */
    float realGridSizeFP[DIM];
    /*! \brief Real-space grid dimensions (padded). The padding as compared to realGridSize includes the (order - 1) overlap. */
    int realGridSizePadded[DIM]; /* Is major dimension of this ever used in kernels? */
    /*! \brief Fourier grid dimensions. This counts the complex numbers! */
    int complexGridSize[DIM];
    /*! \brief Fourier grid dimensions (padded). This counts the complex numbers! */
    int complexGridSizePadded[DIM];
    /*! \brief Ewald solving factor = (M_PI / pme->ewaldcoeff_q)^2 */
    float ewaldFactor;
    /*! \brief Offsets for X/Y/Z components of d_splineModuli */
    int splineValuesOffset[DIM];
    /*! \brief Offsets for X/Y/Z components of d_fractShiftsTable and d_gridlineIndicesTable */
    int tablesOffsets[DIM];
#if CAN_USE_BUFFERS_IN_STRUCTS
    /* Grid pointers */
    /*! \brief Real space grid. */
    DeviceBuffer<float> d_realGrid;
    /*! \brief Complex grid - used in FFT/solve. If inplace cuFFT is used, then it is the same pointer as realGrid. */
    DeviceBuffer<float> d_fourierGrid;
    /*! \brief Grid spline values as in pme->bsp_mod
     * (laid out sequentially (XXX....XYYY......YZZZ.....Z))
     */
    DeviceBuffer<float> d_splineModuli;
    /*! \brief Fractional shifts lookup table as in pme->fshx/fshy/fshz, laid out sequentially (XXX....XYYY......YZZZ.....Z) */
    DeviceBuffer<float> d_fractShiftsTable;
    /*! \brief Gridline indices lookup table
     * (modulo lookup table as in pme->nnx/nny/nnz, laid out sequentially (XXX....XYYY......YZZZ.....Z)) */
    DeviceBuffer<int> d_gridlineIndicesTable;
#else
    char dummy[5 * POINTER_SIZE];
#endif
};

/*! \internal \brief
 * A GPU data structure for storing the PME data of the atoms, local to this process' domain partition.
 * This only has to be updated every DD step.
 */
struct PACKED PmeGpuAtomParams
{
    /*! \brief Number of local atoms */
    int    nAtoms;
    /*! \brief Pointer to the global GPU memory with input rvec atom coordinates.
     * The coordinates themselves change and need to be copied to the GPU for every PME computation,
     * but reallocation happens only at DD.
     */
#if CAN_USE_BUFFERS_IN_STRUCTS
    DeviceBuffer<float> d_coordinates;
    /*! \brief Pointer to the global GPU memory with input atom charges.
     * The charges only need to be reallocated and copied to the GPU at DD step.
     */
    DeviceBuffer<float> d_coefficients;
    /*! \brief Pointer to the global GPU memory with input/output rvec atom forces.
     * The forces change and need to be copied from (and possibly to) the GPU for every PME computation,
     * but reallocation happens only at DD.
     */
    DeviceBuffer<float> d_forces;
    /*! \brief Pointer to the global GPU memory with ivec atom gridline indices.
     * Computed on GPU in the spline calculation part.
     */
    DeviceBuffer<int> d_gridlineIndices;

    /* B-spline parameters are computed entirely on GPU for every PME computation, not copied.
     * Unless we want to try something like GPU spread + CPU gather?
     */
    /*! \brief Pointer to the global GPU memory with B-spline values */
    DeviceBuffer<float> d_theta;
    /*! \brief Pointer to the global GPU memory with B-spline derivative values */
    DeviceBuffer<float> d_dtheta;
#else
    char dummy[6 * POINTER_SIZE];
#endif
};

/*! \internal \brief
 * A GPU data structure for storing the PME data which might change for each new PME computation.
 */
struct PACKED PmeGpuDynamicParams
{
    /* The box parameters. The box only changes size with pressure coupling enabled. */
    /*! \brief
     * Reciprocal (inverted unit cell) box.
     *
     * The box is transposed as compared to the CPU pme->recipbox.
     * Basically, spread uses matrix columns (while solve and gather use rows).
     * This storage format might be not the most optimal since the box is always triangular so there are zeroes.
     */
    float  recipBox[DIM][DIM];
    /*! \brief The unit cell volume for solving. */
    float  boxVolume;
};

/*! \internal \brief
 * A single structure encompassing almost all the PME data used in GPU kernels on device.
 * This is inherited by the GPU framework-specific structure
 * (PmeGpuCudaKernelParams in pme.cuh).
 * This way, most code preparing the kernel parameters can be GPU-agnostic by casting
 * the kernel parameter data pointer to PmeGpuKernelParamsBase.
 */
struct PACKED PmeGpuCudaKernelParams
{
    /*! \brief Constant data that is set once. */
    struct PmeGpuConstParams   constants;
    /*! \brief Data dependent on the grid size/cutoff. */
    struct PmeGpuGridParams    grid;
    /*! \brief Data dependent on the DD and local atoms. */
    struct PmeGpuAtomParams  atoms;
    /*! \brief Data that possibly changes for every new PME computation.
     * This should be kept up-to-date by calling pme_gpu_prepare_computation(...)
     * before launching spreading.
     */
    struct PmeGpuDynamicParams current;
    //FIXME this is criminal
    int fractShiftsTableTexture;
    int gridlineIndicesTableTexture;

};

#endif // PMEOCLDEFINITELYCOMMON_H
