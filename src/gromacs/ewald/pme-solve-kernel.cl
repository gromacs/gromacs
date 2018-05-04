#include "../../ewald/pme-ocl-types-kernel.clh"

/*! \brief
 * PME complex grid solver kernel function.
 *
 * \tparam[in] gridOrdering             Specifies the dimension ordering of the complex grid.
 * \tparam[in] computeEnergyAndVirial   Tells if the reciprocal energy and virial should be computed.
 * \param[in]  kernelParams             Input PME CUDA data in constant memory.
 */
/*
template<
    GridOrderingInternal gridOrdering,
    bool computeEnergyAndVirial
    >
__launch_bounds__(c_solveMaxThreadsPerBlock)
*/
KERNEL_FUNC void CUSTOMIZED_KERNEL_NAME(pme_solve_kernel)(const struct PmeGpuCudaKernelParams kernelParams
#if !CAN_USE_BUFFERS_IN_STRUCTS
                    ,
                    GLOBAL const float * __restrict__ gm_splineModuli,
                    GLOBAL float * __restrict__       gm_virialAndEnergy,
                    GLOBAL float2 * __restrict__      gm_grid
#endif
)
{
    /* This kernel supports 2 different grid dimension orderings: YZX and XYZ */
    int majorDim, middleDim, minorDim;
    if (gridOrdering == YZX)
    {
            majorDim  = YY;
            middleDim = ZZ;
            minorDim  = XX;
    }
        if (gridOrdering == XYZ)
        {
            majorDim  = XX;
            middleDim = YY;
            minorDim  = ZZ;
    }
#if CAN_USE_BUFFERS_IN_STRUCTS
    /* Global memory pointers */
    GLOBAL const float * __restrict__ gm_splineModuli       = kernelParams.grid.d_splineModuli;
    GLOBAL float * __restrict__       gm_virialAndEnergy    = kernelParams.constants.d_virialAndEnergy;
    GLOBAL float2 * __restrict__      gm_grid               = (float2 *)kernelParams.grid.d_fourierGrid;
#endif

    GLOBAL const float * __restrict__ gm_splineValueMajor   = gm_splineModuli + kernelParams.grid.splineValuesOffset[majorDim];
    GLOBAL const float * __restrict__ gm_splineValueMiddle  = gm_splineModuli + kernelParams.grid.splineValuesOffset[middleDim];
    GLOBAL const float * __restrict__ gm_splineValueMinor   = gm_splineModuli + kernelParams.grid.splineValuesOffset[minorDim];




    /* Various grid sizes and indices */
    const int localOffsetMinor = 0, localOffsetMajor = 0, localOffsetMiddle = 0; //unused
    const int localSizeMinor   = kernelParams.grid.complexGridSizePadded[minorDim];
    const int localSizeMiddle  = kernelParams.grid.complexGridSizePadded[middleDim];
    const int localCountMiddle = kernelParams.grid.complexGridSize[middleDim];
    const int localCountMinor  = kernelParams.grid.complexGridSize[minorDim];
    const int nMajor           = kernelParams.grid.realGridSize[majorDim];
    const int nMiddle          = kernelParams.grid.realGridSize[middleDim];
    const int nMinor           = kernelParams.grid.realGridSize[minorDim];
    const int maxkMajor        = (nMajor + 1) / 2;  // X or Y
    const int maxkMiddle       = (nMiddle + 1) / 2; // Y OR Z => only check for !YZX
    const int maxkMinor        = (nMinor + 1) / 2;  // Z or X => only check for YZX

    /* Each thread works on one cell of the Fourier space complex 3D grid (gm_grid).
     * Each block handles up to c_solveMaxThreadsPerBlock cells -
     * depending on the grid contiguous dimension size,
     * that can range from a part of a single gridline to several complete gridlines.
     */
    const int threadLocalId     = getThreadLocalIndex(XX);
    const int gridLineSize      = localCountMinor;
    const int gridLineIndex     = threadLocalId / gridLineSize;
    const int gridLineCellIndex = threadLocalId - gridLineSize * gridLineIndex;
    const int gridLinesPerBlock = getBlockSize(XX) / gridLineSize;
    const int activeWarps       = (getBlockSize(XX) / warp_size);
    const int indexMinor        = getBlockIndex(XX) * getBlockSize(XX) + gridLineCellIndex;
    const int indexMiddle       = getBlockIndex(YY) * gridLinesPerBlock + gridLineIndex;
    const int indexMajor        = getBlockIndex(ZZ);

    /* Optional outputs */
    float energy = 0.0f;
    float virxx  = 0.0f;
    float virxy  = 0.0f;
    float virxz  = 0.0f;
    float viryy  = 0.0f;
    float viryz  = 0.0f;
    float virzz  = 0.0f;

    assert(indexMajor < kernelParams.grid.complexGridSize[majorDim]);
    if ((indexMiddle < localCountMiddle) & (indexMinor < localCountMinor) & (gridLineIndex < gridLinesPerBlock))
    {
        /* The offset should be equal to the global thread index for coalesced access */
        const int             gridIndex     = (indexMajor * localSizeMiddle + indexMiddle) * localSizeMinor + indexMinor;
        GLOBAL float2 * __restrict__ gm_gridCell   = gm_grid + gridIndex;

        const int             kMajor  = indexMajor + localOffsetMajor;
        /* Checking either X in XYZ, or Y in YZX cases */
        const float           mMajor  = (kMajor < maxkMajor) ? kMajor : (kMajor - nMajor);

        const int             kMiddle = indexMiddle + localOffsetMiddle;
        float                 mMiddle = kMiddle;
        /* Checking Y in XYZ case */
        if (gridOrdering == XYZ)
        {
            mMiddle = (kMiddle < maxkMiddle) ? kMiddle : (kMiddle - nMiddle);
        }
        const int             kMinor  = localOffsetMinor + indexMinor;
        float                 mMinor  = kMinor;
        /* Checking X in YZX case */
        if (gridOrdering == YZX)
        {
            mMinor = (kMinor < maxkMinor) ? kMinor : (kMinor - nMinor);
        }
        /* We should skip the k-space point (0,0,0) */
        const bool notZeroPoint  = (kMinor > 0) | (kMajor > 0) | (kMiddle > 0);

        float      mX, mY, mZ;
        if (gridOrdering == YZX)
        {
                mX = mMinor;
                mY = mMajor;
                mZ = mMiddle;
        }
        if (gridOrdering == XYZ)
        {
            mX = mMajor;
            mY = mMiddle;
            mZ = mMinor;
        }

        /* 0.5 correction factor for the first and last components of a Z dimension */
        float corner_fac = 1.0f;
        if (gridOrdering == YZX)
        {
                if ((kMiddle == 0) | (kMiddle == maxkMiddle))
                {
                    corner_fac = 0.5f;
                }
        }
        if (gridOrdering == XYZ)
        {
                if ((kMinor == 0) | (kMinor == maxkMinor))
                {
                    corner_fac = 0.5f;
                }
        }

        if (notZeroPoint)
        {
            const float mhxk = mX * kernelParams.current.recipBox[XX][XX];
            const float mhyk = mX * kernelParams.current.recipBox[XX][YY] + mY * kernelParams.current.recipBox[YY][YY];
            const float mhzk = mX * kernelParams.current.recipBox[XX][ZZ] + mY * kernelParams.current.recipBox[YY][ZZ] + mZ * kernelParams.current.recipBox[ZZ][ZZ];

            const float m2k        = mhxk * mhxk + mhyk * mhyk + mhzk * mhzk;
            assert(m2k != 0.0f);
            //TODO: use LDG/textures for gm_splineValue
            float       denom = m2k * M_PI_F * kernelParams.current.boxVolume * gm_splineValueMajor[kMajor] * gm_splineValueMiddle[kMiddle] * gm_splineValueMinor[kMinor];
            assert(isfinite(denom));
            assert(denom != 0.0f);
            const float   tmp1   = exp(-kernelParams.grid.ewaldFactor * m2k); //FIXME was expf in CUDA
            const float   etermk = kernelParams.constants.elFactor * tmp1 / denom;

            float2        gridValue    = *gm_gridCell;
            const float2  oldGridValue = gridValue;

            bool valueFilter = (fabs(gridValue.x) + fabs(gridValue.y) >= 1e-5);
            //bool indexFilter = (gridIndex < 500);
            bool indexFilter = ((int)mZ == 5) && ((int)mY == 7);
            if (valueFilter && indexFilter)
                ;//printf("grid %d %d %d %f %f\n", (int)mX, (int)mY, (int)mZ, gridValue.x, gridValue.y);
            gridValue.x   *= etermk;
            gridValue.y   *= etermk;
            *gm_gridCell   = gridValue;

            if (computeEnergyAndVirial)
            {
                const float tmp1k = 2.0f * (gridValue.x * oldGridValue.x + gridValue.y * oldGridValue.y);

                float       vfactor = (kernelParams.grid.ewaldFactor + 1.0f / m2k) * 2.0f;
                float       ets2    = corner_fac * tmp1k;
                energy = ets2;

                float ets2vf  = ets2 * vfactor;

                virxx   = ets2vf * mhxk * mhxk - ets2;
                virxy   = ets2vf * mhxk * mhyk;
                virxz   = ets2vf * mhxk * mhzk;
                viryy   = ets2vf * mhyk * mhyk - ets2;
                viryz   = ets2vf * mhyk * mhzk;
                virzz   = ets2vf * mhzk * mhzk - ets2;
            }
        }
    }

    //FIXME thsi is onl for reduction though
    SHARED float sm_virialAndEnergy[c_virialAndEnergyCount * warp_size];

    /* Optional energy/virial reduction */
    if (computeEnergyAndVirial)
    {
#if (GMX_PTX_ARCH >= 300)  && 0 //FIXME
        /* A tricky shuffle reduction inspired by reduce_force_j_warp_shfl.
         * The idea is to reduce 7 energy/virial components into a single variable (aligned by 8).
         * We will reduce everything into virxx.
         */

        /* We can only reduce warp-wise */
        const int          width      = warp_size;
        const unsigned int activeMask = c_fullWarpMask;

        /* Making pair sums */
        virxx  += gmx_shfl_down_sync(activeMask, virxx, 1, width);
        viryy  += gmx_shfl_up_sync  (activeMask, viryy, 1, width);
        virzz  += gmx_shfl_down_sync(activeMask, virzz, 1, width);
        virxy  += gmx_shfl_up_sync  (activeMask, virxy, 1, width);
        virxz  += gmx_shfl_down_sync(activeMask, virxz, 1, width);
        viryz  += gmx_shfl_up_sync  (activeMask, viryz, 1, width);
        energy += gmx_shfl_down_sync(activeMask, energy, 1, width);
        if (threadLocalId & 1)
        {
            virxx = viryy; // virxx now holds virxx and viryy pair sums
            virzz = virxy; // virzz now holds virzz and virxy pair sums
            virxz = viryz; // virxz now holds virxz and viryz pair sums
        }

        /* Making quad sums */
        virxx  += gmx_shfl_down_sync(activeMask, virxx, 2, width);
        virzz  += gmx_shfl_up_sync  (activeMask, virzz, 2, width);
        virxz  += gmx_shfl_down_sync(activeMask, virxz, 2, width);
        energy += gmx_shfl_up_sync  (activeMask, energy, 2, width);
        if (threadLocalId & 2)
        {
            virxx = virzz;  // virxx now holds quad sums of virxx, virxy, virzz and virxy
            virxz = energy; // virxz now holds quad sums of virxz, viryz, energy and unused paddings
        }

        /* Making octet sums */
        virxx += gmx_shfl_down_sync(activeMask, virxx, 4, width);
        virxz += gmx_shfl_up_sync  (activeMask, virxz, 4, width);
        if (threadLocalId & 4)
        {
            virxx = virxz; // virxx now holds all 7 components' octet sums + unused paddings
        }

        /* We only need to reduce virxx now */
#pragma unroll
        for (int delta = 8; delta < width; delta <<= 1)
        {
            virxx += gmx_shfl_down_sync(activeMask, virxx, delta, width);
        }
        /* Now first 7 threads of each warp have the full output contributions in virxx */

        const int        componentIndex      = threadLocalId & (warp_size - 1);
        const bool       validComponentIndex = (componentIndex < c_virialAndEnergyCount);
        /* Reduce 7 outputs per warp in the shared memory */
        const int        stride              = 8; // this is c_virialAndEnergyCount==7 rounded up to power of 2 for convenience, hence the assert
        assert(c_virialAndEnergyCount == 7);
        const int        reductionBufferSize = (c_solveMaxThreadsPerBlock / warp_size) * stride;
        SHARED float sm_virialAndEnergy[reductionBufferSize];

        if (validComponentIndex)
        {
            const int warpIndex = threadLocalId / warp_size;
            sm_virialAndEnergy[warpIndex * stride + componentIndex] = virxx;
        }
        sharedMemoryBarrier();

        /* Reduce to the single warp size */
        const int targetIndex = threadLocalId;
#pragma unroll
        for (int reductionStride = reductionBufferSize >> 1; reductionStride >= warp_size; reductionStride >>= 1)
        {
            const int sourceIndex = targetIndex + reductionStride;
            if ((targetIndex < reductionStride) & (sourceIndex < activeWarps * stride))
            {
                // TODO: the second conditional is only needed on first iteration, actually - see if compiler eliminates it!
                sm_virialAndEnergy[targetIndex] += sm_virialAndEnergy[sourceIndex];
            }
            __syncthreads();
        }

        /* Now use shuffle again */
        if (threadLocalId < warp_size)
        {
            float output = sm_virialAndEnergy[threadLocalId];
#pragma unroll
            for (int delta = stride; delta < warp_size; delta <<= 1)
            {
                output += gmx_shfl_down_sync(activeMask, output, delta, warp_size);
            }
            /* Final output */
            if (validComponentIndex)
            {
                assert(isfinite(output));
                atomicAdd(gm_virialAndEnergy + componentIndex, output);
            }
        }
#else



        /* Shared memory reduction with atomics for compute capability < 3.0.
         * Each component is first reduced into warp_size positions in the shared memory;
         * Then first c_virialAndEnergyCount warps reduce everything further and add to the global memory.
         * This can likely be improved, but is anyway faster than the previous straightforward reduction,
         * which was using too much shared memory (for storing all 7 floats on each thread).
         * [48KB (shared mem limit per SM on CC2.x) / sizeof(float) (4) / c_solveMaxThreadsPerBlock (256) / c_virialAndEnergyCount (7) ==
         * 6 blocks per SM instead of 16 which is maximum on CC2.x].
         */

        const int        lane      = threadLocalId & (warp_size - 1);
        const int        warpIndex = threadLocalId / warp_size;
        const bool       firstWarp = (warpIndex == 0);
        if (firstWarp)
        {
            sm_virialAndEnergy[0 * warp_size + lane] = virxx;
            sm_virialAndEnergy[1 * warp_size + lane] = viryy;
            sm_virialAndEnergy[2 * warp_size + lane] = virzz;
            sm_virialAndEnergy[3 * warp_size + lane] = virxy;
            sm_virialAndEnergy[4 * warp_size + lane] = virxz;
            sm_virialAndEnergy[5 * warp_size + lane] = viryz;
            sm_virialAndEnergy[6 * warp_size + lane] = energy;
        }
        sharedMemoryBarrier();
        if (!firstWarp)
        {
            atomicAdd_l_f(sm_virialAndEnergy + 0 * warp_size + lane, virxx);
            atomicAdd_l_f(sm_virialAndEnergy + 1 * warp_size + lane, viryy);
            atomicAdd_l_f(sm_virialAndEnergy + 2 * warp_size + lane, virzz);
            atomicAdd_l_f(sm_virialAndEnergy + 3 * warp_size + lane, virxy);
            atomicAdd_l_f(sm_virialAndEnergy + 4 * warp_size + lane, virxz);
            atomicAdd_l_f(sm_virialAndEnergy + 5 * warp_size + lane, viryz);
            atomicAdd_l_f(sm_virialAndEnergy + 6 * warp_size + lane, energy);
        }
        sharedMemoryBarrier();

        //FIXME make multiple iterations for wide warps
        assert(activeWarps >= c_virialAndEnergyCount); // we need to cover all components, or have multiple iterations otherwise
        const int componentIndex = warpIndex;
        if (componentIndex < c_virialAndEnergyCount)
        {
            const int targetIndex = threadLocalId;
#pragma unroll
            for (int reductionStride = warp_size >> 1; reductionStride >= 1; reductionStride >>= 1)
            {
                if (lane < reductionStride)
                {
                    sm_virialAndEnergy[targetIndex] += sm_virialAndEnergy[targetIndex + reductionStride];
                    //FIXME thsi added barrrier makes OpenCL correct
                    sharedMemoryBarrier();
                }
            }
            if (lane == 0)
            {
                atomicAdd(gm_virialAndEnergy + componentIndex, sm_virialAndEnergy[targetIndex]);
            }
        }
#endif
    }
}
