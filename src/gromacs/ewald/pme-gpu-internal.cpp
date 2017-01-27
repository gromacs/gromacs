/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

/*! \internal \file
 *
 * \brief This file contains internal function implementations
 * for performing the PME calculations on GPU.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include "pme-gpu-internal.h"

#include "config.h"

#include <cstdlib>
#include <cstring>

#include <string>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "pme-grid.h"
#include "pme-internal.h"

/*! \internal \brief
 * Wrapper for getting a pointer to the plain C++ part of the GPU kernel parameters structure.
 *
 * \param[in] pmeGPU  The PME GPU structure.
 * \returns The pointer to the kernel parameters.
 */
pme_gpu_kernel_params_base_t *pme_gpu_get_kernel_params_base_ptr(const pme_gpu_t *pmeGPU)
{
    // reinterpret_cast is needed because the derived CUDA structure is not known in this file
    pme_gpu_kernel_params_base_t *kernelParamsPtr = reinterpret_cast<pme_gpu_kernel_params_base_t *>(pmeGPU->kernelParams.get());
    return kernelParamsPtr;
}

void pme_gpu_get_energy_virial(const pme_gpu_t *pmeGPU, real *energy, matrix virial)
{
    GMX_ASSERT(energy, "Bad (NULL) energy output in PME GPU");
    size_t j = 0;
    virial[XX][XX] = 0.25f * pmeGPU->staging.h_virialAndEnergy[j++];
    virial[YY][YY] = 0.25f * pmeGPU->staging.h_virialAndEnergy[j++];
    virial[ZZ][ZZ] = 0.25f * pmeGPU->staging.h_virialAndEnergy[j++];
    virial[XX][YY] = virial[YY][XX] = 0.25f * pmeGPU->staging.h_virialAndEnergy[j++];
    virial[XX][ZZ] = virial[ZZ][XX] = 0.25f * pmeGPU->staging.h_virialAndEnergy[j++];
    virial[YY][ZZ] = virial[ZZ][YY] = 0.25f * pmeGPU->staging.h_virialAndEnergy[j++];
    *energy        = 0.5f * pmeGPU->staging.h_virialAndEnergy[j++];
}

void pme_gpu_start_step(pme_gpu_t *pmeGPU, const matrix box, const rvec *h_coordinates)
{
    pme_gpu_copy_input_coordinates(pmeGPU, h_coordinates);

    pme_gpu_kernel_params_base_t *kernelParamsPtr = pme_gpu_get_kernel_params_base_ptr(pmeGPU);

    const size_t                  boxMemorySize        = sizeof(matrix);
    const bool                    haveToUpdateUnitCell = memcmp(pmeGPU->previousBox, box, boxMemorySize);
    /* There could be a pressure coupling check here, but this is more straightforward.
     * This is an exact comparison of float values though.
     */
    if (haveToUpdateUnitCell)
    {
        memcpy(pmeGPU->previousBox, box, boxMemorySize);

        kernelParamsPtr->step.boxVolume = box[XX][XX] * box[YY][YY] * box[ZZ][ZZ];
        GMX_ASSERT(kernelParamsPtr->step.boxVolume != 0.0f, "Zero volume of the unit cell");

#if GMX_DOUBLE
        GMX_THROW(gmx::NotImplementedError("PME is implemented for single-precision only on GPU"));
#else
        matrix recipBox;
        gmx::invertBoxMatrix(box, recipBox);
        /* The GPU recipBox is transposed as compared to the CPU recipBox.
         * Spread uses matrix columns (while solve and gather use rows).
         * There is no particular reason for this; it might be further rethought/optimized for better access patterns.
         */
        const real newRecipBox[DIM][DIM] =
        {
            {recipBox[XX][XX], recipBox[YY][XX], recipBox[ZZ][XX]},
            {             0.0, recipBox[YY][YY], recipBox[ZZ][YY]},
            {             0.0,              0.0, recipBox[ZZ][ZZ]}
        };
        memcpy(kernelParamsPtr->step.recipBox, newRecipBox, boxMemorySize);
#endif
    }
}

/*! \brief \libinternal
 * The PME GPU reinitialization function that is called both at the end of any MD step and on any load balancing step.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
void pme_gpu_reinit_step(const pme_gpu_t *pmeGPU)
{
    pme_gpu_clear_grids(pmeGPU);
    pme_gpu_clear_energy_virial(pmeGPU);
}

void pme_gpu_finish_step(const pme_gpu_t *pmeGPU, const bool bCalcF, const bool bCalcEnerVir)
{
    /* Needed for copy back as well as timing events */
    pme_gpu_synchronize(pmeGPU);

    if (bCalcF)
    {
        pme_gpu_sync_output_forces(pmeGPU);
    }
    if (bCalcEnerVir)
    {
        pme_gpu_sync_output_energy_virial(pmeGPU);
    }
    pme_gpu_update_timings(pmeGPU);
    pme_gpu_reinit_step(pmeGPU);
}

/*! \brief \libinternal
 * (Re-)initializes all the PME GPU data related to the grid size and cut-off.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
void pme_gpu_reinit_grids(pme_gpu_t *pmeGPU)
{
    pme_gpu_kernel_params_base_t *kernelParamsPtr = pme_gpu_get_kernel_params_base_ptr(pmeGPU);
    kernelParamsPtr->grid.ewaldFactor = (M_PI * M_PI) / (pmeGPU->common->ewaldcoeff_q * pmeGPU->common->ewaldcoeff_q);

    /* The grid size variants */
    for (int i = 0; i < DIM; i++)
    {
        kernelParamsPtr->grid.realGridSize[i]       = pmeGPU->common->nk[i];
        kernelParamsPtr->grid.realGridSizeFP[i]     = (float)kernelParamsPtr->grid.realGridSize[i];
        kernelParamsPtr->grid.realGridSizePadded[i] = pmeGPU->common->pmegrid_n[i];

        kernelParamsPtr->grid.complexGridSize[i]       = kernelParamsPtr->grid.realGridSize[i];
        kernelParamsPtr->grid.complexGridSizePadded[i] = kernelParamsPtr->grid.realGridSizePadded[i];
    }
    /* FFT: n real elements correspond to (n / 2 + 1) complex elements in minor dimension */
    kernelParamsPtr->grid.complexGridSize[ZZ] /= 2;
    kernelParamsPtr->grid.complexGridSize[ZZ]++;
    kernelParamsPtr->grid.complexGridSizePadded[ZZ] = kernelParamsPtr->grid.complexGridSize[ZZ];
    //TODO: should this take the X dimension for CPU FFT?

    pme_gpu_realloc_and_copy_fract_shifts(pmeGPU);
    pme_gpu_realloc_and_copy_bspline_values(pmeGPU);
    pme_gpu_realloc_grids(pmeGPU);
    pme_gpu_reinit_3dfft(pmeGPU);
}

/* Several GPU functions that refer to the CPU PME data live here.
 * We would like to keep these away from the GPU-framework specific code for clarity,
 * as well as compilation issues with MPI.
 */

/*! \brief \libinternal
 * Copies everything useful from the PME CPU to the PME GPU structure.
 * The goal is to minimize interaction with the PME CPU structure in the GPU code.
 *
 * \param[in] pme         The PME structure.
 */
void pme_gpu_copy_common_data_from(const gmx_pme_t *pme)
{
    /* TODO: Consider refactoring the CPU PME code to use the same structure,
     * so that this function becomes 2 lines */
    pme_gpu_t *pmeGPU             = pme->gpu;
    pmeGPU->common->ngrids        = pme->ngrids;
    pmeGPU->common->epsilon_r     = pme->epsilon_r;
    pmeGPU->common->ewaldcoeff_q  = pme->ewaldcoeff_q;
    pmeGPU->common->nk[XX]        = pme->nkx;
    pmeGPU->common->nk[YY]        = pme->nky;
    pmeGPU->common->nk[ZZ]        = pme->nkz;
    pmeGPU->common->pmegrid_n[XX] = pme->pmegrid_nx;
    pmeGPU->common->pmegrid_n[YY] = pme->pmegrid_ny;
    pmeGPU->common->pmegrid_n[ZZ] = pme->pmegrid_nz;
    pmeGPU->common->pme_order     = pme->pme_order;
    for (int i = 0; i < DIM; i++)
    {
        pmeGPU->common->bsp_mod[i].assign(pme->bsp_mod[i], pme->bsp_mod[i] + pmeGPU->common->nk[i]);
    }
    const int cellCount = c_pmeNeighborUnitcellCount;
    pmeGPU->common->fsh[XX].assign(pme->fshx, pme->fshx + cellCount * pme->nkx);
    pmeGPU->common->fsh[YY].assign(pme->fshy, pme->fshy + cellCount * pme->nky);
    pmeGPU->common->fsh[ZZ].assign(pme->fshz, pme->fshz + cellCount * pme->nkz);
    pmeGPU->common->nn[XX].assign(pme->nnx, pme->nnx + cellCount * pme->nkx);
    pmeGPU->common->nn[YY].assign(pme->nny, pme->nny + cellCount * pme->nky);
    pmeGPU->common->nn[ZZ].assign(pme->nnz, pme->nnz + cellCount * pme->nkz);
}

/*! \brief \libinternal
 * Finds out if PME with given inputs is possible to run on GPU.
 *
 * \param[in]  pme          The PME structure.
 * \param[out] error        The error message if the input is not supported on GPU.
 * \returns                 True if this PME input is possible to run on GPU, false otherwise.
 */
bool pme_gpu_check_restrictions(const gmx_pme_t *pme,
                                std::string     &error)
{
    if (pme->nnodes != 1)
    {
        error.append("PME is only implemented for a single rank on GPU. ");
    }
    if (pme->pme_order != 4)
    {
        error.append("PME is only implemented for the interpolation order of 4 on GPU. ");
    }
    if (pme->bFEP)
    {
        error.append("PME is only implemented for a single grid on GPU. ");
    }
    if (pme->doLJ)
    {
        error.append("PME LJ is not implemented on GPU. ");
    }
#if GMX_DOUBLE
    {
        error.append("PME is only implemented for single precision on GPU. ");
    }
#endif
#if GMX_GPU != GMX_GPU_CUDA
    {
        error.append("PME is only implemented for CUDA framework on GPU. ");
    }
#endif

    return error.empty();
}

/*! \libinternal \brief
 * Initializes the PME GPU data at the beginning of the run.
 *
 * \param[in] pme     The PME structure.
 * \param[in] hwinfo  The hardware information structure.
 * \param[in] gpu_opt The GPU information structure.
 */
void pme_gpu_init(gmx_pme_t *pme, const gmx_hw_info_t *hwinfo, const gmx_gpu_opt_t *gpu_opt)
{
    std::string error;
    pme->useGPU = pme_gpu_check_restrictions(pme, error);
    if (!pme->useGPU)
    {
        GMX_THROW(gmx::NotImplementedError(error));
    }

    pme->gpu          = new pme_gpu_t();
    pme_gpu_t *pmeGPU = pme->gpu;
    pmeGPU->common = std::shared_ptr<pme_shared_t>(new pme_shared_t());

    /* Some settings are set here at once for the whole run */
    /* A convenience variable. */
    pmeGPU->settings.useDecomposition = (pme->nnodes == 1);
    /* GPU FFT will only get used for a single rank. */
    pmeGPU->settings.performGPUFFT = !pme_gpu_uses_dd(pmeGPU) && !getenv("GMX_PME_GPU_FFTW");
    /* FIXME - CPU solve/FFTW after the GPU spread is definitely broken at the moment - 20160511 */
    pmeGPU->settings.performGPUSolve = true;
    /* FIXME: CPU gather with GPU spread has got to be broken as well due to different theta/dtheta layout. */
    pmeGPU->settings.performGPUGather = true;
    /* This is for the delayed atom data init hack and will get flipped back to false at first step */
    pmeGPU->settings.needToUpdateAtoms = true;

    pme_gpu_set_testing(pmeGPU, false);

    pme_gpu_init_specific(pmeGPU, hwinfo, gpu_opt);
    pme_gpu_init_sync_events(pmeGPU);
    pme_gpu_init_timings(pmeGPU);
    pme_gpu_alloc_energy_virial(pmeGPU);

    pme_gpu_copy_common_data_from(pme);

    GMX_ASSERT(pmeGPU->common->epsilon_r != 0.0f, "PME GPU: bad electrostatic coefficient");

    pme_gpu_kernel_params_base_t *kernelParamsPtr = pme_gpu_get_kernel_params_base_ptr(pmeGPU);
    kernelParamsPtr->constants.elFactor = ONE_4PI_EPS0 / pmeGPU->common->epsilon_r;
}

void pme_gpu_transform_spline_atom_data(const pme_gpu_t *pmeGPU, const pme_atomcomm_t *atc,
                                        PmeSplineDataType type, int dimIndex, PmeLayoutTransform transform)
{
    // The GPU atom spline data is laid out in a different way currently than the CPU one.
    // This function converts the data from GPU to CPU layout (in the host memory).
    // It is only intended for testing purposes so far.
    // Ideally we should use similar layouts on CPU and GPU if we care about mixed modes and their performance
    // (e.g. spreading on GPU, gathering on CPU).
    GMX_RELEASE_ASSERT(atc->nthread == 1, "Only the serial PME data layout is supported");
    const uintmax_t threadIndex  = 0;
    const auto      atomCount    = pme_gpu_get_kernel_params_base_ptr(pmeGPU)->atoms.nAtoms;
    const auto      atomsPerWarp = pme_gpu_get_atom_spline_data_alignment(pmeGPU);
    const auto      pmeOrder     = pmeGPU->common->pme_order;

    real           *cpuSplineBuffer;
    float          *h_splineBuffer;
    switch (type)
    {
        case PmeSplineDataType::Values:
            cpuSplineBuffer = atc->spline[threadIndex].theta[dimIndex];
            h_splineBuffer  = pmeGPU->staging.h_theta;
            break;

        case PmeSplineDataType::Derivatives:
            cpuSplineBuffer = atc->spline[threadIndex].dtheta[dimIndex];
            h_splineBuffer  = pmeGPU->staging.h_dtheta;
            break;

        default:
            GMX_THROW(gmx::InternalError("Unknown spline data type"));
    }

    for (auto atomIndex = 0; atomIndex < atomCount; atomIndex++)
    {
        auto atomWarpIndex = atomIndex % atomsPerWarp;
        auto warpIndex     = atomIndex / atomsPerWarp;
        for (auto orderIndex = 0; orderIndex < pmeOrder; orderIndex++)
        {
            const auto gpuValueIndex = ((pmeOrder * warpIndex + orderIndex) * DIM + dimIndex) * atomsPerWarp + atomWarpIndex;
            const auto cpuValueIndex = atomIndex * pmeOrder + orderIndex;
            GMX_ASSERT(cpuValueIndex < atomCount * pmeOrder, "Atom spline data index out of bounds (while transforming GPU data layout for host)");
            switch (transform)
            {
                case PmeLayoutTransform::GpuToHost:
                    cpuSplineBuffer[cpuValueIndex] = h_splineBuffer[gpuValueIndex];
                    break;

                case PmeLayoutTransform::HostToGpu:
                    h_splineBuffer[gpuValueIndex] = cpuSplineBuffer[cpuValueIndex];
                    break;

                default:
                    GMX_THROW(gmx::InternalError("Unknown layout transform"));
            }
        }
    }
}

void pme_gpu_reinit(gmx_pme_t *pme, const gmx_hw_info_t *hwinfo, const gmx_gpu_opt_t *gpu_opt)
{
    if (!pme_gpu_active(pme))
    {
        return;
    }

    if (!pme->gpu)
    {
        /* First-time initialization */
        pme_gpu_init(pme, hwinfo, gpu_opt);
    }
    else
    {
        /* After this call nothing in the GPU code should refer to the gmx_pme_t *pme itself - until the next pme_gpu_reinit */
        pme_gpu_copy_common_data_from(pme);
    }

    pme_gpu_reinit_grids(pme->gpu);
    pme_gpu_reinit_step(pme->gpu);
}

void pme_gpu_destroy(pme_gpu_t *pmeGPU)
{
    stopGpuProfiler();

    /* Free lots of data */
    pme_gpu_free_energy_virial(pmeGPU);
    pme_gpu_free_bspline_values(pmeGPU);
    pme_gpu_free_forces(pmeGPU);
    pme_gpu_free_coordinates(pmeGPU);
    pme_gpu_free_coefficients(pmeGPU);
    pme_gpu_free_spline_data(pmeGPU);
    pme_gpu_free_grid_indices(pmeGPU);
    pme_gpu_free_fract_shifts(pmeGPU);
    pme_gpu_free_grids(pmeGPU);

    pme_gpu_destroy_3dfft(pmeGPU);
    pme_gpu_destroy_sync_events(pmeGPU);
    pme_gpu_destroy_timings(pmeGPU);

    /* Free the GPU-framework specific data last */
    pme_gpu_destroy_specific(pmeGPU);

    delete pmeGPU;
}

void pme_gpu_reinit_atoms(pme_gpu_t *pmeGPU, const int nAtoms, const real *coefficients)
{
    pme_gpu_kernel_params_base_t *kernelParamsPtr = pme_gpu_get_kernel_params_base_ptr(pmeGPU);
    kernelParamsPtr->atoms.nAtoms = nAtoms;
    const int                     alignment = pme_gpu_get_atom_data_alignment(pmeGPU);
    pmeGPU->nAtomsPadded = ((nAtoms + alignment - 1) / alignment) * alignment;
    int                           nAtomsAlloc   = PME_GPU_USE_PADDING ? pmeGPU->nAtomsPadded : nAtoms;
    const bool                    haveToRealloc = (pmeGPU->nAtomsAlloc < nAtomsAlloc); /* This check might be redundant, but is logical */
    pmeGPU->nAtomsAlloc = nAtomsAlloc;

#if GMX_DOUBLE
    GMX_RELEASE_ASSERT(false, "Only single precision supported");
    GMX_UNUSED_VALUE(coefficients);
#else
    pme_gpu_realloc_and_copy_input_coefficients(pmeGPU, reinterpret_cast<const float *>(coefficients));
    /* Could also be checked for haveToRealloc, but the copy always needs to be performed */
#endif

    if (haveToRealloc)
    {
        pme_gpu_realloc_coordinates(pmeGPU);
        pme_gpu_realloc_forces(pmeGPU);
        pme_gpu_realloc_spline_data(pmeGPU);
        pme_gpu_realloc_grid_indices(pmeGPU);
    }
}
