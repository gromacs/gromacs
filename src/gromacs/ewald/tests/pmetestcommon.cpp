/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * \brief
 * Implements common routines for PME tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "pmetestcommon.h"

#include "config.h"

#include <gtest/gtest.h>

#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme-gpu-internal.h"
#include "gromacs/ewald/pme-grid.h"
#include "gromacs/ewald/pme-internal.h"
#include "gromacs/ewald/pme-spread.h"
#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/loggerbuilder.h"

namespace gmx
{
//! This constructs the test environment
PmeTestEnvironment * const pmeEnv = (PmeTestEnvironment * const)::testing::AddGlobalTestEnvironment(new PmeTestEnvironment);

//! Simple hardware initialization
void PmeTestEnvironment::hardwareInit()
{
    commrec_.reset(init_commrec());
    gmx_init_intranode_counters(commrec_.get());
    //gmx_hw_opt_t *hardwareOptions = new gmx_hw_opt_t();
    gpuOptions_.reset(new gmx_gpu_opt_t {}); // &hardwareOptions->gpu_opt;
    LoggerBuilder builder;
    LoggerOwner   logOwner(builder.build());
    MDLogger      log(logOwner.logger());
    hardwareInfo_.reset(gmx_detect_hardware(log, commrec_.get(), true));
    pick_compatible_gpus(&hardwareInfo_.get()->gpu_info, gpuOptions_.get());
    gmx_select_gpu_ids(log, commrec_.get(), &hardwareInfo_.get()->gpu_info, false, gpuOptions_.get());
    // TODO: iterate through all the GPUs for testing?
}

void PmeTestEnvironment::SetUp()
{
#if GMX_GPU != GMX_GPU_NONE // TODO - remove conditional, when gmx_select_gpu_ids() stops failing in CPU-only reference build
    hardwareInit();
#endif
}

bool PmeSupportsInputForMode(const t_inputrec *inputRec, PmeCodePath mode)
{
    const bool notImplemented = (mode == PmeCodePath::CUDA) && (inputRec->pme_order != 4);
    // This conditional copies some of the logic from pme_gpu_check_restrictions().
    // Duplicating the constraints here emphasizes them (changing them requires changing the test as well).
    return !notImplemented;
}

//! PME initialization - internal
PmeSafePointer PmeInitInternal(const t_inputrec *inputRec, PmeCodePath mode, size_t atomCount)
{
    gmx_pme_t           *pmeDataRaw   = NULL;
    const bool           useGpu       = (mode == PmeCodePath::CPU) ? false : true;

    const gmx_hw_info_t *hardwareInfo = pmeEnv->getHardwareInfo();
    const gmx_gpu_opt_t *gpuOptions   = pmeEnv->getGpuOptions();

    const real           ewaldCoeffQ = 1.0;
    gmx_pme_init(&pmeDataRaw, NULL, 1, 1, inputRec, atomCount, FALSE, FALSE, TRUE,
                 ewaldCoeffQ, 0.0, 1, useGpu, nullptr, hardwareInfo, gpuOptions);
    PmeSafePointer pme(pmeDataRaw, gmx_pme_destroy); // taking ownership
    return pme;
}

//! Simple PME initialization based on input, no atom data
PmeSafePointer PmeInitEmpty(const t_inputrec *inputRec)
{
    return PmeInitInternal(inputRec, PmeCodePath::CPU, 0);
    // hiding the fact that PME actually needs to know the number of atoms in advance
}

//! PME initialization with atom data and system box
PmeSafePointer PmeInitWithAtoms(const t_inputrec        *inputRec,
                                PmeCodePath              mode,
                                const CoordinatesVector &coordinates,
                                const ChargesVector     &charges,
                                const Matrix3x3         &box
                                )
{
    const size_t    atomCount = coordinates.size();
    GMX_RELEASE_ASSERT(atomCount == charges.size(), "Mismatch in atom data");
    PmeSafePointer  pmeSafe = PmeInitInternal(inputRec, mode, atomCount);
    pme_atomcomm_t *atc     = nullptr;

    // TODO get rid of this with proper matrix type
    matrix boxTemp;
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            boxTemp[i][j] = box[i * DIM + j];
        }
    }

    switch (mode)
    {
        case PmeCodePath::CPU:
            gmx::invertBoxMatrix(boxTemp, pmeSafe->recipbox);
            atc              = &(pmeSafe->atc[0]);
            atc->x           = const_cast<rvec *>(as_rvec_array(coordinates.data()));
            atc->coefficient = const_cast<real *>(charges.data());
            /* With decomposition there would be more boilerplate atc code here, e.g. do_redist_pos_coeffs */
            break;

        case PmeCodePath::CUDA:
            pme_gpu_set_testing(pmeSafe->gpu, true);
            gmx_pme_reinit_atoms(pmeSafe.get(), atomCount, charges.data());
            pme_gpu_start_step(pmeSafe->gpu, boxTemp, as_rvec_array(coordinates.data()));
            break;

        default:
            GMX_THROW(gmx::InternalError("Test not implemented for this mode"));
    }

    return pmeSafe;
}

//! PME spline calculation and charge spreading
void PmePerformSplineAndSpread(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                               bool computeSplines, bool spreadCharges)
{
    gmx_pme_t      *pmeUnsafe                    = pmeSafe.get();
    pme_atomcomm_t *atc                          = &(pmeSafe->atc[0]);
    const size_t    gridIndex                    = 0;
    const bool      computeSplinesForZeroCharges = true;
    real           *fftgrid                      = spreadCharges ? pmeSafe->fftgrid[gridIndex] : nullptr;
    real           *h_grid                       = pmeSafe->pmegrid[gridIndex].grid.grid;
    switch (mode)
    {
        case PmeCodePath::CPU:
            spread_on_grid(pmeUnsafe, atc, &pmeSafe->pmegrid[gridIndex], computeSplines, spreadCharges,
                           fftgrid, computeSplinesForZeroCharges, gridIndex);
            if (spreadCharges && !pmeSafe->bUseThreads)
            {
                wrap_periodic_pmegrid(pmeUnsafe, h_grid);
                copy_pmegrid_to_fftgrid(pmeUnsafe, h_grid, fftgrid, gridIndex);
            }
            break;

        case PmeCodePath::CUDA:
            pme_gpu_spread(pmeSafe->gpu, gridIndex, h_grid, computeSplines, spreadCharges);
            break;

        default:
            GMX_THROW(gmx::InternalError("Test not implemented for this mode"));
    }
}

//! Fetching the spline computation outputs of PmePerformSplineAndSpread()
void PmeFetchOutputsSpline(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                           SplineParamsVector &splineValues,
                           SplineParamsVector &splineDerivatives,
                           GridLineIndicesVector &gridLineIndices)

{
    const pme_atomcomm_t *atc         = &(pmeSafe->atc[0]);
    const size_t          atomCount   = atc->n;
    const size_t          pmeOrder    = pmeSafe->pme_order;
    const size_t          threadIndex = 0; // relying on running single threaded on CPU
    switch (mode)
    {
        case PmeCodePath::CUDA:
            pme_gpu_synchronize(pmeSafe->gpu);
            pme_gpu_sync_spline_atom_data(pmeSafe->gpu);
            pme_gpu_transform_spline_atom_data_for_host(pmeSafe->gpu, atc);
        // intentional absence of break - the atom data has been copied into PME CPU buffers with proper layout

        case PmeCodePath::CPU:
            gridLineIndices.assign(atc->idx, atc->idx + atomCount);
            // spline values - XX...XXYY...YYZZ...ZZ
            splineValues.clear();
            splineDerivatives.clear();
            for (int i = 0; i < DIM; i++)
            {
                splineValues.insert(splineValues.end(), atc->spline[threadIndex].theta[i],
                                    atc->spline[threadIndex].theta[i] + atomCount * pmeOrder);
                splineDerivatives.insert(splineDerivatives.end(), atc->spline[threadIndex].dtheta[i],
                                         atc->spline[threadIndex].dtheta[i] + atomCount * pmeOrder);
            }
            break;

        default:
            GMX_THROW(gmx::InternalError("Test not implemented for this mode"));
    }
}

//! Fetching the spreading output of PmePerformSplineAndSpread()
void PmeFetchOutputsSpread(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                           SparseGridValues &gridValues)
{
    const size_t          gridIndex = 0;
    IVec                  gridSize, gridOffsetUnused, paddedGridSize;

    switch (mode)
    {
        case PmeCodePath::CUDA:
            pme_gpu_synchronize(pmeSafe->gpu);
            pme_gpu_sync_spread_grid(pmeSafe->gpu);
        // intentional absence of break - the grid has been copied into pmeSafe->fftgrid, corresponding to the CPU code

        case PmeCodePath::CPU:
            gridValues.clear();
            gmx_parallel_3dfft_real_limits(pmeSafe->pfft_setup[gridIndex], gridSize, gridOffsetUnused, paddedGridSize);
            for (int ix = 0; ix < gridSize[XX]; ix++)
            {
                for (int iy = 0; iy < gridSize[YY]; iy++)
                {
                    for (int iz = 0; iz < gridSize[ZZ]; iz++)
                    {
                        const size_t gridValueIndex = (ix * paddedGridSize[YY] + iy) * paddedGridSize[ZZ] + iz;
                        const real   value          = pmeSafe->fftgrid[gridIndex][gridValueIndex];
                        if (value != 0.0)
                        {
                            IVec key = {ix, iy, iz};
                            gridValues[key] = value;
                        }
                    }
                }
            }
            break;

        default:
            GMX_THROW(gmx::InternalError("Test not implemented for this mode"));
    }
}

}
