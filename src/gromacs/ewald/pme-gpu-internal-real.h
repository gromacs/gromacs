/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
#ifndef PMEGPUINTERNALREAL_H
#define PMEGPUINTERNALREAL_H

// everything should be fully GPU-agnostic here

//! Type of spline data
enum class PmeSplineDataType
{
    Values,      // theta
    Derivatives, // dtheta
};               //TODO move this into new and shiny pme.h (pme-types.h?)

//! PME grid dimension ordering (from major to minor)
enum class GridOrdering
{
    YZX,
    XYZ
};

//! A binary enum for spline data layout transformation
enum class PmeLayoutTransform
{
    GpuToHost,
    HostToGpu
};

//FIXME this is a bunch of duplicate declarations and includes - just to make tests compile.
// Please nuke this and the tests' logic

#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/ewald-utils.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/gpu_macros.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/units.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"

struct PmeGpu;
struct pme_atomcomm_t;
void pme_gpu_set_testing(PmeGpu *pmeGpu, bool testing);
void pme_gpu_update_input_box(PmeGpu      *pmeGpu,
                              const matrix box);
void pme_gpu_get_real_grid_sizes(const PmeGpu *pmeGpu, gmx::IVec *gridSize, gmx::IVec *paddedGridSize);
GPU_FUNC_QUALIFIER void pme_gpu_synchronize(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM
void pme_gpu_transform_spline_atom_data(const PmeGpu *pmeGpu, const pme_atomcomm_t *atc,
                                        PmeSplineDataType type, int dimIndex, PmeLayoutTransform transform);
GPU_FUNC_QUALIFIER void pme_gpu_copy_input_coordinates(const PmeGpu    *GPU_FUNC_ARGUMENT(pmeGpu),
                                                       const rvec      *GPU_FUNC_ARGUMENT(h_coordinates)) GPU_FUNC_TERM
GPU_FUNC_QUALIFIER void pme_gpu_spread(PmeGpu    *GPU_FUNC_ARGUMENT(pmeGpu),
                                       int              GPU_FUNC_ARGUMENT(gridIndex),
                                       real            *GPU_FUNC_ARGUMENT(h_grid),
                                       bool             GPU_FUNC_ARGUMENT(computeSplines),
                                       bool             GPU_FUNC_ARGUMENT(spreadCharges)) GPU_FUNC_TERM
void pme_gpu_get_energy_virial(const PmeGpu *pmeGpu, real *energy, matrix virial);

/*! \libinternal \brief
 * A GPU Fourier space solving function.
 *
 * \param[in]     pmeGpu                  The PME GPU structure.
 * \param[in,out] h_grid                  The host-side input and output Fourier grid buffer (used only with testing or host-side FFT)
 * \param[in]     gridOrdering            Specifies the dimenion ordering of the complex grid. TODO: store this information?
 * \param[in]     computeEnergyAndVirial  Tells if the energy and virial computation should also be performed.
 */
GPU_FUNC_QUALIFIER void pme_gpu_solve(PmeGpu    *GPU_FUNC_ARGUMENT(pmeGpu),
                                      t_complex       *GPU_FUNC_ARGUMENT(h_grid),
                                      GridOrdering     GPU_FUNC_ARGUMENT(gridOrdering),
                                      bool             GPU_FUNC_ARGUMENT(computeEnergyAndVirial)) GPU_FUNC_TERM


/*! \libinternal \brief
 * A GPU force gathering function.
 *
 * \param[in]     pmeGpu           The PME GPU structure.
 * \param[in]     forceTreatment   Tells how data in h_forces should be treated.
 *                                 TODO: determine efficiency/balance of host/device-side reductions.
 * \param[in]     h_grid           The host-side grid buffer (used only in testing mode)
 */
GPU_FUNC_QUALIFIER void pme_gpu_gather(PmeGpu                *GPU_FUNC_ARGUMENT(pmeGpu),
                                       PmeForceOutputHandling GPU_FUNC_ARGUMENT(forceTreatment),
                                       const float           *GPU_FUNC_ARGUMENT(h_grid)
                                       ) GPU_FUNC_TERM


gmx::ArrayRef<gmx::RVec> pme_gpu_get_forces(PmeGpu *pmeGpu);

gmx::ArrayRef<const gmx::IVec> pmeGpuGetGridlineIndices(const PmeGpu *pmeGpu);

void pmeGpuSetGridlineIndices(PmeGpu *pmeGpu, gmx::ArrayRef<const gmx::IVec> gridlineIndices);


GPU_FUNC_QUALIFIER PmePersistentDataHandle pmeGpuAcquirePersistentData(PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM_WITH_RETURN(nullptr)

#endif // PMEGPUINTERNALREAL_H
