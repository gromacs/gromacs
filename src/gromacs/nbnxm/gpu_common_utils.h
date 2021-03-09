/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2019,2020,2021, by the GROMACS development team, led by
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
 * \brief Implements common util routines for different NBNXN GPU implementations
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_GPU_COMMON_UTILS_H
#define GMX_NBNXM_GPU_COMMON_UTILS_H

#include "config.h"

#include "gromacs/listed_forces/gpubonded.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/range.h"
#include "gromacs/nbnxm/nbnxm_gpu.h"

#if GMX_GPU_CUDA
#    include "cuda/nbnxm_cuda_types.h"
#endif

#if GMX_GPU_OPENCL
#    include "opencl/nbnxm_ocl_types.h"
#endif

namespace Nbnxm
{

/*! \brief An early return condition for empty NB GPU workloads
 *
 * This is currently used for non-local kernels/transfers only.
 * Skipping the local kernel is more complicated, since the
 * local part of the force array also depends on the non-local kernel.
 * The skip of the local kernel is taken care of separately.
 */
static inline bool canSkipNonbondedWork(const NbnxmGpu& nb, InteractionLocality iloc)
{
    assert(nb.plist[iloc]);
    return (iloc == InteractionLocality::NonLocal && nb.plist[iloc]->nsci == 0);
}

/*! \brief Check that atom locality values are valid for the GPU module.
 *
 *  In the GPU module atom locality "all" is not supported, the local and
 *  non-local ranges are treated separately.
 *
 *  \param[in] atomLocality atom locality specifier
 */
static inline void validateGpuAtomLocality(const AtomLocality atomLocality)
{
    std::string str = gmx::formatString(
            "Invalid atom locality passed (%d); valid here is only "
            "local (%d) or nonlocal (%d)",
            static_cast<int>(atomLocality),
            static_cast<int>(AtomLocality::Local),
            static_cast<int>(AtomLocality::NonLocal));

    GMX_ASSERT(atomLocality == AtomLocality::Local || atomLocality == AtomLocality::NonLocal, str.c_str());
}

/*! \brief Convert atom locality to interaction locality.
 *
 *  In the current implementation the this is straightforward conversion:
 *  local to local, non-local to non-local.
 *
 *  \param[in] atomLocality Atom locality specifier
 *  \returns                Interaction locality corresponding to the atom locality passed.
 */
static inline InteractionLocality gpuAtomToInteractionLocality(const AtomLocality atomLocality)
{
    validateGpuAtomLocality(atomLocality);

    /* determine interaction locality from atom locality */
    if (atomLocality == AtomLocality::Local)
    {
        return InteractionLocality::Local;
    }
    else if (atomLocality == AtomLocality::NonLocal)
    {
        return InteractionLocality::NonLocal;
    }
    else
    {
        gmx_incons("Wrong locality");
    }
}

/*! \brief Returns true if there is GPU short-range work for the given interaction locality.
 *
 * Note that as, unlike nonbonded tasks, bonded tasks are not split into local/nonlocal,
 * and therefore if there are GPU offloaded bonded interactions, this function will return
 * true for all interaction localities.
 *
 * \param[inout]  nb        Pointer to the nonbonded GPU data structure
 * \param[in]     iLocality Interaction locality identifier
 */
static inline bool haveGpuShortRangeWork(const NbnxmGpu& nb, const gmx::InteractionLocality iLocality)
{
    return nb.haveWork[iLocality];
}

/*! \brief Calculate atom range and return start index and length.
 *
 * \param[in] atomData Atom descriptor data structure
 * \param[in] atomLocality Atom locality specifier
 * \returns Range of indexes for selected locality.
 */
static inline gmx::Range<int> getGpuAtomRange(const NBAtomData* atomData, const AtomLocality atomLocality)
{
    assert(atomData);
    validateGpuAtomLocality(atomLocality);

    /* calculate the atom data index range based on locality */
    if (atomLocality == AtomLocality::Local)
    {
        return gmx::Range<int>(0, atomData->numAtomsLocal);
    }
    else
    {
        return gmx::Range<int>(atomData->numAtomsLocal, atomData->numAtoms);
    }
}

} // namespace Nbnxm

#endif
