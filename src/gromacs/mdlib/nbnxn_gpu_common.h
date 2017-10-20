/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * \brief Implements common routines for different NBNXN GPU implementations
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_mdlib
 */

#ifndef GMX_MDLIB_NBNXN_GPU_COMMON_H
#define GMX_MDLIB_NBNXN_GPU_COMMON_H

#include "config.h"

#include <string>

#if GMX_GPU == GMX_GPU_CUDA
#include "nbnxn_cuda/nbnxn_cuda_types.h"
#endif

#if GMX_GPU == GMX_GPU_OPENCL
#include "nbnxn_ocl/nbnxn_ocl_types.h"
#endif
#include "gromacs/utility/stringutil.h"


/*! \brief An early return condition for empty NB GPU workloads
 *
 * This is currently used for non-local kernels/transfers only.
 * Skipping the local kernel is more complicated, since the
 * local part of the force array also depends on the non-local kernel.
 * The skip of the local kernel is taken care of separately.
 */
static inline bool canSkipWork(const gmx_nbnxn_gpu_t *nb, int iloc)
{
    assert(nb && nb->plist[iloc]);
    return (iloc == eintNonlocal) && (nb->plist[iloc]->nsci == 0);
}

/*! \brief Check that atom locality values are valid for the GPU module.
 *
 *  In the GPU module atom locality "all" is not supported, the local and
 *  non-local ranges are treated separately.
 *
 *  \param[in] atomLocality atom locality specifier
 */
static inline void validateGpuAtomLocality(int atomLocality)
{
    std::string str = gmx::formatString("Invalid atom locality passed (%d); valid here is only "
                                        "local (%d) or nonlocal (%d)", atomLocality, eatLocal, eatNonlocal);

    GMX_ASSERT(LOCAL_OR_NONLOCAL_A(atomLocality), str.c_str());
}

/*! \brief Convert atom locality to interaction locality.
 *
 *  In the current implementation the this is straightforward conversion:
 *  local to local, non-local to non-local.
 *
 *  \param[in] atomLocality Atom locality specifier
 *  \returns                Interaction locality corresponding to the atom locality passed.
 */
static inline int gpuAtomToInteractionLocality(int atomLocality)
{
    validateGpuAtomLocality(atomLocality);

    /* determine interaction locality from atom locality */
    if (LOCAL_A(atomLocality))
    {
        return eintLocal;
    }
    else if (NONLOCAL_A(atomLocality))
    {
        return eintNonlocal;
    }
    else
    {
        // can't be reached
        assert(false);
        return -1;
    }
}

/*! \brief Calculate atom range and return start index and length.
 *
 * \param[in] atomData Atom descriptor data structure
 * \param[in] atomLocality Atom locality specifier
 * \param[out] atomRangeBegin Starting index of the atom range in the atom data array.
 * \param[out] atomRangeLen Atom range length in the atom data array.
 */
template <typename AtomDataT>
static inline void getGpuAtomRange(const AtomDataT *atomData,
                                   int              atomLocality,
                                   int             &atomRangeBegin,
                                   int             &atomRangeLen)
{
    assert(atomData);
    validateGpuAtomLocality(atomLocality);

    /* calculate the atom data index range based on locality */
    if (LOCAL_A(atomLocality))
    {
        atomRangeBegin  = 0;
        atomRangeLen    = atomData->natoms_local;
    }
    else
    {
        atomRangeBegin  = atomData->natoms_local;
        atomRangeLen    = atomData->natoms - atomData->natoms_local;
    }
}

#endif
