/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/*! \internal \file
 *
 * \brief
 * Implements utility functions used by all nbnxm CPU kernels.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#include "gmxpre.h"

#include "kernel_common.h"

#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/gmxassert.h"

//! Clears all elements of buffer
static void clearBufferAll(gmx::ArrayRef<real> buffer)
{
    for (real& elem : buffer)
    {
        elem = 0;
    }
}

/*! \brief Clears elements of size and stride \p numComponentsPerElement
 *
 * Only elements with flags in \p nbat set for index \p outputIndex
 * are cleared.
 */
template<int numComponentsPerElement>
static void clearBufferFlagged(const nbnxn_atomdata_t& nbat, int outputIndex, gmx::ArrayRef<real> buffer)
{
    gmx::ArrayRef<const gmx_bitmask_t> flags = nbat.buffer_flags;
    gmx_bitmask_t                      our_flag; // NOLINT(cppcoreguidelines-init-variables)
    bitmask_init_bit(&our_flag, outputIndex);

    constexpr size_t numComponentsPerBlock = NBNXN_BUFFERFLAG_SIZE * numComponentsPerElement;

    for (size_t b = 0; b < flags.size(); b++)
    {
        if (!bitmask_is_disjoint(flags[b], our_flag))
        {
            clearBufferAll(buffer.subArray(b * numComponentsPerBlock, numComponentsPerBlock));
        }
    }
}

void clearForceBuffer(nbnxn_atomdata_t* nbat, int outputIndex)
{
    if (nbat->bUseBufferFlags)
    {
        GMX_ASSERT(nbat->fstride == DIM, "Only fstride=3 is currently handled here");

        clearBufferFlagged<DIM>(*nbat, outputIndex, nbat->out[outputIndex].f);
    }
    else
    {
        clearBufferAll(nbat->out[outputIndex].f);
    }
}

void clear_fshift(real* fshift)
{
    for (int i = 0; i < gmx::c_numShiftVectors * DIM; i++)
    {
        fshift[i] = 0;
    }
}

void reduce_energies_over_lists(const nbnxn_atomdata_t* nbat, int nlist, real* Vvdw, real* Vc)
{
    const int nenergrp = nbat->params().nenergrp;

    for (int nb = 0; nb < nlist; nb++)
    {
        for (int i = 0; i < nenergrp; i++)
        {
            /* Reduce the diagonal terms */
            int ind = i * nenergrp + i;
            Vvdw[ind] += nbat->out[nb].Vvdw[ind];
            Vc[ind] += nbat->out[nb].Vc[ind];

            /* Reduce the off-diagonal terms */
            for (int j = i + 1; j < nenergrp; j++)
            {
                /* The output should contain only one off-diagonal part */
                int ind  = i * nenergrp + j;
                int indr = j * nenergrp + i;
                Vvdw[ind] += nbat->out[nb].Vvdw[ind] + nbat->out[nb].Vvdw[indr];
                Vc[ind] += nbat->out[nb].Vc[ind] + nbat->out[nb].Vc[indr];
            }
        }
    }
}
