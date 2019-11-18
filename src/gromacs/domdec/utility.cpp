/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * \brief Declares utility functions used in the domain decomposition module.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "utility.h"

#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "domdec_internal.h"

char dim2char(int dim)
{
    char c = '?';

    switch (dim)
    {
        case XX: c = 'X'; break;
        case YY: c = 'Y'; break;
        case ZZ: c = 'Z'; break;
        default: gmx_fatal(FARGS, "Unknown dim %d", dim);
    }

    return c;
}

void make_tric_corr_matrix(int npbcdim, const matrix box, matrix tcm)
{
    if (YY < npbcdim)
    {
        tcm[YY][XX] = -box[YY][XX] / box[YY][YY];
    }
    else
    {
        tcm[YY][XX] = 0;
    }
    if (ZZ < npbcdim)
    {
        tcm[ZZ][XX] = -(box[ZZ][YY] * tcm[YY][XX] + box[ZZ][XX]) / box[ZZ][ZZ];
        tcm[ZZ][YY] = -box[ZZ][YY] / box[ZZ][ZZ];
    }
    else
    {
        tcm[ZZ][XX] = 0;
        tcm[ZZ][YY] = 0;
    }
}

void check_screw_box(const matrix box)
{
    /* Mathematical limitation */
    if (box[YY][XX] != 0 || box[ZZ][XX] != 0)
    {
        gmx_fatal(FARGS,
                  "With screw pbc the unit cell can not have non-zero off-diagonal x-components");
    }

    /* Limitation due to the asymmetry of the eighth shell method */
    if (box[ZZ][YY] != 0)
    {
        gmx_fatal(FARGS, "pbc=screw with non-zero box_zy is not supported");
    }
}
/*! \brief Resize the state and f*/
void dd_resize_state(t_state* state, PaddedHostVector<gmx::RVec>* f, int natoms)
{
    if (debug)
    {
        fprintf(debug, "Resizing state: currently %d, required %d\n", state->natoms, natoms);
    }

    state_change_natoms(state, natoms);

    if (f != nullptr)
    {
        /* We need to allocate one element extra, since we might use
         * (unaligned) 4-wide SIMD loads to access rvec entries.
         */
        f->resizeWithPadding(natoms);
    }
}

/*! \brief Ensure fr, state and f, if != nullptr, can hold numChargeGroups
 *         atoms for the Verlet scheme and charge groups for the group scheme.
 *
 * todo refactor this now that group scheme is removed
 */
void dd_check_alloc_ncg(t_forcerec* fr, t_state* state, PaddedHostVector<gmx::RVec>* f, int numChargeGroups)
{
    fr->cginfo.resize(numChargeGroups);

    /* We use x during the setup of the atom communication */
    dd_resize_state(state, f, numChargeGroups);
}
