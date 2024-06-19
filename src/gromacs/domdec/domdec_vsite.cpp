/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2006- The GROMACS Authors
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
 * \brief This file implements functions for domdec to use
 * while managing inter-atomic constraints.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "domdec_vsite.h"

#include <cassert>
#include <cstdio>

#include <algorithm>
#include <memory>
#include <vector>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/domdec/hashedmap.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#include "domdec_specatomcomm.h"

void dd_move_f_vsites(const gmx_domdec_t& dd, gmx::ArrayRef<gmx::RVec> f, gmx::ArrayRef<gmx::RVec> fshift)
{
    if (dd.vsite_comm)
    {
        dd_move_f_specat(&dd, dd.vsite_comm.get(), f.data(), fshift.data());
    }
}

void dd_clear_f_vsites(const gmx_domdec_t& dd, gmx::ArrayRef<gmx::RVec> f)
{
    if (dd.vsite_comm)
    {
        for (int i = dd.vsite_comm->at_start; i < dd.vsite_comm->at_end; i++)
        {
            clear_rvec(f[i]);
        }
    }
}

void dd_move_x_vsites(const gmx_domdec_t& dd, const matrix box, gmx::ArrayRef<gmx::RVec> x)
{
    if (dd.vsite_comm)
    {
        dd_move_x_specat(&dd, dd.vsite_comm.get(), box, x.data(), nullptr, FALSE);
    }
}

void dd_move_x_and_v_vsites(const gmx_domdec_t&      dd,
                            const matrix             box,
                            gmx::ArrayRef<gmx::RVec> x,
                            gmx::ArrayRef<gmx::RVec> v)
{
    if (dd.vsite_comm)
    {
        dd_move_x_specat(&dd, dd.vsite_comm.get(), box, x.data(), v.data(), FALSE);
    }
}

void dd_clear_local_vsite_indices(gmx_domdec_t* dd)
{
    if (dd->vsite_comm)
    {
        dd->ga2la_vsite->clearAndResizeHashTable();
    }
}

int dd_make_local_vsites(gmx_domdec_t* dd, int at_start, gmx::ArrayRef<InteractionList> lil)
{
    std::vector<int>&    ireq         = dd->vsite_requestedGlobalAtomIndices;
    gmx::HashedMap<int>* ga2la_specat = dd->ga2la_vsite.get();

    ireq.clear();
    /* Loop over all the home vsites */
    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            const int              nral = NRAL(ftype);
            const InteractionList& lilf = lil[ftype];
            for (int i = 0; i < lilf.size(); i += 1 + nral)
            {
                const int* iatoms = lilf.iatoms.data() + i;
                /* Check if we have the other atoms */
                for (int j = 1; j < 1 + nral; j++)
                {
                    if (iatoms[j] < 0)
                    {
                        /* This is not a home atom,
                         * we need to ask our neighbors.
                         */
                        int a = -iatoms[j] - 1;
                        /* Check to not ask for the same atom more than once */
                        if (!dd->ga2la_vsite->find(a))
                        {
                            /* Add this non-home atom to the list */
                            ireq.push_back(a);
                            /* Temporarily mark with -2,
                             * we get the index later.
                             */
                            ga2la_specat->insert(a, -2);
                        }
                    }
                }
            }
        }
    }

    int at_end = setup_specat_communication(
            dd, &ireq, dd->vsite_comm.get(), ga2la_specat, at_start, 2, "vsite", "");

    /* Fill in the missing indices */
    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            const int        nral = NRAL(ftype);
            InteractionList& lilf = lil[ftype];
            for (int i = 0; i < lilf.size(); i += 1 + nral)
            {
                t_iatom* iatoms = lilf.iatoms.data() + i;
                for (int j = 1; j < 1 + nral; j++)
                {
                    if (iatoms[j] < 0)
                    {
                        const int* a = ga2la_specat->find(-iatoms[j] - 1);
                        GMX_ASSERT(a, "We have checked before that this atom index has been set");
                        iatoms[j] = *a;
                    }
                }
            }
        }
    }

    return at_end;
}

void init_domdec_vsites(gmx_domdec_t* dd, int n_intercg_vsite)
{
    if (debug)
    {
        fprintf(debug, "Begin init_domdec_vsites\n");
    }

    /* Use a hash table for the global to local index.
     * The number of keys is a rough estimate, it will be optimized later.
     */
    int numKeysEstimate = std::min(n_intercg_vsite / 20, n_intercg_vsite / (2 * dd->nnodes));
    dd->ga2la_vsite     = std::make_unique<gmx::HashedMap<int>>(numKeysEstimate);

    dd->vsite_comm = std::make_unique<gmx_domdec_specat_comm_t>();
}
