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

#include "gmxpre.h"

#include "MimicUtils.h"

std::vector<int> genQmmmIndices(const gmx_mtop_t &mtop)
{
    std::vector<int> output;
    int              global_at = 0;
    unsigned char   *grpnr     = mtop.groups.grpnr[egcQMMM];
    for (const gmx_molblock_t &molb : mtop.molblock)
    {
        for (int mol = 0; mol < molb.nmol; ++mol)
        {
            for (int n_atom = 0; n_atom < mtop.moltype[molb.type].atoms.nr; ++n_atom)
            {
                if (!grpnr || !grpnr[global_at])
                {
                    output.push_back(global_at);
                }
                ++global_at;
            }
        }
    }
    return output;
}

void addExclusions(t_blocka               *excls,
                   const std::vector<int> &ids)
{
    t_blocka inter_excl {};
    init_blocka(&inter_excl);
    size_t   n_q = ids.size();

    inter_excl.nr  = excls->nr;
    inter_excl.nra = n_q * n_q;

    size_t total_nra = n_q * n_q;

    snew(inter_excl.index, excls->nr + 1);
    snew(inter_excl.a, total_nra);

    for (int i = 0; i < excls->nr; ++i)
    {
        inter_excl.index[i] = 0;
    }

    int prev_index = 0;
    for (int k = 0; k < inter_excl.nr; ++k)
    {
        inter_excl.index[k] = prev_index;
        for (size_t i = 0; i < n_q; ++i)
        {
            if (k != ids[i])
            {
                continue;
            }
            size_t index = n_q * i;
            inter_excl.index[ids[i]]     = index;
            prev_index                   = index + n_q;
            for (size_t j = 0; j < n_q; ++j)
            {
                inter_excl.a[n_q * i + j] = ids[j];
            }
        }
    }
    inter_excl.index[ids[n_q - 1] + 1] = n_q * n_q;

    inter_excl.index[inter_excl.nr] = n_q * n_q;

    t_block2  qmexcl2 {};
    init_block2(&qmexcl2, excls->nr);
    b_to_b2(&inter_excl, &qmexcl2);
    merge_excl(excls, &qmexcl2, nullptr);
    done_block2(&qmexcl2);
}
