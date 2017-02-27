/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "read-conformation.h"

#include <vector>

#include "gromacs/fileio/confio.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/unique_cptr.h"

using gmx::RVec;

std::vector<real>
makeExclusionDistances(const t_atoms *a, gmx_atomprop_t aps,
                       real defaultDistance, real scaleFactor)
{
    std::vector<real> exclusionDistances;

    exclusionDistances.reserve(a->nr);
    for (int i = 0; i < a->nr; ++i)
    {
        real value;
        if (!gmx_atomprop_query(aps, epropVDW,
                                *(a->resinfo[a->atom[i].resind].name),
                                *(a->atomname[i]), &value))
        {
            value = defaultDistance;
        }
        else
        {
            value *= scaleFactor;
        }
        exclusionDistances.push_back(value);
    }
    return exclusionDistances;
}

void readConformation(const char *confin, gmx_mtop_t *top,
                      std::vector<RVec> *x, std::vector<RVec> *v,
                      int *ePBC, matrix box, const char *statusTitle)
{
    fprintf(stderr, "Reading %s configuration%s\n", statusTitle,
            v ? " and velocities" : "");
    rvec                   *x_tmp = nullptr, *v_tmp = nullptr;
    bool                    dummy;
    readConfAndTopology(confin, &dummy, top, ePBC, x ? &x_tmp : nullptr, v ? &v_tmp : nullptr, box);
    const gmx::sfree_guard  xguard(x_tmp);
    const gmx::sfree_guard  vguard(v_tmp);
    if (x && x_tmp)
    {
        *x = std::vector<RVec>(x_tmp, x_tmp + top->natoms);
    }
    if (v && v_tmp)
    {
        *v = std::vector<RVec>(v_tmp, v_tmp + top->natoms);
    }
    fprintf(stderr, "%s\nContaining %d atoms in %d residues\n",
            *top->name, top->natoms, gmx_mtop_nres(top));
}

void readConformation(const char *confin, t_topology *top,
                      std::vector<RVec> *x, std::vector<RVec> *v,
                      int *ePBC, matrix box, const char *statusTitle)
{
    fprintf(stderr, "Reading %s configuration%s\n", statusTitle,
            v ? " and velocities" : "");
    rvec                   *x_tmp = nullptr, *v_tmp = nullptr;
    read_tps_conf(confin, top, ePBC, x ? &x_tmp : nullptr, v ? &v_tmp : nullptr, box, FALSE);
    const gmx::sfree_guard  xguard(x_tmp);
    const gmx::sfree_guard  vguard(v_tmp);
    if (x && x_tmp)
    {
        *x = std::vector<RVec>(x_tmp, x_tmp + top->atoms.nr);
    }
    if (v && v_tmp)
    {
        *v = std::vector<RVec>(v_tmp, v_tmp + top->atoms.nr);
    }
    fprintf(stderr, "%s\nContaining %d atoms in %d residues\n",
            *top->name, top->atoms.nr, top->atoms.nres);
}
