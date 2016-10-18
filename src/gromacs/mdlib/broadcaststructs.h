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
#ifndef GMX_MDLIB_BROADCASTSTRUCTS_H
#define GMX_MDLIB_BROADCASTSTRUCTS_H

#include <cmath>

#include <vector>

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/smalloc.h"

template <typename T>
void block_bc(const t_commrec *cr, T &data)
{
    gmx_bcast(sizeof(T), static_cast<void *>(&data), cr);
}
template <typename T>
void nblock_bc(const t_commrec *cr, int numElements, T *data)
{
    gmx_bcast(numElements * sizeof(T), static_cast<void *>(data), cr);
}
template <typename T>
void snew_bc(const t_commrec *cr, T * &data, int numElements)
{
    if (!MASTER(cr))
    {
        snew(data, numElements);
    }
}
template <typename T>
static void nblock_abc(const t_commrec *cr, int numElements, T **v)
{
    snew_bc(cr, v, numElements);
    nblock_bc(cr, numElements, *v);
}

template <typename T>
static void nblock_abc(const t_commrec *cr, int numElements, std::vector<T> *v)
{
    if (!MASTER(cr))
    {
        v->resize(numElements);
    }
    gmx_bcast(numElements*sizeof(T), v->data(), cr);
}

#endif
