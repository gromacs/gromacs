/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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
#ifndef GMX_GMXPREPROCESS_READ_CONFORMATION_H
#define GMX_GMXPREPROCESS_READ_CONFORMATION_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

struct gmx_atomprop;
struct gmx_mtop_t;
struct t_atoms;
struct t_topology;

/*! \brief Allocate and fill an array of inter-atomic half distances
 *
 * These are either scaled VDW radii taken from vdwradii.dat, or the
 * default value. Used directly and indirectly by solvate and
 * insert-molecules for deciding whether molecules clash. The return
 * pointer should be freed by the caller. */
std::vector<real>
makeExclusionDistances(const t_atoms *a, gmx_atomprop *aps,
                       real defaultDistance, real scaleFactor);

void readConformation(const char *confin, gmx_mtop_t *top,
                      std::vector<gmx::RVec> *x, std::vector<gmx::RVec> *v,
                      int *ePBC, matrix box, const char *statusTitle);
/*! \brief Read a conformation from a file, allocate and fill data structures.
 *
 * Used by solvate and insert-molecules. The returned pointers *x and
 * *v should be freed by the caller. top should have its destructor
 * called. */
void readConformation(const char *confin, t_topology *top,
                      std::vector<gmx::RVec> *x, std::vector<gmx::RVec> *v,
                      int *ePBC, matrix box, const char *statusTitle);

#endif
