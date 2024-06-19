/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
#ifndef GMX_MDLIB_WALLS_H
#define GMX_MDLIB_WALLS_H

#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

struct SimulationGroups;
struct t_forcerec;
struct t_inputrec;
struct t_mdatoms;
struct t_nrnb;

namespace gmx
{
template<typename>
class ArrayRef;
class ForceWithVirial;
} // namespace gmx

void make_wall_tables(FILE*                   fplog,
                      const t_inputrec&       ir,
                      const char*             tabfn,
                      const SimulationGroups* groups,
                      t_forcerec*             fr);

real do_walls(const t_inputrec&                   ir,
              const t_forcerec&                   fr,
              const matrix                        box,
              gmx::ArrayRef<const int>            typeA,
              gmx::ArrayRef<const int>            typeB,
              gmx::ArrayRef<const unsigned short> cENER,
              int                                 homenr,
              int                                 numPerturbedAtoms,
              gmx::ArrayRef<const gmx::RVec>      x,
              gmx::ForceWithVirial*               forceWithVirial,
              real                                lambda,
              gmx::ArrayRef<real>                 Vlj,
              t_nrnb*                             nrnb);

#endif
