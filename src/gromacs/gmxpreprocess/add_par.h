/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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

#ifndef GMX_GMXPREPROCESS_ADD_PAR_H
#define GMX_GMXPREPROCESS_ADD_PAR_H

#include "gromacs/utility/real.h"

struct InteractionsOfType;
struct PreprocessResidue;

namespace gmx
{
template<typename>
class ArrayRef;
}

void add_param(InteractionsOfType* ps, int ai, int aj, gmx::ArrayRef<const real> c, const char* s);

void add_cmap_param(InteractionsOfType* ps, int ai, int aj, int ak, int al, int am, const char* s);

void add_vsite3_atoms(InteractionsOfType* ps, int ai, int aj, int ak, int al, bool bSwapParity);

void add_vsite2_param(InteractionsOfType* ps, int ai, int aj, int ak, real c0);

void add_vsite3_param(InteractionsOfType* ps, int ai, int aj, int ak, int al, real c0, real c1);

void add_vsite4_atoms(InteractionsOfType* ps, int ai, int aj, int ak, int al, int am);

int search_jtype(const PreprocessResidue& localPpResidue, const char* name, bool bFirstRes);

#endif
