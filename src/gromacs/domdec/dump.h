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
/*! \internal \file
 *
 * \brief This file declares functions for DD to write PDB files
 * e.g. when reporting problems.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_DUMP_H
#define GMX_DOMDEC_DUMP_H

#include <cstdint>

#include "gromacs/math/vectypes.h"

struct gmx_ddbox_t;
struct gmx_domdec_t;
struct gmx_mtop_t;
struct t_commrec;

//! Write the DD grid to a PDB file
void write_dd_grid_pdb(const char *fn, int64_t step,
                       gmx_domdec_t *dd, matrix box, gmx_ddbox_t *ddbox);

/*! \brief Dump a pdb file with the current DD home + communicated atoms.
 *
 * When natoms=-1, dump all known atoms.
 */
void write_dd_pdb(const char *fn, int64_t step, const char *title,
                  const gmx_mtop_t *mtop,
                  const t_commrec *cr,
                  int natoms, const rvec x[], const matrix box);

#endif
