/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
#ifndef GMX_PBCUTIL_BOXUTILITIES_H
#define GMX_PBCUTIL_BOXUTILITIES_H

#include <stdio.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct t_inputrec;
struct t_state;

/*! \brief Make sure the relative box shape remains the same
 *
 * This function ensures that the relative box dimensions are
 * preserved, which otherwise might diffuse away due to rounding
 * errors in pressure coupling or the deform option.
 *
 * \param[in] ir      Input record
 * \param[in] box_rel Relative box
 * \param[out] b      The corrected box
 */
void preserve_box_shape(t_inputrec *ir, matrix box_rel, matrix b);

/*! \brief Determine the relative box components
 *
 * Set state->box_rel used in mdrun to preserve the box shape
 * \param[in] ir       Input record
 * \param[inout] state Structure containing the box
 */
void set_box_rel(struct t_inputrec *ir, t_state *state);

#endif
