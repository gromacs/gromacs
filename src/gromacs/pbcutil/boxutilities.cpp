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
/*! \internal \file
 * \brief
 * Implements routines in boxutilities.h.
 *
 * Utility functions for handling boxes.
 */
#include "gmxpre.h"

#include "boxutilities.h"

#include <cmath>

#include <algorithm>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

/*! \brief Change box components to preserve the relative box shape
 *
 * Change box components to box[XX][XX]*box_rel to preserve the relative box shape
 */
static void do_box_rel(t_inputrec *ir, matrix box_rel, matrix b, gmx_bool bInit)
{
    int d, d2;

    for (d = YY; d <= ZZ; d++)
    {
        for (d2 = XX; d2 <= (ir->epct == epctSEMIISOTROPIC ? YY : ZZ); d2++)
        {
            /* We need to check if this box component is deformed
             * or if deformation of another component might cause
             * changes in this component due to box corrections.
             */
            if (ir->deform[d][d2] == 0 &&
                !(d == ZZ && d2 == XX && ir->deform[d][YY] != 0 &&
                  (b[YY][d2] != 0 || ir->deform[YY][d2] != 0)))
            {
                if (bInit)
                {
                    box_rel[d][d2] = b[d][d2]/b[XX][XX];
                }
                else
                {
                    b[d][d2] = b[XX][XX]*box_rel[d][d2];
                }
            }
        }
    }
}

void preserve_box_shape(t_inputrec *ir, matrix box_rel, matrix b)
{
    if (inputrecPreserveShape(ir))
    {
        do_box_rel(ir, box_rel, b, FALSE);
    }
}

void set_box_rel(t_inputrec *ir, t_state *state)
{
    /* Make sure the box obeys the restrictions before we fix the ratios */
    correct_box(NULL, 0, state->box, NULL);

    clear_mat(state->box_rel);

    if (inputrecPreserveShape(ir))
    {
        do_box_rel(ir, state->box_rel, state->box, TRUE);
    }
}
