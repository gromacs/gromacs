/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2014, by the GROMACS development team, led by
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
#ifndef GMX_INPUTREC_H
#define GMX_INPUTREC_H

#include "gromacs/legacyheaders/types/inputrec.h"

#ifdef __cplusplus
extern "C"
{
#endif
#if 0
} /* fixes auto-indentation problems */
#endif



int ir_optimal_nstcalcenergy(const t_inputrec *ir);

int tcouple_min_integration_steps(int etc);

int ir_optimal_nsttcouple(const t_inputrec *ir);

int pcouple_min_integration_steps(int epc);

int ir_optimal_nstpcouple(const t_inputrec *ir);

/* Returns if the Coulomb force or potential is switched to zero */
gmx_bool ir_coulomb_switched(const t_inputrec *ir);

/* Returns if the Coulomb interactions are zero beyond the rcoulomb.
 * Note: always returns TRUE for the Verlet cut-off scheme.
 */
gmx_bool ir_coulomb_is_zero_at_cutoff(const t_inputrec *ir);

/* As ir_coulomb_is_zero_at_cutoff, but also returns TRUE for user tabulated
 * interactions, since these might be zero beyond rcoulomb.
 */
gmx_bool ir_coulomb_might_be_zero_at_cutoff(const t_inputrec *ir);

/* Returns if the Van der Waals force or potential is switched to zero */
gmx_bool ir_vdw_switched(const t_inputrec *ir);

/* Returns if the Van der Waals interactions are zero beyond the rvdw.
 * Note: always returns TRUE for the Verlet cut-off scheme.
 */
gmx_bool ir_vdw_is_zero_at_cutoff(const t_inputrec *ir);

/* As ir_vdw_is_zero_at_cutoff, but also returns TRUE for user tabulated
 * interactions, since these might be zero beyond rvdw.
 */
gmx_bool ir_vdw_might_be_zero_at_cutoff(const t_inputrec *ir);

#ifdef __cplusplus
}
#endif


#endif
