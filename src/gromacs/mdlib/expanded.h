/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_EXPANDED_H
#define GMX_MDLIB_EXPANDED_H

#include <stdio.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

struct df_history_t;
struct gmx_enerdata_t;
struct t_expanded;
struct t_extmass;
struct t_inputrec;
struct t_lambda;
struct t_mdatoms;
struct t_simtemp;
class t_state;

void init_npt_masses(const t_inputrec *ir, t_state *state, t_extmass *MassQ, gmx_bool bInit);

void init_expanded_ensemble(gmx_bool bStateFromCP, const t_inputrec *ir, df_history_t *dfhist);

int ExpandedEnsembleDynamics(FILE *log, const t_inputrec *ir, const gmx_enerdata_t *enerd,
                             t_state *state, t_extmass *MassQ, int fep_state, df_history_t *dfhist,
                             int64_t step,
                             rvec *v, const t_mdatoms *mdatoms);

void PrintFreeEnergyInfoToFile(FILE *outfile, const t_lambda *fep, const t_expanded *expand,
                               const t_simtemp *simtemp, const df_history_t *dfhist,
                               int fep_state, int frequency, int64_t step);

#endif
