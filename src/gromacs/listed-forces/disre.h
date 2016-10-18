/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares functions for handling distance restraints.
 *
 * \inlibraryapi
 * \ingroup module_listed-forces
 */
#ifndef GMX_LISTED_FORCES_DISRE_H
#define GMX_LISTED_FORCES_DISRE_H

#include <cstdio>

#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_mtop_t;
struct history_t;
struct t_commrec;
struct t_inputrec;
struct t_pbc;
struct t_state;

/*! \brief
 * Initiates *fcd data.
 *
 * Must be called once, nbonds is the number
 * of iatoms in the ilist of the idef struct.
 * When time averaging is used, the history is initialized in state,
 * unless it was read before from a checkpoint file.
 * The implementation of distance restraints with -multi
 * must differ according to whether REMD is active.
 */
void init_disres(FILE *fplog, const gmx_mtop_t *mtop,
                 t_inputrec *ir, const t_commrec *cr,
                 t_fcdata *fcd, t_state *state, gmx_bool bIsREMD);

/*! \brief
 * Calculates r and r^-3 (inst. and time averaged) for all pairs
 * and the ensemble averaged r^-6 (inst. and time averaged) for all restraints
 */
void calc_disres_R_6(const t_commrec *cr,
                     int nfa, const t_iatom *fa,
                     const rvec *x, const t_pbc *pbc,
                     t_fcdata *fcd, history_t *hist);

//! Calculates the distance restraint forces, return the potential.
t_ifunc ta_disres;

//! Copies the new time averages that have been calculated in calc_disres_R_6.
void update_disres_history(t_fcdata *fcd, history_t *hist);

#endif
