/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2005,2006,2007,2008,2009,2010,2012,2013,2014,2015,2018, by the GROMACS development team, led by
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
 * \brief This file declares functions for domdec to use
 * while managing inter-atomic constraints.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_DOMDEC_CONSTRAINTS_H
#define GMX_DOMDEC_DOMDEC_CONSTRAINTS_H

namespace gmx
{
class Constraints;
}

struct gmx_domdec_t;
struct gmx_mtop_t;
struct t_ilist;

/*! \brief Clears the local indices for the constraint communication setup */
void dd_clear_local_constraint_indices(gmx_domdec_t *dd);

/*! \brief Sets up communication and atom indices for all local+connected constraints */
int dd_make_local_constraints(struct gmx_domdec_t *dd, int at_start,
                              const struct gmx_mtop_t *mtop,
                              const int *cginfo,
                              gmx::Constraints *constr, int nrec,
                              struct t_ilist *il_local);

/*! \brief Initializes the data structures for constraint communication */
void init_domdec_constraints(gmx_domdec_t     *dd,
                             const gmx_mtop_t *mtop);

#endif
