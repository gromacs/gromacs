/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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

#ifndef GMX_NBNXM_PAIRLISTSET_H
#define GMX_NBNXM_PAIRLISTSET_H

#include "gromacs/math/vectypes.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "locality.h"

/* Initializes a set of pair lists stored in nbnxn_pairlist_set_t
 *
 * TODO: Merge into the constructor
 */
typedef void nbnxn_free_t (void *ptr);

/* Tells if the pair-list corresponding to nb_kernel_type is simple.
 * Returns FALSE for super-sub type pair-list.
 */
gmx_bool nbnxn_kernel_pairlist_simple(int nb_kernel_type);

/* Initializes a set of pair lists stored in nbnxn_pairlist_set_t */
void nbnxn_init_pairlist_set(nbnxn_pairlist_set_t *nbl_list);

/*! \brief Prepare the list-set produced by the search for dynamic pruning
 *
 * \param[in,out] listSet  The list-set to prepare for dynamic pruning.
 */
void nbnxnPrepareListForDynamicPruning(nbnxn_pairlist_set_t *listSet);

#endif
