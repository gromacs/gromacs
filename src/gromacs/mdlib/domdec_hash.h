/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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
#ifndef GMX_DOMDEC_HASH_H
#define GMX_DOMDEC_HASH_H

#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

/* This include file implements the simplest hash table possible.
 * It is limited to integer keys and integer values.
 * The purpose is highest efficiency and lowest memory usage possible.
 */

typedef struct gmx_hash_t gmx_hash_t;

/* Clear all the entries in the hash table.
 * With the current number of keys check if the table size is still good,
 * if not optimize it with the currenr number of keys.
 */
void gmx_hash_clear_and_optimize(gmx_hash_t *hash);

gmx_hash_t *gmx_hash_init(int nkey_used_estimate);

/* Set the hash entry for global atom a_gl to local atom a_loc and cell. */
void gmx_hash_set(gmx_hash_t *hash, int key, int value);

/* Change the hash value if already set, otherwise set the hash value */
void gmx_hash_change_or_set(gmx_hash_t *hash, int key, int value);

/* Returns the value or -1 if the key is not present */
int gmx_hash_get_minone(const gmx_hash_t *hash, int key);

#ifdef __cplusplus
}
#endif

#endif
