/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
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
#ifndef GMX_TOPOLOGY_BLOCK_H
#define GMX_TOPOLOGY_BLOCK_H

#include "gromacs/legacyheaders/types/simple.h"

#ifdef __cplusplus
extern "C" {
#endif

/* the block structure points into an array (usually of atom_ids).
   It is a list of starting indices for objects of consecutive ids, such
   as molecules.
   For example, if this block denotes molecules, then the first molecule
   ranges from index[0] to index[1]-1 in the atom list.

   This makes the mapping from atoms to molecules O(Nmolecules) instead
   of O(Natoms) in size.  */
typedef struct t_block
{
    int      nr;           /* The number of blocks          */
    atom_id *index;        /* Array of indices (dim: nr+1)  */
    int      nalloc_index; /* The allocation size for index */
} t_block;

typedef struct t_blocka
{
    int      nr;    /* The number of blocks              */
    atom_id *index; /* Array of indices in a (dim: nr+1) */
    int      nra;   /* The number of atoms               */
    atom_id *a;     /* Array of atom numbers in each group  */
    /* (dim: nra)                           */
    /* Block i (0<=i<nr) runs from          */
    /* index[i] to index[i+1]-1. There will */
    /* allways be an extra entry in index   */
    /* to terminate the table               */
    int nalloc_index;           /* The allocation size for index        */
    int nalloc_a;               /* The allocation size for a            */
} t_blocka;

void init_block(t_block *block);
void init_blocka(t_blocka *block);
t_blocka *new_blocka(void);
/* allocate new block */

void done_block(t_block *block);
void done_blocka(t_blocka *block);

void copy_blocka(const t_blocka *src, t_blocka *dest);

void stupid_fill_block(t_block *grp, int natom, gmx_bool bOneIndexGroup);
/* Fill a block structure with numbers identical to the index
 * (0, 1, 2, .. natom-1)
 * If bOneIndexGroup, then all atoms are  lumped in one index group,
 * otherwise there is one atom per index entry
 */

void stupid_fill_blocka(t_blocka *grp, int natom);
/* Fill a block structure with numbers identical to the index
 * (0, 1, 2, .. natom-1)
 * There is one atom per index entry
 */

#ifdef __cplusplus
}
#endif

#endif
