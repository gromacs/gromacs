/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#ifndef GMX_GMX_NM2TYPE_H
#define GMX_GMX_NM2TYPE_H

#include <cstdio>

class PreprocessingAtomTypes;
struct t_atoms;
struct InteractionsOfType;
struct t_symtab;

struct t_nm2type
{
    char *  elem, *type;
    double  q, m;
    int     nbonds;
    char**  bond;
    double* blen;
};

t_nm2type* rd_nm2type(const char* ffdir, int* nnm);
/* Read the name 2 type database. nnm is the number of entries
 * ff is the force field.
 */

void dump_nm2type(FILE* fp, int nnm, t_nm2type nm2t[]);
/* Dump the database for debugging. Can be reread by the program */

int nm2type(int nnm, t_nm2type nm2t[], t_atoms* atoms, PreprocessingAtomTypes* atype, int* nbonds, InteractionsOfType* bond);
/* Try to determine the atomtype (force field dependent) for the atoms
 * with help of the bond list
 */

#endif
