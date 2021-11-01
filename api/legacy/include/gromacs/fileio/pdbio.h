/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020,2021, by the GROMACS development team, led by
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

#ifndef GMX_FILEIO_PDBIO_H
#define GMX_FILEIO_PDBIO_H

#include <stdio.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

class AtomProperties;
struct t_atoms;
struct t_symtab;
struct t_topology;
enum class PbcType : int;
enum class PdbRecordType : int;

typedef struct gmx_conect_t* gmx_conect;

/* Write a PDB line with an ATOM or HETATM record directly to a file pointer.
 *
 * Returns the number of characters printed.
 */
int gmx_fprintf_pdb_atomline(FILE*         fp,
                             PdbRecordType record,
                             int           atom_seq_number,
                             const char*   atom_name,
                             char          alternate_location,
                             const char*   res_name,
                             char          chain_id,
                             int           res_seq_number,
                             char          res_insertion_code,
                             real          x,
                             real          y,
                             real          z,
                             real          occupancy,
                             real          b_factor,
                             const char*   element);

/* Enumerated value for indexing an uij entry (anisotropic temperature factors) */
enum
{
    U11,
    U22,
    U33,
    U12,
    U13,
    U23
};

void pdb_use_ter(gmx_bool bSet);
/* set read_pdbatoms to read upto 'TER' or 'ENDMDL' (default, bSet=FALSE).
   This function is fundamentally broken as far as thread-safety is concerned.*/

void gmx_write_pdb_box(FILE* out, PbcType pbcType, const matrix box);
/* write the box in the CRYST1 record,
 * with pbcType=PbcType::Unset the pbc is guessed from the box
 * This function is fundamentally broken as far as thread-safety is concerned.
 */

void write_pdbfile_indexed(FILE*          out,
                           const char*    title,
                           const t_atoms* atoms,
                           const rvec     x[],
                           PbcType        pbcType,
                           const matrix   box,
                           char           chain,
                           int            model_nr,
                           int            nindex,
                           const int      index[],
                           gmx_conect     conect,
                           bool           usePqrFormat,
                           bool           standardCompliantMode = false);
/* REALLY low level */

void write_pdbfile(FILE*          out,
                   const char*    title,
                   const t_atoms* atoms,
                   const rvec     x[],
                   PbcType        pbcType,
                   const matrix   box,
                   char           chain,
                   int            model_nr,
                   gmx_conect     conect);
/* Low level pdb file writing routine.
 *
 *          ONLY FOR SPECIAL PURPOSES,
 *
 *       USE write_sto_conf WHEN YOU CAN.
 *
 * override chain-identifiers with chain when chain>0
 * write ENDMDL when bEndmodel is TRUE.
 *
 * If the gmx_conect structure is not NULL its content is dumped as CONECT records
 * which may be useful for visualization purposes.
 */

void get_pdb_atomnumber(const t_atoms* atoms, AtomProperties* aps);
/* Routine to extract atomic numbers from the atom names */

int read_pdbfile(FILE*            in,
                 char*            title,
                 int*             model_nr,
                 struct t_atoms*  atoms,
                 struct t_symtab* symtab,
                 rvec             x[],
                 PbcType*         pbcType,
                 matrix           box,
                 gmx_conect       conect);
/* Function returns number of atoms found.
 * pbcType and gmx_conect structure may be NULL.
 */

void gmx_pdb_read_conf(const char* infile,
                       t_symtab*   symtab,
                       char**      name,
                       t_atoms*    atoms,
                       rvec        x[],
                       PbcType*    pbcType,
                       matrix      box);
/* Read a pdb file and extract ATOM and HETATM fields.
 * Read a box from the CRYST1 line, return 0 box when no CRYST1 is found.
 * pbcType may be NULL.
 *
 * If name is not nullptr, gmx_strdup the title string into it. */

void get_pdb_coordnum(FILE* in, int* natoms);
/* Read a pdb file and count the ATOM and HETATM fields. */

gmx_bool is_hydrogen(const char* nm);
/* Return whether atom nm is a hydrogen */

gmx_bool is_dummymass(const char* nm);
/* Return whether atom nm is a dummy mass */

/* Routines to handle CONECT records if they have been read in */
void gmx_conect_dump(FILE* fp, gmx_conect conect);

gmx_bool gmx_conect_exist(gmx_conect conect, int ai, int aj);
/* Return TRUE if there is a conection between the atoms */

void gmx_conect_add(gmx_conect conect, int ai, int aj);
/* Add a connection between ai and aj (numbered from 0 to natom-1) */

gmx_conect gmx_conect_generate(const t_topology* top);
/* Generate a conect structure from a topology */

gmx_conect gmx_conect_init();
/* Initiate data structure */

void gmx_conect_done(gmx_conect gc);
/* Free memory */

#endif /* GMX_FILEIO_PDBIO_H */
