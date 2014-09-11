/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#include "gromacs/legacyheaders/types/simple.h"

#ifdef __cplusplus
extern "C" {
#endif

struct gmx_atomprop;
struct t_atoms;
struct t_topology;

typedef struct gmx_conect_t *gmx_conect;

/* Enumerated type for pdb records. The other entries are ignored
 * when reading a pdb file
 */
enum PDB_record {
    epdbATOM,   epdbHETATM, epdbANISOU, epdbCRYST1, epdbCOMPND,
    epdbMODEL,  epdbENDMDL, epdbTER,    epdbHEADER, epdbTITLE, epdbREMARK,
    epdbCONECT, epdbNR
};

/* Write a PDB line with an ATOM or HETATM record directly to a file pointer.
 *
 * Returns the number of characters printed.
 */
int
gmx_fprintf_pdb_atomline(FILE *            fp,
                         enum PDB_record   record,
                         int               atom_seq_number,
                         const char *      atom_name,
                         char              alternate_location,
                         const char *      res_name,
                         char              chain_id,
                         int               res_seq_number,
                         char              res_insertion_code,
                         real              x,
                         real              y,
                         real              z,
                         real              occupancy,
                         real              b_factor,
                         const char *      element);

/* Enumerated value for indexing an uij entry (anisotropic temperature factors) */
enum {
    U11, U22, U33, U12, U13, U23
};

void pdb_use_ter(gmx_bool bSet);
/* set read_pdbatoms to read upto 'TER' or 'ENDMDL' (default, bSet=FALSE).
   This function is fundamentally broken as far as thread-safety is concerned.*/

void gmx_write_pdb_box(FILE *out, int ePBC, matrix box);
/* write the box in the CRYST1 record,
 * with ePBC=-1 the pbc is guessed from the box
 * This function is fundamentally broken as far as thread-safety is concerned.
 */

void write_pdbfile_indexed(FILE *out, const char *title, struct t_atoms *atoms,
                           rvec x[], int ePBC, matrix box, char chain,
                           int model_nr, atom_id nindex, const atom_id index[],
                           gmx_conect conect, gmx_bool bTerSepChains);
/* REALLY low level */

void write_pdbfile(FILE *out, const char *title, struct t_atoms *atoms,
                   rvec x[], int ePBC, matrix box, char chain,
                   int model_nr, gmx_conect conect, gmx_bool bTerSepChains);
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

void get_pdb_atomnumber(struct t_atoms *atoms, struct gmx_atomprop *aps);
/* Routine to extract atomic numbers from the atom names */

int read_pdbfile(FILE *in, char *title, int *model_nr,
                 struct t_atoms *atoms, rvec x[], int *ePBC, matrix box,
                 gmx_bool bChange, gmx_conect conect);
/* Function returns number of atoms found.
 * ePBC and gmx_conect structure may be NULL.
 */

void read_pdb_conf(const char *infile, char *title,
                   struct t_atoms *atoms, rvec x[], int *ePBC, matrix box,
                   gmx_bool bChange, gmx_conect conect);
/* Read a pdb file and extract ATOM and HETATM fields.
 * Read a box from the CRYST1 line, return 0 box when no CRYST1 is found.
 * Change atom names according to protein conventions if wanted.
 * ePBC and gmx_conect structure may be NULL.
 */

void get_pdb_coordnum(FILE *in, int *natoms);
/* Read a pdb file and count the ATOM and HETATM fields. */

gmx_bool is_hydrogen(const char *nm);
/* Return whether atom nm is a hydrogen */

gmx_bool is_dummymass(const char *nm);
/* Return whether atom nm is a dummy mass */

/* Routines to handle CONECT records if they have been read in */
void gmx_conect_dump(FILE *fp, gmx_conect conect);

gmx_bool gmx_conect_exist(gmx_conect conect, int ai, int aj);
/* Return TRUE if there is a conection between the atoms */

void gmx_conect_add(gmx_conect conect, int ai, int aj);
/* Add a connection between ai and aj (numbered from 0 to natom-1) */

gmx_conect gmx_conect_generate(struct t_topology *top);
/* Generate a conect structure from a topology */

gmx_conect gmx_conect_init(void);
/* Initiate data structure */

void gmx_conect_done(gmx_conect gc);
/* Free memory */

#ifdef __cplusplus
}
#endif

#endif  /* GMX_FILEIO_PDBIO_H */
