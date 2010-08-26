/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _pdbio_h
#define _pdbio_h

#include "sysstuff.h"
#include "typedefs.h"
#include "symtab.h"
#include "atomprop.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct gmx_conect_t *gmx_conect;

/* THE pdb format (for ATOM/HETATOM lines) */
static const char *pdbformat ="%-6s%5u  %-4.4s%3.3s %c%4d%c   %8.3f%8.3f%8.3f";
static const char *pdbformat4="%-6s%5u %-4.4s %3.3s %c%4d%c   %8.3f%8.3f%8.3f";

/* Enumerated type for pdb records. The other entries are ignored
 * when reading a pdb file 
 */
enum { epdbATOM,   epdbHETATM, epdbANISOU, epdbCRYST1, epdbCOMPND, 
       epdbMODEL,  epdbENDMDL, epdbTER,    epdbHEADER, epdbTITLE, epdbREMARK, 
       epdbCONECT, epdbNR };

/* Enumerated value for indexing an uij entry (anisotropic temperature factors) */
enum { U11, U22, U33, U12, U13, U23 };
       
void set_pdb_wide_format(gmx_bool bSet);
/* If bSet, use wider format for occupancy and bfactor */

void pdb_use_ter(gmx_bool bSet);
/* set read_pdbatoms to read upto 'TER' or 'ENDMDL' (default, bSet=FALSE).
   This function is fundamentally broken as far as thread-safety is concerned.*/

void gmx_write_pdb_box(FILE *out,int ePBC,matrix box);
/* write the box in the CRYST1 record,
 * with ePBC=-1 the pbc is guessed from the box
 * This function is fundamentally broken as far as thread-safety is concerned.
 */

void write_pdbfile_indexed(FILE *out,const char *title,t_atoms *atoms,
				  rvec x[],int ePBC,matrix box,char chain,
				  int model_nr,atom_id nindex,atom_id index[],
				  gmx_conect conect,gmx_bool bTerSepChains);
/* REALLY low level */

void write_pdbfile(FILE *out,const char *title,t_atoms *atoms,
			  rvec x[],int ePBC,matrix box,char chain,
			  int model_nr,gmx_conect conect,gmx_bool bTerSepChains);
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
  
void get_pdb_atomnumber(t_atoms *atoms,gmx_atomprop_t aps);
/* Routine to extract atomic numbers from the atom names */

int read_pdbfile(FILE *in,char *title,int *model_nr,
			t_atoms *atoms,rvec x[],int *ePBC,matrix box,
			gmx_bool bChange,gmx_conect conect);
/* Function returns number of atoms found.
 * ePBC and gmx_conect structure may be NULL.
 */

void read_pdb_conf(const char *infile,char *title, 
			  t_atoms *atoms,rvec x[],int *ePBC,matrix box,
			  gmx_bool bChange,gmx_conect conect);
/* Read a pdb file and extract ATOM and HETATM fields.
 * Read a box from the CRYST1 line, return 0 box when no CRYST1 is found.
 * Change atom names according to protein conventions if wanted.
 * ePBC and gmx_conect structure may be NULL.
 */

void get_pdb_coordnum(FILE *in,int *natoms);
/* Read a pdb file and count the ATOM and HETATM fields. */

gmx_bool is_hydrogen(const char *nm);
/* Return whether atom nm is a hydrogen */

gmx_bool is_dummymass(const char *nm);
/* Return whether atom nm is a dummy mass */

/* Routines to handle CONECT records if they have been read in */
void gmx_conect_dump(FILE *fp,gmx_conect conect);

gmx_bool gmx_conect_exist(gmx_conect conect,int ai,int aj);
/* Return TRUE if there is a conection between the atoms */

void gmx_conect_add(gmx_conect conect,int ai,int aj);
/* Add a connection between ai and aj (numbered from 0 to natom-1) */ 

gmx_conect gmx_conect_generate(t_topology *top);
/* Generate a conect structure from a topology */

gmx_conect gmx_conect_init();
/* Initiate data structure */

void gmx_conect_done(gmx_conect gc);
/* Free memory */

#ifdef __cplusplus
}
#endif

#endif	/* _pdbio_h */
