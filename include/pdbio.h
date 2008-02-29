/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */

#ifndef _pdbio_h
#define _pdbio_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sysstuff.h"
#include "typedefs.h"
#include "symtab.h"

/* THE pdb format (for ATOM/HETATOM lines) */
static char *pdbformat ="%-6s%5u  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f";
static char *pdbformat4="%-6s%5u %-4.4s %3.3s %c%4d    %8.3f%8.3f%8.3f";

/* Enumerated type for pdb records. The other entries are ignored
 * when reading a pdb file 
 */
enum { epdbATOM,   epdbHETATM, epdbANISOU, epdbCRYST1, epdbCOMPND, 
       epdbMODEL,  epdbENDMDL, epdbTER,    epdbHEADER, epdbTITLE, epdbREMARK, 
       epdbCONECT, epdbNR };

/* Enumerated value for indexing an uij entry (anisotropic temperature factors) */
enum { U11, U22, U33, U12, U13, U23 };
       
extern void set_pdb_wide_format(bool bSet);
/* If bSet, use wider format for occupancy and bfactor */

extern void pdb_use_ter(bool bSet);
/* set read_pdbatoms to read upto 'TER' or 'ENDMDL' (default, bSet=FALSE) */

extern void write_pdbfile_indexed(FILE *out,char *title,t_atoms *atoms,
				  rvec x[],matrix box,char chain,
				  int model_nr,atom_id nindex,atom_id index[]);
/* REALLY low level */

extern void write_pdbfile(FILE *out,char *title,t_atoms *atoms,
			  rvec x[],matrix box,char chain,int model_nr);
/* Low level pdb file writing routine.
 * 
 *          ONLY FOR SPECIAL PURPOSES,
 * 
 *       USE write_sto_conf WHEN YOU CAN.
 *
 * override chain-identifiers with chain when chain>0
 * write ENDMDL when bEndmodel is TRUE */
  
extern int read_pdbfile(FILE *in,char *title,int *model_nr,
			t_atoms *atoms,rvec x[],matrix box,bool bChange);
/* Function returns number of atoms found. */

extern void read_pdb_conf(char *infile,char *title, 
			  t_atoms *atoms,rvec x[],matrix box,bool bChange);
/* Read a pdb file and extract ATOM and HETATM fields.
 * Read a box from the CRYST1 line, return 0 box when no CRYST1 is found.
 * Change atom names according to protein conventions if wanted
 */

extern void get_pdb_coordnum(FILE *in,int *natoms);
/* Read a pdb file and count the ATOM and HETATM fields. */

extern bool is_hydrogen(char *nm);
/* Return whether atom nm is a hydrogen */

extern bool is_dummymass(char *nm);
/* Return whether atom nm is a dummy mass */

#endif	/* _pdbio_h */
