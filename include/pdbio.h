/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _pdbio_h
#define _pdbio_h

static char *SRCID_pdbio_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) pdbio.h 1.12 7/28/97"
#endif /* HAVE_IDENT */
#include "sysstuff.h"
#include "typedefs.h"
#include "symtab.h"

/* THE pdb format (for ATOM/HETATOM lines) */
static char *pdbformat ="%-6s%5u  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f";
static char *pdbformat4="%-6s%5u %-4.4s %3.3s %c%4d    %8.3f%8.3f%8.3f";

/* Enumerated type for pdb records. The other entries are ignored
 * when reading a pdb file 
 */
enum { epdbATOM, epdbHETATM, epdbANISOU, epdbCRYST1, epdbCOMPND, 
       epdbMODEL, epdbENDMDL, epdbTER, epdbHEADER, epdbTITLE, epdbREMARK, 
       epdbNR };

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
