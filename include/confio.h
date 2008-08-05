/*
 * $Id$
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

#ifndef _confio_h
#define _confio_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"

/* For reading coordinate files it is assumed that enough memory
 * has been allocated beforehand.
 */
#ifdef CPLUSPLUS
extern "C" {
#endif
  
extern void init_t_atoms(t_atoms *atoms, int natoms, bool bPdbinfo);
/* allocate memory for the arrays, set nr to natoms and nres to 0
 * set pdbinfo to NULL or allocate memory for it */  

int read_g96_conf(FILE *fp,char *infile,t_trxframe *fr);
/* read a Gromos96 coordinate or trajectory file,                       *
 * returns the number of atoms                                          *
 * sets what's in the frame in info                                     *  
 * read from fp, infile is only needed for error messages               *   
 * nwanted is the number of wanted coordinates,                         *
 * set this to -1 if you want to know the number of atoms in the file   *
 * title, atoms, x, v can all be NULL, in which case they won't be read */

void write_g96_conf(FILE *out,t_trxframe *fr,int nindex,atom_id *index);
/* write a Gromos96 coordinate file or trajectory frame *
 * index can be NULL                                    */

extern bool gro_next_x_or_v(FILE *status,t_trxframe *fr);
extern int gro_first_x_or_v(FILE *status,t_trxframe *fr);
/* read first/next x and/or v frame from gro file */

extern void write_hconf_indexed_p(FILE *out,char *title,t_atoms *atoms,
				  int nx,atom_id index[],int ndec,
				  rvec *x,rvec *v,matrix box);
		
extern void write_hconf_p(FILE *out,char *title,t_atoms *atoms, int ndec,
			  rvec *x,rvec *v,matrix box); 
/* Write a Gromos file with precision ndec: number of decimal places in x,
 * v has one place more. */ 

void write_sto_conf_indexed(char *outfile,char *title,t_atoms *atoms, 
			    rvec x[],rvec *v,int ePBC,matrix box,
			    atom_id nindex,atom_id index[]);
/* like write_sto_conf, but indexed */ 

extern void write_sto_conf(char *outfile, char *title,t_atoms *atoms, 
			   rvec x[],rvec *v,int ePBC,matrix box);
/* write atoms, x, v (if .gro and not NULL) and box (if not NULL) 
 * to an STO (.gro or .pdb) file */ 

extern void write_sto_conf_mtop(char *outfile, char *title,gmx_mtop_t *mtop,
				rvec x[],rvec *v,int ePBC,matrix box);
/* As write_sto_conf, but uses a gmx_mtop_t struct */

extern void get_stx_coordnum (char *infile,int *natoms);
/* read the number of atoms from an STX file */

extern void read_stx_conf(char *infile, char *title,t_atoms *atoms, 
			  rvec x[],rvec *v,int *ePBC,matrix box);
/* Read atoms, x, v and box from an STX file.
 * If ePBC!=NULL return the type of pbc in *ePBC or -1 if unknown.
 */

#ifdef CPLUSPLUS
}
#endif

#endif	/* _confio_h */
