/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */

#ifndef _confio_h
#define _confio_h

static char *SRCID_confio_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) confio.h 1.16 11/23/92"
#endif /* HAVE_IDENT */

#include "typedefs.h"

/* For reading coordinate files it is assumed that enough memory
 * has been allocated beforehand.
 */
#ifdef CPLUSPLUS
extern "C" {
#endif
  
  typedef struct {
    bool bTitle;
    bool bTime;
    bool bAtoms;
    bool bPos;
    bool bVel;
    bool bBox;
    int  step;
    real time;
  } t_g96info;
  
extern void init_t_atoms(t_atoms *atoms, int natoms, bool bPdbinfo);
/* allocate memory for the arrays, set nr to natoms and nres to 0
 * set pdbinfo to NULL or allocate memory for it */  

extern void free_t_atoms(t_atoms *atoms);
/* free all the arrays and set the nr and nres to 0 */

void clear_g96info(t_g96info *info);
/* set all bools in the info struct to FALSE */

int read_g96_conf(FILE *fp,char *infile,int nwanted,t_g96info *info,
		  char *title,t_atoms *atoms,rvec *x, rvec *v,matrix box);
/* read a Gromos96 coordinate or trajectory file,                       *
 * returns the number of atoms                                          *
 * sets what's in the frame in info                                     *  
 * read from fp, infile is only needed for error messages               *   
 * nwanted is the number of wanted coordinates,                         *
 * set this to -1 if you want to know the number of atoms in the file   *
 * title, atoms, x, v can all be NULL, in which case they won't be read */

void write_g96_conf(FILE *out,char *title,t_atoms *atoms,
                    rvec *x,rvec *v,matrix box,
                    int nindex,atom_id *index);
/* write a Gromos96 coordinate file *
 * x, v and index can be NULL       */

extern bool gro_next_x(FILE *status,real *t,int natoms,rvec x[],matrix box);
extern int gro_first_x(FILE *status, real *t, rvec **x, matrix box);
/* read first/next x frame from gro file */

extern bool gro_next_x_or_v(FILE *status,real *t,int natoms,
			    rvec x[],rvec *v,matrix box);
extern int gro_first_x_or_v(FILE *status, real *t, 
			    rvec **x, rvec **v, matrix box);
/* read first/next x and/or v frame from gro file */

extern bool gro_next_v(FILE *status,real *t,int natoms,rvec *v,matrix box);
extern int gro_first_v(FILE *status, real *t, rvec **v, matrix box);
/* read first/next v frame from gro file */

extern bool gro_next_x_v(FILE *status,real *t,int natoms,
			 rvec x[],rvec *v,matrix box);
extern int gro_first_x_v(FILE *status, real *t, 
			 rvec **x, rvec **v, matrix box);
/* read first/next x and v frame from gro file */

extern void write_hconf(FILE *out,char *title,
			t_atoms *atoms,rvec *x, 
			rvec *v,matrix box);
			
extern void write_hconf_indexed(FILE *out,char *title,t_atoms *atoms,
				int nx,atom_id index[],
				rvec *x,rvec *v,matrix box);

extern void write_hconf_p(FILE *out,char *title,t_atoms *atoms, int pr,
			  rvec *x,rvec *v,matrix box); 
/* Write a Gromos file with precision pr: number of decimal places in x,
 * v has one place more. */ 

extern void write_xdr_conf(char *outfile,char *title,t_atoms *atoms,
			   rvec x[],rvec *v,matrix box);

extern void read_xdr_coordnum(char *infile,int *natoms);

extern void read_xdr_conf(char *infile,char *title,t_atoms *atoms,
			  rvec x[],rvec *v,matrix box);

void write_sto_conf_indexed(char *outfile,char *title,t_atoms *atoms, 
			    rvec x[],rvec *v,matrix box,
			    atom_id nindex,atom_id index[]);
/* like write_sto_conf, but indexed */ 

extern void write_sto_conf(char *outfile, char *title,t_atoms *atoms, 
			   rvec x[],rvec *v, matrix box);
/* write atoms, x, v (if .gro and not NULL) and box (if not NULL) 
 * to an STO (.gro or .pdb) file */ 

extern void get_stx_coordnum (char *infile,int *natoms);
/* read the number of atoms from an STX file */

extern void read_stx_conf(char *infile, char *title,t_atoms *atoms, 
			  rvec x[],rvec *v, matrix box);
/* read atoms, x, v and box from an STX file */

#ifdef CPLUSPLUS
}
#endif

#endif	/* _confio_h */
