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

extern void get_coordnum(char *infile, int *natoms);

extern void read_whole_conf(char *infile, char *title,t_atoms *atoms, 
			    rvec x[], rvec v[], matrix box);

extern void read_conf(char *infile,char *title,int *natoms,
		      rvec x[],rvec v[],matrix box);

extern bool gro_next_x(FILE *status,real *t,int natoms,rvec x[],matrix box);
extern int gro_first_x(FILE *status, real *t, rvec **x, matrix box);
/* read first/next x frame from gro file */

extern bool gro_next_x_or_v(FILE *status,real *t,int natoms,
			    rvec x[],rvec v[],matrix box);
extern int gro_first_x_or_v(FILE *status, real *t, 
			    rvec **x, rvec **v, matrix box);
/* read first/next x and/or v frame from gro file */

extern bool gro_next_v(FILE *status,real *t,int natoms,rvec v[],matrix box);
extern int gro_first_v(FILE *status, real *t, rvec **v, matrix box);
/* read first/next v frame from gro file */

extern bool gro_next_x_v(FILE *status,real *t,int natoms,
			 rvec x[],rvec v[],matrix box);
extern int gro_first_x_v(FILE *status, real *t, 
			 rvec **x, rvec **v, matrix box);
/* read first/next x and v frame from gro file */

extern void write_hconf(FILE *out,char *title,
			t_atoms *atoms,rvec *x, 
			rvec *v,matrix box);
			
extern void write_hconf_p(FILE *out,char *title,t_atoms *atoms, int pr,
			  rvec *x,rvec *v,matrix box);
			
extern void write_hconf_indexed(FILE *out,char *title,t_atoms *atoms,
				int nx,atom_id index[],
				rvec *x,rvec *v,matrix box);

extern void write_hconf_p(FILE *out,char *title,t_atoms *atoms, int pr,
			  rvec *x,rvec *v,matrix box); 
/* Write a Gromos file with precision pr: number of decimal places in x,
 * v has one place more. */ 

extern void write_conf(char *outfile,char *title,t_atoms *atoms,
		       rvec *x,rvec *v,matrix box);
/* For three write_conf routines, if v == NULL, it is not written */

extern void write_pdb_conf(char *outfile,t_atoms *atoms,rvec x[],matrix box,
			   bool bChange);
/* Change atom names according to protein conventions if wanted */

extern void write_pdb_confs(char *outfile,t_atoms **atoms,rvec *x[],
			    int number);

extern void hwrite_pdb_conf_indexed(FILE *out,t_atoms *atoms,rvec x[],
				    matrix box, 
				    int gnx,atom_id index[]);
extern void write_pdb_conf_indexed(char *outfile,t_atoms *atoms,rvec x[],
				   matrix box,
				   int gnx,atom_id index[]);
/* Write a pdb file to either FILE *out or first open the file outfile
 * Use an index to only write out selected atoms. */
 
extern void read_pdb_conf(char *infile,t_atoms *atoms, 
			      rvec x[], matrix box);
		   
extern void write_xdr_conf(char *outfile,char *title,t_atoms *atoms,rvec x[],rvec v[],matrix box);

extern void read_xdr_coordnum(char *infile,int *natoms);

extern void read_xdr_conf(char *infile,char *title,t_atoms *atoms,rvec x[],rvec v[],matrix box);

extern void get_stx_coordnum (char *infile,int *natoms);

extern void read_stx_conf(char *infile, char *title,t_atoms *atoms, 
		   rvec x[],rvec v[], matrix box);

#ifdef CPLUSPLUS
}
#endif

#endif	/* _confio_h */
