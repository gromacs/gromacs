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
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_cdist_h = "$Id$";

#define HEAD_LOADED
#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "statutil.h"
#include "fatal.h"
#include "physics.h"
#include "string2.h"
#include "smalloc.h"
#include "macros.h"
#include "futil.h"
#include "copyrite.h"
#include "vec.h"
#include "pdbio.h"

typedef struct {
  real lb,ub,len;
} t_dist;

/* LIST OF FUNCTIONS */

extern void pdih_lengths(int ai,int aj,int ak,int al,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,real *ub,t_atoms *atoms);
			 
extern real idih_lengths(int ai,int aj,int ak,int al,
			 t_ilist ilist[],t_iparams iparams[],t_atoms *atoms);

extern bool dist_set(t_dist *d,int natoms,int ai,int aj);

extern void set_dist(t_dist *d,int natoms,int ai,int aj,real lb,real ub,
		     real len);
		     
extern void read_O_dist(void);

#define NOBOND -666
extern real lookup_bondlength(int ai,int aj,t_ilist ilist[],
			      t_iparams iparams[],bool bFail,t_atoms *atoms);
/* Return the bondlength in Angstrom, or NOBOND if not found.
 * If not found and bFail, a fatal_error will be issued 
 */
 
extern real lookup_angle(int ai,int aj,int ak,t_ilist ilist[],
			 t_iparams iparams[],t_atoms *atoms);
/* Return the angle between atoms in radians */

extern real angle_length(int ai,int aj,int ak,real theta,
			 t_ilist ilist[],t_iparams iparams[],t_atoms *atoms);
/* Return the distance corresponding to atoms in Angstrom 
 */
 			   
extern void gauche(int ai,int aj,int ak,int al,t_ilist ilist[],
		   t_iparams iparams[],real *lb,t_atoms *atoms);
		   
extern void gauche15(int ai,int aj,int ak,int al,int am,real omega1,
		     real omega2,real omega3,
		     t_ilist ilist[],t_iparams iparams[],real *lb,
		     t_atoms *atoms);

extern real d_len(t_dist *d,int natoms,int ai,int aj);
extern real d_ub(t_dist *d,int natoms,int ai,int aj);
extern real d_lb(t_dist *d,int natoms,int ai,int aj);

/* FUNCTIONS cdist.c calls */

void do_smooth (t_dist *d,t_atoms *atoms,real tol);
int do_triangle (t_dist *d,t_atoms *atoms,real tol);
int triangle_upper_bound (t_dist *d,int natoms,real tol);
int triangle_lower_bound (t_dist *d,int natoms,real tol);
/* Return the number of smoothing operations */

/* OTHER FUNCTIONS */
extern int set_virtual (int *ATMS,int N,real margin,t_dist *d,int natoms);

/* Routines from residues.c */
extern void peptide_bonds (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,
			   real weight[], real pep_margin,
			   t_ilist ilist[],t_iparams iparams[],bool bVir);

extern void arg (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
		 real arg_margin,t_ilist ilist[],
		 t_iparams iparams[],bool bVir);

extern void asn (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
		 real end_margin,t_ilist ilist[],
		 t_iparams iparams[],bool bVir);

extern void gln (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
		 real end_margin,t_ilist ilist[],
		 t_iparams iparams[],bool bVir);
		 
extern void phe (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
		 real ring_margin,t_ilist ilist[],
		 t_iparams iparams[],bool bVir);
		 
extern void tyr (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
		 real ring_margin,t_ilist ilist[],
		 t_iparams iparams[],bool bVir);

extern void trp (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
		 real ring_margin,t_ilist ilist[],
		 t_iparams iparams[],bool bVir);

extern void hisb(FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
		 real ring_margin,t_ilist ilist[],
		 t_iparams iparams[],bool bVir);
		  
extern void val (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
		 real val_margin,t_ilist ilist[],t_iparams iparams[]);

extern void leu (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
		 real leu_margin,t_ilist ilist[],t_iparams iparams[]);

extern void ile (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
		 real ile_margin,t_ilist ilist[],t_iparams iparams[]);
		 

