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
 * GRowing Old MAkes el Chrono Sweat
 */

#ifndef _pdb2top_h
#define _pdb2top_h

static char *SRCID_pdb2top_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) pdb2top.h 1.19 2/2/97"
#endif /* HAVE_IDENT */
#include "typedefs.h"
#include "toputil.h"
#include "hackblock.h"

/* this *MUST* correspond to array in pdb2top.c */
enum { ehisA, ehisB, ehisH, ehis1, ehisNR };
extern char *hh[ehisNR];

typedef struct {
  int  res1,res2;
  char *a1,*a2;
} t_ssbond;

typedef struct {
  char *name;
  int  nr;
} t_mols;

extern char *choose_ff(void);
/* Strange place for this function... */

extern void print_top_comment(FILE *out,char *filename,char *title,bool bITP);

extern void print_top_header(FILE *out,char *filename,char *title,bool bITP, 
			     char *ff,real mHmult);

extern void print_top_mols(FILE *out, char *title, int nincl, char **incls,
			   int nmol, t_mols *mols);

extern void pdb2top(FILE *top_file, char *posre_fn, char *molname,
		    t_atoms *atoms,rvec **x,
		    t_atomtype *atype,t_symtab *tab,
		    int bts[],
		    int nrtp,t_restp   rtp[],
		    t_hackblock *ntdb, t_hackblock *ctdb,
		    bool bH14, int nterpairs, int *rn, int *rc, bool bAlldih,
		    bool bDummies, bool bDummyAromatics, real mHmult,
		    int nssbonds, t_ssbond ssbonds[], int nrexcl, 
		    real long_bond_dist, real short_bond_dist,
		    bool bDeuterate);
/* Create a topology ! */
extern bool is_int(double x);
/* Returns TRUE when x is integer */

extern void print_sums(t_atoms *atoms, bool bSystem);

extern void write_top(FILE *out, char *pr,char *molname,
		      t_atoms *at,int bts[],t_params plist[],t_block *excl,
		      t_atomtype *atype,int *cgnr, int nrexcl);
/* NOTE: nrexcl is not the size of *excl! */

#endif	/* _pdb2top_h */
