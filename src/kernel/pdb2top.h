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
 * Good gRace! Old Maple Actually Chews Slate
 */

#ifndef _pdb2top_h
#define _pdb2top_h

static char *SRCID_pdb2top_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) pdb2top.h 1.19 2/2/97"
#endif /* HAVE_IDENT */
#include "typedefs.h"
#include "pdb2gmx.h"
#include "h_db.h"
#include "toputil.h"
#include "gen_dum.h"
#include "ter_db.h"
#include "topexcl.h"

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

extern void print_top_comment(FILE *out, char *title, bool bITP);

extern void print_top_header(FILE *out, char *title, bool bITP, 
			     char *ff, real mHmult);

extern void print_top_mols(FILE *out, int nincl, char **incls,
			   int nmol, t_mols *mols);

extern void pdb2top(char *ff,FILE *top_file,char *posre_fn,char *molname,
		    int nincl, char **incls, int nmol, t_mols *mols,
		    t_atoms *atoms,int nah,t_addh ah[],rvec **x,
		    t_atomtype *atype,t_symtab *tab,
		    int nrb, t_resbond rb[],
		    int nrtp,t_restp rtp[],
		    int nra, t_resang ra[],
		    int nrd, t_resdih rd[],
		    int nid, t_idihres idi[],
		    t_hackblock *ntdb,t_hackblock *ctdb,
		    bool bH14,int rn,int rc,bool bAlldih,
		    bool bDummies,real mHmult,
		    int nssbonds,t_ssbond ssbonds[], int nrexcl);
/* Create a topology ! */

extern void write_top(char *ff,FILE *out,char *pr,char *molname,
		      int nincl, char **incls, int nmol, t_mols *mols,
		      t_atoms *at,t_params plist[],t_block *excl,
		      t_atomtype *atype,int *cgnr, int nrexcl, real mHmult);
/* write a topology 
 * NOTE: nrexcl is not the size of *excl! */

extern void print_sums(t_atoms *atoms, bool bSystem);

#endif	/* _pdb2top_h */
