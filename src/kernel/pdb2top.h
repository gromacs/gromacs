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
  int res1,res2;
  char *a1,*a2;
} t_ssbond;

extern void pdb2top(char *ff,char *fn,char *pr,char *title,char *molname,
		    int nincl, char **incls, int nmol, char **mols,
		    t_atoms *atoms,int nah,t_addh ah[],rvec **x,
		    t_atomtype *atype,t_symtab *tab,
		    int nrb, t_resbond rb[],
		    int nrtp,t_restp rtp[],
		    int nra, t_resang ra[],
		    int nid, t_idihres idi[],
		    t_hackblock *ntdb,t_hackblock *ctdb,
		    bool bH14,int rn,int rc,bool bAlldih,
		    int nddb, t_dumblock *ddb,bool bDummies,real mHmult,
		    int nssbonds,t_ssbond ssbonds[], int nrexcl);
/* Create a topology ! */

extern void write_top(char *ff,char *fn, char *pr,char *title, char *molname,
		      int nincl, char **incls, int nmol, char **mols,
		      t_atoms *at,t_params plist[],t_block *excl,
		      t_atomtype *atype,t_block *cgs, int nrexcl, real mHmult);
/* write a topology 
 * NOTE: nrexcl is not the size of *excl! */

#endif	/* _pdb2top_h */
