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
 * GRoups of Organic Molecules in ACtion for Science
 */
#ifndef _gen_dum_h
#define _gen_dum_h

static char *SRCID_gen_dum_h = "$Id$";

#include "typedefs.h"
#include "grompp.h"
#include "pdb2gmx.h"

typedef real t_dmbp[MAXFORCEPARAM];

typedef struct {
  char *mname; /* name of mass */
  char *mtype; /* type of mass */
  int  nm;     /* number of mass atoms */
  int  tp;     /* type of attachment (must be 1 for now) */
  t_dmbp c;    /* parameters for attachment, depend on type: */
  /*              1 : dCN, dNH, aCNH, mH */
} t_dummass; /* parameters for dummy masses */

typedef struct {
  char *na[4]; /* control atoms i,j,k,l */
  int  tp;     /* type of dummy (1,2,3) */
  real a;      /* parameter a */
  real b;      /* parameter b */
  real c;      /* parameter c */
} t_dumdum;  /* parameters for dummy atoms */

typedef struct {
  char *na[3]; /* control atoms i,j,k */
  t_dmbp c;    /* bond parameter */
} t_dumang;  /* parameters for angle constraints */

typedef struct {
  char *bname; /* atom type to search for */
  int  nmass;  /* number of masses to add */
  t_dummass *mass; /* mass data */
  int  ndum;   /* number of dummies to make */
  t_dumdum  *dum;  /* dummy data */
  int  nang;   /* number of angle constraints to add */
  t_dumang  *ang;  /* angle data */
} t_dumblock;

extern void do_dummies(t_atoms *at,t_atomtype *atype,t_symtab *symtab,
		       int nrtp,t_restp rtp[],rvec **x,
		       int nddb, t_dumblock *ddb, 
		       bool **is_dummy,real mHmult);

extern void clean_dum_angles(t_params *ps, t_params *plist, bool *is_dum);

extern void clean_dum_dihs(t_params *ps, int natom, char dihname[], 
			   t_params *plist, bool *is_dum);

extern void do_dum_top(t_params *psb, t_params *psd2, t_params *psd3, 
		       t_params *psd2FD, t_params *psd2FAD, t_params *psda,
		       int nddb,t_dumblock *ddb,bool is_dum[],
		       t_atoms *at, t_atomtype *atype, 
		       int nrtp, t_restp rtp[], real mHmult);

extern void do_dum_excl(t_block *excl, bool is_dum[]);

#endif	/* _gen_dum_h */
