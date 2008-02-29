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

#ifndef _hackblock_h
#define _hackblock_h

#include "typedefs.h"
#include "pdbio.h"
#include "grompp.h"

/* Used for reading .rtp/.tdb */
/* ebtsBONDS must be the first, new types can be added to the end */
/* these *MUST* correspond to the arrays in hackblock.c */
enum { ebtsBONDS, ebtsANGLES, ebtsPDIHS, ebtsIDIHS, ebtsEXCLS, ebtsNR };
extern const char *btsNames[ebtsNR];
extern const int btsNiatoms[ebtsNR];

/* if changing any of these structs, make sure that all of the
   free/clear/copy/merge_t_* functions stay updated */

/* BONDEDS */
typedef struct {
  char 	*a[MAXATOMLIST]; /* atom names */
  char  *s;              /* optional define string which gets copied from
			    .rtp/.tdb to .top and will be parsed by cpp
			    during grompp */
} t_rbonded;

typedef struct {
  int       nb;          /* number of bondeds */
  t_rbonded *b;          /* bondeds */
} t_rbondeds;

/* RESIDUES (rtp) */
typedef struct {
  char   *resname;
  /* atom data */
  int    natom;
  t_atom *atom;
  char   ***atomname;
  int    *cgnr;
  /* list of bonded interactions to add */
  t_rbondeds rb[ebtsNR];
} t_restp;

/* Block to hack residues */
typedef struct {
  int     nr;       /* Number of atoms to hack    */
  char    *oname;   /* Old name                   */
  char	  *nname;   /* New name                   */
  /* the type of hack depends on the setting of oname and nname:
   * if oname==NULL                we're adding, must have tp>0 also!
   * if oname!=NULL && nname==NULL we're deleting
   * if oname!=NULL && nname!=NULL we're replacing
   */
  t_atom  *atom;    /* New atom data              */
  int     cgnr;     /* chargegroup number. if not read will be NOTSET */
  int     tp;       /* Type of attachment (1..11) */
  int     nctl;     /* How many control atoms there are */
  char 	  *a[4];    /* Control atoms i,j,k,l	  */
  rvec    newx;     /* calculated new position    */
  atom_id newi;     /* new atom index number (after additions) */
} t_hack;

typedef struct {
  char      *name;  /* Name of hack block (residue or terminus) */
  int       nhack;  /* Number of atoms to hack                  */
  int       maxhack;/* used for efficient srenew-ing            */
  t_hack    *hack;  /* Hack list                                */
  /* list of bonded interactions to add */
  t_rbondeds rb[ebtsNR];
} t_hackblock;

/* all libraries and other data to protonate a structure or trajectory */
typedef struct {
  bool        bInit; /* true after init; set false by init_t_protonate */
  /* force field name: */
  char        FF[10];
  /* libarary data: */
  int         *nab;
  t_hack      **ab;
  t_hackblock *ah, *ntdb, *ctdb;
  t_hackblock **sel_ntdb, **sel_ctdb;
  int         nah;
  t_symtab    tab;
  int         *rN, *rC;
  t_atomtype  *atype;
  /* protonated topology: */
  t_atoms     *patoms;
  /* unprotonated topology: */
  t_atoms     *upatoms;
  
} t_protonate;

typedef struct {
  char *res1, *res2;
  char *atom1,*atom2;
  char *newres1,*newres2;
  int  nbond1,nbond2;
  real length;
} t_specbond;

extern t_specbond *get_specbonds(int *nspecbond);
extern void done_specbonds(int nsb,t_specbond sb[]);

extern void free_t_restp(int nrtp, t_restp **rtp);
extern void free_t_hack(int nh, t_hack **h);
extern void free_t_hackblock(int nhb, t_hackblock **hb);
/* free the whole datastructure */

extern void clear_t_hackblock(t_hackblock *hb);
extern void clear_t_hack(t_hack *hack);
/* reset struct */

extern void merge_t_bondeds(t_rbondeds s[], t_rbondeds d[], 
			    bool bMin, bool bPlus);
/* add s[].b[] to d[].b[]
 * If bMin==TRUE, don't copy bondeds with atoms starting with '-'
 * If bPlus==TRUE, don't copy bondeds with atoms starting with '+'
 */
     
extern void copy_t_restp(t_restp *s, t_restp *d);
extern void copy_t_hack(t_hack *s, t_hack *d);
extern void copy_t_hackblock(t_hackblock *s, t_hackblock *d);
/* make copy of whole datastructure */

extern void merge_hacks_lo(int ns, t_hack *s, int *nd, t_hack **d);
/* add s[] to *d[] */

extern void merge_hacks(t_hackblock *s, t_hackblock *d);
/* add s->hacks[] to d->hacks[] */

extern void merge_t_hackblock(t_hackblock *s, t_hackblock *d);
/* add s->hacks[] and s->rb[] to d*/

extern void dump_hb(FILE *out, int nres, t_hackblock hb[]);
/* print out whole datastructure */

extern void init_t_protonate(t_protonate *protonate);
/* initialize t_protein struct */

#endif	/* _hackblock_h */
