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
 * Gyas ROwers Mature At Cryogenic Speed
 */

#ifndef _hackblock_h
#define _hackblock_h

static char *SRCID_hackblock_h = "$Id$";

#include "typedefs.h"
#include "pdbio.h"

/* Used for reading .rtp/.tdb */
/* ebtsBONDS must be the first, new types can be added to the end */
/* these *MUST* correspond to the arrays in hackblock.c */
enum { ebtsBONDS, ebtsANGLES, ebtsPDIHS, ebtsIDIHS, ebtsNR };
extern char *btsNames[ebtsNR];
extern int btsNiatoms[ebtsNR];

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
  int     tp;       /* Type of attachment (1..10) */
  char 	  *a[4];    /* Control atoms i,j,k,l	  */
  rvec    newx;     /* calculated new position    */
} t_hack;

typedef struct {
  char      *name;  /* Name of hack block (residue or terminus) */
  int       nhack;  /* Number of atoms to hack                  */
  int       maxhack;/* used for efficient srenew-ing            */
  t_hack    *hack;  /* Hack list                                */
  /* list of bonded interactions to add */
  t_rbondeds rb[ebtsNR];
} t_hackblock;

extern void free_t_restp(int nrtp, t_restp **rtp);
extern void free_t_hack(int nh, t_hack **h);
extern void free_t_hackblock(int nhb, t_hackblock **hb);
/* free the whole datastructure */

extern void clear_t_hackblock(t_hackblock *hb);
extern void clear_t_hack(t_hack *hack);
/* reset struct */

extern void merge_t_bondeds(t_rbondeds s[], t_rbondeds d[]);
/* add s[].b[] to d[].b[] */
     
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

#endif	/* _hackblock_h */
