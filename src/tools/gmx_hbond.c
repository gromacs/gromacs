/*
 * $Id$
 * 
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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>

#include "statutil.h"
#include "copyrite.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "futil.h"
#include "tpxio.h"
#include "physics.h"
#include "macros.h"
#include "gmx_fatal.h"
#include "index.h"
#include "smalloc.h"
#include "vec.h"
#include "xvgr.h"
#include "gstat.h"
#include "matio.h"
#include "string2.h"
#include "pbc.h"
#include "correl.h"

#define max_hx 7
typedef int t_hx[max_hx];
#define NRHXTYPES max_hx
char *hxtypenames[NRHXTYPES]=
{"n-n","n-n+1","n-n+2","n-n+3","n-n+4","n-n+5","n-n>6"};
#define MAXHH 4

enum { gr0,  gr1,    grI,  grNR };
enum { hbNo, hbDist, hbHB, hbNR, hbR2}; 
  
static char *grpnames[grNR] = {"0","1","I" };

static bool bDebug = FALSE;

#define HB_NO 0
#define HB_YES 1<<0
#define HB_INS 1<<1
#define HB_YESINS HB_YES|HB_INS
#define HB_NR (1<<2)
#define MAXHYDRO 4

#define ISHB(h)   (((h) & 2) == 2)
#define ISDIST(h) (((h) & 1) == 1)
#define ISDIST2(h) (((h) & 4) == 4)


typedef struct {
  int nr;
  int maxnr;
  atom_id *atoms;
} t_ncell;

typedef struct {
  t_ncell d[grNR];
  t_ncell a[grNR];
} t_gridcell;

typedef int     t_icell[grNR];
typedef atom_id h_id[MAXHYDRO];
 
typedef struct {
  int      history[MAXHYDRO]; 
  /* Has this hbond existed ever? If so as hbDist or hbHB or both.
   * Result is stored as a bitmap (1 = hbDist) || (2 = hbHB)
   */
  /* Bitmask array which tells whether a hbond is present
   * at a given time. Either of these may be NULL 
   */
  int      n0;       /* First frame a HB was found     */ 
  int      nframes,maxframes;  /* Amount of frames in this hbond */
  unsigned int **h; 
  unsigned int **g; 
  /* See Xu and Berne, JPCB 105 (2001), p. 11929. We define the
   * function g(t) = [1-h(t)] H(t) where H(t) is one when the donor-
   * acceptor distance is less than the user-specified distance (typically
   * 0.35 nm).
   */
} t_hbond;

typedef struct {
  int     nra,max_nra;
  atom_id *acc;             /* Atom numbers of the acceptors     */
  int     *grp;             /* Group index                       */
  int     *aptr;            /* Map atom number to acceptor index */
} t_acceptors;

typedef struct {
  int      nrd,max_nrd;
  int      *don;               /* Atom numbers of the donors         */
  int      *grp;               /* Group index                        */
  int      *dptr;              /* Map atom number to donor index     */
  int      *nhydro;            /* Number of hydrogens for each donor */
  h_id     *hydro;             /* The atom numbers of the hydrogens  */
  h_id     *nhbonds;           /* The number of HBs per H at current */
} t_donors;

typedef struct {
  bool        bHBmap,bDAnr;
  int         wordlen;
  /* The following arrays are nframes long */
  int         nframes,max_frames,maxhydro;
  int         *nhb,*ndist;
  h_id        *n_bound;
  real        *time;
  t_icell     *danr;
  t_hx        *nhx;
  /* These structures are initialized from the topology at start up */
  t_donors    d;
  t_acceptors a;
  /* This holds a matrix with all possible hydrogen bonds */
  int         nrhb,nrdist;
  t_hbond     ***hbmap;
} t_hbdata;


/* Changed argument 'bMerge' into 'oneHB' below,
 * since -contact should cause maxhydro to be 1,
 * not just -merge.
 * - Erik Marklund May 29, 2006
 */
static t_hbdata *mk_hbdata(bool bHBmap,bool bDAnr,bool oneHB)
{
  t_hbdata *hb;
  
  snew(hb,1);
  hb->wordlen = 8*sizeof(unsigned int);
  hb->bHBmap  = bHBmap;
  hb->bDAnr   = bDAnr;
  if (oneHB)
    hb->maxhydro = 1;
  else
    hb->maxhydro = MAXHYDRO;
  
  return hb;
}

static void mk_hbmap(t_hbdata *hb,bool bTwo,bool bInsert)
{
  int  i,j;

  snew(hb->hbmap,hb->d.nrd);
  for(i=0; (i<hb->d.nrd); i++) {
    snew(hb->hbmap[i],hb->a.nra);
    if (hb->hbmap[i] == NULL)
      gmx_fatal(FARGS,"Could not allocate enough memory for hbmap");
  }
}

static void add_frames(t_hbdata *hb,int nframes)
{
  int  i,j,k,l;
  
  if (nframes >= hb->max_frames) {
    hb->max_frames += 4096;
    srenew(hb->time,hb->max_frames);
    srenew(hb->nhb,hb->max_frames);
    srenew(hb->ndist,hb->max_frames);
    srenew(hb->n_bound,hb->max_frames);
    srenew(hb->nhx,hb->max_frames);
    if (hb->bDAnr)
      srenew(hb->danr,hb->max_frames);
  }
  hb->nframes=nframes;
}

#define OFFSET(frame) (frame / 32)
#define MASK(frame)   (1 << (frame % 32))

static void _set_hb(unsigned int hbexist[],unsigned int frame,bool bValue)
{
  if (bValue)
    hbexist[OFFSET(frame)] |= MASK(frame);
  else
    hbexist[OFFSET(frame)] &= ~MASK(frame);
}

static bool is_hb(unsigned int hbexist[],int frame)
{
  return ((hbexist[OFFSET(frame)] & MASK(frame)) != 0) ? 1 : 0;
}

static void set_hb(t_hbdata *hb,int id,int ih, int ia,int frame,int ihb)
{
  unsigned int *ghptr=NULL;
  
  if (ihb == hbHB)
    ghptr = hb->hbmap[id][ia]->h[ih];
  else if (ihb == hbDist)
    ghptr = hb->hbmap[id][ia]->g[ih];
  else
    gmx_fatal(FARGS,"Incomprehensible iValue %d in set_hb",ihb);

  _set_hb(ghptr,frame-hb->hbmap[id][ia]->n0,TRUE);
}

static void add_ff(t_hbdata *hbd,int id,int h,int ia,int frame,int ihb)
{
  int     i,j,n;
  t_hbond *hb      = hbd->hbmap[id][ia];
  int     maxhydro = hbd->d.nhydro[id];
  int     wlen     = hbd->wordlen;
  int     delta    = 32*wlen;
  
  if (!hb->h[0]) {
    hb->n0        = frame;
    hb->maxframes = delta;
    for(i=0; (i<maxhydro); i++) {
      snew(hb->h[i],hb->maxframes/wlen);
      snew(hb->g[i],hb->maxframes/wlen);
    }
  }
  else {
    hb->nframes = frame-hb->n0;
    /* We need a while loop here because hbonds may be returning
     * after a long time.
     */
    while (hb->nframes >= hb->maxframes) {
      n = hb->maxframes + delta;
      for(i=0; (i<maxhydro); i++) {
	srenew(hb->h[i],n/wlen);
	srenew(hb->g[i],n/wlen);
	for(j=hb->maxframes/wlen; (j<n/wlen); j++) {
	  hb->h[i][j] = 0;
	  hb->g[i][j] = 0;
	}
      }
      hb->maxframes = n;
    }
  }
  if (frame >= 0)
    set_hb(hbd,id,h,ia,frame,ihb);
  /*hb->nframes++;*/
}

static void inc_nhbonds(t_donors *ddd,int d, int h)
{
  int j;
  int dptr = ddd->dptr[d];
  
  for(j=0; (j<ddd->nhydro[dptr]); j++)
    if (ddd->hydro[dptr][j] == h) {
      ddd->nhbonds[dptr][j]++;
      break;
    }
  if (j == ddd->nhydro[dptr])
    gmx_fatal(FARGS,"No such hydrogen %d on donor %d\n",h+1,d+1);
}
/* Added argument bContact. The reason may not be obvious,
   but when using -contacts all contacts are stored in h[],
   while contacts within -r2 (when provided) are stored in g[].
   Therefore bContact needs to be passed in order to prevent
   add_hbond from trying to deal with hydrogens.
 * - Erik Marklund, June 29, 2006 
 */
static void add_hbond(t_hbdata *hb,int d,int a,int h,int grpd,int grpa,
		      int frame,bool bInsert,bool bMerge,int ihb,bool bContact)
{ 
  int k,id,ia,hh;
  
  if ((id = hb->d.dptr[d]) == NOTSET)
    gmx_fatal(FARGS,"No donor atom %d",d+1);
  else if (grpd != hb->d.grp[id])
    gmx_fatal(FARGS,"Inconsistent donor groups, %d iso %d, atom %d",
	      grpd,hb->d.grp[id],d+1);
  if ((ia = hb->a.aptr[a]) == NOTSET)
    gmx_fatal(FARGS,"No acceptor atom %d",a+1);
  else if (grpa != hb->a.grp[ia])
    gmx_fatal(FARGS,"Inconsistent acceptor groups, %d iso %d, atom %d",
	      grpa,hb->a.grp[ia],a+1);

  if (hb->hbmap) {
    /* Loop over hydrogens to find which hydrogen is in this particular HB */
    if ((ihb == hbHB) && !bMerge && !bContact) {
      for(k=0; (k<hb->d.nhydro[id]); k++) 
	if (hb->d.hydro[id][k] == h)
	  break;
      if (k == hb->d.nhydro[id])
	gmx_fatal(FARGS,"Donor %d does not have hydrogen %d (a = %d)",
		  d+1,h+1,a+1);
    }
    else
      k = 0;
    
    if (hb->bHBmap) {
      if (hb->hbmap[id][ia] == NULL) {
	snew(hb->hbmap[id][ia],1);
	snew(hb->hbmap[id][ia]->h,hb->maxhydro);
	snew(hb->hbmap[id][ia]->g,hb->maxhydro);
      }
      add_ff(hb,id,k,ia,frame,ihb);
    }
    
    /* Strange construction with frame >=0 is a relic from old code
     * for selected hbond analysis. It may be necessary again if that
     * is made to work again.
     */
    if (frame >= 0) {
      hh = hb->hbmap[id][ia]->history[k];
      if (ihb == hbHB) {
	hb->nhb[frame]++;
	if (!(ISHB(hh))) {
	  hb->hbmap[id][ia]->history[k] = hh | 2;
	  hb->nrhb++;
	}
      }
      else if (ihb == hbDist) {
	hb->ndist[frame]++;
	if (!(ISDIST(hh))) {
	  hb->hbmap[id][ia]->history[k] = hh | 1;
	  hb->nrdist++;
	}
      }
    }
  } else {
    if (frame >= 0) {
      if (ihb == hbHB) {
	hb->nhb[frame]++;
      } else if (ihb == hbDist) {
	hb->ndist[frame]++;
      }
    }
  }
  /* Increment number if HBonds per H */
  if (ihb == hbHB && !bContact)
    inc_nhbonds(&(hb->d),d,h);
}

static bool in_list(atom_id selection,int isize,atom_id *index)
{
  int i;
  bool bFound;
  
  bFound=FALSE;
  for(i=0; (i<isize) && !bFound; i++)
    if(selection == index[i])
      bFound=TRUE;
  
  return bFound;
}

static char *mkatomname(t_atoms *atoms,int i)
{
  static char buf[32];
  int rnr;
  
  rnr = atoms->atom[i].resnr;
  sprintf(buf,"%4s%d%-4s",*atoms->resname[rnr],rnr+1,*atoms->atomname[i]);
  
  return buf;
}

static void add_acc(t_acceptors *a,int ia,int grp)
{
  if (a->nra >= a->max_nra) {
    a->max_nra += 16;
    srenew(a->acc,a->max_nra);
    srenew(a->grp,a->max_nra);
  }
  a->grp[a->nra]   = grp;
  a->acc[a->nra++] = ia;
}

static void search_acceptors(t_topology *top,int isize, 
			     atom_id *index,t_acceptors *a,int grp,
			     bool bNitAcc,
			     bool bContact,bool bDoIt)
{
  int i;
  
  for (i=0; (i<top->atoms.nr); i++) {
    if (bDoIt) {
      if ((bContact ||
	   (((*top->atoms.atomname[i])[0] == 'O') || 
	    (bNitAcc && ((*top->atoms.atomname[i])[0] == 'N')))) &&
	  in_list(i,isize,index)) {
	add_acc(a,i,grp);
      }
    }
  }
  snew(a->aptr,top->atoms.nr);
  for(i=0; (i<top->atoms.nr); i++)
    a->aptr[i] = NOTSET;
  for(i=0; (i<a->nra); i++)
    a->aptr[a->acc[i]] = i;
}

static int _acceptor_index(t_acceptors *a,int grp,atom_id i,
			   char *file,int line)
{
  int ai = a->aptr[i];

  if (a->grp[ai] != grp) {
    if (debug && bDebug) 
      fprintf(debug,"Acc. group inconsist.. grp[%d] = %d, grp = %d (%s, %d)\n",
	      ai,a->grp[ai],grp,file,line);
    return NOTSET;
  }
  else
    return ai;
}
#define acceptor_index(a,grp,i) _acceptor_index(a,grp,i,__FILE__,__LINE__)

static int _donor_index(t_donors *d,int grp,atom_id i,char *file,int line)
{
  int di = d->dptr[i];
  
  if (d->grp[di] != grp) {
    if (debug && bDebug)
      fprintf(debug,"Don. group inconsist.. grp[%d] = %d, grp = %d (%s, %d)\n",
	      di,d->grp[di],grp,file,line);
    return NOTSET;
  }
  else
    return di;
}
#define donor_index(d,grp,i) _donor_index(d,grp,i,__FILE__,__LINE__)

static void add_h2d(int id,int ih,t_donors *ddd)
{
  int i;
  
  for(i=0; (i<ddd->nhydro[id]); i++) 
    if (ddd->hydro[id][i] == ih) {
      printf("Hm. This isn't first time I find this donor (%d,%d)\n",
	     ddd->don[id],ih);
      break;
    }
  if (i == ddd->nhydro[id]) {
    if (ddd->nhydro[id] >= MAXHYDRO)
      gmx_fatal(FARGS,"Donor %d has more than %d hydrogens!",
		  ddd->don[id],MAXHYDRO);
    ddd->hydro[id][i] = ih;
    ddd->nhydro[id]++;
  }
}
  
static void add_dh(t_donors *ddd,int id,int ih,int grp)
{
  int i;
  
  for(i=0; (i<ddd->nrd); i++) 
    if (ddd->don[i] == id) {
      add_h2d(i,ih,ddd);
      break;
    }
  if (i == ddd->nrd) {
    if (ddd->nrd >= ddd->max_nrd) {
      ddd->max_nrd += 128;
      srenew(ddd->don,ddd->max_nrd);
      srenew(ddd->nhydro,ddd->max_nrd);
      srenew(ddd->hydro,ddd->max_nrd);
      srenew(ddd->nhbonds,ddd->max_nrd);
      srenew(ddd->grp,ddd->max_nrd);
    }
    ddd->don[ddd->nrd] = id;
    ddd->nhydro[ddd->nrd] = 0;
    ddd->grp[ddd->nrd] = grp;
    add_h2d(ddd->nrd,ih,ddd);
    ddd->nrd++;
  }
}

static void search_donors(t_topology *top, int isize, atom_id *index,
			  t_donors *ddd,int grp,bool bContact,bool bDoIt)
{
  int        i,j,nra;
  t_functype func_type;
  t_ilist    *interaction;
  atom_id    nr1,nr2;
  bool       stop;
  
  if (bContact) {
    if (bDoIt)
      for(i=0; (i<isize); i++) 
	add_dh(ddd,index[i],-1,grp);
  }
  else {
    for(func_type=0; (func_type < F_NRE); func_type++) {
      interaction=&(top->idef.il[func_type]);
      for(i=0; i < interaction->nr; 
	  i+=interaction_function[top->idef.functype[interaction->iatoms[i]]].nratoms+1) {
	/* next function */
	if (func_type != top->idef.functype[interaction->iatoms[i]]) {
	  gmx_fatal(FARGS,"Error in func_type %s",
		      interaction_function[func_type].longname);
	}
	
	/* check out this functype */
	if (func_type == F_SETTLE) {
	  nr1=interaction->iatoms[i+1];
	  
	  if (in_list(nr1,  isize,index)) {
	    if (in_list(nr1+1,isize,index))
	      add_dh(ddd,nr1,nr1+1,grp);
	    if (in_list(nr1+2,isize,index))
	      add_dh(ddd,nr1,nr1+2,grp);
	  }
	} 
	else if (IS_CHEMBOND(func_type)) {
	  for (j=0; j<2; j++) {
	    nr1=interaction->iatoms[i+1+j];
	    nr2=interaction->iatoms[i+2-j];
	    if ((*top->atoms.atomname[nr1][0] == 'H') && 
		((*top->atoms.atomname[nr2][0] == 'O') ||
		 (*top->atoms.atomname[nr2][0] == 'N')) &&
		in_list(nr1,isize,index) && in_list(nr2,isize,index))
	      add_dh(ddd,nr2,nr1,grp);
	  }
	}
      }
    }
    for(func_type=0; func_type < F_NRE; func_type++) {
      interaction=&top->idef.il[func_type];
      for(i=0; i < interaction->nr; 
	  i+=interaction_function[top->idef.functype[interaction->iatoms[i]]].nratoms+1) {
	/* next function */
	if (func_type != top->idef.functype[interaction->iatoms[i]])
	  gmx_incons("function type in search_donors");
	
	if ( interaction_function[func_type].flags & IF_VSITE ) {
	  nr1=interaction->iatoms[i+1];
	  if ( *top->atoms.atomname[nr1][0]  == 'H') {
	    nr2=nr1-1;
	    stop=FALSE;
	    while (!stop && ( *top->atoms.atomname[nr2][0] == 'H'))
	      if (nr2)
		nr2--;
	      else
		stop=TRUE;
	    if ( !stop && ( ( *top->atoms.atomname[nr2][0] == 'O') ||
			    ( *top->atoms.atomname[nr2][0] == 'N') ) &&
		 in_list(nr1,isize,index) && in_list(nr2,isize,index) )
	      add_dh(ddd,nr2,nr1,grp);
	  }
	}
      }
    }
  }
  snew(ddd->dptr,top->atoms.nr);
  for(i=0; (i<top->atoms.nr); i++)
    ddd->dptr[i] = NOTSET;
  for(i=0; (i<ddd->nrd); i++)
    ddd->dptr[ddd->don[i]] = i;
}

static t_gridcell ***init_grid(bool bBox,rvec box[],real rcut,ivec ngrid)
{
  t_gridcell ***grid;
  int i,y,z;
  
  if (bBox)
    for(i=0; i<DIM; i++)
      ngrid[i]=(box[i][i]/(1.2*rcut));
  
  if ( !bBox || (ngrid[XX]<3) || (ngrid[YY]<3) || (ngrid[ZZ]<3) )
    for(i=0; i<DIM; i++)
      ngrid[i]=1;
  else 
    printf("\nWill do grid-seach on %dx%dx%d grid, rcut=%g\n",
	   ngrid[XX],ngrid[YY],ngrid[ZZ],rcut);
  snew(grid,ngrid[ZZ]);
  for (z=0; z<ngrid[ZZ]; z++) {
    snew((grid)[z],ngrid[YY]);
    for (y=0; y<ngrid[YY]; y++)
      snew((grid)[z][y],ngrid[XX]);
  }
  return grid;
}

static void reset_nhbonds(t_donors *ddd)
{
  int i,j;
  
  for(i=0; (i<ddd->nrd); i++) 
    for(j=0; (j<MAXHH); j++)
      ddd->nhbonds[i][j] = 0;
}

static void build_grid(t_hbdata *hb,rvec x[], rvec xshell,
		       bool bBox, matrix box, rvec hbox,
		       real rcut, real rshell,
		       ivec ngrid, t_gridcell ***grid)
{
  int     i,m,gr,xi,yi,zi,nr;
  atom_id *ad;
  ivec    grididx;
  rvec    invdelta,dshell;
  t_ncell *newgrid;
  bool    bDoRshell,bInShell,bAcc;
  real    rshell2=0;
  int     gx,gy,gz;
  int     dum = -1;
  
  bDoRshell = (rshell > 0);
  rshell2   = sqr(rshell);
  bInShell  = TRUE;
  
#define DBB(x) if (debug && bDebug) fprintf(debug,"build_grid, line %d, %s = %d\n",__LINE__,#x,x)
  DBB(dum);
  for(m=0; m<DIM; m++) {
    hbox[m]=box[m][m]*0.5;
    if (bBox) {
      invdelta[m]=ngrid[m]/box[m][m];
      if (1/invdelta[m] < rcut)
	gmx_fatal(FARGS,"Your computational box has shrunk too much.\n"
		  "%s can not handle this situation, sorry.\n",
		  ShortProgram());
    } else
      invdelta[m]=0;
  }
  grididx[XX]=0;
  grididx[YY]=0;
  grididx[ZZ]=0;
  DBB(dum);
  /* resetting atom counts */
  for(gr=0; (gr<grNR); gr++) {
    for (zi=0; zi<ngrid[ZZ]; zi++)
      for (yi=0; yi<ngrid[YY]; yi++)
	for (xi=0; xi<ngrid[XX]; xi++) {
	  grid[zi][yi][xi].d[gr].nr=0;
	  grid[zi][yi][xi].a[gr].nr=0;
      }
    DBB(dum);
    
    /* put atoms in grid cells */
    for(bAcc=FALSE; (bAcc<=TRUE); bAcc++) {
      if (bAcc) {
	nr = hb->a.nra;
	ad = hb->a.acc;
      }
      else {
	nr = hb->d.nrd;
	ad = hb->d.don;
      }
      DBB(bAcc);
      for(i=0; (i<nr); i++) {
	/* check if we are inside the shell */
	/* if bDoRshell=FALSE then bInShell=TRUE always */
	DBB(i);
	if ( bDoRshell ) {
	  bInShell=TRUE;
	  rvec_sub(x[ad[i]],xshell,dshell);
	  if (bBox) 
	    for(m=DIM-1; m>=0 && bInShell; m--) {
	      if ( dshell[m] < -hbox[m] )
		rvec_inc(dshell,box[m]);
	      else if ( dshell[m] >= hbox[m] ) 
		dshell[m] -= 2*hbox[m];
	      /* if we're outside the cube, we're outside the sphere also! */
	      if ( (dshell[m]>rshell) || (-dshell[m]>rshell) )
		bInShell=FALSE;
	    }
	  /* if we're inside the cube, check if we're inside the sphere */
	  if (bInShell)
	    bInShell = norm2(dshell) < rshell2;
	}
	DBB(i);
	if ( bInShell ) {
	  if (bBox) 
	    for(m=DIM-1; m>=0; m--) {
	      /* put atom in the box */
	      while( x[ad[i]][m] < 0 ) 
		rvec_inc(x[ad[i]],box[m]);
	      while( x[ad[i]][m] >= box[m][m] ) 
		rvec_dec(x[ad[i]],box[m]);
	      /* determine grid index of atom */
	      grididx[m]=x[ad[i]][m]*invdelta[m];
	      grididx[m] = (grididx[m]+ngrid[m]) % ngrid[m];
	    }
	  gx = grididx[XX];
	  gy = grididx[YY];
	  gz = grididx[ZZ];
	  range_check(gx,0,ngrid[XX]);
	  range_check(gy,0,ngrid[YY]);
	  range_check(gz,0,ngrid[ZZ]);
	  DBB(gx);
	  DBB(gy);
	  DBB(gz);
	  /* add atom to grid cell */
	  if (bAcc)
	    newgrid=&(grid[gz][gy][gx].a[gr]);
	  else
	    newgrid=&(grid[gz][gy][gx].d[gr]);
	  if (newgrid->nr >= newgrid->maxnr) {
	    newgrid->maxnr+=10;
	    DBB(newgrid->maxnr);
	    srenew(newgrid->atoms, newgrid->maxnr);
	  }
	  DBB(newgrid->nr);
	  newgrid->atoms[newgrid->nr]=ad[i];
	  newgrid->nr++;
	}
      }
    }
  }
}

static void count_da_grid(ivec ngrid, t_gridcell ***grid, t_icell danr)
{
  int gr,xi,yi,zi;
  
  for(gr=0; (gr<grNR); gr++) {
    danr[gr]=0;
    for (zi=0; zi<ngrid[ZZ]; zi++)
      for (yi=0; yi<ngrid[YY]; yi++)
	for (xi=0; xi<ngrid[XX]; xi++) {
	  danr[gr] += grid[zi][yi][xi].d[gr].nr;
	}
  }
}

/* The grid loop.
 * Without a box, the grid is 1x1x1, so all loops are 1 long.
 * With a rectangular box (bTric==FALSE) all loops are 3 long.
 * With a triclinic box all loops are 3 long, except when a cell is
 * located next to one of the box edges which is not parallel to the
 * x/y-plane, in that case all cells in a line or layer are searched.
 * This could be implemented slightly more efficient, but the code
 * would get much more complicated.
 */
#define B(n,x,bTric,bEdge) ((n==1) ? x : bTric&&(bEdge) ? 0   : (x-1))
#define E(n,x,bTric,bEdge) ((n==1) ? x : bTric&&(bEdge) ? n-1 : (x+1))
#define GRIDMOD(j,n) (j+n)%(n)
#define LOOPGRIDINNER(x,y,z,xx,yy,zz,xo,yo,zo,n,bTric)\
     for(zz=B(n[ZZ],zo,bTric,FALSE); zz<=E(n[ZZ],zo,bTric,FALSE); zz++) {\
       z=GRIDMOD(zz,n[ZZ]);\
       for(yy=B(n[YY],yo,bTric,z==0||z==n[ZZ]-1); \
	  yy<=E(n[YY],yo,bTric,z==0||z==n[ZZ]-1); yy++) {\
	 y=GRIDMOD(yy,n[YY]);\
	 for(xx=B(n[XX],xo,bTric,y==0||y==n[YY]-1||z==0||z==n[ZZ]-1); \
	     xx<=E(n[XX],xo,bTric,y==0||y==n[YY]-1||z==0||z==n[ZZ]-1); xx++) {\
	   x=GRIDMOD(xx,n[XX]);
#define ENDLOOPGRIDINNER\
	 }\
       }\
     }\

static void dump_grid(FILE *fp, ivec ngrid, t_gridcell ***grid)
{
  int gr,x,y,z,sum[grNR];
  
  fprintf(fp,"grid %dx%dx%d\n",ngrid[XX],ngrid[YY],ngrid[ZZ]);
  for (gr=0; gr<grNR; gr++) {
    sum[gr]=0;
    fprintf(fp,"GROUP %d (%s)\n",gr,grpnames[gr]);
    for (z=0; z<ngrid[ZZ]; z+=2) {
      fprintf(fp,"Z=%d,%d\n",z,z+1);
      for (y=0; y<ngrid[YY]; y++) {
	for (x=0; x<ngrid[XX]; x++) {
	  fprintf(fp,"%3d",grid[x][y][z].d[gr].nr);
	  sum[gr]+=grid[z][y][x].d[gr].nr;
	  fprintf(fp,"%3d",grid[x][y][z].a[gr].nr);
	  sum[gr]+=grid[z][y][x].a[gr].nr;
	  
	}
	fprintf(fp," | ");
	if ( (z+1) < ngrid[ZZ] )
	  for (x=0; x<ngrid[XX]; x++) {
	    fprintf(fp,"%3d",grid[z+1][y][x].d[gr].nr);
	    sum[gr]+=grid[z+1][y][x].d[gr].nr;
	    fprintf(fp,"%3d",grid[z+1][y][x].a[gr].nr);
	    sum[gr]+=grid[z+1][y][x].a[gr].nr;
	  }
	fprintf(fp,"\n");
      }
    }
  }
  fprintf(fp,"TOTALS:");
  for (gr=0; gr<grNR; gr++)
    fprintf(fp," %d=%d",gr,sum[gr]);
  fprintf(fp,"\n");
}

/* New GMX record! 5 * in a row. Congratulations! 
 * Sorry, only four left.
 */
static void free_grid(ivec ngrid, t_gridcell ****grid)
{
  int y,z;
  t_gridcell ***g = *grid;
  
  for (z=0; z<ngrid[ZZ]; z++) {
    for (y=0; y<ngrid[YY]; y++) {
      sfree(g[z][y]);
    }
    sfree(g[z]);
  }
  sfree(g);
  g=NULL;
}

static void pbc_correct(rvec dx,matrix box,rvec hbox)
{
  int m;
  
  for(m=DIM-1; m>=0; m--) {
    if ( dx[m] < -hbox[m] )
      rvec_inc(dx,box[m]);
    else if ( dx[m] >= hbox[m] )
      rvec_dec(dx,box[m]);
  }
}

static int low_is_hbond(int d,int a,int h,
			rvec r_da,real rcut2,real ccut, 
			rvec x[], bool bBox, matrix box,rvec hbox,
			real *d_ha, real *ang,real d2)
{
  rvec r_dh,r_ha;
  real ca;
  
  if (d2 == 0) {
    rvec_sub(x[h],x[a],r_ha);
    if (bBox) 
      pbc_correct(r_ha,box,hbox);    
    d2 = iprod(r_ha,r_ha);
  }
  
  if ( d2 <= rcut2 ) {
    rvec_sub(x[d],x[h],r_dh);
    if (bBox)
      pbc_correct(r_dh,box,hbox);
    
    ca = cos_angle(r_dh,r_da);
    /* if angle is smaller, cos is larger */
    if (ca >= ccut) {
      *d_ha = sqrt(d2);
      *ang = acos(ca);
      return hbHB;
    }
    else
      return hbDist;
  }
  return hbNo;
}

/* Added argument r2cut, changed contact and implemented 
 * use of second cut-off.
 * - Erik Marklund, June 29, 2006
 */
static int is_hbond(t_hbdata *hb,int grpd,int grpa,int d,int a,
		    real rcut, real r2cut, real ccut, 
		    rvec x[], bool bBox, matrix box,rvec hbox,
		    real *d_ha, real *ang,bool bDA,int *hhh,
		    bool bContact)
{
  int  h,hh,id,ja,ihb;
  rvec r_da;
  real rc2,r2c2,d2;
  
  if (d == a)
    return hbNo;

  if (((id = donor_index(&hb->d,grpd,d)) == NOTSET) ||
      ((ja = acceptor_index(&hb->a,grpa,a)) == NOTSET))
    return hbNo;
  
  rvec_sub(x[d],x[a],r_da);
  if (bBox) 
    pbc_correct(r_da,box,hbox);    
  d2 = iprod(r_da,r_da);
  rc2 = rcut*rcut;
  r2c2 = r2cut*r2cut;  
  *hhh = NOTSET;
  if (d2 < rc2) {
    if (bContact) 
      return hbHB; /* Ok, it's not a hb, but it should be stored in h[].*/
      /* return hbDist; */
    
    if (!bDA)
      d2 = 0;
    
    for(h=0; (h < hb->d.nhydro[id]); h++) {
      hh = hb->d.hydro[id][h];
      ihb = low_is_hbond(d,a,hh,r_da,rc2,ccut,x,bBox,box,hbox,d_ha,ang,d2);
      if (ihb == hbHB) {
	*hhh = hh;
	return ihb;
      }
    }
    return hbDist;
  }
  if (d2 < r2c2)
    return hbDist;
  return hbNo;
}

/* Fixed previously undiscovered bug in the merge
   code, where the last frame of each hbond disappears.
   - Erik Marklund, June 1, 2006 */
static void do_merge(t_hbdata *hb,int ntmp,
		     unsigned int htmp[],unsigned int gtmp[],
		     t_hbond *hb0,t_hbond *hb1)
{
  int m,mm,n00,n01,nn0,nnframes;
  /* Decide where to start from when merging */
  n00      = hb0->n0;
  n01      = hb1->n0;
  nn0      = min(n00,n01);
  nnframes = max(n00 + hb0->nframes,n01 + hb1->nframes) - nn0;
  /* Initiate tmp arrays */
  for(m=0; (m<ntmp); m++) {
    htmp[m] = 0;
    gtmp[m] = 0;
  }
  /* Fill tmp arrays with values due to first HB */
  /* Once again '<' had to be replaced with '<='
     to catch the last frame in which the hbond
     appears.
     - Erik Marklund, June 1, 2006 */  
  for(m=0; (m<=hb0->nframes); m++) {
    mm = m+n00-nn0;
    htmp[mm] = is_hb(hb0->h[0],m);
    gtmp[mm] = is_hb(hb0->g[0],m);
  }
  /* Next HB */
  for(m=0; (m<=hb1->nframes); m++) {
    mm = m+n01-nn0;
    htmp[mm] = htmp[mm] || is_hb(hb1->h[0],m);
    gtmp[mm] = gtmp[mm] || is_hb(hb1->g[0],m);
  }
  /* Reallocate target array */
  if (nnframes > hb0->maxframes) {
    srenew(hb0->h[0],4+nnframes/hb->wordlen);
    srenew(hb0->g[0],4+nnframes/hb->wordlen);
  }
  /* Copy temp array to target array */
  for(m=0; (m<=nnframes); m++) {
    _set_hb(hb0->h[0],m,htmp[m]);
    _set_hb(hb0->g[0],m,gtmp[m]);
  }
  /* Set scalar variables */
  hb0->n0        = nn0;
  hb0->maxframes = nnframes;
}

/* Added argument bContact for nicer output.
 * Erik Marklund, June 29, 2006
 */
static void merge_hb(t_hbdata *hb,bool bTwo, bool bContact)
{
  int  i,inrnew,indnew,j,ii,jj,m,id,ia,grp,ogrp,ntmp;
  unsigned int *htmp,*gtmp;
  t_hbond *hb0,*hb1;

  inrnew = hb->nrhb;
  indnew = hb->nrdist;
    
  /* Check whether donors are also acceptors */
  printf("Merging hbonds with Acceptor and Donor swapped\n");

  ntmp = 2*hb->max_frames;
  snew(gtmp,ntmp);
  snew(htmp,ntmp);
  for(i=0; (i<hb->d.nrd); i++) {
    fprintf(stderr,"\r%d/%d",i+1,hb->d.nrd);
    id = hb->d.don[i];
    ii = hb->a.aptr[id];
    for(j=0; (j<hb->a.nra); j++) {
      ia = hb->a.acc[j];
      jj = hb->d.dptr[ia];
      if ((id != ia) && (ii != NOTSET) && (jj != NOTSET) &&
	  (!bTwo || (bTwo && (hb->d.grp[id] != hb->a.grp[ia])))) {
	hb0 = hb->hbmap[i][j];
	hb1 = hb->hbmap[jj][ii];
	if (hb0 && hb1 && ISHB(hb0->history[0]) && ISHB(hb1->history[0])) {
	  do_merge(hb,ntmp,htmp,gtmp,hb0,hb1);
	  if (ISHB(hb1->history[0])) 
	    inrnew--;
	  else if (ISDIST(hb1->history[0])) 
	    indnew--;
	  else
	    if (bContact) 
	      gmx_incons("No contact history");
	    else
	      gmx_incons("Neither hydrogen bond nor distance");
	  sfree(hb1->h[0]);
	  sfree(hb1->g[0]);
	  hb1->h[0] = NULL;
	  hb1->g[0] = NULL;
	  hb1->history[0] = hbNo;
	}
      }
    }
  }
  fprintf(stderr,"\n");
  printf("- Reduced number of hbonds from %d to %d\n",hb->nrhb,inrnew);
  printf("- Reduced number of distances from %d to %d\n",hb->nrdist,indnew);
  hb->nrhb   = inrnew;
  hb->nrdist = indnew;
  sfree(gtmp);
  sfree(htmp);
}

static void do_nhb_dist(FILE *fp,t_hbdata *hb,real t) 
{
  int  i,j,k,n_bound[MAXHH],nbtot;
  h_id nhb;

  
  /* Set array to 0 */
  for(k=0; (k<MAXHH); k++)
    n_bound[k] = 0;
  /* Loop over possible donors */
  for(i=0; (i<hb->d.nrd); i++) {
    for(j=0; (j<hb->d.nhydro[i]); j++)
      n_bound[hb->d.nhbonds[i][j]]++;
  }      
  fprintf(fp,"%12.5e",t);
  nbtot = 0;
  for(k=0; (k<MAXHH); k++) {
    fprintf(fp,"  %8d",n_bound[k]);
    nbtot += n_bound[k]*k;
  }
  fprintf(fp,"  %8d\n",nbtot);
}

/* Added argument bContact in do_hblife(...). Also
 * added support for -contact in function body.
 * - Erik Marklund, May 31, 2006 */
/* Changed the contact code slightly.
 * - Erik Marklund, June 29, 2006
 */
static void do_hblife(char *fn,t_hbdata *hb,bool bMerge,bool bContact)
{
  FILE *fp;
  static char *leg[] = { "p(t)", "t p(t)" };
  int  *histo;
  int  i,j,j0,k,m,nh,ihb,ohb,nhydro,ndump=0;
  int   nframes = hb->nframes;
  unsigned int **h;
  real   t,x1,dt;
  double sum,integral;
  t_hbond *hbh;
  
  snew(h,hb->maxhydro);
  snew(histo,nframes+1);
  /* Total number of hbonds analyzed here */
  for(i=0; (i<hb->d.nrd); i++) {
    for(k=0; (k<hb->a.nra); k++) {
      hbh = hb->hbmap[i][k];
      if (hbh) {
	if (bMerge) {
	  if (hbh->h[0]) {
	    h[0] = hbh->h[0];
	    nhydro = 1;
	  }
	  else
	    nhydro = 0;
	}
	else {
	  nhydro = 0;
	  for(m=0; (m<hb->maxhydro); m++)
	    if (hbh->h[m]) {
	      h[nhydro++] = hbh->h[m];
	    }
	}
	for(nh=0; (nh<nhydro); nh++) {
	  ohb = 0;
	  j0  = 0;

	  /* Changed '<' into '<=' below, just like I
	     did in the hbm-output-loop in the main code.
	     - Erik Marklund, May 31, 2006
	  */
	  for(j=0; (j<=hbh->nframes); j++) {
	    ihb      = is_hb(h[nh],j);
	    if (debug && (ndump < 10))
	      fprintf(debug,"%5d  %5d\n",j,ihb);
	    if (ihb != ohb) {
	      if (ihb) {
		j0 = j;
	      }
	      else {
		histo[j-j0]++;
	      }
	      ohb = ihb;
	    }
	  }
	  ndump++;
	}
      }
    }
  }
  fprintf(stderr,"\n");
  if (bContact)
    fp = xvgropen(fn,"Uninterrupted contact lifetime","Time (ps)","()");
  else
    fp = xvgropen(fn,"Uninterrupted hydrogen bond lifetime","Time (ps)","()");

  xvgr_legend(fp,asize(leg),leg);
  j0 = nframes-1;
  while ((j0 > 0) && (histo[j0] == 0))
    j0--;
  sum = 0;
  for(i=0; (i<=j0); i++)
    sum+=histo[i];
  dt       = hb->time[1]-hb->time[0];
  sum      = dt*sum;
  integral = 0;
  for(i=1; (i<=j0); i++) {
    t  = hb->time[i] - hb->time[0] - 0.5*dt;
    x1 = t*histo[i]/sum;
    fprintf(fp,"%8.3f  %10.3e  %10.3e\n",t,histo[i]/sum,x1);
    integral += x1;
  }
  integral *= dt;
  fclose(fp);
  printf("%s lifetime = %.2f ps\n", bContact?"Contact":"HB", integral);
  sfree(h);
  sfree(histo);
}

/* Changed argument bMerge into oneHB to handle contacts properly.
 * - Erik Marklund, June 29, 2006
 */
static void dump_ac(t_hbdata *hb,bool oneHB,int nDump)
{
  FILE  *fp;
  int   i,j,k,m,nd,ihb,idist;
  int   nframes = hb->nframes;
  bool  bPrint;
  t_hbond *hbh;

  if (nDump <= 0)
    return;
  fp = ffopen("debug-ac.xvg","w");
  for(j=0; (j<nframes); j++) {
    fprintf(fp,"%10.3f",hb->time[j]);
    for(i=nd=0; (i<hb->d.nrd) && (nd < nDump); i++) {
      for(k=0; (k<hb->a.nra) && (nd < nDump); k++) {
	bPrint = FALSE;
	ihb = idist = 0;
	hbh = hb->hbmap[i][k];
	if (oneHB) {
	  if (hbh->h[0]) {
	    ihb   = is_hb(hbh->h[0],j);
	    idist = is_hb(hbh->g[0],j);
	    bPrint = TRUE;
	  }
	} 
	else {
	  for(m=0; (m<hb->maxhydro) && !ihb ; m++) {
	    ihb   = ihb   || ((hbh->h[m]) && is_hb(hbh->h[m],j));
	    idist = idist || ((hbh->g[m]) && is_hb(hbh->g[m],j));
	  }
	  /* This is not correct! */
	  /* What isn't correct? -Erik M */
	  bPrint = TRUE;
	}
	if (bPrint) {
	  fprintf(fp,"  %1d-%1d",ihb,idist);
	  nd++;
	}
      }
    }
    fprintf(fp,"\n");
  }
  ffclose(fp);
}

static real calc_dg(real tau,real temp)
{
  real kbt;
  
  kbt = BOLTZ*temp;
  if (tau <= 0)
    return -666;
  else
    return kbt*log(kbt*tau/PLANCK);  
}

static void analyse_corr(int n,real t[],real ct[],real nt[],real kt[],
			 real fit_start,real temp)
{
  int    i;
  real   k=1,kp=1,dg,dgp,tau_hb,dtau,tau_rlx,e_1,dt;
  double tmp,sn2=0,sc2=0,sk2=0,scn=0,sck=0,snk=0;
  
  for(i=0; (i<n-2) && (t[i] < fit_start); i++)
    ;
  if (i < n-2) { 
    for(; (i<n); i++) {
      sc2 += sqr(ct[i]);
      sn2 += sqr(nt[i]);
      sk2 += sqr(kt[i]);
      sck += ct[i]*kt[i];
      snk += nt[i]*kt[i];
      scn += ct[i]*nt[i];
    }
    printf("   Hydrogen bond thermodynamics at T = %g K\n",temp);
    printf("--------------------------------------------------\n"
	   "Type      Rate (1/ps)  Time (ps)  DG (kJ/mol)\n");
    tmp = (sn2*sc2-sqr(scn));
    if ((tmp > 0) && (sn2 > 0)) {
      k    = (sn2*sck-scn*snk)/tmp;
      kp   = (k*scn-snk)/sn2;
      printf("Forward    %10.3f   %8.3f  %10.3f\n",k,1/k,calc_dg(1/k,temp));
      printf("Backward   %10.3f   %8.3f  %10.3f\n",kp,1/kp,calc_dg(1/kp,temp));
    }
    if (sc2 > 0) {
      k    = 2*sck/sc2;
      printf("One-way    %10.3f   %8.3f  %10.3f\n",k,1/k,calc_dg(1/k,temp));
    }
    else 
      printf(" - Numerical problems computing HB thermodynamics:\n"
	     "sc2 = %g  sn2 = %g  sk2 = %g sck = %g snk = %g scn = %g\n",
	     sc2,sn2,sk2,sck,snk,scn);
    /* Determine integral of the correlation function */
    tau_hb = evaluate_integral(n,t,ct,NULL,(t[n-1]-t[0])/2,&dtau);
    printf("Integral   %10.3f   %8.3f  %10.3f\n",1/tau_hb,tau_hb,
	   calc_dg(tau_hb,temp));
    e_1 = exp(-1);
    for(i=0; (i<n-2); i++) {
      if ((ct[i] > e_1) && (ct[i+1] <= e_1)) {
	break;
      }
    }
    if (i < n-2) {
      /* Determine tau_relax from linear interpolation */
      tau_rlx = t[i]-t[0] + (e_1-ct[i])*(t[i+1]-t[i])/(ct[i+1]-ct[i]);
      printf("Relaxation %10.3f   %8.3f  %10.3f\n",1/tau_rlx,tau_rlx,
	     calc_dg(tau_rlx,temp));
    }
  }
  else 
    printf("Correlation functions too short to compute thermodynamics\n");
}

static void compute_derivative(int nn,real x[],real y[],real dydx[])
{
  int j;
  
  /* Compute k(t) = -dc(t)/dt */
  for(j=1; (j<nn-1); j++)
    dydx[j] = (y[j+1]-y[j-1])/(x[j+1]-x[j-1]);
  /* Extrapolate endpoints */
  dydx[0]    = 2*dydx[1]   -  dydx[2];
  dydx[nn-1] = 2*dydx[nn-2] - dydx[nn-3];
}

/* Added argument bContact in do_hbac(...). Also
 * added support for -contact in the actual code.
 * - Erik Marklund, May 31, 2006 */
/* Changed contact code and added argument R2 
 * - Erik Marklund, June 29, 2006
 */
static void do_hbac(char *fn,t_hbdata *hb,real aver_nhb,real aver_dist,
		    int nDump,bool bMerge,bool bContact,real fit_start,real temp,bool R2)
{
  FILE *fp;
  static char *leg[] = { "Ac\\sfin sys\\v{}\\z{}(t)", "Ac(t)", "Cc\\scontact,hb\\v{}\\z{}(t)", "-dAc\\sfs\\v{}\\z{}/dt" };
  int  i,j,k,m,n,nd,ihb,idist,n2,nn;
  bool bNorm=FALSE;
  double nhb = 0;
  real *rhbex,*ht,*gt,*ght,*dght,*kt;
  real *ct,tail,tail2,dtail,ct_fac,ght_fac,*cct;
  const real tol = 1e-3;
  int   nframes = hb->nframes,nf;
  unsigned int **h,**g;
  int   nh,nhbonds,nhydro,ngh;
  t_hbond *hbh;
#define WITHOUT_FFTW
#ifndef WITHOUT_FFTW
  correl_t *correl;
#endif

  /* build hbexist matrix in reals for autocorr */
  /* Allocate memory for computing ACF (rhbex) and aggregating the ACF (ct) */
  n2=1;
  while (n2 < nframes)
    n2*=2;
  snew(rhbex,2*n2);
  snew(ct,2*n2);
  snew(gt,2*n2);
  snew(ht,2*n2);
  snew(ght,2*n2);
  snew(dght,2*n2);
  
  nn = nframes/2;
  snew(kt,nn);
  snew(cct,nn);
  
  snew(h,hb->maxhydro);
  snew(g,hb->maxhydro);

#ifndef WITHOUT_FFTW
  correl = init_correl(n2);
#endif
  /* Dump hbonds for debugging */
  dump_ac(hb,bMerge||bContact,nDump);
  
  /* Total number of hbonds analyzed here */
  nhbonds = 0;
  ngh     = 0;
  for(i=0; (i<hb->d.nrd); i++) {
    for(k=0; (k<hb->a.nra); k++) {
      nhydro = 0;
      hbh = hb->hbmap[i][k];
      if (hbh) {
	if (bMerge || bContact) {
	  if (ISHB(hbh->history[0])) {
	    h[0] = hbh->h[0];
	    g[0] = hbh->g[0];
	    nhydro = 1;
	  }
	}
	else {
	  for(m=0; (m<hb->maxhydro); m++)
	    if (ISHB(hbh->history[m])) {
	      g[nhydro] = hbh->g[m];
	      h[nhydro] = hbh->h[m];
	      nhydro++;
	    }
	}
	nf = hbh->nframes;
	for(nh=0; (nh<nhydro); nh++) {
	  if ((((nhbonds+1) % 10) == 0) || (nhbonds+1 == hb->nrhb))
	    fprintf(stderr,"\rACF %d/%d",nhbonds+1,hb->nrhb);
	  nhbonds++;
	  for(j=0; (j<nframes); j++) {
	    /* Changed '<' into '<=' below, just like I did in
	       the hbm-output-loop in the gmx_hbond() block.
	       - Erik Marklund, May 31, 2006 */
	    if (j <= nf) {
	      ihb   = is_hb(h[nh],j);
	      idist = is_hb(g[nh],j);
	    }
	    else {
	      ihb = idist = 0;
	    }
	    rhbex[j] = ihb-aver_nhb;
	    /* For contacts: if a second cut-off is provided, use it,
	     * otherwise use g(t) = 1-h(t) */
	    if (!R2 && bContact)
	      gt[j]  = 1-ihb;
	    else
	      gt[j]  = idist*(1-ihb); 
	    ht[j]    = rhbex[j];
	    nhb     += ihb;
	  }
	  
	  /* The autocorrelation function is normalized after summation only */
	  low_do_autocorr(NULL,NULL,nframes,1,-1,&rhbex,hb->time[1]-hb->time[0],
			  eacNormal,1,FALSE,bNorm,FALSE,0,-1,0,1);
	  
	  /* Cross correlation analysis for thermodynamics */
	  for(j=nframes; (j<n2); j++) {
	    ht[j] = 0;
	    gt[j] = 0;
	  }
	  /*correl_fftw(correl,ht,gt,dght);*/
	  cross_corr(n2,ht,gt,dght);
	  
	  for(j=0; (j<nn); j++) {
	    ct[j]  += rhbex[j];
	    ght[j] += dght[j];
	  }
	}
      }
    }
  }
  fprintf(stderr,"\n");
  
  /* Normalize */
  ct_fac  = 1.0/ct[0];
  ght_fac = 1.0/nhb;
  printf("Normalization for c(t) = %g for gh(t) = %g\n",ct_fac,ght_fac);
  for(j=0; (j<nn); j++) {
    ct[j]  *= ct_fac;
    ght[j] *= ght_fac; 
    /* Xu and Berne use the same normalization constant */
  }
  
  /* Determine tail value for statistics */
  tail  = 0;
  tail2 = 0;
  for(j=nn/2; (j<nn); j++) {
    tail  += ct[j];
    tail2 += ct[j]*ct[j];
  }
  tail  /= (nn - nn/2);
  tail2 /= (nn - nn/2);
  dtail  = sqrt(tail2-tail*tail);
  
  /* Check whether the ACF is long enough */
  if (dtail > tol) {
    printf("\nWARNING: Correlation function is probably not long enough\n"
	   "because the standard deviation in the tail of C(t) > %g\n"
	   "Tail value (average C(t) over second half of acf): %g +/- %g\n",
	   tol,tail,dtail);
  }
  for(j=0; (j<nn); j++) {
    cct[j] = ct[j];
    ct[j]  = (cct[j]-tail)/(1-tail); 
  }
  /* Compute negative derivative k(t) = -dc(t)/dt */
  compute_derivative(nn,hb->time,ct,kt);
  for(j=0; (j<nn); j++)
    kt[j] = -kt[j];

  if (bContact)
    fp = xvgropen(fn, "Contact Autocorrelation","Time (ps)","C(t)");
  else
    fp = xvgropen(fn, "Hydrogen Bond Autocorrelation","Time (ps)","C(t)");
  xvgr_legend(fp,asize(leg),leg);

  for(j=0; (j<nn); j++)
    fprintf(fp,"%10g  %10g  %10g  %10g  %10g\n",
	    hb->time[j]-hb->time[0],ct[j],cct[j],ght[j],kt[j]);
  fclose(fp);
  
  analyse_corr(nn,hb->time,ct,ght,kt,fit_start,temp);
  
  do_view(fn,NULL);
  sfree(rhbex);
  sfree(ct);
  sfree(gt);
  sfree(ht);
  sfree(ght);
  sfree(dght);
  sfree(cct);
  sfree(kt);
  sfree(h);
  sfree(g);
#ifndef WITHOUT_FFTW
  /*done_correl(correl);*/
#endif
}

static void init_hbframe(t_hbdata *hb,int nframes,real t)
{
  int i,j,m;
  
  hb->time[nframes]   = t;
  hb->nhb[nframes]    = 0;
  hb->ndist[nframes]  = 0;
  for (i=0; (i<max_hx); i++)
    hb->nhx[nframes][i]=0;
  /* Loop invalidated */
  if (hb->bHBmap && 0)
    for (i=0; (i<hb->d.nrd); i++)
      for (j=0; (j<hb->a.nra); j++)
	for (m=0; (m<hb->maxhydro); m++)
	  if (hb->hbmap[i][j] && hb->hbmap[i][j]->h[m])
	    set_hb(hb,i,m,j,nframes,HB_NO);
  /*set_hb(hb->hbmap[i][j]->h[m],nframes-hb->hbmap[i][j]->n0,HB_NO);*/
}

static void analyse_donor_props(char *fn,t_hbdata *hb,int nframes,real t)
{
  static FILE *fp = NULL;
  static char *leg[] = { "Nbound", "Nfree" };
  int i,j,k,nbound,nb,nhtot;
  
  if (!fn)
    return;
  if (!fp) {
    fp = xvgropen(fn,"Donor properties","Time (ps)","Number");
    xvgr_legend(fp,asize(leg),leg);
  }
  nbound = 0;
  nhtot  = 0;
  for(i=0; (i<hb->d.nrd); i++) {
    for(k=0; (k<hb->d.nhydro[i]); k++) {
      nb = 0;
      nhtot ++;
      for(j=0; (j<hb->a.nra) && (nb == 0); j++) {
	if (hb->hbmap[i][j] && hb->hbmap[i][j]->h[k] && 
	    is_hb(hb->hbmap[i][j]->h[k],nframes)) 
	  nb = 1;
      }
      nbound += nb;
    }
  }
  fprintf(fp,"%10.3e  %6d  %6d\n",t,nbound,nhtot-nbound);
}

static void dump_hbmap(t_hbdata *hb,
		       int nfile,t_filenm fnm[],bool bTwo,bool bInsert,
		       bool bContact, int isize[],int *index[],char *grpnames[],
		       t_atoms *atoms)
{
  FILE *fp,*fplog;
  int  ddd,hhh,aaa,i,j,k,m,grp;
  char ds[32],hs[32],as[32];
  
  fp = opt2FILE("-hbn",nfile,fnm,"w");
  if (opt2bSet("-g",nfile,fnm)) {
    fplog = ffopen(opt2fn("-g",nfile,fnm),"w");
    if (bContact)
      fprintf(fplog,"# %10s  %12s  %12s\n","Donor","Hydrogen","Acceptor");
    else
      fprintf(fplog,"# %10s  %12s  %12s\n","Donor","Hydrogen","Acceptor");
  }
  else
    fplog = NULL;
  for (grp=gr0; grp<=(bTwo?gr1:gr0); grp++) {
    fprintf(fp,"[ %s ]",grpnames[grp]);
    for (i=0; i<isize[grp]; i++) {
      fprintf(fp,(i%15)?" ":"\n");
      fprintf(fp," %4u",index[grp][i]+1);
    }
    fprintf(fp,"\n");
    /*
      Added -contact support below.
      - Erik Marklund, May 29, 2006
     */
    if (!bContact) {
      fprintf(fp,"[ donors_hydrogens_%s ]",grpnames[grp]);
      for (i=0; (i<hb->d.nrd); i++) {
	for(j=0; (j<hb->d.nhydro[i]); j++)
	  fprintf(fp," %4u %4u",hb->d.don[i]+1,
		  hb->d.hydro[i][j]+1);
	fprintf(fp,"\n");
      }
      fprintf(fp,"[ acceptors_%s ]",grpnames[grp]);
      for (i=0; (i<hb->a.nra); i++) {
	fprintf(fp,(i%15)?" ":"\n");
	fprintf(fp," %4u",hb->a.acc[i]+1);
      }
      fprintf(fp,"\n");
    }
  }
  if (bTwo)
    fprintf(fp,bContact ? "[ contacts_%s-%s ]\n" :
	    "[ hbonds_%s-%s ]\n",grpnames[0],grpnames[1]);
  else
    fprintf(fp,bContact ? "[ contacts_%s ]" : "[ hbonds_%s ]\n",grpnames[0]);
  
  for(i=0; (i<hb->d.nrd); i++) {
    ddd = hb->d.don[i];
    for(k=0; (k<hb->a.nra); k++) {
      aaa = hb->a.acc[k];
      for(m=0; (m<hb->d.nhydro[i]); m++) {
	if (hb->hbmap[i][k] && ISHB(hb->hbmap[i][k]->history[m])) {
	  sprintf(ds,"%s",mkatomname(atoms,ddd));
	  sprintf(as,"%s",mkatomname(atoms,aaa));
	  if (bContact) {
	    fprintf(fp," %6u %6u\n",ddd+1,aaa+1);
	    if (fplog) 
	      fprintf(fplog,"%12s  %12s\n",ds,as);
	  } else {
	    hhh = hb->d.hydro[i][m];
	    sprintf(hs,"%s",mkatomname(atoms,hhh));
	    fprintf(fp," %6u %6u %6u\n",ddd+1,hhh+1,aaa+1);
	    if (fplog) 
	      fprintf(fplog,"%12s  %12s  %12s\n",ds,hs,as);
	  }
	}
      }
    }
  }
  if (bInsert) {
    if (bTwo)
      fprintf(fp,"[ insert_%s->%s-%s ]",
	      grpnames[2],grpnames[0],grpnames[1]);
    else
      fprintf(fp,"[ insert_%s->%s ]",grpnames[2],grpnames[0]);
    j=0;
    
    /*	for(i=0; (i<hb->nrhb); i++) {
	if (hb->hb[i].bInsert) {
	fprintf(fp,(j%15)?" ":"\n");
	fprintf(fp,"%4d",i+1);
	j++;
	}
	}*/
    fprintf(fp,"\n");
  }
  fclose(fp);
  if (fplog)
    fclose(fplog);
}

int gmx_hbond(int argc,char *argv[])
{
  static char *desc[] = {
    "g_hbond computes and analyzes hydrogen bonds. Hydrogen bonds are",
    "determined based on cutoffs for the angle Acceptor - Donor - Hydrogen",
    "(zero is extended) and the distance Hydrogen - Acceptor.",
    "OH and NH groups are regarded as donors, O is an acceptor always,",
    "N is an acceptor by default, but this can be switched using",
    "[TT]-nitacc[tt]. Dummy hydrogen atoms are assumed to be connected",
    "to the first preceding non-hydrogen atom.[PAR]",
    
    "You need to specify two groups for analysis, which must be either",
    "identical or non-overlapping. All hydrogen bonds between the two",
    "groups are analyzed.[PAR]",
    
    "If you set -shell, you will be asked for an additional index group",
    "which should contain exactly one atom. In this case, only hydrogen",
    "bonds between atoms within the shell distance from the one atom are",
    "considered.[PAR]"
    
    /*    "It is also possible to analyse specific hydrogen bonds with",
    "[TT]-sel[tt]. This index file must contain a group of atom triplets",
    "Donor Hydrogen Acceptor, in the following way:[PAR]",
    */
    "[TT]",
    "[ selected ][BR]",
    "     20    21    24[BR]",
    "     25    26    29[BR]",
    "      1     3     6[BR]",
    "[tt][BR]",
    "Note that the triplets need not be on separate lines.",
    "Each atom triplet specifies a hydrogen bond to be analyzed,",
    "note also that no check is made for the types of atoms.[PAR]",
    
    "[TT]-ins[tt] turns on computing solvent insertion into hydrogen bonds.",
    "In this case an additional group must be selected, specifying the",
    "solvent molecules.[PAR]",
    
    "[BB]Output:[bb][BR]",
    "[TT]-num[tt]:  number of hydrogen bonds as a function of time.[BR]",
    "[TT]-ac[tt]:   average over all autocorrelations of the existence",
    "functions (either 0 or 1) of all hydrogen bonds.[BR]",
    "[TT]-dist[tt]: distance distribution of all hydrogen bonds.[BR]",
    "[TT]-ang[tt]:  angle distribution of all hydrogen bonds.[BR]",
    "[TT]-hx[tt]:   the number of n-n+i hydrogen bonds as a function of time",
    "where n and n+i stand for residue numbers and i ranges from 0 to 6.",
    "This includes the n-n+3, n-n+4 and n-n+5 hydrogen bonds associated",
    "with helices in proteins.[BR]",
    "[TT]-hbn[tt]:  all selected groups, donors, hydrogens and acceptors",
    "for selected groups, all hydrogen bonded atoms from all groups and",
    "all solvent atoms involved in insertion.[BR]",
    "[TT]-hbm[tt]:  existence matrix for all hydrogen bonds over all",
    "frames, this also contains information on solvent insertion",
    "into hydrogen bonds. Ordering is identical to that in [TT]-hbn[tt]",
    "index file.[BR]",
    "[TT]-dan[tt]: write out the number of donors and acceptors analyzed for",
    "each timeframe. This is especially usefull when using [TT]-shell[tt].[BR]",
    "[TT]-nhbdist[tt]: compute the number of HBonds per hydrogen in order to",
    "compare results to Raman Spectroscopy.",
    "[PAR]",
    "Note: options [TT]-ac[tt], [TT]-life[tt], [TT]-hbn[tt] and [TT]-hbm[tt]",
    "require an amount of memory proportional to the total numbers of donors",
    "times the total number of acceptors in the selected group(s)."
  };
  
  static real acut=30, abin=1, rcut=0.35, r2cut=0, rbin=0.005, rshell=-1;
  static real maxnhb=0,fit_start=1,temp=298.15;
  static bool bNitAcc=TRUE,bInsert=FALSE,bDA=TRUE,bMerge=TRUE;
  static int  nDump=0;
  static bool bContact=FALSE;
  /* options */
  t_pargs pa [] = {
    { "-ins",  FALSE,  etBOOL, {&bInsert},
      "Analyze solvent insertion" },
    { "-a",    FALSE,  etREAL, {&acut},
      "Cutoff angle (degrees, Acceptor - Donor - Hydrogen)" },
    { "-r",    FALSE,  etREAL, {&rcut},
      "Cutoff radius (nm, X - Acceptor, see next option)" },
    { "-da",   FALSE,  etBOOL, {&bDA},
      "Use distance Donor-Acceptor (if TRUE) or Hydrogen-Acceptor (FALSE)" },
    { "-r2",   FALSE,  etREAL, {&r2cut},
      "Second cutoff radius. Mainly useful with -contact and -ac"},
    { "-abin", FALSE,  etREAL, {&abin},
      "Binwidth angle distribution (degrees)" },
    { "-rbin", FALSE,  etREAL, {&rbin},
      "Binwidth distance distribution (nm)" },
    { "-nitacc",FALSE, etBOOL, {&bNitAcc},
      "Regard nitrogen atoms as acceptors" },
    { "-contact",FALSE,etBOOL, {&bContact},
      "Do not look for hydrogen bonds, but merely for contacts within the cut-off distance" },
    { "-shell", FALSE, etREAL, {&rshell},
      "when > 0, only calculate hydrogen bonds within # nm shell around "
      "one particle" },
    { "-fitstart", FALSE, etREAL, {&fit_start},
      "Time from which to start fitting the correlation functions in order to obtain the forward and backward rate constants for HB breaking and formation" }, 
    { "-temp",  FALSE, etREAL, {&temp},
      "Temperature (K) for computing the Gibbs energy corresponding to HB breaking and reforming" },
    { "-dump",  FALSE, etINT, {&nDump},
      "Dump the first N hydrogen bond ACFs in a single xvg file for debugging" },
    { "-max_hb",FALSE, etREAL, {&maxnhb},
      "Theoretical maximum number of hydrogen bonds used for normalizing HB autocorrelation function. Can be useful in case the program estimates it wrongly" },
    { "-merge", FALSE, etBOOL, {&bMerge},
      "H-bonds between the same donor and acceptor, but with different hydrogen are treated as a single H-bond. Mainly important for the ACF." }
  };
  static char *bugs[] = {
    "The option [TT]-sel[tt] that used to work on selected hbonds is out of order, and therefore not available for the time being."
  };
  t_filenm fnm[] = {
    { efTRX, "-f",   NULL,     ffREAD  },
    { efTPX, NULL,   NULL,     ffREAD  },
    { efNDX, NULL,   NULL,     ffOPTRD },
    /*    { efNDX, "-sel", "select", ffOPTRD },*/
    { efXVG, "-num", "hbnum",  ffWRITE },
    { efXVG, "-ac",  "hbac",   ffOPTWR },
    { efXVG, "-dist","hbdist", ffOPTWR },
    { efXVG, "-ang", "hbang",  ffOPTWR },
    { efXVG, "-hx",  "hbhelix",ffOPTWR },
    { efNDX, "-hbn", "hbond",  ffOPTWR },
    { efXPM, "-hbm", "hbmap",  ffOPTWR },
    { efXVG, "-don", "donor",  ffOPTWR },
    { efXVG, "-dan", "danum",  ffOPTWR },
    { efXVG, "-life","hblife", ffOPTWR },
    { efXVG, "-nhbdist", "nhbdist",ffOPTWR }
    
  };
#define NFILE asize(fnm)
  
  char  hbmap [HB_NR]={ ' ',    'o',      '-',       '*' };
  char *hbdesc[HB_NR]={ "None", "Present", "Inserted", "Present & Inserted" };
  t_rgb hbrgb [HB_NR]={ {1,1,1},{1,0,0},   {0,0,1},    {1,0,1} };
  
  int     status;
  t_topology top;
  t_inputrec ir;
  t_pargs *ppa;
  int     npargs,natoms,nframes=0,shatom;
  int     *isize;
  char    **grpnames;
  atom_id **index;
  rvec    *x,hbox;
  matrix  box;
  real    t,ccut,dist,ang;
  double  max_nhb,aver_nhb,aver_dist;
  int     h,i,j,k,l,start,end,id,ja,ogrp,nsel;
  int     xi,yi,zi,ai;
  int     xj,yj,zj,aj,xjj,yjj,zjj;
  int     xk,yk,zk,ak,xkk,ykk,zkk;
  bool    bSelected,bHBmap,bStop,bTwo,new,was,bBox,bTric;
  int     *adist,*rdist;
  int        grp,nabin,nrbin,bin,resdist,ihb;
  char       **leg;
  t_hbdata   *hb;
  FILE       *fp,*fpins=NULL,*fplog,*fpnhb=NULL;
  t_gridcell ***grid;
  t_ncell    *icell,*jcell,*kcell;
  ivec       ngrid;
    
  CopyRight(stdout,argv[0]);

  npargs = asize(pa);  
  ppa    = add_acf_pargs(&npargs,pa);
  
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_BE_NICE,NFILE,fnm,npargs,
		    ppa,asize(desc),desc,asize(bugs),bugs);

  /* process input */
  bSelected = opt2bSet("-sel",NFILE,fnm);
  ccut = cos(acut*DEG2RAD);
  
  if (bContact) {
    if (bSelected)
      gmx_fatal(FARGS,"Can not analyze selected contacts: turn off -sel");
    if (bInsert)
      gmx_fatal(FARGS,"Can not analyze inserted contacts: turn off -ins");
  }
  
  /* Initiate main data structure! */
  bHBmap = (opt2bSet("-ac",NFILE,fnm) ||
	    opt2bSet("-life",NFILE,fnm) ||
	    opt2bSet("-hbn",NFILE,fnm) || 
	    opt2bSet("-hbm",NFILE,fnm));
  
  if (opt2bSet("-nhbdist",NFILE,fnm)) {
    char *leg[MAXHH+1] = { "0 HBs", "1 HB", "2 HBs", "3 HBs", "Total" };
    fpnhb = xvgropen(opt2fn("-nhbdist",NFILE,fnm),
		     "Number of donor-H with N HBs","Time (ps)","N");
    xvgr_legend(fpnhb,asize(leg),leg);
  }
  
  hb = mk_hbdata(bHBmap,opt2bSet("-dan",NFILE,fnm),bMerge || bContact);
  
  /* get topology */
  read_tpx(ftp2fn(efTPX,NFILE,fnm),&i,&t,&t,
	   &ir,box,&natoms,NULL,NULL,NULL,&top);
  
  snew(grpnames,grNR);
  snew(index,grNR);
  snew(isize,grNR);
  if (bSelected) {
    /* analyze selected hydrogen bonds */
    printf("Select group with selected atoms:\n");
    get_index(&(top.atoms),opt2fn("-sel",NFILE,fnm),
	      1,&nsel,index,grpnames);
    if (nsel % 3)
      gmx_fatal(FARGS,"Number of atoms in group '%s' not a multiple of 3\n"
		  "and therefore cannot contain triplets of "
		  "Donor-Hydrogen-Acceptor",grpnames[0]);
    bTwo=FALSE;
    
    for(i=0; (i<nsel); i+=3) {
      int dd = index[0][i];
      int hh = index[0][i+1];
      int aa = index[0][i+2];
      add_dh (&hb->d,dd,hh,i);
      add_acc(&hb->a,aa,i);
      /* Should this be here ? */
      snew(hb->d.dptr,top.atoms.nr);
      snew(hb->a.aptr,top.atoms.nr);
      add_hbond(hb,dd,aa,hh,gr0,gr0,0,FALSE,bMerge,0,bContact);
    }
    printf("Analyzing %d selected hydrogen bonds from '%s'\n",
	   isize[0],grpnames[0]);
  } 
  else {
    /* analyze all hydrogen bonds: get group(s) */
    printf("Specify 2 groups to analyze:\n");
    get_index(&(top.atoms),ftp2fn_null(efNDX,NFILE,fnm),
	      2,isize,index,grpnames);
    
    /* check if we have two identical or two non-overlapping groups */
    bTwo = isize[0] != isize[1];
    for (i=0; (i<isize[0]) && !bTwo; i++)
      bTwo = index[0][i] != index[1][i];
    if (bTwo) {
      printf("Checking for overlap in atoms between %s and %s\n",
	     grpnames[0],grpnames[1]);
      for (i=0; i<isize[0]; i++)
	for (j=0; j<isize[1]; j++)
	  if (index[0][i] == index[1][j]) 
	    gmx_fatal(FARGS,"Partial overlap between groups '%s' and '%s'",
			grpnames[0],grpnames[1]);
    }
    if (bTwo)
      printf("Calculating %s "
	     "between %s (%d atoms) and %s (%d atoms)\n",
	     bContact ? "contacts" : "hydrogen bonds",
	     grpnames[0],isize[0],grpnames[1],isize[1]);
    else
      fprintf(stderr,"Calculating %s in %s (%d atoms)\n",
	      bContact?"contacts":"hydrogen bonds",grpnames[0],isize[0]);
  }
  if (bInsert) {
    printf("Specify group for insertion analysis:\n");
    get_index(&(top.atoms),ftp2fn_null(efNDX,NFILE,fnm),
	      1,&(isize[grI]),&(index[grI]),&(grpnames[grI]));
    printf("Checking for overlap...\n");
    for (i=0; i<isize[grI]; i++)
      for (grp=0; grp<(bTwo?2:1); grp++)
	for (j=0; j<isize[grp]; j++)
	  if (index[grI][i] == index[grp][j]) 
	    gmx_fatal(FARGS,"Partial overlap between groups '%s' and '%s'",
			grpnames[grp],grpnames[grI]);
    fpins=ffopen("insert.dat","w");
    fprintf(fpins,"%4s: %15s -> %15s (%7s) - %15s (%7s)\n",
	    "time","insert","donor","distang","acceptor","distang");
  }
  
  /* search donors and acceptors in groups */
  for (i=0; (i<grNR); i++)
    if ( ((i==gr0) && !bSelected ) ||
	 ((i==gr1) && bTwo ) ||
	 ((i==grI) && bInsert ) ) {
      if (bContact) {
	search_acceptors(&top,isize[i],index[i],&hb->a,i,
			 bNitAcc,TRUE,(bTwo && (i==gr0)) || !bTwo);
	search_donors   (&top,isize[i],index[i],&hb->d,i,
			 TRUE,(bTwo && (i==gr1)) || !bTwo);
      }
      else {
	search_acceptors(&top,isize[i],index[i],&hb->a,i,bNitAcc,FALSE,TRUE);
	search_donors   (&top,isize[i],index[i],&hb->d,i,FALSE,TRUE);
      }
    }
  printf("Found %d donors and %d acceptors\n",hb->d.nrd,hb->a.nra);
  /*if (bSelected)
    snew(donors[gr0D], dons[gr0D].nrd);*/

  if (bHBmap) {
    /* Generate hbond data structure */
    mk_hbmap(hb,bTwo,bInsert);
  }
  
  /* check input */
  bStop=FALSE;
  if (hb->d.nrd + hb->a.nra == 0) {
    printf("No Donors or Acceptors found\n");
    bStop=TRUE;
  }
  if (!bStop) {
    if (hb->d.nrd == 0) {
      printf("No Donors found\n");
      bStop=TRUE;
    }
    if (hb->a.nra == 0) {
      printf("No Acceptors found\n");
      bStop=TRUE;
    }
  }
  if (bStop)
    gmx_fatal(FARGS,"Nothing to be done");

  shatom=0;
  if (rshell > 0) {
    int shisz;
    atom_id *shidx;
    char *shgrpnm;
    /* get index group with atom for shell */
    do {
      printf("Select atom for shell (1 atom):\n");
      get_index(&(top.atoms),ftp2fn_null(efNDX,NFILE,fnm),
		1,&shisz,&shidx,&shgrpnm);
      if (shisz != 1)
	printf("group contains %d atoms, should be 1 (one)\n",shisz);
    } while(shisz != 1);
    shatom = shidx[0];
    printf("Will calculate hydrogen bonds within a shell "
	   "of %g nm around atom %i\n",rshell,shatom+1);
  }

  /* Analyze trajectory */
  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  if ( natoms > top.atoms.nr )
    gmx_fatal(FARGS,"Topology (%d atoms) does not match trajectory (%d atoms)",
	      top.atoms.nr,natoms);
		
  bBox  = ir.ePBC!=epbcNONE;
  grid  = init_grid(bBox, box, (rcut>r2cut)?rcut:r2cut, ngrid);
  nabin = acut/abin;
  nrbin = rcut/rbin;
  snew(adist,nabin+1);
  snew(rdist,nrbin+1);
  
  do {
    bTric = bBox && TRICLINIC(box);
    build_grid(hb,x,x[shatom], bBox,box,hbox, (rcut>r2cut)?rcut:r2cut, rshell, ngrid,grid);
    
    reset_nhbonds(&(hb->d));

    if (debug && bDebug)
      dump_grid(debug, ngrid, grid);
    
    add_frames(hb,nframes);
    init_hbframe(hb,nframes,t);

    if (hb->bDAnr)
      count_da_grid(ngrid, grid, hb->danr[nframes]);
      
    if (bSelected) {
      int ii;
      
      for(ii=0; (ii<nsel); ii++) {
	int dd = index[0][i];
	int hh = index[0][i+1];
	int aa = index[0][i+2];
	ihb = is_hbond(hb,ii,ii,dd,aa,rcut,r2cut,ccut,x,bBox,box,
		       hbox,&dist,&ang,bDA,&h,bContact);
	
	if (ihb) {
	  /* add to index if not already there */
	  /* Add a hbond */
	  add_hbond(hb,dd,aa,hh,ii,ii,nframes,FALSE,bMerge,ihb,bContact);
	}
      }
    }
    else {
      /* loop over all gridcells (xi,yi,zi)      */
      /* Removed confusing macro, DvdS 27/12/98  */
      for(xi=0; (xi<ngrid[XX]); xi++)
	for(yi=0; (yi<ngrid[YY]); yi++)
	  for(zi=0; (zi<ngrid[ZZ]); zi++) {
	    
	    /* loop over donor groups gr0 (always) and gr1 (if necessary) */
	    for (grp=gr0; (grp <= (bTwo?gr1:gr0)); grp++) {
	      icell=&(grid[zi][yi][xi].d[grp]);
	      
	      if (bTwo)
		ogrp = 1-grp;
	      else
		ogrp = grp;
	      
	      /* loop over all hydrogen atoms from group (grp) 
	       * in this gridcell (icell) 
	       */
	      for (ai=0; (ai<icell->nr); ai++) {
		i  = icell->atoms[ai];
		
		/* loop over all adjacent gridcells (xj,yj,zj) */
		/* This is a macro!!! */
		LOOPGRIDINNER(xj,yj,zj,xjj,yjj,zjj,xi,yi,zi,ngrid,bTric) {
		  jcell=&(grid[zj][yj][xj].a[ogrp]);
		  /* loop over acceptor atoms from other group (ogrp) 
		   * in this adjacent gridcell (jcell) 
		   */
		  for (aj=0; (aj<jcell->nr); aj++) {
		    j = jcell->atoms[aj];
		  
		    /* check if this once was a h-bond */
		    ihb = is_hbond(hb,grp,ogrp,i,j,rcut,r2cut,ccut,x,bBox,box,
				   hbox,&dist,&ang,bDA,&h,bContact);
		    
		    if (ihb) {
		      /* add to index if not already there */
		      /* Add a hbond */
		      add_hbond(hb,i,j,h,grp,ogrp,nframes,FALSE,bMerge,ihb,bContact);
		      
		      /* make angle and distance distributions */
		      if (ihb == hbHB && !bContact) {
			ang*=RAD2DEG;
			adist[(int)( ang/abin)]++;
			rdist[(int)(dist/rbin)]++;
			
			if (!bTwo) {
			  int id,ia;
			  if ((id = donor_index(&hb->d,grp,i)) == NOTSET)
			    gmx_fatal(FARGS,"Invalid donor %d",i);
			  if ((ia = acceptor_index(&hb->a,ogrp,j)) == NOTSET)
			    gmx_fatal(FARGS,"Invalid acceptor %d",j);
			  resdist=abs(top.atoms.atom[i].resnr-
				      top.atoms.atom[j].resnr);
			  if (resdist >= max_hx) 
			    resdist = max_hx-1;
			  hb->nhx[nframes][resdist]++;
			}
		      }
		      if (bInsert && bSelected) {
			/* this has been a h-bond, or we are analyzing 
			   selected bonds: check for inserted */
			bool ins_d, ins_a;
			real ins_d_dist, ins_d_ang, ins_a_dist, ins_a_ang;
			int  ins_d_k=0,ins_a_k=0;
			
			ins_d=ins_a=FALSE;
			ins_d_dist=ins_d_ang=ins_a_dist=ins_a_ang=1e6;
			
			/* loop over gridcells adjacent to i (xk,yk,zk) */
			LOOPGRIDINNER(xk,yk,zk,xkk,ykk,zkk,xi,yi,zi,ngrid,bTric){
			  kcell=&(grid[zk][yk][xk].a[grI]);
			  /* loop over acceptor atoms from ins group 
			     in this adjacent gridcell (kcell) */
			  for (ak=0; (ak<kcell->nr); ak++) {
			    k=kcell->atoms[ak];
			    ihb = is_hbond(hb,grp,grI,i,k,rcut,r2cut,ccut,x,
					   bBox,box,hbox,&dist,&ang,bDA,&h,
					   bContact);
			    if (ihb == hbHB) {
			      if (dist < ins_d_dist) {
				ins_d=TRUE;
				ins_d_dist=dist;
				ins_d_ang =ang ;
				ins_d_k   =k   ;
			      }
			    }
			  }
			}
			ENDLOOPGRIDINNER;
			/* loop over gridcells adjacent to j (xk,yk,zk) */
			LOOPGRIDINNER(xk,yk,zk,xkk,ykk,zkk,xj,yj,zj,ngrid,bTric){
			  kcell=&grid[zk][yk][xk].d[grI];
			  /* loop over hydrogen atoms from ins group 
			     in this adjacent gridcell (kcell) */
			  for (ak=0; ak<kcell->nr; ak++) {
			    k   = kcell->atoms[ak];
			    ihb = is_hbond(hb,grI,ogrp,k,j,rcut,r2cut,ccut,x,
					   bBox,box,hbox,&dist,&ang,bDA,&h,
					   bContact);
			    if (ihb == hbHB) {
			      if (dist<ins_a_dist) {
				ins_a=TRUE;
				ins_a_dist=dist;
				ins_a_ang =ang ;
				ins_a_k   =k   ;
			      }
			    }
			  }
			}
			ENDLOOPGRIDINNER;
			
			{
			  ihb = is_hbond(hb,grI,grI,ins_d_k,ins_a_k,rcut,r2cut,ccut,x,
					 bBox,box,hbox,&dist,&ang,bDA,&h,bContact);
			  if (ins_d && ins_a && ihb) {
			    /* add to hbond index if not already there */
			    add_hbond(hb,ins_d_k,ins_a_k,h,grI,ogrp,
				      nframes,TRUE,bMerge,ihb,bContact);
			    
			    /* print insertion info to file */
			    /*fprintf(fpins,
			      "%4g: %4u:%3.3s%4d%3.3s -> "
			      "%4u:%3.3s%4d%3.3s (%4.2f,%2.0f) - "
			      "%4u:%3.3s%4d%3.3s (%4.2f,%2.0f)\n",t,
			      a[grIA][ins_d_k]+1,
			      *top.atoms.resname[top.atoms.atom[a[grIA][ins_d_k]].resnr],
			      top.atoms.atom[a[grIA][ins_d_k]].resnr+1,
			      *top.atoms.atomname[a[grIA][ins_d_k]],
			      a[grp+grD][i]+1,
			      *top.atoms.resname[top.atoms.atom[a[grp+grD][i]].resnr],
			      top.atoms.atom[a[grp+grD][i]].resnr+1,
			      *top.atoms.atomname[a[grp+grD][i]],
			      ins_d_dist,ins_d_ang*RAD2DEG,
			      a[ogrp+grA][j]+1,
			      *top.atoms.resname[top.atoms.atom[a[ogrp+grA][j]].resnr],
			      top.atoms.atom[a[ogrp+grA][j]].resnr+1,
			      *top.atoms.atomname[a[ogrp+grA][j]],
			      ins_a_dist,ins_a_ang*RAD2DEG);*/
			  }
			}
		      }
		    }
		  } /* for aj  */
		}
		ENDLOOPGRIDINNER;
	      } /* for ai  */
	    } /* for grp */
	  } /* for xi,yi,zi */
    }
    analyse_donor_props(opt2fn_null("-don",NFILE,fnm),hb,nframes,t);
    if (fpnhb)
      do_nhb_dist(fpnhb,hb,t);
    
    nframes++;
  } while (read_next_x(status,&t,natoms,x,box));
  
  free_grid(ngrid,&grid);
  
  close_trj(status);
  if (bInsert)
    fclose(fpins);
  if (fpnhb)
    fclose(fpnhb);
    
  /* Compute maximum possible number of different hbonds */
  if (maxnhb > 0)
    max_nhb = maxnhb;
  else {
    max_nhb = 0.5*(hb->d.nrd*hb->a.nra);
  }
  /* Added support for -contact below.
   * - Erik Marklund, May 29-31, 2006 */
  /* Changed contact code.
   * - Erik Marklund, June 29, 2006 */
  if (bHBmap) {
    if (hb->nrhb==0) {
      printf("No %s found!!\n", bContact ? "contacts" : "hydrogen bonds");
    } else {
	printf("Found %d different %s in trajectory\n"
	       "Found %d different atom-pairs within %s distance\n",
	       hb->nrhb, bContact?"contacts":"hydrogen bonds",
	       hb->nrdist,(r2cut>0)?"second cut-off":"hydrogen bonding");

	if (bMerge)
	  merge_hb(hb,bTwo,bContact);

	if (opt2bSet("-hbn",NFILE,fnm)) 
	  dump_hbmap(hb,NFILE,fnm,bTwo,bInsert,bContact,isize,index,grpnames,&top.atoms);

	/* Moved the call to merge_hb() to a line BEFORE dump_hbmap
	 * to make the -hbn and -hmb output match eachother. 
	 * - Erik Marklund, May 30, 2006 */
    }
  }
  /* Print out number of hbonds and distances */
  aver_nhb  = 0;    
  aver_dist = 0;
  fp = xvgropen(opt2fn("-num",NFILE,fnm),bContact ? "Contacts" :
		"Hydrogen Bonds","Time","Number");
  snew(leg,2);
  snew(leg[0],STRLEN);
  snew(leg[1],STRLEN);
  sprintf(leg[0],"%s",bContact?"Contacts":"Hydrogen bonds");
  sprintf(leg[1],"Pairs within %g nm",(r2cut>0)?r2cut:rcut);
  xvgr_legend(fp,2,leg);
  sfree(leg[1]);
  sfree(leg[0]);
  sfree(leg);
  for(i=0; (i<nframes); i++) {
    fprintf(fp,"%10g  %10d  %10d\n",hb->time[i],hb->nhb[i],hb->ndist[i]);
    aver_nhb  += hb->nhb[i];
    aver_dist += hb->ndist[i];
  }
  fclose(fp);
  aver_nhb  /= nframes;
  aver_dist /= nframes;
  /* Print HB distance distribution */
  if (opt2bSet("-dist",NFILE,fnm)) {
    long sum;
    
    sum=0;
    for(i=0; i<nrbin; i++)
      sum+=rdist[i];
    
    fp = xvgropen(opt2fn("-dist",NFILE,fnm),
		  "Hydrogen Bond Distribution",
		  "Hydrogen - Acceptor Distance (nm)","");
    for(i=0; i<nrbin; i++)
      fprintf(fp,"%10g %10g\n",(i+0.5)*rbin,rdist[i]/(rbin*(real)sum));
    fclose(fp);
  }
  
  /* Print HB angle distribution */
  if (opt2bSet("-ang",NFILE,fnm)) {
    long sum;
    
    sum=0;
    for(i=0; i<nabin; i++)
      sum+=adist[i];
    
    fp = xvgropen(opt2fn("-ang",NFILE,fnm),
		  "Hydrogen Bond Distribution",
		  "Donor - Hydrogen - Acceptor Angle (\\SO\\N)","");
    for(i=0; i<nabin; i++)
      fprintf(fp,"%10g %10g\n",(i+0.5)*abin,adist[i]/(abin*(real)sum));
    fclose(fp);
  }
    
  /* Print HB in alpha-helix */
  if (opt2bSet("-hx",NFILE,fnm)) {
    fp = xvgropen(opt2fn("-hx",NFILE,fnm),
		  "Hydrogen Bonds","Time(ps)","Count");
    xvgr_legend(fp,NRHXTYPES,hxtypenames);
    for(i=0; i<nframes; i++) {
      fprintf(fp,"%10g",hb->time[i]);
      for (j=0; j<max_hx; j++)
	fprintf(fp," %6d",hb->nhx[i][j]);
      fprintf(fp,"\n");
    }
    fclose(fp);
  }
  printf("Average number of %s per timeframe %.3f out of %g possible\n",
	 bContact ? "contacts" : "hbonds",
	 bContact ? aver_dist : aver_nhb, max_nhb);
	 
  /* Do Autocorrelation etc. */
  if (hb->bHBmap) {
    /* 
       Added support for -contact in ac and hbm calculations below.
       - Erik Marklund, May 29, 2006
     */

    if (opt2bSet("-ac",NFILE,fnm) || opt2bSet("-life",NFILE,fnm))
      please_cite(stdout,"Spoel2006b");
    if (opt2bSet("-ac",NFILE,fnm)) 
      do_hbac(opt2fn("-ac",NFILE,fnm),hb,aver_nhb/max_nhb,aver_dist,nDump,
	      bMerge,bContact,fit_start,temp,r2cut>0);
    if (opt2bSet("-life",NFILE,fnm))
      do_hblife(opt2fn("-life",NFILE,fnm),hb,bMerge,bContact);
    if (opt2bSet("-hbm",NFILE,fnm)) {
      t_matrix mat;
      int id,ia,hh,x,y;
      
      mat.nx=nframes;
      mat.ny=(bContact ? hb->nrdist : hb->nrhb);

      snew(mat.matrix,mat.nx);
      for(x=0; (x<mat.nx); x++) 
	snew(mat.matrix[x],mat.ny);
      y=0;
      for(id=0; (id<hb->d.nrd); id++) 
	for(ia=0; (ia<hb->a.nra); ia++) {
	  for(hh=0; (hh<hb->maxhydro); hh++) {
	    if (hb->hbmap[id][ia]) {
	      if (ISHB(hb->hbmap[id][ia]->history[hh])) {
		/* Changed '<' into '<=' in the for-statement below.
		 * It fixed the previously undiscovered bug that caused
		 * the last occurance of an hbond/contact to not be
		 * set in mat.matrix. Have a look at any old -hbm-output
		 * and you will notice that the last column is allways empty.
		 * - Erik Marklund May 30, 2006
		 */
		for(x=0; (x<=hb->hbmap[id][ia]->nframes); x++) {
		  int nn0 = hb->hbmap[id][ia]->n0;
		  range_check(y,0,mat.ny);
		  mat.matrix[x+nn0][y] = is_hb(hb->hbmap[id][ia]->h[hh],x);
		}
		y++;
	      }
	    }
	  }
	}
      mat.axis_x=hb->time;
      snew(mat.axis_y,mat.ny);
      for(j=0; j<mat.ny; j++)
	mat.axis_y[j]=j;
      sprintf(mat.title,bContact ? "Contact Existence Map":
	      "Hydrogen Bond Existence Map");
      sprintf(mat.legend,bContact ? "Contacts" : "Hydrogen Bonds");
      sprintf(mat.label_x,"Time (ps)");
      sprintf(mat.label_y, bContact ? "Contact Index" : "Hydrogen Bond Index");
      mat.bDiscrete=TRUE;
      if (bInsert)
	mat.nmap=HB_NR;
      else
	mat.nmap=2;
      snew(mat.map,mat.nmap);
      for(i=0; i<mat.nmap; i++) {
	mat.map[i].code.c1=hbmap[i];
	mat.map[i].desc=hbdesc[i];
	mat.map[i].rgb=hbrgb[i];
      }
      fp = opt2FILE("-hbm",NFILE,fnm,"w");
      write_xpm_m(fp, mat);
      fclose(fp);
      for(x=0; x<mat.nx; x++)
	sfree(mat.matrix[x]);
      sfree(mat.axis_y);
      sfree(mat.matrix);
      sfree(mat.map);
    }
  }
    
  if (hb->bDAnr) {
    int  i,j,nleg;
    char **legnames;
    char buf[STRLEN];
    
#define USE_THIS_GROUP(j) ( (j == gr0) || (bTwo && (j == gr1)) || (bInsert && (j == grI)) )
    
    fp = xvgropen(opt2fn("-dan",NFILE,fnm),
		  "Donors and Acceptors","Time(ps)","Count");
    nleg = (bTwo?2:1)*2 + (bInsert?2:0);
    snew(legnames,nleg);
    i=0;
    for(j=0; j<grNR; j++)
      if ( USE_THIS_GROUP(j) ) {
	sprintf(buf,"Donors %s",grpnames[j]);
	legnames[i++] = strdup(buf);
	sprintf(buf,"Acceptors %s",grpnames[j]);
	legnames[i++] = strdup(buf);
      }
    if (i != nleg)
      gmx_incons("number of legend entries");
    xvgr_legend(fp,nleg,legnames);
    for(i=0; i<nframes; i++) {
      fprintf(fp,"%10g",hb->time[i]);
      for (j=0; (j<grNR); j++)
	if ( USE_THIS_GROUP(j) )
	  fprintf(fp," %6d",hb->danr[i][j]);
      fprintf(fp,"\n");
    }
    fclose(fp);
  }
  
  thanx(stdout);
  
  return 0;
}
