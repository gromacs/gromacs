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
 * Good ROcking Metal Altar for Chronical Sinners */
static char *SRCID_ns_c = "$Id$";

#include <math.h>
#include <string.h>
#include "sysstuff.h"
#include "assert.h"
#include "smalloc.h"
#include "macros.h"
#include "maths.h"
#include "vec.h"
#include "network.h"
#include "nsgrid.h"
#include "force.h"
#include "ns.h"
#include "pbc.h"
#include "names.h"
#include "fatal.h"
#include "nrnb.h"
#include "ns.h"
#include "force.h"
#include "txtdump.h"

#define MAX_CG 1024

typedef struct {
  int     ncg;
  int     nj;
  atom_id jcg[MAX_CG];
} t_ns_buf;

/* 
 *    E X C L U S I O N   H A N D L I N G
 */
typedef unsigned long t_excl;

#ifdef DEBUG
static void SETEXCL_(t_excl e[],atom_id i,atom_id j)
{   e[j] = e[j] | (1<<i); }
static void RMEXCL_(t_excl e[],atom_id i,atom_id j) 
{ e[j]=e[j] & ~(1<<i); }
static bool ISEXCL_(t_excl e[],atom_id i,atom_id j) 
{ return (bool)(e[j] & (1<<i)); }
static bool NOTEXCL_(t_excl e[],atom_id i,atom_id j)
{  return !(ISEXCL(e,i,j)); }
#else
#define SETEXCL(e,i,j) (e)[((atom_id) j)] |= (1<<((atom_id) i))
#define RMEXCL(e,i,j)  (e)[((atom_id) j)] &= (~(1<<((atom_id) i)))
#define ISEXCL(e,i,j)  (bool) ((e)[((atom_id) j)] & (1<<((atom_id) i)))
#define NOTEXCL(e,i,j) !(ISEXCL(e,i,j))
#endif

/************************************************
 *
 *  U T I L I T I E S    F O R    N S
 *
 ************************************************/

static int NLJ_INC = 16384;

static void init_nblist(t_nblist *nl,int homenr,int il_code)
{
  nl->il_code = il_code;
  nl->maxnri  = homenr*MAXSHIFT;
  nl->maxnrj  = 0;
  nl->nri     = 0;
  nl->nrj     = 0;
  snew(nl->iinr,   nl->maxnri+2);
  snew(nl->gid,    nl->maxnri+2);
  snew(nl->shift,  nl->maxnri+2);
  snew(nl->jindex, nl->maxnri+2);
  if (nl->maxnri > 0)
    nl->iinr[0] = -1;
}

static unsigned int nbf_index(bool bCoul,bool bRF,bool bBham,
			      bool bTab,bool bWater,bool bEwald)
{
  /* lot of redundancy here,
   * since we cant have RF and EWALD simultaneously...
   */

  int inloop[64] = { 
    eNR_LJC,        eNR_QQ,         eNR_BHAM,       eNR_QQ,
    eNR_LJCRF,      eNR_QQRF,       eNR_BHAMRF,     eNR_QQRF,
    eNR_TAB,        eNR_COULTAB,    eNR_BHAMTAB,    eNR_COULTAB,
    eNR_TAB,        eNR_COULTAB,    eNR_BHAMTAB,    eNR_COULTAB,
    eNR_LJC_WAT,    eNR_QQ_WAT,     eNR_BHAM_WAT,   eNR_QQ_WAT,
    eNR_LJCRF_WAT,  eNR_QQRF_WAT,   eNR_BHAMRF_WAT, eNR_QQRF_WAT,
    eNR_TAB_WAT,    eNR_COULTAB_WAT,eNR_BHAMTAB_WAT,eNR_COULTAB_WAT,
    eNR_TAB_WAT,    eNR_COULTAB_WAT,eNR_BHAMTAB_WAT,eNR_COULTAB_WAT,
    eNR_LJC_EW,     eNR_QQ_EW,      eNR_BHAM_EW,    eNR_QQ_EW,
    eNR_LJC_EW,     eNR_QQ_EW,      eNR_BHAM_EW,    eNR_QQ_EW,
    eNR_TAB,        eNR_COULTAB,    eNR_BHAMTAB,    eNR_COULTAB,
    eNR_TAB,        eNR_COULTAB,    eNR_BHAMTAB,    eNR_COULTAB,
    eNR_LJC_WAT_EW, eNR_QQ_WAT_EW,  eNR_BHAM_WAT_EW,eNR_QQ_WAT_EW,
    eNR_LJC_WAT_EW, eNR_QQ_WAT_EW,  eNR_BHAM_WAT_EW,eNR_QQ_WAT_EW,
    eNR_TAB_WAT,    eNR_COULTAB_WAT,eNR_BHAMTAB_WAT,eNR_COULTAB_WAT,
    eNR_TAB_WAT,    eNR_COULTAB_WAT,eNR_BHAMTAB_WAT,eNR_COULTAB_WAT
  };

  unsigned int ni;
  
  ni = bCoul | (bBham << 1) | (bRF << 2) | (bTab << 3) | (bWater << 4)
      | (bEwald << 5);
  
  return inloop[ni];
}

void init_neighbor_list(FILE *log,t_forcerec *fr,int homenr)
{
  /* Make this tunable! (does not seem to be a big difference though) 
   * This parameter determines the number of i particles in a long range 
   * neighbourlist. Too few means many function calls, too many means
   * cache trashing.
   */
  int maxlr=max(homenr/20,50);
  
  init_nblist(&fr->nlist_sr[eNL_VDW],(homenr-fr->nWatMol*3),
	      nbf_index(FALSE,fr->bRF,fr->bBHAM,fr->bTab,FALSE,fr->bEwald));
  init_nblist(&fr->nlist_sr[eNL_QQ],(homenr-fr->nWatMol*3),
	      nbf_index(TRUE,fr->bRF,fr->bBHAM,fr->bTab,FALSE,fr->bEwald));
  if (fr->bPert)
    init_nblist(&fr->nlist_sr[eNL_FREE],homenr,
		fr->bBHAM ? eNR_BHAM_FREE : eNR_LJC_FREE);
  if (fr->bWaterOpt) {
    init_nblist(&fr->nlist_sr[eNL_VDW_WAT],fr->nWatMol,
		nbf_index(FALSE,fr->bRF,fr->bBHAM,fr->bTab,TRUE,fr->bEwald));
    init_nblist(&fr->nlist_sr[eNL_QQ_WAT],fr->nWatMol,
		nbf_index(TRUE,fr->bRF,fr->bBHAM,fr->bTab,TRUE,fr->bEwald));
  }
  if (fr->bTwinRange) {
    fprintf(log,"Allocating space for long range neighbour list of %d atoms\n",
	    maxlr);
    init_nblist(&fr->nlist_lr[eNL_VDW],maxlr,
		nbf_index(FALSE,fr->bRF,fr->bBHAM,fr->bTab,FALSE,fr->bEwald));
    init_nblist(&fr->nlist_lr[eNL_QQ],maxlr,
		nbf_index(TRUE,fr->bRF,fr->bBHAM,fr->bTab,FALSE,fr->bEwald));
    if (fr->bPert)
      init_nblist(&fr->nlist_lr[eNL_FREE],maxlr,
		  fr->bBHAM ? eNR_BHAM_FREE : eNR_LJC_FREE);
    if (fr->bWaterOpt) {
      init_nblist(&fr->nlist_lr[eNL_VDW_WAT],maxlr,
		  nbf_index(FALSE,fr->bRF,fr->bBHAM,fr->bTab,TRUE,fr->bEwald));
      init_nblist(&fr->nlist_lr[eNL_QQ_WAT],maxlr,
		  nbf_index(TRUE,fr->bRF,fr->bBHAM,fr->bTab,TRUE,fr->bEwald));
    }
  }
}

static void reset_nblist(t_nblist *nl)
{
  nl->nri       = 0;
  nl->nrj       = 0;
  if (nl->maxnri > 0) {
    nl->iinr[0]   = -1;
    nl->jindex[0] = 0;
    nl->jindex[1] = 0;
  }
}

static void reset_neighbor_list(t_forcerec *fr,bool bLR,int eNL)
{
  int i;
  
  if (bLR) 
    reset_nblist(&(fr->nlist_lr[eNL]));
  else {
    for(i=0; (i<eNL_NR); i++)
      reset_nblist(&(fr->nlist_sr[i]));
  }
}

static gmx_inline void new_i_nblist(t_nblist *nlist,
				    int ftype,int i_atom,int shift,int gid)
{
  int    i,nri,nshift;
    
  if (nlist->maxnrj <= nlist->nrj + NLJ_INC-1) {
    if (debug)
      fprintf(debug,"Adding %5d J particles for nblist %s\n",NLJ_INC,
	      interaction_function[ftype].longname);

    nlist->maxnrj += NLJ_INC;
    srenew(nlist->jjnr,nlist->maxnrj);
  }

  nri = nlist->nri;
  
  /* Check whether we have to increase the i counter */
  if ((nlist->iinr[nri]  != i_atom) || 
      (nlist->shift[nri] != shift) || 
      (nlist->gid[nri]   != gid)) {
    /* This is something else. Now see if any entries have 
     * been added in the list of the previous atom.
     */
    if ((nlist->jindex[nri+1] > nlist->jindex[nri]) && 
	(nlist->iinr[nri] != -1)) {
      
      /* If so increase the counter */
      if (nlist->nri < nlist->maxnri) {
	nlist->nri++;
	nri++;
      }
      else
	fatal_error(0,"Too many i-atoms for %s (i_atom = %d, maxnri = %d)",
		    interaction_function[ftype].longname,i_atom,nlist->maxnri);
    }
    /* Set the number of neighbours and the atom number */
    nlist->jindex[nri+1] = nlist->jindex[nri];
    nlist->iinr[nri]     = i_atom;
    nlist->gid[nri]      = gid;
    nlist->shift[nri]    = shift;
  }
}

#ifdef SORTNLIST
#define aswap(v,i,j) {  \
  atom_id temp;         \
                        \
  temp=v[i];            \
  v[i]=v[j];            \
  v[j]=temp;            \
}

static void quicksort(atom_id v[], int left, int right)
{
  int i,last;

  if (left >= right)                    /* Do nothing if array contains */
    return;                             /* fewer than two elements      */
  aswap(v,left,(left+right)/2);         /* Move partition element       */
  last=left;                            /* to v[0]                      */
  for(i=left+1; (i<=right); i++)        /* partition                    */
    if (v[i] < v[left]) {
      last++;
      aswap(v,last,i);                  /* watch out for macro trick    */
    }
  aswap(v,left,last);                   /* restore partition element    */
  quicksort(v,left,last-1);
  quicksort(v,last+1,right);
}
#endif

static gmx_inline void close_i_nblist(t_nblist *nlist) 
{
  int nri = nlist->nri;
  
  nlist->jindex[nri+1] = nlist->nrj;
}

static gmx_inline void close_nblist(t_nblist *nlist)
{
  if (nlist->maxnri > 0) {
    int nri = nlist->nri;
    
    if ((nlist->jindex[nri+1] > nlist->jindex[nri]) && 
	(nlist->iinr[nri] != -1)) {
      nlist->nri++;
      nlist->jindex[nri+2] = nlist->nrj;
    }
  }
}

static gmx_inline void close_neighbor_list(t_forcerec *fr,bool bLR,int eNL)
{
  int i;

  if (bLR)
    close_nblist(&(fr->nlist_lr[eNL]));
  else {
    for(i=0; (i<eNL_NR); i++) 
      close_nblist(&(fr->nlist_sr[i]));
  }
}

static void add_j_to_nblist(t_nblist *nlist,int j_atom)
{
  int nrj=nlist->nrj;
  
  nlist->jjnr[nrj] = j_atom;
  nlist->nrj ++;
}

static gmx_inline void put_in_list(bool bHaveLJ[],
				   int nWater,int ngid,t_mdatoms *md,
				   int icg,int jgid,int nj,atom_id jjcg[],
				   atom_id index[],atom_id a[],
				   t_excl bExcl[],int shift,
				   t_forcerec *fr,bool bLR,bool bCoulOnly)
{
  t_nblist  *vdw,*coul,*free=NULL;
  
  int 	    i,j,jcg,igid,gid,ind_ij;
  atom_id   jj,jj0,jj1,i_atom,j_atom;
  int       i0,nicg;
  
  int       *type;
  ushort    *cENER;
  real      *charge;
  real      qi,qq,rlj;
  bool      bWater,bFree,bFreeJ,bNotEx,*bPert;
  
#ifdef SORTNLIST
  /* Quicksort the charge groups in the neighbourlist to obtain
   * better caching properties. We do this only for the short range, 
   * i.e. when we use the nlist more than once
   */
  if (!bLR) 
    quicksort(jjcg,0,nj-1);
#endif
  /* Copy some pointers */
  charge  = md->chargeA;
  type    = md->typeA;
  cENER   = md->cENER;
  bPert   = md->bPerturbed;
  
  /* Check whether this molecule is a water molecule */
  i0     = index[icg];
  nicg   = index[icg+1]-i0;
  bWater = (((type[a[i0]] == nWater) && (nicg == 3)) && 
	    (!bPert[a[i0]]) && (!bPert[a[i0+1]]) && (!bPert[a[i0+2]]));
  if (bWater)
    nicg = 1;
    
  if (bLR) {
    if (bWater) {
      vdw  = &fr->nlist_lr[eNL_VDW_WAT];
      coul = &fr->nlist_lr[eNL_QQ_WAT];
    }
    else {
      vdw  = &fr->nlist_lr[eNL_VDW];
      coul = &fr->nlist_lr[eNL_QQ];
    }
    if (fr->bPert)
      free = &fr->nlist_lr[eNL_FREE];
  }
  else {
    if (bWater) {
      vdw  = &fr->nlist_sr[eNL_VDW_WAT];
      coul = &fr->nlist_sr[eNL_QQ_WAT];
    }
    else {
      vdw  = &fr->nlist_sr[eNL_VDW];
      coul = &fr->nlist_sr[eNL_QQ];
    }
    if (fr->bPert)
      free = &fr->nlist_sr[eNL_FREE];
  }

  /* Loop over the atoms in the i charge group */    
  for(i=0; (i<nicg); i++) {
    i_atom  = a[i0+i];
    igid    = cENER[i_atom];
    gid     = GID(igid,jgid,ngid);
    
    /* Create new i_atom for each energy group */
    if (!bCoulOnly)
      new_i_nblist(vdw,bLR ? F_LJLR : F_LJ,i_atom,shift,gid);
    new_i_nblist(coul,bLR ? F_LR : F_SR,i_atom,shift,gid);
    if (fr->bPert)
      new_i_nblist(free,F_DVDL,i_atom,shift,gid);

    /* Loop over the j charge groups */
    for(j=0; (j<nj); j++) {
      jcg=jjcg[j];
      
      if (bWater && (jcg==icg))
	continue;

      /* Check for large charge groups */
      if (jcg == icg) 
	jj0 = i0 + i + 1;
      else 
	jj0 = index[jcg];
      
      jj1=index[jcg+1];

      if (bWater && !fr->bPert) {
	if (bCoulOnly) {
	  for(jj=jj0; (jj<jj1); jj++) {
	    j_atom = a[jj];
	
	    add_j_to_nblist(coul,j_atom);
	  }
	}
	else {
	  for(jj=jj0; (jj<jj1); jj++) {
	    j_atom = a[jj];
	    
	    if (bHaveLJ[type[j_atom]])
	      add_j_to_nblist(vdw,j_atom);
	    else if (charge[j_atom] != 0)
	      add_j_to_nblist(coul,j_atom);
	  }
	}
      }
      else {
	/* Finally loop over the atoms in the j-charge group */	
	bFree = bPert[i_atom];
	qi    = charge[i_atom];
	for(jj=jj0; (jj<jj1); jj++) {
	  j_atom = a[jj];
	  bFreeJ = bFree || bPert[j_atom];
	  bNotEx = NOTEXCL(bExcl,i,j_atom);
	
	  if (bNotEx) {
	    if (bFreeJ) 
	      add_j_to_nblist(free,j_atom);
	    else if (bCoulOnly) 
	      /* This is done whether or  not bWater is set */
	      add_j_to_nblist(coul,j_atom);
	    else {
	      if (bHaveLJ[type[j_atom]])
		add_j_to_nblist(vdw,j_atom);
	      else if (qi*charge[j_atom] != 0)
		add_j_to_nblist(coul,j_atom);
	    }
	  }
	}
      }
    }
    close_i_nblist(coul);
    close_i_nblist(vdw);
    if (fr->bPert)
      close_i_nblist(free);
  }
}  

void setexcl(int nri,atom_id ia[],t_block *excl,bool b,
	     t_excl bexcl[])
{
  atom_id k;
  int     i,inr;
  
  if (b) {
    for(i=0; (i<nri); i++) {
      inr = ia[i];
      for(k=excl->index[inr]; (k<excl->index[inr+1]); k++) {
	SETEXCL(bexcl,i,excl->a[k]);
      }
    }
  }
  else {
    for(i=0; (i<nri); i++) {
      inr = ia[i];
      for(k=excl->index[inr]; (k<excl->index[inr+1]); k++) {
	RMEXCL(bexcl,i,excl->a[k]);
      }
    }
  }
}

int calc_naaj(int icg,int cgtot)
{
  int naaj;
  
  if ((cgtot % 2) == 1) {
    /* Odd number of charge groups, easy */
    naaj = 1+(cgtot/2);
  }
  else if ((cgtot % 4) == 0) {
    /* Multiple of four is hard */
    if (icg < cgtot/2) {
      if ((icg % 2) == 0)
	naaj=1+(cgtot/2);
      else
	naaj=cgtot/2;
    }
    else {
      if ((icg % 2) == 1)
	naaj=1+(cgtot/2);
      else
	naaj=cgtot/2;
    }
  }
  else {
    /* cgtot/2 = odd */
    if ((icg % 2) == 0)
      naaj=1+(cgtot/2);
    else
      naaj=cgtot/2;
  }
#ifdef DEBUG
  fprintf(log,"naaj=%d\n",naaj);
#endif
  return naaj;
}

/************************************************
 *
 *  S I M P L E      C O R E     S T U F F
 *
 ************************************************/

static real calc_image_rect(rvec xi,rvec xj,rvec box_size,
			    rvec b_inv,int *shift)
{
  const real h15=1.5;
  real ddx,ddy,ddz;
  real dx,dy,dz;
  real r2;
  int  tx,ty,tz;
  
  /* Compute diff vector */
  dx=xj[XX]-xi[XX];
  dy=xj[YY]-xi[YY];
  dz=xj[ZZ]-xi[ZZ];
  
  /* Perform NINT operation, using trunc operation, therefore
   * we first add 1.5 then subtract 1 again
   */
  tx=dx*b_inv[XX]+h15;
  ty=dy*b_inv[YY]+h15;
  tz=dz*b_inv[ZZ]+h15;
  tx--;
  ty--;
  tz--;
  
  /* Correct diff vector for translation */
  ddx=tx*box_size[XX]-dx;
  ddy=ty*box_size[YY]-dy;
  ddz=tz*box_size[ZZ]-dz;
  
  /* Distance squared */
  r2=(ddx*ddx)+(ddy*ddy)+(ddz*ddz);
  
  *shift=XYZ2IS(tx,ty,tz);

  return r2;
}

static void ns_inner_rect(rvec x[],int icg,int njcg,atom_id jcg[],
			  bool bBox,rvec box_size,rvec b_inv,real rcut2,
			  t_block *cgs,t_ns_buf **ns_buf,ushort gid[])
{
  int      shift;
  int      j,nrj,jgid;
  atom_id  cg_j,*cgindex,*cga;
  t_ns_buf *nsbuf;

  cgindex = cgs->index;
  cga     = cgs->a;
  if (bBox) {
    shift = CENTRAL;
    for(j=0; (j<njcg); j++) {
      cg_j   = jcg[j];
      nrj    = cgindex[cg_j+1]-cgindex[cg_j];
      if (calc_image_rect(x[icg],x[cg_j],box_size,b_inv,&shift) < rcut2) {
	jgid  = gid[cga[cgindex[cg_j]]];
	nsbuf = &ns_buf[jgid][shift];
	if (nsbuf->ncg >= MAX_CG) 
	  fatal_error(0,"Too many charge groups (%d) in buffer",nsbuf->ncg);
	nsbuf->jcg[nsbuf->ncg++]=cg_j;
	nsbuf->nj += nrj;
      }
    }
  } else {
    for(j=0; (j<njcg); j++) {
      cg_j   = jcg[j];
      nrj    = cgindex[cg_j+1]-cgindex[cg_j];
      if ((rcut2 == 0) || (distance2(x[icg],x[cg_j]) < rcut2)) {
	jgid  = gid[cga[cgindex[cg_j]]];
	nsbuf = &ns_buf[jgid][CENTRAL];
	if (nsbuf->ncg >= MAX_CG) 
	  fatal_error(0,"Too many charge groups (%d) in buffer",nsbuf->ncg);
	nsbuf->jcg[nsbuf->ncg++]=cg_j;
	nsbuf->nj += nrj;
      }
    }
  }
}

static int ns_simple_core(t_forcerec *fr,
			  t_topology *top,
			  t_mdatoms *md,
			  rvec box_size,t_excl bexcl[],
			  int ngid,t_ns_buf **ns_buf,
			  bool bHaveLJ[])
{
  static   atom_id  *aaj=NULL;
  int      naaj,k;
  real     rlist2;
  int      nsearch,icg,jcg,i0,nri,nn;
  t_ns_buf *nsbuf;
  atom_id  *i_atoms;
  t_block  *cgs=&(top->blocks[ebCGS]);
  t_block  *excl=&(top->atoms.excl);
  rvec     b_inv;
  int      m;
  bool     bBox;
  
  if (aaj==NULL) {
    snew(aaj,2*cgs->nr);
    for(jcg=0; (jcg<cgs->nr); jcg++) {
      aaj[jcg]=jcg;
      aaj[jcg+cgs->nr]=jcg;
    }
  }
  rlist2 = sqr(fr->rlist);

  bBox = (fr->eBox != ebtNONE);
  if (bBox)
    for(m=0; (m<DIM); m++)
      b_inv[m]=divide(1.0,box_size[m]);

  nsearch=0;
  for (icg=fr->cg0; (icg<fr->hcg); icg++) {
    i0      = cgs->index[icg];
    nri     = cgs->index[icg+1]-i0;
    i_atoms = &(cgs->a[i0]);
    setexcl(nri,i_atoms,excl,TRUE,bexcl);
    
    naaj=calc_naaj(icg,cgs->nr);
    ns_inner_rect(fr->cg_cm,icg,naaj,&(aaj[icg]),
		  bBox,box_size,b_inv,rlist2,cgs,ns_buf,md->cENER);
    nsearch += naaj;
    
    for(nn=0; (nn<ngid); nn++) {
      for(k=0; (k<SHIFTS); k++) {
	nsbuf = &(ns_buf[nn][k]);
	if (nsbuf->ncg > 0) { 
	  put_in_list(bHaveLJ,fr->nWater,
		      ngid,md,icg,nn,nsbuf->ncg,nsbuf->jcg,
		      cgs->index,cgs->a,bexcl,k,fr,FALSE,FALSE);
	  nsbuf->ncg=nsbuf->nj=0;
	}
      }
    }
    setexcl(nri,i_atoms,excl,FALSE,bexcl);
  }
  close_neighbor_list(fr,FALSE,-1);
  
  return nsearch;
}

/************************************************
 *
 *    N S 5     G R I D     S T U F F
 *
 ************************************************/
static bool get_dx(int cx,int Nx,int tx,int delta,int *dx0,int *dx1)
{
  if (tx == 1) {
    if (cx >= delta)
      return TRUE;
    else {
      *dx0=Nx+cx-delta;
      *dx1=Nx-1;
    }
  } else if (tx == 0) {
    if (cx < delta)
      *dx0=0;
    else
      *dx0=cx-delta;
    if (cx < Nx-delta)
      *dx1=cx+delta;
    else
      *dx1=Nx-1;
  }
  else {
    if (cx < Nx-delta)
      return TRUE;
    else {
      *dx0=0;
      *dx1=cx+delta-Nx;
    }
  }
 
  return FALSE; 
}

#define sqr(x) ((x)*(x))
#define calc_dx2(XI,YI,ZI,y) (sqr(XI-y[XX])+sqr(YI-y[YY])+sqr(ZI-y[ZZ]))
#define calc_cyl_dx2(XI,YI,y) (sqr(XI-y[XX])+sqr(YI-y[YY]))
/****************************************************
 *
 *    F A S T   N E I G H B O R  S E A R C H I N G
 *
 *    Optimized neighboursearching routine using grid 
 *    of at least 5x5x5, see GROMACS manual
 *
 ****************************************************/

static void do_longrange(FILE *log,t_topology *top,t_forcerec *fr,
			 int ngid,t_mdatoms *md,int icg,
			 int jgid,int nlr,
			 atom_id lr[],t_excl bexcl[],int shift,
			 rvec x[],rvec box_size,t_nrnb *nrnb,
			 real lambda,real *dvdlambda,
			 t_groups *grps,bool bCoulOnly,
			 bool bDoForces,bool bHaveLJ[])
{
  int i;

  for(i=0; (i<eNL_NR); i++) {
    if ((fr->nlist_lr[i].nri > fr->nlist_lr[i].maxnri-32) || bDoForces) {
      close_neighbor_list(fr,TRUE,i);

      /* Evaluate the forces */
      do_fnbf(log,fr,x,fr->flr,md,
	      grps->estat.ee[egLJLR],grps->estat.ee[egLR],box_size,
	      nrnb,lambda,dvdlambda,TRUE,i);
      
      reset_neighbor_list(fr,TRUE,i);
    }
  }

  if (!bDoForces) {  
    /* Put the long range particles in a list */
    put_in_list(bHaveLJ,fr->nWater,
		ngid,md,icg,jgid,nlr,lr,top->blocks[ebCGS].index,
		top->blocks[ebCGS].a,bexcl,shift,fr,TRUE,bCoulOnly);
  }
}

static int ns5_core(FILE *log,t_forcerec *fr,int cg_index[],
		    rvec box_size,int ngid,
		    t_topology *top,t_groups *grps,
		    t_grid *grid,rvec x[],t_excl bexcl[],
		    t_nrnb *nrnb,t_mdatoms *md,
		    real lambda,real *dvdlambda,
		    bool bHaveLJ[])
{
  static atom_id **nl_lr_ljc,**nl_lr_coul,**nl_sr=NULL;
  static int     *nlr_ljc,*nlr_coul,*nsr;
  
  t_block *cgs=&(top->blocks[ebCGS]);
  ushort  *gid=md->cENER;
  atom_id *i_atoms,*cgsatoms=cgs->a,*cgsindex=cgs->index;
  int     tx,ty,tz,cx,cy,cz,dx,dy,dz,cj;
  int     dx0,dx1,dy0,dy1,dz0,dz1;
  int     Nx,Ny,Nz,delta,shift=-1,j,nrj,nns,nn=-1;
  int     icg=-1,iicg,cgsnr,i0,nri,naaj,min_icg,icg_naaj,jjcg,cgj0,jgid;
  int     *grida,*gridnra,*gridind;
  rvec    xi,*cgcm,*svec;
  real    r2,rs2,rvdw2,rcoul2,XI,YI,ZI;
  
  cgsnr    = cgs->nr;
  rs2      = sqr(fr->rlist);
  rvdw2    = sqr(fr->rvdw);
  rcoul2   = sqr(fr->rcoulomb);
  
  if (nl_sr == NULL) {
    /* Short range buffers */
    snew(nl_sr,ngid);
    /* Counters */
    snew(nsr,ngid);
    snew(nlr_ljc,ngid);
    snew(nlr_coul,ngid);
    
    if (rvdw2 > rs2) 
      /* Long range LJC buffers */
      snew(nl_lr_ljc,ngid);
    
    if (rcoul2 > rvdw2) 
      /* Long range Coul buffers */
      snew(nl_lr_coul,ngid);
      
    for(j=0; (j<ngid); j++) {
      snew(nl_sr[j],MAX_CG);
      if (rvdw2 > rs2)
	snew(nl_lr_ljc[j],MAX_CG);
      if (rcoul2 > rvdw2)
	snew(nl_lr_coul[j],MAX_CG);
    }
    if (debug)
      fprintf(debug,"ns5_core: rs2 = %g, rvdw2 = %g, rcoul2 = %g (nm^2)\n",
	      rs2,rvdw2,rcoul2);
  }
  
  /* Unpack arrays */
  cgcm    = fr->cg_cm;
  svec    = fr->shift_vec;
  Nx      = grid->nrx;
  Ny      = grid->nry;
  Nz      = grid->nrz;
  grida   = grid->a;
  gridind = grid->index;
  gridnra = grid->nra;
  delta   = grid->delta;
  nns     = 0;
  
  /* Loop over charge groups */
  for(iicg=fr->cg0; (iicg < fr->hcg); iicg++) {
    icg      = cg_index[iicg];
    /* Consistency check */
    if (icg != iicg)
      fatal_error(0,"icg = %d, iicg = %d, file %s, line %d",icg,iicg,__FILE__,
		  __LINE__);
    i0       = cgsindex[icg];
    nri      = cgsindex[icg+1]-i0;
    i_atoms  = &(cgsatoms[i0]);

    /* Set the exclusions for the atoms in charge group icg using
     * a bitmask
     */    
    setexcl(nri,i_atoms,&top->atoms.excl,TRUE,bexcl);
    
    /* Compute the number of charge groups that fall within the control
     * of this one (icg)
     */
    naaj     = calc_naaj(icg,cgsnr);
    icg_naaj = icg+naaj;
    min_icg  = icg_naaj-cgsnr;
    
    /* Changed iicg to icg, DvdS 990115 
     * (but see consistency check above, DvdS 990330) 
     */
    ci2xyz(grid,icg,&cx,&cy,&cz);
#ifdef NS5DB
    fprintf(log,"icg=%5d, naaj=%5d, cx=%2d, cy=%2d, cz=%2d\n",
	    icg,naaj,cx,cy,cz);
#endif
    /* Loop over shift vectors in three dimensions */
    for (tx=-1; (tx<=1); tx++) {
      /* Calculate range of cells in X direction that have the shift tx */
      if (get_dx(cx,Nx,tx,delta,&dx0,&dx1))
	continue;
      for (ty=-1; (ty<=1); ty++) {
	/* Calculate range of cells in Y direction that have the shift ty */
	if (get_dx(cy,Ny,ty,delta,&dy0,&dy1))
	  continue;
	for (tz=-1; (tz<=1); tz++) {
	  /* Calculate range of cells in Z direction that have the shift tz */
	  if (get_dx(cz,Nz,tz,delta,&dz0,&dz1))
	    continue;

	  /* Get shift vector */	  
	  shift=XYZ2IS(tx,ty,tz);
#ifdef NS5DB
	  assert(shift >= 0);
	  assert(shift < SHIFTS);
#endif
	  XI       = cgcm[icg][XX]+svec[shift][XX];
	  YI       = cgcm[icg][YY]+svec[shift][YY];
	  ZI       = cgcm[icg][ZZ]+svec[shift][ZZ];
	  for(nn=0; (nn<ngid); nn++) {
	    nsr[nn]      = 0;
	    nlr_ljc[nn]  = 0;
	    nlr_coul[nn] = 0;
	  }
#ifdef NS5DB
	  fprintf(log,"shift: %2d, dx0,1: %2d,%2d, dy0,1: %2d,%2d, dz0,1: %2d,%2d\n",
		  shift,dx0,dx1,dy0,dy1,dz0,dz1);
	  fprintf(log,"cgcm: %8.3f  %8.3f  %8.3f\n",cgcm[icg][XX],
		  cgcm[icg][YY],cgcm[icg][ZZ]);
	  fprintf(log,"xi:   %8.3f  %8.3f  %8.3f\n",XI,YI,ZI);
#endif
	  for (dx=dx0; (dx<=dx1); dx++) {
	    for (dy=dy0; (dy<=dy1); dy++) {
	      for (dz=dz0; (dz<=dz1); dz++) {
		/* Find the grid-cell cj in which possible neighbours are */
		cj   = xyz2ci(Ny,Nz,dx,dy,dz);
		
		/* Check out how many cgs (nrj) there in this cell */
		nrj  = gridnra[cj];
		
		/* Find the offset in the cg list */
		cgj0 = gridind[cj];
		
		/* Loop over cgs */
		for (j=0; (j<nrj); j++) {
		  jjcg = grida[cgj0+j];

		  /* check whether this guy is in range! */
		  if (((jjcg >= icg) && (jjcg < icg_naaj)) ||
		      ((jjcg < min_icg))) {
		    r2=calc_dx2(XI,YI,ZI,cgcm[jjcg]);
		    if (r2 < rcoul2) {
		      jgid = gid[cgsatoms[cgsindex[jjcg]]];
		      if (r2 < rs2) {
			if (nsr[jgid] >= MAX_CG) {
			  put_in_list(bHaveLJ,fr->nWater,
				      ngid,md,icg,jgid,nsr[jgid],nl_sr[jgid],
				      cgsindex,cgsatoms,bexcl,
				      shift,fr,FALSE,FALSE);
			  nsr[jgid]=0;
			}
			nl_sr[jgid][nsr[jgid]++]=jjcg;
		      } 
		      else if (r2 < rvdw2) {
			if (nlr_ljc[jgid] >= MAX_CG) {
			  do_longrange(log,top,fr,ngid,md,icg,jgid,
				       nlr_ljc[jgid],
				       nl_lr_ljc[jgid],bexcl,shift,x,
				       box_size,nrnb,
				       lambda,dvdlambda,grps,FALSE,FALSE,
				       bHaveLJ);
			  nlr_ljc[jgid]=0;
			}
			nl_lr_ljc[jgid][nlr_ljc[jgid]++]=jjcg;
		      } 
		      else {
			if (nlr_coul[jgid] >= MAX_CG) {
			  do_longrange(log,top,fr,ngid,md,icg,jgid,
				       nlr_coul[jgid],
				       nl_lr_coul[jgid],bexcl,shift,x,
				       box_size,nrnb,
				       lambda,dvdlambda,grps,TRUE,FALSE,
				       bHaveLJ);
			  nlr_coul[jgid]=0;
			}
			nl_lr_coul[jgid][nlr_coul[jgid]++]=jjcg;
		      }
		    }
		    nns++;
		  }
		}
	      }
	    }
	  }
	  /* CHECK whether there is anything left in the buffers */
	  for(nn=0; (nn<ngid); nn++) {
	    if (nsr[nn] > 0)
	      put_in_list(bHaveLJ,fr->nWater,
			  ngid,md,icg,nn,nsr[nn],nl_sr[nn],
			  cgsindex,cgsatoms,bexcl,shift,fr,FALSE,FALSE);
	    
	    if (nlr_ljc[nn] > 0) 
	      do_longrange(log,top,fr,ngid,md,icg,nn,nlr_ljc[nn],
			   nl_lr_ljc[nn],bexcl,shift,x,box_size,nrnb,
			   lambda,dvdlambda,grps,FALSE,FALSE,bHaveLJ);
	    
	    if (nlr_coul[nn] > 0) 
	      do_longrange(log,top,fr,ngid,md,icg,nn,nlr_coul[nn],
			   nl_lr_coul[nn],bexcl,shift,x,box_size,nrnb,
			   lambda,dvdlambda,grps,TRUE,FALSE,bHaveLJ);
	  }
	}
      }
    }
    setexcl(nri,i_atoms,&top->atoms.excl,FALSE,bexcl);
  }
  /* Perform any left over force calculations */
  if (rvdw2 > rs2)
    do_longrange(log,top,fr,0,md,icg,nn,nlr_ljc[nn],
		 nl_lr_ljc[nn],bexcl,shift,x,box_size,nrnb,
		 lambda,dvdlambda,grps,FALSE,TRUE,bHaveLJ);
  if (rcoul2 > rvdw2)
    do_longrange(log,top,fr,0,md,icg,nn,nlr_coul[nn],
		 nl_lr_coul[nn],bexcl,shift,x,box_size,nrnb,
		 lambda,dvdlambda,grps,TRUE,TRUE,bHaveLJ);
		 
  /* Close off short range neighbourlists */
  close_neighbor_list(fr,FALSE,-1);
  
  return nns;
}

static rvec *sptr;
static int  sdim;

static int  rv_comp(const void *a,const void *b)
{
  int  ia = *(int *)a;
  int  ib = *(int *)b;
  real diff;
  
  diff = sptr[ia][sdim] - sptr[ib][sdim];
  if (diff < 0)
    return -1;
  else if (diff == 0)
    return 0;
  else
    return 1;
}

static void sort_charge_groups(t_commrec *cr,int cg_index[],int slab_index[],
			       rvec cg_cm[],int Dimension)
{
  int i,nrcg,cgind;
  
  nrcg = slab_index[cr->nprocs];
  sptr = cg_cm;
  sdim = Dimension;
  qsort(cg_index,nrcg,sizeof(cg_index[0]),rv_comp);
  
  if (debug) {
    fprintf(debug,"Just sorted the cg_cm array on dimension %d\n",Dimension);
    fprintf(debug,"Index:  Coordinates of cg_cm\n");
    for(i=0; (i<nrcg); i++) {
      cgind = cg_index[i];
      fprintf(debug,"%8d%10.3f%10.3f%10.3f\n",cgind,
	      cg_cm[cgind][XX],cg_cm[cgind][YY],cg_cm[cgind][ZZ]);
    }
  }
  sptr = NULL;
  sdim = -1;
}

int search_neighbours(FILE *log,t_forcerec *fr,
		      rvec x[],matrix box,
		      t_topology *top,t_groups *grps,
		      t_commrec *cr,t_nsborder *nsb,
		      t_nrnb *nrnb,t_mdatoms *md,
		      real lambda,real *dvdlambda)
{
  static   bool        bFirst=TRUE;
  static   t_grid      *grid=NULL;
  static   t_excl      *bexcl;
  static   bool        *bHaveLJ;
  static   t_ns_buf    **ns_buf=NULL;
  static   int         *cg_index=NULL,*slab_index=NULL;
  static   bool        bSwitched=FALSE;
  
  t_block  *cgs=&(top->blocks[ebCGS]);
  rvec     box_size;
  int      i,j,m,ngid;

  int      nsearch;
  bool     bGrid;
  char     *ptr;
  
  /* Set some local variables */
  bGrid=fr->bGrid;
  ngid=top->atoms.grps[egcENER].nr;
  
  for(m=0; (m<DIM); m++)
    box_size[m]=box[m][m];

  /* First time initiation of arrays etc. */  
  if (bFirst) {
    int icg,nr_in_cg,maxcg;
    
    /* Compute largest charge groups size (# atoms) */
    nr_in_cg=1;
    for (icg=0; (icg < cgs->nr); icg++)
      nr_in_cg=max(nr_in_cg,(int)(cgs->index[icg+1]-cgs->index[icg]));

    /* Verify whether largest charge group is <= max cg.
     * This is determined by the type of the local exclusion type 
     * Exclusions are stored in bits. (If the type is not large
     * enough, enlarge it, unsigned char -> unsigned short -> unsigned long)
     */
    maxcg=sizeof(t_excl)*8;
    if (nr_in_cg > maxcg)
      fatal_error(0,"Max #atoms in a charge group: %d > %d\n",
		  nr_in_cg,maxcg);
      
    snew(bexcl,cgs->nra);
    debug_gmx();

    if ((ptr=getenv("NLIST")) != NULL) {
      sscanf(ptr,"%d",&NLJ_INC);
      
      fprintf(log,"%s: I will increment J-lists by %d\n",
	      __FILE__,NLJ_INC);
    }

    /* Check whether we have to do domain decomposition,
     * if so set local variables for the charge group index and the
     * slab index.
     */
    if (fr->bDomDecomp) {
      snew(slab_index,cr->nprocs+1);
      for(i=0; (i<=cr->nprocs); i++)
	slab_index[i] = i*((real) cgs->nr/((real) cr->nprocs));
      fr->cg0 = slab_index[cr->pid];
      fr->hcg = slab_index[cr->pid+1] - fr->cg0;
      if (debug)
	fprintf(debug,"Will use DOMAIN DECOMPOSITION, "
		"from charge group index %d to %d on cpu %d\n",
		fr->cg0,fr->cg0+fr->hcg,cr->pid);
    }
    snew(cg_index,cgs->nr+1);
    for(i=0; (i<=cgs->nr);  i++)
      cg_index[i] = i;
    
    if (bGrid) {
      snew(grid,1);
      init_grid(log,grid,fr->ndelta,box,fr->rlistlong,cgs->nr);
    }
    
    /* Create array that determines whether or not atoms have LJ */
    snew(bHaveLJ,fr->ntype);
    for(i=0; (i<fr->ntype); i++) {
      for(j=0; (j<fr->ntype); j++) {
	bHaveLJ[i] = (bHaveLJ[i] || 
		      (fr->bBHAM ? 
		       ((BHAMA(fr->nbfp,fr->ntype,i,j) != 0) ||
			(BHAMB(fr->nbfp,fr->ntype,i,j) != 0) ||
			(BHAMC(fr->nbfp,fr->ntype,i,j) != 0)) :
		       ((C6(fr->nbfp,fr->ntype,i,j) != 0) ||
			(C12(fr->nbfp,fr->ntype,i,j) != 0))));
      }
    }
    if (debug) 
      pr_ivec(debug,0,"bHaveLJ",bHaveLJ,fr->ntype);
    
    bFirst=FALSE;
  }
  debug_gmx();
  
  /* Reset the neighbourlists */
  reset_neighbor_list(fr,FALSE,-1);
  
  if (bGrid) {
    grid_first(log,grid,box,fr->rlistlong);
    /* Check if box is big enough to do grid searching... */
    if ( !( (grid->nrx >= 2*grid->delta+1) && 
	    (grid->nry >= 2*grid->delta+1) && 
	    (grid->nrz >= 2*grid->delta+1) ) ) {
      if (!bSwitched)
	fprintf(log,"WARNING: Box too small for grid-search, "
		"switching to simple neighboursearch.\n");
      if (fr->bTwinRange)
	fatal_error(0,"TWIN-RANGE cut-off with Simple "
		    "neighboursearching not implemented.\n"
		    "Use grid neighboursearching, and make (rlong < 0.4 box)");
      bGrid=FALSE;
      bSwitched=TRUE;
    } else {
      if (bSwitched)
	fprintf(log,"WARNING: Box large enough again for grid-search\n");
      bSwitched=FALSE;
    }
  }
  debug_gmx();
  
  if (bGrid) {
    /* Don't know why this all is... (DvdS 3/99) */
#ifndef SEGV
    int start = 0;
    int end   = cgs->nr;
#else
    int start = fr->cg0;
    int end   = (cgs->nr+1)/2;
#endif

    if (fr->bDomDecomp)
      sort_charge_groups(cr,cg_index,slab_index,fr->cg_cm,fr->Dimension);
    debug_gmx();
    
    fill_grid(log,fr->bDomDecomp,cg_index,
	      grid,box,cgs->nr,fr->cg0,fr->hcg,fr->cg_cm);
    debug_gmx();

    if (PAR(cr))
      mv_grid(cr,fr->bDomDecomp,cg_index,grid,nsb->workload);
    debug_gmx();
      
    calc_elemnr(log,fr->bDomDecomp,cg_index,grid,start,end,cgs->nr);
    calc_ptrs(grid);
    grid_last(log,fr->bDomDecomp,cg_index,grid,start,end,cgs->nr);

    if (debug) {
      check_grid(debug,grid);
      print_grid(debug,grid,fr->bDomDecomp,cg_index);
    }
  }
  debug_gmx();
  
  /* Do the core! */
  if (bGrid)
    nsearch = ns5_core(log,fr,cg_index,box_size,ngid,top,grps,
		       grid,x,bexcl,nrnb,md,lambda,dvdlambda,bHaveLJ);
  else {
    /* Only allocate this when necessary, saves 100 kb */
    if (ns_buf == NULL) {
      snew(ns_buf,ngid);
      for(i=0; (i<ngid); i++)
	snew(ns_buf[i],SHIFTS);
    }
    nsearch = ns_simple_core(fr,top,md,box_size,
			     bexcl,ngid,ns_buf,bHaveLJ);
  }
  debug_gmx();
  
#ifdef DEBUG
  pr_nsblock(log);
#endif

  inc_nrnb(nrnb,eNR_NS,nsearch);
  /* inc_nrnb(nrnb,eNR_LR,fr->nlr); */

  return nsearch;
}


