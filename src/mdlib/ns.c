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

static int NLI_INC = 1000;
static int NLJ_INC = 16384;

static void reallocate_nblist(t_nblist *nl)
{
  if (debug)
    fprintf(debug,"reallocating neigborlist il_code=%d, maxnri=%d\n",
	    nl->il_code,nl->maxnri); 
  srenew(nl->iinr,   nl->maxnri+2);
  srenew(nl->gid,    nl->maxnri+2);
  srenew(nl->shift,  nl->maxnri+2);
  srenew(nl->jindex, nl->maxnri+2);
}

static void init_nblist(t_nblist *nl,int homenr,int il_code)
{
  nl->il_code = il_code;
  /* maxnri is influenced by the number of shifts (maximum is 8)
   * and the number of energy groups.
   * If it is not enough, nl memory will be reallocated during the run.
   * 4 seems to be a reasonable factor, which only causes reallocation
   * during runs with tiny and many energygroups.
   */
  nl->maxnri  = homenr*4;
  nl->maxnrj  = 0;
  nl->nri     = 0;
  nl->nrj     = 0;
  nl->iinr    = NULL;
  nl->gid     = NULL;
  nl->shift   = NULL;
  nl->jindex  = NULL;
  reallocate_nblist(nl);
  nl->jindex[0] = 0;
  nl->jindex[1] = 0;
  if (nl->maxnri > 0)
    nl->iinr[0] = -1;
}

static unsigned int nbf_index(bool bVDWOnly,bool bCoulOnly,bool bRF,bool bBham,
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
  
  ni = bCoulOnly | (bBham << 1) | (bRF << 2) | (bTab << 3) | (bWater << 4)
      | (bEwald << 5);
  
  return inloop[ni];
}

void init_neighbor_list(FILE *log,t_forcerec *fr,int homenr)
{
  /* Make maxlr tunable! (does not seem to be a big difference though) 
   * This parameter determines the number of i particles in a long range 
   * neighbourlist. Too few means many function calls, too many means
   * cache trashing.
   */
  int maxsr,maxsr_wat,maxlr,maxlr_wat;
  
  maxsr     = homenr-fr->nWatMol*3;
  if (maxsr < 0)
    fatal_error(0,"%s, %d: Negative number of short range atoms.\n"
		"Call your Gromacs dealer for assistance.",__FILE__,__LINE__);
  maxsr_wat = fr->nWatMol; 
  maxlr     = max(maxsr/10,50);
  maxlr_wat = max(maxsr_wat/10,50);

  init_nblist(&fr->nlist_sr[eNL_VDWQQ],maxsr,
	   nbf_index(FALSE,FALSE,fr->bRF,fr->bBHAM,fr->bTab,
		     FALSE,fr->bEwald));
  init_nblist(&fr->nlist_sr[eNL_VDW],maxsr,
	   nbf_index(TRUE,FALSE,fr->bRF,fr->bBHAM,fr->bTab,
		     FALSE,fr->bEwald));
  init_nblist(&fr->nlist_sr[eNL_QQ],maxsr,
	   nbf_index(FALSE,TRUE,fr->bRF,fr->bBHAM,fr->bTab,
		     FALSE,fr->bEwald));
  if (fr->efep != efepNO)
    init_nblist(&fr->nlist_sr[eNL_FREE],maxsr,
	     fr->bBHAM ? eNR_BHAM_FREE : eNR_LJC_FREE);
  if (fr->bWaterOpt) {
    init_nblist(&fr->nlist_sr[eNL_VDWQQ_WAT],maxsr_wat,
	     nbf_index(FALSE,FALSE,fr->bRF,fr->bBHAM,fr->bTab,
		       TRUE,fr->bEwald));
    init_nblist(&fr->nlist_sr[eNL_QQ_WAT],maxsr_wat,
	     nbf_index(FALSE,TRUE,fr->bRF,fr->bBHAM,fr->bTab,
		       TRUE,fr->bEwald));
  }
  if (fr->bTwinRange) {
    fprintf(log,"Allocating space for long range neighbour list of %d atoms\n",
	    maxlr);
    init_nblist(&fr->nlist_lr[eNL_VDWQQ],maxlr,
	     nbf_index(FALSE,FALSE,fr->bRF,fr->bBHAM,fr->bTab,
		       FALSE,fr->bEwald));
    init_nblist(&fr->nlist_lr[eNL_VDW],maxlr,
	     nbf_index(TRUE,FALSE,fr->bRF,fr->bBHAM,fr->bTab,
		       FALSE,fr->bEwald));
    init_nblist(&fr->nlist_lr[eNL_QQ],maxlr,
	     nbf_index(FALSE,TRUE,fr->bRF,fr->bBHAM,fr->bTab,
		       FALSE,fr->bEwald));
    if (fr->efep != efepNO)
      init_nblist(&fr->nlist_lr[eNL_FREE],maxlr,
	       fr->bBHAM ? eNR_BHAM_FREE : eNR_LJC_FREE);
    if (fr->bWaterOpt) {
      init_nblist(&fr->nlist_lr[eNL_VDWQQ_WAT],maxlr_wat,
	       nbf_index(FALSE,FALSE,fr->bRF,fr->bBHAM,fr->bTab,
			 TRUE,fr->bEwald));
      init_nblist(&fr->nlist_lr[eNL_QQ_WAT],maxlr_wat,
	       nbf_index(FALSE,TRUE,fr->bRF,fr->bBHAM,fr->bTab,
			 TRUE,fr->bEwald));
    }
  }
}

static void reset_nblist(t_nblist *nl)
{
  nl->nri       = 0;
  nl->nrj       = 0;
  if (nl->maxnri > 0) {
    nl->iinr[0]   = -1;
    if (nl->maxnrj > 1) {
      nl->jindex[0] = 0;
      nl->jindex[1] = 0;
    }
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
      nlist->nri++;
      nri++;
      if (nlist->nri >= nlist->maxnri) {
	nlist->maxnri += NLI_INC;
	reallocate_nblist(nlist);
      }
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
				   int ngid,t_mdatoms *md,
				   int icg,int jgid,int nj,atom_id jjcg[],
				   atom_id index[],atom_id a[],
				   t_excl bExcl[],int shift,
				   t_forcerec *fr,bool bLR,
				   bool bVDWOnly,bool bCoulOnly)
{
  /* The a[] index has been removed,
     to put it back in i_atom should be a[i0] and jj should be a[jj].
     */
  t_nblist  *vdwc,*vdw,*coul,*free=NULL;
  
  int 	    i,j,jcg,igid,gid,ind_ij;
  atom_id   jj,jj0,jj1,i_atom;
  int       i0,nicg;
  
  int       *type;
  unsigned short    *cENER;
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
  bWater = (((type[i0] == fr->nWater) && (nicg == 3)) && 
	    (!bPert[i0]) && (!bPert[i0+1]) && (!bPert[i0+2]));
  if (bWater && fr->efep==efepNO)
    nicg = 1;
    
  /* Unpack pointers to neighbourlist structs */
  if (bLR) {
    /* Long range */
    if (bWater) {
      vdwc = &fr->nlist_lr[eNL_VDWQQ_WAT];
      vdw  = &fr->nlist_lr[eNL_VDW];
      coul = &fr->nlist_lr[eNL_QQ_WAT];
    }
    else {
      vdwc = &fr->nlist_lr[eNL_VDWQQ];
      vdw  = &fr->nlist_lr[eNL_VDW];
      coul = &fr->nlist_lr[eNL_QQ];
    }
    if (fr->efep != efepNO)
      free = &fr->nlist_lr[eNL_FREE];
  }
  else {
    /* Short range */
    if (bWater) {
      vdwc = &fr->nlist_sr[eNL_VDWQQ_WAT];
      vdw  = &fr->nlist_sr[eNL_VDW];
      coul = &fr->nlist_sr[eNL_QQ_WAT];
    }
    else {
      vdwc = &fr->nlist_sr[eNL_VDWQQ];
      vdw  = &fr->nlist_sr[eNL_VDW];
      coul = &fr->nlist_sr[eNL_QQ];
    }
    if (fr->efep != efepNO)
      free = &fr->nlist_sr[eNL_FREE];
  }

  if (bWater && fr->efep==efepNO) {
    /* Loop over the atoms in the i charge group */    
    for(i=0; (i<nicg); i++) {
      i_atom  = i0;
      igid    = cENER[i_atom];
      gid     = GID(igid,jgid,ngid);
      
      /* Create new i_atom for each energy group */
      if (!bCoulOnly && !bVDWOnly)
	new_i_nblist(vdwc,bLR ? F_LJLR : F_LJ,i_atom,shift,gid);
      if (!bCoulOnly)
	new_i_nblist(vdw,bLR ? F_LJLR : F_LJ,i_atom,shift,gid);
      if (!bVDWOnly)
	new_i_nblist(coul,bLR ? F_LR : F_SR,i_atom,shift,gid);
      
      /* Loop over the j charge groups */
      for(j=0; (j<nj); j++) {
	jcg=jjcg[j];
	
	if (jcg==icg)
	  continue;
	
	jj0 = index[jcg];
	if (type[jj0] == fr->nWater) {
	  if (bVDWOnly)
	    add_j_to_nblist(vdw,jj0);
	  else {
	    if (bCoulOnly)
	      add_j_to_nblist(coul,jj0);
	    else
	      add_j_to_nblist(vdwc,jj0);
	    add_j_to_nblist(coul,jj0+1);
	    add_j_to_nblist(coul,jj0+2);
	  }
	} else {
	  jj1 = index[jcg+1];
	  
	  if (bCoulOnly) {
	    for(jj=jj0; (jj<jj1); jj++) {
	      if (charge[jj] != 0)
		add_j_to_nblist(coul,jj);
	    }
	  } else if (bVDWOnly) {
	    for(jj=jj0; (jj<jj1); jj++) 
	      if (bHaveLJ[type[jj]])
		add_j_to_nblist(vdw,jj);
	  } else {
	    for(jj=jj0; (jj<jj1); jj++) {
	      if (bHaveLJ[type[jj]]) {
		if (charge[jj] != 0)
		  add_j_to_nblist(vdwc,jj);
		else
		  add_j_to_nblist(vdw,jj);
	      } else if (charge[jj] != 0)
		add_j_to_nblist(coul,jj);
	    }
	  }
	}
      }
      close_i_nblist(vdw);
      close_i_nblist(coul);
      close_i_nblist(vdwc);
    }
  } else {
    /* Loop over the atoms in the i charge group */    
    for(i=0; (i<nicg); i++) {
      i_atom  = i0+i;
      igid    = cENER[i_atom];
      gid     = GID(igid,jgid,ngid);
      qi      = charge[i_atom];
      
      /* Create new i_atom for each energy group */
      if (!bCoulOnly && !bVDWOnly)
	new_i_nblist(vdwc,bLR ? F_LJLR : F_LJ,i_atom,shift,gid);
      if (!bCoulOnly)
	new_i_nblist(vdw,bLR ? F_LR : F_SR,i_atom,shift,gid);
      if (!bVDWOnly)
	new_i_nblist(coul,bLR ? F_LR : F_SR,i_atom,shift,gid);
      if (fr->efep != efepNO)
	new_i_nblist(free,F_DVDL,i_atom,shift,gid);
      
      if (!(bVDWOnly || qi==0) || !(bCoulOnly || !bHaveLJ[type[i_atom]])) {
	/* Loop over the j charge groups */
	for(j=0; (j<nj); j++) {
	  jcg=jjcg[j];
	  
	  /* Check for large charge groups */
	  if (jcg == icg) 
	    jj0 = i0 + i + 1;
	  else 
	    jj0 = index[jcg];
	  
	  jj1=index[jcg+1];
	  
	  /* Finally loop over the atoms in the j-charge group */	
	  bFree = bPert[i_atom];
	  for(jj=jj0; (jj<jj1); jj++) {
	    bFreeJ = bFree || bPert[jj];
	    /* Complicated if, because the water H's should also
	     * see perturbed j-particles
	     */
	    if (!bWater || i==0 || bFreeJ) {
	      bNotEx = NOTEXCL(bExcl,i,jj);
	      
	      if (bNotEx) {
		if (bFreeJ) 
		    add_j_to_nblist(free,jj);
		else if (bCoulOnly) { 
		  /* This is done whether or  not bWater is set */
		  if (charge[jj] != 0)
		    add_j_to_nblist(coul,jj);
		} else if (bVDWOnly) { 
		  if (bHaveLJ[type[jj]])
		    add_j_to_nblist(vdw,jj);
		} else {
		  if (bHaveLJ[type[jj]]) {
		    if (qi != 0 && charge[jj] != 0)
		      add_j_to_nblist(vdwc,jj);
		    else
		      add_j_to_nblist(vdw,jj);
		  } else if (qi != 0 && charge[jj] != 0)
		    add_j_to_nblist(coul,jj);
		}
	      }
	    }
	  }
	}
      }
      close_i_nblist(vdw);
      close_i_nblist(coul);
      close_i_nblist(vdwc);
      if (fr->efep != efepNO)
	close_i_nblist(free);
    }
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

static real calc_image_tric(rvec xi,rvec xj,matrix box,
			    rvec b_inv,int *shift)
{
  /* This code assumes that the cut-off is smaller than
   * a half times the smallest diagonal element of the box.
   */
  const real h15=1.5;
  real dx,dy,dz;
  real r2;
  int  tx,ty,tz;
  
  /* Compute diff vector */
  dz=xj[ZZ]-xi[ZZ];
  dy=xj[YY]-xi[YY];
  dx=xj[XX]-xi[XX];
  
  /* Perform NINT operation, using trunc operation, therefore
   * we first add 1.5 then subtract 1 again
   */
  tz=dz*b_inv[ZZ]+h15;
  tz--;
  dz-=tz*box[ZZ][ZZ];
  dy-=tz*box[ZZ][YY];
  dx-=tz*box[ZZ][XX];

  ty=dy*b_inv[YY]+h15;
  ty--;
  dy-=ty*box[YY][YY];
  dx-=ty*box[YY][XX];

  tx=dx*b_inv[XX]+h15;
  tx--;
  dx-=tx*box[XX][XX];
  
  /* Distance squared */
  r2=(dx*dx)+(dy*dy)+(dz*dz);
  
  *shift=XYZ2IS(tx,ty,tz);

  return r2;
}

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

static void add_simple(t_ns_buf *nsbuf,int nrj,atom_id cg_j,
		       bool bHaveLJ[],int ngid,t_mdatoms *md,
		       int icg,int jgid,t_block *cgs,t_excl bexcl[],
		       int shift,t_forcerec *fr)
{
  if (nsbuf->ncg >= MAX_CG) {
    put_in_list(bHaveLJ,ngid,md,icg,jgid,nsbuf->ncg,nsbuf->jcg,
		cgs->index,cgs->a,bexcl,shift,fr,FALSE,FALSE,FALSE);
    /* Reset buffer contents */
    nsbuf->ncg = nsbuf->nj = 0;
  }
  nsbuf->jcg[nsbuf->ncg++]=cg_j;
  nsbuf->nj += nrj;
}

static void ns_inner_tric(rvec x[],int icg,bool *i_eg_excl,
			  int njcg,atom_id jcg[],
			  matrix box,rvec b_inv,real rcut2,
			  t_block *cgs,t_ns_buf **ns_buf,ushort gid[],
			  bool bHaveLJ[],int ngid,t_mdatoms *md,
			  t_excl bexcl[],t_forcerec *fr)
{
  int      shift;
  int      j,nrj,jgid;
  atom_id  cg_j,*cgindex /* ,*cga */;
  t_ns_buf *nsbuf;
  
  cgindex = cgs->index;
  /* cga     = cgs->a; */
  shift   = CENTRAL;
  for(j=0; (j<njcg); j++) {
    cg_j   = jcg[j];
    nrj    = cgindex[cg_j+1]-cgindex[cg_j];
    if (calc_image_tric(x[icg],x[cg_j],box,b_inv,&shift) < rcut2) {
      /* jgid  = gid[cga[cgindex[cg_j]]]; */
      jgid  = gid[cgindex[cg_j]];
      if (!i_eg_excl[jgid]) {
	add_simple(&ns_buf[jgid][shift],nrj,cg_j,
		   bHaveLJ,ngid,md,icg,jgid,cgs,bexcl,shift,fr);
      }
    }
  }
}

static void ns_inner_rect(rvec x[],int icg,bool *i_eg_excl,
			  int njcg,atom_id jcg[],
			  bool bBox,rvec box_size,rvec b_inv,real rcut2,
			  t_block *cgs,t_ns_buf **ns_buf,unsigned short gid[],
			  bool bHaveLJ[],int ngid,t_mdatoms *md,
			  t_excl bexcl[],t_forcerec *fr)
{
  int      shift;
  int      j,nrj,jgid;
  atom_id  cg_j,*cgindex /* ,*cga */;
  t_ns_buf *nsbuf;

  cgindex = cgs->index;
  /* cga     = cgs->a; */
  if (bBox) {
    shift = CENTRAL;
    for(j=0; (j<njcg); j++) {
      cg_j   = jcg[j];
      nrj    = cgindex[cg_j+1]-cgindex[cg_j];
      if (calc_image_rect(x[icg],x[cg_j],box_size,b_inv,&shift) < rcut2) {
	/* jgid  = gid[cga[cgindex[cg_j]]]; */
	jgid  = gid[cgindex[cg_j]];
	if (!i_eg_excl[jgid]) {
	  add_simple(&ns_buf[jgid][shift],nrj,cg_j,
		     bHaveLJ,ngid,md,icg,jgid,cgs,bexcl,shift,fr);
	}
      }
    }
  } 
  else {
    for(j=0; (j<njcg); j++) {
      cg_j   = jcg[j];
      nrj    = cgindex[cg_j+1]-cgindex[cg_j];
      if ((rcut2 == 0) || (distance2(x[icg],x[cg_j]) < rcut2)) {
	/* jgid  = gid[cga[cgindex[cg_j]]]; */
	jgid  = gid[cgindex[cg_j]];
	if (!i_eg_excl[jgid]) {
	  add_simple(&ns_buf[jgid][CENTRAL],nrj,cg_j,
		     bHaveLJ,ngid,md,icg,jgid,cgs,bexcl,CENTRAL,fr);
	}
      }
    }
  }
}

static int ns_simple_core(t_forcerec *fr,
			  t_topology *top,
			  t_mdatoms *md,
			  matrix box,rvec box_size,
			  t_excl bexcl[],
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
  bool     bBox,bTriclinic,*i_eg_excl;
  
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
  bTriclinic = TRICLINIC(box);

  nsearch=0;
  for (icg=fr->cg0; (icg<fr->hcg); icg++) {
    i0        = cgs->index[icg];
    nri       = cgs->index[icg+1]-i0;
    i_atoms   = &(cgs->a[i0]);
    i_eg_excl = fr->eg_excl + ngid*md->cENER[*i_atoms];
    setexcl(nri,i_atoms,excl,TRUE,bexcl);
    
    naaj=calc_naaj(icg,cgs->nr);
    if (bTriclinic)
      ns_inner_tric(fr->cg_cm,icg,i_eg_excl,naaj,&(aaj[icg]),
		    box,b_inv,rlist2,cgs,ns_buf,md->cENER,
		    bHaveLJ,ngid,md,bexcl,fr);
    else
      ns_inner_rect(fr->cg_cm,icg,i_eg_excl,naaj,&(aaj[icg]),
		    bBox,box_size,b_inv,rlist2,cgs,ns_buf,md->cENER,
		    bHaveLJ,ngid,md,bexcl,fr);
    nsearch += naaj;
    
    for(nn=0; (nn<ngid); nn++) {
      for(k=0; (k<SHIFTS); k++) {
	nsbuf = &(ns_buf[nn][k]);
	if (nsbuf->ncg > 0) {
	  put_in_list(bHaveLJ,ngid,md,icg,nn,nsbuf->ncg,nsbuf->jcg,
		      cgs->index,cgs->a,bexcl,k,fr,FALSE,FALSE,FALSE);
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

static gmx_inline void get_dx(int Nx,real gridx,real grid_x,real rc2,real x,
			      int *dx0,int *dx1,real *dcx2)
{
  real dcx,tmp;
  int  xgi,xgi0,xgi1,i;

  xgi = (int)(Nx+x*grid_x)-Nx;
  if (xgi < 0) {
    *dx0 = 0;
    xgi0 = -1;
    *dx1 = -1;
    xgi1 = 0;
  } else if (xgi >= Nx) {
    *dx0 = Nx;
    xgi0 = Nx-1;
    *dx1 = Nx-1;
    xgi1 = Nx;
  } else {
    dcx2[xgi] = 0;
    *dx0 = xgi;
    xgi0 = xgi-1;
    *dx1 = xgi;
    xgi1 = xgi+1;
  }

  for(i=xgi0; i>=0; i--) {
     dcx = (i+1)*gridx-x;
     tmp = dcx*dcx;
     if (tmp >= rc2)
       break;
     *dx0 = i;
     dcx2[i] = tmp;
  }
  for(i=xgi1; i<Nx; i++) {
     dcx = i*gridx-x;
     tmp = dcx*dcx;
     if (tmp >= rc2)
       break;
     *dx1 = i;
     dcx2[i] = tmp;
  }
}

#define sqr(x) ((x)*(x))
#define calc_dx2(XI,YI,ZI,y) (sqr(XI-y[XX])+sqr(YI-y[YY])+sqr(ZI-y[ZZ]))
#define calc_cyl_dx2(XI,YI,y) (sqr(XI-y[XX])+sqr(YI-y[YY]))
/****************************************************
 *
 *    F A S T   N E I G H B O R  S E A R C H I N G
 *
 *    Optimized neighboursearching routine using grid 
 *    at least 1x1x1, see GROMACS manual
 *
 ****************************************************/

static void do_longrange(FILE *log,t_topology *top,t_forcerec *fr,
			 int ngid,t_mdatoms *md,int icg,
			 int jgid,int nlr,
			 atom_id lr[],t_excl bexcl[],int shift,
			 rvec x[],rvec box_size,t_nrnb *nrnb,
			 real lambda,real *dvdlambda,
			 t_groups *grps,int bVDWOnly,bool bCoulOnly,
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
    put_in_list(bHaveLJ,ngid,md,icg,jgid,nlr,lr,top->blocks[ebCGS].index,
		top->blocks[ebCGS].a,bexcl,shift,fr,
		TRUE,bVDWOnly,bCoulOnly);
  }
}

static int ns5_core(FILE *log,t_forcerec *fr,int cg_index[],
		    matrix box,rvec box_size,int ngid,
		    t_topology *top,t_groups *grps,
		    t_grid *grid,rvec x[],t_excl bexcl[],
		    t_nrnb *nrnb,t_mdatoms *md,
		    real lambda,real *dvdlambda,
		    bool bHaveLJ[])
{
  static atom_id **nl_lr_ljc,**nl_lr_one,**nl_sr=NULL;
  static int     *nlr_ljc,*nlr_one,*nsr;
  static real *dcx2=NULL,*dcy2=NULL,*dcz2=NULL;
  
  t_block *cgs=&(top->blocks[ebCGS]);
  unsigned short  *gid=md->cENER;
  atom_id *i_atoms,*cgsatoms=cgs->a,*cgsindex=cgs->index;
  int     tx,ty,tz,dx,dy,dz,cj;
  int     dx0,dx1,dy0,dy1,dz0,dz1;
  int     Nx,Ny,Nz,shift=-1,j,nrj,nns,nn=-1;
  real    gridx,gridy,gridz,grid_x,grid_y,grid_z;
  real    margin_x,margin_y,sh_zx,sh_zy,sh_x;
  int     icg=-1,iicg,cgsnr,i0,nri,naaj,min_icg,icg_naaj,jjcg,cgj0,jgid;
  int     *grida,*gridnra,*gridind;
  bool    bVDWOnly,bCoulOnly;
  rvec    xi,*cgcm;
  real    r2,rs2,rvdw2,rcoul2,rm2,rl2,XI,YI,ZI,dcx,dcy,dcz,tmp1,tmp2;
  bool    *i_eg_excl;
  
  cgsnr    = cgs->nr;
  rs2      = sqr(fr->rlist);
  if (fr->bTwinRange) {
    rvdw2  = sqr(fr->rvdw);
    rcoul2 = sqr(fr->rcoulomb);
  } else {
    /* Workaround for a gcc -O3 or -ffast-math problem */
    rvdw2  = rs2;
    rcoul2 = rs2;
  }
  rm2 = min(rvdw2,rcoul2);
  rl2 = max(rvdw2,rcoul2);
  if (rvdw2 > rcoul2) {
    bVDWOnly  = TRUE;
    bCoulOnly = FALSE;
  } else {
    bVDWOnly  = FALSE;
    bCoulOnly = TRUE;
  }

  if (nl_sr == NULL) {
    /* Short range buffers */
    snew(nl_sr,ngid);
    /* Counters */
    snew(nsr,ngid);
    snew(nlr_ljc,ngid);
    snew(nlr_one,ngid);

    if (rm2 > rs2) 
      /* Long range VdW and Coul buffers */
      snew(nl_lr_ljc,ngid);
    
    if (rl2 > rm2)
      /* Long range VdW or Coul only buffers */
      snew(nl_lr_one,ngid);
    
    for(j=0; (j<ngid); j++) {
      snew(nl_sr[j],MAX_CG);
      if (rm2 > rs2)
	snew(nl_lr_ljc[j],MAX_CG);
      if (rl2 > rm2)
	snew(nl_lr_one[j],MAX_CG);
    }
    if (debug)
      fprintf(debug,"ns5_core: rs2 = %g, rvdw2 = %g, rcoul2 = %g (nm^2)\n",
	      rs2,rvdw2,rcoul2);
  }

  /* Unpack arrays */
  cgcm    = fr->cg_cm;
  Nx      = grid->nrx;
  Ny      = grid->nry;
  Nz      = grid->nrz;
  grida   = grid->a;
  gridind = grid->index;
  gridnra = grid->nra;
  nns     = 0;

  if (dcx2 == NULL) {
    /* Allocate tmp arrays */
    snew(dcx2,Nx*2);
    snew(dcy2,Ny*2);
    snew(dcz2,Nz*2);
  }

  gridx      = box[XX][XX]/grid->nrx;
  gridy      = box[YY][YY]/grid->nry;
  gridz      = box[ZZ][ZZ]/grid->nrz;
  grid_x     = 1/gridx;
  grid_y     = 1/gridy;
  grid_z     = 1/gridz;

  debug_gmx();
  
  /* Loop over charge groups */
  for(iicg=fr->cg0; (iicg < fr->hcg); iicg++) {
    icg      = cg_index[iicg];
    /* Consistency check */
    if (icg != iicg)
      fatal_error(0,"icg = %d, iicg = %d, file %s, line %d",icg,iicg,__FILE__,
		  __LINE__);
    i0        = cgsindex[icg];
    nri       = cgsindex[icg+1]-i0;
    i_atoms   = &(cgsatoms[i0]);
    i_eg_excl = fr->eg_excl + ngid*gid[*i_atoms];
    
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
#ifdef NS5DB
    fprintf(log,"icg=%5d, naaj=%5d\n",icg,naaj);
#endif
    /* Loop over shift vectors in three dimensions */
    for (tz=-1; tz<=1; tz++) {
      ZI = cgcm[icg][ZZ]+tz*box[ZZ][ZZ];
      /* Calculate range of cells in Z direction that have the shift tz */
      get_dx(Nz,gridz,grid_z,rcoul2,ZI,&dz0,&dz1,dcz2);
      if (dz0 > dz1)
	continue;
      for (ty=-1; ty<=1; ty++) {
	YI = cgcm[icg][YY]+ty*box[YY][YY]+tz*box[ZZ][YY];
	/* Calculate range of cells in Y direction that have the shift ty */
	get_dx(Ny,gridy,grid_y,rcoul2,YI,&dy0,&dy1,dcy2);
	if (dy0 > dy1)
	  continue;
	for (tx=-1; tx<=1; tx++) {
	  XI = cgcm[icg][XX]+tx*box[XX][XX]+ty*box[YY][XX]+tz*box[ZZ][XX];
	  /* Calculate range of cells in X direction that have the shift tx */
	  get_dx(Nx,gridx,grid_x,rcoul2,XI,&dx0,&dx1,dcx2);
	  if (dx0 > dx1)
	    continue;
	  /* Get shift vector */	  
	  shift=XYZ2IS(tx,ty,tz);
#ifdef NS5DB
	  assert(shift >= 0);
	  assert(shift < SHIFTS);
#endif
	  for(nn=0; (nn<ngid); nn++) {
	    nsr[nn]      = 0;
	    nlr_ljc[nn]  = 0;
	    nlr_one[nn] = 0;
	  }
#ifdef NS5DB
	  fprintf(log,"shift: %2d, dx0,1: %2d,%2d, dy0,1: %2d,%2d, dz0,1: %2d,%2d\n",
		  shift,dx0,dx1,dy0,dy1,dz0,dz1);
	  fprintf(log,"cgcm: %8.3f  %8.3f  %8.3f\n",cgcm[icg][XX],
		  cgcm[icg][YY],cgcm[icg][ZZ]);
	  fprintf(log,"xi:   %8.3f  %8.3f  %8.3f\n",XI,YI,ZI);
#endif
	  for (dx=dx0; (dx<=dx1); dx++) {
	    tmp1 = rcoul2 - dcx2[dx];
	    for (dy=dy0; (dy<=dy1); dy++) {
	      tmp2 = tmp1 - dcy2[dy];
	      if (tmp2 > 0)
		for (dz=dz0; (dz<=dz1); dz++) {
		  if (tmp2 > dcz2[dz]) {
		    /* Find grid-cell cj in which possible neighbours are */
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
			  /* jgid = gid[cgsatoms[cgsindex[jjcg]]]; */
			  jgid = gid[cgsindex[jjcg]];
			  /* check energy group exclusions */
			  if (!i_eg_excl[jgid]) {
			    if (r2 < rs2) {
			      if (nsr[jgid] >= MAX_CG) {
				put_in_list(bHaveLJ,ngid,md,icg,jgid,
					    nsr[jgid],nl_sr[jgid],
					    cgsindex,cgsatoms,bexcl,
					    shift,fr,FALSE,FALSE,FALSE);
				nsr[jgid]=0;
			      }
			      nl_sr[jgid][nsr[jgid]++]=jjcg;
			    } else if (r2 < rm2) {
			      if (nlr_ljc[jgid] >= MAX_CG) {
				do_longrange(log,top,fr,ngid,md,icg,jgid,
					     nlr_ljc[jgid],
					     nl_lr_ljc[jgid],bexcl,shift,x,
					     box_size,nrnb,
					     lambda,dvdlambda,grps,
					     FALSE,FALSE,FALSE,
					     bHaveLJ);
				nlr_ljc[jgid]=0;
			      }
			      nl_lr_ljc[jgid][nlr_ljc[jgid]++]=jjcg;
			    } else {
			      if (nlr_one[jgid] >= MAX_CG) {
				do_longrange(log,top,fr,ngid,md,icg,jgid,
					     nlr_one[jgid],
					     nl_lr_one[jgid],bexcl,shift,x,
					     box_size,nrnb,
					     lambda,dvdlambda,grps,
					     bVDWOnly,bCoulOnly,FALSE,
					     bHaveLJ);
				nlr_one[jgid]=0;
			      }
			      nl_lr_one[jgid][nlr_one[jgid]++]=jjcg;
			    }
			  }
			}
			nns++;
		      }
		    }
		  }
		}
	    }
	  }
	  /* CHECK whether there is anything left in the buffers */
	  for(nn=0; (nn<ngid); nn++) {
	    if (nsr[nn] > 0)
	      put_in_list(bHaveLJ,ngid,md,icg,nn,nsr[nn],nl_sr[nn],
			  cgsindex,cgsatoms,bexcl,shift,fr,FALSE,FALSE,FALSE);
	    
	    if (nlr_ljc[nn] > 0) 
	      do_longrange(log,top,fr,ngid,md,icg,nn,nlr_ljc[nn],
			   nl_lr_ljc[nn],bexcl,shift,x,box_size,nrnb,
			   lambda,dvdlambda,grps,FALSE,FALSE,FALSE,bHaveLJ);
	    
	    if (nlr_one[nn] > 0) 
	      do_longrange(log,top,fr,ngid,md,icg,nn,nlr_one[nn],
			   nl_lr_one[nn],bexcl,shift,x,box_size,nrnb,
			   lambda,dvdlambda,grps,
			   bVDWOnly,bCoulOnly,FALSE,bHaveLJ);
	  }
	}
      }
    }
    setexcl(nri,i_atoms,&top->atoms.excl,FALSE,bexcl);
  }
  /* Perform any left over force calculations */
  for (nn=0; (nn<ngid); nn++) {
    if (rm2 > rs2)
      do_longrange(log,top,fr,0,md,icg,nn,nlr_ljc[nn],
		   nl_lr_ljc[nn],bexcl,shift,x,box_size,nrnb,
		   lambda,dvdlambda,grps,FALSE,FALSE,TRUE,bHaveLJ);
    if (rl2 > rm2)
      do_longrange(log,top,fr,0,md,icg,nn,nlr_one[nn],
		   nl_lr_one[nn],bexcl,shift,x,box_size,nrnb,
		   lambda,dvdlambda,grps,bVDWOnly,bCoulOnly,TRUE,bHaveLJ);
  }
  debug_gmx();
  
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
  
  t_block  *cgs=&(top->blocks[ebCGS]);
  rvec     box_size;
  int      i,j,m,ngid;
  real     min_size;

  int      nsearch;
  bool     bGrid;
  char     *ptr;
  
  /* Set some local variables */
  bGrid=fr->bGrid;
  ngid=top->atoms.grps[egcENER].nr;
  
  for(m=0; (m<DIM); m++)
    box_size[m]=box[m][m];
  
  if (fr->eBox != ebtNONE) {
    if (bGrid) {
      min_size = min(norm2(box[XX]),min(norm2(box[YY]),norm2(box[ZZ])));
      if (sqr(2*fr->rlistlong) >= min_size)
	fatal_error(0,"One of the box vectors has become shorter than twice the cut-off length.");
    } else {
      min_size = min(box_size[XX],min(box_size[YY],box_size[ZZ]));
      if (2*fr->rlistlong >= min_size)
	fatal_error(0,"One of the box lengths has become smaller than twice the cut-off length.");
    }
  }

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
  
  if (bGrid)
    grid_first(log,grid,box,fr->rlistlong);
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
    nsearch = ns5_core(log,fr,cg_index,box,box_size,ngid,top,grps,
		       grid,x,bexcl,nrnb,md,lambda,dvdlambda,bHaveLJ);
  else {
    /* Only allocate this when necessary, saves 100 kb */
    if (ns_buf == NULL) {
      snew(ns_buf,ngid);
      for(i=0; (i<ngid); i++)
	snew(ns_buf[i],SHIFTS);
    }
    nsearch = ns_simple_core(fr,top,md,box,box_size,
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


