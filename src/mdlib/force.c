/*
 * $Id$
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
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "smalloc.h"
#include "macros.h"
#include "force.h"
#include "invblock.h"
#include "confio.h"
#include "nsb.h"
#include "names.h"
#include "network.h"
#include "wnblist.h"
#include "pbc.h"
#include "ns.h"
#include "nrnb.h"
#include "bondf.h"
#include "mshift.h"
#include "txtdump.h"
#include "ewald_util.h"
#include "shift_util.h"
#include "pppm.h"
#include "ewald.h"
#include "pme.h"
#include "mdrun.h"

t_forcerec *mk_forcerec(void)
{
  t_forcerec *fr;
  
  snew(fr,1);
  
  return fr;
}

#ifdef DEBUG
static void pr_nbfp(FILE *fp,real *nbfp,bool bBHAM,int atnr)
{
  int i,j;
  
  for(i=0; (i<atnr); i++) {
    for(j=0; (j<atnr); j++) {
      fprintf(fp,"%2d - %2d",i,j);
      if (bBHAM)
	fprintf(fp,"  a=%10g, b=%10g, c=%10g\n",BHAMA(nbfp,atnr,i,j),
		BHAMB(nbfp,atnr,i,j),BHAMC(nbfp,atnr,i,j));
      else
	fprintf(fp,"  c6=%10g, c12=%10g\n",C6(nbfp,atnr,i,j),
		C12(nbfp,atnr,i,j));
    }
  }
}
#endif

static real *mk_nbfp(const t_idef *idef,bool bBHAM)
{
  real *nbfp;
  int  i,j,k,atnr;
  
  atnr=idef->atnr;
  if (bBHAM) {
    snew(nbfp,3*atnr*atnr);
    for(i=k=0; (i<atnr); i++) {
      for(j=0; (j<atnr); j++,k++) {
	BHAMA(nbfp,atnr,i,j) = idef->iparams[k].bham.a;
	BHAMB(nbfp,atnr,i,j) = idef->iparams[k].bham.b;
	BHAMC(nbfp,atnr,i,j) = idef->iparams[k].bham.c;
      }
    }
  }
  else {
    snew(nbfp,2*atnr*atnr);
    for(i=k=0; (i<atnr); i++) {
      for(j=0; (j<atnr); j++,k++) {
	C6(nbfp,atnr,i,j)   = idef->iparams[k].lj.c6;
	C12(nbfp,atnr,i,j)  = idef->iparams[k].lj.c12;
      }
    }
  }
  return nbfp;
}

static void check_solvent(FILE *fp,const t_topology *top,t_forcerec *fr,
			  const t_mdatoms *md,const t_nsborder *nsb)
{
  /* This routine finds out whether a charge group can be used as
   * solvent in the innerloops. The routine should be called once
   * at the beginning of the MD program.
   */
  const t_block *cgs,*excl,*mols;
  atom_id *cgid;
  int     i,j,m,j0,j1,nj,k,aj,ak,tjA,tjB,nl_m,nl_n,nl_o;
  int     warncount;
  bool    bOneCG;
  bool    *bAllExcl,bAE,bOrder;
  bool    *bHaveLJ,*bHaveCoul;
  
  cgs  = &(top->blocks[ebCGS]);
  excl = &(top->atoms.excl);
  mols = &(top->blocks[ebMOLS]);

  if (fp)
    fprintf(fp,"Going to determine what solvent types we have.\n");
  snew(fr->solvent_type,cgs->nr+1);
  snew(fr->mno_index,(cgs->nr+1)*3);
  
  /* Generate charge group number for all atoms */
  cgid = make_invblock(cgs,cgs->nra);

  /* Statistics for MFlops counting */  
  fr->nMNOMol = 0;
  fr->nWatMol = 0;
  for(m=0; m<3; m++)
    fr->nMNOav[m] = 0;
  
  warncount=0;

  /* Loop over molecules */
  if (fp)
    fprintf(fp,"There are %d molecules, %d charge groups and %d atoms\n",
	    mols->nr,cgs->nr,cgs->nra);
  for(i=0; (i<mols->nr); i++) {
    j = mols->a[mols->index[i]];
    if (j>=START(nsb) && j<START(nsb)+HOMENR(nsb)) {
      /* Set boolean that determines whether the molecules consists of one CG */
      bOneCG = TRUE;
      /* Set some counters */
      j0     = mols->index[i];
      j1     = mols->index[i+1];
      nj     = j1-j0;
      for(j=j0+1; (j<j1); j++) {
	bOneCG = bOneCG && (cgid[mols->a[j]] == cgid[mols->a[j-1]]);
      }
      if (fr->bSolvOpt && bOneCG && nj>1) {
	/* Check whether everything is excluded */
	snew(bAllExcl,nj);
	bAE = TRUE;
	/* Loop over all atoms in molecule */
	for(j=j0; (j<j1) && bAE; j++) {
	  /* Set a flag for each atom in the molecule that determines whether
	   * it is excluded or not 
	   */
	  for(k=0; (k<nj); k++)
	    bAllExcl[k] = FALSE;
	  /* Now check all the exclusions of this atom */
	  for(k=excl->index[j]; (k<excl->index[j+1]); k++) {
	    ak = excl->a[k];
	    /* Consistency and range check */
	    if ((ak < j0) || (ak >= j1)) 
	      gmx_fatal(FARGS,"Exclusion outside molecule? ak = %d, j0 = %d, j1 = 5d, mol is %d",ak,j0,j1,i);
	    bAllExcl[ak-j0] = TRUE;
	  }
	  /* Now sum up the booleans */
	  for(k=0; (k<nj); k++)
	    bAE = bAE && bAllExcl[k];
	}
	if (bAE) {
	  snew(bHaveCoul,nj);
	  snew(bHaveLJ,nj);
	  for(j=j0; (j<j1); j++) {
	    /* Check for coulomb */
	    aj = mols->a[j];
	    bHaveCoul[j-j0] = ((top->atoms.atom[aj].q != 0.0) ||
			       (top->atoms.atom[aj].qB != 0.0));
	    /* Check for LJ. */
	    tjA = top->atoms.atom[aj].type;
	    tjB = top->atoms.atom[aj].typeB;
	    bHaveLJ[j-j0] = FALSE;
	    for(k=0; (k<fr->ntype); k++) {
	      if (fr->bBHAM) 
		bHaveLJ[j-j0] = (bHaveLJ[j-j0] || 
				 (BHAMA(fr->nbfp,fr->ntype,tjA,k) != 0.0) ||
				 (BHAMB(fr->nbfp,fr->ntype,tjA,k) != 0.0) ||
				 (BHAMC(fr->nbfp,fr->ntype,tjA,k) != 0.0) ||
				 (BHAMA(fr->nbfp,fr->ntype,tjB,k) != 0.0) ||
				 (BHAMB(fr->nbfp,fr->ntype,tjB,k) != 0.0) ||
				 (BHAMC(fr->nbfp,fr->ntype,tjB,k) != 0.0));
	      else
		bHaveLJ[j-j0] = (bHaveLJ[j-j0] || 
				 (C6(fr->nbfp,fr->ntype,tjA,k)  != 0.0) ||
				 (C12(fr->nbfp,fr->ntype,tjA,k) != 0.0) ||
				 (C6(fr->nbfp,fr->ntype,tjB,k)  != 0.0) ||
				 (C12(fr->nbfp,fr->ntype,tjB,k) != 0.0));
	    }
	  }
	  /* Now we have determined what particles have which interactions 
	   * In the case of water-like molecules we only check for the number
	   * of particles and the LJ, not for the Coulomb. Let's just assume
	   * that the water loops are faster than the MNO loops anyway. DvdS
	   */
	  /* No - there's another problem: To optimize the water
	   * innerloop assumes the charge of the first i atom is constant
	   * qO, and charge on atoms 2/3 is constant qH. /EL
	   */
	  /* I won't write any altivec versions of the general solvent inner 
	   * loops. Thus, when USE_PPC_ALTIVEC is defined it is faster 
	   * to use the normal loops instead of the MNO solvent version. /EL
	   */
	  aj=mols->a[j0];
	  if((nj==3) && bHaveCoul[0] && bHaveLJ[0] &&
	     !bHaveLJ[1] && !bHaveLJ[2] &&
	     (top->atoms.atom[aj+1].q == top->atoms.atom[aj+2].q)) {
	    fr->solvent_type[cgid[aj]] = esolWATER;
	    fr->nWatMol++;
	  }
	  else {
#ifdef USE_PPC_ALTIVEC
	    fr->solvent_type[cgid[aj]] = esolNO;
#else
	    /* Time to compute M & N & O */
	    for(k=0; (k<nj) && (bHaveLJ[k] && bHaveCoul[k]); k++)
	      ;
	    nl_n = k;
	    for(; (k<nj) && (!bHaveLJ[k] && bHaveCoul[k]); k++)
	      ;
	    nl_o = k;
	    for(; (k<nj) && (bHaveLJ[k] && !bHaveCoul[k]); k++)
	      ;
	    nl_m = k;
	    /* Now check whether we're at the end of the pack */
	    bOrder = FALSE;
	    for(; (k<nj); k++)
	      bOrder = bOrder || (bHaveLJ[k] || bHaveCoul[k]);
	    if (bOrder) {
	      /* If we have a solvent molecule with LJC everywhere, then
	       * we shouldn't issue a warning. Only if we suspect something
	       * could be better.
	       */
	      if (nl_n != nj) {
		warncount++;
		if(warncount<11) 
		  fprintf(fp,"The order in molecule %d could be optimized"
			  " for better performance\n",i);
		if(warncount==10)
		  fprintf(fp,"(More than 10 molecules where the order can be optimized)\n");
	      }
	      nl_m = nl_n = nl_o = nj;
	    }
	    fr->mno_index[cgid[aj]*3]   = nl_m;
	    fr->mno_index[cgid[aj]*3+1] = nl_n;
	    fr->mno_index[cgid[aj]*3+2] = nl_o;
	    fr->solvent_type[cgid[aj]]  = esolMNO;
	    /* Now compute the number of solvent molecules for MFlops accounting.
	     * The base for this counting is not determined by the order in the
	     * algorithm but by the number of actual interactions.
	     */
	    
	    fr->nMNOMol++;
	    for(k=0; (k<nj); k++)
	      if (bHaveLJ[k] || bHaveCoul[k])
		fr->nMNOav[0]++;
	    for(k=0; (k<nj); k++)
	      if (!bHaveLJ[k] && bHaveCoul[k])
		fr->nMNOav[1]++;
	    for(k=nl_m=0; (k<nj); k++)
	      if (bHaveLJ[k] && !bHaveCoul[k])
		fr->nMNOav[2]++;
#endif /* MNO solvent if not using altivec */
	  }    
	  /* Last check for perturbed atoms */
	  for(j=j0; (j<j1); j++)
	    if (md->bPerturbed[mols->a[j]])
	      fr->solvent_type[cgid[mols->a[j0]]] = esolNO;
	  
	  sfree(bHaveLJ);
	  sfree(bHaveCoul);
	}
	else {
	  /* Turn off solvent optimization for all cg's in the molecule,
	   * here there is only one.
	   */
	  fr->solvent_type[cgid[mols->a[j0]]] = esolNO;
	}
	sfree(bAllExcl);
      }
      else {
	/* Turn off solvent optimization for all cg's in the molecule */
	for(j=mols->index[i]; (j<mols->index[i+1]); j++) {
	  fr->solvent_type[cgid[mols->a[j]]] = esolNO;
	}
      }
    }      
  }
  if (debug) {
    for(i=0; (i<cgs->nr); i++) 
      fprintf(debug,"MNO: cg = %5d, m = %2d, n = %2d, o = %2d\n",
	      i,fr->mno_index[3*i],fr->mno_index[3*i+1],fr->mno_index[3*i+2]);
  }
  
  if (fr->nMNOMol > 0)
    for(m=0; (m<3); m++)
      fr->nMNOav[m] /= fr->nMNOMol;
  
  sfree(cgid);

  if (fp) {
    fprintf(fp,"There are %d optimized solvent molecules on node %d\n",
	    fr->nMNOMol,nsb->nodeid);
    if (fr->nMNOMol > 0)
      fprintf(fp,"  aver. nr. of atoms per molecule: vdwc %.1f coul %.1f vdw %.1f\n",
	      fr->nMNOav[0],fr->nMNOav[1],fr->nMNOav[2]);
    fprintf(fp,"There are %d optimized water molecules on node %d\n",
	    fr->nWatMol,nsb->nodeid);
  }
}

void set_chargesum(FILE *log,t_forcerec *fr,const t_mdatoms *mdatoms)
{
  double qsum;
  int    i;

  qsum = 0;
  for(i=0; i<mdatoms->nr; i++)
    qsum += mdatoms->chargeA[i];
  fr->qsum[0] = qsum;
  if (fr->efep != efepNO) {
    qsum = 0;
    for(i=0; i<mdatoms->nr; i++)
    qsum += mdatoms->chargeB[i];
    fr->qsum[1] = qsum;
  } else {
    fr->qsum[1] = fr->qsum[0];
  }
  if (log) {
    if (fr->efep == efepNO)
      fprintf(log,"System total charge: %.3f\n",fr->qsum[0]);
    else
      fprintf(log,"System total charge, top. A: %.3f top. B: %.3f\n",
	      fr->qsum[0],fr->qsum[1]);
  }
}

void update_forcerec(FILE *log,t_forcerec *fr,matrix box)
{
  fr->epsfac = ONE_4PI_EPS0;
  if (EEL_RF(fr->eeltype))
    calc_rffac(log,fr->eeltype,
	       fr->epsilon_r,fr->rcoulomb,fr->temp,fr->zsquare,box,
	       &fr->kappa,&fr->k_rf,&fr->c_rf);
  else {
    if (fr->epsilon_r != 0)
      /* multiply the dielectric constant by 1/eps */
      fr->epsfac /= fr->epsilon_r;
    else
      /* eps = 0 is infinite dieletric: no coulomb interactions */
      fr->epsfac = 0;
  }
}

void set_avcsixtwelve(FILE *log,t_forcerec *fr,
		      const t_mdatoms *mdatoms,const t_block *excl)
{
  int    i,j,tpi,tpj,j1,j2,k,nexcl;
  double csix,ctwelve;
  int    natoms,ntp,*type;
  bool   bBHAM;
  real   *nbfp;
  atom_id *AA;

  natoms = mdatoms->nr;
  ntp = fr->ntype;
  type = mdatoms->typeA;
  bBHAM = fr->bBHAM;
  nbfp = fr->nbfp;
  AA = excl->a;

  csix = 0;
  ctwelve = 0;
  /* We loop over all the atom pairs and subtract the excluded pairs.
   * The main reason for substracting exclusions is that in some cases some
   * combinations might never occur and the parameters could have any value.
   * These unused values should not influence the dispersion correction.
   */
  nexcl = 0;
  for(i=0; (i<natoms); i++) {
    tpi = type[i];
#ifdef DEBUG
    if (tpi >= ntp)
      gmx_fatal(FARGS,"Atomtype[%d] = %d, maximum = %d",i,tpi,ntp);
#endif
    for(j=i+1; (j<natoms); j++) {
      tpj   = type[j];
#ifdef DEBUG
      if (tpj >= ntp)
	gmx_fatal(FARGS,"Atomtype[%d] = %d, maximum = %d",j,tpj,ntp);
#endif
      if (bBHAM) {
	csix += BHAMC(nbfp,ntp,tpi,tpj);
      } else {
	csix    += C6 (nbfp,ntp,tpi,tpj);
	ctwelve += C12(nbfp,ntp,tpi,tpj);
      }
    }
    /* Subtract the exclusions */
    j1  = excl->index[i];
    j2  = excl->index[i+1];
    for(j=j1; j<j2; j++) {
      k = AA[j];
      if (k > i) {
	tpj   = type[k];
	if (bBHAM) {
	  csix -= BHAMC(nbfp,ntp,tpi,tpj);
	} else {
	  csix    -= C6 (nbfp,ntp,tpi,tpj);
	  ctwelve -= C12(nbfp,ntp,tpi,tpj);
	}
	nexcl++;
      }
    }
  }
  csix    /= 0.5*natoms*(natoms - 1) - nexcl;
  ctwelve /= 0.5*natoms*(natoms - 1) - nexcl;
  if (debug) {
    fprintf(debug,"Counted %d exclusions\n",nexcl);
    fprintf(debug,"Average C6 parameter is: %10g\n",csix);
    fprintf(debug,"Average C12 parameter is: %10g\n",ctwelve);
  }
  fr->avcsix = csix;
  fr->avctwelve = ctwelve;
}

static void set_bham_b_max(FILE *log,t_forcerec *fr,const t_mdatoms *mdatoms)
{
  int  i,j,tpi,tpj,ntypes,natoms,*type;
  real b,bmin;
  real *nbfp;

  fprintf(log,"Determining largest Buckingham b parameter for table\n");
  nbfp   = fr->nbfp;
  ntypes = fr->ntype;
  type   = mdatoms->typeA;
  natoms = mdatoms->nr;

  bmin           = -1;
  fr->bham_b_max = 0;
  for(i=0; (i<natoms); i++) {
    tpi = type[i];
    if (tpi >= ntypes)
      gmx_fatal(FARGS,"Atomtype[%d] = %d, maximum = %d",i,tpi,ntypes);
    
    for(j=0; (j<natoms); j++) {
      tpj   = type[j];
      if (tpj >= ntypes)
	gmx_fatal(FARGS,"Atomtype[%d] = %d, maximum = %d",j,tpj,ntypes);
      b = BHAMB(nbfp,ntypes,tpi,tpj);
      if (b > fr->bham_b_max)
	fr->bham_b_max = b;
      if ((b < bmin) || (bmin==-1))
	bmin = b;
    }
  }
  fprintf(log,"Buckingham b parameters, min: %g, max: %g\n",
	  bmin,fr->bham_b_max);
}

void init_forcerec(FILE *fp,
		   t_forcerec *fr,
		   const t_inputrec *ir,
		   const t_topology *top,
		   const t_commrec  *cr,
		   const t_mdatoms  *mdatoms,
		   const t_nsborder *nsb,
		   matrix     box,
		   bool       bMolEpot,
		   const char *tabfn,
		   bool       bNoSolvOpt)
{
  int     i,j,m,natoms,ngrp;
  real    q,zsq,nrdf,T,rtab;
  rvec    box_size;
  const t_block *mols,*cgs;
  const t_idef  *idef;
  bool    bTab,bSep14tab;
  
  if (check_box(box))
    gmx_fatal(FARGS,check_box(box));

  cgs            = &(top->blocks[ebCGS]);
  mols           = &(top->blocks[ebMOLS]);
  idef           = &(top->idef);
  
  natoms         = mdatoms->nr;

  /* Copy the user determined parameters */
  fr->userint1 = ir->userint1;
  fr->userint2 = ir->userint2;
  fr->userint3 = ir->userint3;
  fr->userint4 = ir->userint4;
  fr->userreal1 = ir->userreal1;
  fr->userreal2 = ir->userreal2;
  fr->userreal3 = ir->userreal3;
  fr->userreal4 = ir->userreal4;

  /* Shell stuff */
  fr->fc_stepsize = ir->fc_stepsize;

  /* Free energy */
  fr->efep       = ir->efep;
  fr->sc_alpha   = ir->sc_alpha;
  fr->sc_sigma6  = pow(ir->sc_sigma,6);

  /* Neighbour searching stuff */
  fr->bGrid      = (ir->ns_type == ensGRID);
  fr->ndelta     = ir->ndelta;
  fr->ePBC       = ir->ePBC;
  fr->rlist      = ir->rlist;
  fr->rlistlong  = max(ir->rlist,max(ir->rcoulomb,ir->rvdw));
  fr->eeltype    = ir->coulombtype;
  fr->vdwtype    = ir->vdwtype;

  fr->bTwinRange = fr->rlistlong > fr->rlist;
  fr->bEwald     = fr->eeltype==eelPME || fr->eeltype==eelEWALD;
  fr->bvdwtab    = fr->vdwtype != evdwCUT;
  fr->bcoultab   = fr->eeltype != eelCUT && !(EEL_RF(fr->eeltype) &&
					      fr->vdwtype == evdwCUT);
  if (getenv("GMX_FORCE_TABLES")) {
    fr->bvdwtab  = TRUE;
    fr->bcoultab = TRUE;
  }

  if (fp) {
    fprintf(fp,"Table routines are used for coulomb: %s\n",bool_names[fr->bcoultab]);
    fprintf(fp,"Table routines are used for vdw:     %s\n",bool_names[fr->bvdwtab ]);
  }
  
  /* Tables are used for direct ewald sum */
  if(fr->bEwald) {
    fr->ewaldcoeff=calc_ewaldcoeff(ir->rcoulomb, ir->ewald_rtol);
    if (fp)
      fprintf(fp,"Using a Gaussian width (1/beta) of %g nm for Ewald\n",
	      1/fr->ewaldcoeff);
  }

  /* Domain decomposition parallellism... */
  fr->bDomDecomp = ir->bDomDecomp;
  fr->Dimension  = ir->decomp_dir;
  
  /* Electrostatics */
  fr->epsilon_r  = ir->epsilon_r;
  fr->fudgeQQ    = ir->fudgeQQ;
  fr->rcoulomb_switch = ir->rcoulomb_switch;
  fr->rcoulomb        = ir->rcoulomb;
  
  if (bNoSolvOpt || getenv("GMX_NO_SOLV_OPT"))
    fr->bSolvOpt = FALSE;
  else
    fr->bSolvOpt = TRUE;
  
  /* Parameters for generalized RF */
  fr->zsquare = 0.0;
  fr->temp    = 0.0;
  
  if (fr->eeltype == eelGRF) {
    if (ir->efep!=efepNO)
      fprintf(fp,"\nWARNING: the generalized reaction field constants are determined from topology A only\n\n");
    zsq = 0.0;
    for (i=0; (i<cgs->nr); i++) {
      q = 0;
      for(j=cgs->index[i]; (j<cgs->index[i+1]); j++)
	q+=mdatoms->chargeA[cgs->a[j]];
      if (q != 0.0)
	/* Changed from square to fabs 990314 DvdS 
	 * Does not make a difference for monovalent ions, but doe for 
	 * divalent ions (Ca2+!!)
	 */
	zsq += fabs(q);
    }
    fr->zsquare = zsq;
    
    T    = 0.0;
    nrdf = 0.0;
    for(i=0; (i<ir->opts.ngtc); i++) {
      nrdf += ir->opts.nrdf[i];
      T    += (ir->opts.nrdf[i] * ir->opts.ref_t[i]);
    }
    if (nrdf == 0) 
      gmx_fatal(FARGS,"No degrees of freedom!");
    fr->temp   = T/nrdf;
  }
  else if (EEL_FULL(fr->eeltype) || (fr->eeltype == eelSHIFT) || 
	   (fr->eeltype == eelUSER) || (fr->eeltype == eelSWITCH)) {
    /* We must use the long range cut-off for neighboursearching...
     * An extra range of e.g. 0.1 nm (half the size of a charge group)
     * is necessary for neighboursearching. This allows diffusion 
     * into the cut-off range (between neighborlist updates), 
     * and gives more accurate forces because all atoms within the short-range
     * cut-off rc must be taken into account, while the ns criterium takes
     * only those with the center of geometry within the cut-off.
     * (therefore we have to add half the size of a charge group, plus
     * something to account for diffusion if we have nstlist > 1)
     */
    for(m=0; (m<DIM); m++)
      box_size[m]=box[m][m];

    if (fr->phi == NULL)
      snew(fr->phi,mdatoms->nr);
    
    if ((fr->eeltype==eelPPPM) || (fr->eeltype==eelPOISSON) || 
	(fr->eeltype == eelSHIFT && fr->rcoulomb > fr->rcoulomb_switch))
	set_shift_consts(fp,fr->rcoulomb_switch,fr->rcoulomb,box_size,fr);
  }

  /* Initiate arrays */
  if (fr->bTwinRange) {
    snew(fr->f_twin,natoms);
    snew(fr->fshift_twin,SHIFTS);
  }
  
  if (EEL_FULL(fr->eeltype)) {
    if (ir->efep != efepNO) {
      fprintf(fp,"\nWARNING: With %s the reciprocal part only uses the charges from topology A\n\n",eel_names[fr->eeltype]);
    }
    snew(fr->f_el_recip,natoms);
  }
  
  /* Mask that says whether or not this NBF list should be computed */
  /*  if (fr->bMask == NULL) {
      ngrp = ir->opts.ngener*ir->opts.ngener;
      snew(fr->bMask,ngrp);*/
  /* Defaults to always */
  /*    for(i=0; (i<ngrp); i++)
	fr->bMask[i] = TRUE;
	}*/
  
  if (fr->cg_cm == NULL)
    snew(fr->cg_cm,cgs->nr);
  if (fr->shift_vec == NULL)
    snew(fr->shift_vec,SHIFTS);
    
  if (fr->fshift == NULL)
    snew(fr->fshift,SHIFTS);
  
  if (bMolEpot && (fr->nmol==0)) {
    fr->nmol=mols->nr;
    fr->mol_nr=make_invblock(mols,natoms);
    snew(fr->mol_epot,fr->nmol);
    fr->nstcalc=ir->nstenergy;
  }
  
  if (fr->nbfp == NULL) {
    fr->ntype = idef->atnr;
    fr->bBHAM = (idef->functype[0] == F_BHAM);
    fr->nbfp  = mk_nbfp(idef,fr->bBHAM);
  }
  /* Copy the energy group exclusions */
  fr->eg_excl = ir->opts.eg_excl;

  /* Van der Waals stuff */
  fr->rvdw        = ir->rvdw;
  fr->rvdw_switch = ir->rvdw_switch;
  if ((fr->vdwtype != evdwCUT) && (fr->vdwtype != evdwUSER) && !fr->bBHAM) {
    if (fr->rvdw_switch >= fr->rvdw)
      gmx_fatal(FARGS,"rvdw_switch (%g) must be < rvdw (%g)",
		  fr->rvdw_switch,fr->rvdw);
    if (fp)
      fprintf(fp,"Using %s Lennard-Jones, switch between %g and %g nm\n",
	      (fr->eeltype==eelSWITCH) ? "switched":"shifted",
	      fr->rvdw_switch,fr->rvdw);
  } 

  if (fr->bBHAM && (fr->vdwtype == evdwSHIFT || fr->vdwtype == evdwSWITCH))
    gmx_fatal(FARGS,"Switch/shift interaction not supported with Buckingham");
  
  if (fp)
    fprintf(fp,"Cut-off's:   NS: %g   Coulomb: %g   %s: %g\n",
	    fr->rlist,fr->rcoulomb,fr->bBHAM ? "BHAM":"LJ",fr->rvdw);
  
  if (ir->eDispCorr != edispcNO)
    set_avcsixtwelve(fp,fr,mdatoms,&top->atoms.excl);

  if (fr->bBHAM)
    set_bham_b_max(fp,fr,mdatoms);

  /* Copy the GBSA data (radius, volume and surftens for each
   * atomtype) from the topology atomtype section to forcerec.
   */
  snew(fr->atype_radius,fr->ntype);
  snew(fr->atype_vol,fr->ntype);
  snew(fr->atype_surftens,fr->ntype);
  if (top->atomtypes.nr > 0) {
    for(i=0;i<fr->ntype;i++)
      fr->atype_radius[i]=top->atomtypes.radius[i];
    for(i=0;i<fr->ntype;i++)
      fr->atype_vol[i]=top->atomtypes.vol[i];
    for(i=0;i<fr->ntype;i++)
      fr->atype_surftens[i]=top->atomtypes.surftens[i];
  }    

  /* Now update the rest of the vars */
  update_forcerec(fp,fr,box);
  
  set_chargesum(fp,fr,mdatoms);

  /* if we are using LR electrostatics, and they are tabulated,
   * the tables will contain shifted coulomb interactions.
   * Since we want to use the non-shifted ones for 1-4
   * coulombic interactions, we must have an extra set of
   * tables. This should be done in tables.c, instead of this
   * ugly hack, but it works for now...
   */

  /* Construct tables.
   * A little unnecessary to make both vdw and coul tables sometimes,
   * but what the heck... */

  bTab = fr->bcoultab || fr->bvdwtab || (fr->efep != efepNO);
  bSep14tab = !bTab || EEL_FULL(fr->eeltype) || (EEL_RF(fr->eeltype) &&
						 fr->eeltype != eelRF_OLD);
  if (bTab) {
    rtab = fr->rlistlong + ir->tabext;
    if (fr->rlistlong == 0 && fr->efep != efepNO && rtab < 5) {
      rtab = 5;
      if (fp)
	fprintf(fp,
		"\nWARNING: Increasing the free energy table length from\n"
		"         table extension = %f nm to %g nm,\n"
		"         you can set the table extension in the mdp file\n\n",
		ir->tabext,rtab);
    }
    /* make tables for ordinary interactions */
    fr->tab = make_tables(fp,fr,MASTER(cr),tabfn,rtab,FALSE);
    /* Copy the contents of the table to separate coulomb and LJ
     * tables too, to improve cache performance.
     */
    snew(fr->coultab,4*(fr->tab.n+1));
    snew(fr->vdwtab,8*(fr->tab.n+1));  
    for(i=0; i<=fr->tab.n; i++) {
      for(j=0; j<4; j++) 
	fr->coultab[4*i+j]=fr->tab.tab[12*i+j];
      for(j=0; j<8; j++) 
	fr->vdwtab[8*i+j]=fr->tab.tab[12*i+4+j];
    }
    if (!bSep14tab)
      fr->tab14 = fr->tab;
  }
  if (bSep14tab)
    /* generate extra tables with plain Coulomb for 1-4 interactions only */
    fr->tab14 = make_tables(fp,fr,MASTER(cr),tabfn,ir->tabext,TRUE);

  if (!fr->mno_index)
    check_solvent(fp,top,fr,mdatoms,nsb);
}
 
#define pr_real(fp,r) fprintf(fp,"%s: %e\n",#r,r)
#define pr_int(fp,i)  fprintf((fp),"%s: %d\n",#i,i)
#define pr_bool(fp,b) fprintf((fp),"%s: %s\n",#b,bool_names[b])

void pr_forcerec(FILE *fp,t_forcerec *fr,t_commrec *cr)
{
  pr_real(fp,fr->rlist);
  pr_real(fp,fr->rcoulomb);
  pr_real(fp,fr->fudgeQQ);
  pr_int(fp,fr->ndelta);
  pr_bool(fp,fr->bGrid);
  pr_bool(fp,fr->bTwinRange);
  /*pr_int(fp,fr->cg0);
    pr_int(fp,fr->hcg);*/
  pr_int(fp,fr->tab.n);
  if (fr->tab.n > 0) {
    pr_real(fp,fr->rcoulomb_switch);
    pr_real(fp,fr->rcoulomb);
  }
  
  pr_int(fp,fr->nmol);
  pr_int(fp,fr->nstcalc);
  
  fflush(fp);
}

void ns(FILE *fp,
	t_forcerec *fr,
	rvec       x[],
	rvec       f[],
	matrix     box,
	t_groups   *grps,
	t_grpopts  *opts,
	t_topology *top,
	t_mdatoms  *md,
	t_commrec  *cr,
	t_nrnb     *nrnb,
	t_nsborder *nsb,
	int        step,
	real       lambda,
	real       *dvdlambda)
{
  static bool bFirst=TRUE;
  static int  nDNL;
  char   *ptr;
  int    nsearch;
  
  if (bFirst) {
    ptr=getenv("DUMPNL");
    if (ptr) {
      nDNL=atoi(ptr);
      fprintf(fp,"nDNL = %d\n",nDNL);  
    } else
      nDNL=0;
    /* Allocate memory for the neighbor lists */
    init_neighbor_list(fp,fr,HOMENR(nsb));
      
    bFirst=FALSE;
  }
    
  if (fr->bTwinRange) 
    fr->nlr=0;

  /* Whether or not we do dynamic load balancing,
   * workload contains the proper numbers of charge groups
   * to be searched.
   */
  if (cr->nodeid == 0)
    fr->cg0=0;
  else
    fr->cg0=nsb->workload[cr->nodeid-1];
  fr->hcg=nsb->workload[cr->nodeid];

  nsearch = search_neighbours(fp,fr,x,box,top,grps,cr,nsb,nrnb,md,
			      lambda,dvdlambda);
  if (debug)
    fprintf(debug,"nsearch = %d\n",nsearch);
    
  /* Check whether we have to do dynamic load balancing */
  /*if ((nsb->nstDlb > 0) && (mod(step,nsb->nstDlb) == 0))
    count_nb(cr,nsb,&(top->blocks[ebCGS]),nns,fr->nlr,
    &(top->idef),opts->ngener);
  */
  if (nDNL > 0)
    dump_nblist(fp,fr,nDNL);
}

void force(FILE       *fplog,   int        step,
	   t_forcerec *fr,      t_inputrec *ir,
	   t_idef     *idef,    t_nsborder *nsb,
	   t_commrec  *cr,      t_commrec *mcr,
	   t_nrnb     *nrnb,
	   t_groups   *grps,    t_mdatoms  *md,
	   int        ngener,   t_grpopts  *opts,
	   rvec       x[],      rvec       f[],
	   real       epot[],   t_fcdata   *fcd,
	   bool       bVerbose, matrix     box,
	   real       lambda,   t_graph    *graph,
	   t_block    *excl,    bool       bNBFonly,
	   rvec       mu_tot[],
	   bool       bGatherOnly)
{
  int     i,nit;
  bool    bDoEpot,bSepDVDL;
  rvec    box_size;
  real    dvdlambda,Vlr,Vcorr=0;
  t_pbc   pbc;
#define PRINT_SEPDVDL(s,v,dvdl) if (bSepDVDL) fprintf(fplog,"  %-30s V %12.5e  dVdl %12.5e\n",s,v,dvdl);

  /* Reset box */
  for(i=0; (i<DIM); i++)
    box_size[i]=box[i][i];
    
  bDoEpot=((fr->nmol > 0) && (fr->nstcalc > 0) && (mod(step,fr->nstcalc)==0));
  bSepDVDL=(fr->bSepDVDL && do_per_step(step,ir->nstlog));
  /* Reset epot... */
  if (bDoEpot) 
    for(i=0; (i<fr->nmol); i++)
      fr->mol_epot[i]=0.0;
  debug_gmx();

  if (bSepDVDL)
    fprintf(fplog,"Step %d: non-bonded V and dVdl for node %d:\n",
	    step,cr->nodeid);
  
  /* Call the short range functions all in one go. */
  dvdlambda = 0;
  do_fnbf(fplog,cr,fr,x,f,md,
	  fr->bBHAM ? grps->estat.ee[egBHAMSR] : grps->estat.ee[egLJSR],
	  grps->estat.ee[egCOULSR],box_size,nrnb,
	  lambda,&dvdlambda,FALSE,-1);
  epot[F_DVDL] += dvdlambda;
  PRINT_SEPDVDL("VdW and Coulomb particle-p.",0.0,dvdlambda);
  debug_gmx();

  if (debug) 
    pr_rvecs(debug,0,"fshift after SR",fr->fshift,SHIFTS);
  
  /* Shift the coordinates. Must be done before bonded forces and PPPM, 
   * but is also necessary for SHAKE and update, therefore it can NOT 
   * go when no bonded forces have to be evaluated.
   */
  
  /* Check whether we need to do bondeds or correct for exclusions */
  if (!bNBFonly || EEL_RF(fr->eeltype) || EEL_FULL(fr->eeltype)) {
    if (graph) {
      shift_self(graph,box,x);
      if (TRICLINIC(box))
	inc_nrnb(nrnb,eNR_SHIFTX,2*graph->nnodes);
      else
	inc_nrnb(nrnb,eNR_SHIFTX,graph->nnodes);
    }
    /* We may need the pbc structure for e.g. position restraints
     * or when pbc=full
     */
    set_pbc(&pbc,box);
    debug_gmx();
  }
  
  if (EEL_FULL(fr->eeltype)) {
    dvdlambda = 0;
    switch (fr->eeltype) {
    case eelPPPM:
      Vlr = do_pppm(fplog,FALSE,x,fr->f_el_recip,md->chargeA,
		    box_size,fr->phi,cr,nsb,nrnb);
      break;
    case eelPME:
      Vlr = do_pme(fplog,FALSE,ir,x,fr->f_el_recip,md->chargeA,
		   box,cr,nsb,nrnb,fr->vir_el_recip,fr->ewaldcoeff,
		   bGatherOnly);
      PRINT_SEPDVDL("PME mesh",Vlr,dvdlambda);
      break;
    case eelEWALD:
      Vlr = do_ewald(fplog,FALSE,ir,x,fr->f_el_recip,md->chargeA,
		     box_size,cr,nsb,fr->vir_el_recip,fr->ewaldcoeff);
      PRINT_SEPDVDL("Ewald long-range",Vlr,dvdlambda);
      break;
    default:
      Vlr = 0;
      gmx_fatal(FARGS,"No such electrostatics method implemented %s",
		  eel_names[fr->eeltype]);
    }
    epot[F_DVDL] += dvdlambda;
    if(fr->bEwald) {
      dvdlambda = 0;
      Vcorr = ewald_LRcorrection(fplog,nsb,cr,fr,md->chargeA,excl,
				 x,box,mu_tot[0],fr->qsum[0],
				 ir->ewald_geometry,ir->epsilon_surface,
				 fr->vir_el_recip);
      PRINT_SEPDVDL("Ewald excl./charge/dip. corr.",Vcorr,dvdlambda);
      epot[F_DVDL] += dvdlambda;
    } else {
      Vcorr = shift_LRcorrection(fplog,nsb,cr,fr,md->chargeA,excl,x,TRUE,box,
				 fr->vir_el_recip);
    }
    epot[F_COUL_RECIP] = Vlr + Vcorr;
    if (debug)
      fprintf(debug,"Vlr = %g, Vcorr = %g, Vlr_corr = %g\n",
	      Vlr,Vcorr,epot[F_COUL_RECIP]);
    if (debug) {
      pr_rvecs(debug,0,"vir_el_recip after corr",fr->vir_el_recip,DIM);
      pr_rvecs(debug,0,"fshift after LR Corrections",fr->fshift,SHIFTS);
    }
  } else if (EEL_RF(fr->eeltype)) {
    dvdlambda = 0;
    if (fr->eeltype != eelRF_OLD)
      epot[F_RF_EXCL] = RF_excl_correction(fplog,nsb,fr,graph,md,excl,x,f,
					   fr->fshift,&pbc,lambda,&dvdlambda);
    epot[F_DVDL] += dvdlambda;
    PRINT_SEPDVDL("RF exclusion correction",epot[F_RF_EXCL],dvdlambda);
  }
  debug_gmx();
  
  if (debug)    
    print_nrnb(debug,nrnb); 
  debug_gmx();
  
  if (!bNBFonly) {
    calc_bonds(fplog,cr,mcr,
	       idef,x,f,fr,&pbc,graph,epot,nrnb,lambda,md,
	       opts->ngener,grps->estat.ee[egLJ14],grps->estat.ee[egCOUL14],
	       fcd,step,fr->bSepDVDL && do_per_step(step,ir->nstlog));    
    debug_gmx();
  }
  if (debug) 
    pr_rvecs(debug,0,"fshift after bondeds",fr->fshift,SHIFTS);
}
