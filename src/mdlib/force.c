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
#include "assert.h"
#include "physics.h"
#include "macros.h"
#include "vec.h"
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
#include "copyrite.h"

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

static real *mk_nbfp(t_idef *idef,bool bBHAM)
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

static void check_solvent(FILE *fp,t_topology *top,t_forcerec *fr,
			  t_mdatoms *md,t_nsborder *nsb)
{
  /* This routine finds out whether a charge group can be used as
   * solvent in the innerloops. The routine should be called once
   * at the beginning of the MD program.
   */
  t_block *cgs,*excl,*mols;
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
  
  warncount=0;

  /* Loop over molecules */
  if (fp)
    fprintf(fp,"There are %d molecules, %d charge groups and %d atoms\n",
	    mols->nr,cgs->nr,cgs->nra);
  for(i=0; (i<mols->nr); i++) {
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
	    fatal_error(0,"Exclusion outside molecule? ak = %d, j0 = %d, j1 = 5d, mol is %d",ak,j0,j1,i);
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
	   (top->atoms.atom[aj+1].q == top->atoms.atom[aj+2].q))
	  fr->solvent_type[cgid[aj]] = esolWATER;
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
  if (debug) {
    for(i=0; (i<cgs->nr); i++) 
      fprintf(debug,"MNO: cg = %5d, m = %2d, n = %2d, o = %2d\n",
	      i,fr->mno_index[3*i],fr->mno_index[3*i+1],fr->mno_index[3*i+2]);
  }

  /* Now compute the number of solvent molecules, could be merged with code above */
  fr->nMNOMol = 0;
  fr->nWatMol = 0;
  for(m=0; m<3; m++)
    fr->nMNOav[m] = 0;
  for(i=0; i<mols->nr; i++) {
    j = mols->a[mols->index[i]];
    if (j>=START(nsb) && j<START(nsb)+HOMENR(nsb)) {
	if (fr->solvent_type[cgid[j]] == esolMNO) {
	  fr->nMNOMol++;
	  for(m=0; m<3; m++)
	    fr->nMNOav[m] += fr->mno_index[3*cgid[j]+m];
	}
	else if (fr->solvent_type[cgid[j]] == esolWATER)
	  fr->nWatMol++;
    }
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
	      fr->nMNOav[1],fr->nMNOav[2]-fr->nMNOav[1],fr->nMNOav[0]-fr->nMNOav[2]);
    fprintf(fp,"There are %d optimized water molecules on node %d\n",
	    fr->nWatMol,nsb->nodeid);
  }
}

static void calc_rffac(FILE *log,int eel,real eps,real Rc,real Temp,
		       real zsq,matrix box,
		       real *kappa,real *epsfac,real *krf,real *crf)
{
  /* Compute constants for Generalized reaction field */
  static bool bFirst=TRUE;
  real   k1,k2,I,vol,rmin;
  
  if ((eel == eelRF) || (eel == eelGRF)) {
    vol     = det(box);
    I       = zsq/vol;
    if (eel == eelGRF) {
      /* Consistency check */
      if (Temp <= 0.0)
	fatal_error(0,"Temperature is %f while using"
		    " Generalized Reaction Field\n",Temp);
      /* Ionic strength (only needed for eelGRF */
      *kappa  = sqrt(2*I/(EPSILON0*eps*BOLTZ*Temp));
    }
    else
      *kappa = 0;

    /* eps == 0 signals infinite dielectric */
    if (eps == 0) {
      *krf = 1/(2*Rc*Rc*Rc);
      *crf = 0;
    }
    else {
      k1      = (1+*kappa*Rc);
      k2      = eps*sqr((real)(*kappa*Rc));
      
      *krf    = (((eps-1)*k1+k2)/((2*eps+1)*k1+2*k2)/(Rc*Rc*Rc));
      *crf    = 1/Rc + *krf*Rc*Rc;
    }
    *epsfac = ONE_4PI_EPS0;
    rmin    = pow(*krf*2.0,-1.0/3.0);
    
    if (bFirst) {
      if (eel == eelGRF)
	please_cite(log,"Tironi95a");
      fprintf(log,"%s:\n"
	      "epsRF = %10g, I   = %10g, volume = %10g, kappa  = %10g\n"
	      "rc    = %10g, krf = %10g, crf    = %10g, epsfac = %10g\n",
	      eel_names[eel],eps,I,vol,*kappa,Rc,*krf,*crf,*epsfac);
      fprintf(log,
	      "The electrostatics potential has its minimum at rc = %g\n",
	      rmin);
      
      bFirst=FALSE;
    }
  }
  else {
    /* If we're not using a reaction field, set the factor to 0
     * and multiply the dielectric constant by 1/eps
     */
    *kappa  = 0.0;
    *krf    = 0.0;
    *crf    = 0.0;
    if (eps == 0)
      eps = 1;
    *epsfac = ONE_4PI_EPS0/eps;
  } 
}

void update_forcerec(FILE *log,t_forcerec *fr,matrix box)
{
  calc_rffac(log,fr->eeltype,
	     fr->epsilon_r,fr->rcoulomb,fr->temp,fr->zsquare,box,
	     &fr->kappa,&fr->epsfac,&fr->k_rf,&fr->c_rf);
}

static double calc_avcsix(FILE *log,real *nbfp,int atnr,
			  int natoms,int type[],bool bBHAM)
{
  int    i,j,tpi,tpj;
  double csix;

  /* Check this code: do we really need a double loop? */  
  csix = 0;
  for(i=0; (i<natoms); i++) {
    tpi = type[i];
#ifdef DEBUG
    if (tpi >= atnr)
      fatal_error(0,"Atomtype[%d] = %d, maximum = %d",i,tpi,atnr);
#endif
    for(j=0; (j<natoms); j++) {
      tpj   = type[j];
#ifdef DEBUG
      if (tpj >= atnr)
	fatal_error(0,"Atomtype[%d] = %d, maximum = %d",j,tpj,atnr);
#endif
      if (bBHAM)
	csix += BHAMC(nbfp,atnr,tpi,tpj);
      else
	csix += C6(nbfp,atnr,tpi,tpj);
    }
  }
  csix /= ((double)natoms*(double)natoms);
  if (debug)
    fprintf(debug,"Average C6 parameter is: %10g\n",csix);
  
  return csix;
}

void set_avcsix(FILE *log,t_forcerec *fr,t_mdatoms *mdatoms)
{
  fr->avcsix=calc_avcsix(log,fr->nbfp,fr->ntype,mdatoms->nr,
			 mdatoms->typeA,fr->bBHAM);
}

static double calc_avctwelve(FILE *log,real *nbfp,int atnr,
			     int natoms,int type[],bool bBHAM)
{
  int    i,j,tpi,tpj;
  double ctwel;
  
  /* Check this code: do we really need a double loop? */  
  ctwel = 0;
  for(i=0; (i<natoms); i++) {
    tpi = type[i];
#ifdef DEBUG
    if (tpi >= atnr)
      fatal_error(0,"Atomtype[%d] = %d, maximum = %d",i,tpi,atnr);
#endif
    for(j=0; (j<natoms); j++) {
      tpj   = type[j];
#ifdef DEBUG
      if (tpj >= atnr)
        fatal_error(0,"Atomtype[%d] = %d, maximum = %d",j,tpj,atnr);
#endif
      if (bBHAM)
        /*no repulsion correction for Buckingham for now*/
        ctwel += 0;
      else
        ctwel += C12(nbfp,atnr,tpi,tpj);
    }
  }
  ctwel /= (natoms*natoms);
  if (debug)
    fprintf(debug,"Average C12 parameter is: %10g\n",ctwel);
  
  return ctwel;
}

void set_avctwelve(FILE *log,t_forcerec *fr,t_mdatoms *mdatoms)
{
  fr->avctwelve=calc_avctwelve(log,fr->nbfp,fr->ntype,mdatoms->nr,
			       mdatoms->typeA,fr->bBHAM);
}


static void set_bham_b_max(FILE *log,t_forcerec *fr,t_mdatoms *mdatoms)
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
      fatal_error(0,"Atomtype[%d] = %d, maximum = %d",i,tpi,ntypes);
    
    for(j=0; (j<natoms); j++) {
      tpj   = type[j];
      if (tpj >= ntypes)
	fatal_error(0,"Atomtype[%d] = %d, maximum = %d",j,tpj,ntypes);
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
		   t_inputrec *ir,
		   t_topology *top,
		   t_commrec  *cr,
		   t_mdatoms  *mdatoms,
		   t_nsborder *nsb,
		   matrix     box,
		   bool       bMolEpot,
		   char       *tabfn,
		   bool       bNoSolvOpt)
{
  int     i,j,m,natoms,ngrp;
  real    q,zsq,nrdf,T;
  rvec    box_size;
  double  rtab;
  t_block *mols,*cgs;
  t_idef  *idef;

  if (check_box(box))
    fatal_error(0,check_box(box));

  cgs            = &(top->blocks[ebCGS]);
  mols           = &(top->blocks[ebMOLS]);
  idef           = &(top->idef);
  
  natoms         = mdatoms->nr;

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
  fr->bRF        = (fr->eeltype==eelRF || fr->eeltype==eelGRF) &&
		    fr->vdwtype==evdwCUT;
  fr->bcoultab   = (fr->eeltype!=eelCUT && !fr->bRF) || fr->bEwald;

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
    zsq = 0.0;
    for (i=0; (i<cgs->nr); i++) {
      q = 0;
      for(j=cgs->index[i]; (j<cgs->index[i+1]); j++)
	q+=mdatoms->chargeT[cgs->a[j]];
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
      fatal_error(0,"No degrees of freedom!");
    fr->temp   = T/nrdf;
  }
  else if (EEL_LR(fr->eeltype) || (fr->eeltype == eelSHIFT) || 
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
  
  if (EEL_LR(fr->eeltype)) {
    snew(fr->f_pme,natoms);
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
      fatal_error(0,"rvdw_switch (%g) must be < rvdw (%g)",
		  fr->rvdw_switch,fr->rvdw);
    if (fp)
      fprintf(fp,"Using %s Lennard-Jones, switch between %g and %g nm\n",
	      (fr->eeltype==eelSWITCH) ? "switched":"shifted",
	      fr->rvdw_switch,fr->rvdw);
  } 

  if (fr->bBHAM && (fr->vdwtype == evdwSHIFT || fr->vdwtype == evdwSWITCH))
    fatal_error(0,"Switch/shift interaction not supported with Buckingham");
  
  if (fp)
    fprintf(fp,"Cut-off's:   NS: %g   Coulomb: %g   %s: %g\n",
	    fr->rlist,fr->rcoulomb,fr->bBHAM ? "BHAM":"LJ",fr->rvdw);
  
  if (ir->eDispCorr != edispcNO) {
    set_avcsix(fp,fr,mdatoms);
    set_avctwelve(fp,fr,mdatoms);
  }
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

  if (fr->bcoultab || fr->bvdwtab) {
    if (EEL_LR(fr->eeltype)) {
      bool bcoulsave,bvdwsave;
      /* generate extra tables for 1-4 interactions only
       * fake the forcerec so make_tables thinks it should
       * just create the non shifted version 
       */
      bcoulsave=fr->bcoultab;
      bvdwsave=fr->bvdwtab;
      fr->bcoultab=FALSE;
      fr->bvdwtab=FALSE;
      fr->rtab=ir->tabext;
      make_tables(fp,fr,MASTER(cr),tabfn);
      fr->bcoultab=bcoulsave;
      fr->bvdwtab=bvdwsave;
      fr->coulvdw14tab=fr->coulvdwtab;
      fr->coulvdwtab=NULL;
    }
    fr->rtab = fr->rlistlong + ir->tabext;
  }
  else if (fr->efep != efepNO) {
    if (fr->rlistlong == 0) {
      char *ptr,*envvar="FEP_TABLE_LENGTH";
      fr->rtab = 5;
      ptr = getenv(envvar);
      if (ptr) {
	sscanf(ptr,"%lf",&rtab);
	fr->rtab = rtab;
      }
      if (fp)
	fprintf(fp,"\nNote: Setting the free energy table length to %g nm\n"
		"      You can set this value with the environment variable %s"
		"\n\n",fr->rtab,envvar);
    } 
    else
      fr->rtab = fr->rlistlong + ir->tabext;
  } 
  else
    fr->rtab = ir->tabext;
  
  /* make tables for ordinary interactions */
  make_tables(fp,fr,MASTER(cr),tabfn);
  if(!(EEL_LR(fr->eeltype) && (fr->bcoultab || fr->bvdwtab)))
    fr->coulvdw14tab=fr->coulvdwtab;

  /* Copy the contents of the table to separate coulomb and LJ
   * tables too, to improve cache performance.
   */
  snew(fr->coultab,4*(fr->ntab+1));
  snew(fr->vdwtab,8*(fr->ntab+1));  
  for(i=0; i<=fr->ntab; i++) {
    for(j=0; j<4; j++) 
      fr->coultab[4*i+j]=fr->coulvdwtab[12*i+j];
    for(j=0; j<8; j++) 
      fr->vdwtab[8*i+j]=fr->coulvdwtab[12*i+4+j];
  }
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
  pr_int(fp,fr->ntab);
  if (fr->ntab > 0) {
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

void force(FILE       *fp,     int        step,
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
	   matrix lr_vir,       rvec       mu_tot,
	   real       qsum,     bool       bGatherOnly)
{
  int     i,nit;
  bool    bDoEpot;
  rvec    box_size;
  real    Vlr,Vcorr=0;
  
  /* Reset box */
  for(i=0; (i<DIM); i++)
    box_size[i]=box[i][i];
    
  bDoEpot=((fr->nmol > 0) && (fr->nstcalc > 0) && (mod(step,fr->nstcalc)==0));
  /* Reset epot... */
  if (bDoEpot) 
    for(i=0; (i<fr->nmol); i++)
      fr->mol_epot[i]=0.0;
  debug_gmx();
  
  /* Call the short range functions all in one go. */
  do_fnbf(fp,cr,fr,x,f,md,
	  fr->bBHAM ? grps->estat.ee[egBHAM] : grps->estat.ee[egLJ],
	  grps->estat.ee[egCOUL],box_size,nrnb,
	  lambda,&epot[F_DVDL],FALSE,-1);
  debug_gmx();

  if (debug) 
    pr_rvecs(debug,0,"fshift after SR",fr->fshift,SHIFTS);
  
  /* Shift the coordinates. Must be done before bonded forces and PPPM, 
   * but is also necessary for SHAKE and update, therefore it can NOT 
   * go when no bonded forces have to be evaluated.
   */
  if (debug && graph && 0)
    p_graph(debug,"DeBUGGGG",graph);
  
  /* Check whether we need to do bondeds */
  if (!bNBFonly) {
    if (graph) {
      shift_self(graph,box,x);
      if (TRICLINIC(box))
	inc_nrnb(nrnb,eNR_SHIFTX,2*graph->nnodes);
      else
	inc_nrnb(nrnb,eNR_SHIFTX,graph->nnodes);
    }
    debug_gmx();
  }
  
  if (EEL_LR(fr->eeltype)) {
    switch (fr->eeltype) {
    case eelPPPM:
      Vlr = do_pppm(fp,FALSE,x,fr->f_pme,md->chargeT,
		    box_size,fr->phi,cr,nsb,nrnb);
      break;
    case eelPME:
      Vlr = do_pme(fp,FALSE,ir,x,fr->f_pme,md->chargeT,
		   box,cr,nsb,nrnb,lr_vir,fr->ewaldcoeff,bGatherOnly);
      break;
    case eelEWALD:
      Vlr = do_ewald(fp,FALSE,ir,x,fr->f_pme,md->chargeT,
		     box_size,cr,nsb,lr_vir,fr->ewaldcoeff);
      break;
    default:
      Vlr = 0;
      fatal_error(0,"No such electrostatics method implemented %s",
		  eel_names[fr->eeltype]);
    }
    if(fr->bEwald)
      Vcorr =
	ewald_LRcorrection(fp,nsb,cr,fr,md->chargeT,excl,x,box,mu_tot,qsum,
			   ir->ewald_geometry,ir->epsilon_surface,lr_vir);
    else
      Vcorr = shift_LRcorrection(fp,nsb,cr,fr,md->chargeT,excl,x,TRUE,box,lr_vir);
    epot[F_LR] = Vlr + Vcorr;
    if (debug)
      fprintf(debug,"Vlr = %g, Vcorr = %g, Vlr_corr = %g\n",
	      Vlr,Vcorr,epot[F_LR]);
    if (debug) {
      pr_rvecs(debug,0,"lr_vir after corr",lr_vir,DIM);
      pr_rvecs(debug,0,"fshift after LR Corrections",fr->fshift,SHIFTS);
    }
  }
  debug_gmx();
  
  if (debug)    
    print_nrnb(debug,nrnb); 
  debug_gmx();
  
  if (!bNBFonly) {
    calc_bonds(fp,cr,mcr,
	       idef,x,f,fr,graph,epot,nrnb,box,lambda,md,
	       opts->ngener,grps->estat.ee[egLJ14],grps->estat.ee[egCOUL14],
	       fcd,step,fr->bSepDVDL && do_per_step(step,ir->nstlog));    
    debug_gmx();
  }
  if (debug) 
    pr_rvecs(debug,0,"fshift after bondeds",fr->fshift,SHIFTS);
  
  for(i=0; (i<F_EPOT); i++)
    if (i != F_DISRESVIOL && i != F_ORIRESDEV && i != F_DIHRESVIOL)
      epot[F_EPOT]+=epot[i];
}
