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

#include "confio.h"
#include "constr.h"
#include "copyrite.h"
#include "invblock.h"
#include "main.h"
#include "mdrun.h"
#include "nrnb.h"
#include "smalloc.h"
#include "vec.h"
#include "physics.h"
#include "names.h"
#include "txtdump.h"

typedef struct {
  atom_id iatom[3];
  atom_id blocknr;
} t_sortblock;

static int pcount=0;

static int pcomp(const void *p1, const void *p2)
{
  int     db;
  atom_id min1,min2,max1,max2;
  t_sortblock *a1=(t_sortblock *)p1;
  t_sortblock *a2=(t_sortblock *)p2;

  pcount++;
  
  db=a1->blocknr-a2->blocknr;
  
  if (db != 0)
    return db;
    
  min1=min(a1->iatom[1],a1->iatom[2]);
  max1=max(a1->iatom[1],a1->iatom[2]);
  min2=min(a2->iatom[1],a2->iatom[2]);
  max2=max(a2->iatom[1],a2->iatom[2]);
  
  if (min1 == min2)
    return max1-max2;
  else
    return min1-min2;
}

static int icomp(const void *p1, const void *p2)
{
  atom_id *a1=(atom_id *)p1;
  atom_id *a2=(atom_id *)p2;

  return (*a1)-(*a2);
}

static void dump_confs(int step,t_atoms *atoms,
		       rvec x[],rvec xprime[],matrix box)
{
  char buf[256];
  
  sprintf(buf,"step%d.pdb",step-1);
  write_sto_conf(buf,"one step before crash",atoms,x,NULL,box);
  sprintf(buf,"step%d.pdb",step);
  write_sto_conf(buf,"crashed",atoms,xprime,NULL,box);
  fprintf(stdlog,"Wrote pdb files with previous and current coordinates\n");
  fprintf(stderr,"Wrote pdb files with previous and current coordinates\n");
}

static void pr_sortblock(FILE *fp,char *title,int nsb,t_sortblock sb[])
{
  int i;
  
  fprintf(fp,"%s\n",title);
  for(i=0; (i<nsb); i++)
    fprintf(fp,"i: %5d, iatom: (%5d %5d %5d), blocknr: %5d\n",
	    i,sb[i].iatom[0],sb[i].iatom[1],sb[i].iatom[2],
	    sb[i].blocknr);
}

static bool low_constrain(FILE *log,t_topology *top,t_inputrec *ir,
			  int step,t_mdatoms *md,int start,int homenr,
			  rvec *x,rvec *xprime,rvec *min_proj,matrix box,
			  real lambda,real *dvdlambda,tensor *vir,
			  t_nrnb *nrnb,bool bCoordinates,bool bInit)
{
  static bool      bFirst=TRUE;
  static int       nblocks=0;
  static int       *sblock=NULL;
  static int       nsettle,nsettle_alloc=0,settle_type;
  static int       *owptr=NULL;
  static t_lincsdata *lincsd=NULL;
  static bool      bDumpOnError = TRUE;
  
  char        buf[STRLEN];
  bool        bOK;
  t_sortblock *sb;
  t_block     *blocks=&(top->blocks[ebSBLOCKS]);
  t_idef      *idef=&(top->idef);
  t_iatom     *iatom;
  atom_id     *inv_sblock;
  int         i,j,m,bnr;
  int         ncons,bstart,error;
  tensor      rmdr;
  real        hdt_2;
  
  bOK = TRUE;
  if (bInit) {
    /* Output variables, initiate them right away */

    if (bFirst) {
      if ((ir->etc==etcBERENDSEN) || (ir->epc==epcBERENDSEN))
	please_cite(log,"Berendsen84a");
      
      bDumpOnError = (getenv("NO_SHAKE_ERROR") == NULL);
    }
    
    /* Put the oxygen atoms in the owptr array */
    nsettle=idef->il[F_SETTLE].nr/2;
    if (nsettle > 0) {
      if (nsettle > nsettle_alloc) {
	nsettle_alloc = over_alloc(nsettle);
	srenew(owptr,nsettle_alloc);
      }
      settle_type=idef->il[F_SETTLE].iatoms[0];
      for (j=0; (j<idef->il[F_SETTLE].nr); j+=2) {
	if (idef->il[F_SETTLE].iatoms[j] != settle_type)
	  gmx_fatal(FARGS,"More than one settle type (%d and %d)",
		      settle_type,idef->il[F_SETTLE].iatoms[j]);
	owptr[j/2]=idef->il[F_SETTLE].iatoms[j+1];
#ifdef DEBUG
	fprintf(log,"owptr[%d]=%d\n",j/2,owptr[j/2]);
#endif
      }
      /* We used to free this memory, but ED sampling needs it later on 
       *  sfree(idef->il[F_SETTLE].iatoms);
       */
      
      if (bFirst)
	please_cite(log,"Miyamoto92a");
    }
    
    ncons=idef->il[F_SHAKE].nr/3;
    if (ncons > 0) {
      if (!bFirst)
	gmx_fatal(FARGS,
		  "Constraint reinitialization only implemented for settle");

      bstart=(idef->nodeid > 0) ? blocks->multinr[idef->nodeid-1] : 0;
      nblocks=blocks->multinr[idef->nodeid] - bstart;
      if (debug) 
	fprintf(debug,"ncons: %d, bstart: %d, nblocks: %d\n",
		ncons,bstart,nblocks);
      
      /* Calculate block number for each atom */
      inv_sblock=make_invblock(blocks,md->nr);
      
      /* Store the block number in temp array and
       * sort the constraints in order of the sblock number 
       * and the atom numbers, really sorting a segment of the array!
       */
#ifdef DEBUGIDEF 
      pr_idef(stdlog,0,"Before Sort",idef);
#endif
      iatom=idef->il[F_SHAKE].iatoms;
      snew(sb,ncons);
      for(i=0; (i<ncons); i++,iatom+=3) {
	for(m=0; (m<3); m++)
	  sb[i].iatom[m]=iatom[m];
	sb[i].blocknr=inv_sblock[iatom[1]];
      }
      
      /* Now sort the blocks */
      if (debug) {
	pr_sortblock(debug,"Before sorting",ncons,sb);
	fprintf(debug,"Going to sort constraints\n");
      }
      
      qsort(sb,ncons,(size_t)sizeof(*sb),pcomp);
      
      if (debug) {
	fprintf(debug,"I used %d calls to pcomp\n",pcount);
	pr_sortblock(debug,"After sorting",ncons,sb);
      }
      
      iatom=idef->il[F_SHAKE].iatoms;
      for(i=0; (i<ncons); i++,iatom+=3) 
	for(m=0; (m<DIM); m++)
	  iatom[m]=sb[i].iatom[m];
#ifdef DEBUGIDEF
      pr_idef(stdlog,0,"After Sort",idef);
#endif
      
      j=0;
      snew(sblock,nblocks+1);
      bnr=-2;
      for(i=0; (i<ncons); i++) {
	if (sb[i].blocknr != bnr) {
	  bnr=sb[i].blocknr;
	  sblock[j++]=3*i;
	}
      }
      /* Last block... */
      sblock[j++]=3*ncons;
      
      if (j != (nblocks+1)) {
	fprintf(log,"bstart: %d\n",bstart);
	fprintf(log,"j: %d, nblocks: %d, ncons: %d\n",
		j,nblocks,ncons);
	for(i=0; (i<ncons); i++)
	  fprintf(log,"i: %5d  sb[i].blocknr: %5u\n",i,sb[i].blocknr);
	for(j=0; (j<=nblocks); j++)
	  fprintf(log,"sblock[%3d]=%5d\n",j,(int) sblock[j]);
	gmx_fatal(FARGS,"DEATH HORROR: "
		    "top->blocks[ebSBLOCKS] does not match idef->il[F_SHAKE]");
      }
      sfree(sb);
      sfree(inv_sblock);
    }
    
    if (idef->il[F_SHAKE].nr) {
      if (ir->eConstrAlg == estLINCS || !bCoordinates) {
	if (bFirst)
	  please_cite(stdlog,"Hess97a");
	lincsd = init_lincs(stdlog,&top->idef,start,homenr,
			    EI_DYNAMICS(ir->eI));
	set_lincs_matrix(lincsd,md->invmass);
	lincsd->matlam = lambda;
      } 
      else {
	if (bFirst)
	  please_cite(stdlog,"Ryckaert77a");
      }
    }
    
    bFirst = FALSE;
  } 
  else {
    /* !bInit */
    if (vir != NULL)
      clear_mat(rmdr);

    if (nblocks > 0) {
      where();

      if (ir->eConstrAlg == estLINCS || !bCoordinates)
	bOK = constrain_lincs(stdlog,ir,step,lincsd,md,
			      x,xprime,min_proj,box,lambda,dvdlambda,
			      vir!=NULL,rmdr,
			      bCoordinates,nrnb,bDumpOnError);
      else if (ir->eConstrAlg == estSHAKE)
	bOK = bshakef(stdlog,homenr,md->invmass,nblocks,sblock,idef,
		      ir,box,x,xprime,nrnb,lambda,dvdlambda,
		      vir!=NULL,rmdr,bDumpOnError);

      if (!bOK && bDumpOnError)
	fprintf(stdlog,"Constraint error in algorithm %s at step %d\n",
		eshake_names[ir->eConstrAlg],step);
    }
    if (nsettle > 0) {
      int  ow1;
      real mO,mH,dOH,dHH;
      
      ow1  = owptr[0];
      mO   = md->massT[ow1];
      mH   = md->massT[ow1+1];
      dOH  = top->idef.iparams[settle_type].settle.doh;
      dHH  = top->idef.iparams[settle_type].settle.dhh;
      csettle(stdlog,nsettle,owptr,x[0],xprime[0],dOH,dHH,mO,mH,
	      vir!=NULL,rmdr,&error);
      inc_nrnb(nrnb,eNR_SETTLE,nsettle);
      if (vir != NULL)
	inc_nrnb(nrnb,eNR_CONSTR_VIR,nsettle*3);
      
      bOK = (error < 0);
      if (!bOK && bDumpOnError)
	fprintf(stdlog,"\nt = %.3f ps: Water molecule starting at atom %d can not be "
		"settled.\nCheck for bad contacts and/or reduce the timestep.",
		ir->init_t+step*ir->delta_t,owptr[error]+1);
    }
    if (vir != NULL) {
      hdt_2 = 0.5/(ir->delta_t*ir->delta_t);
      for(i=0; i<DIM; i++)
	for(j=0; j<DIM; j++)
	  (*vir)[i][j] = hdt_2*rmdr[i][j];
    }
    if (!bOK && bDumpOnError) 
      dump_confs(step,&(top->atoms),x,xprime,box);
  }
  return bOK;
}

bool constrain(FILE *log,t_topology *top,t_inputrec *ir,int step,
	       t_mdatoms *md,int start,int homenr,
	       rvec *x,rvec *xprime,rvec *min_proj,matrix box,
	       real lambda,real *dvdlambda,tensor *vir,
	       t_nrnb *nrnb,bool bCoordinates)
{
  return low_constrain(log,top,ir,step,md,start,homenr,x,xprime,min_proj,box,
		       lambda,dvdlambda,vir,nrnb,bCoordinates,FALSE);
}

int count_constraints(t_topology *top,t_commrec *cr)
{
  int nc;
  
  nc = top->idef.il[F_SETTLE].nr*3/2 + top->idef.il[F_SHAKE].nr/3;
  if (PAR(cr))
    gmx_sumi(1,&nc,cr);
  
  return nc;
}

int init_constraints(FILE *log,t_topology *top,t_inputrec *ir,
		      t_mdatoms *md,int start,int homenr,bool bOnlyCoords,
		      t_commrec *cr)
{
  int count;

  low_constrain(log,top,ir,0,md,start,homenr,NULL,NULL,NULL,NULL,
		0,NULL,NULL,NULL,bOnlyCoords,TRUE);

  if (cr)
    count = count_constraints(top,cr);
  else
    count = -1;
  
  return count;
}
