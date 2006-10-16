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
#include "domdec.h"
#include "pdbio.h"

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

static void write_constr_pdb(char *fn,char *title,t_atoms *atoms,
			     int start,int homenr,gmx_domdec_t *dd,
			     rvec x[],matrix box)
{
  char fname[STRLEN],format[STRLEN];
  FILE *out;
  int  i,ii,resnr;

  if (dd) {
    sprintf(fname,"%s_n%d.pdb",fn,dd->nodeid);
    start = 0;
    homenr = dd->nat_tot_con;
  } else {
    sprintf(fname,"%s.pdb",fn);
  }
  sprintf(format,"%s\n",pdbformat);

  out = ffopen(fname,"w");

  fprintf(out,"TITLE     %s\n",title);
  gmx_write_pdb_box(out,box);
  for(i=start; i<start+homenr; i++) {
    if (dd) {
      if (i >= dd->nat_home && i < dd->nat_tot_vsite)
	continue;
      ii = dd->gatindex[i];
    } else {
      ii = i;
    }
    resnr = atoms->atom[ii].resnr;
    fprintf(out,format,"ATOM",(ii+1)%100000,
	    *atoms->atomname[ii],*atoms->resname[resnr],' ',(resnr+1)%10000,
	    10*x[i][XX],10*x[i][YY],10*x[i][ZZ]);
  }
  fprintf(out,"TER\n");

  fclose(out);
}
			     
static void dump_confs(int step,t_atoms *atoms,
		       int start,int homenr,gmx_domdec_t *dd,
		       rvec x[],rvec xprime[],matrix box)
{
  char buf[256];
  
  sprintf(buf,"step%db",step);
  write_constr_pdb(buf,"initial coordinates",atoms,start,homenr,dd,x,box);
  sprintf(buf,"step%dc",step);
  write_constr_pdb(buf,"coordinates after constraining",atoms,start,homenr,dd,xprime,box);
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

static bool low_constrain(FILE *log,t_topology *top,t_ilist *settle,
			  t_inputrec *ir,
			  gmx_domdec_t *dd,
			  int step,t_mdatoms *md,int start,int homenr,
			  rvec *x,rvec *xprime,rvec *min_proj,matrix box,
			  real lambda,real *dvdlambda,
			  real invdt,rvec *v,tensor *vir,
			  t_nrnb *nrnb,bool bCoordinates,bool bInit)
{
  static bool      bFirst=TRUE;
  static int       nblocks=0;
  static int       *sblock=NULL;
  static t_lincsdata *lincsd=NULL;
  static bool      bDumpOnError = TRUE;
  
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
  int         settle_type;

  bOK = TRUE;
  if (bInit) {
    if (bFirst) {
      bDumpOnError = (getenv("NO_SHAKE_ERROR") == NULL);

      if (idef->il[F_SETTLE].nr > 0) {
	/* Check that we have only one settle type */
	settle_type=idef->il[F_SETTLE].iatoms[0];
	for (j=0; j<idef->il[F_SETTLE].nr; j+=2) {
	  if (idef->il[F_SETTLE].iatoms[j] != settle_type)
	    gmx_fatal(FARGS,"More than one settle type (%d and %d)",
		      settle_type,idef->il[F_SETTLE].iatoms[j]);
	}
	please_cite(log,"Miyamoto92a");
      }
    }

    if (dd == NULL) {
      ncons = idef->il[F_CONSTR].nr/3;
    } else {
      if (dd->constraints)
	ncons = dd->constraints->ncon;
      else
	ncons = 0;
    }
    if (ncons > 0 || (dd && dd->constraints)) {
      if (ir->eConstrAlg == estLINCS || !bCoordinates) {
	if (bFirst) {
	  please_cite(stdlog,"Hess97a");
	  snew(lincsd,1);
	}
	init_lincs(stdlog,&top->idef,start,homenr,EI_DYNAMICS(ir->eI),dd,
		   lincsd);
	set_lincs_matrix(lincsd,md->invmass);
	lincsd->matlam = lambda;
      } 
      else {
	if (bFirst)
	  please_cite(stdlog,"Ryckaert77a");
	else
	  gmx_fatal(FARGS,
		  "Constraint reinitialization not implemented for shake");
	
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
	iatom=idef->il[F_CONSTR].iatoms;
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
	
	iatom=idef->il[F_CONSTR].iatoms;
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
		    "top->blocks[ebSBLOCKS] does not match idef->il[F_CONSTR]");
	}
	sfree(sb);
	sfree(inv_sblock);
      }
    }
    
    bFirst = FALSE;
  } 
  else {
    /* !bInit */
    if (vir != NULL)
      clear_mat(rmdr);

    if (lincsd || nblocks > 0) {
      where();

      if (ir->eConstrAlg == estLINCS || !bCoordinates)
	bOK = constrain_lincs(stdlog,ir,step,lincsd,md,dd,
			      x,xprime,min_proj,box,lambda,dvdlambda,
			      invdt,v,vir!=NULL,rmdr,
			      bCoordinates,nrnb,bDumpOnError);
      else if (ir->eConstrAlg == estSHAKE)
	bOK = bshakef(stdlog,homenr,md->invmass,nblocks,sblock,idef,
		      ir,box,x,xprime,nrnb,lambda,dvdlambda,
		      invdt,v,vir!=NULL,rmdr,bDumpOnError);

      if (!bOK && bDumpOnError)
	fprintf(stdlog,"Constraint error in algorithm %s at step %d\n",
		eshake_names[ir->eConstrAlg],step);
    }
    if (settle->nr > 0) {
      int nsettle;
      real mO,mH,dOH,dHH;
      
      nsettle = settle->nr/2;
      mO   = md->massT[settle->iatoms[1]];
      mH   = md->massT[settle->iatoms[1]+1];
      dOH  = top->idef.iparams[settle->iatoms[0]].settle.doh;
      dHH  = top->idef.iparams[settle->iatoms[0]].settle.dhh;
      csettle(stdlog,nsettle,settle->iatoms,x[0],xprime[0],dOH,dHH,mO,mH,
	      invdt,v[0],vir!=NULL,rmdr,&error);
      inc_nrnb(nrnb,eNR_SETTLE,nsettle);
      if (v != NULL)
	inc_nrnb(nrnb,eNR_CONSTR_V,nsettle*3);
      if (vir != NULL)
	inc_nrnb(nrnb,eNR_CONSTR_VIR,nsettle*3);
      
      bOK = (error < 0);
      if (!bOK && bDumpOnError)
	fprintf(stdlog,"\nt = %.3f ps: Water molecule starting at atom %d can not be "
		"settled.\nCheck for bad contacts and/or reduce the timestep.",
		ir->init_t+step*ir->delta_t,
		glatnr(dd,settle->iatoms[error*2+1]));
    }
    if (vir != NULL) {
      hdt_2 = 0.5/(ir->delta_t*ir->delta_t);
      for(i=0; i<DIM; i++)
	for(j=0; j<DIM; j++)
	  (*vir)[i][j] = hdt_2*rmdr[i][j];
    }
    if (!bOK && bDumpOnError) 
      dump_confs(step,&(top->atoms),start,homenr,dd,x,xprime,box);
  }
  return bOK;
}

bool constrain(FILE *log,t_topology *top,t_ilist *settle,
	       t_inputrec *ir,
	       gmx_domdec_t *dd,
	       int step,t_mdatoms *md,int start,int homenr,
	       rvec *x,rvec *xprime,rvec *min_proj,matrix box,
	       real lambda,real *dvdlambda,
	       real dt,rvec *v,tensor *vir,
	       t_nrnb *nrnb,bool bCoordinates)
{
  return low_constrain(log,top,settle,ir,dd,
		       step,md,start,homenr,x,xprime,min_proj,box,
		       lambda,dvdlambda,
		       dt==0 ? 0 : 1/dt,v,vir,nrnb,bCoordinates,FALSE);
}

int count_constraints(t_topology *top,t_commrec *cr)
{
  int nc;
  
  if (DOMAINDECOMP(cr)) {
    nc = top->idef.il[F_SETTLE].nr*3/2;
    if (cr->dd->constraints)
      nc += cr->dd->constraints->ncon_global;
  } else {
    nc = top->idef.il[F_SETTLE].nr*3/2 + top->idef.il[F_CONSTR].nr/3;
    if (PAR(cr))
      gmx_sumi(1,&nc,cr);
  }

  return nc;
}

int init_constraints(FILE *log,t_topology *top,t_ilist *settle,t_inputrec *ir,
		      t_mdatoms *md,int start,int homenr,bool bOnlyCoords,
		      t_commrec *cr,gmx_domdec_t *dd)
{
  int count;

  low_constrain(log,top,settle,ir,dd,0,md,start,homenr,NULL,NULL,NULL,NULL,
		0,NULL,0,NULL,NULL,NULL,bOnlyCoords,TRUE);
  
  if (cr)
    count = count_constraints(top,cr);
  else
    count = -1;
  
  return count;
}
