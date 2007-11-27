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
#include "partdec.h"

typedef struct gmx_constr {
  int             nflexcon;     /* The number of flexible constraints */
  t_block         at2con;       /* A list of atoms to constraints     */
  bool            bInterCGcons; /* Are there inter charge-group con's */
  gmx_lincsdata_t lincsd;       /* LINCS data                         */
  int             nblocks;      /* The number of SHAKE blocks         */
  int             *sblock;      /* The SHAKE blocks                   */
  int             maxwarn;      /* The maximum number of warnings     */
  int             warncount_lincs;
  int             warncount_settle;
} t_gmx_constr;

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

int n_flexible_constraints(struct gmx_constr *constr)
{
  int nflexcon;

  if (constr)
    nflexcon = constr->nflexcon;
  else
    nflexcon = 0;

  return nflexcon;
}

void too_many_constraint_warnings(int eConstrAlg,int warncount)
{
  char *abort="- aborting to avoid logfile runaway.\n"
    "This normally happens when your system is not sufficiently equilibrated,"
    "or if you are changing lambda too fast in free energy simulations.\n";
  
  gmx_fatal(FARGS,
	    "Too many %s warnings (%d)\n"
	    "If you know what you are doing you can %s"
	    "set the environment variable GMX_MAXCONSTRWARN to -1,\n"
	    "but normally it is better to fix the problem",
	    (eConstrAlg == estLINCS) ? "LINCS" : "SETTLE",warncount,
	    (eConstrAlg == estLINCS) ?
	    "adjust the lincs warning threshold in your mdp file\nor " : "\n");
}

static void write_constr_pdb(char *fn,char *title,t_atoms *atoms,
			     int start,int homenr,gmx_domdec_t *dd,
			     rvec x[],matrix box)
{
  char fname[STRLEN],format[STRLEN];
  FILE *out;
  int  dd_ac0=0,dd_ac1=0,i,ii,resnr;

  if (dd) {
    sprintf(fname,"%s_n%d.pdb",fn,dd->sim_nodeid);
    dd_get_constraint_range(dd,&dd_ac0,&dd_ac1);
    start = 0;
    homenr = dd_ac1;
  } else {
    sprintf(fname,"%s.pdb",fn);
  }
  sprintf(format,"%s\n",pdbformat);

  out = ffopen(fname,"w");

  fprintf(out,"TITLE     %s\n",title);
  gmx_write_pdb_box(out,box);
  for(i=start; i<start+homenr; i++) {
    if (dd) {
      if (i >= dd->nat_home && i < dd_ac0)
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
			     
static void dump_confs(FILE *fplog,int step,t_atoms *atoms,
		       int start,int homenr,gmx_domdec_t *dd,
		       rvec x[],rvec xprime[],matrix box)
{
  char buf[256];
  
  sprintf(buf,"step%db",step);
  write_constr_pdb(buf,"initial coordinates",atoms,start,homenr,dd,x,box);
  sprintf(buf,"step%dc",step);
  write_constr_pdb(buf,"coordinates after constraining",atoms,start,homenr,dd,xprime,box);
  if (fplog)
    fprintf(fplog,"Wrote pdb files with previous and current coordinates\n");
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

bool constrain(FILE *fplog,bool bLog,bool bEner,
	       struct gmx_constr *constr,
	       t_topology *top,t_inputrec *ir,
	       gmx_domdec_t *dd,
	       int step,t_mdatoms *md,
	       rvec *x,rvec *xprime,rvec *min_proj,matrix box,
	       real lambda,real *dvdlambda,
	       real dt,rvec *v,tensor *vir,
	       t_nrnb *nrnb,int econq)
{
  bool    bOK;
  int     start,homenr;
  int     i,j;
  int     ncons,error;
  tensor  rmdr;
  real    invdt,hdt_2;
  t_ilist *settle;
  int     nsettle;
  real    mO,mH,dOH,dHH;
    
  bOK = TRUE;

  start  = md->start;
  homenr = md->homenr;
  if (ir->delta_t == 0)
    invdt = 0;
  else
    invdt  = 1/ir->delta_t;

  if (vir != NULL)
    clear_mat(rmdr);
    
  where();
  if (constr->lincsd) {
    bOK = constrain_lincs(fplog,bLog,bEner,ir,step,constr->lincsd,md,dd,
			  x,xprime,min_proj,box,lambda,dvdlambda,
			  invdt,v,vir!=NULL,rmdr,
			  econq,nrnb,
			  constr->maxwarn,&constr->warncount_lincs);
    if (!bOK && constr->maxwarn >= 0 && fplog)
      fprintf(fplog,"Constraint error in algorithm %s at step %d\n",
	      eshake_names[estLINCS],step);
  }
  
  if (constr->nblocks > 0) {
    if (econq != econqCoord)
      gmx_fatal(FARGS,"Internal error, SHAKE called for constraining something else than coordinates");

    bOK = bshakef(fplog,homenr,md->invmass,constr->nblocks,constr->sblock,
		  &top->idef,ir,box,x,xprime,nrnb,lambda,dvdlambda,
		  invdt,v,vir!=NULL,rmdr,constr->maxwarn>=0);
    if (!bOK && constr->maxwarn >= 0 && fplog)
      fprintf(fplog,"Constraint error in algorithm %s at step %d\n",
	      eshake_names[estSHAKE],step);
  }
  
  settle  = &top->idef.il[F_SETTLE];
  if (settle->nr > 0) {
    if (econq != econqCoord)
	gmx_fatal(FARGS,"For this system also velocities and/or forces need to be constrained, this can not be done with SETTLE");

    nsettle = settle->nr/2;
    mO   = md->massT[settle->iatoms[1]];
    mH   = md->massT[settle->iatoms[1]+1];
    dOH  = top->idef.iparams[settle->iatoms[0]].settle.doh;
    dHH  = top->idef.iparams[settle->iatoms[0]].settle.dhh;
    csettle(fplog,nsettle,settle->iatoms,x[0],xprime[0],dOH,dHH,mO,mH,
	    invdt,v[0],vir!=NULL,rmdr,&error);
    inc_nrnb(nrnb,eNR_SETTLE,nsettle);
    if (v != NULL)
      inc_nrnb(nrnb,eNR_CONSTR_V,nsettle*3);
    if (vir != NULL)
      inc_nrnb(nrnb,eNR_CONSTR_VIR,nsettle*3);
    
    bOK = (error < 0);
    if (!bOK && constr->maxwarn >= 0) {
      char buf[256];
      sprintf(buf,
	      "\nt = %.3f ps: Water molecule starting at atom %d can not be "
	      "settled.\nCheck for bad contacts and/or reduce the timestep.\n",
	      ir->init_t+step*ir->delta_t,
	      glatnr(dd,settle->iatoms[error*2+1]));
      if (fplog)
	fprintf(fplog,"%s",buf);
      fprintf(stderr,"%s",buf);
      constr->warncount_settle++;
      if (constr->warncount_settle > constr->maxwarn)
	too_many_constraint_warnings(-1,constr->warncount_settle);
    }
  }

  if (vir != NULL) {
    hdt_2 = 0.5/(ir->delta_t*ir->delta_t);
    for(i=0; i<DIM; i++)
      for(j=0; j<DIM; j++)
	(*vir)[i][j] = hdt_2*rmdr[i][j];
  }

  if (!bOK && constr->maxwarn >= 0) 
    dump_confs(fplog,step,&(top->atoms),start,homenr,dd,x,xprime,box);

  return bOK;
}

real *constr_rmsd_data(struct gmx_constr *constr)
{
  if (constr->lincsd)
    return lincs_rmsd_data(constr->lincsd);
  else
    return NULL;
}

real constr_rmsd(struct gmx_constr *constr,bool bSD2)
{
  if (constr->lincsd)
    return lincs_rmsd(constr->lincsd,bSD2);
  else
    return 0;
}

void set_constraints(struct gmx_constr *constr,
		     t_topology *top,t_inputrec *ir,
		     t_mdatoms *md,gmx_domdec_t *dd)
{
  int  i,j,m,ncons;
  t_idef *idef=&(top->idef);
  int  bstart,bnr;
  t_sortblock *sb;
  t_block     *blocks=&(top->blocks[ebSBLOCKS]);
  t_iatom     *iatom;
  atom_id     *inv_sblock;
  int  settle_type;

  if (dd == NULL) {
    ncons = idef->il[F_CONSTR].nr/3;
  } else {
    if (dd->constraints)
      ncons = dd->constraints->ncon;
    else
      ncons = 0;
  }
  if (ncons > 0 || (dd && dd->constraints)) {
    if (ir->eConstrAlg == estLINCS) {
      set_lincs(&top->idef,md->start,md->homenr,&constr->at2con,
		EI_DYNAMICS(ir->eI),dd,constr->lincsd);
      if (dd == NULL) {
	/* The atom to constraints list is no longer needed */
	done_block(&constr->at2con);
      }
      set_lincs_matrix(constr->lincsd,md->invmass,md->lambda);
    } 
    if (ir->eConstrAlg == estSHAKE) {
      if (constr->nblocks > 0)
	gmx_fatal(FARGS,
		  "Constraint reinitialization not implemented for shake");
      
      /*
	bstart=(idef->nodeid > 0) ? blocks->multinr[idef->nodeid-1] : 0;
	nblocks=blocks->multinr[idef->nodeid] - bstart;
      */
      bstart  = 0;
      constr->nblocks = blocks->nr;
      if (debug) 
	fprintf(debug,"ncons: %d, bstart: %d, nblocks: %d\n",
		ncons,bstart,constr->nblocks);
      
      /* Calculate block number for each atom */
      inv_sblock = make_invblock(blocks,md->nr);
      
      /* Store the block number in temp array and
       * sort the constraints in order of the sblock number 
       * and the atom numbers, really sorting a segment of the array!
       */
#ifdef DEBUGIDEF 
      pr_idef(fplog,0,"Before Sort",idef);
#endif
      iatom=idef->il[F_CONSTR].iatoms;
      snew(sb,ncons);
      for(i=0; (i<ncons); i++,iatom+=3) {
	for(m=0; (m<3); m++)
	  sb[i].iatom[m] = iatom[m];
	sb[i].blocknr = inv_sblock[iatom[1]];
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
      pr_idef(fplog,0,"After Sort",idef);
#endif
      
      j=0;
      snew(constr->sblock,constr->nblocks+1);
      bnr=-2;
      for(i=0; (i<ncons); i++) {
	if (sb[i].blocknr != bnr) {
	  bnr=sb[i].blocknr;
	  constr->sblock[j++]=3*i;
	}
      }
      /* Last block... */
      constr->sblock[j++] = 3*ncons;
      
      if (j != (constr->nblocks+1)) {
	fprintf(stderr,"bstart: %d\n",bstart);
	fprintf(stderr,"j: %d, nblocks: %d, ncons: %d\n",
		j,constr->nblocks,ncons);
	for(i=0; (i<ncons); i++)
	  fprintf(stderr,"i: %5d  sb[i].blocknr: %5u\n",i,sb[i].blocknr);
	for(j=0; (j<=constr->nblocks); j++)
	  fprintf(stderr,"sblock[%3d]=%5d\n",j,(int)constr->sblock[j]);
	gmx_fatal(FARGS,"DEATH HORROR: "
		  "top->blocks[ebSBLOCKS] does not match idef->il[F_CONSTR]");
      }
      sfree(sb);
      sfree(inv_sblock);
    }
  }
}

static bool interCGconstraints(t_block *cgs,t_block *at2con,t_iatom *iatoms)
{
  bool bInterCGConstraints;
  int  cg,a0,a1,a,c,a2;

  bInterCGConstraints = FALSE;
  for(cg=0; cg<cgs->nr && !bInterCGConstraints; cg++) {
    a0 = cgs->index[cg];
    a1 = cgs->index[cg+1];
    for(a=a0; a<a1; a++) {
      for(c=at2con->index[a]; c<at2con->index[a+1]; c++) {
	/* Single all constraints are stored for both atoms
	 * we only need to check the second atom.
	 */
	a2 = iatoms[3*at2con->a[c]+2];
	if (a2 < a0 || a2 >= a1)
	  bInterCGConstraints = TRUE;
      }
    }
  }

  return bInterCGConstraints;
}

static t_block make_at2con(int start,int natoms,t_idef *idef,
			   bool bDynamics,int *nflexiblecons)
{
  int *count,ncon,con,nflexcon,i,a;
  t_iatom *ia;
  t_block at2con;
  bool bFlexCon;
  
  snew(count,natoms);
  ncon = idef->il[F_CONSTR].nr/3;
  ia   = idef->il[F_CONSTR].iatoms;
  nflexcon = 0;
  for(con=0; con<ncon; con++) {
    bFlexCon = (idef->iparams[ia[0]].constr.dA == 0 &&
		idef->iparams[ia[0]].constr.dB == 0);
    if (bFlexCon)
      nflexcon++;
    if (bDynamics || !bFlexCon) {
      for(i=1; i<3; i++) {
	a = ia[i] - start;
	count[a]++;
      }
    }
    ia += 3;
  }
  *nflexiblecons = nflexcon;

  at2con.nr = natoms;
  at2con.nalloc_index = at2con.nr+1;
  snew(at2con.index,at2con.nalloc_index);
  at2con.index[0] = 0;
  for(a=0; a<natoms; a++) {
    at2con.index[a+1] = at2con.index[a] + count[a];
    count[a] = 0;
  }
  at2con.nra = at2con.index[natoms];
  at2con.nalloc_a = at2con.nra;
  snew(at2con.a,at2con.nalloc_a);

  ia   = idef->il[F_CONSTR].iatoms;
  for(con=0; con<ncon; con++) {
    bFlexCon = (idef->iparams[ia[0]].constr.dA == 0 &&
		idef->iparams[ia[0]].constr.dB == 0);
    if (bDynamics || !bFlexCon) {
      for(i=1; i<3; i++) {
	a = ia[i] - start;
	at2con.a[at2con.index[a]+count[a]++] = con;
      }
    }
    ia += 3;
  }
  
  sfree(count);

  return at2con;
}

gmx_constr_t init_constraints(FILE *fplog,t_commrec *cr,
			      t_topology *top,t_inputrec *ir)
{
  int  nc[2],settle_type,j,start,natoms;
  struct gmx_constr *constr;
  char *env;

  /* The number of normal constraints */
  nc[0] = top->idef.il[F_CONSTR].nr/3;
  /* The number of constraints handled by SETTLE */
  nc[1] = top->idef.il[F_SETTLE].nr*3/2;
  if (PARTDECOMP(cr))
    gmx_sumi(2,nc,cr);

  if (nc[0]+nc[1] > 0) {
    snew(constr,1);

    if (nc[0] == 0) {
      constr->nflexcon = 0;
    } else {
      if (PARTDECOMP(cr)) {
	pd_at_range(cr,&start,&j);
	natoms = start + j;
      } else {
	start  = 0;
	natoms = top->atoms.nr;
      }

      constr->at2con = make_at2con(start,natoms,&top->idef,
				   EI_DYNAMICS(ir->eI),&constr->nflexcon);

      if (DOMAINDECOMP(cr)) {
	constr->bInterCGcons =
	  interCGconstraints(&top->blocks[ebCGS],&constr->at2con,
			     top->idef.il[F_CONSTR].iatoms);
      }

      if (PARTDECOMP(cr))
	gmx_sumi(1,&constr->nflexcon,cr);
      
      if (constr->nflexcon > 0) {
	if (fplog) {
	  fprintf(fplog,"There are %d flexible constraints\n",
		  constr->nflexcon);
	  if (ir->fc_stepsize == 0) {
	    fprintf(fplog,"\n"
		    "WARNING: step size for flexible constraining = 0\n"
		    "         All flexible constraints will be rigid.\n"
		    "         Will try to keep all flexible constraints at their original length,\n"
		    "         but the lengths may exhibit some drift.\n\n");
	    constr->nflexcon = 0;
	  }
	}
	if (constr->nflexcon > 0)
	  please_cite(fplog,"Hess2002");
      }
      
      if (ir->eConstrAlg == estLINCS) {
	constr->lincsd = init_lincs(fplog,&top->idef,
				    constr->nflexcon,&constr->at2con,
				    DOMAINDECOMP(cr) && constr->bInterCGcons,
				    ir->nLincsIter,ir->nProjOrder);
      }
      
      if (ir->eConstrAlg == estSHAKE) {
	if (constr->nflexcon)
	  gmx_fatal(FARGS,"For this system also velocities and/or forces need to be constrained, this can not be done with SHAKE, you should select LINCS");
	please_cite(fplog,"Ryckaert77a");
      }
    }
    
    if (nc[1] > 0) {
      please_cite(fplog,"Miyamoto92a");
      if (top->idef.il[F_SETTLE].nr > 0) {
	/* Check that we have only one settle type */
	settle_type=top->idef.il[F_SETTLE].iatoms[0];
	for (j=0; j<top->idef.il[F_SETTLE].nr; j+=2) {
	  if (top->idef.il[F_SETTLE].iatoms[j] != settle_type)
	    gmx_fatal(FARGS,"More than one settle type (%d and %d)",
		    settle_type,top->idef.il[F_SETTLE].iatoms[j]);
	}
      }
    }

    constr->maxwarn = 999;
    env = getenv("GMX_MAXCONSTRWARN");
    if (env) {
      constr->maxwarn = 0;
      sscanf(env,"%d",&constr->maxwarn);
      if (fplog)
	fprintf(fplog,
		"Setting the maximum number of constraint warnings to %d\n",
		constr->maxwarn);
      if (MASTER(cr))
	fprintf(stderr,
		"Setting the maximum number of constraint warnings to %d\n",
		constr->maxwarn);
    }
    if (constr->maxwarn < 0 && fplog)
      fprintf(fplog,"maxwarn < 0, will not stop on constraint errors\n");
    constr->warncount_lincs  = 0;
    constr->warncount_settle = 0;
  } else {
    constr = NULL;
  }
  
  return constr;
}

t_block *atom2constraints(gmx_constr_t constr)
{
  return &constr->at2con;
}

bool inter_charge_group_constraints(gmx_constr_t constr)
{
  return constr->bInterCGcons;
}
