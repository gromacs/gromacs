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
#include "callf77.h"
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

static int int_comp(const void *a,const void *b)
{
  return (*(int *)a) - (*(int *)b);
}

static void init_lincs(FILE *log,t_topology *top,t_inputrec *ir,
		       t_mdatoms *md,int start,int homenr,
		       int *nrtot,
		       rvec **r,int **bla1,int **bla2,int **blnr,int **blbnb,
		       real **bllen,real **blc,real **blcc,real **blm,
		       real **tmp1,real **tmp2,real **tmp3,
		       real **lincslam,real **bllen0,real **ddist)
{
  t_idef      *idef=&(top->idef);
  t_iatom     *iatom;
  int         i,j,k,n,b1,b,cen;
  int         ncons,nZeroLen;
  int         type,a1,a2,b2,nr,n1,n2,nc4;
  real        len=0,len1,sign;
  real        im1,im2;
  int         **at_c,*at_cn,*at_cm;
  
  ncons  = idef->il[F_SHAKE].nr/3;
  *nrtot = 0;
  
  if (ncons > 0) {

    iatom=idef->il[F_SHAKE].iatoms;

    /* Make atom-constraint connection list for temporary use */
    snew(at_c,homenr);
    snew(at_cn,homenr);
    snew(at_cm,homenr);

    for(i=0; i<ncons; i++) {
      a1=iatom[3*i+1]-start;
      a2=iatom[3*i+2]-start;
      if (at_cn[a1] >= at_cm[a1]) {
	at_cm[a1] += 4;
	srenew(at_c[a1],at_cm[a1]);
      }
      at_c[a1][at_cn[a1]] = i;
      at_cn[a1]++;
      if (at_cn[a2] >= at_cm[a2]) {
	at_cm[a2] += 4;
	srenew(at_c[a2],at_cm[a2]);
      }
      at_c[a2][at_cn[a2]] = i;
      at_cn[a2]++;
    }
    sfree(at_cm);
    
    for(i=0; i<ncons; i++) {
      a1=iatom[3*i+1]-start;
      a2=iatom[3*i+2]-start;
      *nrtot += at_cn[a1] + at_cn[a2] - 2;
    }      

    snew(*r,ncons);
    snew(*bla1,ncons);
    snew(*bla2,ncons);
    snew(*blnr,ncons+1);
    snew(*bllen,ncons);
    snew(*blc,ncons);
    snew(*tmp1,ncons);
    snew(*tmp2,ncons);
    snew(*tmp3,ncons);
    snew(*lincslam,ncons);
    snew(*bllen0,ncons);
    snew(*ddist,ncons);
    snew(*blbnb,*nrtot);
    snew(*blcc,*nrtot);
    snew(*blm,*nrtot);
    
    /* Make constraint-neighbor list */
    (*blnr)[0] = 0;
    nZeroLen = 0;
    for(i=0; (i<ncons); i++) {
      j=3*i;
      a1=iatom[j+1];
      a2=iatom[j+2];
      /* (*blnr)[i+1] = (*blnr)[i] + at_cn[a1] + at_cn[a2] - 2; */
      type=iatom[j];
      len =idef->iparams[type].shake.dA;
      len1=idef->iparams[type].shake.dB;
      if ((len==0 || len1==0) && len1!=len)
	gmx_fatal(FARGS,"It is not allowed to have a constraint length "
		    "zero and non-zero in the A and B topology");
      if (len == 0)
	nZeroLen++;
      (*bla1)[i]=a1;
      (*bla2)[i]=a2;
      (*bllen)[i]=len;
      (*bllen0)[i]=len;
      (*ddist)[i]=len1-len;
      im1=md->invmass[a1];
      im2=md->invmass[a2];
      (*blc)[i]=invsqrt(im1+im2);
      /* Construct the constraint connection matrix blbnb */
      (*blnr)[i+1]=(*blnr)[i];
      for(k=0; k<at_cn[a1-start]; k++)
	if (at_c[a1-start][k] != i)
	  (*blbnb)[((*blnr)[i+1])++]=at_c[a1-start][k];
      for(k=0; k<at_cn[a2-start]; k++)
	if (at_c[a2-start][k] != i)
	  (*blbnb)[((*blnr)[i+1])++]=at_c[a2-start][k];
      /* Order the blbnb matrix to optimize memory access */
      qsort(&((*blbnb)[(*blnr)[i]]),(*blnr)[i+1]-(*blnr)[i],
            sizeof((*blbnb)[0]),int_comp);
    }

    sfree(at_cn);
    for(i=0; i<homenr; i++)
      sfree(at_c[i]);
    sfree(at_c);
    
    fprintf(log,"\nInitializing LINear Constraint Solver\n");
    fprintf(log,"  number of constraints is %d\n",ncons);
    fprintf(log,"  average number of constraints coupled to one constraint is %.1f\n",
	    (real)(*nrtot)/ncons);
    if (nZeroLen)
      fprintf(log,"  found %d constraints with zero length\n",nZeroLen);
    fprintf(log,"\n");
    fflush(log);

    /* Construct the coupling coefficient matrix blcc */
    for(b=0; (b<ncons); b++) {
      i=(*bla1)[b];
      j=(*bla2)[b];
      for(n=(*blnr)[b]; (n<(*blnr)[b+1]);n++) {
	k = (*blbnb)[n];
	if (i==(*bla1)[k] || j==(*bla2)[k])
	  sign=-1;
	else
	  sign=1;
	if (i==(*bla1)[k] || i==(*bla2)[k])
	  cen=i;
	else
	  cen=j;
	(*blcc)[n]=sign*md->invmass[cen]*(*blc)[b]*(*blc)[k];
      }
    }
#ifdef DEBUG
    for(i=0; i<ncons; i++) {
      fprintf(log,"%d  %d %d  %g  %g  %d\n",i,(*bla1)[i],(*bla2)[i],
      (*bllen)[i],(*blc)[i],(*blnr)[i+1]-(*blnr)[i]);
      for(n=(*blnr)[i]; (n<(*blnr)[i+1]);n++) {
	k = (*blbnb)[n];
	fprintf(log,"  %d  %g\n",k,(*blcc)[n]);
      }
    }
#endif
  }
}

static bool constrain_lincs(FILE *log,t_topology *top,t_inputrec *ir,
			    int step,t_mdatoms *md,int start,int homenr,
			    int *nbl,int **sbl,
			    rvec *x,rvec *xprime,rvec *min_proj,matrix box,
			    real lambda,real *dvdlambda,bool bCoordinates,
			    bool bInit,t_nrnb *nrnb,bool bDumpOnError)
{
  static int       *bla1,*bla2,*blnr,*blbnb,nrtot=0;
  static rvec      *r;
  static real      *bllen,*blc,*blcc,*blm,*tmp1,*tmp2,*tmp3,*lincslam,
                   *bllen0,*ddist;
  static int       nc;

  char             buf[STRLEN];
  int              b,i,j,nit,warn,p_imax,error;
  real             wang,p_max,p_rms;
  real             dt,dt_2,tmp;
  bool             bOK;
  
  bOK = TRUE;
  if (bInit) {
    nc = top->idef.il[F_SHAKE].nr/3;
    init_lincs(stdlog,top,ir,md,start,homenr,
	       &nrtot,
	       &r,&bla1,&bla2,&blnr,&blbnb,
	       &bllen,&blc,&blcc,&blm,&tmp1,&tmp2,&tmp3,&lincslam,
	       &bllen0,&ddist);
  } 
  else if (nc != 0) {
    /* If there are any constraints */
    if (bCoordinates) {
      dt   = ir->delta_t;
      dt_2 = 1.0/(dt*dt);
      
      if (ir->efep != efepNO)
	for(i=0;i<nc;i++)
	  bllen[i]=bllen0[i]+lambda*ddist[i];

      /* Set the zero lengths to the old lengths */
      for(b=0; b<nc; b++)
	if (bllen0[b]==0) {
	  i = bla1[b];
	  j = bla2[b];
	  tmp =   sqr(x[i][XX]-x[j][XX])+
			  sqr(x[i][YY]-x[j][YY])+
			  sqr(x[i][ZZ]-x[j][ZZ]);
	  bllen[b] = tmp * invsqrt(tmp);	  
	}

      wang=ir->LincsWarnAngle;
      
      if (do_per_step(step,ir->nstlog) || step<0)
	cconerr(&p_max,&p_rms,&p_imax,xprime,nc,bla1,bla2,bllen);

      nit = ir->nLincsIter;
      
#ifdef USE_FORTRAN
      F77_FUNC(flincs,FLINCS)(x[0],xprime[0],&nc,bla1,bla2,blnr,blbnb,
			      bllen,blc,blcc,blm,&nit,&ir->nProjOrder,
			      md->invmass,r[0],tmp1,tmp2,tmp3,&wang,&warn,
			      lincslam);
#else
      clincs(x,xprime,nc,bla1,bla2,blnr,blbnb,
	     bllen,blc,blcc,blm,nit,ir->nProjOrder,
	     md->invmass,r,tmp1,tmp2,tmp3,wang,&warn,lincslam);
#endif
      if (ir->efep != efepNO) {
	real dvdl=0;
	
	for(i=0; (i<nc); i++)
	  dvdl+=lincslam[i]*dt_2*ddist[i];
	*dvdlambda+=dvdl;
      }
      
      if (do_per_step(step,ir->nstlog) || (step<0)) {
	fprintf(stdlog,"   Rel. Constraint Deviation:  Max    between atoms     RMS\n");
	fprintf(stdlog,"       Before LINCS         %.6f %6d %6d   %.6f\n",
		p_max,bla1[p_imax]+1,bla2[p_imax]+1,p_rms);
	cconerr(&p_max,&p_rms,&p_imax,xprime,nc,bla1,bla2,bllen);
	fprintf(stdlog,"        After LINCS         %.6f %6d %6d   %.6f\n\n",
		p_max,bla1[p_imax]+1,bla2[p_imax]+1,p_rms);
      }
      
      if (warn > 0) {
	if (bDumpOnError) {
	  cconerr(&p_max,&p_rms,&p_imax,xprime,nc,bla1,bla2,bllen);
	  sprintf(buf,"\nStep %d, time %g (ps)  LINCS WARNING\n"
		  "relative constraint deviation after LINCS:\n"
		  "max %.6f (between atoms %d and %d) rms %.6f\n",
		  step,ir->init_t+step*ir->delta_t,
		  p_max,bla1[p_imax]+1,bla2[p_imax]+1,p_rms);
	  fprintf(stdlog,"%s",buf);
	  fprintf(stderr,"%s",buf);
	  lincs_warning(x,xprime,nc,bla1,bla2,bllen,wang);
	}
	bOK = (p_max < 0.5);
      }
      for(b=0; (b<nc); b++)
	if (bllen0[b] == 0)
	  bllen[b] = 0;
    } 
    else {
#ifdef USE_FORTRAN
      F77_FUNC(flincsp,FLINCSP)(x[0],xprime[0],min_proj[0],&nc,bla1,bla2,blnr,blbnb,
				blc,blcc,blm,&ir->nProjOrder,
				md->invmass,r[0],tmp1,tmp2,tmp3);
#else
      clincsp(x,xprime,min_proj,nc,bla1,bla2,blnr,blbnb,
	      blc,blcc,blm,ir->nProjOrder,
	      md->invmass,r,tmp1,tmp2,tmp3);
#endif
    }

    /* count assuming nit=1 */
    inc_nrnb(nrnb,eNR_LINCS,nc);
    inc_nrnb(nrnb,eNR_LINCSMAT,(2+ir->nProjOrder)*nrtot);
  }
  return bOK;
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
			  real lambda,real *dvdlambda,t_nrnb *nrnb,
			  bool bCoordinates,bool bInit)
{
  static int       nblocks=0;
  static int       *sblock=NULL;
  static int       nsettle,settle_type;
  static int       *owptr;
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
  
  bOK = TRUE;
  if (bInit) {
    /* Output variables, initiate them right away */
    
    if ((ir->etc==etcBERENDSEN) || (ir->epc==epcBERENDSEN))
      please_cite(log,"Berendsen84a");
    
    bDumpOnError = (getenv("NO_SHAKE_ERROR") == NULL);
    
    /* Put the oxygen atoms in the owptr array */
    nsettle=idef->il[F_SETTLE].nr/2;
    if (nsettle > 0) {
      snew(owptr,nsettle);
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
      
      please_cite(log,"Miyamoto92a");
    }
    
    ncons=idef->il[F_SHAKE].nr/3;
    if (ncons > 0) {
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
	please_cite(stdlog,"Hess97a");
	bOK = constrain_lincs(stdlog,top,ir,0,md,start,homenr,&nblocks,&sblock,
			      NULL,NULL,NULL,NULL,0,NULL,bCoordinates,TRUE,nrnb,
			      bDumpOnError);
      } 
      else
	please_cite(stdlog,"Ryckaert77a");
    }
  } 
  else {
    /* !bInit */
    if (nblocks > 0) {
      where();
      
      if (ir->eConstrAlg == estSHAKE)
	bOK = bshakef(stdlog,homenr,md->invmass,nblocks,sblock,idef,
		      ir,box,x,xprime,nrnb,lambda,dvdlambda,bDumpOnError);
      else if (ir->eConstrAlg == estLINCS)
	bOK = constrain_lincs(stdlog,top,ir,step,md,
			      start,homenr,&nblocks,&sblock,
			      x,xprime,min_proj,box,lambda,dvdlambda,
			      bCoordinates,FALSE,nrnb,bDumpOnError);
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
#ifdef USE_FORTRAN
      F77_FUNC(fsettle,FSETTLE)(&nsettle,owptr,x[0],xprime[0],
				&dOH,&dHH,&mO,&mH,&error);
#else
      csettle(stdlog,nsettle,owptr,x[0],xprime[0],dOH,dHH,mO,mH,&error);
#endif
      inc_nrnb(nrnb,eNR_SETTLE,nsettle);
      bOK = (error < 0);
      if (!bOK && bDumpOnError)
	fprintf(stdlog,"\nt = %.3f ps: Water molecule starting at atom %d can not be "
		"settled.\nCheck for bad contacts and/or reduce the timestep.",
		ir->init_t+step*ir->delta_t,owptr[error]+1);
    }
    if (!bOK && bDumpOnError) 
      dump_confs(step,&(top->atoms),x,xprime,box);
  }
  return bOK;
}

bool constrain(FILE *log,t_topology *top,t_inputrec *ir,int step,
	       t_mdatoms *md,int start,int homenr,
	       rvec *x,rvec *xprime,rvec *min_proj,matrix box,
	       real lambda,real *dvdlambda,t_nrnb *nrnb,bool bCoordinates)
{
  return low_constrain(log,top,ir,step,md,start,homenr,x,xprime,min_proj,box,
		       lambda,dvdlambda,nrnb,bCoordinates,FALSE);
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
  low_constrain(log,top,ir,0,md,start,homenr,NULL,NULL,NULL,NULL,
		0,NULL,NULL,bOnlyCoords,TRUE);
  
  return count_constraints(top,cr);
}

void lincs_warning(rvec *x,rvec *xprime,
		   int ncons,int *bla1,int *bla2,real *bllen,real wangle)
{
  int b,i,j;
  rvec v0,v1;
  real wfac,d0,d1,cosine;
  char buf[STRLEN];
  
  wfac=cos(DEG2RAD*wangle);
  
  sprintf(buf,"bonds that rotated more than %g degrees:\n"
	  " atom 1 atom 2  angle  previous, current, constraint length\n",
	  wangle);
  fprintf(stderr,buf);
  fprintf(stdlog,buf); 

  for(b=0;b<ncons;b++) {
    i=bla1[b];
    j=bla2[b];
    rvec_sub(x[i],x[j],v0);
    rvec_sub(xprime[i],xprime[j],v1);
    d0=norm(v0);
    d1=norm(v1);
    cosine=iprod(v0,v1)/(d0*d1);
    if (cosine<wfac) {
      sprintf(buf," %6d %6d  %5.1f  %8.4f %8.4f    %8.4f\n",
	      i+1,j+1,RAD2DEG*acos(cosine),d0,d1,bllen[b]);
      fprintf(stderr,buf);
      fprintf(stdlog,buf);
    }
  }
}

void cconerr(real *max,real *rms,int *imax,rvec *xprime,
	     int ncons,int *bla1,int *bla2,real *bllen)
     
{
  real      len,d,ma,ms,tmp0,tmp1,tmp2,r2;
  int       b,i,j,im;
  
  ma=0;
  ms=0;
  im=0;
  for(b=0;b<ncons;b++) {
    i=bla1[b];
    j=bla2[b];
    tmp0=xprime[i][0]-xprime[j][0];
    tmp1=xprime[i][1]-xprime[j][1];
    tmp2=xprime[i][2]-xprime[j][2];
    r2=tmp0*tmp0+tmp1*tmp1+tmp2*tmp2;
    len=r2*invsqrt(r2);
    d=fabs(len/bllen[b]-1);
    if (d > ma) {
      ma=d;
      im=b;
    }
    ms=ms+d*d;
  }
  *max=ma;
  *rms=sqrt(ms/ncons);
  *imax=im;
}
