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
static char *SRCID_constr_c = "$Id$";

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

static void init_lincs(FILE *log,t_topology *top,t_inputrec *ir,
		       t_mdatoms *md,
		       int *nrtot,
		       rvec **r,int **bla1,int **bla2,int **blnr,int **blbnb,
		       real **bllen,real **blc,real **blcc,real **blm,
		       real **tmp1,real **tmp2,real **tmp3,
		       real **lincslam,real **bllen0,real **ddist)
{
  t_idef      *idef=&(top->idef);
  t_iatom     *iatom;
  int         i,j,k,n,b1,b,cen;
  int         ncons;
  int         type,a1,a2,b2,nr,n1,n2,nc4;
  real        len=0,len1,sign;
  real        im1,im2;
  
  ncons  = idef->il[F_SHAKE].nr/3;
  *nrtot = 0;
  
  if (ncons > 0) {

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
    
    iatom=idef->il[F_SHAKE].iatoms;

    /* Make constraint-neighbor list */
    (*blnr)[0] = 0;
    for(i=0; (i<ncons); i++) {
      j=3*i;
      a1=iatom[j+1];
      a2=iatom[j+2];
      nr=0;
      for(k=0; (k<ncons); k++) {
	b1=iatom[3*k+1];
	b2=iatom[3*k+2];
	if ((a1==b1 || a1==b2) || (a2==b1 || a2==b2)) 
	  if (i != k) nr++;
      }
      *nrtot += nr;
      (*blnr)[i+1] = *nrtot;
      type=iatom[j];
      len =idef->iparams[type].shake.dA;
      len1=idef->iparams[type].shake.dB;
      (*bla1)[i]=a1;
      (*bla2)[i]=a2;
      (*bllen)[i]=len;
      (*bllen0)[i]=len;
      (*ddist)[i]=len1-len;
    }

    fprintf(log,"\nInitializing LINear Constraint Solver\n");
    fprintf(log,"  number of constraint is %d\n",ncons);
    fprintf(log,"  average number of constraints coupled to one constraint is %.1f\n\n",
	    (real)(*nrtot)/ncons);
    fflush(log);

    snew(*blbnb,*nrtot); 
    snew(*blcc,*nrtot);
    snew(*blm,*nrtot); 

    /* Construct the constraint connection matrix blbnb */
    for(i=0; (i<ncons); i++) {
      a1=(*bla1)[i];
      a2=(*bla2)[i];
      im1=md->invmass[a1];
      im2=md->invmass[a2];
      (*blc)[i]=invsqrt(im1+im2);
      nr=0;
      for(k=0; (k<ncons); k++) {
	b1=(*bla1)[k];
	b2=(*bla2)[k];
	if ((a1==b1 || a1==b2) || (a2==b1 || a2==b2)) 
	  if (i != k) {
	    (*blbnb)[(*blnr)[i]+nr]=k;
	    nr++;
	  }
      }
    }
    
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

static void constrain_lincs(FILE *log,t_topology *top,t_inputrec *ir,
			    int step,t_mdatoms *md,int start,int homenr,
			    int *nbl,int **sbl,
			    rvec *x,rvec *xprime,matrix box,
			    real lambda,real *dvdlambda,bool bInit,
			    t_nrnb *nrnb)
{
  static int       *bla1,*bla2,*blnr,*blbnb,nrtot=0;
  static rvec      *r;
  static real      *bllen,*blc,*blcc,*blm,*tmp1,*tmp2,*tmp3,*lincslam,
                   *bllen0,*ddist;
  static int       nc;

  char             buf[STRLEN];
  int              i,nit,warn,p_imax,error;
  real             wang,p_max,p_rms;
  real             dt,dt_2;

  if (bInit) {
    nc = top->idef.il[F_SHAKE].nr/3;
    init_lincs(stdlog,top,ir,md,
	       &nrtot,
	       &r,&bla1,&bla2,&blnr,&blbnb,
	       &bllen,&blc,&blcc,&blm,&tmp1,&tmp2,&tmp3,&lincslam,
	       &bllen0,&ddist);
  } else {
    if (nc == 0)
      return;
    dt   = ir->delta_t;
    dt_2 = 1.0/(dt*dt);
    
    if (ir->efep != efepNO)
      for(i=0;i<nc;i++)
	bllen[i]=bllen0[i]+lambda*ddist[i];
    
    wang=ir->LincsWarnAngle;
    
    if (do_per_step(step,ir->nstLincsout))
      cconerr(&p_max,&p_rms,&p_imax,xprime,nc,bla1,bla2,bllen);

    if ((ir->eI == eiSteep) || (ir->eI == eiCG))
      /* Use more iterations when doing energy minimization, *
       * because we need very accurate positions and forces. */
      nit = ir->nProjOrder;
    else
      nit = 1;
    
#ifdef USEF77
    flincs(x[0],xprime[0],&nc,bla1,bla2,blnr,blbnb,
	   bllen,blc,blcc,blm,&nit,&ir->nProjOrder,
	   md->invmass,r[0],tmp1,tmp2,tmp3,&wang,&warn,
	   lincslam);
#else
    clincs(x,xprime,nc,bla1,bla2,blnr,blbnb,
	   bllen,blc,blcc,blm,nit,ir->nProjOrder,
	   md->invmass,r,tmp1,tmp2,tmp3,wang,&warn,lincslam);
#endif
    /* count assuming nit=1 */
    inc_nrnb(nrnb,eNR_LINCS,nc);
    inc_nrnb(nrnb,eNR_LINCSMAT,(2+ir->nProjOrder)*nrtot);
    if (ir->efep != efepNO) {
      real dvdl=0;
      
      for(i=0; (i<nc); i++)
	dvdl+=lincslam[i]*dt_2*ddist[i];
      *dvdlambda+=dvdl;
    }
    
    if (do_per_step(step,ir->nstLincsout)) {
      fprintf(stdlog,"Step %d\nRel. Constraint Deviation:  Max    between atoms     RMS\n",step);
      fprintf(stdlog,"    Before LINCS         %.6f %6d %6d   %.6f\n",
	      p_max,bla1[p_imax]+1,bla2[p_imax]+1,p_rms);
      cconerr(&p_max,&p_rms,&p_imax,xprime,nc,bla1,bla2,bllen);
      fprintf(stdlog,"     After LINCS         %.6f %6d %6d   %.6f\n\n",
	      p_max,bla1[p_imax]+1,bla2[p_imax]+1,p_rms);
    }
    
    if (warn > 0) {
      cconerr(&p_max,&p_rms,&p_imax,xprime,nc,bla1,bla2,bllen);
      sprintf(buf,"\nStep %d, time %g (ps)  LINCS WARNING\n"
	      "relative constraint deviation after LINCS:\n"
	      "max %.6f (between atoms %d and %d) rms %.6f\n",
	      step,ir->init_t+step*ir->delta_t,
	      p_max,bla1[p_imax]+1,bla2[p_imax]+1,p_rms);
      fprintf(stdlog,"%s",buf);
      fprintf(stderr,"%s",buf);
      lincs_warning(x,xprime,nc,bla1,bla2,bllen,wang);
      if (p_max > 0.5) {
	dump_confs(step,&(top->atoms),x,xprime,box);
	fatal_error(0,"Bond deviates more than half its own length");
      }
    }
  }
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
			  rvec *x,rvec *xprime,matrix box,
			  real lambda,real *dvdlambda,t_nrnb *nrnb,
			  bool bInit)
{
  static int       nblocks=0;
  static int       *sblock=NULL;
  static int       nsettle,settle_type;
  static int       *owptr;

  char        buf[STRLEN];
  t_sortblock *sb;
  t_block     *blocks=&(top->blocks[ebSBLOCKS]);
  t_idef      *idef=&(top->idef);
  t_iatom     *iatom;
  atom_id     *inv_sblock;
  int         i,j,m,bnr;
  int         ncons,bstart,error;
  
  if (bInit) {
    /* Output variables, initiate them right away */
    
    if ((ir->btc) || (ir->epc != epcNO))
      please_cite(log,"Berendsen84a");
    
    /* Put the oxygen atoms in the owptr array */
    nsettle=idef->il[F_SETTLE].nr/2;
    if (nsettle > 0) {
      snew(owptr,nsettle);
      settle_type=idef->il[F_SETTLE].iatoms[0];
      for (j=0; (j<idef->il[F_SETTLE].nr); j+=2) {
	if (idef->il[F_SETTLE].iatoms[j] != settle_type)
	  fatal_error(0,"More than one settle type (%d and %d)",
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
      bstart=(idef->pid > 0) ? blocks->multinr[idef->pid-1] : 0;
      nblocks=blocks->multinr[idef->pid] - bstart;
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
	fatal_error(0,"DEATH HORROR: "
		    "top->blocks[ebSBLOCKS] does not match idef->il[F_SHAKE]");
      }
      sfree(sb);
      sfree(inv_sblock);
    }
    
    if (idef->il[F_SHAKE].nr) {
      if (ir->eConstrAlg == estLINCS) {
	please_cite(stdlog,"Hess97a");
	constrain_lincs(stdlog,top,ir,0,md,
			start,homenr,
			&nblocks,&sblock,
			NULL,NULL,NULL,0,NULL,TRUE,nrnb);
      } else
	please_cite(stdlog,"Ryckaert77a");
    }
  } else {
    if (nblocks > 0) {
      where();
      
      if (ir->eConstrAlg == estSHAKE) {
	error=bshakef(stdlog,homenr,md->invmass,nblocks,sblock,idef,
		      ir,box,x,xprime,nrnb,dvdlambda);
	if (error == -1) {
	  dump_confs(step,&(top->atoms),x,xprime,box);
	  fprintf(stdlog,"SHAKE ERROR at step %d\n",step);
	  fatal_error(0,"SHAKE ERROR at step %d\n",step);
	}
      }
      
      if (ir->eConstrAlg == estLINCS)
	constrain_lincs(stdlog,top,ir,step,md,
			start,homenr,
			&nblocks,&sblock,
			x,xprime,box,lambda,dvdlambda,FALSE,nrnb);
    }
    if (nsettle > 0) {
      int  ow1;
      real mO,mH,dOH,dHH;
      
      ow1  = owptr[0];
      mO   = md->massA[ow1];
      mH   = md->massA[ow1+1];
      dOH  = top->idef.iparams[settle_type].settle.doh;
      dHH  = top->idef.iparams[settle_type].settle.dhh;
#ifdef USEF77
      fsettle(&nsettle,owptr,x[0],xprime[0],&dOH,&dHH,&mO,&mH,&error);
#else
      csettle(stdlog,nsettle,owptr,x[0],xprime[0],dOH,dHH,mO,mH,&error);
#endif
      inc_nrnb(nrnb,eNR_SETTLE,nsettle);
      if (error>=0) {
	dump_confs(step,&(top->atoms),x,xprime,box);
	sprintf(buf,"\nMolecule starting at atomnr. %d can not be settled, "
		"step %d, time %g (ps)",
		owptr[error]+1,step,ir->init_t+step*ir->delta_t);
	fprintf(stdlog,buf);
	fatal_error(0,buf);
      }
    }
  }
  
  return (nsettle+idef->il[F_SHAKE].nr > 0);
}

void constrain(FILE *log,t_topology *top,t_inputrec *ir,int step,
	       t_mdatoms *md,int start,int homenr,
	       rvec *x,rvec *xprime,matrix box,
	       real lambda,real *dvdlambda,t_nrnb *nrnb)
{
  low_constrain(log,top,ir,step,md,start,homenr,x,xprime,box,
		lambda,dvdlambda,nrnb,FALSE);
}

bool init_constraints(FILE *log,t_topology *top,t_inputrec *ir,
		      t_mdatoms *md,int start,int homenr)
{
  return low_constrain(log,top,ir,0,md,start,homenr,NULL,NULL,NULL,
		       0,NULL,NULL,TRUE);
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
  real      len,d,ma,ms,tmp0,tmp1,tmp2;
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
    len=sqrt(tmp0*tmp0+tmp1*tmp1+tmp2*tmp2);
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
