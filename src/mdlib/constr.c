#include "typedefs.h"
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

typedef struct {
  atom_id iatom[3];
  atom_id blocknr;
} t_sortblock;

int pcomp(const void *p1, const void *p2)
{
  int     db;
  atom_id min1,min2,max1,max2;
  t_sortblock *a1=(t_sortblock *)p1;
  t_sortblock *a2=(t_sortblock *)p2;
#ifdef DEBUG
  pcount++;
#endif
  
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
		       t_mdatoms *md,int start,int homenr,
		       int *nbl,int **sbl,
		       int *ncm,int *cmax,
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
  real        im1,im2,imcen;
  
  ncons=idef->il[F_SHAKE].nr/3;

  if (ncons > 0) {

    /* Make constraint-neighbour list */

    snew(*r,ncons);
    snew(*bla1,ncons);
    snew(*bla2,ncons);
    snew(*blnr,ncons);
    snew(*bllen,ncons);
    snew(*blc,ncons);
    snew(*tmp1,ncons);
    snew(*tmp2,ncons);
    snew(*tmp3,ncons);
    snew(*lincslam,ncons);
    snew(*bllen0,ncons);
    snew(*ddist,ncons);
    
    iatom=idef->il[F_SHAKE].iatoms;

    /* Number of coupling of a bond is defined as the number of
       bonds directly connected to that bond (not to an atom!).
       The constraint are divided into two groups, the first
       group consists of bonds with 4 or less couplings, the
       second group consists of bonds with more than 4 couplings
       (in proteins most of the bonds have 2 or 4 couplings). 
       
       cmax: maximum number of bonds coupled to one bond 
       ncm:  number of bonds with mor than 4 couplings 
       */

    *cmax=0;
    n1=0;
    n2=ncons-1;

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
      if (nr > *cmax) *cmax=nr;
      type=iatom[j];
      len =idef->iparams[type].shake.dA;
      len1=idef->iparams[type].shake.dB;
      if (nr <=4) {
	(*bla1)[n1]=a1;
	(*bla2)[n1]=a2;
	(*bllen)[n1]=len;
	(*bllen0)[n1]=len;
	(*ddist)[n1]=len1-len;
	n1++;
      }
      else {
	(*bla1)[n2]=a1;
	(*bla2)[n2]=a2;
	(*bllen)[n2]=len;
	(*bllen0)[n2]=len;
	(*ddist)[n2]=len1-len;
	n2--;
      }
    }
    
    *ncm=ncons-n1;
    nc4=(*cmax-4)*n1;

    i=4*n1+(*cmax)*(*ncm);
    snew(*blbnb,i); 
    snew(*blcc,i);
    snew(*blm,i); 

    for(i=0; (i<ncons); i++) {
      a1=(*bla1)[i];
      a2=(*bla2)[i];
      im1=md->invmass[a1];
      im2=md->invmass[a2];
      (*blc)[i]=invsqrt(im1+im2);
      /* printf("%d %d %f\n",a1,a2,blist[i].len); */
      nr=0;
      /* printf("\n%d",i); */
      for(k=0; (k<ncons); k++) {
	b1=(*bla1)[k];
	b2=(*bla2)[k];
	if ((a1==b1 || a1==b2) || (a2==b1 || a2==b2)) 
	  if (i != k) {
	    if (i<n1)
	      (*blbnb)[4*i+nr]=k;
	    else
	      (*blbnb)[(*cmax)*i-nc4+nr]=k;
	    nr++;
	    /* printf(" %d",k); */
	  }
      }
      (*blnr)[i]=nr;
      /* printf("\n"); */
    }
    fprintf(stdlog,"\nInitializing LINear Constraint Solver\n");
    fprintf(stdlog,"%d constraints\nof which %d with more than 4 neighbours\n",ncons,*ncm);
    fprintf(stdlog,"maximum number of bonds coupled to one bond is %d\n\n",*cmax);
    fflush(stdlog);
 
    for(b=0; (b<ncons); b++) {
      i=(*bla1)[b];
      j=(*bla2)[b];
      nr=(*blnr)[b];
      if (nr) 
	for(n=0; (n<nr);n++) {
	  if (b < n1) 
	    k=(*blbnb)[4*b+n];
	  else
	    k=(*blbnb)[(*cmax)*b-nc4+n];
	  if (i==(*bla1)[k] || j==(*bla2)[k])
	    sign=-1;
	    else
	    sign=1;
	  if (i==(*bla1)[k] || i==(*bla2)[k])
	    cen=i;
	  else
	    cen=j;
	  if (ir->eI==eiMD) {
	    imcen=md->invmass[cen];
	    len=sign*imcen*(*blc)[b]*(*blc)[k];
	  }
	  if (ir->eI==eiLD) 
	    len=sign*0.5;
	  if (b<n1) 
	    (*blcc)[4*b+n]=len;
	  else
	    (*blcc)[(*cmax)*b-nc4+n]=len;
	}
    }
    
  }
}

static void constrain_lincs(FILE *log,t_topology *top,t_inputrec *ir,
			    int step,t_mdatoms *md,int start,int homenr,
			    int *nbl,int **sbl,
			    rvec *x,rvec *xprime,matrix box,
			    real lambda,real *dvdlambda,bool bInit)
{
  static int       *bla1,*bla2,*blnr,*blbnb;
  static rvec      *r;
  static real      *bllen,*blc,*blcc,*blm,*tmp1,*tmp2,*tmp3,*lincslam,
                   *bllen0,*ddist;
  static int       nc,ncm,cmax;

  char             buf[STRLEN];
  int              i;
  int              warn,p_imax,error;
  real             wang,p_max,p_rms;
  real             dt,dt_2;

  if (bInit) {
    nc = top->idef.il[F_SHAKE].nr/3;
    init_lincs(stdlog,top,ir,md,
	       start,homenr,
	       nbl,sbl,
	       &ncm,&cmax,
	       &r,&bla1,&bla2,&blnr,&blbnb,
	       &bllen,&blc,&blcc,&blm,&tmp1,&tmp2,&tmp3,&lincslam,
	       &bllen0,&ddist);
  } else {
    dt   = ir->delta_t;
    dt_2 = 1.0/(dt*dt);
    
    if (ir->bPert)
      for(i=0;i<nc;i++)
	bllen[i]=bllen0[i]+lambda*ddist[i];
    
    wang=ir->LincsWarnAngle;
    
    if (do_per_step(step,ir->nstLincsout))
      cconerr(&p_max,&p_rms,&p_imax,xprime,nc,bla1,bla2,bllen);
    
    if (ir->eI != eiLD) {
#ifdef USEF77
      flincs(x[0],xprime[0],&nc,&ncm,&cmax,bla1,bla2,blnr,blbnb,
	     bllen,blc,blcc,blm,&ir->nProjOrder,
	     md->invmass,r[0],tmp1,tmp2,tmp3,&wang,&warn,
	     lincslam);
#else
      clincs(x,xprime,nc,ncm,cmax,bla1,bla2,blnr,blbnb,
	     bllen,blc,blcc,blm,ir->nProjOrder,
	     md->invmass,r,tmp1,tmp2,tmp3,wang,&warn,lincslam);
#endif
      if (ir->bPert) {
	real dvdl=0;
	
	for(i=0; (i<nc); i++)
	  dvdl+=lincslam[i]*dt_2*ddist[i];
	*dvdlambda+=dvdl;
      }
    }
    
    if (ir->eI==eiLD) {
#ifdef USEF77
      flincsld(x[0],xprime[0],&nc,&ncm,&cmax,bla1,bla2,blnr,
	       blbnb,bllen,blcc,blm,&ir->nProjOrder,
	       r[0],tmp1,tmp2,tmp3,&wang,&warn);
#else
      clincsld(x,xprime,nc,ncm,cmax,bla1,bla2,blnr,
	       blbnb,bllen,blcc,blm,ir->nProjOrder,
	       r,tmp1,tmp2,tmp3,wang,&warn);
#endif
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
#ifdef DEBUGIDEF
      fprintf(stdlog,"ncons: %d, bstart: %d, nblocks: %d\n",
	      ncons,bstart,nblocks);
      fflush(stdlog);
#endif
      
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
#ifdef DEBUG
      pr_sortblock(log,"Before sorting",ncons,sb);
      fprintf(log,"Going to sort constraints\n");
#endif
      
      qsort(sb,ncons,(size_t)sizeof(*sb),pcomp);
      
#ifdef DEBUG
      fprintf(log,"I used %d calls to pcomp\n",pcount);
      pr_sortblock(log,"After sorting",ncons,sb);
#endif
      
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
	fprintf(stdlog,"bstart: %d\n",bstart);
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
    
    if (idef->il[F_SHAKE].nr)
      if (ir->eConstrAlg == estLINCS) {
	please_cite(stdlog,"Hess97a");
	constrain_lincs(stdlog,top,ir,0,md,
			start,homenr,
			&nblocks,&sblock,
			NULL,NULL,NULL,0,NULL,TRUE);
      } else
	please_cite(stdlog,"Ryckaert77a");
  } else {
    if (nblocks > 0) {
      where();
      
      if (ir->eConstrAlg == estSHAKE) {
	ncons=bshakef(stdlog,homenr,md->invmass,nblocks,sblock,idef,
		      ir,box,x,xprime,nrnb);
	if (ncons == -1) {
	  dump_confs(step,&(top->atoms),x,xprime,box);
	  fprintf(stdlog,"SHAKE ERROR at step %d\n",step);
	  fatal_error(0,"SHAKE ERROR at step %d\n",step);
	}
      }
      
      if (ir->eConstrAlg == estLINCS)
	constrain_lincs(stdlog,top,ir,step,md,
			start,homenr,
			&nblocks,&sblock,
			x,xprime,box,lambda,dvdlambda,FALSE);
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
