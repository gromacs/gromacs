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
 * GRoups of Organic Molecules in ACtion for Science
 */
static char *SRCID_do_gct_c = "$Id$";

#include "typedefs.h"
#include "do_gct.h"
#include "block_tx.h"
#include "futil.h"
#include "xvgr.h"
#include "macros.h"
#include "physics.h"
#include "network.h"
#include "smalloc.h"
#include "string2.h"
#include "readinp.h"
#include "filenm.h"
#include "mdrun.h"
#include "update.h"

t_coupl_rec *init_coupling(FILE *log,int nfile,t_filenm fnm[],
			   t_commrec *cr,t_forcerec *fr,
			   t_mdatoms *md,t_idef *idef)
{
  int         i,nc,index,j;
  int         ati,atj;
  t_coupl_rec *tcr;
  
  snew(tcr,1);
  read_gct (opt2fn("-j",nfile,fnm), tcr);
  write_gct(opt2fn("-jo",nfile,fnm),tcr,idef);
  
  if ((tcr->dipole != 0.0) && (!ftp2bSet(efNDX,nfile,fnm))) {
    fatal_error(0,"Trying to use polarization correction to energy without specifying an index file for molecule numbers. This will generate huge induced dipoles, and hence not work!");
  }
  
  copy_ff(tcr,fr,md,idef);
    
  /* Update all processors with coupling info */
  if (PAR(cr))
    comm_tcr(log,cr,&tcr);
  
  return tcr;
}

real Ecouple(t_coupl_rec *tcr,real ener[])
{
  if (tcr->bInter)
    return ener[F_SR]+ener[F_LJ]+ener[F_LR]+ener[F_LJLR];
  else
    return ener[F_EPOT];
}

char *mk_gct_nm(char *fn,int ftp,int ati,int atj)
{
  static char buf[256];
  
  strcpy(buf,fn);
  if (atj == -1)
    sprintf(buf+strlen(fn)-4,"%d.%s",ati,ftp2ext(ftp));
  else
    sprintf(buf+strlen(fn)-4,"%d_%d.%s",ati,atj,ftp2ext(ftp));
  
  return buf;
}

static void pr_ff(t_coupl_rec *tcr,real time,t_idef *idef,real ener[],t_commrec *cr,
		  int nfile,t_filenm fnm[])
{
  static FILE *prop;
  static FILE **out=NULL;
  static FILE **qq=NULL;
  static FILE **ip=NULL;
  t_coupl_LJ  *tclj;
  t_coupl_BU  *tcbu;
  char        buf[256];
  char        *leg[] = { "C12", "C6" };
  char        *bleg[] = { "A", "B", "C" };
  char        *pleg[] = { "RA-Pres", "Pres", "RA-Epot", "Epot" };
  int         i,index;
  
  if ((prop == NULL) && (out == NULL) && (qq == NULL) && (ip == NULL)) {
    prop=xvgropen(opt2fn("-runav",nfile,fnm),
		  "Properties and Running Averages","Time (ps)","");
    xvgr_legend(prop,asize(pleg),pleg);
    if (tcr->nLJ) {
      snew(out,tcr->nLJ);
      for(i=0; (i<tcr->nLJ); i++) {
	if (tcr->tcLJ[i].bPrint) {
	  tclj   = &(tcr->tcLJ[i]);
	  out[i] = 
	    xvgropen(mk_gct_nm(opt2fn("-ffout",nfile,fnm),
			       efXVG,tclj->at_i,tclj->at_j),
		     "General Coupling Lennard Jones","Time (ps)",
		     "Force constant (units)");
	  fprintf(out[i],"@ subtitle \"Interaction between types %d and %d\"\n",
		  tclj->at_i,tclj->at_j);
	  xvgr_legend(out[i],asize(leg),leg);
	  fflush(out[i]);
	}
      }
    }
    else if (tcr->nBU) {
      snew(out,tcr->nBU);
      for(i=0; (i<tcr->nBU); i++) {
	if (tcr->tcBU[i].bPrint) {
	  tcbu=&(tcr->tcBU[i]);
	  out[i] = 
	    xvgropen(mk_gct_nm(opt2fn("-ffout",nfile,fnm),efXVG,
			       tcbu->at_i,tcbu->at_j),
		     "General Coupling Buckingham","Time (ps)",
		     "Force constant (units)");
	  fprintf(out[i],"@ subtitle \"Interaction between types %d and %d\"\n",
		  tcbu->at_i,tcbu->at_j);
	  xvgr_legend(out[i],asize(bleg),bleg);
	  fflush(out[i]);
	}
      }
    }
    snew(qq,tcr->nQ);
    for(i=0; (i<tcr->nQ); i++) {
      if (tcr->tcQ[i].bPrint) {
	qq[i] = xvgropen(mk_gct_nm(opt2fn("-ffout",nfile,fnm),efXVG,
				   tcr->tcQ[i].at_i,-1),
			 "General Coupling Charge","Time (ps)","Charge (e)");
	fprintf(qq[i],"@ subtitle \"Type %d\"\n",tcr->tcQ[i].at_i);
	fflush(qq[i]);
      }
    }
    snew(ip,tcr->nIP);
    for(i=0; (i<tcr->nIP); i++) {
      sprintf(buf,"gctIP%d",tcr->tIP[i].type);
      ip[i]=xvgropen(mk_gct_nm(opt2fn("-ffout",nfile,fnm),efXVG,0,-1),
		     "General Coupling iparams","Time (ps)","ip ()");
      index=tcr->tIP[i].type;
      fprintf(ip[i],"@ subtitle \"Coupling to %s\"\n",
	      interaction_function[idef->functype[index]].longname);
      fflush(ip[i]);
    }
  }
  fprintf(prop,"%10g  %12.5e  %12.5e  %12.5e  %12.5e\n",
	  time,tcr->pres,ener[F_PRES],tcr->epot,Ecouple(tcr,ener));
  fflush(prop);
  for(i=0; (i<tcr->nLJ); i++) {
    tclj=&(tcr->tcLJ[i]);
    if (tclj->bPrint) {
      fprintf(out[i],"%14.7e  %14.7e  %14.7e\n",
	      time,tclj->c12,tclj->c6);
      fflush(out[i]);
    }
  }
  for(i=0; (i<tcr->nBU); i++) {
    tcbu=&(tcr->tcBU[i]);
    if (tcbu->bPrint) {
      fprintf(out[i],"%14.7e  %14.7e  %14.7e  %14.7e\n",
	      time,tcbu->a,tcbu->b,tcbu->c);
      fflush(out[i]);
    }
  }
  for(i=0; (i<tcr->nQ); i++) {
    if (tcr->tcQ[i].bPrint) {
      fprintf(qq[i],"%14.7e  %14.7e\n",time,tcr->tcQ[i].Q);
      fflush(qq[i]);
    }
  }
  for(i=0; (i<tcr->nIP); i++) {
    fprintf(ip[i],"%10g  ",time);
    index=tcr->tIP[i].type;
    switch(idef->functype[index]) {
    case F_BONDS:
      fprintf(ip[i],"%10g  %10g\n",tcr->tIP[i].iprint.harmonic.krA,
	      tcr->tIP[i].iprint.harmonic.rA);
      break;
    default:
      break;
    }
    fflush(ip[i]);
  }
}

static void pr_dev(real t,real dev[eoNR],t_commrec *cr,int nfile,t_filenm fnm[])
{
  static FILE *fp=NULL;
  int    i;
  
  if (!fp) {
    fp=xvgropen(opt2fn("-devout",nfile,fnm),
		"Deviations from target value","Pres","Epot");
  }
  fprintf(fp,"%12.5e  %12.5e\n",dev[eoPres],dev[eoEpot]);
  fflush(fp);
}

static void upd_nbfplj(FILE *log,t_coupl_LJ *tclj,
		       real **nbfp,int ati,int atj,int atnr,real f6,real f12)
{
  int n;
  
  /* update the nonbonded force parameters */
  if (atj != -1) {
    C6 (nbfp,ati,atj) *= f6;
    C12(nbfp,ati,atj) *= f12;
    C6 (nbfp,atj,ati) *= f6;
    C12(nbfp,atj,ati) *= f12;
    
    /* Save for printing */
    tclj->c6  = C6(nbfp,ati,atj);
    tclj->c12 = C12(nbfp,ati,atj);
  }
  else {
    /* To implement *GEOMETRIC* mixing rules, the factors on the force
     * constants must be rooted first. This way the parameters
     * on the diagonal will be multiplied by the full constant.
     */
    f6  = sqrt(f6);
    f12 = sqrt(f12);
    
    if (debug)
      fprintf(debug,"Updating LJ-FF for %d entries, f6=%g, f12=%g\n",atnr,f6,f12);
    
    for(n=0; (n<atnr); n++) {
      C6 (nbfp,ati,n) *= f6;
      C12(nbfp,ati,n) *= f12;
      C6 (nbfp,n,ati) *= f6;
      C12(nbfp,n,ati) *= f12;
    }
    /* Save diagonal elements for printing */
    tclj->c6  = C6 (nbfp,ati,ati);
    tclj->c12 = C12(nbfp,ati,ati);
  }
}

static void upd_nbfplj2(FILE *log,real **nbfp,int atnr,real f6[],real f12[])
{
  int n,m,k;
  
  /* Update the nonbonded force parameters */
  for(k=n=0; (n<atnr); n++) {
    for(m=0; (m<atnr); m++,k++) {
      C6 (nbfp,n,m) *= f6[k];
      C12(nbfp,n,m) *= f12[k];
    }
  }
}

static void upd_nbfpbu2(FILE *log,real **nbfp,int atnr,real fa[],real fb[],real fc[])
{
  int n,m,k;
  
  /* Update the nonbonded force parameters */
  for(k=n=0; (n<atnr); n++) {
    for(m=0; (m<atnr); m++) {
      (nbfp)[n][3*m]   *= fa[k];
      (nbfp)[n][3*m+1] *= fb[k];
      (nbfp)[n][3*m+2] *= fc[k];
    }
  }
}

static void upd_nbfpbu(FILE *log,t_coupl_BU *tcbu,
		       real **nbfp,int ati,int atj,int atnr,
		       real fa,real fb,real fc)
{
  int n;
  
  /* update the nonbonded force parameters */
  if (atj != -1) {
    (nbfp)[(ati)][3*(atj)]   *= fa;
    (nbfp)[(ati)][3*(atj)+1] *= fb;
    (nbfp)[(ati)][3*(atj)+2] *= fc;
    (nbfp)[(atj)][3*(ati)]   *= fa;
    (nbfp)[(atj)][3*(ati)+1] *= fb;
    (nbfp)[(atj)][3*(ati)+2] *= fc;
    
    /* Save for printing */
    tcbu->a = (nbfp)[(ati)][3*(atj)];
    tcbu->b = (nbfp)[(ati)][3*(atj)+1];
    tcbu->c = (nbfp)[(ati)][3*(atj)+2];
  }
  else {
    /* To implement *GEOMETRIC* mixing rules, the factors on the force
     * constants must be rooted first. This way the parameters
     * on the diagonal will be multiplied by the full constant.
     */
    fa  = sqrt(fa);
    fb  = sqrt(fb);
    fc  = sqrt(fc);
    if (debug)
      fprintf(debug,"Updating Buck-FF for %d entries, fa=%g, fb=%g, fc=%g\n",
	      atnr,fa,fb,fc);
    
    for(n=0; (n<atnr); n++) {
      (nbfp)[(ati)][3*(n)]   *= fa;
      (nbfp)[(ati)][3*(n)+1] *= fb;
      (nbfp)[(ati)][3*(n)+2] *= fc;
      (nbfp)[(n)][3*(ati)]   *= fa;
      (nbfp)[(n)][3*(ati)+1] *= fb;
      (nbfp)[(n)][3*(ati)+2] *= fc;
    }
    /* Save diagonal elements for printing */
    tcbu->a = (nbfp)[(ati)][3*(ati)];
    tcbu->b = (nbfp)[(ati)][3*(ati)+1];
    tcbu->c = (nbfp)[(ati)][3*(ati)+2];
  }
}

void gprod(t_commrec *cr,int n,real f[])
{
  /* Compute the global product of all elements in an array 
   * such that after gprod f[i] = PROD_j=1,nprocs f[i][j]
   */
  int  i,j;
  real *buf;
  
  snew(buf,n);
  
  for(i=0; (i<=cr->nprocs); i++) {
    gmx_tx(cr->left,array(f,n));
    gmx_rx(cr->right,array(buf,n));
    gmx_wait(cr->left,cr->right);
    for(j=0; (j<n); j++)
      f[j] *= buf[j];
  }
  sfree(buf);
}

void set_factor_matrix(int ntypes,real f[],real fmult,int ati,int atj)
{
#define FMIN 0.95
#define FMAX 1.05
  int i;

  fmult = min(FMAX,max(FMIN,fmult));  
  if (atj != -1) {
    f[ntypes*ati+atj] *= fmult;
    f[ntypes*atj+ati] *= fmult;
  }
  else {
    for(i=0; (i<ntypes); i++) {
      f[ntypes*ati+i] *= fmult;
      f[ntypes*i+ati] *= fmult;
    }
  }
#undef FMIN
#undef FMAX
}

void do_coupling(FILE *log,int nfile,t_filenm fnm[],
		 t_coupl_rec *tcr,real t,int step,real ener[],
		 t_forcerec *fr,t_inputrec *ir,bool bMaster,
		 t_mdatoms *md,t_idef *idef,real mu_aver,int nmols,
		 t_commrec *cr)
{

#define enm2Debye 48.0321
#define d2e(x) (x)/enm2Debye
#define enm2kjmol(x) (x)*0.0143952 /* = 2.0*4.0*M_PI*EPSILON0 */

  static real *f6,*f12,*fa,*fb,*fc,*fq;
  static bool bFirst = TRUE;
  
  int         i,j,ati,atj,atnr2,type,ftype;
  real        deviation[eoNR];
  real        ff6,ff12,ffa,ffb,ffc,ffq,factor,dt,mu_ind,Epol,Eintern;
  bool        bTest,bPrint;
  t_coupl_LJ  *tclj;
  t_coupl_BU  *tcbu;
  t_coupl_Q   *tcq;
  t_coupl_iparams *tip;
  
  atnr2 = idef->atnr * idef->atnr;
  if (bFirst) {
    if (PAR(cr))
      fprintf(log,"DOGCT: this is parallel\n");
    else
      fprintf(log,"DOGCT: this is not parallel\n");
    snew(f6, atnr2);
    snew(f12,atnr2);
    snew(fa, atnr2);
    snew(fb, atnr2);
    snew(fc, atnr2);
    snew(fq, idef->atnr);
    bFirst = FALSE;
  }
  
  bPrint = do_per_step(step,ir->nstprint);
  dt     = ir->delta_t;

  /* Initiate coupling to the reference pressure and temperature to start
   * coupling slowly.
   */
  if (step == 0) {
    tcr->pres = tcr->pres0;
    tcr->epot = tcr->epot0*nmols;
  }

  /* We want to optimize the LJ params, usually to the Vaporization energy 
   * therefore we only count intermolecular degrees of freedom.
   * Note that this is now optional. switch UseEinter to yes in your gct file
   * if you want this.
   */
  Eintern = Ecouple(tcr,ener);
  
  /* Use a memory of tcr->nmemory steps, so we actually couple to the
   * average observable over the last tcr->nmemory steps. This may help
   * in avoiding local minima in parameter space.
   */
  tcr->pres = run_aver(tcr->pres,ener[F_PRES],step,tcr->nmemory);
  tcr->epot = run_aver(tcr->epot,Eintern,     step,tcr->nmemory);
  
  if (bPrint)
    pr_ff(tcr,t,idef,ener,cr,nfile,fnm);
  
  /* Calculate the deviation of average value from the target value */
  deviation[eoPres] = tcr->pres0 - tcr->pres;

  /* if dipole != 0.0 assume we want to use polarization corrected coupling */
  if ((tcr->dipole) != 0.0) {
    mu_ind = mu_aver - d2e(tcr->dipole); /* in e nm */
    /* Epol = mu_ind*mu_ind/(2.0*(tcr->polarizability)*4.0*M_PI*EPSILON0); */
    Epol = mu_ind*mu_ind/(enm2kjmol(tcr->polarizability));

    deviation[eoEpot] = (tcr->epot0 - Epol)*nmols - tcr->epot;
    if (debug) {
      fprintf(debug,"mu_ind: %g (%g D) mu_aver: %g (%g D)\n",
	      mu_ind,mu_ind*enm2Debye,mu_aver,mu_aver*enm2Debye);
      fprintf(debug,"Eref %g Epol %g Erunav %g Dev %g Eact %g\n",
	      (tcr->epot0)*nmols, Epol*nmols, tcr->epot,
	      deviation[eoEpot],ener[F_EPOT]);
    }
  }
  else 
    deviation[eoEpot] = tcr->epot0*nmols - tcr->epot;
  
  if (bPrint)
    pr_dev(t,deviation,cr,nfile,fnm);
  
  /* First set all factors to 1 */
  for(i=0; (i<atnr2); i++) {
    f6[i] = f12[i] = fa[i] = fb[i] = fc[i] = 1.0;
  }
  for(i=0; (i<idef->atnr); i++) 
    fq[i] = 1.0;
  
  
  if (!fr->bBHAM) {
    for(i=0; (i<tcr->nLJ); i++) {
      tclj=&(tcr->tcLJ[i]);
      
      factor=deviation[tclj->eObs];
      
      ati=tclj->at_i;
      atj=tclj->at_j;
      
      ff6 = ff12 = 1.0;	
      if (tclj->xi_6)      
	ff6  += (dt/tclj->xi_6)  * factor;
      if (tclj->xi_12)     
	ff12 += (dt/tclj->xi_12) * factor;
      
      set_factor_matrix(idef->atnr,f6, sqrt(ff6), ati,atj);
      set_factor_matrix(idef->atnr,f12,sqrt(ff12),ati,atj);
      
      /*f6  = min(max(f6, FMIN),FMAX);
	f12 = min(max(f12,FMIN),FMAX);
	
	  upd_nbfplj(log,tclj,fr->nbfp,ati,atj,idef->atnr,f6,f12);
      */
    }
    if (PAR(cr)) {
      gprod(cr,atnr2,f6);
      gprod(cr,atnr2,f12);
    }
    upd_nbfplj2(log,fr->nbfp,idef->atnr,f6,f12);
    
    /* Copy for printing */
    for(i=0; (i<tcr->nLJ); i++) {
      tclj=&(tcr->tcLJ[i]);
      tclj->c6  =  C6(fr->nbfp,tclj->at_i,tclj->at_i);
      tclj->c12 = C12(fr->nbfp,tclj->at_i,tclj->at_i);
    }
  }
  else {
    for(i=0; (i<tcr->nBU); i++) {
      tcbu=&(tcr->tcBU[i]);
      
      factor=deviation[tcbu->eObs];
      
      if (tcbu->xi_a)      
	ffa = 1 + (dt/tcbu->xi_a)  * factor;
      if (tcbu->xi_b)      
	  ffb = 1 + (dt/tcbu->xi_b)  * factor;
      if (tcbu->xi_c)      
	ffc = 1 + (dt/tcbu->xi_c)  * factor;
      
      ati=tcbu->at_i;
      atj=tcbu->at_j;
      set_factor_matrix(idef->atnr,fa,sqrt(ffa),ati,atj);
      set_factor_matrix(idef->atnr,fa,sqrt(ffb),ati,atj);
      set_factor_matrix(idef->atnr,fc,sqrt(ffc),ati,atj);
      /*fa  = min(max(fa, FMIN),FMAX);
	fb  = min(max(fb, FMIN),FMAX);
	fc  = min(max(fc, FMIN),FMAX);
	upd_nbfpbu(log,tcbu,fr->nbfp,ati,atj,idef->atnr,fa,fb,fc);
      */
    }
    if (PAR(cr)) {
      gprod(cr,atnr2,fa);
      gprod(cr,atnr2,fb);
      gprod(cr,atnr2,fc);
    }
    upd_nbfpbu2(log,fr->nbfp,idef->atnr,fa,fb,fc);
  }
  
  for(i=0; (i<tcr->nQ); i++) {
    tcq=&(tcr->tcQ[i]);
    
    if (tcq->xi_Q)     
      ffq = 1.0 + (dt/tcq->xi_Q) * deviation[tcq->eObs];
    else
      ffq=1.0;
    fq[tcq->at_i] *= ffq;
    
  }
  if (PAR(cr))
    gprod(cr,idef->atnr,fq);
  
  for(j=0; (j<md->nr); j++) {
    md->chargeA[j] *= fq[md->typeA[j]];
  }
  
  for(i=0; (i<tcr->nIP); i++) {
    tip    = &(tcr->tIP[i]);
    type   = tip->type;
    ftype  = idef->functype[type];
    factor = dt*deviation[tip->eObs];
    
    /* Time for a killer macro */
#define DOIP(ip) if (tip->xi.##ip) idef->iparams[type].##ip *= (1+factor/tip->xi.##ip)
      
    switch(ftype) {
    case F_BONDS:
      DOIP(harmonic.krA);
      DOIP(harmonic.rA);
	break;
    default:
      break;
    }
#undef DOIP
    tip->iprint=idef->iparams[type];
  }
}

