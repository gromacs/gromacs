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

t_coupl_rec *init_coupling(FILE *log,int nfile,t_filenm fnm[],
			   t_commrec *cr,t_forcerec *fr,
			   t_mdatoms *md,t_idef *idef)
{
  int         i,nc,index,j;
  int         ati,atj;
  t_coupl_rec *tcr;
  t_coupl_LJ  *tclj;
  t_coupl_BU  *tcbu;
  t_coupl_Q   *tcq;
  
  snew(tcr,1);
  read_gct(opt2fn("-j",nfile,fnm),tcr);
  write_gct(opt2fn("-jo",nfile,fnm),tcr,idef);

  copy_ff(tcr,fr,md,idef);
    
  /* Update all processors with coupling info */
  if (PAR(cr))
    comm_tcr(log,cr,&tcr);
  
  return tcr;
}

static void pr_ff(t_coupl_rec *tcr,real time,t_idef *idef,real ener[])
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
    prop=xvgropen("runaver.xvg","Properties and Running Averages",
		  "Time (ps)","");
    xvgr_legend(prop,asize(pleg),pleg);
    if (tcr->nLJ) {
      snew(out,tcr->nLJ);
      for(i=0; (i<tcr->nLJ); i++) {
	if (tcr->tcLJ[i].bPrint) {
	  tclj=&(tcr->tcLJ[i]);
	  sprintf(buf,"gctLJ%d-%d.xvg",tclj->at_i,tclj->at_j);
	  out[i]=xvgropen(buf,"General Coupling Lennard Jones","Time (ps)",
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
	  sprintf(buf,"gctBU%d-%d.xvg",tcbu->at_i,tcbu->at_j);
	  out[i]=xvgropen(buf,"General Coupling Buckingham","Time (ps)",
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
	sprintf(buf,"gctQ%d.xvg",tcr->tcQ[i].at_i);
	qq[i]=xvgropen(buf,"General Coupling Charge","Time (ps)","Charge (e)");
	fprintf(qq[i],"@ subtitle \"Type %d\"\n",tcr->tcQ[i].at_i);
	fflush(qq[i]);
      }
    }
    snew(ip,tcr->nIP);
    for(i=0; (i<tcr->nIP); i++) {
      sprintf(buf,"gctIP%d.xvg",tcr->tIP[i].type);
      ip[i]=xvgropen(buf,"General Coupling iparams","Time (ps)","ip ()");
      index=tcr->tIP[i].type;
      fprintf(ip[i],"@ subtitle \"Coupling to %s\"\n",
	      interaction_function[idef->functype[index]].longname);
      fflush(ip[i]);
    }
  }
  fprintf(prop,"%10g  %12.5e  %12.5e  %12.5e  %12.5e\n",
	  time,tcr->pres,ener[F_PRES],tcr->epot,ener[F_EPOT]);
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

static void pr_dev(real t,real dev[eoNR])
{
  static FILE *fp=NULL;
  int    i;
  
  if (!fp) {
    fp=xvgropen("deviatie.xvg","Deviations from target value",
		"Pres","Epot");
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
    /* 
       fprintf(log,"Updating FF for %d entries, f6=%g, f12=%g\n",atnr,f6,f12);
       */
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
    /* 
       fprintf(log,"Updating FF for %d entries, f6=%g, f12=%g\n",atnr,f6,f12);
       */
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

void do_coupling(FILE *log,t_coupl_rec *tcr,real t,int step,real ener[],
		 t_forcerec *fr,t_inputrec *ir,bool bMaster,
		 t_mdatoms *md,t_idef *idef,real mu_aver)
{
#define FMIN 0.95
#define FMAX 1.05

#define enm2Debye 48.0321
#define d2e(x) (x)/enm2Debye
#define enm2kjmol(x) (x)*0.0143952 /* = 2.0*4.0*M_PI*EPSILON0 */

  int        i,j,ati,atj,type,ftype,nmem,nmol;
  real       deviation[eoNR];
  real       f6,f12,fa,fb,fc,fq,factor,dt,mu_ind,Epol;
  bool       bTest,bPrint;
  t_coupl_LJ *tclj;
  t_coupl_BU *tcbu;
  t_coupl_Q  *tcq;
  t_coupl_iparams *tip;
  
  bPrint = do_per_step(step,ir->nstprint);
  dt     = ir->delta_t;
  nmol   = ir->userint3; /* Should contain the number of molecules */

  /* Use a memory of nmem steps, so we actually couple to the
   * average observable over the last nmem steps. This may help
   * in avoiding local minima in parameter space.
   */
  nmem=ir->userint2;

  if (step == 0) {
    tcr->pres=ener[F_PRES];
    tcr->epot=ener[F_EPOT];
  }

  tcr->pres=run_aver(tcr->pres,ener[F_PRES],step,nmem);
  tcr->epot=run_aver(tcr->epot,ener[F_EPOT],step,nmem);
  
  if (bMaster && bPrint)
    pr_ff(tcr,t,idef,ener);
  
  /* Calculate the deviation of average value from the target value */
  deviation[eoPres] = tcr->pres0 - tcr->pres;

  /* if dipole != 0.0 assume we want to use polarization corrected coupling */
  if ((tcr->dipole) != 0.0) {
    mu_ind = mu_aver - d2e(tcr->dipole); /* in e nm */
    /* Epol = mu_ind*mu_ind/(2.0*(tcr->polarizability)*4.0*M_PI*EPSILON0); */
    Epol = mu_ind*mu_ind/(enm2kjmol(tcr->polarizability));

    deviation[eoEpot] = (tcr->epot0)*nmol - Epol*nmol - tcr->epot;
    /*
    fprintf(stderr,":  %g %g %g %g %g \n",mu_ind,mu_aver, d2e(tcr->dipole),
	    d2e((2.27 - tcr->dipole)),(2.27 - tcr->dipole)
	   );
    fprintf(stderr,"Eref %g Epol %g Erunav %g Dev %g Eact %g\n",
	    (tcr->epot0)*nmol, Epol*nmol, tcr->epot,
	    deviation[eoEpot],ener[F_EPOT]);
	    */
  }
  else 
    deviation[eoEpot] = (tcr->epot0)*nmol - tcr->epot;
  
  /* Start coupling only after we have the correct average over 
   * nmem points.
   */
  if (step < nmem) 
    return;

  if (bPrint)
    pr_dev(t,deviation);
      
  for(i=0; (i<tcr->nLJ); i++) {
    tclj=&(tcr->tcLJ[i]);

    factor=deviation[tclj->eObs];
    
    f6=f12=1.0;
    if (tclj->xi_6)      
      f6  += (dt/tclj->xi_6)  * factor;
    if (tclj->xi_12)     
      f12 += (dt/tclj->xi_12) * factor;
      
    ati=tclj->at_i;
    atj=tclj->at_j;
    
    f6  = min(max(f6, FMIN),FMAX);
    f12 = min(max(f12,FMIN),FMAX);
    
    upd_nbfplj(log,tclj,fr->nbfp,ati,atj,idef->atnr,f6,f12);
  }
  
  for(i=0; (i<tcr->nBU); i++) {
    tcbu=&(tcr->tcBU[i]);

    factor=deviation[tcbu->eObs];
    
    fa=fb=fc=1.0;
    if (tcbu->xi_a)      
      fa  += (dt/tcbu->xi_a)  * factor;
    if (tcbu->xi_b)      
      fb  += (dt/tcbu->xi_b)  * factor;
    if (tcbu->xi_c)      
      fc  += (dt/tcbu->xi_c)  * factor;
      
    ati=tcbu->at_i;
    atj=tcbu->at_j;
    
    fa  = min(max(fa, FMIN),FMAX);
    fb  = min(max(fb, FMIN),FMAX);
    fc  = min(max(fc, FMIN),FMAX);
    
    upd_nbfpbu(log,tcbu,fr->nbfp,ati,atj,idef->atnr,fa,fb,fc);
  }
  
  for(i=0; (i<tcr->nQ); i++) {
    tcq=&(tcr->tcQ[i]);

    if (tcq->xi_Q)     
      fq = 1.0 + (dt/tcq->xi_Q) * deviation[tcq->eObs];
    else
      fq=1.0;
    fq  = min(max(fq, FMIN),FMAX);
    
    for(j=0; (j<md->nr); j++) 
      if (md->typeA[j] == tcq->at_i) {
	md->chargeA[j] *= fq;
	tcq->Q=md->chargeA[j];
      }
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


