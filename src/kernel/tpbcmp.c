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
static char *SRCID_tpbcmp_c = "$Id$";

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "main.h"
#include "macros.h"
#include "smalloc.h"
#include "futil.h"
#include "statutil.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "fatal.h"
#include "names.h"
#include "tpxio.h"
#include "enxio.h"

static void cmp_int(FILE *fp,char *s,int index,int i1,int i2)
{
  if (i1 != i2) {
    if (index != -1)
      fprintf(fp,"%s[%d] (%d - %d)\n",s,index,i1,i2);
    else
      fprintf(fp,"%s (%d - %d)\n",s,i1,i2);
  }
}

static void cmp_us(FILE *fp,char *s,int index,unsigned short i1,unsigned short i2)
{
  if (i1 != i2) {
    if (index != -1)
      fprintf(fp,"%s[%d] (%d - %d)\n",s,index,i1,i2);
    else
      fprintf(fp,"%s (%d - %d)\n",s,i1,i2);
  }
}

static void cmp_uc(FILE *fp,char *s,int index,unsigned char i1,unsigned char i2)
{
  if (i1 != i2) {
    if (index != -1)
      fprintf(fp,"%s[%d] (%d - %d)\n",s,index,i1,i2);
    else
      fprintf(fp,"%s (%d - %d)\n",s,i1,i2);
  }
}

static void cmp_real(FILE *fp,char *s,int index,real i1,real i2)
{
  if (i1 != i2) {
    if (index != -1)
      fprintf(fp,"%s[%d] (%e - %e)\n",s,index,i1,i2);
    else
      fprintf(fp,"%s (%e - %e)\n",s,i1,i2);
  }
}

static void cmp_rvec(FILE *fp,char *s,int index,rvec i1,rvec i2)
{
  if ((i1[XX] != i2[XX]) || (i1[YY] != i2[YY]) || (i1[ZZ] != i2[ZZ])) {
    if (index != -1)
      fprintf(fp,"%s[%d] (%e,%e,%e - %e,%e,%e)\n",s,index,
	      i1[XX],i1[YY],i1[ZZ],i2[XX],i2[YY],i2[ZZ]);
    else
      fprintf(fp,"%s (%e,%e,%e - %e,%e,%e)\n",s,
	      i1[XX],i1[YY],i1[ZZ],i2[XX],i2[YY],i2[ZZ]);
  }
}

static void cmp_ivec(FILE *fp,char *s,int index,ivec i1,ivec i2)
{
  if ((i1[XX] != i2[XX]) || (i1[YY] != i2[YY]) || (i1[ZZ] != i2[ZZ])) {
    if (index != -1)
      fprintf(fp,"%s[%d] (%d,%d,%d - %d,%d,%d)\n",s,index,
	      i1[XX],i1[YY],i1[ZZ],i2[XX],i2[YY],i2[ZZ]);
    else
      fprintf(fp,"%s (%d,%d,%d - %d,%d,%d)\n",s,
	      i1[XX],i1[YY],i1[ZZ],i2[XX],i2[YY],i2[ZZ]);
  }
}

static void cmp_ilist(FILE *fp,int ftype,t_ilist *il1,t_ilist *il2)
{
  int i;
  char buf[256];
 
  fprintf(fp,"comparing ilist %s\n",interaction_function[ftype].name);
  sprintf(buf,"%s->nr",interaction_function[ftype].name);
  cmp_int(fp,buf,0,il1->nr,il2->nr);
  sprintf(buf,"%s->multinr",interaction_function[ftype].name);
  for(i=0; (i<MAXNODES); i++)
    cmp_int(fp,buf,i,il1->multinr[i],il2->multinr[i]);
  sprintf(buf,"%s->iatoms",interaction_function[ftype].name);
  for(i=0; (i<il1->nr); i++) 
    cmp_int(fp,buf,i,il1->iatoms[i],il2->iatoms[i]);
}

void cmp_iparm(FILE *fp,char *s,int index,t_functype ft,
	       t_iparams ip1,t_iparams ip2) 
{
  if (memcmp(&ip1,&ip2,(size_t)sizeof(ip1)) != 0) {
    fprintf(fp,"%s1: ",s);
    pr_iparams(fp,ft,&ip1);
    fprintf(fp,"%s2: ",s);
    pr_iparams(fp,ft,&ip2);
  }
}

static void cmp_idef(FILE *fp,t_idef *id1,t_idef *id2)
{
  int i;
 
  fprintf(fp,"comparing idef\n");
  cmp_int(fp,"idef->ntypes",-1,id1->ntypes,id2->ntypes);
  cmp_int(fp,"idef->nodeid",   -1,id1->nodeid,id2->nodeid);
  cmp_int(fp,"idef->atnr",  -1,id1->atnr,id2->atnr);
  for(i=0; (i<id1->ntypes); i++) {
    cmp_int(fp,"idef->functype",i,(int)id1->functype[i],(int)id2->functype[i]);
    cmp_iparm(fp,"idef->iparam",i,id1->functype[i],
	      id1->iparams[i],id2->iparams[i]);
  }
  for(i=0; (i<F_NRE); i++)
    cmp_ilist(fp,i,&(id1->il[i]),&(id2->il[i]));
}

static void cmp_block(FILE *fp,t_block *b1,t_block *b2,char *s)
{
  int i,j,k;
  char buf[32];
  
  fprintf(fp,"comparing block %s\n",s);
  sprintf(buf,"%s.nr",s);
  cmp_int(fp,buf,-1,b1->nr,b2->nr);
  sprintf(buf,"%s.nra",s);
  cmp_int(fp,buf,-1,b1->nra,b2->nra);
  sprintf(buf,"%s.multinr",s);
  for(i=0; (i<MAXNODES); i++)
    cmp_int(fp,buf,i,b1->multinr[i],b2->multinr[i]);
} 

static void cmp_atom(FILE *fp,int index,t_atom *a1,t_atom *a2)
{
  int  i;
  char buf[256];
  
  cmp_us(fp,"atom.type",index,a1->type,a2->type);
  cmp_us(fp,"atom.ptype",index,a1->ptype,a2->ptype);
  cmp_int(fp,"atom.resnr",index,a1->resnr,a2->resnr);
  cmp_real(fp,"atom.m",index,a1->m,a2->m);
  cmp_real(fp,"atom.q",index,a1->q,a2->q);
  cmp_us(fp,"atom.typeB",index,a1->typeB,a2->typeB);
  cmp_real(fp,"atom.mB",index,a1->mB,a2->mB);
  cmp_real(fp,"atom.qB",index,a1->qB,a2->qB);
  for(i=0; (i<egcNR); i++) {
    sprintf(buf,"atom.grpnr(%d)",i);
    cmp_uc(fp,buf,index,a1->grpnr[i],a2->grpnr[i]);
  }
}

static void cmp_atoms(FILE *fp,t_atoms *a1,t_atoms *a2)
{
  int i;
  
  fprintf(fp,"comparing atoms\n");
  cmp_int(fp,"atoms->nr",-1,a1->nr,a2->nr);
  for(i=0; (i<a1->nr); i++)
    cmp_atom(fp,i,&(a1->atom[i]),&(a2->atom[i]));
  cmp_block(fp,&a1->excl,&a2->excl,"excl");
}

static void cmp_top(FILE *fp,t_topology *t1,t_topology *t2)
{
  int i;
  
  fprintf(fp,"comparing top\n");
  cmp_idef(fp,&(t1->idef),&(t2->idef));
  cmp_atoms(fp,&(t1->atoms),&(t2->atoms));
  for(i=0; (i<ebNR); i++)
    cmp_block(fp,&t1->blocks[i],&t2->blocks[i],EBLOCKS(i));
}

static void cmp_rvecs(FILE *fp,char *title,int n,rvec x1[],rvec x2[])
{
  int i;
  
  fprintf(fp,"comparing %s\n",title);
  for(i=0; (i<n); i++)
    cmp_rvec(fp,title,i,x1[i],x2[i]);
}

static void cmp_grpopts(FILE *fp,t_grpopts *opt1,t_grpopts *opt2)
{
  int i;
  
  cmp_int(fp,"inputrec->grpopts.ngtc",0,  opt1->ngtc,opt2->ngtc);
  cmp_int(fp,"inputrec->grpopts.ngacc",0, opt1->ngacc,opt2->ngacc);
  cmp_int(fp,"inputrec->grpopts.ngfrz",0, opt1->ngfrz,opt2->ngfrz);
  cmp_int(fp,"inputrec->grpopts.ngener",0,opt1->ngener,opt2->ngener);
  for(i=0; (i<min(opt1->ngtc,opt2->ngtc)); i++) {
    cmp_real(fp,"inputrec->grpopts.nrdf",i,opt1->nrdf[i],opt2->nrdf[i]);
    cmp_real(fp,"inputrec->grpopts.ref_t",i,opt1->ref_t[i],opt2->ref_t[i]);
    cmp_real(fp,"inputrec->grpopts.tau_t",i,opt1->tau_t[i],opt2->tau_t[i]);
  }
  for(i=0; (i<min(opt1->ngacc,opt2->ngacc)); i++)
    cmp_rvec(fp,"inputrec->grpopts.acc",i,opt1->acc[i],opt2->acc[i]);
  for(i=0; (i<min(opt1->ngfrz,opt2->ngfrz)); i++)
    cmp_ivec(fp,"inputrec->grpopts.nFreeze",i,opt1->nFreeze[i],opt2->nFreeze[i]);
}

static void cmp_cosines(FILE *fp,char *s,t_cosines c1[DIM],t_cosines c2[DIM])
{
  int i,m;
  char buf[256];
  
  for(m=0; (m<DIM); m++) {
    sprintf(buf,"inputrec->%s[%d]",s,m);
    cmp_int(fp,buf,0,c1->n,c2->n);
    for(i=0; (i<min(c1->n,c2->n)); i++) {
      cmp_real(fp,buf,i,c1->a[i],c2->a[i]);
      cmp_real(fp,buf,i,c1->phi[i],c2->phi[i]);
    }
  }
}

static void cmp_inputrec(FILE *fp,t_inputrec *ir1,t_inputrec *ir2)
{
  fprintf(fp,"comparing inputrec\n");

  /* gcc 2.96 doesnt like these defines at all, but issues a huge list
   * of warnings. Maybe it will change in future versions, but for the
   * moment I've spelled them out instead. /EL 000820 
   * #define CIB(s) cmp_int(fp,"inputrec->"#s,0,ir1->##s,ir2->##s)
   * #define CII(s) cmp_int(fp,"inputrec->"#s,0,ir1->##s,ir2->##s)
   * #define CIR(s) cmp_real(fp,"inputrec->"#s,0,ir1->##s,ir2->##s)
   */
  cmp_int(fp,"inputrec->eI",0,ir1->eI,ir2->eI);
  cmp_int(fp,"inputrec->nsteps",0,ir1->nsteps,ir2->nsteps);
  cmp_int(fp,"inputrec->ePBC",0,ir1->ePBC,ir2->ePBC);
  cmp_int(fp,"inputrec->ns_type",0,ir1->ns_type,ir2->ns_type);
  cmp_int(fp,"inputrec->nstlist",0,ir1->nstlist,ir2->nstlist);
  cmp_int(fp,"inputrec->ndelta",0,ir1->ndelta,ir2->ndelta);
  cmp_int(fp,"inputrec->bDomDecomp",0,ir1->bDomDecomp,ir2->bDomDecomp);
  cmp_int(fp,"inputrec->decomp_dir",0,ir1->decomp_dir,ir2->decomp_dir);
  cmp_int(fp,"inputrec->nstcomm",0,ir1->nstcomm,ir2->nstcomm);
  cmp_int(fp,"inputrec->nstlog",0,ir1->nstlog,ir2->nstlog);
  cmp_int(fp,"inputrec->nstxout",0,ir1->nstxout,ir2->nstxout);
  cmp_int(fp,"inputrec->nstvout",0,ir1->nstvout,ir2->nstvout);
  cmp_int(fp,"inputrec->nstfout",0,ir1->nstfout,ir2->nstfout);
  cmp_int(fp,"inputrec->nstenergy",0,ir1->nstenergy,ir2->nstenergy);
  cmp_int(fp,"inputrec->nstxtcout",0,ir1->nstxtcout,ir2->nstxtcout);
  cmp_real(fp,"inputrec->init_t",0,ir1->init_t,ir2->init_t);
  cmp_real(fp,"inputrec->delta_t",0,ir1->delta_t,ir2->delta_t);
  cmp_real(fp,"inputrec->xtcprec",0,ir1->xtcprec,ir2->xtcprec);
  cmp_int(fp,"inputrec->nkx",0,ir1->nkx,ir2->nkx);
  cmp_int(fp,"inputrec->nky",0,ir1->nky,ir2->nky);
  cmp_int(fp,"inputrec->nkz",0,ir1->nkz,ir2->nkz);
  cmp_int(fp,"inputrec->pme_order",0,ir1->pme_order,ir2->pme_order);
  cmp_real(fp,"inputrec->ewald_rtol",0,ir1->ewald_rtol,ir2->ewald_rtol);
  cmp_real(fp,"inputrec->epsilon_surface",0,ir1->epsilon_surface,ir2->epsilon_surface);
  cmp_int(fp,"inputrec->bOptFFT",0,ir1->bOptFFT,ir2->bOptFFT);
  cmp_int(fp,"inputrec->bUncStart",0,ir1->bUncStart,ir2->bUncStart);
  cmp_int(fp,"inputrec->etc",0,ir1->etc,ir2->etc);
  cmp_int(fp,"inputrec->epc",0,ir1->epc,ir2->epc);
  cmp_int(fp,"inputrec->epct",0,ir1->epct,ir2->epct);
  cmp_real(fp,"inputrec->tau_p",0,ir1->tau_p,ir2->tau_p);
  cmp_rvec(fp,"inputrec->ref_p(x)",0,ir1->ref_p[XX],ir2->ref_p[XX]);
  cmp_rvec(fp,"inputrec->ref_p(y)",0,ir1->ref_p[YY],ir2->ref_p[YY]);
  cmp_rvec(fp,"inputrec->ref_p(z)",0,ir1->ref_p[ZZ],ir2->ref_p[ZZ]);
  cmp_rvec(fp,"inputrec->compress(x)",0,ir1->compress[XX],ir2->compress[XX]);
  cmp_rvec(fp,"inputrec->compress(y)",0,ir1->compress[YY],ir2->compress[YY]);
  cmp_rvec(fp,"inputrec->compress(z)",0,ir1->compress[ZZ],ir2->compress[ZZ]);
  cmp_int(fp,"inputrec->bSimAnn",0,ir1->bSimAnn,ir2->bSimAnn);
  cmp_real(fp,"inputrec->zero_temp_time",0,ir1->zero_temp_time,ir2->zero_temp_time);
  cmp_real(fp,"inputrec->rlist",0,ir1->rlist,ir2->rlist);
  cmp_int(fp,"inputrec->coulombtype",0,ir1->coulombtype,ir2->coulombtype);
  cmp_real(fp,"inputrec->rcoulomb_switch",0,ir1->rcoulomb_switch,ir2->rcoulomb_switch);
  cmp_real(fp,"inputrec->rcoulomb",0,ir1->rcoulomb,ir2->rcoulomb);
  cmp_int(fp,"inputrec->vdwtype",0,ir1->vdwtype,ir2->vdwtype);
  cmp_real(fp,"inputrec->rvdw_switch",0,ir1->rvdw_switch,ir2->rvdw_switch);
  cmp_real(fp,"inputrec->rvdw",0,ir1->rvdw,ir2->rvdw);
  cmp_real(fp,"inputrec->epsilon_r",0,ir1->epsilon_r,ir2->epsilon_r);
  cmp_int(fp,"inputrec->bDispCorr",0,ir1->bDispCorr,ir2->bDispCorr);
  cmp_real(fp,"inputrec->shake_tol",0,ir1->shake_tol,ir2->shake_tol);
  cmp_real(fp,"inputrec->fudgeQQ",0,ir1->fudgeQQ,ir2->fudgeQQ);
  cmp_int(fp,"inputrec->efep",0,ir1->efep,ir2->efep);
  cmp_real(fp,"inputrec->init_lambda",0,ir1->init_lambda,ir2->init_lambda);
  cmp_real(fp,"inputrec->delta_lambda",0,ir1->delta_lambda,ir2->delta_lambda);
  cmp_real(fp,"inputrec->sc_alpha",0,ir1->sc_alpha,ir2->sc_alpha);
  cmp_real(fp,"inputrec->sc_sigma",0,ir1->sc_sigma,ir2->sc_sigma);
  cmp_real(fp,"inputrec->dr_fc",0,ir1->dr_fc,ir2->dr_fc);
  cmp_int(fp,"inputrec->eDisreWeighting",0,ir1->eDisreWeighting,ir2->eDisreWeighting);
  cmp_int(fp,"inputrec->bDisreMixed",0,ir1->bDisreMixed,ir2->bDisreMixed);
  cmp_int(fp,"inputrec->nstdisreout",0,ir1->nstdisreout,ir2->nstdisreout);
  cmp_real(fp,"inputrec->dr_tau",0,ir1->dr_tau,ir2->dr_tau);
  cmp_real(fp,"inputrec->em_stepsize",0,ir1->em_stepsize,ir2->em_stepsize);
  cmp_real(fp,"inputrec->em_tol",0,ir1->em_tol,ir2->em_tol);
  cmp_int(fp,"inputrec->nstcgsteep",0,ir1->nstcgsteep,ir2->nstcgsteep);
  cmp_int(fp,"inputrec->eConstrAlg",0,ir1->eConstrAlg,ir2->eConstrAlg);
  cmp_int(fp,"inputrec->nProjOrder",0,ir1->nProjOrder,ir2->nProjOrder);
  cmp_real(fp,"inputrec->LincsWarnAngle",0,ir1->LincsWarnAngle,ir2->LincsWarnAngle);
  cmp_real(fp,"inputrec->ld_temp",0,ir1->ld_temp,ir2->ld_temp);
  cmp_real(fp,"inputrec->ld_fric",0,ir1->ld_fric,ir2->ld_fric);
  cmp_int(fp,"inputrec->ld_seed",0,ir1->ld_seed,ir2->ld_seed);
  cmp_real(fp,"inputrec->cos_accel",0,ir1->cos_accel,ir2->cos_accel);
  cmp_int(fp,"inputrec->userint1",0,ir1->userint1,ir2->userint1);
  cmp_int(fp,"inputrec->userint2",0,ir1->userint2,ir2->userint2);
  cmp_int(fp,"inputrec->userint3",0,ir1->userint3,ir2->userint3);
  cmp_int(fp,"inputrec->userint4",0,ir1->userint4,ir2->userint4);
  cmp_real(fp,"inputrec->userreal1",0,ir1->userreal1,ir2->userreal1);
  cmp_real(fp,"inputrec->userreal2",0,ir1->userreal2,ir2->userreal2);
  cmp_real(fp,"inputrec->userreal3",0,ir1->userreal3,ir2->userreal3);
  cmp_real(fp,"inputrec->userreal4",0,ir1->userreal4,ir2->userreal4);
  cmp_grpopts(fp,&(ir1->opts),&(ir2->opts));
  cmp_cosines(fp,"ex",ir1->ex,ir2->ex);
  cmp_cosines(fp,"et",ir1->et,ir2->et);
}


void comp_tpx(char *fn1,char *fn2)
{
  char        *ff[2];
  t_tpxheader sh[2];
  t_inputrec  ir[2];
  rvec        *xx[2],*vv[2];
  t_topology  top[2];
  matrix      box[2];
  int         i,step,natoms;
  real        t,lambda;

  ff[0]=fn1;
  ff[1]=fn2;
  for(i=0; (i<2); i++) {
    read_tpxheader(ff[i],&(sh[i]));
    snew(xx[i],sh[i].natoms);
    snew(vv[i],sh[i].natoms);
    read_tpx(ff[i],&step,&t,&lambda,&(ir[i]),box[i],&natoms,
	     xx[i],vv[i],NULL,&(top[i]));
  }
  cmp_inputrec(stdout,&ir[0],&ir[1]);
  cmp_top(stdout,&top[0],&top[1]);
  cmp_rvecs(stdout,"box",DIM,box[0],box[1]);
  cmp_rvecs(stdout,"x",natoms,xx[0],xx[1]);
  cmp_rvecs(stdout,"v",natoms,vv[0],vv[1]);
}

static void cmp_energies(FILE *fp,int frame1,int frame2,int nre,
			 t_energy e1[],t_energy e2[],char *enm1[],char *enm2[],
			 real ftol)
{
  int  i;
  
  for(i=0; (i<nre); i++)
    if (fabs(e1[i].e - e2[i].e) > ftol)
      fprintf(fp,"%-15s  frame %3d:  %12g,  frame %3d: %12g\n",
	      enm1[i],frame1,e1[i].e,frame2,e2[i].e);
}

void comp_enx(char *fn1,char *fn2,real ftol)
{
  int       in1,in2,nre,nre1,nre2,frame1,frame2,ndr1,ndr2;
  int       i;
  char      **enm1=NULL,**enm2=NULL;
  t_energy  *ee1,*ee2;
  t_drblock dr1,dr2;
  bool      b1,b2;
  real      t1,t2;
  
  fprintf(stdout,"comparing energy file %s and %s\n\n",fn1,fn2);

  in1 = open_enx(fn1,"r");
  in2 = open_enx(fn2,"r");
  do_enxnms(in1,&nre1,&enm1);
  do_enxnms(in2,&nre2,&enm2);
  if (nre1 != nre2) {
    fprintf(stdout,"%s: nre=%d, %s: nre=%d\n",fn1,nre1,fn2,nre2);
    return;
  }
  nre = nre1;
  fprintf(stdout,"There are %d terms in the energy files\n\n",nre);
  
  for(i=0; (i<nre); i++) {
    if (strcmp(enm1[i],enm2[i]) != 0)
      fprintf(stdout,"%s: %16s - %s: %16s\n",fn1,enm1[i],fn2,enm2[i]);
  }
  
  snew(ee1,nre);
  snew(ee2,nre);
  do { 
    b1 = do_enx(in1,&t1,&frame1,&nre1,ee1,&ndr1,&dr1);
    b2 = do_enx(in2,&t2,&frame2,&nre2,ee2,&ndr2,&dr2);
    if (b1 && !b2)
      fprintf(stdout,"\nEnd of file on %s but not on %s\n",fn2,fn1);
    else if (!b1 && b2) 
      fprintf(stdout,"\nEnd of file on %s but not on %s\n",fn1,fn2);
    else if (!b1 && !b2)
      fprintf(stdout,"\nFiles read succesfully\n");
    else {
      cmp_real(stdout,"t",-1,t1,t2);
      cmp_int(stdout,"frame",-1,frame1,frame2);
      cmp_int(stdout,"nre",-1,nre1,nre2);
      cmp_int(stdout,"ndr",-1,ndr1,ndr2);

      if (nre1 != nre) {
	fprintf(stdout,
		"file %s internally inconsistent: nre changed from %d to %d\n",
		fn1,nre,nre1);
	break;
      }
      if (nre2 != nre) {
	fprintf(stdout,
		"file %s internally inconsistent: nre changed from %d to %d\n",
		fn2,nre,nre2);
	break;
      }
      cmp_energies(stdout,frame1,frame2,nre,ee1,ee2,enm1,enm2,ftol);
      /*if ((ndr1 == ndr2) && (ndr1 > 0))
	cmp_disres(stdout,frame1,frame2,ndr1,dr1,dr2);*/
    }
  } while (b1 && b2);
    
  close_enx(in1);
  close_enx(in2);

}
