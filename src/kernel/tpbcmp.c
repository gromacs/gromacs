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
 * S  C  A  M  O  R  G
 */
static char *SRCID_tpbcmp_c = "$Id$";

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

static void cmp_int(FILE *fp,char *s,int index,int i1,int i2)
{
  if (i1 != i2) {
    if (index != -1)
      fprintf(fp,"%s[%d] (%d - %d)\n",s,index,i1,i2);
    else
      fprintf(fp,"%s (%d - %d)\n",s,i1,i2);
  }
}

static void cmp_us(FILE *fp,char *s,int index,ushort i1,ushort i2)
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
  
  sprintf(buf,"%s->nr",interaction_function[ftype].name);
  cmp_int(fp,buf,0,il1->nr,il2->nr);
  sprintf(buf,"%s->multinr",interaction_function[ftype].name);
  for(i=0; (i<MAXPROC); i++)
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
  cmp_int(fp,"idef->pid",   -1,id1->pid,id2->pid);
  cmp_int(fp,"idef->atnr",  -1,id1->atnr,id2->atnr);
  for(i=0; (i<id1->ntypes); i++) {
    cmp_int(fp,"idef->functype",i,(int)id1->functype[i],(int)id2->functype[i]);
    cmp_iparm(fp,"idef->iparam",i,id1->functype[i],
	      id1->iparams[i],id2->iparams[i]);
  }
  for(i=0; (i<F_NRE); i++)
    cmp_ilist(fp,i,&(id1->il[i]),&(id2->il[i]));
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
}

static void cmp_top(FILE *fp,t_topology *t1,t_topology *t2)
{
  fprintf(fp,"comparing top\n");
  cmp_idef(fp,&(t1->idef),&(t2->idef));
  cmp_atoms(fp,&(t1->atoms),&(t2->atoms));
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
    cmp_int(fp,"inputrec->grpopts.nrdf",i,opt1->nrdf[i],opt2->nrdf[i]);
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
#define CIB(s) cmp_int(fp,"inputrec->"#s,0,ir1->##s,ir2->##s)
#define CII(s) cmp_int (fp,"inputrec->"#s,0,ir1->##s,ir2->##s)
#define CIR(s) cmp_real(fp,"inputrec->"#s,0,ir1->##s,ir2->##s)
  CII(eI);
  CII(nsteps);
  CII(ns_type);
  CII(nstlist);
  CII(ndelta);
  CII(nstcomm);
  CII(nstlog);
  CII(nstxout);
  CII(nstvout);
  CII(nstfout);
  CII(nstenergy);
  CII(nstxtcout);
  CIR(init_t);
  CIR(delta_t);
  CIR(xtcprec);
  CII(niter);
  CII(solvent_opt);
  CIR(gausswidth);
  CII(nkx);
  CII(nky);
  CII(nkz);
  CII(eBox);
  CIB(bUncStart);
  CIB(btc);
  CII(epc);
  CII(ntcmemory);
  CII(npcmemory);
  CIR(tau_p);
  cmp_rvec(fp,"inputrec->ref_p",0,ir1->ref_p,ir2->ref_p);
  cmp_rvec(fp,"inputrec->compress",0,ir1->compress,ir2->compress);
  CIB(bSimAnn);
  CIR(zero_temp_time);
  CII(eeltype);
  CIR(rshort);
  CIR(rlong);
  CIR(ns_dr);
  CII(bLJshift);
  CII(bLJcorr);
  CIR(epsilon_r);
  CIR(shake_tol);
  CIR(fudgeLJ);
  CIR(fudgeQQ);
  CIB(bPert);
  CIR(init_lambda);
  CIR(delta_lambda);
  CII(eDisreWeighting);
  CIB(bDisreMixed);
  CIR(dr_fc);
  CIR(dr_tau);
  CII(nstdisreout);
  CIR(em_stepsize);
  CIR(em_tol);
  CII(nstcgsteep);
  CII(eConstrAlg);
  CII(nProjOrder);
  CIR(LincsWarnAngle);
  CII(nstLincsout);
  CIR(ld_temp);
  CIR(ld_fric);
  CII(ld_seed);
  CII(userint1);
  CII(userint2);
  CII(userint3);
  CII(userint4);
  CIR(userreal1);
  CIR(userreal2);
  CIR(userreal3);
  CIR(userreal4);
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
