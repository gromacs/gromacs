/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * Great Red Oystrich Makes All Chemists Sane
 */

#include <stdio.h>
#include "smalloc.h"
#include "typedefs.h"
#include "names.h"
#include "txtdump.h"
#include "string2.h"

int available(FILE *fp,void *p,char *title)
{
  if (!p) (void) fprintf(fp,"not available: %s\n",title);
  return (p!=NULL);
}

int pr_indent(FILE *fp,int n)
{
  int i;

  for (i=0; i<n; i++) (void) fprintf(fp," ");
  return n;
}

int pr_title(FILE *fp,int indent,char *title)
{
  (void) pr_indent(fp,indent);
  (void) fprintf(fp,"%s:\n",title);
  return (indent+INDENT);
}

int pr_title_n(FILE *fp,int indent,char *title,int n)
{
  (void) pr_indent(fp,indent);
  (void) fprintf(fp,"%s (%d):\n",title,n);
  return (indent+INDENT);
}

int pr_title_nxn(FILE *fp,int indent,char *title,int n1,int n2)
{
  (void) pr_indent(fp,indent);
  (void) fprintf(fp,"%s (%dx%d):\n",title,n1,n2);
  return (indent+INDENT);
}

void pr_ivec(FILE *fp,int indent,char *title,int vec[],int n)
{
  int i;

  if (available(fp,vec,title))
    {
      indent=pr_title_n(fp,indent,title,n);
      for (i=0; i<n; i++)
        {
          (void) pr_indent(fp,indent);
          (void) fprintf(fp,"%s[%d]=%d\n",title,i,vec[i]);
        }
    }
}

void pr_ivecs(FILE *fp,int indent,char *title,ivec vec[],int n)
{
  int i,j;

  if (available(fp,vec,title))
    {  
      indent=pr_title_nxn(fp,indent,title,n,DIM);
      for (i=0; i<n; i++)
        {
          (void) pr_indent(fp,indent);
          (void) fprintf(fp,"%s[%d]={",title,i);
          for (j=0; j<DIM; j++)
            {
              if (j!=0) (void) fprintf(fp,", ");
              fprintf(fp,"%d",vec[i][j]);
            }
          (void) fprintf(fp,"}\n");
        }
    }
}

void pr_rvec(FILE *fp,int indent,char *title,real vec[],int n)
{
  int i;

  if (available(fp,vec,title))
    {  
      indent=pr_title_n(fp,indent,title,n);
      for (i=0; i<n; i++)
        {
          (void) pr_indent(fp,indent);
          (void) fprintf(fp,"%s[%d]=%12.5e\n",title,i,vec[i]);
        }
    }
}

void pr_rvecs(FILE *fp,int indent,char *title,rvec vec[],int n)
{
  int i,j;

  if (available(fp,vec,title))
    {  
      indent=pr_title_nxn(fp,indent,title,n,DIM);
      for (i=0; i<n; i++)
        {
          (void) pr_indent(fp,indent);
          (void) fprintf(fp,"%s[%5d]={",title,i);
          for (j=0; j<DIM; j++)
            {
              if (j!=0) (void) fprintf(fp,", ");
              fprintf(fp,"%12.5e",vec[i][j]);
            }
          (void) fprintf(fp,"}\n");
        }
    }
}

void pr_energies(FILE *fp,int indent,char *title,t_energy *e,int n)
{
  int i;

  if (available(fp,e,title)) {
    indent=pr_title_n(fp,indent,title,n);
    for (i=0; i<n; i++) {
      (void) pr_indent(fp,indent);
      fprintf(fp,"%s[%2d]={e=%10.3e, eav=%10.3e, esum=%10.3e, e2sum=%10.3e}\n",
	      title,i,e[i].e,e[i].eav,e[i].esum,e[i].e2sum);
    }
  }
}

void pr_grp_opts(FILE *out,int indent,char *title,t_grpopts *opts)
{
  int i,m;

  fprintf(out,"%s:\n",title);
  
  pr_indent(out,indent);
  fprintf(out,"nrdf:\t");
  for(i=0; (i<opts->ngtc); i++)
    fprintf(out,"  %10d",opts->nrdf[i]);
  fprintf(out,"\n");
  
  pr_indent(out,indent);
  fprintf(out,"ref_t:\t");
  for(i=0; (i<opts->ngtc); i++)
    fprintf(out,"  %10g",opts->ref_t[i]);
  fprintf(out,"\n");
  
  pr_indent(out,indent);
  fprintf(out,"tau_t:\t");
  for(i=0; (i<opts->ngtc); i++)
    fprintf(out,"  %10g",opts->tau_t[i]);
  fprintf(out,"\n");
  
  pr_indent(out,indent);
  fprintf(out,"acc:\t");
  for(i=0; (i<opts->ngacc); i++)
    for(m=0; (m<DIM); m++)
      fprintf(out,"  %10g",opts->acc[i][m]);
  fprintf(out,"\n");

  pr_indent(out,indent);
  fprintf(out,"nfreeze:");
  for(i=0; (i<opts->ngfrz); i++)
    for(m=0; (m<DIM); m++)
      fprintf(out,"  %10s",opts->nFreeze[i][m] ? "Y" : "N");
  fprintf(out,"\n");

  fflush(out);
}

static void pr_cosine(FILE *fp,int indent,char *title,t_cosines *cos)
{
  int j;
  
  indent=pr_title(fp,indent,title);
  (void) pr_indent(fp,indent);
  fprintf(fp,"n = %d\n",cos->n);
  if (cos->n > 0) {
    (void) pr_indent(fp,indent+2);
    fprintf(fp,"a =");
    for(j=0; (j<cos->n); j++)
      fprintf(fp," %e",cos->a[j]);
    fprintf(fp,"\n");
    (void) pr_indent(fp,indent+2);
    fprintf(fp,"phi =");
    for(j=0; (j<cos->n); j++)
      fprintf(fp," %e",cos->phi[j]);
    fprintf(fp,"\n");
  }
}

static void pr_int(FILE *fp,int indent,char *title,int i)
{
  pr_indent(fp,indent);
  fprintf(fp,"%-20s = %d\n",title,i);
}

static void pr_real(FILE *fp,int indent,char *title,real r)
{
  pr_indent(fp,indent);
  fprintf(fp,"%-20s = %g\n",title,r);
}

static void pr_str(FILE *fp,int indent,char *title,char *s)
{
  pr_indent(fp,indent);
  fprintf(fp,"%-20s = %s\n",title,s);
}

void pr_inputrec(FILE *fp,int indent,char *title,t_inputrec *ir)
{
  if (available(fp,ir,title)) {
    indent=pr_title(fp,indent,title);
#define PS(t,s) pr_str(fp,indent,t,s)
#define PI(t,s) pr_int(fp,indent,t,s)
#define PR(t,s) pr_real(fp,indent,t,s)
    PS("integrator",EI(ir->eI));
    PI("nsteps",ir->nsteps);
    PS("ns_type",ENS(ir->ns_type));
    PI("nstlist",ir->nstlist);
    PI("ndelta",ir->ndelta);
    PI("nstcomm",ir->nstcomm);
    PI("nstprint",ir->nstprint);
    PI("nstxout",ir->nstxout);
    PI("nstvout",ir->nstvout);
    PI("nstfout",ir->nstfout);
    PI("nstgrp",ir->nstgrp);
    PI("nstxtcout",ir->nstxtcout);
    PR("init_t",ir->init_t);
    PR("delta_t",ir->delta_t);
    PR("xtcprec",ir->xtcprec);
    PI("niter",ir->niter);
    PI("watertype",ir->watertype);
    PI("nwatoms",ir->nwatoms);
    PR("gausswidth",ir->gausswidth);
    PI("nkx",ir->nkx);
    PI("nky",ir->nky);
    PI("nkz",ir->nkz);
    PS("eBox",EBOXTYPE(ir->eBox));
    PS("bShakeFirst",BOOL(ir->bShakeFirst));
    PS("btc",BOOL(ir->btc));
    PI("ntcmemory",ir->ntcmemory);
    PS("bpc",BOOL(ir->bpc));
    PI("npcmemory",ir->npcmemory);
    PR("tau_p",ir->tau_p);
    pr_rvec(fp,indent,"ref_p",ir->ref_p,DIM);
    pr_rvec(fp,indent,"compress",ir->compress,DIM);
    PS("bSimAnn",BOOL(ir->bSimAnn)); 
    PR("zero_temp_time",ir->zero_temp_time); 
    PS("eeltype",EELTYPE(ir->eeltype));
    PR("rshort",ir->rshort);
    PR("rlong",ir->rlong);
    PR("epsilon_r",ir->epsilon_r);
    PS("bLJcorr",BOOL(ir->bLJcorr));
    PR("tol",ir->tol);
    PR("fudgeLJ",ir->fudgeLJ);
    PR("fudgeQQ",ir->fudgeQQ);
    PS("free_energy",BOOL(ir->bPert));
    PR("init_lambda",ir->init_lambda);
    PR("delta_lambda",ir->delta_lambda);
    PR("dr_fc",ir->dr_fc);
    PR("dihr_fc",ir->dihr_fc);
    PR("dr_tau",ir->dr_tau);
    PR("em_stepsize",ir->em_stepsize);
    PR("em_tol",ir->em_tol);
    PS("shake_type",ESHAKETYPE(ir->eShakeType));
    PR("shake-tol",ir->tol);
    PI("lincs_order",ir->nProjOrder);
    PR("lincs_warnangle",ir->LincsWarnAngle);
    PI("nstLincsout",ir->nstLincsout);
    PR("ld_temp",ir->ld_temp);
    PR("ld_fric",ir->ld_fric);
    PI("ld_seed",ir->ld_seed);
    PI("userint1",ir->userint1);
    PI("userint2",ir->userint2);
    PI("userint3",ir->userint3);
    PI("userint4",ir->userint4);
    PR("userreal1",ir->userreal1);
    PR("userreal2",ir->userreal2);
    PR("userreal3",ir->userreal3);
    PR("userreal4",ir->userreal4);
#undef PS
#undef PR
#undef PI
    pr_grp_opts(fp,indent,"grpopts",&(ir->opts));
    pr_cosine(fp,indent,"efield-x",&(ir->ex[XX]));
    pr_cosine(fp,indent,"efield-xt",&(ir->et[XX]));
    pr_cosine(fp,indent,"efield-y",&(ir->ex[YY]));
    pr_cosine(fp,indent,"efield-yt",&(ir->et[YY]));
    pr_cosine(fp,indent,"efield-z",&(ir->ex[ZZ]));
    pr_cosine(fp,indent,"efield-zt",&(ir->et[ZZ]));
  }
}

static void pr_harm(FILE *fp,t_iparams *iparams,char *r,char *kr)
{
  fprintf(fp,"%sA=%12.5e, %sA=%12.5e, %sB=%12.5e, %sB=%12.5e\n",
	  r,iparams->harmonic.rA,kr,iparams->harmonic.krA,
	  r,iparams->harmonic.rB,kr,iparams->harmonic.krB);
}

void pr_iparams(FILE *fp,int ftype,t_iparams *iparams)
{
  int i;
  
  switch (ftype) {
  case F_ANGLES:
    pr_harm(fp,iparams,"th","ct");
    break;
  case F_BHAM:
    fprintf(fp,"a=%15.8e, b=%15.8e, c=%15.8e\n",
	    iparams->bham.a,iparams->bham.b,iparams->bham.c);
    break;
  case F_BONDS:
    pr_harm(fp,iparams,"b0","cb");
    break;
  case F_IDIHS:
    pr_harm(fp,iparams,"xi","cx");
    break;
  case F_MORSE:
    fprintf(fp,"b0=%15.8e, cb=%15.8e, beta=%15.8e\n",
	    iparams->morse.b0,iparams->morse.cb,iparams->morse.beta);
    break;
  case F_LJ:
    fprintf(fp,"c6=%15.8e, c12=%15.8e\n",iparams->lj.c6,iparams->lj.c12);
    break;
  case F_LJ14:
    fprintf(fp,"c6A=%15.8e, c12A=%15.8e, c6B=%15.8e, c12B=%15.8e\n",
	    iparams->lj14.c6A,iparams->lj14.c12A,
	    iparams->lj14.c6B,iparams->lj14.c12B);
    break;
  case F_PDIHS:
    fprintf(fp,"phiA=%15.8e, cpA=%15.8e, phiB=%15.8e, cpB=%15.8e, mult=%d\n",
	    iparams->pdihs.phiA,iparams->pdihs.cpA,
	    iparams->pdihs.phiB,iparams->pdihs.cpB,
	    iparams->pdihs.mult);
    break;
  case F_DISRES:
    fprintf(fp,"index=%4d, type=%1d, rx=(%15.8e,%15.8e,%15.8e,%15.8e)\n",
	    iparams->disres.index,iparams->disres.type,
	    iparams->disres.rx0,iparams->disres.rx1,
	    iparams->disres.rx2,iparams->disres.rx3);
    break;
  case F_POSRES:
    fprintf(fp,"pos0=(%15.8e,%15.8e,%15.8e), fc=(%15.8e,%15.8e,%15.8e)\n",
	    iparams->posres.pos0[XX],iparams->posres.pos0[YY],
	    iparams->posres.pos0[ZZ],iparams->posres.fc[XX],
	    iparams->posres.fc[YY],iparams->posres.fc[ZZ]);
    break;
  case F_RBDIHS:
    for (i=0; i<NR_RBDIHS; i++)
      fprintf(fp,"%srbc[%d]=%15.8e",i==0?"":", ",i,iparams->rbdihs.rbc[i]);
    fprintf(fp,"\n");
    break;
  case F_SHAKE:
    fprintf(fp,"dA=%15.8e, dB=%15.8e\n",
	    iparams->shake.dA,iparams->shake.dB);
    break;
  case F_SETTLE:
    fprintf(fp,"doh=%15.8e, dhh=%15.8e\n",iparams->settle.doh,
	    iparams->settle.dhh);
    break;
  case F_DUMMY1:
    fprintf(fp,"a=%15.8e, b=%15.8e\n",
	    iparams->dummy.a,iparams->dummy.b);
    break;
  case F_DUMMY2:
  case F_DUMMY3:
    fprintf(fp,"a=%15.8e, b=%15.8e, c=%15.8e\n",
	    iparams->dummy.a,iparams->dummy.b,iparams->dummy.c);
    break;
  default:
    fatal_error(0,"unknown function type %d (%s) in %s line %d",
		ftype,interaction_function[ftype].name,__FILE__,__LINE__);
  }
}

static void pr_ilist(FILE *fp,int indent,char *title,
                     t_idef *idef,t_ilist *ilist)
{
  int i,j,k,type,ftype;
  t_iatom *iatoms;

  if (available(fp,ilist,title))
    {  
      indent=pr_title(fp,indent,title);
      (void) pr_indent(fp,indent);
      fprintf(fp,"multinr:");
      for (i=0; i<MAXPROC; i++) 
	(void) fprintf(fp," %d",ilist->multinr[i]);
      fprintf(fp,"\n");
      (void) pr_indent(fp,indent);
      fprintf(fp,"iatoms(%d):\n",ilist->nr);
      iatoms=ilist->iatoms;
      for (i=j=0; i<ilist->nr;) {
#ifndef DEBUG
	(void) pr_indent(fp,indent+INDENT);
	type=*(iatoms++);
	ftype=idef->functype[type];
	(void) fprintf(fp,"%d type=%d (%s)",
		       j++,type,interaction_function[ftype].name);
	for (k=0; k<interaction_function[ftype].nratoms; k++)
	  (void) fprintf(fp," %d",*(iatoms++));
	(void) fprintf(fp,"\n");
	i+=1+interaction_function[ftype].nratoms;
#else
	fprintf(fp,"%5d%5d\n",i,iatoms[i]);
	i++;
#endif
      }
    }
}

void pr_idef(FILE *fp,int indent,char *title,t_idef *idef)
{
  int i,j;
  
  if (available(fp,idef,title)) {  
    indent=pr_title(fp,indent,title);
    (void) pr_indent(fp,indent);
    (void) fprintf(fp,"atnr=%d\n",idef->atnr);
    (void) pr_indent(fp,indent);
    (void) fprintf(fp,"pid=%d\n",idef->pid);
    (void) pr_indent(fp,indent);
    (void) fprintf(fp,"ntypes=%d\n",idef->ntypes);
    for (i=0; i<idef->ntypes; i++) {
      (void) pr_indent(fp,indent+INDENT);
      (void) fprintf(fp,"functype[%d]=%s, ",
		     i,interaction_function[idef->functype[i]].name);
      pr_iparams(fp,idef->functype[i],&idef->iparams[i]);
    }
    for(j=0; (j<F_NRE); j++)
      pr_ilist(fp,indent,interaction_function[j].longname,
	       idef,&idef->il[j]);
  }
}

static int pr_block_title(FILE *fp,int indent,char *title,t_block *block)
{
  int i;

  if (available(fp,block,title))
    {
      indent=pr_title(fp,indent,title);
      (void) pr_indent(fp,indent);
      fprintf(fp,"multinr:");
      for (i=0; i<MAXPROC; i++) (void) fprintf(fp," %d",block->multinr[i]);
      fprintf(fp,"\n");
      (void) pr_indent(fp,indent);
      (void) fprintf(fp,"nr=%d\n",block->nr);
      (void) pr_indent(fp,indent);
      (void) fprintf(fp,"nra=%d\n",block->nra);
    }
  return indent;
}

static void low_pr_block(FILE *fp,int indent,char *title,t_block *block)
{
  int i;
  
  if (available(fp,block,title))
    {
      indent=pr_block_title(fp,indent,title,block);
      for (i=0; i<=block->nr; i++)
        {
          (void) pr_indent(fp,indent+INDENT);
          (void) fprintf(fp,"%s->index[%d]=%d\n",title,i,block->index[i]);
        }
      for (i=0; i<block->nra; i++)
        {
          (void) pr_indent(fp,indent+INDENT);
          (void) fprintf(fp,"%s->a[%d]=%d\n",title,i,block->a[i]);
        }
    }
}

void pr_block(FILE *fp,int indent,char *title,t_block *block)
{
  int i,j,ok,size,start,end;
  
  if (available(fp,block,title))
    {
      indent=pr_block_title(fp,indent,title,block);
      start=0;
      end=start;
      if ((ok=(block->index[start]==0))==0)
        (void) fprintf(fp,"block->index[%d] should be 0\n",start);
      else
        for (i=0; i<block->nr; i++)
          {
            end=block->index[i+1];
            size=pr_indent(fp,indent);
            if (end<=start)
              size+=fprintf(fp,"%s[%d]={",title,i);
            else
              size+=fprintf(fp,"%s[%d][%d..%d]={",title,i,start,end-1);
            for (j=start; j<end; j++)
              {
                if (j>start) size+=fprintf(fp,", ");
                if ((size)>(USE_WIDTH))
                  {
                    (void) fprintf(fp,"\n");
                    size=pr_indent(fp,indent+INDENT);
                  }
                size+=fprintf(fp,"%d",block->a[j]);
              }
            (void) fprintf(fp,"}\n");
            start=end;
          }
      if ((end!=block->nra)||(!ok)) 
        {
          (void) pr_indent(fp,indent);
          (void) fprintf(fp,"tables inconsistent, dumping complete tables:\n");
          low_pr_block(fp,indent,title,block);
        }
    }
}

static void pr_blocks(FILE *fp,int indent,char *title,
                      t_block block[],int n,char *block_names[])
{
  int i;
  char s[STRLEN];
  
  if (available(fp,block,title))
    {
      indent=pr_title_n(fp,indent,title,n);
      for (i=0; i<n; i++)
        {
          if (block_names==NULL)
            sprintf(s,"%s[%d]",title,i);
          else
            sprintf(s,"%s[%s]",title,block_names[i]);
          pr_block(fp,indent,s,&(block[i]));
        }
    }
}

static void pr_atom(FILE *fp,int indent,char *title,t_atom *atom,int n)
{
  int i,j;
  
  if (available(fp,atom,title)) {  
    indent=pr_title_n(fp,indent,title,n);
    for (i=0; i<n; i++) {
      (void) pr_indent(fp,indent);
      fprintf(fp,"%s[%6d]={type=%3d, typeB=%3d, ptype=%8s, m=%12.5e, "
	      "q=%12.5e, mB=%12.5e, qB=%12.5e, resnr=%5d} grpnrs=[",
	      title,i,atom[i].type,atom[i].typeB,ptype_str[atom[i].ptype],
	      atom[i].m,atom[i].q,atom[i].mB,atom[i].qB,atom[i].resnr);
      for(j=0; (j<egcNR); j++)
	fprintf(fp," %d",(int)atom[i].grpnr[j]);
      fprintf(fp," ]}\n");
    }
  }
}

static void pr_grps(FILE *fp,int indent,char *title,t_grps grps[],int ngrp,
		    char **grpname[])
{
  int i,j;
  
  for(i=0; (i<ngrp); i++) {
    fprintf(fp,"%s[%d] nr=%d, name=[",title,i,grps[i].nr);
    for(j=0; (j<grps[i].nr); j++)
      fprintf(fp," %s",*(grpname[grps[i].nm_ind[j]]));
    fprintf(fp,"]\n");
  }
}

static void pr_strings(FILE *fp,int indent,char *title,char ***nm,int n)
{
  int i;

  if (available(fp,nm,title))
    {  
      indent=pr_title_n(fp,indent,title,n);
      for (i=0; i<n; i++)
        {
          (void) pr_indent(fp,indent);
          (void) fprintf(fp,"%s[%d]={name=\"%s\"}\n",title,i,*(nm[i]));
        }
    }
}

static void pr_atoms(FILE *fp,int indent,char *title,t_atoms *atoms)
{
  if (available(fp,atoms,title))
    {
      indent=pr_title(fp,indent,title);
      pr_atom(fp,indent,"atom",atoms->atom,atoms->nr);
      pr_grps(fp,indent,"grp",atoms->grps,egcNR,atoms->grpname);
      pr_strings(fp,indent,"atom",atoms->atomname,atoms->nr);
      pr_strings(fp,indent,"residue",atoms->resname,atoms->nres);
      pr_strings(fp,indent,"grpname",atoms->grpname,atoms->ngrpname);
      pr_block(fp,indent,"excl",&atoms->excl);
    }
}

void pr_top(FILE *fp,int indent,char *title,t_topology *top)
{
  if (available(fp,top,title)) {
    indent=pr_title(fp,indent,title);
    (void) pr_indent(fp,indent);
    (void) fprintf(fp,"name=\"%s\"\n",*(top->name));
    pr_atoms(fp,indent,"atoms",&(top->atoms));
    pr_blocks(fp,indent,"blocks",top->blocks,ebNR,eblock_names);
    pr_idef(fp,indent,"idef",&top->idef);
  }
}

