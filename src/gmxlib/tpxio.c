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
 * Good ROcking Metal Altar for Chronical Sinners
 */
 
#include <ctype.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "string2.h"
#include "fatal.h"
#include "macros.h"
#include "assert.h"
#include "names.h"
#include "symtab.h"
#include "futil.h"
#include "filenm.h"
#include "gmxfio.h"
#include "tpxio.h"
#include "confio.h"
#include "mass.h"
#include "copyrite.h"

/* This number should be increased whenever the file format changes! */
static int tpx_version  = 2;
/* This is the version of the file we are reading */
static int file_version = 0;

void _do_section(int fp,int key,bool bRead,char *src,int line)
{
  char buf[STRLEN];
  bool bDbg;

  if (fio_getftp(fp) == efTPA) {
    if (!bRead) {
      do_string(itemstr[key]);
      bDbg       = fio_getdebug(fp);
      fio_setdebug(fp,FALSE);
      do_string(comment_str[key]);
      fio_setdebug(fp,bDbg);
    }
    else {
      if (fio_getdebug(fp))
	fprintf(stderr,"Looking for section %s (%s, %d)",
		itemstr[key],src,line);
      
      do {
	do_string(buf);
      } while ((strcasecmp(buf,itemstr[key]) != 0));
      
      if (strcasecmp(buf,itemstr[key]) != 0) 
	fatal_error(0,"\nCould not find section heading %s",itemstr[key]);
      else if (fio_getdebug(fp))
	fprintf(stderr," and found it\n");
    }
  }
}

#define do_section(key,bRead) _do_section(fp,key,bRead,__FILE__,__LINE__)

/**************************************************************
 *
 * Now the higer level routines that do io of the structures and arrays
 *
 **************************************************************/
static void do_inputrec(t_inputrec *ir,bool bRead)
{
  int i,j; 
  bool bDum=TRUE;

  if (file_version >= 1) {  
    /* Basic inputrec stuff */  
    do_int(ir->eI); 
    do_int(ir->nsteps); 
    do_int(ir->ns_type); 
    do_int(ir->nstlist); 
    do_int(ir->ndelta); 
    do_int(ir->nstcomm); 
    do_int(ir->nstcgsteep); 
    do_int(ir->nstprint); 
    do_int(ir->nstxout); 
    do_int(ir->nstvout); 
    do_int(ir->nstfout); 
    do_int(ir->nstgrp); 
    do_int(ir->nstxtcout); 
    do_real(ir->init_t); 
    do_real(ir->delta_t); 
    do_real(ir->xtcprec); 
    do_int(ir->niter); 
    do_int(ir->solvent_opt); 
    do_int(ir->nsatoms); 
    do_real(ir->gausswidth); 
    do_int(ir->nkx); 
    do_int(ir->nky); 
    do_int(ir->nkz); 
    do_int(ir->eBox); 
    do_int(ir->bUncStart); 
    do_int(ir->btc); 
    do_int(ir->ntcmemory); 
    do_int(ir->epc); 
    do_int(ir->npcmemory); 
    do_real(ir->tau_p); 
    do_rvec(ir->ref_p); 
    do_rvec(ir->compress); 
    do_int(ir->bSimAnn); 
    do_real(ir->zero_temp_time); 
    do_int(ir->eeltype); 
    do_real(ir->rshort); 
    do_real(ir->rlong); 
    do_int(ir->bLJcorr); 
    do_real(ir->epsilon_r); 
    do_real(ir->shake_tol); 
    do_real(ir->fudgeLJ); 
    do_real(ir->fudgeQQ); 
    do_int(ir->bPert); 
    do_real(ir->init_lambda); 
    do_real(ir->delta_lambda); 
    do_real(ir->dr_fc); 
    do_real(ir->dr_tau); 
    do_real(ir->dihr_fc); 
    do_real(ir->em_stepsize); 
    do_real(ir->em_tol); 
    do_int(ir->eConstrAlg); 
    do_int(ir->nProjOrder); 
    do_real(ir->LincsWarnAngle); 
    do_int(ir->nstLincsout); 
    do_real(ir->ld_temp); 
    do_real(ir->ld_fric); 
    do_int(ir->ld_seed); 
    do_int(ir->userint1); 
    do_int(ir->userint2); 
    do_int(ir->userint3); 
    do_int(ir->userint4); 
    do_real(ir->userreal1); 
    do_real(ir->userreal2); 
    do_real(ir->userreal3); 
    do_real(ir->userreal4); 
    
    /* grpopts stuff */
    do_int(ir->opts.ngtc); 
    do_int(ir->opts.ngacc); 
    do_int(ir->opts.ngfrz); 
    do_int(ir->opts.ngener); 
    if (bRead) {
      snew(ir->opts.nrdf,   ir->opts.ngtc); 
      snew(ir->opts.ref_t,  ir->opts.ngtc); 
      snew(ir->opts.tau_t,  ir->opts.ngtc); 
      snew(ir->opts.nFreeze,ir->opts.ngfrz); 
      snew(ir->opts.acc,    ir->opts.ngacc); 
    } 
    if (ir->opts.ngtc > 0) {
      ndo_int (ir->opts.nrdf, ir->opts.ngtc,bDum); 
      ndo_real(ir->opts.ref_t,ir->opts.ngtc,bDum); 
      ndo_real(ir->opts.tau_t,ir->opts.ngtc,bDum); 
    }
    if (ir->opts.ngfrz > 0) 
      ndo_ivec(ir->opts.nFreeze,ir->opts.ngfrz,bDum);
    if (ir->opts.ngacc > 0) 
      ndo_rvec(ir->opts.acc,ir->opts.ngacc); 
    
    /* Cosine stuff for electric fields */
    for(j=0; (j<DIM); j++) {
      do_int  (ir->ex[j].n);
      do_int  (ir->et[j].n);
      if (bRead) {
	snew(ir->ex[j].a,  ir->ex[j].n);
	snew(ir->ex[j].phi,ir->ex[j].n);
	snew(ir->et[j].a,  ir->et[j].n);
	snew(ir->et[j].phi,ir->et[j].n);
      }
      ndo_real(ir->ex[j].a,  ir->ex[j].n,bDum);
      ndo_real(ir->ex[j].phi,ir->ex[j].n,bDum);
      ndo_real(ir->et[j].a,  ir->et[j].n,bDum);
      ndo_real(ir->et[j].phi,ir->et[j].n,bDum);
    }
  }  
  /* New version of the inputrec will work like this */
  if (file_version >= 2) {
    /* Do version 2 specific things */
  }  
  if (file_version != tpx_version) {
    /* Give a warning about features that are not accessible */
    fprintf(stderr,"WARNING: tpx file_version %d, software version %d\n",
	    tpx_version,file_version);
  }
}

static void do_harm(t_iparams *iparams,bool bRead)
{
  do_real(iparams->harmonic.rA);
  do_real(iparams->harmonic.krA);
  do_real(iparams->harmonic.rB);
  do_real(iparams->harmonic.krB);
}

void do_iparams(t_functype ftype,t_iparams *iparams,bool bRead)
{
  int i;
  bool bDum;
  
  if (!bRead)
    set_comment(interaction_function[ftype].name);
  switch (ftype) {
  case F_ANGLES:
    do_harm(iparams,bRead);
    break;
  case F_BHAM:
    do_real(iparams->bham.a);
    do_real(iparams->bham.b);
    do_real(iparams->bham.c);
    break;
  case F_BONDS:
    do_harm(iparams,bRead);
    break;
  case F_IDIHS:
    do_harm(iparams,bRead);
    break;
  case F_MORSE:
    do_real(iparams->morse.b0);
    do_real(iparams->morse.cb);
    do_real(iparams->morse.beta);
    break;
  case F_LJ:
    do_real(iparams->lj.c6);
    do_real(iparams->lj.c12);
    break;
  case F_LJ14:
    do_real(iparams->lj14.c6A);
    do_real(iparams->lj14.c12A);
    do_real(iparams->lj14.c6B);
    do_real(iparams->lj14.c12B);
    break;
  case F_PDIHS:
    do_real(iparams->pdihs.phiA);
    do_real(iparams->pdihs.cpA);
    do_real(iparams->pdihs.phiB);
    do_real(iparams->pdihs.cpB);
    do_int (iparams->pdihs.mult);
    break;
  case F_DISRES:
    do_int (iparams->disres.index);
    do_int (iparams->disres.type);
    do_real(iparams->disres.rx0);
    do_real(iparams->disres.rx1);
    do_real(iparams->disres.rx2);
    do_real(iparams->disres.rx3);
    break;
  case F_POSRES:
    do_rvec(iparams->posres.pos0);
    do_rvec(iparams->posres.fc);
    break;
  case F_RBDIHS:
    ndo_real(iparams->rbdihs.rbc,NR_RBDIHS,bDum);
    break;
  case F_SHAKE:
    do_real(iparams->shake.dA);
    do_real(iparams->shake.dB);
    break;
  case F_SETTLE:
    do_real(iparams->settle.doh);
    do_real(iparams->settle.dhh);
    break;
  case F_DUMMY1:
    do_real(iparams->dummy.a);
    break;
  case F_DUMMY2:
  case F_DUMMY2FD:
  case F_DUMMY2FAD:
    do_real(iparams->dummy.a);
    do_real(iparams->dummy.b);
    break;
  case F_DUMMY3:
    do_real(iparams->dummy.a);
    do_real(iparams->dummy.b);
    do_real(iparams->dummy.c);
    break;
  default:
    fatal_error(0,"unknown function type %d (%s) in %s line %d",
		ftype,interaction_function[ftype].name,__FILE__,__LINE__);
  }
  if (!bRead)
    unset_comment();
}

static void do_ilist(t_ilist *ilist,bool bRead,char *name)
{
  int i;
  bool bDum=TRUE;
  
  if (!bRead)
    set_comment(name);
  ndo_int(ilist->multinr,MAXPROC,bDum);
  do_int (ilist->nr);
  if (bRead)
    snew(ilist->iatoms,ilist->nr);
  ndo_int(ilist->iatoms,ilist->nr,bDum);
  if (!bRead)
    unset_comment();
}

static void do_idef(t_idef *idef,bool bRead)
{
  int i,j;
  bool bDum=TRUE;
  
  do_int(idef->atnr);
  do_int(idef->pid);
  do_int(idef->ntypes);
  if (bRead) {
    snew(idef->functype,idef->ntypes);
    snew(idef->iparams,idef->ntypes);
  }
  ndo_int(idef->functype,idef->ntypes,bDum);
  
  for (i=0; (i<idef->ntypes); i++) 
    do_iparams(idef->functype[i],&idef->iparams[i],bRead);
  
  for(j=0; (j<F_NRE); j++)
    do_ilist(&idef->il[j],bRead,interaction_function[j].name);
}

static void do_block(t_block *block,bool bRead)
{
  int i;
  bool bDum=TRUE;

  ndo_int(block->multinr,MAXPROC,bDum);
  do_int (block->nr);
  do_int (block->nra);
  if (bRead) {
    snew(block->index,block->nr+1);
    snew(block->a,block->nra);
  }
  ndo_int(block->index,block->nr+1,bDum);
  ndo_int(block->a,block->nra,bDum);
}

static void do_atom(t_atom *atom,bool bRead)
{
  int ngrp=egcNR;
  
  do_real (atom->m);
  do_real (atom->q);
  do_real (atom->mB);
  do_real (atom->qB);
  do_ushort(atom->type);
  do_ushort(atom->typeB);
  do_int (atom->ptype);
  do_int (atom->resnr);
  do_nuchar(atom->grpnr,ngrp);
}

static void do_grps(int ngrp,t_grps grps[],bool bRead)
{
  int i,j;
  bool bDum=TRUE;
  
  for(j=0; (j<ngrp); j++) {
    do_int (grps[j].nr);
    if (bRead)
      snew(grps[j].nm_ind,grps[j].nr);
    ndo_int(grps[j].nm_ind,grps[j].nr,bDum);
  }
}

static void do_symstr(char ***nm,bool bRead,t_symtab *symtab)
{
  int ls;
  
  if (bRead) {
    do_int(ls);
    *nm = get_symtab_handle(symtab,ls);
  }
  else {
    ls = lookup_symtab(symtab,*nm);
    do_int(ls);
  }
}

static void do_strstr(int nstr,char ***nm,bool bRead,t_symtab *symtab)
{
  int  j;
  
  for (j=0; (j<nstr); j++) 
    do_symstr(&(nm[j]),bRead,symtab);
}

static void do_atoms(t_atoms *atoms,bool bRead,t_symtab *symtab)
{
  int i;
  
  do_int(atoms->nr);
  do_int(atoms->nres);
  do_int(atoms->ngrpname);
  if (bRead) {
    snew(atoms->atom,atoms->nr);
    snew(atoms->atomname,atoms->nr);
    snew(atoms->resname,atoms->nres);
    snew(atoms->grpname,atoms->ngrpname);
  }
  for(i=0; (i<atoms->nr); i++)
    do_atom(&atoms->atom[i],bRead);
  do_strstr(atoms->nr,atoms->atomname,bRead,symtab);
  do_strstr(atoms->nres,atoms->resname,bRead,symtab);
  do_strstr(atoms->ngrpname,atoms->grpname,bRead,symtab);
  
  do_grps(egcNR,atoms->grps,bRead);
  
  do_block(&(atoms->excl),bRead);
}

static void do_symtab(t_symtab *symtab,bool bRead)
{
  int i,nr;
  t_symbuf *symbuf;
  char buf[STRLEN];
  
  do_int(symtab->nr);
  nr     = symtab->nr;
  if (bRead) {
    snew(symtab->symbuf,1);
    symbuf = symtab->symbuf;
    symbuf->bufsize = nr;
    snew(symbuf->buf,nr);
    for (i=0; (i<nr); i++) {
      do_string(buf);
      symbuf->buf[i]=strdup(buf);
    }
  }
  else {
    symbuf = symtab->symbuf;
    while (symbuf!=NULL) {
      for (i=0; (i<symbuf->bufsize) && (i<nr); i++) 
	do_string(symbuf->buf[i]);
      nr-=i;
      symbuf=symbuf->next;
    }
    assert(nr==0);
  }
}

static void make_chain_identifiers(t_atoms *atoms,t_block *mols)
{
  int m,a,a0,a1;
  char c,chain;
#define CHAIN_MIN_ATOMS 15

  chain='A';
  for(m=0; m<mols->nr; m++) {
    a0=mols->index[m];
    a1=mols->index[m+1];
    if ((a1-a0 >= CHAIN_MIN_ATOMS) && (chain <= 'Z')) {
      c=chain;
      chain++;
    } else
      c=' ';
    for(a=a0; a<a1; a++)
      atoms->atom[mols->a[a]].chain=c;  
  }
  if (chain == 'B')
    for(a=0; a<atoms->nr; a++)
      atoms->atom[a].chain=' ';
}
  
static void do_top(t_topology *top,bool bRead)
{
  int  i;
  
  do_symtab(&(top->symtab),bRead);
  do_symstr(&(top->name),bRead,&(top->symtab));
  do_atoms (&(top->atoms),bRead,&(top->symtab));
  do_idef  (&(top->idef),bRead);
  for(i=0; (i<ebNR); i++)
    do_block(&(top->blocks[i]),bRead);
  if (bRead) {
    close_symtab(&(top->symtab));
    make_chain_identifiers(&(top->atoms),&(top->blocks[ebMOLS]));
  }
}

static void do_tpxheader(int fp,bool bRead,t_tpxheader *tpx)
{
  char  buf[STRLEN];
  bool  bDouble;
  int   i,precision;
 
  fio_select(fp);
  fio_setdebug(fp,bDebugMode());
  
  /* NEW! XDR tpb file */
  precision = sizeof(real);
  if (bRead) {
    do_string(buf);
    do_int(precision);
    bDouble = (precision == sizeof(double));
    if ((precision != sizeof(float)) && !bDouble)
      fatal_error(0,"Unknown precision in file %s: real is %d bytes "
		  "instead of %d or %d",
		  fio_getname(fp),precision,sizeof(float),sizeof(double));
    fio_setprecision(fp,bDouble);
    fprintf(stderr,"Reading file %s, %s (%s precision)\n",
	    fio_getname(fp),buf,bDouble ? "double" : "single");
  }
  else {
    do_string(GromacsVersion());
    bDouble = (precision == sizeof(double));
    fio_setprecision(fp,bDouble);
    do_int(precision);
    file_version = tpx_version;
  }
  
  /* Check versions! */
  do_int(file_version);
  if (file_version == 1) 
    fatal_error(0,"Reading tpx file (%s) version %d with version %d program.\n",
		fio_getname(fp),file_version,tpx_version);
  else if (file_version != tpx_version) 
    fprintf(stderr,"WARNING: reading tpx file (%s) version %d with version %d"
	    " program. Some options may not work.\n",
	    fio_getname(fp),file_version,tpx_version);
    
  do_section(eitemHEADER,bRead);
  do_int (tpx->natoms);
  do_int (tpx->step);
  do_real(tpx->t);
  do_real(tpx->lambda);
  do_int (tpx->bIr);
  do_int (tpx->bTop);
  do_int (tpx->bX);
  do_int (tpx->bV);
  do_int (tpx->bF);
  do_int (tpx->bBox);
}

static void do_tpx(int fp,bool bRead,int *step,real *t,real *lambda,
		   t_inputrec *ir,rvec *box,int *natoms,
		   rvec *x,rvec *v,rvec *f,t_topology *top)
{
  t_tpxheader tpx;
  t_inputrec  dum_ir;
  t_topology  dum_top;
  rvec        dum_x;
   
  if (!bRead) {
    tpx.natoms = *natoms;
    tpx.step   = *step;
    tpx.t      = *t;
    tpx.lambda = *lambda;
    tpx.bIr  = (ir  != NULL);
    tpx.bTop = (top != NULL);
    tpx.bX   = (x   != NULL);
    tpx.bV   = (v   != NULL);
    tpx.bF   = (f   != NULL);
    tpx.bBox = (box != NULL);
  }

  do_tpxheader(fp,bRead,&tpx);

  if (bRead) {
    *natoms = tpx.natoms;
    *step   = tpx.step;
    *t      = tpx.t;
    *lambda = tpx.lambda;
  }
    
#define do_test(b,p) if (bRead && (p!=NULL) && !b) fatal_error(0,"No %s in %s",#p,fio_getname(fp)) 

  do_test(tpx.bBox,box);
  do_section(eitemBOX,bRead);
  if (tpx.bBox) ndo_rvec(box,DIM);

  do_test(tpx.bIr,ir);
  do_section(eitemIR,bRead);
  if (tpx.bIr) {
    if (ir)
      do_inputrec(ir,bRead);
    else {
      init_inputrec(&dum_ir);
      do_inputrec  (&dum_ir,bRead);
      done_inputrec(&dum_ir);
    }
  }
  do_test(tpx.bTop,top);
  do_section(eitemTOP,bRead);
  if (tpx.bTop) {
    if (top)
      do_top(top,bRead);
    else {
      do_top(&dum_top,bRead);
      done_top(&dum_top);
    }
  }
  do_test(tpx.bX,x);  
  do_section(eitemX,bRead);
  if (tpx.bX) ndo_rvec(x,*natoms);
  
  do_test(tpx.bV,v);
  do_section(eitemV,bRead);
  if (tpx.bV) ndo_rvec(v,*natoms);

  do_test(tpx.bF,f);
  do_section(eitemF,bRead);
  if (tpx.bF) ndo_rvec(f,*natoms);
}

/************************************************************
 *
 *  The following routines are the exported ones
 *
 ************************************************************/

int open_tpx(char *fn,char *mode)
{
  return fio_open(fn,mode);
}    
 
void close_tpx(int fp)
{
  fio_close(fp);
}

void read_tpxheader(char *fn,t_tpxheader *tpx)
{
  int fp;

  fp = open_tpx(fn,"r");
  do_tpxheader(fp,TRUE,tpx);
  close_tpx(fp);
}

void write_tpx(char *fn,int step,real t,real lambda,
	       t_inputrec *ir,rvec *box,int natoms,
	       rvec *x,rvec *v,rvec *f,t_topology *top)
{
  int fp;

  fp = open_tpx(fn,"w");
  do_tpx(fp,FALSE,&step,&t,&lambda,ir,box,&natoms,x,v,f,top);
  close_tpx(fp);
}

void read_tpx(char *fn,int *step,real *t,real *lambda,
	      t_inputrec *ir,rvec *box,int *natoms,
	      rvec *x,rvec *v,rvec *f,t_topology *top)
{
  int fp;
  
  fp = open_tpx(fn,"r");
  do_tpx(fp,TRUE,step,t,lambda,ir,box,natoms,x,v,f,top);
  close_tpx(fp);
}

void fwrite_tpx(int fp,int step,real t,real lambda,
		t_inputrec *ir,rvec *box,int natoms,
		rvec *x,rvec *v,rvec *f,t_topology *top)
{
  do_tpx(fp,FALSE,&step,&t,&lambda,ir,box,&natoms,x,v,f,top);
}


void fread_tpx(int fp,int *step,real *t,real *lambda,
	       t_inputrec *ir,rvec *box,int *natoms,
	       rvec *x,rvec *v,rvec *f,t_topology *top)
{
  do_tpx(fp,TRUE,step,t,lambda,ir,box,natoms,x,v,f,top);
}

bool fn_bTPX(char *file)
{
  switch (fn2ftp(file)) {
  case efTPR:
  case efTPB:
  case efTPA:
    return TRUE;
  default:
    return FALSE;
  }
}

int read_tps_conf(char *infile,char *title,t_topology *top,t_atoms **atoms, 
		  rvec **x,rvec **v,matrix box,bool bMass)
{
  t_tpxheader  header;
  int          natoms,i;
  
  if (fn_bTPX(infile)) {
    read_tpxheader(infile,&header);
    snew(*x,header.natoms);
    if (v)
      snew(*v,header.natoms);
    read_tpx(infile,NULL,NULL,NULL,NULL,box,&natoms,
	     *x,(v==NULL) ? NULL : *v,NULL,top);
    strcpy(title,*top->name);
    *atoms = &top->atoms;
  }
  else {
    get_stx_coordnum(infile,&natoms);
    snew(*atoms,1);
    init_t_atoms(*atoms,natoms,FALSE);
    snew(*x,natoms);
    if (v)
      snew(*v,natoms);
    read_stx_conf(infile,title,*atoms,*x,(v==NULL) ? NULL : *v,box);
    if (bMass)
      for(i=0; i<natoms; i++)
	(*atoms)->atom[i].m = 
	  get_mass(*(*atoms)->resname[(*atoms)->atom[i].resnr],
		   *(*atoms)->atomname[i]);
    top = NULL;
  }

  return natoms;
}
