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
 * GRowing Old MAkes el Chrono Sweat
 */

#ifdef HAVE_IDENT
#ident "@(#) toputil.c 1.68 9/30/97"
#endif

#include "assert.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "macros.h"
#include "string2.h"
#include "topdirs.h"
#include "toputil.h"
#include "topdirs.h"
#include "toputil.h"
#include "symtab.h"
#include "fatal.h"

/* UTILITIES */

int at2type(char *str,t_atomtype *at)
{
  int i;

  for (i=0; (i<at->nr) && strcasecmp(str,*(at->atomname[i])); i++)
    ;
  if (i == at->nr) 
    fatal_error(0,"Atomtype '%s' not found!",str);
  
  return i;
}

char *type2nm(int nt, t_atomtype *at)
{
  if ((nt < 0) || (nt >= at->nr))
    fatal_error(0,"nt out of range: %d",nt);
  
  return *(at->atomname[nt]);
}

void pr_alloc (int extra, t_params *pr)
{
  /* get new space for arrays */
  if (pr->nr+extra < 0) {
    fprintf(stderr,"Trying to make array < 0 bytes\n");
    exit(1);
  }
  if (extra) {
    if (pr->nr+extra == 0) {
      sfree(pr->param);
      pr->param=NULL;
    }
    else {
      int i,j;
      srenew(pr->param,pr->nr+extra);
      for(i=pr->nr; (i<pr->nr+extra); i++) {
	for(j=0; (j<MAXATOMLIST); j++)
	  pr->param[i].a[j]=0;
	for(j=0; (j<MAXFORCEPARAM); j++)
	  pr->param[i].c[j]=0;
      }
    }
  }
  if (pr->nr+extra < 0) {
    fprintf(stderr,"Trying to make array < 0 bytes\n");
    exit(1);
  }
  if (extra) {
    if (pr->nr+extra == 0) {
      sfree(pr->param);
      pr->param=NULL;
    }
    else {
      int i,j;
      srenew(pr->param,pr->nr+extra);
      for(i=pr->nr; (i<pr->nr+extra); i++) {
	for(j=0; (j<MAXATOMLIST); j++)
	  pr->param[i].a[j]=0;
	for(j=0; (j<MAXFORCEPARAM); j++)
	  pr->param[i].c[j]=0;
      }
    }
  }
}

/* INIT STRUCTURES */

void init_block(t_block *block)
{
  int i;

  block->nr    = 0;
  block->nra   = 0;
  snew(block->index,1);
  block->index[0] = 0;
  block->a     = NULL;
  for(i=0; (i<MAXPROC); i++)
    block->multinr[i]=0;
}

void init_atom (t_atoms *at)
{
  int i;
  
  init_block(&(at->excl));
  at->nr       = 0;
  at->nres     = 0;
  at->ngrpname = 0;
  at->atom     = NULL;
  at->resname  = NULL;
  at->atomname = NULL;
  at->grpname  = NULL;
  for(i=0; (i<egcNR); i++) {
    at->grps[i].nr=0;
    at->grps[i].nm_ind=NULL;
  }
}

void init_atomtype (t_atomtype *at)
{
  at->nr   = 0;
  at->atom     = NULL;
  at->atomname = NULL;
  at->nb       = NULL;
}

void init_plist(t_params plist[])
{
  int i;
  
  for(i=0; (i<F_NRE); i++) {
    plist[i].nr    = 0;
    plist[i].param = NULL;
  }
}

void init_molinfo(t_molinfo *mol)
{
  mol->nrexcl = 0;
  mol->excl_set = FALSE;
  mol->bProcessed = FALSE;
  init_plist(mol->plist);
  init_block(&mol->cgs);
  init_block(&mol->mols);
  init_atom(&mol->atoms);
}

void init_top (t_topology *top)
{
  int i;
  
  top->name = NULL;
  init_atom (&(top->atoms));
  for (i=0; (i<ebNR); i++)
    init_block(&(top->blocks[i]));
}

/* FREEING MEMORY */

void done_bt (t_params *pl)
{
  sfree(pl->param);
}

void done_block(t_block *block)
{
  block->nr    = 0;
  block->nra   = 0;
  sfree(block->index);
  sfree(block->a);
}

void done_atom (t_atoms *at)
{
  done_block(&(at->excl));
  at->nr       = 0;
  at->nres     = 0;
  sfree(at->atom);
  sfree(at->resname);
  sfree(at->atomname);
}

void done_top(t_topology *top)
{
  int i;
  
  done_atom (&(top->atoms));
  for (i=0; (i<ebNR); i++)
    done_block(&(top->blocks[i]));
}

void done_mi(t_molinfo *mi)
{
  int i;
  
  done_atom (&(mi->atoms));
  done_block(&(mi->cgs));
  done_block(&(mi->mols));
  for(i=0; (i<F_NRE); i++)
    done_bt(&(mi->plist[i]));
}

/* PRINTING STRUCTURES */

void print_at (FILE * out, t_atomtype *at)
{
  int i;
  t_atom     *atom = at->atom;
  t_param    *nb   = at->nb;
  
  fprintf (out,"[ %s ]\n",dir2str(d_atomtypes));
  fprintf (out,"; %6s  %8s  %8s  %8s  %12s  %12s\n",
	   "type","mass","charge","particle","c6","c12");
  for (i=0; (i<at->nr); i++) 
    fprintf(out,"%8s  %8.3f  %8.3f  %8s  %12e  %12e\n",
	    *(at->atomname[i]),atom[i].m,atom[i].q,"A",
	    nb[i].C0,nb[i].C1);
  
  fprintf (out,"\n");
}

static void print_nbt (FILE *out,char *title,t_atomtype *at,
		       int ftype,t_params *nbt)
{
  int f,i,j,k,l,nrfp;
  
  if (ftype == F_LJ)
    f=1;
  else
    f=2;
  nrfp=NRFP(ftype);
  
  if (nbt->nr) {
    /* header */
    fprintf (out,"; %s\n",title);
    fprintf (out,"[ %s ]\n",dir2str(d_nonbond_params));
    fprintf (out,"; %6s  %8s","ai","aj");
    fprintf (out,"  %8s","funct");
    for (j=0; (j<nrfp); j++)
      fprintf (out,"  %11c%1d",'c',j);
    fprintf (out,"\n");
    
    /* print non-bondtypes */
    for (i=k=0; (i<at->nr); i++)
      for(j=0; (j<at->nr); j++,k++) {
	fprintf (out,"%8s  %8s  %8d",
		 *(at->atomname[i]),*(at->atomname[j]),f);
	for(l=0; (l<nrfp); l++)
	  fprintf (out,"  %12.5e",nbt->param[k].c[l]);
	fprintf (out,"\n");
      }
  }
  fprintf (out,"\n");
}

void print_bt(FILE *out, directive d, t_atomtype *at,
	      int ftype,t_params plist[],bool bConsts,
	      bool bFullDih)
{
  /* This dihp is a DIRTY patch because the dih-types do not use
   * all four atoms to determine the type.
   */
  static int dihp[2][2] = { { 1,2 }, { 0,3 } };
  t_params *bt;
  int      i,j,f,nral,nrfp;
  bool     bDih;
  
  bt=&(plist[ftype]);
  
  if (!bt->nr) 
    return;
  
  f = 0;
  switch(ftype) {
  case F_MORSE:
    f = 1;
    break;
  case F_PDIHS:
  case F_RBDIHS:
    bDih=TRUE;
    break;
  case F_IDIHS:
    f=1;
    bDih=TRUE;
    break;
  default:
    bDih=FALSE;
  }
  if (bFullDih)
    bDih=FALSE;
    
  nral = NRAL(ftype);
  nrfp = NRFP(ftype);
  
  /* header */
  fprintf(out,"[ %s ]\n",dir2str(d));
  fprintf(out,"; ");
  if (!bDih) {
    fprintf (out,"%3s  %4s","ai","aj");
    for (j=2; (j<nral); j++)
      fprintf (out,"  %3c%c",'a','i'+j);
  }
  else 
    for (j=0; (j<2); j++)
      fprintf (out,"%3c%c",'a','i'+dihp[f][j]);
  
  fprintf (out," funct");
  for (j=0; (j<nrfp); j++)
    fprintf (out," %11c%1d",'c',j);
  fprintf (out,"\n");
  
  /* print bondtypes */
  for (i=0; (i<bt->nr); i++) {
    if (!bDih)
      for (j=0; (j<nral); j++)
	fprintf (out,"%5s ",*(at->atomname[bt->param[i].a[j]]));
    else 
      for(j=0; (j<2); j++)
	fprintf (out,"%5s ",
		 *(at->atomname[bt->param[i].a[dihp[f][j]]]));
    fprintf (out,"%5d ",f+1);
    for (j=0; (j<nrfp && (bConsts || (bt->param[i].c[j] != NOTSET))); j++)
      fprintf (out,"%12e ",bt->param[i].c[j]);
    
    fprintf (out,"\n");
  }
  fprintf (out,"\n");
  fflush (out);
}

void print_block (FILE *out, char *szName, 
		  char *szIndex, char *szA,
		  t_block *block)
{
  int i,j;
  
  fprintf (out,"; %s\n",szName);
  fprintf (out,"; %4s    %s\n",szIndex,szA);
  for (i=0; (i < block->nr); i++) {
    for (i=0; (i < block->nr); i++) {
      fprintf (out,"%6d",i+1);
      for (j=block->index[i]; (j < ((int)block->index[i+1])); j++)
	fprintf (out,"%5d",block->a[j]+1);
      fprintf (out,"\n");
    }
    fprintf (out,"\n");
  }
}

void print_excl(FILE *out, t_block *excl)
{
  int     i;
  atom_id j;
  
  fprintf (out,"[ %s ]\n",dir2str(d_exclusions));
  fprintf (out,"; %4s    %s\n","i","excluded from i");
  for (i=0; (i < excl->nr); i++) {
    fprintf (out,"%6d",i+1);
    for (j=excl->index[i]; (j < excl->index[i+1]); j++)
      fprintf (out,"%5d",excl->a[j]+1);
    fprintf (out,"\n");
  }
  fprintf (out,"\n");
  fflush(out);
}

void print_atoms(FILE *out,t_atomtype *atype,t_atoms *at,t_block *cgs)
{
  int  i,j,k;
  int  itype,nrdum;
  char *as;
  bool bFound;
  real qtot;
  
  as=dir2str(d_atoms);
  fprintf(out,"[ %s ]\n",as);
  fprintf(out,"; %4s  %6s  %6s  %6s  %6s  %6s  %12s  %12s  %6s  %12s  %12s\n",
	  "nr","type","resnr","residu","atom","cgnr","charge","mass","typeB","chargeB","massB");
    
  /* A little shortcut (at, nrdum) */
  nrdum = at->nr;
  qtot  = 0;
  
  if (debug)
    fprintf(debug,"This molecule has %d atoms and %d residues\n",
	    nrdum,at->nres);
  
  if (at->nres) {
    /* if the information is present... */
    for (i=0; (i < nrdum); i++) {
      /* quadratic algorithm to find the cg-nr */
      bFound=FALSE;
      for (j=0; ((j<cgs->nr) && !bFound); j++)
	for(k=cgs->index[j]; ((k<((int)cgs->index[j+1])) && !bFound); k++)
	  if (cgs->a[k]==i) {
	    bFound=TRUE;
	    break;
	  }
      if (!bFound)
	fatal_error(0,"Atom %d not in any charge group!",i);
	
      itype=at->atom[i].type;
      if ((itype < 0) || (itype > atype->nr))
	fatal_error(0,"itype = %d, i= %d in print_atoms",itype,i);
	
      fprintf(out,"%6d  %6s  %6d  %6s  %6s  %6d  %12g  %12g",
	      i+1,*(atype->atomname[itype]),
	      at->atom[i].resnr+1,  
	      *(at->resname[at->atom[i].resnr]),
	      *(at->atomname[i]),j,
	      at->atom[i].q,at->atom[i].m);
      if (PERTURBED(at->atom[i])) {
	fprintf(out,"  %6s  %12g  %12g",
		*(atype->atomname[at->atom[i].typeB]),
		at->atom[i].qB,at->atom[i].mB);
      }
      qtot+=at->atom[i].q;
      fprintf(out,"\t; qtot: %g\n",qtot);
    }
  }
  fprintf(out,"\n");
  fflush(out);
}

void print_bonds(FILE *out,int natoms,directive d,
		 int ftype,t_params plist[],bool bConsts)
{
  t_symtab   dumtab;
  t_atomtype dumatype;
  int i;

  snew(dumatype.atom,natoms);
  snew(dumatype.atomname,natoms);
  open_symtab(&dumtab);
  for (i=0; (i < natoms); i++) {
    char buf[12];
    sprintf(buf,"%4d",(i+1));
    dumatype.atomname[i]=put_symtab(&dumtab,buf);
  }
  print_bt(out,d,&dumatype,ftype,plist,bConsts,TRUE);
    
  rm_symtab(&dumtab);
  sfree(dumatype.atom);
  sfree(dumatype.atomname);
}

