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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include <string.h>
#include "sysstuff.h"
#include "strdb.h"
#include "futil.h"
#include "macros.h"
#include "names.h"
#include "string2.h"
#include "statutil.h"
#include "confio.h"
#include "main.h"
#include "copyrite.h"
#include "typedefs.h"
#include "smalloc.h"
#include "invblock.h"
#include "macros.h"
#include "index.h"
#include "txtdump.h"

typedef enum { etOther, etProt, etDNA, erestNR } eRestp;
static const char *ResTP[erestNR] = { "OTHER", "PROTEIN", "DNA" };

static const char   *Sugars[]     = { "A", "T", "G", "C", "U" };
#define  NDNA asize(Sugars)

static bool yn(bool bASK)
{
  char c;

  if (bASK) {
    do {
      c=toupper(fgetc(stdin));
    } while ((c != 'Y') && (c != 'N'));

    return (c == 'Y');
  }
  else
    return FALSE;
}

t_block *new_block(void)
{
  t_block *block;

  snew(block,1);
  snew(block->index,1);

  return block;
}

void write_index(char *outf, t_block *b,char **gnames)
{
  FILE *out;
  int  i,j,k;

  out=ffopen(outf,"w");
  /* fprintf(out,"%5d  %5d\n",b->nr,b->nra); */
  for(i=0; (i<b->nr); i++) {
    fprintf(out,"[ %s ]\n",gnames[i]);
    for(k=0,j=b->index[i]; j<b->index[i+1]; j++,k++) {
      fprintf(out,"%4d ",b->a[j]+1);
      if ((k % 15) == 14)
	fprintf(out,"\n");
    }
    fprintf(out,"\n");
  }
  fclose(out);
}

void add_grp(t_block *b,char ***gnames,int nra,atom_id a[],const char *name)
{
  int i;

  srenew(b->index,b->nr+2);
  srenew(*gnames,b->nr+1);
  (*gnames)[b->nr]=strdup(name);

  srenew(b->a,b->nra+nra);
  for(i=0; (i<nra); i++)
    b->a[b->nra++]=a[i];
  b->nr++;
  b->index[b->nr]=b->nra;
}

/* compare index in `a' with group in `b' at `index', 
   when `index'<0 it is relative to end of `b' */
static bool grp_cmp(t_block *b, int nra, atom_id a[], int index)
{
  int i;
  
  if (index < 0)
    index = b->nr-1+index;
  if (index >= b->nr)
    gmx_fatal(FARGS,"no such index group %d in t_block (nr=%d)",index,b->nr);
  /* compare sizes */
  if ( nra != b->index[index+1] - b->index[index] )
    return FALSE;
  for(i=0; i<nra; i++)
    if ( a[i] != b->a[b->index[index]+i] )
      return FALSE;
  return TRUE;
}

static void p_status(int nres,eRestp restp[],bool bVerb)
{
  int i,j,ntp[erestNR];

  for(i=0; (i<erestNR); i++)
    ntp[i]=0;
  for(j=0; (j<nres); j++)
    ntp[restp[j]]++;
  
  if (bVerb)
    for(i=0; (i<erestNR); i++) 
      printf("There are: %5d %10s residues\n",ntp[i],ResTP[i]);
}

atom_id *mk_aid(t_atoms *atoms,eRestp restp[],eRestp res,int *nra,
		bool bTrue)
/* Make an array of atom_ids for all atoms with:
 * (restp[i] == res) == bTrue
 */
{
  atom_id *a;
  int     i;

  snew(a,atoms->nr);
  *nra=0;
  for(i=0; (i<atoms->nr); i++) 
    if ((restp[atoms->atom[i].resnr] == res) == bTrue)
      a[(*nra)++]=i;
  
  return a;
}

static void analyse_other(eRestp Restp[],t_atoms *atoms,
			  t_block *gb,char ***gn,bool bASK,bool bVerb)
{
  char **restp=NULL;
  char **attp=NULL;
  char *rname,*aname;
  atom_id *other_ndx,*aid,*aaid;
  int  i,j,k,l,resnr,naid,naaid,natp,nrestp=0;
  
  for(i=0; (i<atoms->nres); i++)
    if (Restp[i] == etOther)
      break;
  if (i < atoms->nres) {
    /* we have others */
    if (bVerb)
      printf("Analysing Other...\n");
    snew(other_ndx,atoms->nr);
    for(k=0; (k<atoms->nr); k++) {
      resnr=atoms->atom[k].resnr;
      rname=*atoms->resname[resnr];
      if (Restp[resnr] ==  etOther) {
	for(l=0; (l<nrestp); l++)
	  if (strcmp(restp[l],rname) == 0)
	    break;
	if (l==nrestp) {
	  srenew(restp,++nrestp);
	  restp[nrestp-1]=strdup(rname);
	}
      }
    }
    for(i=0; (i<nrestp); i++) {
      snew(aid,atoms->nr);
      naid=0;
      for(j=0; (j<atoms->nr); j++) {
	rname=*atoms->resname[atoms->atom[j].resnr];
	if (strcmp(restp[i],rname) == 0) 
	  aid[naid++] = j;
      }
      add_grp(gb,gn,naid,aid,restp[i]);
      if (bASK) {
	printf("split %s into atoms (y/n) ? ",restp[i]);
	fflush(stdout);
	if (yn(bASK)) {
	  natp=0;
	  for(k=0; (k<naid); k++) {
	    aname=*atoms->atomname[aid[k]];
	    for(l=0; (l<natp); l++)
	      if (strcmp(aname,attp[l]) == 0)
		break;
	    if (l == natp) {
	      srenew(attp,++natp);
	      attp[natp-1]=aname;
	    }
	  }
	  if (natp > 1) {
	    for(l=0; (l<natp); l++) {
	      snew(aaid,naid);
	      naaid=0;
	      for(k=0; (k<naid); k++) {
		aname=*atoms->atomname[aid[k]];
		if (strcmp(aname,attp[l])==0) 
		  aaid[naaid++]=aid[k];
	      }
	      add_grp(gb,gn,naaid,aaid,attp[l]);
	      sfree(aaid);
	    }
	  }
	  sfree(attp);
	  attp=NULL;
	}
	sfree(aid);
      }
    }
    sfree(other_ndx);
  }
}

static void analyse_prot(eRestp restp[],t_atoms *atoms,
			 t_block *gb,char ***gn,bool bASK,bool bVerb)
{
  /* atomnames to be used in constructing index groups: */
  static const char *pnoh[]    = { "H" };
  static const char *pnodum[]  = { "MN1",  "MN2",  "MCB1", "MCB2", "MCG1", "MCG2", 
			     "MCD1", "MCD2", "MCE1", "MCE2", "MNZ1", "MNZ2" };
  static const char *calpha[]  = { "CA" };
  static const char *bb[]      = { "N","CA","C" };
  static const char *mc[]      = { "N","CA","C","O","O1","O2","OXT" };
  static const char *mcb[]     = { "N","CA","CB","C","O","O1","O2","OT","OXT" };
  static const char *mch[]     = { "N","CA","C","O","O1","O2","OT","OXT",
			     "H1","H2","H3","H" };
  /* array of arrays of atomnames: */
  static const char **chains[] = { NULL,pnoh,calpha,bb,mc,mcb,mch,mch,mch,pnodum };
#define NCH asize(chains)
  /* array of sizes of arrays of atomnames: */
  const int       sizes[NCH] = { 
    0, asize(pnoh), asize(calpha), asize(bb), 
    asize(mc), asize(mcb), asize(mch), asize(mch), asize(mch), asize(pnodum)
  };
  /* descriptive names of index groups */
  const char   *ch_name[NCH] = { 
    "Protein", "Protein-H", "C-alpha", "Backbone", 
    "MainChain", "MainChain+Cb", "MainChain+H", "SideChain", "SideChain-H", 
    "Prot-Masses"
  };
  /* construct index group containing (TRUE) or excluding (FALSE)
     given atom names */
  const bool complement[NCH] = { 
    TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE
  };
  const int  wholename[NCH]  = { -1, 0,-1,-1,-1,-1,-1,-1, 11,-1 };
  /* the index in wholename gives the first item in the arrays of 
   * atomtypes that should be tested with 'strncasecmp' in stead of
   * strcasecmp, or -1 if all items should be tested with strcasecmp
   * This is comparable to using a '*' wildcard at the end of specific
   * atom names, but that is more involved to implement...
   */
  /* only add index group if it differs from the specified one, 
     specify -1 to always add group */
  const int compareto[NCH] = { -1,-1,-1,-1,-1,-1,-1,-1,-1, 0 };

  int     i,n,j;
  atom_id *aid;
  int     nra,nnpres,npres;
  bool    match;
  char    ndx_name[STRLEN],*atnm;

  if (bVerb)
    printf("Analysing Protein...\n");
  snew(aid,atoms->nr);

  /* calculate the number of protein residues */
  npres=0;
  for(i=0; (i<atoms->nres); i++)
    if (restp[i] == etProt)
      npres++;

  /* find matching or complement atoms */
  for(i=0; (i<NCH); i++) {
    nra=0;
    for(n=0; (n<atoms->nr); n++) {
      if (restp[atoms->atom[n].resnr] == etProt) {
	match=FALSE;
	for(j=0; (j<sizes[i]); j++) {
	  /* skip digits at beginning of atomname, e.g. 1H */
	  atnm=*atoms->atomname[n];
	  while (isdigit(atnm[0]))
	    atnm++;
	  if ( (wholename[i]==-1) || (j<wholename[i]) ) {
	    if (strcasecmp(chains[i][j],atnm) == 0)
	      match=TRUE;
	  } else {
	    if (strncasecmp(chains[i][j],atnm,strlen(chains[i][j])) == 0)
	      match=TRUE;
	  }
	}
	if (match != complement[i])
	  aid[nra++]=n;
      }
    }
    /* if we want to add this group always or it differs from previous 
       group, add it: */
    if ( compareto[i] == -1 || !grp_cmp(gb,nra,aid,compareto[i]-i) )
      add_grp(gb,gn,nra,aid,ch_name[i]); 
  }
  
  if (bASK) {
    for(i=0; (i<NCH); i++) {
      printf("Split %12s into %5d residues (y/n) ? ",ch_name[i],npres);
      if (yn(bASK)) {
	int resnr;
	nra = 0;
	for(n=0;((atoms->atom[n].resnr<npres) && (n<atoms->nr));) {
	  resnr = atoms->atom[n].resnr;
	  for(;((atoms->atom[n].resnr==resnr) && (n<atoms->nr));n++) {
	    match=FALSE;
	    for(j=0;(j<sizes[i]); j++) 
	      if (strcasecmp(chains[i][j],*atoms->atomname[n]) == 0)
		match=TRUE;
	    if (match != complement[i])
	      aid[nra++]=n;
	  }
	  /* copy the residuename to the tail of the groupname */
	  if (nra > 0) {
	    sprintf(ndx_name,"%s_%s%d",
		    ch_name[i],*atoms->resname[resnr],resnr+1);
	    add_grp(gb,gn,nra,aid,ndx_name);
	    nra = 0;
	  }
	}
      } 
    }
    printf("Make group with sidechain and C=O swapped (y/n) ? ");
    if (yn(bASK)) {
      /* Make swap sidechain C=O index */
      int resnr,hold;
      nra = 0;
      for(n=0;((atoms->atom[n].resnr<npres) && (n<atoms->nr));) {
	resnr = atoms->atom[n].resnr;
	hold  = -1;
	for(;((atoms->atom[n].resnr==resnr) && (n<atoms->nr));n++)
	  if (strcmp("CA",*atoms->atomname[n]) == 0) {
	    aid[nra++]=n;
	    hold=nra;
	    nra+=2;
	  } else if (strcmp("C",*atoms->atomname[n]) == 0) {
	    if (hold == -1)
	      gmx_incons("Atom naming problem");
	    aid[hold]=n;
	  } else if (strcmp("O",*atoms->atomname[n]) == 0) {
	    if (hold == -1)
	      gmx_incons("Atom naming problem");
	    aid[hold+1]=n;
	  } else if (strcmp("O1",*atoms->atomname[n]) == 0) {
	    if (hold == -1)
	      gmx_incons("Atom naming problem");
	    aid[hold+1]=n;
	  } else 
	    aid[nra++]=n;
      }
      /* copy the residuename to the tail of the groupname */
      if (nra > 0) {
	add_grp(gb,gn,nra,aid,"SwapSC-CO");
	nra = 0;
      } 
    }
  }
  sfree(aid);
}

static void analyse_dna(eRestp restp[],t_atoms *atoms,
			t_block *gb,char ***gn,bool bASK,bool bVerb)
{
  if (bVerb)
    printf("Analysing DNA... (not really)\n");
  if (debug)
    printf("eRestp %p; atoms %p; gb %p; gn %p; bASK %s; bASK %s",
	   (void *)restp, (void *)atoms, (void *)gb, (void *)gn, 
	   bool_names[bASK], bool_names[bVerb]);
}

t_aa_names *get_aa_names(void) 
{
  /* Read the database in aminoacids.dat */
  t_aa_names *aan;
  
  snew(aan,1);
  
  aan->n = get_strings("aminoacids.dat",&(aan->aa));
  /* qsort(aan->aa,aan->n,sizeof(aan->aa[0]),strcmp);*/
  
  return aan;
}

bool is_protein(t_aa_names *aan,char *resnm)
{
  /* gives true if resnm occurs in aminoacids.dat */
  int i;
  
  for(i=0; (i<aan->n); i++)
    if (strcasecmp(aan->aa[i],resnm) == 0)
      return TRUE;
  
  return FALSE;
}

void done_aa_names(t_aa_names **aan)
{
  /* Free memory */
  int i;
  
  for(i=0; (i<(*aan)->n); i++)
    sfree((*aan)->aa[i]);
  sfree((*aan)->aa);
  sfree(*aan);
  *aan = NULL;
}

void analyse(t_atoms *atoms,t_block *gb,char ***gn,bool bASK,bool bVerb)
{
  t_aa_names *aan;
  eRestp  *restp;
  char    *resnm;
  atom_id *aid;
  int     nra;
  int     i,j;

  if (bVerb)
    printf("Analysing residue names:\n");
  snew(restp,atoms->nres);
  aid=mk_aid(atoms,restp,etOther,&nra,TRUE);
  add_grp(gb,gn,nra,aid,"System"); 
  sfree(aid);

  aan = get_aa_names();
  for(i=0; (i<atoms->nres); i++) {
    resnm=*atoms->resname[i];
    if ((restp[i] == etOther) && is_protein(aan,resnm))
      restp[i] = etProt;
    if (restp[i] == etOther)
      for(j=0; (j<NDNA);  j++) {
	if (strcasecmp(Sugars[j],resnm) == 0)
	  restp[i] = etDNA;
      }
  }
  p_status(atoms->nres,restp,bVerb);
  done_aa_names(&aan);
  
  /* Protein */
  aid=mk_aid(atoms,restp,etProt,&nra,TRUE);
  if (nra > 0) 
    analyse_prot(restp,atoms,gb,gn,bASK,bVerb);
  
  sfree(aid);

  /* Non-Protein */
  aid=mk_aid(atoms,restp,etProt,&nra,FALSE);
  if ((nra > 0) && (nra < atoms->nr))
    add_grp(gb,gn,nra,aid,"Non-Protein"); 
  sfree(aid);

  /* DNA */
  aid=mk_aid(atoms,restp,etDNA,&nra,TRUE);
  if (nra > 0) {
    add_grp(gb,gn,nra,aid,"DNA"); 
    analyse_dna(restp,atoms,gb,gn,bASK,bVerb);
  }
  sfree(aid);

  /* Other */
  analyse_other(restp,atoms,gb,gn,bASK,bVerb);
  aid=mk_aid(atoms,restp,etOther,&nra,TRUE);
  if ((nra > 0) && (nra < atoms->nr))
    add_grp(gb,gn,nra,aid,"Other"); 
  sfree(aid);
  sfree(restp);
}

void check_index(char *gname,int n,atom_id index[],char *traj,int natoms)
{
  int i;
  
  for(i=0; i<n; i++)
    if (index[i] >= natoms)
      gmx_fatal(FARGS,"%s atom number (index[%d]=%d) is larger than the number of atoms in %s (%d)",
		  gname ? gname : "Index",i+1, index[i]+1,
		  traj ? traj : "the trajectory",natoms);
}

t_block *init_index(char *gfile, char ***grpname)
{
  FILE     *in;
  t_block  *b;
  int      a,maxentries;
  int      i,j,ng;
  char     line[STRLEN],*pt,str[STRLEN];

  in=ffopen(gfile,"r");
  snew(b,1);
  get_a_line(in,line,STRLEN);
  if ( line[0]=='[' ) {
    /* new format */
    b->nr=0;
    b->index=NULL;
    b->nra=0;
    b->a=NULL;
    *grpname=NULL;
    maxentries=0;
    do {
      if (get_header(line,str)) {
	b->nr++;
	srenew(b->index,b->nr+1);
	srenew(*grpname,b->nr);
	if (b->nr==1)
	  b->index[0]=0;
	b->index[b->nr]=b->index[b->nr-1];
	(*grpname)[b->nr-1]=strdup(str);
      } else {
	pt=line;
	while ((i=sscanf(pt,"%s",str)) == 1) {
	  i=b->index[b->nr];
	  if (i>=maxentries) {
	    maxentries+=1024;
	    srenew(b->a,maxentries);
	  }
	  b->a[i]=atoi(str)-1;
	  b->index[b->nr]++;
	  (b->nra)++;
	  pt=strstr(pt,str)+strlen(str);
	}
      }
    } while (get_a_line(in,line,STRLEN));
  } 
  else {
    /* old format */
    sscanf(line,"%d%d",&b->nr,&b->nra);
    snew(b->index,b->nr+1);
    snew(*grpname,b->nr);
    b->index[0]=0;
    snew(b->a,b->nra);
    for (i=0; (i<b->nr); i++) {
      fscanf(in,"%s%d",str,&ng);
      (*grpname)[i]=strdup(str);
      b->index[i+1]=b->index[i]+ng;
      if (b->index[i+1] > b->nra)
	gmx_fatal(FARGS,"Something wrong in your indexfile at group %s",str);
      for(j=0; (j<ng); j++) {
	fscanf(in,"%d",&a);
	b->a[b->index[i]+j]=a;
      }
    }
  }
  ffclose(in);

  return b;
}

static void minstring(char *str)
{
  int i;

  for (i=0; (i < (int)strlen(str)); i++) 
    if (str[i]=='-')
      str[i]='_';
}

int find_group(char s[], int ngrps, char **grpname)
{
  int aa, i, n;
  char string[STRLEN];
  bool bMultiple;
  
  bMultiple = FALSE;
  n = strlen(s);
  aa=NOTSET;
  /* first look for whole name match */
  if (aa==NOTSET)
    for(i=0; i<ngrps; i++)
      if (strcasecmp_min(s,grpname[i])==0) {
	if(aa!=NOTSET)
	  bMultiple = TRUE;
	aa=i;
      }
  /* second look for first string match */
  if (aa==NOTSET)
    for(i=0; i<ngrps; i++)
      if (strncasecmp_min(s,grpname[i],n)==0) {
	if(aa!=NOTSET)
	  bMultiple = TRUE;
	aa=i;
      }
  /* last look for arbitrary substring match */
  if (aa==NOTSET) {
    upstring(s);
    minstring(s);
    for(i=0; i<ngrps; i++) {
      strcpy(string, grpname[i]);
      upstring(string);
      minstring(string);
      if (strstr(s,string)!=NULL) {
	if(aa!=NOTSET)
	  bMultiple = TRUE;
	aa=i;
      }
    }
  }
  if (bMultiple) {
    printf("Error: Multiple groups '%s' selected\n", s);
    aa=NOTSET;
  }
  return aa;
}

static int qgroup(int *a, int ngrps, char **grpname)
{
  char s[STRLEN];
  int  aa;
  bool bInRange;
  
  do {
  fprintf(stderr,"Select a group: ");
  do {
    if ( scanf("%s",s)!=1 ) 
      gmx_fatal(FARGS,"Cannot read from input");
  trim(s); /* remove spaces */
  } while (strlen(s)==0);
  aa = atoi(s);
    if (aa==0 && strcmp(s,"0")!=0 ) /* string entered */
      aa = find_group(s, ngrps, grpname);
    bInRange = aa>=0 && aa<ngrps;
    if (!bInRange)
      printf("Error: No such group '%s'\n", s);
  } while (!bInRange);
  printf("Selected %d: '%s'\n", aa, grpname[aa]);
  *a = aa;
  return aa;
}

static void rd_groups(t_block *grps,char **grpname,char *gnames[],
		      int ngrps,int isize[],atom_id *index[],int grpnr[])
{
  int i,j,gnr1;

  if (grps->nr==0)
    gmx_fatal(FARGS,"Error: no groups in indexfile");
  for(i=0; (i<grps->nr); i++)
    fprintf(stderr,"Group %5d (%12s) has %5d elements\n",i,grpname[i],
	   grps->index[i+1]-grps->index[i]);
  for(i=0; (i<ngrps); i++) {
    if (grps->nr > 1)
      do {
	gnr1=qgroup(&grpnr[i], grps->nr, grpname);
	if ((gnr1<0) || (gnr1>=grps->nr))
	  fprintf(stderr,"Select between %d and %d.\n",0,grps->nr-1);
      }	while ((gnr1<0) || (gnr1>=grps->nr));
    else {
      fprintf(stderr,"There is one group in the index\n");
      gnr1=0;
    }
    gnames[i]=strdup(grpname[gnr1]);
    isize[i]=grps->index[gnr1+1]-grps->index[gnr1];
    snew(index[i],isize[i]);
    for(j=0; (j<isize[i]); j++)
      index[i][j]=grps->a[grps->index[gnr1]+j];
  }
}

void rd_index(char *statfile,int ngrps,int isize[],
	      atom_id *index[],char *grpnames[])
{
  char    **gnames;
  t_block *grps;
  int     *grpnr;
  
  snew(grpnr,ngrps);
  if (!statfile)
    gmx_fatal(FARGS,"No index file specified");
  grps=init_index(statfile,&gnames);
  rd_groups(grps,gnames,grpnames,ngrps,isize,index,grpnr);
}

void rd_index_nrs(char *statfile,int ngrps,int isize[],
		  atom_id *index[],char *grpnames[],int grpnr[])
{
  char    **gnames;
  t_block *grps;
  
  if (!statfile)
    gmx_fatal(FARGS,"No index file specified");
  grps=init_index(statfile,&gnames);
  
  rd_groups(grps,gnames,grpnames,ngrps,isize,index,grpnr);
}

void get_index(t_atoms *atoms, char *fnm, int ngrps,
	       int isize[], atom_id *index[],char *grpnames[])
{
  char    ***gnames;
  t_block *grps;
  int     *grpnr;
  
  snew(grpnr,ngrps);
  snew(gnames,1);
  if (fnm != NULL) {
    grps=init_index(fnm,gnames);
  }
  else {
    snew(grps,1);
    snew(grps->index,1);
    analyse(atoms,grps,gnames,FALSE,FALSE);
  } 
  rd_groups(grps,*gnames,grpnames,ngrps,isize,index,grpnr);
}

t_cluster_ndx *cluster_index(char *ndx)
{
  t_cluster_ndx *c;
  int i;
  
  snew(c,1);
  c->clust     = init_index(ndx,&c->grpname);
  c->maxframe = -1;
  for(i=0; (i<c->clust->nra); i++)
    c->maxframe = max(c->maxframe,c->clust->a[i]);
  fprintf(stdlog ? stdlog : stdout,
	  "There are %d clusters containing %d structures, highest framenr is %d\n",
	  c->clust->nr,c->clust->nra,c->maxframe);
  if (debug) {
    pr_block(debug,0,"clust",c->clust,TRUE);
    for(i=0; (i<c->clust->nra); i++)
      if ((c->clust->a[i] < 0) || (c->clust->a[i] > c->maxframe))
	gmx_fatal(FARGS,"Ramge check error for c->clust->a[%d] = %d\n"
		  "should be within 0 and %d",i,c->clust->a[i],c->maxframe+1);
  }
  c->inv_clust=make_invblock(c->clust,c->maxframe);
	  
  return c;
}

