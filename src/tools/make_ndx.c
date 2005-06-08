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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include "sysstuff.h"
#include "strdb.h"
#include "futil.h"
#include "macros.h"
#include "string2.h"
#include "statutil.h"
#include "confio.h"
#include "copyrite.h"
#include "typedefs.h"
#include "index.h"
#include "smalloc.h"
#include "vec.h"
#include "index.h"

#define MAXNAMES 20

bool bCase=FALSE;

static int or_groups(atom_id nr1,atom_id *at1,atom_id nr2,atom_id *at2,
		     atom_id *nr,atom_id *at)
{
  atom_id i1,i2,max=0;
  bool bNotIncr;

  *nr=0;
  
  bNotIncr=FALSE;
  for(i1=0; i1<nr1; i1++) {
    if ((i1>0) && (at1[i1] <= max))
      bNotIncr=TRUE;
    max=at1[i1];
  }
  for(i1=0; i1<nr2; i1++) {
    if ((i1>0) && (at2[i1] <= max))
      bNotIncr=TRUE;
    max=at2[i1];
  }

  if (bNotIncr)
    printf("One of your groups is not ascending\n");
  else {
    i1=0;
    i2=0;
    *nr=0;
    while ((i1 < nr1) || (i2 < nr2)) {
      if ((i2 == nr2) || ((i1<nr1) && (at1[i1] < at2[i2]))) {
	at[*nr]=at1[i1];
	(*nr)++;
	i1++;
      }
      else {
	if ((i2<nr2) && ((i1==nr1) || (at1[i1] > at2[i2]))) {
	  at[*nr]=at2[i2];
	  (*nr)++;
	}
	i2++;
      }
    }
    
    printf("Merged two groups with OR: %u %u -> %u\n",nr1,nr2,*nr);
  }

  return *nr;
}

static int and_groups(atom_id nr1,atom_id *at1,atom_id nr2,atom_id *at2,
		      atom_id *nr,atom_id *at)
{
  atom_id i1,i2;
  
  *nr=0;
  for (i1=0; i1<nr1; i1++) {
    for (i2=0; i2<nr2; i2++) {
      if (at1[i1]==at2[i2]) {
	at[*nr]=at1[i1];
	(*nr)++;
      }
    }
  }

  printf("Merged two groups with AND: %u %u -> %u\n",nr1,nr2,*nr);

  return *nr;
}

static bool isalnum_star(char c)
{
  return (isalnum(c) || (c=='*'));
}

static int parse_names(char **string,int *n_names,char **names)
{
  int i;
  
  *n_names=0;
  while ((isalnum_star((*string)[0]) || ((*string)[0]==' '))) {
    if (isalnum_star((*string)[0])) {
      if (*n_names >= MAXNAMES) 
	gmx_fatal(FARGS,"To many names: %d\n",*n_names+1);
      i=0;
      while (isalnum_star((*string)[i])) {
	names[*n_names][i]=(*string)[i];
	i++;
	if (i==5) {
	  printf("Name to long: %d characters\n",i);	  
	  return 0;
	}
      }
      names[*n_names][i]='\0';
      if (!bCase)
	upstring(names[*n_names]);
      *string += i;
      (*n_names)++;
    }
    else
      (*string)++;
  }  

  return *n_names;
}

static bool parse_int(char **string,int *nr)
{
  bool bRet;

  while ((*string)[0]==' ')
    (*string)++;

  bRet=FALSE;

  if (isdigit((*string)[0])) {
    *nr=(*string)[0]-'0';
    (*string)++;
    while (isdigit((*string)[0])) {
      *nr = (*nr)*10+(*string)[0]-'0';
      (*string)++;
    }
    bRet=TRUE;
  }
  else
    *nr=NOTSET;

  return bRet;
}

static bool isquote(char c)
{
  char quotes[] = "'\"";
  
  return strchr(quotes, c)!=NULL;
}

static bool parse_string(char **string,int *nr, int ngrps, char **grpname)
{
  char *s, *sp;
  char c;
  
  while ((*string)[0]==' ')
    (*string)++;

  (*nr) = NOTSET;
  if (isquote((*string)[0])) {
    c=(*string)[0];
    (*string)++;
    s = strdup((*string));
    sp = strchr(s, c);
    if (sp!=NULL) {
      (*string) += sp-s + 1;
      sp[0]='\0';
      (*nr) = find_group(s, ngrps, grpname);
    }
  }
  
  return (*nr) != NOTSET;
}

static int select_atomnumbers(char **string,t_atoms *atoms,atom_id n1,
			      atom_id *nr,atom_id *index,char *gname)
{
  char    buf[STRLEN];
  int     i,up;

  *nr=0;
  while ((*string)[0]==' ')
    (*string)++;
  if ((*string)[0]=='-') {
    (*string)++;
    parse_int(string,&up);
    if ((n1<1) || (n1>atoms->nr) || (up<1) || (up>atoms->nr))
      printf("Invalid atom range\n");
    else {
      for(i=n1-1; i<=up-1; i++) {
	index[*nr]=i;
	(*nr)++;
      }  
      printf("Found %u atom%s in range %u-%d\n",*nr,(*nr==1)?"":"s",n1,up);
      if (n1==up)
	sprintf(buf,"a_%u",n1);
      else
	sprintf(buf,"a_%u-%d",n1,up);
      strcpy(gname,buf);
    }
  }
  else {
    i=n1;
    sprintf(gname,"a");
    do {
      if ((i-1>=0) && (i-1<atoms->nr)) {
	index[*nr] = i-1;
	(*nr)++;
	sprintf(buf,"_%d",i);
	strcat(gname,buf);
      } else {
	printf("Invalid atom number %d\n",i);
	*nr = 0;
      }
    } while ((*nr!=0) && (parse_int(string,&i)));
  }
  
  return *nr;
}

static int select_residuenumbers(char **string,t_atoms *atoms,atom_id n1,
				 atom_id *nr,atom_id *index,char *gname)
{
  char    buf[STRLEN];
  int     j,resnr;
  int     i,up;

  *nr=0;
  while ((*string)[0]==' ')
    (*string)++;
  if ((*string)[0]=='-') {
    (*string)++;
    parse_int(string,&up);

    for(i=0; i<atoms->nr; i++) {
      resnr=atoms->atom[i].resnr;
      for(j=n1; (j<=up); j++) {
	if (resnr==j-1) {
	  index[*nr]=i;
	  (*nr)++;
	}
      }
    }
    printf("Found %u atom%s with res.nr. in range %u-%d\n",
	   *nr,(*nr==1)?"":"s",n1,up);
    if (n1==up)
      sprintf(buf,"r_%u",n1);
    else
      sprintf(buf,"r_%u-%d",n1,up);
    strcpy(gname,buf);
  }
  else {
    j=n1;
    sprintf(gname,"r");
    do {
      for(i=0; i<atoms->nr; i++) {
	if (atoms->atom[i].resnr==j-1) {
	index[*nr]=i;
	(*nr)++;
	}
      }
      sprintf(buf,"_%d",j);
      strcat(gname,buf);
    } while (parse_int(string,&j));
  }
  
  return *nr;
}

static bool atoms_from_residuenumbers(t_atoms *atoms,int group,t_block *block,
				      atom_id *nr,atom_id *index,char *gname)
{
  int i,j,j0,j1,resnr,nres;
  
  j0=block->index[group];
  j1=block->index[group+1];
  nres = atoms->atom[atoms->nr-1].resnr;
  for(j=j0; j<j1; j++)
    if (block->a[j]>=nres) {
      printf("Index %s contains number>nres (%d>%d)\n",
	     gname,block->a[j]+1,nres);
      return FALSE;
    }
  for(i=0; i<atoms->nr; i++) {
    resnr=atoms->atom[i].resnr;
    for (j=j0; j<j1; j++)
      if (block->a[j]==resnr) {
	index[*nr]=i;
	(*nr)++;
	break;
      }
  }
  printf("Found %u atom%s in %d residues from group %s\n",
	 *nr,(*nr==1)?"":"s",j1-j0,gname);
  return *nr;
}

static bool comp_name(char *name,char *search)
{
  int n;

  n=strlen(search)-1;

  return (((search[n]!='*') && 
	   (bCase ? strcmp(name,search) : strcasecmp(name,search))) ||
	  ((search[n]=='*') && 
	   (bCase ? strncmp(name,search,n) : strncasecmp(name,search,n))));
}

static int select_chainnames(t_atoms *atoms,int n_names,char **names,
			     atom_id *nr,atom_id *index)
{
  char name[2];
  int j;
  atom_id i;
  
  name[1]=0;
  *nr=0;
  for(i=0; i<atoms->nr; i++) {
    name[0]=atoms->atom[i].chain;
    j=0; 
    while (j<n_names && comp_name(name,names[j])) 
      j++;
    if (j<n_names) {
      index[*nr]=i;
      (*nr)++;
    }
  }
  printf("Found %u atom%s with chain identifier%s",
	 *nr,(*nr==1)?"":"s",(n_names==1)?"":"s");
  for(j=0; (j<n_names); j++)
    printf(" %s",names[j]);
  printf("\n");

  return *nr;
}

static int select_atomnames(t_atoms *atoms,int n_names,char **names,
			    atom_id *nr,atom_id *index,bool bType)
{
  char *name;
  int j;
  atom_id i;
  
  *nr=0;
  for(i=0; i<atoms->nr; i++) {
    if (bType)
      name=*(atoms->atomtype[i]);
    else
      name=*(atoms->atomname[i]);
    j=0; 
    while (j<n_names && comp_name(name,names[j])) 
      j++;
    if (j<n_names) {
      index[*nr]=i;
      (*nr)++;
    }
  }
  printf("Found %u atoms with %s%s",
	 *nr,bType ? "type" : "name",(n_names==1)?"":"s");
  for(j=0; (j<n_names); j++)
    printf(" %s",names[j]);
  printf("\n");

  return *nr;
}

static int select_residuenames(t_atoms *atoms,int n_names,char **names,
			       atom_id *nr,atom_id *index)
{
  char *name;
  int j;
  atom_id i;

  *nr=0;
  for(i=0; i<atoms->nr; i++) {
    name=*(atoms->resname[atoms->atom[i].resnr]);
    j=0; 
    while (j<n_names && comp_name(name,names[j])) 
      j++;
    if (j<n_names) {
      index[*nr]=i;
      (*nr)++;
    }
  }
  printf("Found %u atoms with residue name%s",*nr,(n_names==1)?"":"s");
  for(j=0; (j<n_names); j++)
    printf(" %s",names[j]);
  printf("\n");

  return *nr;
}

static void copy2block(int n,atom_id *index,t_block *block)
{
  int i,n0;

  block->nr++;
  n0=block->nra;
  block->nra=n0+n;
  srenew(block->index,block->nr+1);
  block->index[block->nr]=n0+n;
  srenew(block->a,n0+n);
  for(i=0; (i<n); i++)
    block->a[n0+i]=index[i];
}

static void make_gname(int n,char **names,char *gname)
{
  int i;
  
  strcpy(gname,names[0]);
  for (i=1; i<n; i++) {
    strcat(gname,"_");
    strcat(gname,names[i]);
  }
}

static void copy_group(int g,t_block *block,atom_id *nr,atom_id *index)
{
  int i,i0;
  
  i0=block->index[g];
  *nr=block->index[g+1]-i0;
  for (i=0; i<=*nr; i++)
    index[i]=block->a[i0+i];
}

static void remove_group(int nr,int nr2,t_block *block,char ***gn)
{
  int i,j,shift;
  char *name;
  
  if (nr2==NOTSET)
    nr2=nr;
  
  for(j=0; j<=nr2-nr; j++) {
    if ((nr<0) || (nr>=block->nr))
      printf("Group %d does not exist\n",nr+j);
    else {
      shift=block->index[nr+1]-block->index[nr];
      for(i=block->index[nr+1]; i<block->nra; i++)
	block->a[i-shift]=block->a[i];
      
      for(i=nr; i<block->nr; i++)
	block->index[i]=block->index[i+1]-shift;
      name = strdup((*gn)[nr]);
      sfree((*gn)[nr]);
      for(i=nr; i<block->nr-1; i++) {  
	(*gn)[i]=(*gn)[i+1];
      }
      block->nr--;
      block->nra=block->index[block->nr];
      printf("Removed group %d '%s'\n",nr+j,name);
      sfree(name);
    }
  }
}

static void split_group(t_atoms *atoms,int sel_nr,t_block *block,char ***gn,
			bool bAtom)
{
  char buf[STRLEN],*name;
  int i,nr;
  atom_id a,n0,n1,max;

  if (bAtom)
    printf("Splitting group %d '%s' into residues\n",sel_nr,(*gn)[sel_nr]);
  else
    printf("Splitting group %d '%s' into atoms\n",sel_nr,(*gn)[sel_nr]);
  
  max=block->nra;
  n0=block->index[sel_nr];
  n1=block->index[sel_nr+1];
  for (i=n0; i<n1; i++) {
    a=block->a[i];
    nr=atoms->atom[a].resnr;
    name=*(atoms->resname[nr]);
    if (bAtom || (i==n0) || (atoms->atom[block->a[i-1]].resnr!=nr)) { 
      if (i>n0)
	block->index[block->nr]=block->nra;
      block->nr++;
      srenew(block->index,block->nr+1);
      srenew(*gn,block->nr);
      if (bAtom)
	sprintf(buf,"%s_%s_%u",(*gn)[sel_nr],*atoms->atomname[a],a+1);
      else
	sprintf(buf,"%s_%s_%d",(*gn)[sel_nr],name,nr+1);
      (*gn)[block->nr-1]=strdup(buf);
    }
    if (block->nra == max) {
      max+=20;
      srenew(block->a,max);
    }
    block->a[block->nra]=a;
    block->nra++;
  }
  block->index[block->nr]=block->nra;
}

static int split_chain(t_atoms *atoms,rvec *x,
			int sel_nr,t_block *block,char ***gn)
{
  char    buf[STRLEN];
  int     j,nchain;
  atom_id i,a,max,natoms,*start=NULL,*end=NULL,ca_start,ca_end;
  rvec    vec;

  natoms=atoms->nr;
  nchain=0;
  ca_start=0;

  while (ca_start<natoms) {
    while((ca_start<natoms) && strcmp(*atoms->atomname[ca_start],"CA"))
      ca_start++;
    if (ca_start<natoms) {
      srenew(start,nchain+1);
      srenew(end,nchain+1);
      start[nchain]=ca_start;
      while ((start[nchain]>0) && 
	     (atoms->atom[start[nchain]-1].resnr==atoms->atom[ca_start].resnr))
	start[nchain]--;

      i=ca_start;
      do {
	ca_end=i;
	do {
	  i++;
	} while ((i<natoms) && strcmp(*atoms->atomname[i],"CA"));
	if (i<natoms)
	  rvec_sub(x[ca_end],x[i],vec);
      } while ((i<natoms) && (norm(vec)<0.45));
      
      end[nchain]=ca_end;
      while ((end[nchain]+1<natoms) && 
	    (atoms->atom[end[nchain]+1].resnr==atoms->atom[ca_end].resnr)) 
	end[nchain]++;
      ca_start=end[nchain]+1;
      nchain++;
    }
  }
  if (nchain==1)
    printf("Found 1 chain, will not split\n");
  else
    printf("Found %d chains\n",nchain);
  for (j=0; j<nchain; j++)
    printf("%d:%6u atoms (%u to %u)\n",
	   j+1,end[j]-start[j]+1,start[j]+1,end[j]+1);

  if (nchain>1) {
    for (j=0; j<nchain; j++) {
      block->nr++;
      srenew(block->index,block->nr+1);
      srenew(*gn,block->nr);
      sprintf(buf,"%s_chain%d",(*gn)[sel_nr],j+1);
      (*gn)[block->nr-1]=strdup(buf);
      max=block->nra;
      for (i=block->index[sel_nr]; i<block->index[sel_nr+1]; i++) {
	a=block->a[i];
	if ((a>=start[j]) && (a<=end[j])) {
	  if (block->nra==max) {
	    max+=20;
	    srenew(block->a,max);
	  }
	  block->a[block->nra]=a;
	  block->nra++;
	}
      }
      block->index[block->nr]=block->nra;
      if (block->index[block->nr-1]==block->index[block->nr])
	remove_group(block->nr-1,NOTSET,block,gn);
    }
  }
  sfree(start);
  sfree(end);

  return nchain;
}

static bool check_have_atoms(t_atoms *atoms, char *string)
{
  if ( atoms==NULL ) {
    printf("Can not process '%s' without atoms info\n", string);
    return FALSE;
  } else
    return TRUE;
}

static bool parse_entry(char **string,int natoms,t_atoms *atoms,
			t_block *block,char ***gn,
			atom_id *nr,atom_id *index,char *gname)
{
  static char **names, *ostring;
  static int  namelen=5;
  static bool bFirst=TRUE;
  int         j,n_names,sel_nr1;
  atom_id     i,nr1,*index1;
  bool        bRet,bCompl;

  if (bFirst) {
    bFirst=FALSE;
    snew(names,MAXNAMES);
    for (i=0; i<MAXNAMES; i++)
      snew(names[i],namelen);
  }

  bRet=FALSE;
  sel_nr1=NOTSET;

  while(*string[0]==' ')
    (*string)++;

  if ((*string)[0]=='!') {
    bCompl=TRUE;
    (*string)++;
     while(*string[0]==' ')
       (*string)++;
  } else 
    bCompl=FALSE;

  ostring = *string;
  
  if (parse_int(string,&sel_nr1) || 
      parse_string(string,&sel_nr1,block->nr,*gn)) {
    if ((sel_nr1>=0) && (sel_nr1<block->nr)) {
      copy_group(sel_nr1,block,nr,index);
      strcpy(gname,(*gn)[sel_nr1]);
      printf("Copied index group %d '%s'\n",sel_nr1,(*gn)[sel_nr1]);
      bRet=TRUE;
    } else
      printf("Group %d does not exist\n",sel_nr1); 
  } 
  else if ((*string)[0]=='a') {
    (*string)++;
    if (check_have_atoms(atoms, ostring)) {
    if (parse_int(string,&sel_nr1)) {
      bRet=select_atomnumbers(string,atoms,sel_nr1,nr,index,gname);
    } 
    else if (parse_names(string,&n_names,names)) {
      bRet=select_atomnames(atoms,n_names,names,nr,index,FALSE);
      make_gname(n_names,names,gname);
      }
    }
  }
  else if ((*string)[0]=='t') {
    (*string)++;
    if (check_have_atoms(atoms, ostring) &&
        parse_names(string,&n_names,names)) {
      if (atoms->atomtype == NULL)
	printf("Need a run input file to select atom types\n");
      else {
	bRet=select_atomnames(atoms,n_names,names,nr,index,TRUE);
	make_gname(n_names,names,gname);
      }
    }
  }
  else if (strncmp(*string,"res",3)==0) {
    (*string)+=3;
    if ( check_have_atoms(atoms, ostring) &&
	 parse_int(string,&sel_nr1) &&
	 (sel_nr1>=0) && (sel_nr1<block->nr) ) {
      bRet=atoms_from_residuenumbers(atoms,
				     sel_nr1,block,nr,index,(*gn)[sel_nr1]);
      sprintf(gname,"atom_%s",(*gn)[sel_nr1]);
    }
  }
  else if ((*string)[0]=='r') {
    (*string)++;
    if (check_have_atoms(atoms, ostring)) {
    if (parse_int(string,&sel_nr1)) {
      bRet=select_residuenumbers(string,atoms,sel_nr1,nr,index,gname);
    } 
    else if (parse_names(string,&n_names,names)) {
      bRet=select_residuenames(atoms,n_names,names,nr,index);
      make_gname(n_names,names,gname);
    }
  }
  }
  else if (strncmp(*string,"chain",5)==0) {
    (*string)+=5;
    if (check_have_atoms(atoms, ostring) &&
        parse_names(string,&n_names,names)) {
      bRet=select_chainnames(atoms,n_names,names,nr,index);
      sprintf(gname,"ch%s",names[0]);
      for (i=1; i<n_names; i++)
	strcat(gname,names[i]);
    }
  }
  if (bRet && bCompl) {
    snew(index1,natoms-*nr);
    nr1=0;
    for(i=0; i<natoms; i++) {
      j=0;
      while ((j<*nr) && (index[j] != i))
	j++;
      if (j==*nr) {
	if (nr1 >= natoms-*nr) {
	  printf("There are double atoms in your index group\n");
	  break;
	}
	index1[nr1]=i;
	nr1++;
      }
    }
    *nr=nr1;
    for(i=0; i<nr1; i++) 
      index[i]=index1[i];
    sfree(index1);
    
    for (i=strlen(gname)+1; i>0; i--)
      gname[i]=gname[i-1];
    gname[0]='!';
    printf("Complemented group: %u atoms\n",*nr);
  }
  
  return bRet;
}

static void list_residues(t_atoms *atoms)
{
  int i,j,start,end,prev_resnr,resnr;
  bool bDiff;

  /* Print all the residues, assuming continuous resnr count */ 
  start = atoms->atom[0].resnr;
  prev_resnr = start;
  for(i=0; i<atoms->nr; i++) {
    resnr = atoms->atom[i].resnr;
    if ((resnr != prev_resnr) || (i==atoms->nr-1)) {
      if ((bDiff=strcmp(*atoms->resname[resnr],*atoms->resname[start])) || 
	  (i==atoms->nr-1)) {
	if (bDiff)
	  end = prev_resnr;
	else
	  end = resnr;
	if (end < start+3)
	  for(j=start; j<=end; j++)
	    printf("%4d %-5s",j+1,*(atoms->resname[j]));
	else
	  printf(" %4d - %4d %-5s  ",start+1,end+1,*(atoms->resname[start]));
	start = resnr;
      }
    }
    prev_resnr = resnr;
  }   
  printf("\n");
}

static void edit_index(int natoms, t_atoms *atoms,rvec *x,t_block *block, char ***gn, bool bVerbose)
{
  static char **atnames, *ostring;
  static bool bFirst=TRUE;
  static int  namelen=5;
  char inp_string[STRLEN],*string;
  char gname[STRLEN],gname1[STRLEN],gname2[STRLEN];
  int  i,i0,i1,sel_nr,sel_nr2,newgroup;
  atom_id nr,nr1,nr2,*index,*index1,*index2;
  bool bAnd,bOr,bPrintOnce;
  
  if (bFirst) {
    bFirst=FALSE;
    snew(atnames,MAXNAMES);
    for (i=0; i<MAXNAMES; i++)
      snew(atnames[i],namelen);
  }

  string=NULL;

  snew(index,natoms);
  snew(index1,natoms);
  snew(index2,natoms);

  newgroup=NOTSET;
  bPrintOnce = TRUE;
  do { 
    gname1[0]='\0';
    if (bVerbose || bPrintOnce || newgroup!=NOTSET) {
      printf("\n");
      if (bVerbose || bPrintOnce || newgroup==NOTSET) {
	i0=0;
	i1=block->nr;
      } else {
	i0=newgroup;
	i1=newgroup+1;
      }
      for(i=i0; i<i1; i++)
      printf("%3d %-20s: %5u atoms\n",i,(*gn)[i],
	     block->index[i+1]-block->index[i]);
      newgroup=NOTSET;
    }
    if (bVerbose || bPrintOnce) {
      printf("\n");
      printf(" nr : group       !   'name' nr name   'splitch' nr    Enter: list groups\n");
      printf(" 'a': atom        &   'del' nr         'splitres' nr   'l': list residues\n");
      printf(" 't': atom type   |   'keep' nr        'splitat' nr    'h': help\n");
      printf(" 'r': residue         'res' nr         'chain' char\n");
      printf(" \"name\": group        'case': case %s         'q': save and quit\n",
	     bCase ? "insensitive" : "sensitive  "); 
      bPrintOnce = FALSE;
    }
    printf("\n");
    printf("> ");
    fgets(inp_string,STRLEN,stdin);
    inp_string[strlen(inp_string)-1]=0;
    printf("\n");
    string=inp_string;
    while (string[0]==' ')
      string++;

    ostring = string;
    nr=0;
    if (string[0] == 'h') {
      printf(" nr                : selects an index group by number or quoted string.\n");
      printf("                     The string is first matched against the whole group name,\n");
      printf("                     then against the beginning and finally against an\n");
      printf("                     arbitrary substring. A multiple match is an error.\n");
      
      printf(" 'a' nr1 [nr2 ...] : selects atoms, atom numbering starts at 1.\n");
      printf(" 'a' nr1 - nr2     : selects atoms in the range from nr1 to nr2.\n"); 
      printf(" 'a' name1[*] [name2[*] ...] : selects atoms by name(s), wildcard allowed\n"); 
      printf("                               at the end of a name.\n");
      printf(" 't' type1[*] [type2[*] ...] : as 'a' but for type, run input file required.\n");
      printf(" 'r'               : analogous to 'a', but for residues.\n");
      printf(" 'chain' ch1 [ch2 ...] : selects atoms by chain identifier(s),\n");
      printf("                         not available with a .gro file as input.\n");
      printf(" !                 : takes the complement of a group with respect to all\n");
      printf("                     the atoms in the input file.\n");
      printf(" & |               : AND and OR, can be placed between any of the options\n");
      printf("                     above, the input is processed from left to right.\n");
      printf(" 'name' nr name    : rename group nr to name.\n");
      printf(" 'del' nr1 [- nr2] : deletes one group or groups in the range from nr1 to nr2.\n");
      printf(" 'keep' nr         : deletes all groups except nr.\n");
      printf(" 'case'            : make all name compares case (in)sensitive.\n"); 
      printf(" 'splitch' nr      : split group into chains using CA distances.\n");
      printf(" 'splitres' nr     : split group into residues.\n");
      printf(" 'splitat' nr      : split group into atoms.\n");
      printf(" 'res' nr          : interpret numbers in group as residue numbers\n");
      printf(" Enter             : list the currently defined groups\n");
      printf(" 'l'               : list the residues.\n");
      printf(" 'h'               : show this help.\n");
      printf(" 'q'               : save and quit.\n");
      printf("\n");
      printf(" Examples:\n");
      printf(" > 2 | 4 & r 3-5\n");
      printf(" will select all atoms from group 2 and 4 which have residue numbers\n 3, 4 or 5\n");
      printf(" > a C* & !a C CA\n");
      printf(" will select all atoms starting with 'C' but not the atoms 'C' and 'CA'\n");  
      printf(" > \"protein\" & ! \"backb\"\n");
      printf(" from the default groups, this selects all protein atoms but not the backbone\n");
      printf("\npress Enter");
      getchar();
    }
    else if (strncmp(string,"del",3)==0) {
      string+=3;
      if (parse_int(&string,&sel_nr)) {
	while(string[0]==' ')
	  string++;
	if (string[0]=='-') {
	  string++;
	  parse_int(&string,&sel_nr2);
	} else
	  sel_nr2=NOTSET;
	while(string[0]==' ')
	  string++;
	if (string[0]=='\0')
	  remove_group(sel_nr,sel_nr2,block,gn);
	else
	  printf("\nSyntax error: \"%s\"\n",string);
      }
    }
    else if (strncmp(string,"keep",4)==0) {
      string+=4;
      if (parse_int(&string,&sel_nr)) {
	remove_group(sel_nr+1,block->nr-1,block,gn);
	remove_group(0,sel_nr-1,block,gn);
      }
    }
    else if (strncmp(string,"name",4)==0) {
      string+=4;
      if (parse_int(&string,&sel_nr)) {
	if ((sel_nr>=0) && (sel_nr<block->nr)) { 
	  sscanf(string,"%s",gname);
	  sfree((*gn)[sel_nr]);
	  (*gn)[sel_nr]=strdup(gname);
	}
      }
    }
    else if (strncmp(string,"case",4)==0) {
      bCase=!bCase;
      printf("  Switched to case %s\n",bCase ? "sensitive" : "insensitive");
    }
    else if (string[0] == 'v' ) {
      bVerbose=!bVerbose;
      printf("  Turned verbose %s\n",bVerbose ? "on" : "off");
    }
    else if (string[0] == 'l') {
      if ( check_have_atoms(atoms, ostring) )
      list_residues(atoms);
    }
    else if (strncmp(string,"splitch",7)==0) {
      string+=7;
      if ( check_have_atoms(atoms, ostring) && 
	   parse_int(&string,&sel_nr) &&
	   (sel_nr>=0) && (sel_nr<block->nr)) 
	  split_chain(atoms,x,sel_nr,block,gn);
    }
    else if (strncmp(string,"splitres",8)==0 ) {
      string+=8;
      if ( check_have_atoms(atoms, ostring) &&
	   parse_int(&string,&sel_nr) &&
	   (sel_nr>=0) && (sel_nr<block->nr)) 
	  split_group(atoms,sel_nr,block,gn,FALSE);
    }
    else if (strncmp(string,"splitat",7)==0 ) {
     string+=7;
      if ( check_have_atoms(atoms, ostring) &&
	   parse_int(&string,&sel_nr) &&
	   (sel_nr>=0) && (sel_nr<block->nr)) 
	  split_group(atoms,sel_nr,block,gn,TRUE);   
    }
    else if (string[0] == '\0') {
      bPrintOnce = TRUE;
    }
    else if (string[0] != 'q') {
      nr1=-1;
      nr2=-1;
      if (parse_entry(&string,natoms,atoms,block,gn,&nr,index,gname)) {
	do {
	  while (string[0]==' ')
	    string++;
	  
	  bAnd=FALSE;
	  bOr=FALSE;
	  if (string[0]=='&')
	    bAnd=TRUE;
	  else if (string[0]=='|')
	    bOr=TRUE;
	  
	  if (bAnd || bOr) {
	    string++;
	    nr1=nr;
	    for(i=0; i<nr; i++)
	      index1[i]=index[i];
	    strcpy(gname1,gname);
	    if (parse_entry(&string,natoms,atoms,block,gn,&nr2,index2,gname2)) {
	      if (bOr) {
		or_groups(nr1,index1,nr2,index2,&nr,index);
		sprintf(gname,"%s_%s",gname1,gname2);
	      }
	      else {
		and_groups(nr1,index1,nr2,index2,&nr,index);
		sprintf(gname,"%s_&_%s",gname1,gname2);
	      }
	    }
	  }
	} while (bAnd || bOr);
      }
      while(string[0]==' ')
	string++;
      if (string[0])
	printf("\nSyntax error: \"%s\"\n",string);
      else if (nr>0) {
	copy2block(nr,index,block);
	srenew(*gn,block->nr);
	newgroup = block->nr-1;
	(*gn)[newgroup]=strdup(gname);
      }
      else
	printf("Group is empty\n");
    }
  } while (string[0]!='q');

  sfree(index);
  sfree(index1);
  sfree(index2);
}

static int block2natoms(t_block *block)
{
  int i, natoms;
  
  natoms = 0;
  for(i=0; i<block->nra; i++)
    natoms = max(natoms, block->a[i]);
  natoms++;
  
  return natoms;
}

void merge_blocks(t_block *dest, t_block *source)
{
  int     i,nra0,i0;
  
  /* count groups, srenew and fill */
  i0 = dest->nr;
  nra0 = dest->nra;
  dest->nr+=source->nr;
  srenew(dest->index, dest->nr+1);
  for(i=0; i<source->nr; i++)
    dest->index[i0+i] = nra0 + source->index[i];
  /* count atoms, srenew and fill */
  dest->nra+=source->nra;
  srenew(dest->a, dest->nra);
  for(i=0; i<source->nra; i++)
    dest->a[nra0+i] = source->a[i];
  
  /* terminate list */
  dest->index[dest->nr]=dest->nra;
  
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "Index groups are necessary for almost every gromacs program.",
    "All these programs can generate default index groups. You ONLY",
    "have to use make_ndx when you need SPECIAL index groups.",
    "There is a default index group for the whole system, 9 default", 
    "index groups are generated for proteins, a default index group",
    "is generated for every other residue name.[PAR]"
    "When no index file is supplied, also make_ndx will generate the",
    "default groups.",
    "With the index editor you can select on atom, residue and chain names",
    "and numbers.",
    "When a run input file is supplied you can also select on atom type.",
    "You can use NOT, AND and OR, you can split groups",
    "into chains, residues or atoms. You can delete and rename groups.[PAR]",
    "The atom numbering in the editor and the index file starts at 1."
  };

  static int natoms=0;
  static bool bVerbose=FALSE;
  t_pargs pa[] = {
    { "-natoms",  FALSE, etINT, {&natoms}, 
      "set number of atoms (default: read from coordinate or index file)" },
    { "-verbose", FALSE, etBOOL, {&bVerbose}, 
      "HIDDENVerbose output" }
  };
#define NPA asize(pa)

  char     title[STRLEN];
  int      nndxin;   
  char     *stxfile, **ndxinfiles, *ndxoutfile;
  bool     bNatoms;
  int      i,j;
  t_atoms  *atoms;
  rvec     *x,*v;
  matrix   box;
  t_block  *block,*block2;
  char     **gnames,**gnames2;
  t_filenm fnm[] = {
    { efSTX, "-f", NULL,     ffOPTRD  },
    { efNDX, "-n", NULL,     ffOPTRDMULT },
    { efNDX, "-o", NULL,     ffWRITE }
  };
#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  
  parse_common_args(&argc,argv,0,NFILE,fnm,NPA,pa,asize(desc),desc,
		    0,NULL);
  
  stxfile = ftp2fn_null(efSTX,NFILE,fnm);
  if (opt2bSet("-n",NFILE,fnm)) {
    nndxin = opt2fns(&ndxinfiles,"-n",NFILE,fnm);
  } else {
    nndxin = 0;
  }
  ndxoutfile = opt2fn("-o",NFILE,fnm);
  bNatoms = opt2parg_bSet("-natoms",NPA,pa);
  
  if (!stxfile && !nndxin)
    gmx_fatal(FARGS,"No input files (structure or index)");
  
  if (stxfile) {
    snew(atoms,1);
    get_stx_coordnum(stxfile,&(atoms->nr));
    init_t_atoms(atoms,atoms->nr,TRUE);
    snew(x,atoms->nr);
    snew(v,atoms->nr);
    fprintf(stderr,"\nReading structure file\n");
    read_stx_conf(stxfile,title,atoms,x,v,box);
    natoms = atoms->nr;
    bNatoms=TRUE;
  } else {
    atoms = NULL;
    x = NULL;
  }

  /* read input file(s) */
  block = new_block();
  gnames = NULL;
  printf("Going to read %d old index file(s)\n",nndxin);
  if (nndxin) {
    for(i=0; i<nndxin; i++) {
      block2 = init_index(ndxinfiles[i],&gnames2);
      srenew(gnames, block->nr+block2->nr);
      for(j=0; j<block2->nr; j++)
	gnames[block->nr+j]=gnames2[j];
      sfree(gnames2);
      merge_blocks(block, block2);
      sfree(block2->a);
      sfree(block2->index);
/*       done_block(block2); */
      sfree(block2);
    }
  } 
  else {
    snew(gnames,1);
    analyse(atoms,block,&gnames,FALSE,TRUE);
  }
  
  if (!bNatoms) {
    natoms = block2natoms(block);
    printf("Counted atom numbers up to %d in index file\n", natoms);
  }

  edit_index(natoms,atoms,x,block,&gnames,bVerbose);

  write_index(ndxoutfile,block,gnames);

  thanx(stderr);
    
  return 0;
}
