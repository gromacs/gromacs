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
 * Gnomes, ROck Monsters And Chili Sauce
 */
#include "sysstuff.h"
#include "smalloc.h"
#include "typedefs.h"
#include "string2.h"
#include "rdgroup.h"
#include "futil.h"
#include "fatal.h"
#include "index.h"

t_block *init_index(char *gfile, char ***grpname)
{
  FILE     *in;
  t_block  *b;
  int      a;
  int      i,j,ng,nratot;
  char     name[128];

  in=ffopen(gfile,"r");
  snew(b,1);
  fscanf(in,"%d%d",&b->nr,&b->nra);
  snew(b->index,b->nr+1);
  snew(*grpname,b->nr);
  b->index[0]=0;
  snew(b->a,b->nra);
  for (i=0; (i<b->nr); i++) {
    fscanf(in,"%s%d",name,&ng);
    (*grpname)[i]=strdup(name);
    b->index[i+1]=b->index[i]+ng;
    if (b->index[i+1] > b->nra)
      fatal_error(0,"Something wrong in your indexfile at group %s",name);
    for(j=0; (j<ng); j++) {
      fscanf(in,"%d",&a);
      b->a[b->index[i]+j]=a;
    }
    /* while (fgetc(in)!='\n'); */
  }
  ffclose(in);

  return b;
}

static int qgroup(int *a)
{
  printf("Select a group: ");
  scanf("%d",a);
  return *a;
}

static int rd_pascal_set(int **set)
{
  char buf[STRLEN+1];
  char *bb[STRLEN],*ptr;
  int   i,j,nb,ns,ne;
  
  printf("Give a set of numbers (eg.: 1..8,13,17..21)\n");
  fgets(buf,STRLEN,stdin);
  
  nb=0;
  bb[nb++]=buf;
  for(i=0; (i<(int)strlen(buf)); i++)
    if (buf[i] == ',') {
      buf[i] = '\0';
      if (buf[i+1] != '\0')
	bb[nb++] = &(buf[i+1]);
    }
  ns=0;
  for(i=0; (i<nb); i++) {
    if ((ptr=strstr(bb[i],"..")) != NULL) {
      ptr[0]=' ', ptr[1]=' ';
      sscanf(bb[i],"%d %d",&nb,&ne);
      for(j=nb; (j<=ne); j++)
	(*set)[ns++]=j;
    }
    else {
      sscanf(bb[i],"%d",&nb);
      (*set)[ns++]=nb;
    }
  }
  printf("\nI found the following set:\n");
  for(i=0; (i<ns); i++)
    printf("%4d ",(*set)[i]);
  printf("\n\n");
  
  return ns;
}

static void rd_groups(t_block *grps,char **grpname,char *gnames[],
		      int ngrps,int isize[],atom_id *index[],int grpnr[])
{
  int i,j,gnr1;

  fprintf(stderr,"Reading groups\n");
  if (grps->nr==0) {
    fprintf(stderr,"Error: no groups in indexfile\n");
    exit(1);
  }
  printf("There are %d groups in the index\n",grps->nr);
  for(i=0; (i<grps->nr); i++)
    printf("Group %5d (%12s) has %5d elements\n",i,grpname[i],
	   grps->index[i+1]-grps->index[i]);
  for(i=0; (i<ngrps); i++) {
    if (grps->nr > 1)
      do {
	gnr1=qgroup(&grpnr[i]);
	if ((gnr1<0) || (gnr1>=grps->nr))
	  printf("Select between %d and %d.\n",0,grps->nr-1);
      }	while ((gnr1<0) || (gnr1>=grps->nr));
    else 
      gnr1=0;
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
  if (!statfile) {
    fprintf(stderr,"No index file specified\n");
    exit(1);
  }
  grps=init_index(statfile,&gnames);
  rd_groups(grps,gnames,grpnames,ngrps,isize,index,grpnr);
}

void rd_index_nrs(char *statfile,int ngrps,int isize[],
	        atom_id *index[],char *grpnames[],int grpnr[])
{
  char    **gnames;
  t_block *grps;
  
  if (!statfile) {
    fprintf(stderr,"No index file specified\n");
    exit(1);
  }
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
    fprintf(stderr,"Making index\n");
    snew(grps,1);
    snew(grps->index,1);
    analyse(atoms,grps,gnames,FALSE,FALSE);
  } 
  rd_groups(grps,*gnames,grpnames,ngrps,isize,index,grpnr);
}






