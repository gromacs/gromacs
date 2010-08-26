/*
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

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "string2.h"
#include "strdb.h"
#include "typedefs.h"
#include "macros.h"
#include "copyrite.h"
#include "smalloc.h"
#include "statutil.h"
#include "confio.h"
#include "calch.h"

typedef struct {
  char *key;
  int  nexp;
  char **exp;
} t_expansion;

t_expansion *read_expansion_map(char *fn,int *nexpand)
{
  char        ibuf[12],buf[12][10];
  char        **ptr;
  t_expansion *exp;
  int         i,k,nexp,nn;
  
  nexp=get_file(fn,&ptr);
  
  snew(exp,nexp);
  for(i=0; (i<nexp); i++) {
    /* Let scanf do the counting... */
    nn=sscanf(ptr[i],"%s%s%s%s%s%s%s%s%s%s%s",
	      ibuf,buf[0],buf[1],buf[2],buf[3],buf[4],
	      buf[5],buf[6],buf[7],buf[8],buf[9]);
    if (nn <= 1)
      break;
    exp[i].key=strdup(ibuf);
    exp[i].nexp=nn-1;
    snew(exp[i].exp,nn-1);
    for(k=0; (k<nn-1); k++)
      exp[i].exp[k]=strdup(buf[k]);
  }
  fprintf(stderr,"I found %d expansion mapping entries!\n",i);
  
  /* Clean up */
  for(i=0; (i<nexp); i++)
    sfree(ptr[i]);
  sfree(ptr);
  
  *nexpand=nexp;
  return exp;  
}

char **get_exp(int NEXP,t_expansion expansion[],char **ai,int *nexp)
{
  int  i;

  for(i=0; (i<NEXP); i++)
    if (strcmp(*ai,expansion[i].key) == 0) {
      *nexp=expansion[i].nexp;
      return expansion[i].exp;
    }
  *nexp=1;

  return ai;
}

int find_atom(char *ai,char *ri,
	      int resi,int r0,
	      int natoms,char ***aname,t_atom atom[],
	      int linec,gmx_bool bVerbose)
{
  int i;

  /* Locate residue */
  for(i=0; (i<natoms) && (atom[i].resnr != resi); i++)
    ;
  if (i == natoms)
    return -1;
    
  /* Compare atom names */
  for(   ; (i<natoms) && (atom[i].resnr == resi); i++)
    if (strcmp(*(aname[i]),ai) == 0)
      return i;
      
  /* Not found?! */
  if (bVerbose)
    fprintf(stderr,"Warning: atom %s not found in res %s%d (line %d)\n",
	    ai,ri ? ri : "",resi+r0,linec);
  
  return -1;
}

void conv_dr(FILE *in,FILE *out,char *map,t_atoms *atoms,int r0,gmx_bool bXplor,
	     gmx_bool bVerbose)
{
  static char *format="%s%d%s%s%d%s%lf%lf";
  static char *xplorformat="%d%s%d%s";
  gmx_bool   bOK;
  int    i,j,nc,nindex,ni,nj,nunres;
  int    atomi,atomj,resi,resj;
  char   **aiexp,**ajexp;
  char   *ai,*aj;
  char   *ri,*rj;
  char   buf[1024];
  double ub,lb;
  int    linec;
  int    NEXP;
  t_expansion *exp;
  
  exp=read_expansion_map(map,&NEXP);
  
  nc=0;
  nindex=0;
  nunres=0;
  snew(ai,10);
  snew(aj,10);
  fprintf(out,"[ distance_restraints ]\n");
  linec=1;
  
  if (bXplor) {
    ri = rj = NULL;
  }
  else {
    snew(ri,16);
    snew(rj,16);
  }
  while (fgets2(buf,1023,in) != NULL) {
    /* Parse the input string. If your file format is different,
     * modify it here...
     * If your file contains no spaces but colon (:) for separators
     * it may be easier to just replace those by a space using a
     * text editor.
     */
    if (bXplor) {
      bOK = (sscanf(buf,xplorformat,&resi,ai,&resj,aj) == 4);
      /* Cut atomnames at 4 characters */
      if (strlen(ai) >= 4)
	ai[4] = '\0';
      if (strlen(aj) >= 4)
	aj[4] = '\0';
      ub = 5.0;
      lb = 2.0;
    }
    else {
      bOK = (sscanf(buf,format,ri,&resi,ai,rj,&resj,aj,&lb,&ub) == 8);
    }
    if (bOK) {
      aiexp=get_exp(NEXP,exp,&ai,&ni);
      ajexp=get_exp(NEXP,exp,&aj,&nj);
      
      /* Turn bounds into nm */
      ub*=0.1;
      lb*=0.1;
      
      /* Subtract starting residue to match topology */
      resi-=r0;
      resj-=r0;
      
      /* Test whether residue names match 
       * Again, if there is a mismatch between GROMACS names
       * and your file (eg. HIS vs. HISH) it may be easiest to
       * use your text editor...
       */
       
      if (!bXplor) {
	bOK = (strcmp(*atoms->resname[resi],ri) == 0);
	if (!bOK) {
	  fprintf(stderr,"Warning resname in disres file %s%d, in tpx file %s%d\n",
		  ri,resi+r0,*atoms->resname[resi],resi+r0);
	  nunres++;
	}
	else {
	  /* Residue j */
	  bOK = (strcmp(*atoms->resname[resj],rj) != 0);
	  if (!bOK) {
	    fprintf(stderr,"Warning resname in disres file %s%d, in tpx file %s%d\n",
		    rj,resj+r0,*atoms->resname[resj],resj+r0);
	    nunres++;
	  }
	}
      }
      if (bOK) {
	/* Here, both residue names match ! */
	for(i=0; (i<ni); i++) {
	  if ((atomi=find_atom(aiexp[i],ri,resi,r0,atoms->nr,
			       atoms->atomname,atoms->atom,linec,bVerbose)) == -1)
	    nunres++;
	  else {
	    /* First atom is found now... */
	    for(j=0; (j<nj); j++) {
	      if ((atomj=find_atom(ajexp[j],rj,resj,r0,atoms->nr,
				   atoms->atomname,atoms->atom,linec,bVerbose)) == -1)
		nunres++;
	      else {
		/* BOTH atoms are found now! */
		fprintf(out,"%5d  %5d  1  %5d  1  %8.3f  %8.3f  %8.3f  %8.3f\n",
			1+atomi,1+atomj,nindex,lb,ub,0.0,0.0);
		nc++;
	      }
	    }
	  }
	}
      }
      nindex++;
    }
    linec++;
  }
  fprintf(stderr,"Total number of NOES: %d\n",nindex);
  fprintf(stderr,"Total number of restraints: %d\n",nc);
  fprintf(stderr,"Total number of unresolved atoms: %d\n",nunres);
  if (nunres+nc != nindex) 
    fprintf(stderr,"Holy Cow! some lines have disappeared.\n");
}

int main (int argc,char *argv[])
{
  const char *desc[] = {
    "gendr generates a distance restraint entry for a gromacs topology",
    "from another format. The format of the input file must be:[BR]",
    "resnr-i resname-i atomnm-i resnr-j resname-j atomnm-j lower upper[BR]"  ,
    "where lower and upper are the distance bounds.",
    "The entries must be separated by spaces, but may be otherwise in",
    "free format. Some expansion of templates like MB -> HB1, HB2 is done",
    "but this is not really well tested."
  };
  const char *bugs[] = {
    "This program is not well tested. Use at your own risk."
  };
  
  static int  r0       = 1;
  static gmx_bool bXplor   = FALSE;
  static gmx_bool bVerbose = FALSE;
  t_pargs pa[] = {
    { "-r",     FALSE, etINT,  {&r0},       "starting residue number" },
    { "-xplor", FALSE, etBOOL, {&bXplor},   "Use xplor format for input" },
    { "-v",     FALSE, etBOOL, {&bVerbose}, "Be loud and noisy" }
  };

  FILE        *in,*out;
  t_topology  *top;
  
  t_filenm fnm[] = {
    { efTPX, "-s", NULL, ffREAD  },
    { efDAT, "-d", NULL, ffREAD  },
    { efITP, "-o", NULL, ffWRITE },
    { efDAT, "-m", "expmap", ffREAD }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,asize(desc),desc,
		    asize(bugs),bugs);

  fprintf(stderr,"******************* WARNING *****************************\n");
  fprintf(stderr,"*** Use at your own risk. When in doubt check the source.\n");
  fprintf(stderr,"*** Hang on: check the source anyway.\n");
  fprintf(stderr,"******************* WARNING *****************************\n");
		    
  fprintf(stderr,"Will subtract %d from res numbers in %s\n",
	  r0,ftp2fn(efDAT,NFILE,fnm));
    
  top=read_top(ftp2fn(efTPX,NFILE,fnm));

  in  = opt2FILE("-d",NFILE,fnm,"r");
  out = ftp2FILE(efITP,NFILE,fnm,"w");
  conv_dr(in,out,opt2fn("-m",NFILE,fnm),&(top->atoms),r0,bXplor,bVerbose);
  ffclose(in);
  ffclose(out);
  
  thanx(stderr);
  
  return 0;
}


