/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
static char *SRCID_pdb2gmx_c = "$Id$";
#include <time.h>
#include <ctype.h>
#include "assert.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "copyrite.h"
#include "string2.h"
#include "confio.h"
#include "symtab.h"
#include "vec.h"
#include "statutil.h"
#include "futil.h"
#include "fatal.h"
#include "pdbio.h"
#include "toputil.h"
#include "h_db.h"
#include "physics.h"
#include "pgutil.h"
#include "calch.h"
#include "resall.h"
#include "pdb2top.h"
#include "ter_db.h"
#include "strdb.h"
#include "txtdump.h"
#include "gbutil.h"
#include "genhydro.h"
#include "readinp.h"
#include "xlate.h"
#include "specbond.h"
#include "index.h"
#include "hizzie.h"

static void dump_res(FILE *fp,int indent,t_restp *restp,int nah,
		     t_hackblock *ah)
{
  char *rtype[ebtsNR] = { "rbond", "rangle", "rdihedral", "rimproper" };
  int ntype[ebtsNR] = { 2, 3, 4, 4 };
  int  i,j,k;
  
  pr_indent(fp,indent);
  fprintf(fp,"<residue restype=\"%s\">\n",restp->resname);
  indent += 2;

  for(i=0; (i<restp->natom); i++) {
    pr_indent(fp,indent);
    fprintf(fp,"<ratom name=\"%s\">\n",*restp->atomname[i]);
  }
  for(i=0; (i<ebtsNR); i++) {
    for(j=0; (j<restp->rb[i].nb); j++) {
      pr_indent(fp,indent);
      fprintf(fp,"<%s",rtype[i]);
      for(k=0; (k<ntype[i]); k++)
	fprintf(fp," a%d=\"%s\"",k+1,restp->rb[i].b[j].a[k]);
      fprintf(fp,"/>\n");
    }
  }
  for(i=0; (i<nah); i++) {
    if (strcmp(restp->resname,ah[i].name) == 0) {
      for(j=0; (j<ah[i].nhack); j++) {
	pr_indent(fp,indent);
	if ((ah[i].hack[j].a[0][0] == 'O') || (ah[i].hack[j].a[0][0] == 'N'))
	  fprintf(fp,"<raddh hclass=\"polar\"    ");
	else
	  fprintf(fp,"<raddh hclass=\"aliphatic\"");
	fprintf(fp," addgeom=\"%d\" addnum=\"%d\"",
		ah[i].hack[j].tp,ah[i].hack[j].nr);
	fprintf(fp," addto=\"%s %s %s",ah[i].hack[j].a[0],ah[i].hack[j].a[1],
		ah[i].hack[j].a[2]);
	if (ah[i].hack[j].a[3] > 0)
	  fprintf(fp," %s",ah[i].hack[j].a[3]);
	fprintf(fp,"\"/>\n");
      }
      break;
    }
  }
  indent -= 2;
  pr_indent(fp,indent);
  fprintf(fp,"</residue>\n");
}

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "This program reads an rtp file and dumps an xml file."
  };

  typedef struct {
    char chain;
    int  start;
    int  natom;
    bool bAllWat;
    int  nterpairs;
    int  *chainstart;
  } t_pdbchain;

  typedef struct {
    char chain;
    bool bAllWat;
    int nterpairs;
    int *chainstart;
    t_hackblock **ntdb;
    t_hackblock **ctdb;
    int *rN;
    int *rC;
    t_atoms *pdba;
    rvec *x;
  } t_chain;
  
  FILE       *fp;
  int        natom,nres;
  t_atoms    pdba_all,*pdba;
  t_atoms    *atoms;
  t_block    *block;
  int        chain,nch,maxch,nwaterchain;
  t_pdbchain *pdb_ch;
  t_chain    *chains,*cc;
  char       pchain,select[STRLEN];
  int        nincl,nmol;
  char       **incls;
  t_mols     *mols;
  char       **gnames;
  matrix     box;
  rvec       box_space;
  char       *ff;
  int        i,j,k,l,nrtp;
  int        *swap_index,si;
  int        bts[ebtsNR];
  t_restp    *restp;
  t_hackblock *ah;
  t_symtab   symtab;
  t_atomtype *atype;
  char       fn[256],*top_fn,itp_fn[STRLEN],posre_fn[STRLEN],buf_fn[STRLEN];
  char       molname[STRLEN],title[STRLEN];
  char       *c;
  int        nah,nNtdb,nCtdb;
  t_hackblock *ntdb,*ctdb;
  int        nssbonds;
  t_ssbond   *ssbonds;
  rvec       *pdbx,*x;
  bool       bUsed,bDummies=FALSE,bWat,bPrevWat=FALSE,bITP,bDummyAromatics=FALSE;
  real       mHmult=0;
  
  CopyRight(stderr,argv[0]);
	
  ff = strdup("ffgmx2");
  
  /* Read atomtypes... */
  atype=read_atype(ff,&symtab);
    
  /* read residue database */
  printf("Reading residue database... (%s)\n",ff);
  nrtp=read_resall(ff,bts,&restp,atype,&symtab);
  
  /* read hydrogen database */
  nah=read_h_db(ff,&ah);
  
  /* Read Termini database... */
  nNtdb=read_ter_db(ff,'n',&ntdb,atype);
  nCtdb=read_ter_db(ff,'c',&ctdb,atype);
  
  fp=fopen("residues.xml","w");
  fprintf(fp,"<?xml version=\"1.0\"?>\n");
  fprintf(fp,"<!DOCTYPE residues SYSTEM \"residues.dtd\">\n");
  fprintf(fp,"\n<residues>\n");
  for(i=0; (i<nrtp); i++) {
    dump_res(fp,2,&(restp[i]),nah,ah);
  }
  fprintf(fp,"</residues>\n");
  fclose(fp);  

  thanx(stderr);
  
  return 0;
}
