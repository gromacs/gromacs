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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <time.h>
#include <ctype.h>
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
#include "gmx_fatal.h"
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
		     t_hackblock *ah,int nres_long,char **res_long)
{
  char *rtype[ebtsNR] = { "rbond", "rangle", "rdihedral", "rimproper" };
  int ntype[ebtsNR] = { 2, 3, 4, 4 };
  int  i,j,k,nn;
  char *tmp,descr[128];
  
  descr[0] = '\0';
  for(nn=0; (nn < nres_long) && (strstr(res_long[nn],restp->resname) == NULL); nn++) 
    ;
  if (nn < nres_long) {
    tmp = res_long[nn] + strlen(restp->resname);
    while (*tmp && isspace(*tmp)) {
      tmp++;
    }
    if (strlen(tmp) > 0)
      strcpy(descr,tmp);
  }
  
  pr_indent(fp,indent);
  fprintf(fp,"<residue restype=\"%s\" longname=\"%s\">\n",restp->resname,descr);
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
	fprintf(fp," addto=\"%s",ah[i].hack[j].a[0]);
	for(k=1; ((k <= 3) && (ah[i].hack[j].a[k] > 0)); k++)
	  fprintf(fp," %s",ah[i].hack[j].a[k]);
	fprintf(fp,"\"/>\n");
      }
      break;
    }
  }
  indent -= 2;
  pr_indent(fp,indent);
  fprintf(fp,"</residue>\n\n");
}

static void dump_hack_add(FILE *fp,int indent,t_hack *hack)
{
  pr_indent(fp,indent);
  fprintf(fp,"<modadd addname=\"%s\" addgeom=\"%d\" addto=\"%s %s %s\"/>\n",
	  hack->nname,hack->tp,hack->a[0],hack->a[1],hack->a[2]);
}

static void dump_hack_del(FILE *fp,int indent,t_hack *hack)
{
  pr_indent(fp,indent);
  fprintf(fp,"<moddelete delname=\"%s\"/>\n",hack->oname);
}

static void dump_hack_rep(FILE *fp,int indent,t_hack *hack)
{
  pr_indent(fp,indent);
  fprintf(fp,"<modreplace oldname=\"%s\" newname=\"%s\"/>\n",hack->oname,hack->nname);
}

static void dump_mod(FILE *fp,int indent,t_hackblock *tdb,t_atomtype *atype)
{
  int i,j,k;

  pr_indent(fp,indent);
  fprintf(fp,"<moddef modtype=\"%s\"\n",tdb->name);
  indent += 2;
  for(j=0; (j<tdb->nhack); j++) {
    if (tdb->hack[j].oname == NULL)
      dump_hack_add(fp,indent,&(tdb->hack[j]));
    else if (tdb->hack[j].nname == NULL)
      dump_hack_del(fp,indent,&(tdb->hack[j]));
    else 
      dump_hack_rep(fp,indent,&(tdb->hack[j]));
  }
  indent -= 2;
  pr_indent(fp,indent);
  fprintf(fp,"</moddef>\n\n");
}

static void dump_mods(FILE *fp,int indent,
 		      int n,t_hackblock tdb[],t_atomtype *atype)
{
  int i,j;
  
  for(i=0; (i<n); i++)
    dump_mod(fp,indent,&(tdb[i]),atype);
}

static void do_specbonds(FILE *fp,int indent)
{
  t_specbond *sb;
  int i,nsb;
  
  pr_indent(fp,indent);
  fprintf(fp,"<linkdef linktype=\"peptide\" restype=\"* *\" atomprev=\"C\" atomnext=\"N\" refdist=\"0.133\"/>\n");
  sb = get_specbonds(&nsb);
  for(i=0; (i<nsb); i++) {
    pr_indent(fp,indent);
    fprintf(fp,"<linkdef linktype=\"unknown\" restype=\"%s %s\" atomprev=\"%s\" atomnext=\"%s\" refdist=\"%g\"/>\n",
	    sb[i].res1,sb[i].res2,sb[i].atom1,sb[i].atom2,sb[i].length);
  }
  done_specbonds(nsb,sb);
  sfree(sb);
  fprintf(fp,"\n");
}

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "This program reads an [TT].rtp[tt] file and dumps an [TT].xml[tt] file."
  };

  typedef struct {
    char chain;
    int  start;
    int  natom;
    gmx_bool bAllWat;
    int  nterpairs;
    int  *chainstart;
  } t_pdbchain;

  typedef struct {
    char chain;
    gmx_bool bAllWat;
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
  char       **res_long;
  int        nres_long;
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
  gmx_bool       bUsed,bDummies=FALSE,bWat,bPrevWat=FALSE,bITP,bDummyAromatics=FALSE;
  real       mHmult=0;
  int        nrexcl;
  gmx_bool       bAlldih,bH14,bRemoveDih;

  CopyRight(stderr,argv[0]);
	
  ff = strdup("ffgmx2");

  /* Read long residue names */
  nres_long = get_file("res-long.dat",&res_long);
    
  /* Read atomtypes... */
  atype=read_atype(ff,&symtab);
    
  /* read residue database */
  printf("Reading residue database... (%s)\n",ff);
  nrtp=read_resall(ff,bts,&restp,atype,&symtab,&bAlldih,&nrexcl,&bH14,&bRemoveDih);
  
  /* read hydrogen database */
  nah=read_h_db(ff,&ah);
  
  /* Read Termini database... */
  nNtdb=read_ter_db(ff,'n',&ntdb,atype);
  nCtdb=read_ter_db(ff,'c',&ctdb,atype);
  
  fp=gmx_fio_fopen("residues.xml","w");
  /* fprintf(fp,"<?xml version=\"1.0\"?>\n");
     fprintf(fp,"<!DOCTYPE residues SYSTEM \"residues.dtd\">\n"); */
  fprintf(fp,"\n<residues>\n");
  for(i=0; (i<nrtp); i++) {
    dump_res(fp,2,&(restp[i]),nah,ah,nres_long,res_long);
  }
  do_specbonds(fp,2);
  dump_mods(fp,2,nNtdb,ntdb,atype);
  dump_mods(fp,2,nCtdb,ctdb,atype);
  fprintf(fp,"</residues>\n");
  gmx_fio_fclose(fp);  

  thanx(stderr);
  
  return 0;
}
