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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sysstuff.h"
#include <ctype.h>
#include "string2.h"
#include "strdb.h"
#include "futil.h"
#include "smalloc.h"
#include "gmx_fatal.h"
#include "symtab.h"
#include "macros.h"
#include "resall.h"
#include "pgutil.h"

t_atomtype *read_atype(char *adb,t_symtab *tab)
{
#define MAXAT 5000
  FILE       *in;
  char       aadb[STRLEN];
  char       buf[STRLEN],name[STRLEN];
  double     m;
  int        nratt;
  t_atomtype *at;

  sprintf(aadb,"%s.atp",adb);
  in=libopen(aadb);
  
  snew(at,1);
  snew(at->atom,MAXAT);
  snew(at->atomname,MAXAT);
  
  for(nratt=0; ; nratt++) {
    if (nratt >= MAXAT)
      gmx_fatal(FARGS,"nratt >= MAXAT(%d). Increase the latter",MAXAT);
    if (feof(in))
      break;

    /* Skip blank or comment-only lines */
    do {
      fgets2(buf,STRLEN,in);
      if(buf) {
	strip_comment(buf);
	trim(buf);
      }
    } while (buf && strlen(buf)==0);

    if(buf==NULL)
      break;
    
    if (sscanf(buf,"%s%lf",name,&m) != 2)
      break;
    set_at(&(at->atom[nratt]),m,0.0,nratt,0);
    at->atomname[nratt]=put_symtab(tab,name);
    fprintf(stderr,"\rAtomtype %d",nratt+1);
  }
  fclose(in);
  fprintf(stderr,"\n");
  at->nr=nratt;
  
  return at;
}

static void print_resatoms(FILE *out,t_atomtype *atype,t_restp *rtp)
{
  int j,tp;
  
  /* fprintf(out,"%5s\n",rtp->resname);
     fprintf(out,"%5d\n",rtp->natom); */
  fprintf(out,"[ %s ]\n",rtp->resname);
  fprintf(out," [ atoms ]\n");
  
  for(j=0; (j<rtp->natom); j++) {
    tp=rtp->atom[j].type;
    if ((tp < 0) || (tp >= atype->nr))
      gmx_fatal(FARGS,"tp (%d) out of range (0 .. %d)",tp,atype->nr);
    fprintf(out,"%6s%6s%8.3f%6d\n",
	    *(rtp->atomname[j]),*(atype->atomname[tp]),
	    rtp->atom[j].q,rtp->cgnr[j]);
  }
}

static bool read_atoms(FILE *in,char *line,
		       t_restp *r0,t_symtab *tab,t_atomtype *atype)
{
  int    i,j,cg,maxentries;
  char   buf[256],buf1[256];
  double q;

  /* Read Atoms */
  maxentries=0;
  r0->atom=     NULL;
  r0->atomname= NULL;
  r0->cgnr=     NULL;
  i=0;
  while (get_a_line(in,line,STRLEN) && (strchr(line,'[')==NULL)) { 
    if (sscanf(line,"%s%s%lf%d",buf,buf1,&q,&cg) != 4)
      return FALSE;
    if (i>=maxentries) {
      maxentries+=100;
      srenew(r0->atom,     maxentries);
      srenew(r0->atomname, maxentries);
      srenew(r0->cgnr,     maxentries);
    }
    r0->atomname[i]=put_symtab(tab,buf);
    r0->atom[i].q=q;
    r0->cgnr[i]=cg;
    for(j=0; (j<atype->nr); j++)
      if (strcasecmp(buf1,*(atype->atomname[j])) == 0)
	break;
    if (j == atype->nr)
      gmx_fatal(FARGS,"Atom type %s (residue %s) not found in atomtype "
		  "database",buf1,r0->resname);
    r0->atom[i].type=j;
    r0->atom[i].m=atype->atom[j].m;
    i++;
  }
  r0->natom=i;
  srenew(r0->atom,i);
  srenew(r0->atomname,i);
  srenew(r0->cgnr,i);

  return TRUE;
}

bool read_bondeds(int bt, FILE *in, char *line, t_restp *rtp)
{
  char str[STRLEN];
  int  j,n,ni,maxrb;
  
  maxrb = rtp->rb[bt].nb;
  while (get_a_line(in,line,STRLEN) && (strchr(line,'[')==NULL)) {
    if ( rtp->rb[bt].nb >= maxrb ) {
      maxrb+=100;
      srenew(rtp->rb[bt].b,maxrb);
    }
    n=0;
    for(j=0; j < btsNiatoms[bt]; j++) {
      if ( sscanf(line+n,"%s%n",str,&ni)==1 )
	rtp->rb[bt].b[rtp->rb[bt].nb].a[j]=strdup(str);
      else
	return FALSE;
      n+=ni;
    }
    for(  ; j < MAXATOMLIST; j++)
      rtp->rb[bt].b[rtp->rb[bt].nb].a[j]=NULL;
    while (isspace(line[n]))
      n++;
    rtrim(line+n);
    rtp->rb[bt].b[rtp->rb[bt].nb].s=strdup(line+n);
    rtp->rb[bt].nb++;
  }
  /* give back unused memory */
  srenew(rtp->rb[bt].b,rtp->rb[bt].nb);
  
  return TRUE;
}

static void print_resbondeds(FILE *out, int bt, t_restp *rtp)
{
  int i,j;
  
  if (rtp->rb[bt].nb) {
    fprintf(out," [ %s ]\n",btsNames[bt]);
    
    for(i=0; i < rtp->rb[bt].nb; i++) {
      for(j=0; j<btsNiatoms[bt]; j++)
	fprintf(out,"%6s ",rtp->rb[bt].b[i].a[j]);
      if (rtp->rb[bt].b[i].s[0])
	fprintf(out,"    %s",rtp->rb[bt].b[i].s);
      fprintf(out,"\n");
    }
  }
}

static void check_rtp(int nrtp,t_restp rtp[],char *libfn)
{
  int i;

  /* check for double entries, assuming list is already sorted */
  for(i=1; (i<nrtp); i++) {
    if (strcasecmp(rtp[i-1].resname,rtp[i].resname) == 0)
      fprintf(stderr,"WARNING double entry %s in file %s\n",
	      rtp[i].resname,libfn);
  }
}

static int comprtp(const void *a,const void *b)
{
  t_restp *ra,*rb;

  ra=(t_restp *)a;
  rb=(t_restp *)b;

  return strcasecmp(ra->resname,rb->resname);
}

int get_bt(char* header)
{
  int i;

  for(i=0; i<ebtsNR; i++)
    if ( strcasecmp(btsNames[i],header)==0 )
      return i;
  return NOTSET;
}

void clear_t_restp(t_restp *rrtp)
{
  memset((void *)rrtp, 0, sizeof(t_restp));
}

int read_resall(char *ff, int bts[], t_restp **rtp, 
		t_atomtype *atype, t_symtab *tab, bool *bAlldih, int *nrexcl,
		bool *HH14, bool *bRemoveDih)
{
  FILE      *in;
  char      rrdb[STRLEN],line[STRLEN],header[STRLEN];
  int       i,nrtp,maxrtp,bt,nparam;
  t_restp   *rrtp;
  bool      bNextResidue,bError;
  
  sprintf(rrdb,"%s.rtp",ff);
  in=libopen(rrdb);
  rrtp=NULL;
  if (debug) {
    fprintf(debug,"%9s %5s", "Residue", "atoms");
    for(i=0; i<ebtsNR; i++)
      fprintf(debug," %10s",btsNames[i]);
    fprintf(debug,"\n");
  }

  /* these bonded parameters will overwritten be when  *
   * there is a [ bondedtypes ] entry in the .rtp file */
  bts[0] = 1; /* normal bonds     */
  bts[1] = 1; /* normal angles    */
  bts[2] = 1; /* normal dihedrals */
  bts[3] = 2; /* normal impropers */
  
  
  /* Column 5 & 6 aren't really bonded types, but we include
   * them here to avoid introducing a new section:
   * Column 5: 1 means generate all dihedrals, 0 not.
   * Column 6: Number of bonded neighbors to exclude.
   * Coulmn 7: Generate 1,4 interactions between pairs of hydrogens
   * Column 8: Remove impropers over the same bond as a proper dihedral
   */
  
  nrtp=0;
  maxrtp=0;
  get_a_line(in,line,STRLEN);
  if (!get_header(line,header))
    gmx_fatal(FARGS,"in .rtp file at line:\n%s\n",line);
  if (strncasecmp("bondedtypes",header,5)==0) {
    get_a_line(in,line,STRLEN);
    if ((nparam=sscanf(line,"%d %d %d %d %d %d %d %d",&bts[0],&bts[1],&bts[2],&bts[3],bAlldih,nrexcl,HH14,bRemoveDih)) < 4 )
      gmx_fatal(FARGS,"need at least 4 (up to 8) parameters in .rtp file at line:\n%s\n",line);
    get_a_line(in,line,STRLEN);
    if(nparam<5) {
      fprintf(stderr,"Using default: not generating all possible dihedrals\n");
      *bAlldih=FALSE;
    }
    if(nparam<6) {
      fprintf(stderr,"Using default: excluding 3 bonded neighbors\n");
      *nrexcl=3;
    }
    if(nparam<7) {
      fprintf(stderr,"Using default: generating 1,4 H--H interactions\n");
      *HH14=TRUE;
    }
    if(nparam<8) {
      fprintf(stderr,"Using default: removing impropers on same bond as a proper\n");
      *bRemoveDih=TRUE;
    }
  } else {
    fprintf(stderr,
	    "Reading .rtp file without '[ bondedtypes ]' directive,\n"
	    "Will proceed as if the entry\n"
	    "\n"
	    "\n[ bondedtypes ]"
	    "\n; bonds  angles  dihedrals  impropers all_dihedrals nr_exclusions HH14 remove_dih"
	    "\n   %3d     %3d        %3d        %3d           %3d           %3d   %3d    %3d"
	    "\n"
	    "was present at the beginning of %s",
	    bts[0],bts[1],bts[2],bts[3], (*bAlldih) ? 1 : 0,*nrexcl,*HH14,*bRemoveDih,rrdb);
  }
  nrtp=0;
  while (!feof(in)) {
    if (nrtp >= maxrtp) {
      maxrtp+=100;
      srenew(rrtp,maxrtp);
    }
    clear_t_restp(&rrtp[nrtp]);
    if (!get_header(line,header))
      gmx_fatal(FARGS,"in .rtp file at line:\n%s\n",line);
    rrtp[nrtp].resname=strdup(header);
    
    get_a_line(in,line,STRLEN);
    bError=FALSE;
    bNextResidue=FALSE;
    do {
      if (!get_header(line,header)) {
	bError = TRUE;
      } else {
	bt = get_bt(header);
	if (bt != NOTSET) {
	  /* header is an bonded directive */
	  bError = !read_bondeds(bt,in,line,&rrtp[nrtp]);
	} else if (strncasecmp("atoms",header,5) == 0) {
	  /* header is the atoms directive */
	  bError = !read_atoms(in,line,&(rrtp[nrtp]),tab,atype);
	} else {
	  /* else header must be a residue name */
	  bNextResidue = TRUE;
	}
      }
      if (bError)
	gmx_fatal(FARGS,"in .rtp file in residue %s at line:\n%s\n",
		    rrtp[nrtp].resname,line);
    } while (!feof(in) && !bNextResidue);

    if (rrtp[nrtp].natom == 0)
      gmx_fatal(FARGS,"No atoms found in .rtp file in residue %s\n",
		  rrtp[nrtp].resname);
    if (debug) {
      fprintf(debug,"%3d %5s %5d",
	      nrtp+1,rrtp[nrtp].resname,rrtp[nrtp].natom);
      for(i=0; i<ebtsNR; i++)
	fprintf(debug," %10d",rrtp[nrtp].rb[i].nb);
      fprintf(debug,"\n");
    }
    nrtp++;
    fprintf(stderr,"\rResidue %d",nrtp);
  }
  fclose(in);
  /* give back unused memory */
  srenew(rrtp,nrtp);
  
  fprintf(stderr,"\nSorting it all out...\n");
  qsort(rrtp,nrtp,(size_t)sizeof(rrtp[0]),comprtp);
  
  check_rtp(nrtp,rrtp,rrdb);
  
  *rtp  = rrtp;
  
  return nrtp;
}

void print_resall(FILE *out, int bts[], int nrtp, t_restp rtp[],
		  t_atomtype *atype, bool bAlldih, int nrexcl, 
		  bool HH14, bool bRemoveDih)
{
  int i,bt;

  /* print all the ebtsNR type numbers */
  fprintf(out,"[ bondedtypes ]\n");
  fprintf(out,"; bonds  angles  dihedrals  impropers all_dihedrals nr_exclusions  HH14  remove_dih\n");
  fprintf(out," %5d  %6d  %9d  %9d  %14d  %14d %14d %14d\n\n",bts[0],bts[1],bts[2],bts[3],bAlldih,nrexcl,HH14,bRemoveDih);

  for(i=0; i<nrtp; i++) {
    if (rtp[i].natom > 0) {
      print_resatoms(out,atype,&rtp[i]);
      for(bt=0; bt<ebtsNR; bt++)
	print_resbondeds(out,bt,&rtp[i]);
    }
  }
}

/************************************************************
 *
 *                  SEARCH   ROUTINES
 * 
 ***********************************************************/
int neq_str(char *a1,char *a2)
{
  int j,l;
  
  l=min((int)strlen(a1),(int)strlen(a2));
  j=0;
  while ( (j<l) && (toupper(a1[j]) == toupper(a2[j])) )
    j++;
  
  return j;
}

t_restp *search_rtp(char *key,int nrtp,t_restp rtp[])
{
  int i,n,best,besti;

  besti=-1;
  best=1;
  for(i=0; (i<nrtp); i++) {
    n=neq_str(key,rtp[i].resname);
    if (n > best) {
      besti=i;
      best=n;
    }
  }
  if (besti == -1)
    gmx_fatal(FARGS,"Residue '%s' not found in residue topology database\n",key);
  if (strlen(rtp[besti].resname) != strlen(key))
    fprintf(stderr,"Warning: '%s' not found in residue topology database, "
	    "trying to use '%s'\n", key, rtp[besti].resname);
  
  return &rtp[besti];
}
