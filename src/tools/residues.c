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
#include <math.h>

#include "cdist.h"
#include "names.h"

/* Some Macros that apply in each routine */
#define MAXLOGG asize(sa)
#define NPD     asize(pd)
#define NGAU    asize(gau)
#define M_PI43  (4*M_PI/3)
#define M_PI23  (2*M_PI/3) 

enum { CIS, TRANS, GAUCHE };

typedef struct {
  char anm[12];
  int  anr;
} t_cdatom;

typedef struct {
  char *ai,*aj,*ak,*al;
  int  type;
} t_cdpdih;

typedef struct {
  char *ai,*aj,*ak,*al,*am;
  int  type;
} t_cd15;

typedef struct {
  char *ai,*aj,*ak,*al,*am;
  int  type,bCis;
} t_cd15a;

typedef struct {
  char *ai,*aj,*ak,*al,*am;
  real om1,om2,om3;
} t_cd15g;

typedef struct {
  char *ai,*aj,*ak,*al,*am,*an;
  int  type;
} t_cd16;

/* Debug output */
static void pr_logg(FILE *fp,int nlogg,t_cdatom cda[],bool bVir,char *label)
{
  int i;
  
  if (fp) {
    fprintf(fp,"%5s  virtual = %3s",label,yesno_names[bVir]);
    if (bVir)
      nlogg++;
    for(i=0; (i<nlogg); i++)
      fprintf(fp,"  %s:%d",cda[i].anm,cda[i].anr);
    fprintf(fp,"\n");
  }
}

static void pr_ndist(FILE *fp,char *res,int n14,int n15,int n16,int n17,int nV)
{
  static bool bFirst = TRUE;
  int    ntot;
  
  if ((ntot = n14+n15+n16+n17+nV) == 0)
    return;
  if (bFirst) {
    fprintf(fp,"Res     #1-4    #1-5    #1-6    #1-7   #Virt  #Total\n");
    bFirst = FALSE;
  }
  fprintf(fp,"%-4s%8d%8d%8d%8d%8d%8d\n",res,n14,n15,n16,n17,nV,ntot);
}

static t_cdatom *init_cdatoms(int nsa,char *sa[])
{
  int      k;
  t_cdatom *cda;
  
  snew(cda,nsa);
  for(k=0; (k<nsa); k++) {
    strcpy(cda[k].anm,sa[k]);
    cda[k].anr = -1;
  }
  return cda;
}

static int set_cdatoms(t_atoms *atoms,int j,int residnr,int ncda,t_cdatom cda[])
{
  int i,k,nlogg=0,resj,len=0;
  char *anmj;
  bool bNext;

  resj = atoms->atom[j].resnr;
  anmj = *atoms->atomname[j];
  
  for(i=0; (i<ncda); i++) 
    cda[i].anr=-1;
  while ((resj <= residnr+1) && (j<atoms->nr)) {
    k = -1;
    for(i=0; ((i<ncda) && (k == -1)); i++) {
      bNext = (strchr(cda[i].anm,'+') != NULL);
      if (bNext)
	len = strlen(cda[i].anm)-1;

      if ((bNext  && (resj == residnr+1) && (strlen(anmj) == len) &&
	   (strncmp(anmj,cda[i].anm,len) == 0)) ||
	  (!bNext && (resj == residnr) && 
	   (strcmp(anmj,cda[i].anm) == 0)))
	k = i;
    }
    if (k != -1) {
      if (cda[k].anr != -1)
	gmx_fatal(FARGS,"Overwriting cda entry for %s [%s], was %d now %d\n"
		    "resj=%d, residnr=%d, k=%d, ncda=%d",
		    cda[k].anm,anmj,cda[k].anr+1,j+1,resj,residnr,k,ncda);
      cda[k].anr = j;
      nlogg++;
    }
    j++;
    if (j < atoms->nr) {
      resj = atoms->atom[j].resnr;
      anmj = *atoms->atomname[j];
    }
  }
  return nlogg;
}

static int logg_index(char *anm,int ncda,t_cdatom cda[],char *file,int line)
{
  int i,li=-1;
  
  for(i=0; (i<ncda) && (li == -1); i++)
    if (strcmp(anm,cda[i].anm) == 0)
      li = cda[i].anr;
  if (li == -1) 
    gmx_fatal(FARGS,"Can not find atom %s (file %s, line %d)",anm,file,line);
  return li;
}
#define logger(sss) logg_index(sss,MAXLOGG,cda,__FILE__,__LINE__)

static int do_a_bond(int ai,int aj,t_ilist ilist[],t_iparams iparams[],
		     bool bFlag,t_atoms *atoms,real margin,
		     real weight[],t_dist *d)
{
  real blen,lb,ub;
  int  natoms = atoms->nr;
  
  blen = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  lb   = (1.0-margin)*blen;
  ub   = (1.0+margin)*blen;
  if (((weight[ai] != 0.0) || (weight[aj] != 0.0)) &&
      (!dist_set(d,natoms,ai,aj))) {
    set_dist(d,natoms,ai,aj,lb,ub,blen);
    return 1;
  }
  return 0;
}

static int do_a_dist(int ai,int aj,int natoms,real margin,
		     real weight[],t_dist *d,real blen)
{
  real lb,ub;
  
  lb   = (1.0-margin)*blen;
  ub   = (1.0+margin)*blen;
  if (((weight[ai] != 0.0) || (weight[aj] != 0.0)) &&
      (!dist_set(d,natoms,ai,aj))) {
    set_dist(d,natoms,ai,aj,lb,ub,blen);
    return 1;
  }
  return 0;
}

static int do_an_angle(int ai,int aj,int ak,
		       t_ilist ilist[],t_iparams iparams[],
		       bool bFlag,t_atoms *atoms,real margin,
		       real weight[],t_dist *d)
{
  real angle,blen,lb,ub;
  int  natoms = atoms->nr;
  
  angle = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  blen  = angle_length(ai,aj,ak,RAD2DEG*angle,ilist,iparams,atoms);
  lb    = (1.0-margin)*blen;
  ub    = (1.0+margin)*blen;
  if (((weight[ai] != 0.0) || (weight[ak] != 0.0)) &&
      (!dist_set(d,natoms,ai,ak))) {
    set_dist(d,natoms,ai,ak,lb,ub,blen);
    return 1;
  }
  return 0;
}

static int do_a_pdih_(int ai,int aj,int ak,int al,int type,
		      t_ilist ilist[],t_iparams iparams[],
		      bool bFlag,t_atoms *atoms,real margin,
		      real weight[],t_dist *d,char *file,int line)
{
  real angle,blen,lb,ub;
  int  natoms = atoms->nr;
  
  if (type == GAUCHE)
    gauche_(ai,aj,ak,al,ilist,iparams,&blen,atoms,file,line);
  else {
    pdih_lengths_(ai,aj,ak,al,ilist,iparams,&lb,&ub,atoms,file,line);
    if (type == CIS)
      blen = lb;
    else
      blen = ub;
  }
  lb=(1.0-margin)*blen;
  ub=(1.0+margin)*blen;
  if (((weight[ai] != 0.0) || (weight[al] != 0.0)) &&
      (!dist_set(d,natoms,ai,al))) {
    set_dist(d,natoms,ai,al,lb,ub,blen);
    return 1;
  }
  return 0;
}
#define do_a_pdih(ai,aj,ak,al,type,ilist,iparams,bFlag,atoms,margin,weight,d)\
  do_a_pdih_(ai,aj,ak,al,type,ilist,iparams,bFlag,atoms,margin,weight,d,__FILE__,__LINE__)

int set_virtual (t_cdatom cda[],int N,real margin,t_dist *d,int natoms)
{
  /* Routine to add distances to virtual particle. 
     The virtual particle is placed 10A outside
     the plane common to all atoms, straight above
     the center of mass, which is calculated assuming 
     unit mass for all atoms.

     Adam Kirrander 990211                        */
  
  int ndist=0,i,j;
  real CONST=0.0,Ki,tmp,len,lb,ub;
  
  /* Calculate the constant part */
  for (i=1 ; i<N ; i++ ) {
    for (j=0 ; j<i ; j++) {
      tmp = d_len(d,natoms,cda[i].anr,cda[j].anr);
      CONST += tmp*tmp;
    }
  }
  CONST = CONST/N;
  
  /* Calculate and set distances */
  for (i=0 ; i<N ; i++ ) {
    Ki = 0.0;
    for (j=0 ; j<N ; j++) {                      /*Calc. variable part*/
      if (i == j) continue;
      tmp = d_len(d,natoms,cda[i].anr,cda[j].anr); 
      Ki += tmp*tmp;
    }
    len = sqrt(64.0+((Ki-CONST)/N));              /*Pythagoras*/
    lb  = (1.0-margin)*len;
    ub  = (1.0+margin)*len; 
    set_dist(d,natoms,cda[i].anr,cda[N].anr,lb,ub,len); 
    ndist++;
  }
  
  /* If number of virtual dist. should correspond to nr. atoms */
  if (ndist != N) 
    fprintf(stderr,"Check routine set_virtual!\n");
  
  return ndist;
} 

/**********************************************************
 *
 *     A R G A N I N E
 *
 **********************************************************/
/* Notice for both 15-routines: */
/* CD-HHX1 corresponds to lb and CD-HHX2 corresponds to ub */ 

static void arg_15_CDHH1(int ai,int aj,int ak,int al,int am,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,real *ub,t_atoms *atoms)
{
  /* Returns the length corresponding to a 1-5 distance */
  real rij,rjk,rkl,rlm,rik,rkm,rim;
  real th1,th2,th3,thikj,thmkl;
  real half_pi = M_PI*0.5;
  
  rij = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  rjk = lookup_bondlength(aj,ak,ilist,iparams,TRUE,atoms);
  rkl = lookup_bondlength(ak,al,ilist,iparams,TRUE,atoms);
  rlm = lookup_bondlength(al,am,ilist,iparams,TRUE,atoms);
  th1 = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  th2 = lookup_angle(aj,ak,al,ilist,iparams,atoms);
  th3 = lookup_angle(ak,al,am,ilist,iparams,atoms);

  /* Compute distance from i to k using law of cosines */
  rik = sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(th1));

  /* Compute distance from k to m using law of cosines */
  rkm = sqrt(rkl*rkl+rlm*rlm-2.0*rkl*rlm*cos(th3));
  
  /* Compute angle th21 using law of sines */
  thikj = asin(rij*sin(th1)/rik);
  
  /* Compute th99 using law of sines */
  thmkl = asin(rlm*sin(th3)/rkm);

  /* Compute cis length using law of cosines */
  *lb = sqrt(rik*rik+rkm*rkm-2.0*rik*rkm*cos(th2-thikj-thmkl));
  
  /* Compute trans length using law of cosines */
  *ub = sqrt(rik*rik+rkm*rkm-2.0*rik*rkm*cos(th2-thikj+thmkl));
}

static void arg_15_CDHH2(int ai,int aj,int ak,int al,int am,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,real *ub,t_atoms *atoms)
{
  /* Returns the length corresponding to a 1-5 distance */
  real rij,rjk,rkl,rlm,rik,rkm,rim;
  real th1,th2,th3,thikj,thmkl;
  real half_pi = M_PI*0.5;
  
  rij = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  rjk = lookup_bondlength(aj,ak,ilist,iparams,TRUE,atoms);
  rkl = lookup_bondlength(ak,al,ilist,iparams,TRUE,atoms);
  rlm = lookup_bondlength(al,am,ilist,iparams,TRUE,atoms);
  th1 = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  th2 = lookup_angle(aj,ak,al,ilist,iparams,atoms);
  th3 = lookup_angle(ak,al,am,ilist,iparams,atoms);

  /* Compute distance from i to k using law of cosines */
  rik = sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(th1));

  /* Compute distance from k to m using law of cosines */
  rkm = sqrt(rkl*rkl+rlm*rlm-2.0*rkl*rlm*cos(th3));
  
  /* Compute angle th99 using law of sines */
  thikj = asin(rij*sin(th1)/rik);
  
  /* Compute th21 using law of sines */
  thmkl = asin(rlm*sin(th3)/rkm);

  /* Compute cis length using law of cosines */
  *lb = sqrt(rik*rik+rkm*rkm-2.0*rik*rkm*cos(th2+thikj-thmkl));
  
  /* Compute trans length using law of cosines */
  *ub = sqrt(rik*rik+rkm*rkm-2.0*rik*rkm*cos(th2+thikj+thmkl));
}

typedef void arg_15_func(int ai,int aj,int ak,int al,int am,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,real *ub,t_atoms *atoms);

void arg (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	  real arg_margin,t_ilist ilist[],t_iparams iparams[],bool bVir)
{
  static char *sa[] = { "NE", "HE", "CZ", "NH1", "HH11", "HH12", "NH2",
			"HH21", "HH22", "CD", "VF" };
  t_cdatom *cda;
  t_cd15a  cd15a[] = {
    { "HE",   "NE",  "CZ", "NH1", "HH12", 2, 0 },
    { "HE",   "NE",  "CZ", "NH1", "HH11", 2, 1 },
    { "HE",   "NE",  "CZ", "NH2", "HH22", 1, 0 },
    { "HE",   "NE",  "CZ", "NH2", "HH21", 1, 1 },
    { "HH11", "NH1", "CZ", "NH2", "HH21", 2, 0 },
    { "HH11", "NH1", "CZ", "NH2", "HH22", 2, 1 },
    { "HH12", "NH1", "CZ", "NH2", "HH21", 1, 0 },
    { "HH12", "NH1", "CZ", "NH2", "HH22", 1, 1 },
    { "CD",   "NE",  "CZ", "NH1", "HH11", 1, 1 },
    { "CD",   "NE",  "CZ", "NH1", "HH12", 1, 0 },
    { "CD",   "NE",  "CZ", "NH2", "HH21", 2, 1 },
    { "CD",   "NE",  "CZ", "NH2", "HH22", 2, 0 }
  };
  arg_15_func *arg_15[2] = { arg_15_CDHH1, arg_15_CDHH2 };
  t_cdpdih pd[] = {
    { "NE",   "CZ",  "NH1", "HH12", TRANS }, 
    { "NE",   "CZ",  "NH2", "HH22", TRANS },
    { "HE",   "NE",  "CZ",  "NH1",  TRANS }, 
    { "HH21", "NH2", "CZ",  "NH1",  TRANS },
    { "HH11", "NH1", "CZ",  "NH2",  TRANS }, 
    { "CD",   "NE",  "CZ",  "NH2",  TRANS },
    { "CD",   "NE",  "CZ",  "NH1",  CIS }, 
    { "HH12", "NH1", "CZ",  "NH2",  CIS },
    { "HH22", "NH2", "CZ",  "NH1",  CIS }, 
    { "HE",   "NE",  "CZ",  "NH2",  CIS }, 
    { "NE",   "CZ",  "NH2", "HH21", CIS }, 
    { "NE",   "CZ",  "NH1", "HH11", CIS }
  };
  int    natoms,i,j,kk,q,residnr,oldresidnr;
  int    n14dist,n15dist,nVdist,atom_index,nlogg;
  real   blen,lb,ub,angle;

  cda        = init_cdatoms(MAXLOGG,sa);
  natoms     = atoms->nr;
  n14dist    = 0;
  n15dist    = 0;
  nVdist     = 0;
  residnr    = -1;
  oldresidnr = -1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),"ARG") == 0) {
      oldresidnr = residnr;
      residnr    = atoms->atom[i].resnr;
      if (oldresidnr == residnr)
	continue;
      
      nlogg = set_cdatoms(atoms,i,residnr,MAXLOGG,cda);

      if ( ((nlogg == MAXLOGG-1) && !bVir) || ((nlogg == MAXLOGG) && bVir) ) {

	pr_logg(debug,MAXLOGG-1,cda,bVir,"Arg");
	
	/* Proper dihedrals */
	for(kk=0; (kk<NPD); kk++)
	  n14dist += do_a_pdih(logger(pd[kk].ai),logger(pd[kk].aj),
			       logger(pd[kk].ak),logger(pd[kk].al),
			       pd[kk].type,ilist,iparams,TRUE,atoms,arg_margin,
			       weight,d);
	for(kk=0; (kk<asize(cd15a)); kk++) {
	  arg_15[cd15a[kk].type-1](logger(cd15a[kk].ai),logger(cd15a[kk].aj),
				   logger(cd15a[kk].ak),logger(cd15a[kk].al),
				   logger(cd15a[kk].am),
				   ilist,iparams,&lb,&ub,atoms);
	  n15dist += do_a_dist(logger(cd15a[kk].ai),logger(cd15a[kk].am),
			       natoms,arg_margin,
			       weight,d,cd15a[kk].bCis ? lb : ub);
	}
	/* VIRTUAL DISTANCES */
	if (bVir) 
	  nVdist += set_virtual (cda,MAXLOGG-1,arg_margin,d,natoms);

      }
    }
  }
  pr_ndist(log,"ARG",n14dist,n15dist,0,0,nVdist);
  
  sfree(cda);
}

/**********************************************************
 *
 *     A S P A R A G I N E
 *
 **********************************************************/
void asn (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	  real end_margin,t_ilist ilist[],t_iparams iparams[],bool bVir)
{
  static char *sa[] = { "CB", "CG", "OD1", "ND2", "HD21", "HD22", "VF" };
  t_cdatom *cda;
  t_cdpdih pd[] = {
    { "CB",  "CG", "ND2", "HD22", TRANS }, 
    { "OD1", "CG", "ND2", "HD21", TRANS },
    { "CB",  "CG", "ND2", "HD21", CIS }, 
    { "OD1", "CG", "ND2", "HD22", CIS }
  };
  int    natoms,i,j,kk,q,residnr,oldresidnr;
  int    n14dist,nVdist,atom_index,nlogg;
  real   blen,lb,ub,angle;

  cda        = init_cdatoms(MAXLOGG,sa);
  natoms     = atoms->nr;
  n14dist    = 0;
  nVdist     = 0;
  residnr    = -1;
  oldresidnr = -1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),"ASN") == 0) {
      oldresidnr = residnr;
      residnr    = atoms->atom[i].resnr;
      if (oldresidnr == residnr)
	continue;

      nlogg = set_cdatoms(atoms,i,residnr,MAXLOGG,cda);

      if ( ((nlogg == MAXLOGG-1) && !bVir) || ((nlogg == MAXLOGG) && bVir) ) {
	
	pr_logg(debug,MAXLOGG-1,cda,bVir,"Asn");

	/* Proper dihedrals */
	for(kk=0; (kk<NPD); kk++)
	  n14dist += do_a_pdih(logger(pd[kk].ai),logger(pd[kk].aj),
			       logger(pd[kk].ak),logger(pd[kk].al),
			       pd[kk].type,ilist,iparams,TRUE,atoms,end_margin,
			       weight,d);
			       
	/* VIRTUAL DISTANCES */
	if (bVir) 
	  nVdist += set_virtual (cda,MAXLOGG-1,end_margin,d,natoms);
      }
      nlogg=0;
    }
  }
  pr_ndist(log,"ASN",n14dist,0,0,0,nVdist);

  sfree(cda);
}

/**********************************************************
 *
 *     G L U T A M I N E
 *
 **********************************************************/
void gln (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	  real end_margin,t_ilist ilist[],t_iparams iparams[],bool bVir)
{
  static char *sa[] = { "CG", "CD", "OE1", "NE2", "HE21", "HE22", "VF" };
  t_cdatom *cda;
  t_cdpdih pd[] = { 
    { "CG",  "CD", "NE2", "HE22", TRANS }, 
    { "OE1", "CD", "NE2", "HE21", TRANS },
    { "CG",  "CD", "NE2", "HE21", CIS }, 
    { "OE1", "CD", "NE2", "HE22", CIS }
  };
  int    natoms,i,j,kk,q,residnr,oldresidnr;
  int    n14dist,nVdist,atom_index,nlogg;
  real   blen,lb,ub,angle;

  cda        = init_cdatoms(MAXLOGG,sa);
  natoms     = atoms->nr;
  n14dist    = 0;
  nVdist     = 0;
  residnr    = -1;
  oldresidnr = -1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),"GLN") == 0) {
      oldresidnr = residnr;
      residnr    = atoms->atom[i].resnr;
      if (oldresidnr == residnr)
	continue;
      
      nlogg = set_cdatoms(atoms,i,residnr,MAXLOGG,cda);      
      
      if ( ((nlogg == MAXLOGG-1) && !bVir) || ((nlogg == MAXLOGG) && bVir) ) {

	pr_logg(debug,MAXLOGG-1,cda,bVir,"Gln");

	/* Proper dihedrals */
	for(kk=0; (kk<NPD); kk++)
	  n14dist += do_a_pdih(logger(pd[kk].ai),logger(pd[kk].aj),
			       logger(pd[kk].ak),logger(pd[kk].al),
			       pd[kk].type,ilist,iparams,TRUE,atoms,end_margin,
			       weight,d);

	/* VIRTUAL DISTANCES */
	if (bVir) 
	  nVdist += set_virtual (cda,MAXLOGG-1,end_margin,d,natoms);
      }
      nlogg=0;
    }
  }
  pr_ndist(log,"GLN",n14dist,0,0,0,nVdist);

  sfree(cda);
}

/**********************************************************
 *
 *     H I S T I D I N E
 *
 **********************************************************/
static void hisb_15_type2(int ai,int aj,int ak,int al,int am,
			  t_ilist ilist[],t_iparams iparams[],
			  real *lb,t_atoms *atoms)
{
  /* Returns the length corresponding to a 1-5 distance */
  real rij,rjk,rkl,rlm,rik,rkm,rim;
  real thijk,thjkl,thklm,thikj,thmkl;
  real real_pi = M_PI*2.0;
  real half_pi = M_PI*0.5;

  /*  fprintf(stderr,"Got past initialisations\n");*/
  rij   = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  rjk   = lookup_bondlength(aj,ak,ilist,iparams,TRUE,atoms);
  rkl   = lookup_bondlength(ak,al,ilist,iparams,TRUE,atoms);
  rlm   = lookup_bondlength(al,am,ilist,iparams,TRUE,atoms);
  /*fprintf(stderr,"Got past lookup_bondlength\n");*/
  thijk = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  thjkl = lookup_angle(aj,ak,al,ilist,iparams,atoms);
  thklm = lookup_angle(ak,al,am,ilist,iparams,atoms);
  /*fprintf(stderr,"Got past lookup_angle\n");*/
  /*fprintf(stderr,"%g %g %g %g %g %g %g\n",rij,rjk,rkl,rlm,RAD2DEG*thijk,
    RAD2DEG*thjkl,RAD2DEG*thklm);*/
  rik   =  sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(thijk));
  rkm   =  sqrt(rkl*rkl+rlm*rlm-2.0*rkl*rlm*cos(thklm));
  /*fprintf(stderr,"Got past angle_length\n");*/

  /* Compute angle thikj using law of sines */
  thikj = asin(rij*sin(thijk)/rik);
  
  /* Compute thmkl using law of sines */
  thmkl = asin(rlm*sin(thklm)/rkm);

  /* Compute cis length using law of cosines */
  *lb = sqrt(rik*rik+rkm*rkm-2.0*rik*rkm*cos(real_pi-thikj-thmkl-thjkl));
}

void hisb (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	   real ring_margin,t_ilist ilist[],t_iparams iparams[],bool bVir)
{
  static char *sa[] = { "CB", "CG", "ND1", "CE1", "NE2", "CD2", 
			"HE1", "HE2", "HD2", "VF" };
  t_cdatom *cda;
  t_cdpdih pd[] = {
    { "CB",  "CG",  "ND1", "CE1", TRANS }, 
    { "CB",  "CG",  "CD2", "NE2", TRANS }, 
    { "HE1", "CE1", "ND1", "CG",  TRANS }, 
    { "HE1", "CE1", "NE2", "CD2", TRANS },
    { "HE2", "NE2", "CE1", "ND1", TRANS }, 
    { "HE2", "NE2", "CD2", "CG",  TRANS },
    { "HE2", "NE2", "CE1", "HE1", CIS }, 
    { "HE2", "NE2", "CD2", "HD2", CIS },
    { "CB",  "CG",  "CD2", "HD2", CIS }, 
    { "ND1", "CG",  "CD2", "HD2", TRANS },
    { "CE1", "NE2", "CD2", "HD2", TRANS }
  };
  int    natoms,i,j,kk,q,residnr,oldresidnr;
  int    n14dist,n15dist,nVdist,atom_index,nlogg;
  real   blen,lb,ub,angle;

  cda     = init_cdatoms(MAXLOGG,sa);
  natoms  = atoms->nr;
  n14dist = 0;
  n15dist = 0;
  nVdist  = 0;
  residnr = -1;
  oldresidnr=-1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),"HISB") == 0) {
      oldresidnr=residnr;
      residnr=atoms->atom[i].resnr;
      if (oldresidnr == residnr) 
	continue;
      
      nlogg = set_cdatoms(atoms,i,residnr,MAXLOGG,cda);      
      
      if ( ((nlogg == MAXLOGG-1) && !bVir) || ((nlogg == MAXLOGG) && bVir) ) {

	pr_logg(debug,MAXLOGG-1,cda,bVir,"Hisb");
	
	/* Proper dihedrals */
	for(kk=0; (kk<NPD); kk++)
	  n14dist += do_a_pdih(logger(pd[kk].ai),logger(pd[kk].aj),
			       logger(pd[kk].ak),logger(pd[kk].al),
			       pd[kk].type,ilist,iparams,TRUE,atoms,ring_margin,
			       weight,d);
	
	/*SETDISTANCE for HE1 and HD2 (1-5) */
	hisb_15_type2(logger("HE1"),logger("CE1"),logger("NE2"),logger("CD2"),
		      logger("HD2"),ilist,iparams,&blen,atoms);
	n15dist += do_a_dist(logger("HE1"),logger("HD2"),natoms,ring_margin,
			     weight,d,blen);

	/*SETDISTANCE for CB and HE1 (1-5) */
	hisb_15_type2(logger("CB"),logger("CG"),logger("ND1"),logger("CE1"),
		      logger("HE1"),ilist,iparams,&blen,atoms);
	n15dist += do_a_dist(logger("CB"),logger("HE1"),natoms,ring_margin,
			     weight,d,blen);
		      
	/*SETDISTANCE for CB and HE2 (1-5) */
	hisb_15_type2(logger("CB"),logger("CG"),logger("CD2"),logger("NE2"),
		      logger("HE2"),ilist,iparams,&blen,atoms);
	n15dist += do_a_dist(logger("CB"),logger("HE2"),natoms,ring_margin,
			     weight,d,blen);
	
	/* VIRTUAL DISTANCES */
	if (bVir) 
	  nVdist = set_virtual (cda,MAXLOGG,ring_margin,d,natoms);
      }
    }
  }
  pr_ndist(log,"HIS",n14dist,n15dist,0,0,nVdist);

  sfree(cda);
}

/**********************************************************
 *
 *     I S O L E U C I N E
 *
 **********************************************************/
void ile (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	  real ile_margin,t_ilist ilist[],t_iparams iparams[])
{
  static char *sa[] = { "CA", "CB", "HB", "CD", "HD1", "HD2", "HD3",
			"CG1", "HG11", "HG12", 
			"CG2", "HG21", "HG22", "HG23" };
  t_cdatom *cda;
  t_cdpdih pd[] = {
    { "HD1",  "CD",  "CG1", "HG11", GAUCHE }, 
    { "HD1",  "CD",  "CG1", "CB",   TRANS }, 
    { "HD1",  "CD",  "CG1", "HG12", GAUCHE },
    { "HD2",  "CD",  "CG1", "HG11", TRANS },
    { "HD2",  "CD",  "CG1", "CB",   GAUCHE }, 
    { "HD2",  "CD",  "CG1", "HG12", GAUCHE },
    { "HD3",  "CD",  "CG1", "HG12", TRANS }, 
    { "HD3",  "CD",  "CG1", "HG11", GAUCHE }, 
    { "HD3",  "CD",  "CG1", "CB",   GAUCHE },
    { "HG21", "CG2", "CB",  "CA",   TRANS },
    { "HG21", "CG2", "CB",  "CG1",  GAUCHE }, 
    { "HG21", "CG2", "CB",  "HB",   GAUCHE },
    { "HG22", "CG2", "CB",  "HB",   TRANS }, 
    { "HG22", "CG2", "CB",  "CG1",  GAUCHE }, 
    { "HG22", "CG2", "CB",  "CA",   GAUCHE },
    { "HG23", "CG2", "CB",  "CG1",  TRANS },
    { "HG23", "CG2", "CB",  "CA",   GAUCHE }, 
    { "HG23", "CG2", "CB",  "HB",   GAUCHE }
  };
  
  int   natoms,i,j,kk,residnr,oldresidnr;
  int   n14dist,n15dist,atom_index,nlogg;
  real  blen,lb,ub,angle;
  real  pi = M_PI;

  cda        = init_cdatoms(MAXLOGG,sa);
  natoms     = atoms->nr;
  n14dist    = 0;
  n15dist    = 0;
  residnr    = -1;
  oldresidnr = -1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),"ILE") == 0) {
      oldresidnr=residnr;
      residnr=atoms->atom[i].resnr;
      if (oldresidnr == residnr)
	continue;

      nlogg = set_cdatoms(atoms,i,residnr,MAXLOGG,cda);
	
      if ( nlogg == MAXLOGG ) {
	pr_logg(debug,MAXLOGG,cda,FALSE,"Ile");

	/* Proper dihedrals */
	for(kk=0; (kk<NPD); kk++)
	  n14dist += do_a_pdih(logger(pd[kk].ai),logger(pd[kk].aj),
			       logger(pd[kk].ak),logger(pd[kk].al),
			       pd[kk].type,ilist,iparams,TRUE,atoms,
			       ile_margin,weight,d);
      }
    }
  }
  pr_ndist(log,"ILE",n14dist,0,0,0,0);

  sfree(cda);
}

/**********************************************************
 *
 *     L E U C I N E
 *
 **********************************************************/
void leu (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	  real leu_margin,t_ilist ilist[],t_iparams iparams[])
{
  static char *sa[] = { "CB", "CG", "HG", "CD1", "HD11", "HD12", "HD13",
			"CD2", "HD21", "HD22", "HD23" };
  t_cdatom *cda;
  t_cdpdih pd[] = {
    { "HD11", "CD1", "CG", "HG",  TRANS }, 
    { "HD12", "CD1", "CG", "CB",  TRANS },
    { "HD13", "CD1", "CG", "CD2", TRANS }, 
    { "HD21", "CD2", "CG", "CB",  TRANS },
    { "HD22", "CD2", "CG", "HG",  TRANS }, 
    { "HD23", "CD2", "CG", "CD1", TRANS },
    { "HD11", "CD1", "CG", "CB",  GAUCHE }, 
    { "HD11", "CD1", "CG", "CD2", GAUCHE },
    { "HD12", "CD1", "CG", "CD2", GAUCHE }, 
    { "HD12", "CD1", "CG", "HG",  GAUCHE },
    { "HD13", "CD1", "CG", "HG",  GAUCHE }, 
    { "HD13", "CD1", "CG", "CB",  GAUCHE },
    { "HD21", "CD2", "CG", "HG",  GAUCHE }, 
    { "HD21", "CD2", "CG", "CD1", GAUCHE },
    { "HD22", "CD2", "CG", "CD1", GAUCHE }, 
    { "HD22", "CD2", "CG", "CB",  GAUCHE },
    { "HD23", "CD2", "CG", "CB",  GAUCHE }, 
    { "HD23", "CD2", "CG", "HG",  GAUCHE }
  };
  t_cd15g cd15g[] = {
    { "HD11", "CD1", "CG", "CD2", "HD21", M_PI, M_PI43, M_PI43 },
    { "HD11", "CD1", "CG", "CD2", "HD22", M_PI, M_PI43, M_PI23 },
    { "HD11", "CD1", "CG", "CD2", "HD23", M_PI, M_PI43, 0 },
    { "HD12", "CD1", "CG", "CD2", "HD21", M_PI, M_PI23, M_PI43 },
    { "HD12", "CD1", "CG", "CD2", "HD22", M_PI, M_PI23, M_PI23 },
    { "HD12", "CD1", "CG", "CD2", "HD23", M_PI, M_PI23, 0 },
    { "HD13", "CD1", "CG", "CD2", "HD21", M_PI, 0,      M_PI43 },
    { "HD13", "CD1", "CG", "CD2", "HD22", M_PI, 0,      M_PI23 },
    { "HD13", "CD1", "CG", "CD2", "HD23", M_PI, 0,      0 }
  };
  int   natoms,i,j,kk,residnr,oldresidnr;
  int   n14dist,n15dist,atom_index,nlogg;
  real  blen,lb,ub,angle;
  real  pi = M_PI;

  cda        =  init_cdatoms(MAXLOGG,sa);  
  natoms     =  atoms->nr;
  n14dist    =  0;
  n15dist    =  0;
  residnr    = -1;
  oldresidnr = -1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),"LEU") == 0) {
      oldresidnr = residnr;
      residnr    = atoms->atom[i].resnr;
      if (oldresidnr == residnr)
	continue;
      
      nlogg = set_cdatoms(atoms,i,residnr,MAXLOGG,cda);

      if ( nlogg == MAXLOGG ) {
	pr_logg(debug,MAXLOGG,cda,FALSE,"Leu");

	/* Proper dihedrals */
	for(kk=0; (kk<NPD); kk++)
	  n14dist += do_a_pdih(logger(pd[kk].ai),logger(pd[kk].aj),
			       logger(pd[kk].ak),logger(pd[kk].al),
			       pd[kk].type,ilist,iparams,TRUE,atoms,
			       leu_margin,weight,d);
	
	/* 1-5 dihedrals */
	for(kk=0; (kk<asize(cd15g)); kk++) {
	  gauche15(logger(cd15g[kk].ai),logger(cd15g[kk].aj),
		   logger(cd15g[kk].ak),logger(cd15g[kk].al),
		   logger(cd15g[kk].am),
		   cd15g[kk].om1,cd15g[kk].om2,cd15g[kk].om3,
		   ilist,iparams,&blen,atoms);
	  n15dist += do_a_dist(logger(cd15g[kk].ai),logger(cd15g[kk].am),
			       natoms,leu_margin,weight,d,blen);
	}			
      }
    }
  }
  pr_ndist(log,"LEU",n14dist,n15dist,0,0,0);

  sfree(cda);
}

/**********************************************************
 *
 *     P H E N Y L A L A N I N E
 *
 **********************************************************/
static void phe_15_CBHE(int ai,int aj,int ak,int al,int am,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,t_atoms *atoms)
{
  /* Returns the length corresponding to a 1-5 distance */
  real rij,rjk,rkl,rlm,rik,rkm,rim;
  real thijk,thjkl,thklm,thikj,thmkl;
  real real_pi = M_PI*2.0;
  real half_pi = M_PI*0.5;

  rij   = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  rjk   = lookup_bondlength(aj,ak,ilist,iparams,TRUE,atoms);
  rkl   = lookup_bondlength(ak,al,ilist,iparams,TRUE,atoms);
  rlm   = lookup_bondlength(al,am,ilist,iparams,TRUE,atoms);

  thijk = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  thjkl = lookup_angle(aj,ak,al,ilist,iparams,atoms);
  thklm = lookup_angle(ak,al,am,ilist,iparams,atoms);

  rik   =  sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(thijk));
  rkm   =  sqrt(rkl*rkl+rlm*rlm-2.0*rkl*rlm*cos(thklm));

  /* Compute angle thikj using law of sines */
  thikj = asin(rij*sin(thijk)/rik);
  
  /* Compute thmkl using law of sines */
  thmkl = asin(rlm*sin(thklm)/rkm);

  /* Compute cis length using law of cosines */
  *lb = sqrt(rik*rik+rkm*rkm-2.0*rik*rkm*cos(real_pi-thikj-thmkl-thjkl));
}

static void phe_15_CBCZ(int ai,int aj,int ak,int al,int am,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,t_atoms *atoms)
{
  /* Returns the length corresponding to a 1-5 distance */
  real rij,rjk,rkl,rlm,rik,rkm;
  real thijk,thjkl,thklm,thikj,thmkl;
  real real_pi = M_PI*2.0;
  real half_pi = M_PI*0.5;
  
  rij   = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  rjk   = lookup_bondlength(aj,ak,ilist,iparams,TRUE,atoms);
  rkl   = lookup_bondlength(ak,al,ilist,iparams,TRUE,atoms);
  rlm   = lookup_bondlength(al,am,ilist,iparams,TRUE,atoms);
  thijk = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  thjkl = lookup_angle(aj,ak,al,ilist,iparams,atoms);
  thklm = lookup_angle(ak,al,am,ilist,iparams,atoms);
  rik   = sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(thijk));
  rkm   = sqrt(rkl*rkl+rlm*rlm-2.0*rkl*rlm*cos(thklm));
  
  /* Compute angle thikj using law of sines */
  thikj = asin(rij*sin(thijk)/rik);
  
  /* Compute thmkl using law of sines */
  thmkl = asin(rlm*sin(thklm)/rkm);

  /* Compute cis length using law of cosines */
  *lb = sqrt(rik*rik+rkm*rkm-2.0*rik*rkm*cos(thjkl-thmkl+thikj));
}

static void phe_15_CGHZ(int ai,int aj,int ak,int al,int am,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,t_atoms *atoms)
{
  /* Returns the length corresponding to a 1-5 distance */
  real rij,rjk,rkl,rlm,rik,rkm;
  real thijk,thjkl,thklm,thikj,thmkl;
  real real_pi = M_PI*2.0;

  rij   = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  rjk   = lookup_bondlength(aj,ak,ilist,iparams,TRUE,atoms);
  rkl   = lookup_bondlength(ak,al,ilist,iparams,TRUE,atoms);
  rlm   = lookup_bondlength(al,am,ilist,iparams,TRUE,atoms);
  thijk = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  thjkl = lookup_angle(aj,ak,al,ilist,iparams,atoms);
  thklm = lookup_angle(ak,al,am,ilist,iparams,atoms);
  rik   = sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(thijk));
  rkm   = sqrt(rkl*rkl+rlm*rlm-2.0*rkl*rlm*cos(thklm));

  /* Compute angle thikj using law of sines */
  thikj = asin(rij*sin(thijk)/rik);
  
  /* Compute thmkl using law of sines */
  thmkl = asin(rlm*sin(thklm)/rkm);

  /* Compute cis length using law of cosines */
  *lb = sqrt(rik*rik+rkm*rkm-2.0*rik*rkm*cos(thjkl-thikj+thmkl));
}

static void phe_16_type2(int ai,int aj,int ak,int al,int am,int an,
		  t_ilist ilist[],t_iparams iparams[],
		  real *lb,t_atoms *atoms)
{
  /* Returns the length corresponding to a 1-6 distance */
  real rij,rjk,rkl,rlm,rmn,rik,rln,ril;
  real thijk,thjkl,thklm,thlmn,thikj,thikl,thilk,thnlm,thiln;
  real real_pi = M_PI*2.0;    
  real half_pi = M_PI*0.5;
  
  rij   = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  rjk   = lookup_bondlength(aj,ak,ilist,iparams,TRUE,atoms);
  rkl   = lookup_bondlength(ak,al,ilist,iparams,TRUE,atoms);
  rlm   = lookup_bondlength(al,am,ilist,iparams,TRUE,atoms);
  rmn   = lookup_bondlength(am,an,ilist,iparams,TRUE,atoms);
  thijk = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  thjkl = lookup_angle(aj,ak,al,ilist,iparams,atoms);
  thklm = lookup_angle(ak,al,am,ilist,iparams,atoms);
  thlmn = lookup_angle(al,am,an,ilist,iparams,atoms);

  /* Compute rik and rlm */
  rik   = sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(thijk));
  rln   = sqrt(rlm*rlm+rmn*rmn-2.0*rlm*rmn*cos(thlmn));

  /* Compute thikj */
  thikj = asin(rij*sin(thijk)/rik);
  thikl = thikj+thjkl;
  
  /* Compute ril */
  ril   = sqrt(rik*rik+rkl*rkl-2.0*rik*rkl*cos(thikl));

  /* Compute thilk, thlnm and thiln */
  thilk = asin(rik*sin(thikl)/ril);
  thnlm = asin(rmn*sin(thlmn)/rln);
  thiln = thklm-thilk+thnlm;

  *lb = sqrt(ril*ril+rln*rln-2.0*ril*rln*cos(thiln));
}

typedef void phe_15_func(int ai,int aj,int ak,int al,int am,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,t_atoms *atoms);
			 
void phe_tyr(FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	     real ring_margin,t_ilist ilist[],t_iparams iparams[],bool bVir,
	     bool bPhe)
{
  static char *sa[] = { "CB", "CG", "CD2", "HD2", "CE2", "HE2", "CZ",
			"HZ", "CE1", "HE1", "CD1", "HD1", "VF" };
  t_cdatom *cda;
  t_cd15 cd15[] = {
    { "HD2", "CD2", "CG", "CD1", "HD1", 1},
    { "HE2", "CE2", "CZ", "CE1", "HE1", 1},
    { "HE1", "CE1", "CZ", "CE2", "CD2", 2},
    { "HD2", "CD2", "CE2", "CZ", "CE1", 2},
    { "HD1", "CD1", "CG", "CD2", "CE2", 2},
    { "HE2", "CE2", "CZ", "CE1", "CD1", 2},
    { "CB",  "CG", "CD2", "CE2", "HE2", 1},
    { "CB",  "CG", "CD1", "CE1", "HE1", 1},
    { "CB",  "CG", "CD2", "CE2", "CZ",  2}
  };
  phe_15_func *phe15[3] = {
    phe_15_CBHE, phe_15_CBCZ, phe_15_CGHZ
  };
  t_cdpdih pd[] = {
    { "CD2", "CG",  "CD1", "HD1", TRANS }, 
    { "CZ",  "CE1", "CD1", "HD1", TRANS },
    { "CD2", "CG",  "CD1", "CE1", CIS }, 
    { "CE2", "CZ",  "CE1", "CD1", CIS },
    { "CB",  "CG",  "CD1", "CE1", TRANS }, 
    { "CB",  "CG",  "CD2", "CE2", TRANS },
    { "CD1", "CG",  "CD2", "HD2", TRANS }, 
    { "CZ",  "CE2", "CD2", "HD2", TRANS },
    { "CG",  "CD2", "CE2", "HE2", TRANS }, 
    { "CE1", "CZ",  "CE2", "HE2", TRANS },
    { "CD1", "CE1", "CZ",  "HZ",  TRANS }, 
    { "CD2", "CE2", "CZ",  "HZ",  TRANS },
    { "CE2", "CZ",  "CE1", "HE1", TRANS }, 
    { "CG",  "CD1", "CE1", "HE1", TRANS },
    { "CB",  "CG",  "CD1", "HD1", CIS }, 
    { "CB",  "CG",  "CD2", "HD2", CIS },
    { "CG",  "CD2", "CE2", "CZ",  CIS }, 
    { "HD2", "CD2", "CE2", "HE2", CIS },
    { "HE2", "CE2", "CZ",  "HZ",  CIS }, 
    { "HZ",  "CZ",  "CE1", "HE1", CIS },
    { "HE1", "CE1", "CD1", "HD1", CIS }
  };
  int    natoms,i,j,k,kk,q,residnr,oldresidnr,nhz;
  int    n14dist,n15dist,n16dist,nVdist,atom_index,nlogg;
  real   blen,lb,ub,angle;
  char   NHZ[8],RES[8];

  cda        = init_cdatoms(MAXLOGG,sa);
  nhz        = 7;
  /* Consistency check */
  if (strcmp(cda[nhz].anm,"HZ") != 0)
    gmx_fatal(FARGS,"wrong atom name %s, expeced HZ",cda[nhz].anm);
  if (!bPhe) {
    strcpy(NHZ,"OH");
    strcpy(cda[nhz].anm,"OH");
    strcpy(RES,"TYR");
  }
  else {
    strcpy(NHZ,"HZ");
    strcpy(RES,"PHE");
  }
  natoms     = atoms->nr;
  n14dist    = 0;
  n15dist    = 0;
  n16dist    = 0;
  nVdist     = 0;
  residnr    = -1;
  oldresidnr =-1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),RES) == 0) {
      oldresidnr = residnr;
      residnr    = atoms->atom[i].resnr;
      if (oldresidnr == residnr)
	continue;

      nlogg = set_cdatoms(atoms,i,residnr,MAXLOGG,cda);
      
      if ( ((nlogg == MAXLOGG-1) && !bVir) || ((nlogg == MAXLOGG) && bVir) ) {

	pr_logg(debug,MAXLOGG-1,cda,bVir,RES);

#define tp(s) (strcmp(s,"HZ") == 0) ? NHZ : s
	/* Proper dihedrals */
	for(kk=0; (kk<NPD); kk++)
	  n14dist += do_a_pdih(logger(tp(pd[kk].ai)),logger(tp(pd[kk].aj)),
			       logger(tp(pd[kk].ak)),logger(tp(pd[kk].al)),
			       pd[kk].type,ilist,iparams,TRUE,atoms,
			       ring_margin,weight,d);
	
	for(k=0; (k<asize(cd15)); k++) {
	  phe15[cd15[k].type-1](logger(cd15[k].ai),logger(cd15[k].aj),
				logger(cd15[k].ak),logger(cd15[k].al),
				logger(cd15[k].am),
				ilist,iparams,&blen,atoms);
	  n15dist += do_a_dist(logger(cd15[k].ai),
			       logger(cd15[k].am),natoms,ring_margin,
			       weight,d,blen);
	}

	/*SETDISTANCE for HD2 and HZ (1-5) */
	phe_15_CBHE(logger("HD2"),logger("CD2"),logger("CE2"),logger("CZ"),
		    logger(NHZ),ilist,iparams,&blen,atoms);
	n15dist += do_a_dist(logger("HD2"),logger(NHZ),natoms,ring_margin,
			     weight,d,blen);
			     
	/*SETDISTANCE for HZ and HD1 (1-5) */
	phe_15_CBHE(logger(NHZ),logger("CZ"),logger("CE1"),logger("CD1"),
		    logger("HD1"),ilist,iparams,&blen,atoms);
	n15dist += do_a_dist(logger(NHZ),logger("HD1"),natoms,ring_margin,
			     weight,d,blen);

	/*SETDISTANCE for CG and HZ (1-5) */
	phe_15_CGHZ(logger("CG"),logger("CD2"),logger("CE2"),logger("CZ"),
		    logger(NHZ),ilist,iparams,&blen,atoms);
	n15dist += do_a_dist(logger("CG"),logger(NHZ),natoms,ring_margin,
			     weight,d,blen);
	
	/* 1111111111-------------666666666666 */
	/*SETDISTANCE for CB and HZ (1-6) */
	phe_16_type2(logger("CB"),logger("CG"),logger("CD2"),logger("CE2"),
		     logger("CZ"),logger(NHZ),ilist,iparams,&blen,atoms);
	n16dist += do_a_dist(logger("CB"),logger(NHZ),natoms,ring_margin,
			     weight,d,blen);

	/*SETDISTANCE for HD1 and HE2 (1-6) */
	phe_16_type2(logger("HD1"),logger("CD1"),logger("CG"),logger("CD2"),
		     logger("CE2"),logger("HE2"),ilist,iparams,&blen,atoms);
	n16dist += do_a_dist(logger("HD1"),logger("HE2"),natoms,ring_margin,
			     weight,d,blen);

	/*SETDISTANCE for HD2 and HE1 (1-6) */
	phe_16_type2(logger("HD2"),logger("CD2"),logger("CE2"),logger("CZ"),
		     logger("CE1"),logger("HE1"),ilist,iparams,&blen,atoms);
	n16dist += do_a_dist(logger("HD2"),logger("HE1"),natoms,ring_margin,
			     weight,d,blen);
	
	/* VIRTUAL DISTANCES */
	if (bVir) 
	  nVdist += set_virtual (cda,MAXLOGG-1,ring_margin,d,natoms);
      }
    }
  }
  pr_ndist(log,RES,n14dist,n15dist,n16dist,0,nVdist);

  sfree(cda);
}

void phe(FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	 real ring_margin,t_ilist ilist[],t_iparams iparams[],bool bVir)
{
  phe_tyr(log,d,idef,atoms,weight,ring_margin,ilist,iparams,bVir,TRUE);
}

/**********************************************************
 *
 *      T Y R O S I N E
 *
 **********************************************************/
void tyr(FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	 real ring_margin,t_ilist ilist[],t_iparams iparams[],bool bVir)
{
  phe_tyr(log,d,idef,atoms,weight,ring_margin,ilist,iparams,bVir,FALSE);
}

/**********************************************************
 *
 *     T R Y P T O P H A N
 *
 **********************************************************/
/* Remember that the order of atoms sent in to the 1-5-routines */
/* is important */

static void trp_15_type1(int ai,int aj,int ak,int al,int am,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,t_atoms *atoms)
{
  /* Returns the length corresponding to a 1-5 distance */
  real rij,rjk,rkl,rlm,rik,rkm;
  real thijk,thjkl,thklm,thikj,thmkl;
  real real_pi = M_PI*2.0;
  real half_pi = M_PI*0.5;
  
  rij   = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  rjk   = lookup_bondlength(aj,ak,ilist,iparams,TRUE,atoms);
  rkl   = lookup_bondlength(ak,al,ilist,iparams,TRUE,atoms);
  rlm   = lookup_bondlength(al,am,ilist,iparams,TRUE,atoms);
  thijk = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  thjkl = lookup_angle(aj,ak,al,ilist,iparams,atoms);
  thklm = lookup_angle(ak,al,am,ilist,iparams,atoms);
  rik   = sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(thijk));
  rkm   = sqrt(rkl*rkl+rlm*rlm-2.0*rkl*rlm*cos(thklm));
  
  /* Compute angle thikj using law of sines */
  thikj = asin(rij*sin(thijk)/rik);
  
  /* Compute thmkl using law of sines */
  thmkl = asin(rlm*sin(thklm)/rkm);

  /* Compute cis length using law of cosines */
  *lb = sqrt(rik*rik+rkm*rkm-2.0*rik*rkm*cos(thjkl-thmkl+thikj));
}

static void trp_15_type2(int ai,int aj,int ak,int al,int am,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,t_atoms *atoms)
{
  /* Returns the length corresponding to a 1-5 distance */
  real rij,rjk,rkl,rlm,rik,rkm,rim;
  real thijk,thjkl,thklm,thikj,thmkl;
  real real_pi = M_PI*2.0;
  real half_pi = M_PI*0.5;

  rij   = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  rjk   = lookup_bondlength(aj,ak,ilist,iparams,TRUE,atoms);
  rkl   = lookup_bondlength(ak,al,ilist,iparams,TRUE,atoms);
  rlm   = lookup_bondlength(al,am,ilist,iparams,TRUE,atoms);
  thijk = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  thjkl = lookup_angle(aj,ak,al,ilist,iparams,atoms);
  thklm = lookup_angle(ak,al,am,ilist,iparams,atoms);

  rik   =  sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(thijk));
  rkm   =  sqrt(rkl*rkl+rlm*rlm-2.0*rkl*rlm*cos(thklm));

  /* Compute angle thikj using law of sines */
  thikj = asin(rij*sin(thijk)/rik);
  
  /* Compute thmkl using law of sines */
  thmkl = asin(rlm*sin(thklm)/rkm);

  /* Compute cis length using law of cosines */
  *lb = sqrt(rik*rik+rkm*rkm-2.0*rik*rkm*cos(real_pi-thikj-thmkl-thjkl));
}

static void trp_15_type3(int ai,int aj,int ak,int al,int am,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,t_atoms *atoms)
{
  /* Returns the length corresponding to a 1-5 distance */
  real rij,rjk,rkl,rlm,rik,rkm;
  real thijk,thjkl,thklm,thikj,thmkl;
  real real_pi = M_PI*2.0;
  real half_pi = M_PI*0.5;
  
  rij   = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  rjk   = lookup_bondlength(aj,ak,ilist,iparams,TRUE,atoms);
  rkl   = lookup_bondlength(ak,al,ilist,iparams,TRUE,atoms);
  rlm   = lookup_bondlength(al,am,ilist,iparams,TRUE,atoms);
  thijk = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  thjkl = lookup_angle(aj,ak,al,ilist,iparams,atoms);
  thklm = lookup_angle(ak,al,am,ilist,iparams,atoms);
  rik   = sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(thijk));
  rkm   = sqrt(rkl*rkl+rlm*rlm-2.0*rkl*rlm*cos(thklm));
  
  /* Compute angle thikj using law of sines */
  thikj = asin(rij*sin(thijk)/rik);
  
  /* Compute thmkl using law of sines */
  thmkl = asin(rlm*sin(thklm)/rkm);

  /* Compute cis length using law of cosines */
  *lb = sqrt(rik*rik+rkm*rkm-2.0*rik*rkm*cos(thjkl-thikj-thmkl));
}

static void trp_16_type1(int ai,int aj,int ak,int al,int am,int an,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,t_atoms *atoms)
{

  real rij,rjk,rkl,rlm,rmn,rik,ril,rim;
  real thijk,thjkl,thklm,thlmn,thikj,thilk,thiml,thimn,thikl,thilm;
  real real_pi = M_PI*2.0;
  real half_pi = M_PI*0.5;
  
  rij   = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  rjk   = lookup_bondlength(aj,ak,ilist,iparams,TRUE,atoms);
  rkl   = lookup_bondlength(ak,al,ilist,iparams,TRUE,atoms);
  rlm   = lookup_bondlength(al,am,ilist,iparams,TRUE,atoms);
  rmn   = lookup_bondlength(am,an,ilist,iparams,TRUE,atoms);
  thijk = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  thjkl = lookup_angle(aj,ak,al,ilist,iparams,atoms);
  thklm = lookup_angle(ak,al,am,ilist,iparams,atoms);
  thlmn = lookup_angle(al,am,an,ilist,iparams,atoms);

  rik   = sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(thijk));
  
  /* Compute angle thikl */
  thikj = asin(rij*sin(thijk)/rik);
  thikl = thjkl-thikj;
  
  /* Compute ril */
  ril   = sqrt(rik*rik+rkl*rkl-2.0*rik*rkl*cos(thikl));

  /* Compute thilk */
  thilk = asin(rik*sin(thikl)/ril);

  /* Compute rim */
  thilm = thilk+thklm;
  rim   = sqrt(ril*ril+rlm*rlm-2.0*ril*rlm*cos(thilm));

  /* Compute rin */
  thiml = asin(ril*sin(thilm)/rim);
  thimn = thlmn-thiml;
  *lb = sqrt(rim*rim+rmn*rmn-2.0*rim*rmn*cos(thimn));
}

static void trp_16_type2(int ai,int aj,int ak,int al,int am,int an,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,t_atoms *atoms)
{
  /* Returns the length corresponding to a 1-6 distance */

  real rij,rjk,rkl,rlm,rmn,rik,rln,ril;
  real thijk,thjkl,thklm,thlmn,thikj,thikl,thilk,thnlm,thiln;
  real real_pi = M_PI*2.0;    
  real half_pi = M_PI*0.5;
  
  rij   = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  rjk   = lookup_bondlength(aj,ak,ilist,iparams,TRUE,atoms);
  rkl   = lookup_bondlength(ak,al,ilist,iparams,TRUE,atoms);
  rlm   = lookup_bondlength(al,am,ilist,iparams,TRUE,atoms);
  rmn   = lookup_bondlength(am,an,ilist,iparams,TRUE,atoms);
  thijk = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  thjkl = lookup_angle(aj,ak,al,ilist,iparams,atoms);
  thklm = lookup_angle(ak,al,am,ilist,iparams,atoms);
  thlmn = lookup_angle(al,am,an,ilist,iparams,atoms);

  /* Compute rik and rlm */
  rik   = sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(thijk));
  rln   = sqrt(rlm*rlm+rmn*rmn-2.0*rlm*rmn*cos(thlmn));

  /* Compute thikj */
  thikj = asin(rij*sin(thijk)/rik);
  thikl = thikj+thjkl;
  
  /* Compute ril */
  ril   = sqrt(rik*rik+rkl*rkl-2.0*rik*rkl*cos(thikl));

  /* Compute thilk, thlnm and thiln */
  thilk = asin(rik*sin(thikl)/ril);
  thnlm = asin(rmn*sin(thlmn)/rln);
  thiln = thklm-thilk+thnlm;

  *lb = sqrt(ril*ril+rln*rln-2.0*ril*rln*cos(thiln));
}

static void trp_16_type3(int ai,int aj,int ak,int al,int am,int an,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,t_atoms *atoms)
{
  real rij,rjk,rkl,rlm,rmn,rik,ril,rim;
  real thijk,thjkl,thklm,thlmn,thikj,thilk,thiml,thimn,thikl,thilm;
  real real_pi = M_PI*2.0;
  real half_pi = M_PI*0.5;
  
  rij   = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  rjk   = lookup_bondlength(aj,ak,ilist,iparams,TRUE,atoms);
  rkl   = lookup_bondlength(ak,al,ilist,iparams,TRUE,atoms);
  rlm   = lookup_bondlength(al,am,ilist,iparams,TRUE,atoms);
  rmn   = lookup_bondlength(am,an,ilist,iparams,TRUE,atoms);
  thijk = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  thjkl = lookup_angle(aj,ak,al,ilist,iparams,atoms);
  thklm = lookup_angle(ak,al,am,ilist,iparams,atoms);
  thlmn = lookup_angle(al,am,an,ilist,iparams,atoms);

  rik   = sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(thijk));
  
  /* Compute angle thikl */
  thikj = asin(rij*sin(thijk)/rik);
  thikl = thjkl-thikj;
  
  /* Compute ril */
  ril   = sqrt(rik*rik+rkl*rkl-2.0*rik*rkl*cos(thikl));

  /* Compute thilk */
  thilk = asin(rik*sin(thikl)/ril);

  /* Compute rim */
  thilm = thilk+thklm;
  rim   = sqrt(ril*ril+rlm*rlm-2.0*ril*rlm*cos(thilm));

  /* Compute rin */
  thiml = asin(ril*sin(thilm)/rim);
  thimn = thlmn+thiml;
  *lb   = sqrt(rim*rim+rmn*rmn-2.0*rim*rmn*cos(thimn));
}


static void trp_16_type4(int ai,int aj,int ak,int al,int am,int an,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,t_atoms *atoms)
{
  real rij,rjk,rkl,rlm,rmn,rik,rln,ril;
  real thijk,thjkl,thklm,thlmn,thikj,thikl,thilk,thnlm,thiln;
  real real_pi = M_PI*2.0;
  real half_pi = M_PI*0.5;

  rij   = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  rjk   = lookup_bondlength(aj,ak,ilist,iparams,TRUE,atoms);
  rkl   = lookup_bondlength(ak,al,ilist,iparams,TRUE,atoms);
  rlm   = lookup_bondlength(al,am,ilist,iparams,TRUE,atoms);
  rmn   = lookup_bondlength(am,an,ilist,iparams,TRUE,atoms);
  thijk = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  thjkl = lookup_angle(aj,ak,al,ilist,iparams,atoms);
  thklm = lookup_angle(ak,al,am,ilist,iparams,atoms);
  thlmn = lookup_angle(al,am,an,ilist,iparams,atoms);
  
  /* Compute rik and rln */
  rik   = sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(thijk));
  rln   = sqrt(rlm*rlm+rmn*rmn-2.0*rlm*rmn*cos(thlmn));
  
  /* Compute thikj and thikl */
  thikj = asin(rij*sin(thijk)/rik);
  thikl = thikj+thjkl;

  /* Compute ril */
  ril = sqrt(rik*rik+rkl*rkl-2.0*rik*rkl*cos(thikl));

  /* Compute thilk, thnlm and thiln */
  thilk = asin(rik*sin(thikl)/ril);
  thnlm = asin(rmn*sin(thlmn)/rln);
  thiln = thilk+thklm+thnlm;

  /* Compute rin */
  *lb = sqrt(ril*ril+rln*rln-2.0*ril*rln*cos(thiln));
}

static void trp_17_type1(int ai,int aj,int ak,int al,int am,int an,int ao,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,t_atoms *atoms)
{
  real rij,rjk,rkl,rlm,rmn,rno,rik,rkm,rmo,rim;
  real thijk,thjkl,thklm,thlmn,thmno,thikj,thmkl,thikm,thomn,
    thkml,thimk,thimo;
  real real_pi = M_PI*2.0;
  real half_pi = M_PI*0.5;

  rij   = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  rjk   = lookup_bondlength(aj,ak,ilist,iparams,TRUE,atoms);
  rkl   = lookup_bondlength(ak,al,ilist,iparams,TRUE,atoms);
  rlm   = lookup_bondlength(al,am,ilist,iparams,TRUE,atoms);
  rmn   = lookup_bondlength(am,an,ilist,iparams,TRUE,atoms);
  rno   = lookup_bondlength(an,ao,ilist,iparams,TRUE,atoms);
  thijk = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  thjkl = lookup_angle(aj,ak,al,ilist,iparams,atoms);
  thklm = lookup_angle(ak,al,am,ilist,iparams,atoms);
  thlmn = lookup_angle(al,am,an,ilist,iparams,atoms);
  thmno = lookup_angle(am,an,ao,ilist,iparams,atoms);

  /* Compute rik, rkm, rmo */
  rik   = sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(thijk));
  rkm   = sqrt(rkl*rkl+rlm*rlm-2.0*rkl*rlm*cos(thklm));
  rmo   = sqrt(rmn*rmn+rno*rno-2.0*rmn*rno*cos(thmno));

  /* Compute thikj,thkml,thomn,thikm */
  thikj = asin(rij*sin(thijk)/rik);
  thmkl = asin(rlm*sin(thklm)/rkm);
  thikm = thikj+thjkl-thmkl;

  /* Compute rim */
  rim   = sqrt(rik*rik+rkm*rkm-2.0*rik*rkm*cos(thikm));

  /* Compute thimk,thmkl,thimo */
  thomn = asin(rno*sin(thmno)/rmo);
  thkml = asin(rkl*sin(thklm)/rkm);
  thimk = asin(rik*sin(thikm)/rim);
  thimo = thimk+thkml+thlmn+thomn;

  /* Compute rio */
  *lb = sqrt(rim*rim+rmo*rmo-2.0*rim*rmo*cos(thimo));
}

static void trp_17_type2(int ai,int aj,int ak,int al,int am,int an,int ao,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,t_atoms *atoms)
{
  real rij,rjk,rkl,rlm,rmn,rno,rik,rmo,ril,rlo;
  real thijk,thjkl,thklm,thlmn,thmno,thikj,thomn,thikl,thlmo,thkli,
    tholm,thilo;
  real real_pi = M_PI*2.0;
  real half_pi = M_PI*0.5;


  rij   = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  rjk   = lookup_bondlength(aj,ak,ilist,iparams,TRUE,atoms);
  rkl   = lookup_bondlength(ak,al,ilist,iparams,TRUE,atoms);
  rlm   = lookup_bondlength(al,am,ilist,iparams,TRUE,atoms);
  rmn   = lookup_bondlength(am,an,ilist,iparams,TRUE,atoms);
  rno   = lookup_bondlength(an,ao,ilist,iparams,TRUE,atoms);
  thijk = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  thjkl = lookup_angle(aj,ak,al,ilist,iparams,atoms);
  thklm = lookup_angle(ak,al,am,ilist,iparams,atoms);
  thlmn = lookup_angle(al,am,an,ilist,iparams,atoms);
  thmno = lookup_angle(am,an,ao,ilist,iparams,atoms);

  /* Compute rik, rmo */
  rik   = sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(thijk));
  rmo   = sqrt(rmn*rmn+rno*rno-2.0*rmn*rno*cos(thmno));
  
  /* Compute thikj,thomn, thikl and thlmo */
  thikj = asin(rij*sin(thijk)/rik);
  thomn = asin(rno*sin(thmno)/rmo);
  thikl = thikj+thjkl;
  thlmo = thlmn+thomn;
 
  /* Compute ril and rlo */
  ril  = sqrt(rik*rik+rkl*rkl-2.0*rik*rkl*cos(thikl));
  rlo  = sqrt(rlm*rlm+rmo*rmo-2.0*rlm*rmo*cos(thlmo));

  /* Compute thkli, tholm and thilo*/
  thkli = asin(rik*sin(thikl)/ril);
  tholm = asin(rmo*sin(thlmo)/rlo);
  thilo = thkli+tholm+thklm;

  /* Compute rio */
  *lb = sqrt(ril*ril+rlo*rlo-2.0*ril*rlo*cos(thilo));
}

/* FUNCTION TYPES */
typedef void trp_15_func(int ai,int aj,int ak,int al,int am,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,t_atoms *atoms);
typedef void trp_16_func(int ai,int aj,int ak,int al,int am,int an,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,t_atoms *atoms);
			 
void trp(FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	 real ring_margin,t_ilist ilist[],t_iparams iparams[],bool bVir)
{
  static char *sa[] = { "CB",  "CG",  "CD2", "CE3", "CZ3", "CH2", 
			"CZ2", "CE2", "NE1", "CD1", "HE3", "HZ3", 
			"HH2", "HZ2", "HE1", "HD1", "VF" };
  t_cdatom *cda;
  t_cd16 cd16[] = {
    { "CB",  "CG",  "CD2", "CE3", "CZ3", "CH2", 1 },
    { "CZ3", "CH2", "CZ2", "CE2", "NE1", "HE1", 1 },
    { "HH2", "CH2", "CZ2", "CE2", "CD2", "CG",  2 },
    { "HZ3", "CZ3", "CE3", "CD2", "CE2", "NE1", 2 },
    { "HZ2", "CZ2", "CH2", "CZ3", "CE3", "HE3", 2 },
    { "CB",  "CG",  "CD2", "CE3", "CZ3", "HZ3", 3 },
    { "HE1", "NE1", "CE2", "CZ2", "CH2", "HH2", 3 },
    { "HZ2", "CZ2", "CE2", "NE1", "CD1", "HD1", 3 },
    { "HE3", "CE3", "CD2", "CG",  "CD1", "HD1", 3 },
    { "CB",  "CG",  "CD2", "CE2", "CZ2", "HZ2", 4 },
    { "CZ3", "CE3", "CD2", "CG",  "CD1", "HD1", 4 },
    { "CH2", "CZ2", "CE2", "NE1", "CD1", "HD1", 4 },
    { "HZ3", "CZ3", "CE3", "CD2", "CG",  "CD1", 4 },
    { "HH2", "CH2", "CZ2", "CE2", "NE1", "CD1", 4 },
    { "HE3", "CE3", "CD2", "CE2", "NE1", "HE1", 4 }
  };
  trp_16_func *trp16[4] = {
    trp_16_type1, trp_16_type2, trp_16_type3, trp_16_type4
  };
  t_cd15 cd15[] = {
    { "HE1", "NE1", "CD1", "CG",  "CB",  2 },
    { "CB",  "CG",  "CD2", "CE2", "CZ2", 2 },
    { "CD1", "CG",  "CD2", "CE3", "CZ3", 2 },
    { "CB",  "CG",  "CD2", "CE3", "HE3", 3 },
    { "HD1", "CD1", "CG",  "CD2", "CE3", 2 },
    { "HD1", "CD1", "NE1", "CE2", "CZ2", 2 },
    { "HE1", "NE1", "CE2", "CZ2", "HZ2", 3 },
    { "CH2", "CZ2", "CE2", "NE1", "HE1", 1 },
    { "HE1", "NE1", "CE2", "CD2", "CE3", 2 },
    { "HZ2", "CZ2", "CE2", "CD2", "CG",  2 },
    { "HH2", "CH2", "CZ2", "CE2", "CD2", 1 },
    { "HZ3", "CZ3", "CE3", "CD2", "CE2", 1 },
    { "HE3", "CE3", "CD2", "CE2", "CZ2", 1 },
    { "HZ2", "CZ2", "CE2", "CD2", "CE3", 1 }, 
    { "CZ3", "CE3", "CD2", "CG",  "CB",  1 },
    { "CG",  "CD2", "CE2", "CZ2", "CH2", 1 },
    { "NE1", "CE2", "CD2", "CE3", "CZ3", 1 },
    { "CD1", "CG",  "CD2", "CE3", "HE3", 1 },
    { "CD1", "NE1", "CE2", "CZ2", "HZ2", 1 },
    { "NE1", "CE2", "CZ2", "CH2", "HH2", 2 },
    { "HE3", "CE3", "CZ3", "CH2", "HH2", 2 },
    { "HZ3", "CZ3", "CH2", "CZ2", "HZ2", 2 },
    { "CG",  "CD2", "CE3", "CZ3", "HZ3", 2 },
    { "CH2", "CZ2", "CE2", "NE1", "CD1", 2 },
    { "NE1", "CE2", "CD2", "CE3", "HE3", 2 }
  };
  trp_15_func *trp15[3] = {
    trp_15_type1, trp_15_type2, trp_15_type3
  };
  t_cdpdih pd[] = {
    { "CG",  "CD2", "CE3", "CZ3", TRANS }, 
    { "CD1", "NE1", "CE2", "CZ2", TRANS },
    { "NE1", "CE2", "CZ2", "CH2", TRANS }, 
    { "NE1", "CE2", "CD2", "CE3", TRANS },
    { "CD2", "CE3", "CZ3", "HZ3", TRANS }, 
    { "CE3", "CZ3", "CH2", "HH2", TRANS },
    { "CZ3", "CH2", "CZ2", "HZ2", TRANS }, 
    { "CG",  "CD2", "CE2", "CZ2", TRANS },
    { "CE2", "NE1", "CD1", "HD1", TRANS }, 
    { "CD2", "CG",  "CD1", "HD1", TRANS },
    { "CB",  "CG",  "CD2", "CE3", CIS }, 
    { "CB",  "CG",  "CD1", "HD1", CIS },
    { "CG",  "CD2", "CE3", "HE3", CIS }, 
    { "HE3", "CE3", "CZ3", "HZ3", CIS },
    { "HZ3", "CZ3", "CH2", "HH2", CIS }, 
    { "HH2", "CH2", "CZ2", "HZ2", CIS },
    { "HZ2", "CZ2", "CE2", "NE1", CIS }, 
    { "CZ2", "CE2", "NE1", "HE1", CIS },
    { "HE1", "NE1", "CD1", "HD1", CIS }, 
    { "CD2", "CE3", "CZ3", "CH2", CIS },
    { "CE3", "CZ3", "CH2", "CZ2", CIS }, 
    { "CZ3", "CH2", "CZ2", "CE2", CIS },
    { "CB",  "CG",  "CD2", "CE2", TRANS }, 
    { "CB",  "CG",  "CD1", "NE1", TRANS },
    { "CG",  "CD1", "NE1", "HE1", TRANS }, 
    { "CD2", "CE2", "CZ2", "HZ2", TRANS },
    { "CD2", "CE2", "NE1", "HE1", TRANS }, 
    { "CE3", "CD2", "CG",  "CD1", TRANS },
    { "CH2", "CZ3", "CE3", "HE3", TRANS }, 
    { "CZ2", "CH2", "CZ3", "HZ3", TRANS },
    { "CE2", "CD2", "CE3", "HE3", TRANS }, 
    { "CE2", "CZ2", "CH2", "HH2", TRANS }
  };
  int    natoms,i,j,k,kk,q,residnr,oldresidnr,nlogg;
  int    n14dist,n15dist,n16dist,n17dist,nVdist,atom_index;
  real   blen,lb,ub,angle;

  cda        = init_cdatoms(MAXLOGG,sa);
  natoms     = atoms->nr;
  n14dist    = 0;
  n15dist    = 0;
  n16dist    = 0;
  n17dist    = 0;
  nVdist     = 0;
  residnr    = -1;
  oldresidnr = -1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),"TRP") == 0) {
      oldresidnr = residnr;
      residnr    = atoms->atom[i].resnr;
      if (oldresidnr == residnr) 
	continue;
      
      nlogg = set_cdatoms(atoms,i,residnr,MAXLOGG,cda);
      
      if ( ((nlogg == MAXLOGG-1) && !bVir) || ((nlogg == MAXLOGG) && bVir) ) {
	
	pr_logg(debug,MAXLOGG-1,cda,bVir,"Trp");

	for(kk=0; (kk<NPD); kk++)
	  n14dist += do_a_pdih(logger(pd[kk].ai),logger(pd[kk].aj),
			       logger(pd[kk].ak),logger(pd[kk].al),
			       pd[kk].type,ilist,iparams,TRUE,atoms,
			       ring_margin,weight,d);
	
	for(k=0; (k<asize(cd15)); k++) {
	  trp15[cd15[k].type-1](logger(cd15[k].ai),logger(cd15[k].aj),
				logger(cd15[k].ak),logger(cd15[k].al),
				logger(cd15[k].am),
				ilist,iparams,&blen,atoms);
	  n15dist += do_a_dist(logger(cd15[k].ai),
			       logger(cd15[k].am),natoms,0.5*ring_margin,
			       weight,d,blen);
	}

	for(k=0; (k<asize(cd16)); k++) {
	  trp16[cd16[k].type-1](logger(cd16[k].ai),logger(cd16[k].aj),
				logger(cd16[k].ak),logger(cd16[k].al),
				logger(cd16[k].am),logger(cd16[k].an),
				ilist,iparams,&blen,atoms);
	  n16dist += do_a_dist(logger(cd16[k].ai),
			       logger(cd16[k].an),natoms,0.5*ring_margin,
			       weight,d,blen);
	}
	
	/* SETDISTANCE for HH2 and HD1 (1-7) */
	trp_17_type2(logger("HH2"),logger("CH2"),logger("CZ2"),logger("CE2"),
		     logger("NE1"),logger("CD1"),logger("HD1"),
		     ilist,iparams,&blen,atoms);
	n17dist += do_a_dist(logger("HH2"),logger("HD1"),
			     natoms,0.5*ring_margin,weight,d,blen);
	
	/* SETDISTANCE for HZ3 and HE1 (1-7) */
	trp_17_type1(logger("HZ3"),logger("CZ3"),logger("CE3"),logger("CD2"),
		     logger("CE2"),logger("NE1"),logger("HE1"),
		     ilist,iparams,&blen,atoms);
	n17dist += do_a_dist(logger("HZ3"),logger("HE1"),
			     natoms,0.5*ring_margin,weight,d,blen);
		     
	/* SETDISTANCE for HH2 and CB (1-7) */
	trp_17_type1(logger("HH2"),logger("CH2"),logger("CZ2"),logger("CE2"),
		     logger("CD2"),logger("CG"),logger("CB"),
		     ilist,iparams,&blen,atoms);
	n17dist += do_a_dist(logger("HH2"),logger("CB"),
			     natoms,0.5*ring_margin,weight,d,blen);

	/* SETDISTANCE for HZ3 and HD1 (1-7) */
	trp_17_type2(logger("HZ3"),logger("CZ3"),logger("CE3"),logger("CD2"),
		     logger("CG"),logger("CD1"),logger("HD1"),
		     ilist,iparams,&blen,atoms);
	n17dist += do_a_dist(logger("HZ3"),logger("HD1"),
			     natoms,0.5*ring_margin,weight,d,blen);
		     
	/* VIRTUAL DISTANCES */
	if (bVir) 
	  nVdist += set_virtual (cda,MAXLOGG,ring_margin,d,natoms);
      }
    }
  }
  pr_ndist(log,"TRP",n14dist,n15dist,n16dist,n17dist,nVdist);
  
  sfree(cda);
}
 
/**********************************************************
 *
 *     V A L I N E
 *
 **********************************************************/

void val (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	  real val_margin,t_ilist ilist[],t_iparams iparams[])
{
  static char *sa[] = { "CA",  "CB",   "HB", 
			"CG1", "HG11", "HG12", "HG13",
			"CG2", "HG21", "HG22", "HG23" };
  t_cdatom *cda;
  t_cd15g  cd15g[] = {
    { "HG11", "CG1", "CB", "CG2", "HG21", M_PI, M_PI43, M_PI43 },
    { "HG11", "CG1", "CB", "CG2", "HG22", M_PI, M_PI43, M_PI23 },
    { "HG11", "CG1", "CB", "CG2", "HG23", M_PI, M_PI43, 0 },
    { "HG12", "CG1", "CB", "CG2", "HG21", M_PI, M_PI23, M_PI43 },
    { "HG12", "CG1", "CB", "CG2", "HG22", M_PI, M_PI23, M_PI23 },
    { "HG12", "CG1", "CB", "CG2", "HG23", M_PI, M_PI23, 0 },
    { "HG13", "CG1", "CB", "CG2", "HG21", M_PI, 0,      M_PI43 },
    { "HG13", "CG1", "CB", "CG2", "HG22", M_PI, 0,      M_PI23 },
    { "HG13", "CG1", "CB", "CG2", "HG23", M_PI, 0,      0 }
  };
  t_cdpdih pd[] = { 
    { "HG11", "CG1", "CB", "CA",  GAUCHE }, 
    { "HG11", "CG1", "CB", "CG2", GAUCHE },
    { "HG12", "CG1", "CB", "CG2", GAUCHE }, 
    { "HG12", "CG1", "CB", "HB",  GAUCHE }, 
    { "HG13", "CG1", "CB", "HB",  GAUCHE }, 
    { "HG13", "CG1", "CB", "CA",  GAUCHE },
    { "HG21", "CG2", "CB", "HB",  GAUCHE }, 
    { "HG21", "CG2", "CB", "CG1", GAUCHE }, 
    { "HG22", "CG2", "CB", "CG1", GAUCHE }, 
    { "HG22", "CG2", "CB", "CA",  GAUCHE }, 
    { "HG23", "CG2", "CB", "CA",  GAUCHE }, 
    { "HG23", "CG2", "CB", "HB",  GAUCHE },
    { "HG11", "CG1", "CB", "HB",  TRANS }, 
    { "HG12", "CG1", "CB", "CA",  TRANS },
    { "HG13", "CG1", "CB", "CG2", TRANS }, 
    { "HG21", "CG2", "CB", "CA",  TRANS },
    { "HG22", "CG2", "CB", "HB",  TRANS }, 
    { "HG23", "CG2", "CB", "CG1", TRANS }
  };
  int  natoms,i,j,kk,residnr,oldresidnr,n14dist,n15dist,nlogg;
  real blen,lb,ub,angle;
  real pi = M_PI;

  cda        =  init_cdatoms(MAXLOGG,sa);
  natoms     =  atoms->nr;
  n14dist    =  0;
  n15dist    =  0;
  residnr    = -1;
  oldresidnr = -1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),"VAL") == 0) {
      oldresidnr = residnr;
      residnr    = atoms->atom[i].resnr;
      if (oldresidnr == residnr)
	continue;
      
      nlogg = set_cdatoms(atoms,i,residnr,MAXLOGG,cda);
      
      if (nlogg == MAXLOGG) {
	pr_logg(debug,MAXLOGG-1,cda,FALSE,"Val");

	/* Proper dihedrals */
	for(kk=0; (kk<NPD); kk++)
	  n14dist += do_a_pdih(logger(pd[kk].ai),logger(pd[kk].aj),
			       logger(pd[kk].ak),logger(pd[kk].al),
			       pd[kk].type,ilist,iparams,TRUE,atoms,val_margin,
			       weight,d);

	/* 1-5 dihedrals */
	for(kk=0; (kk<asize(cd15g)); kk++) {
	  gauche15(logger(cd15g[kk].ai),logger(cd15g[kk].aj),
		   logger(cd15g[kk].ak),logger(cd15g[kk].al),
		   logger(cd15g[kk].am),
		   cd15g[kk].om1,cd15g[kk].om2,cd15g[kk].om3,
		   ilist,iparams,&blen,atoms);
	  n15dist += do_a_dist(logger(cd15g[kk].ai),logger(cd15g[kk].am),
			       natoms,val_margin,weight,d,blen);
	}			
      }
    }
  }
  pr_ndist(log,"VAL",n14dist,n15dist,0,0,0);
  
  sfree(cda);
}

/**********************************************************
 *
 *     A L L   B O N D S   a n d    A N G L E S
 *
 **********************************************************/

/* Hacked by David and probably full of errors. */
void simple_bonds_and_angles(FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,
			     real weight[],real bond_margin,real angle_margin)
{
  static int bbb[] = { F_BONDS, F_SHAKE, F_G96BONDS,
		       F_CUBICBONDS, F_CONNBONDS,
		       F_MORSE };
  static int aaa[] = { F_ANGLES, F_G96ANGLES };
  
  int     i,j,ftype,nratoms;
  int     ndist12,ndist13;
  t_iatom *ia;
  
  ndist12 = ndist13 = 0;
  for(i=0; (i<asize(bbb)); i++) {
    ftype   = bbb[i];
    nratoms = interaction_function[ftype].nratoms+1;
    ia      = idef->il[ftype].iatoms;
    for(j=0; (j<idef->il[ftype].nr); ) {
      ndist12 += do_a_bond(ia[1],ia[2],idef->il,idef->iparams,
			   TRUE,atoms,bond_margin,weight,d);
      j  += nratoms;
      ia += nratoms;
    }
  }
  for(i=0; (i<asize(aaa)); i++) {
    ftype   = aaa[i];
    nratoms = interaction_function[ftype].nratoms+1;
    ia      = idef->il[ftype].iatoms;
    for(j=0; (j<idef->il[ftype].nr); ) {
      ndist13 += do_an_angle(ia[1],ia[2],ia[3],idef->il,idef->iparams,
			     TRUE,atoms,angle_margin,weight,d);
      j  += nratoms;
      ia += nratoms;
    }
  }
  fprintf(stderr,"Added %d generic bonds and %d generic angles\n",
	  ndist12,ndist13);
}

/**********************************************************
 *
 *     P E P T I D E   B O N D S
 *
 **********************************************************/

/* Hacked by Adam and probably full of errors. */
void peptide_bonds (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,
		    real weight[], real pep_margin,
		    t_ilist ilist[],t_iparams iparams[],bool bVir)
     
{
  /* Note that N, H(CD), CA are from next residue (indicated by +) */
  static char *sa[]    = { "CA", "C", "O", "N+", "H+",  "CA+", "VP" };
  static char *sapro[] = { "CA", "C", "O", "N+", "CD+", "CA+", "VP" };
  t_cdatom *cdanorm,*cdapro,*cda;
  t_cdpdih pdnorm[] = {
    { "O",  "C", "N+", "H+",   TRANS }, 
    { "O",  "C", "N+", "CA+",  CIS },
    { "CA", "C", "N+", "CA+",  TRANS }, 
    { "CA", "C", "N+", "H+",   CIS }
  };  
  t_cdpdih pdpro[] = {
    { "O",  "C", "N+", "CD+",  TRANS }, 
    { "O",  "C", "N+", "CA+",  CIS },
    { "CA", "C", "N+", "CA+",  TRANS }, 
    { "CA", "C", "N+", "CD+",  CIS }
  };  
  t_cdpdih *pd;
  int  natoms,odist,i,j,kk,q,nlogg,nloggpro,maxlogg;
  int  n14dist,nVdist,residnr,oldresidnr;
  real blen,lb,ub,angle;
  bool bPro;
  
  cdanorm  =  init_cdatoms(MAXLOGG,sa);
  cdapro   =  init_cdatoms(asize(sapro),sapro);  
  natoms   =  atoms->nr;
  n14dist  =  0;
  nVdist   =  0;
  residnr  = -1;

  for (i=0; (i<natoms); i++) {
    oldresidnr = residnr;
    residnr    = atoms->atom[i].resnr;
    if (oldresidnr == residnr)
      continue;
  
    /* Just search for atoms twice, if the next residue is a proline
     * the search should find one more atom... 
     */
    nloggpro = set_cdatoms(atoms,i,residnr,MAXLOGG,cdapro);
    nlogg    = set_cdatoms(atoms,i,residnr,MAXLOGG,cdanorm);
    maxlogg  = max(nlogg,nloggpro);
    bPro     = (nloggpro > nlogg);
      
    if ((bVir && (maxlogg == MAXLOGG)) || (!bVir && (maxlogg == MAXLOGG-1))) {
      if (bPro) {
	pd  = pdpro;
	cda = cdapro;
      }
      else {
	pd  = pdnorm;
	cda = cdanorm;
      }
      pr_logg(debug,MAXLOGG-1,cda,bVir,"PEP");
      
      for(kk=0; (kk<asize(pdnorm)); kk++)
	n14dist += do_a_pdih(logger(pd[kk].ai),logger(pd[kk].aj),
			     logger(pd[kk].ak),logger(pd[kk].al),
			     pd[kk].type,ilist,iparams,TRUE,atoms,
			     pep_margin,weight,d);
      
      /* VIRTUAL DISTANCES */
      if (bVir) 
	nVdist += set_virtual (cda,MAXLOGG-1,pep_margin,d,natoms);
    }
  }
  pr_ndist(log,"PEP",n14dist,0,0,0,nVdist);
  
  sfree(cdanorm);
  sfree(cdapro);
}
