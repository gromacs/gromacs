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
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include "typedefs.h"
#include "pdbio.h"
#include "smalloc.h"
#include "vec.h"
#include "physics.h"
#include "toputil.h"
#include "pdb2top.h"
#include "string2.h"

static int in_strings(char *key,int nstr,const char **str)
{
  int j;
  
  for(j=0; (j<nstr); j++)
    if (strcmp(str[j],key) == 0)
      return j;
      
  return -1;
}

static gmx_bool hbond(rvec x[],int i,int j,real distance)
{
  real tol = distance*distance;
  rvec   tmp;
  
  rvec_sub(x[i],x[j],tmp);
  
  return (iprod(tmp,tmp) < tol);
}

static void chk_allhb(t_atoms *pdba,rvec x[],t_blocka *hb,
		      gmx_bool donor[],gmx_bool accept[],real dist)
{
  int i,j,k,ii,natom;
  
  natom=pdba->nr;
  snew(hb->index,natom+1);
  snew(hb->a,6*natom);
  hb->nr  = natom;
  hb->nra = 6*natom;
  
  k = ii = 0;
  hb->index[ii++] = 0;
  for(i=0; (i<natom); i++) {
    if (donor[i]) {
      for(j=i+1; (j<natom); j++) 
	if ((accept[j]) && (hbond(x,i,j,dist))) 
	  hb->a[k++] = j;
    }
    else if (accept[i]) {
      for(j=i+1; (j<natom); j++) 
	if ((donor[j]) && (hbond(x,i,j,dist))) 
	  hb->a[k++] = j;
    }
    hb->index[ii++] = k;
  }
  hb->nra = k;
}

static void pr_hbonds(FILE *fp,t_blocka *hb,t_atoms *pdba)
{
  int i,j,k,j0,j1;
  
  fprintf(fp,"Dumping all hydrogen bonds!\n");
  for(i=0; (i<hb->nr); i++) {
    j0=hb->index[i];
    j1=hb->index[i+1];
    for(j=j0; (j<j1); j++) {
      k=hb->a[j];
      fprintf(fp,"%5s%4d%5s - %5s%4d%5s\n",
	      *pdba->resinfo[pdba->atom[i].resind].name,
	      pdba->resinfo[pdba->atom[i].resind].nr,*pdba->atomname[i],
	      *pdba->resinfo[pdba->atom[k].resind].name,
	      pdba->resinfo[pdba->atom[k].resind].nr,*pdba->atomname[k]);
    }
  }
}

static gmx_bool chk_hbonds(int i,t_atoms *pdba, rvec x[],
		       gmx_bool ad[],gmx_bool hbond[],rvec xh,
		       real angle,real dist)
{
  gmx_bool bHB;
  int  j,aj,ri,natom;
  real d2,dist2,a;
  rvec nh,oh;
  
  natom=pdba->nr;
  bHB = FALSE;
  ri = pdba->atom[i].resind;
  dist2=sqr(dist);
  for(j=0; (j<natom); j++) {
    /* Check whether the other atom is a donor/acceptor and not i */
    if ((ad[j]) && (j != i)) {
      /* Check whether the other atom is on the same ring as well */
      if ((pdba->atom[j].resind != ri) ||
	  ((strcmp(*pdba->atomname[j],"ND1") != 0) &&
	   (strcmp(*pdba->atomname[j],"NE2") != 0))) {
	aj = j;
	d2  = distance2(x[i],x[j]);
	rvec_sub(x[i],xh,nh);
	rvec_sub(x[aj],xh,oh);
	a  = RAD2DEG * acos(cos_angle(nh,oh));
	if ((d2 < dist2) && (a > angle)) {
	  if (debug)
	    fprintf(debug,
		    "HBOND between %s%d-%s and %s%d-%s is %g nm, %g deg\n",
		    *pdba->resinfo[pdba->atom[i].resind].name, 
		    pdba->resinfo[pdba->atom[i].resind].nr,*pdba->atomname[i],
		    *pdba->resinfo[pdba->atom[aj].resind].name,
		    pdba->resinfo[pdba->atom[aj].resind].nr,*pdba->atomname[aj],
		    sqrt(d2),a);
	  hbond[i] = TRUE;
	  bHB      = TRUE;
	}
      }
    }
  }
  return bHB;
}

static void calc_ringh(rvec xattach,rvec xb,rvec xc,rvec xh)
{
  rvec tab,tac;
  real n;
 
  /* Add a proton on a ring to atom attach at distance 0.1 nm */ 
  rvec_sub(xattach,xb,tab);
  rvec_sub(xattach,xc,tac);
  rvec_add(tab,tac,xh);
  n=0.1/norm(xh);
  svmul(n,xh,xh);
  rvec_inc(xh,xattach);
}

void set_histp(t_atoms *pdba,rvec *x,real angle,real dist){
  static const char *prot_acc[] = {
    "O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "OW"
  };
#define NPA asize(prot_acc)
  static const char *prot_don[] = {
    "N", "NH1", "NH2", "NE", "ND1", "ND2", "NE2", "NZ", "OG", "OG1", "OH", "NE1", "OW"
  };
#define NPD asize(prot_don)
  
  gmx_bool *donor,*acceptor;
  gmx_bool *hbond,bHaveH=FALSE;
  gmx_bool bHDd,bHEd;
  rvec xh1,xh2;
  int  natom;
  int  i,nd,na,aj,hisind,his0,type=-1;
  int  nd1,ne2,cg,cd2,ce1;
  t_blocka *hb;
  real d;
  char *atomnm;
  
  natom=pdba->nr;
  snew(donor,natom);
  snew(acceptor,natom);
  snew(hbond,natom);
  snew(hb,1);
  
  nd=na=0;
  for(i=0; (i<natom); i++) {
    if (in_strings(*pdba->atomname[i],NPA,prot_acc) != -1) {
      acceptor[i] = TRUE;
      na++;
    }
    if (in_strings(*pdba->atomname[i],NPD,prot_don) != -1) {
      donor[i] = TRUE;
      nd++;
    }
  }
  fprintf(stderr,"There are %d donors and %d acceptors\n",nd,na);
  chk_allhb(pdba,x,hb,donor,acceptor,dist);
  if (debug)
    pr_hbonds(debug,hb,pdba);
  fprintf(stderr,"There are %d hydrogen bonds\n",hb->nra);
  
  /* Now do the HIS stuff */
  hisind=-1;
  for(i=0; (i<natom); ) {
    if (gmx_strcasecmp(*pdba->resinfo[pdba->atom[i].resind].name,"HIS") != 0) 
      i++;
    else {
      if (pdba->atom[i].resind != hisind) {
	hisind=pdba->atom[i].resind;
	
	/* Find the  atoms in the ring */
	nd1=ne2=cg=cd2=ce1=-1;
	while (i<natom && pdba->atom[i].resind==hisind) {
	  atomnm = *pdba->atomname[i];
	  if (strcmp(atomnm,"CD2") == 0)
	    cd2 = i;
	  else if (strcmp(atomnm,"CG") == 0)
	    cg  = i;
	  else if (strcmp(atomnm,"CE1") == 0)
	    ce1 = i;
	  else if (strcmp(atomnm,"ND1") == 0)
	    nd1 = i;
	  else if (strcmp(atomnm,"NE2") == 0)
	    ne2 = i;

	  i++;
	}
	
	if (!((cg == -1 ) || (cd2 == -1) || (ce1 == -1) ||
	      (nd1 == -1) || (ne2 == -1))) {
	  calc_ringh(x[nd1],x[cg],x[ce1],xh1);
	  calc_ringh(x[ne2],x[ce1],x[cd2],xh2);
	  
	  bHDd = chk_hbonds(nd1,pdba,x,acceptor,hbond,xh1,angle,dist);
	  chk_hbonds(nd1,pdba,x,donor,hbond,xh1,angle,dist);
	  bHEd = chk_hbonds(ne2,pdba,x,acceptor,hbond,xh2,angle,dist);
	  chk_hbonds(ne2,pdba,x,donor,hbond,xh2,angle,dist);
	  
	  if (bHDd) {
	    if (bHEd)
	      type = ehisH;
	    else 
	      type = ehisA;
	  }
	  else 
	    type = ehisB;
	  fprintf(stderr,"Will use %s for residue %d\n",
		  hh[type],pdba->resinfo[hisind].nr);
	}
	else 
	  gmx_fatal(FARGS,"Incomplete ring in HIS%d",
		    pdba->resinfo[hisind].nr);
	
	snew(pdba->resinfo[hisind].rtp,1);
	*pdba->resinfo[hisind].rtp = strdup(hh[type]);
      }
    }
  }
  done_blocka(hb);
  sfree(hb);
  sfree(donor);
  sfree(acceptor);
  sfree(hbond);
}
