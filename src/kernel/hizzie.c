/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
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
 * Giant Rising Ordinary Mutants for A Clerical Setup
 */
static char *SRCID_hizzie_c = "$Id$";

#include <stdio.h>
#include <string.h>
#include "typedefs.h"
#include "pdbio.h"
#include "smalloc.h"
#include "vec.h"
#include "physics.h"
#include "toputil.h"
#include "pdb2gmx.h"
#include "pdb2top.h"

static int in_strings(char *key,int nstr,char **str)
{
  int j;
  
  for(j=0; (j<nstr); j++)
    if (strcmp(str[j],key) == 0)
      return j;
      
  return -1;
}

static bool hbond(rvec x[],int i,int j,real distance)
{
  real tol = distance*distance;
  rvec   tmp;
  
  rvec_sub(x[i],x[j],tmp);
  
  return (iprod(tmp,tmp) < tol);
}

static void chk_allhb(t_atoms *pdba,rvec x[],t_block *hb,
		      bool donor[],bool accept[],real dist)
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

static void pr_hbonds(FILE *fp,t_block *hb,t_atoms *pdba)
{
  int i,j,k,j0,j1;
  
  fprintf(fp,"Dumping all hydrogen bonds!\n");
  for(i=0; (i<hb->nr); i++) {
    j0=hb->index[i];
    j1=hb->index[i+1];
    for(j=j0; (j<j1); j++) {
      k=hb->a[j];
      fprintf(fp,"%5s%4d%5s - %5s%4d%5s\n",
	      *pdba->resname[pdba->atom[i].resnr],
	      pdba->atom[i].resnr+1,*pdba->atomname[i],
	      *pdba->resname[pdba->atom[k].resnr],
	      pdba->atom[k].resnr+1,*pdba->atomname[k]);
    }
  }
}

static bool chk_hbonds(int i,t_atoms *pdba, rvec x[],
		       bool ad[],bool hbond[],rvec xh,
		       real angle,real dist)
{
  bool bHB;
  int  j,aj,ri,natom;
  real d2,dist2,a;
  rvec nh,oh;
  
  natom=pdba->nr;
  bHB = FALSE;
  ri = pdba->atom[i].resnr;
  dist2=sqr(dist);
  for(j=0; (j<natom); j++) {
    /* Check whether the other atom is a donor/acceptor and not i */
    if ((ad[j]) && (j != i)) {
      /* Check whether the other atom is on the same ring as well */
      if ((pdba->atom[j].resnr != ri) ||
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
		    *pdba->resname[pdba->atom[i].resnr], 
		    pdba->atom[i].resnr+1, *pdba->atomname[i],
		    *pdba->resname[pdba->atom[aj].resnr],
		    pdba->atom[aj].resnr+1,*pdba->atomname[aj],sqrt(d2),a);
	  hbond[i] = TRUE;
	  bHB      = TRUE;
	}
      }
    }
  }
  return bHB;
}

void calc_ringh(rvec xattach,rvec xb,rvec xc,rvec xh)
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
  static char *prot_acc[] = {
    "O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "OW"
  };
#define NPA asize(prot_acc)
  static char *prot_don[] = {
    "N", "NH1", "NH2", "NE", "ND1", "ND2", "NE2", "NZ", "OG", "OG1", "OH", "NE1", "OW"
  };
#define NPD asize(prot_don)
  
  bool *donor,*acceptor;
  bool *hbond,bHaveH=FALSE;
  bool bHDd,bHEd;
  rvec xh1,xh2;
  int  natom;
  int  i,j,nd,na,aj,hisnr,his0,type;
  int  nd1,ne2,cg,cd2,ce1;
  t_block *hb;
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
  hisnr=-1;
  for(i=0; (i<natom); ) {
    if (strcasecmp(*pdba->resname[pdba->atom[i].resnr],"HIS") != 0) 
      i++;
    else {
      if (pdba->atom[i].resnr != hisnr) {
	hisnr=pdba->atom[i].resnr;
	
	/* Find the  atoms in the ring */
	nd1=ne2=cg=cd2=ce1=-1;
	for(j=i; (pdba->atom[j].resnr==hisnr) && (j<natom); j++) {
	  atomnm=*pdba->atomname[j];
	  if (strcmp(atomnm,"CD2") == 0)
	    cd2=j;
	  else if (strcmp(atomnm,"CG") == 0)
	    cg=j;
	  else if (strcmp(atomnm,"CE1") == 0)
	    ce1=j;
	  else if (strcmp(atomnm,"ND1") == 0)
	    nd1=j;
	  else if (strcmp(atomnm,"NE2") == 0)
	    ne2=j;
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
	  fprintf(stderr,"Will use %s for residue %d\n",hh[type],hisnr+1);
	}
	else 
	  fatal_error(0,"Incomplete ring in HIS%d",hisnr+1);
	
	sfree(*pdba->resname[hisnr]);
	*pdba->resname[hisnr]=strdup(hh[type]);
      }
    }
  }
  done_block(hb);
  sfree(hb);
  sfree(donor);
  sfree(acceptor);
  sfree(hbond);
}
