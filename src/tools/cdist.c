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

#include <stdlib.h>
#include <ctype.h>
#include "macros.h"
#include "vec.h"
#include "fatal.h"
#include "txtdump.h"
#include "cdist.h"
#include "invblock.h"
#include "futil.h"
#include "tpxio.h"
#include "index.h"

#define NINDEX(ai,aj,natom) (natom*(ai)-((ai)*(ai+1))/2+(aj)-(ai))
#define INDEX(ai,aj,natom) ((ai) < (aj))?NINDEX((ai),(aj),natom):NINDEX(aj,ai,natom) 
#define CHECK(ai,natom)    if (((ai)<0) || ((ai)>=natom)) gmx_fatal(FARGS,"Invalid atom number %d",ai+1)

static real vdwlen(t_atoms *atoms,int i,int j)
{
#define NAT 5
  char anm[NAT] = "HCNOS";
  /* For time being I give S the same vdw-parameters as O */
  real dist[NAT][NAT] = { 
    {  2.0,  2.4,  2.4,  2.3,  2.3 },
    {  2.4,  3.0,  2.9, 2.75, 2.75 },
    {  2.4,  2.9,  2.7,  2.7,  2.7 },
    {  2.3, 2.75,  2.7,  2.8,  2.8 },
    {  2.3, 2.75,  2.7,  2.8,  2.8 }
  };
  char ati,atj;
  int  ai,aj;
  
  ati=toupper((*atoms->atomname[i])[0]);
  atj=toupper((*atoms->atomname[j])[0]);
  
  for(ai=0; (ai<NAT) && (ati != anm[ai]); ai++)
    ;
  for(aj=0; (aj<NAT) && (atj != anm[aj]); aj++)
    ;
  if ((ai < NAT) && (aj < NAT))
    return dist[ai][aj];
  else {
    fprintf(stderr,"No combination for nbs (%c %c) using 2.0 A\n",ati,atj);
    return 2.0;
  }
}
  
void set_dist(t_dist *d,int natoms,int ai,int aj,real lb,
	      real ub,real len)
{
  int index;

  if (ub <= lb)
    fprintf(stderr,"set_dist: lb = %f, ub = %f, len = %f, atoms %d,%d\n",
	    lb,ub,len,ai,aj); 
  else {
    CHECK(ai,natoms);
    CHECK(aj,natoms);
    index = INDEX(ai,aj,natoms);
    d[index].lb =  lb;
    d[index].ub =  ub;
    d[index].len = len;
  }
}

void set_ideal(t_dist *d,int natoms,int ai,int aj,real len)
{
  int index;

  CHECK(ai,natoms);
  CHECK(aj,natoms);
  index = INDEX(ai,aj,natoms);
  d[index].len = len;
}

bool dist_set(t_dist *d,int natoms,int ai,int aj)
{
  int index;
  
  index = INDEX(ai,aj,natoms);
  
  return (d[index].lb != 0.0);
}

static t_dist *new_dist(int natom)
{
  t_dist *d;
  
  snew(d,(natom*(natom+1))/2);

  return d;
}

real d_lb(t_dist *d,int natoms,int ai,int aj)
{
  int index;
  
  index = INDEX(ai,aj,natoms);
  
  return d[index].lb;
}

real d_ub(t_dist *d,int natoms,int ai,int aj)
{
  int index;
  
  index = INDEX(ai,aj,natoms);
  
  return d[index].ub;
}

real d_len(t_dist *d,int natoms,int ai,int aj)
{
  int index;
  
  index = INDEX(ai,aj,natoms);
  
  return d[index].len;
}

static t_dist *read_dist(FILE *log,char *fn,int natom,real weight[])
{
#define  BLEN 255
  FILE   *fp;
  char   buf[BLEN+1];
  int    i,ai,aj,www,ndist;
  double lb,ub,len;
  t_dist *d;
  
  fprintf(log,"Going to read %s\n",fn);
  d     = new_dist(natom);
  fp    = ffopen(fn,"r");
  ndist = 0;
  www   = 0;
  while (fgets2(buf,BLEN,fp) != NULL) {
    if (buf[0] != '#') {
      if (sscanf(buf,"%d%d%lf%lf%lf",&ai,&aj,&len,&lb,&ub) != 5)
	if (sscanf(buf,"%d%d%lf%lf",&ai,&aj,&lb,&ub) != 4)
	  gmx_fatal(FARGS,"Invalid dist format in %s",fn);
      ai--;
      aj--;
      if ((weight[ai] != 0) || (weight[aj] != 0)) {
	if (!dist_set(d,natom,ai,aj)) {
	  if (debug)
	    fprintf(debug,"\r%d %d %d",ndist,ai,aj);
	  set_dist(d,natom,ai,aj,lb,ub,-1);
	  ndist++;
	}
      }
      else
	www++;
    }
  }
  fclose(fp);
  
  fprintf(stderr,
	  "Read %d distances from %s (discarded %d due to zero weight)\n",
	  ndist,fn,www);
  fprintf(log,
	  "Read %d distances from %s (discarded %d due to zero weight)\n",
	  ndist,fn,www);
  fflush(log);
  
  return d;
}

static void measure_report(FILE *fp,int nm,int natom,int nover,int nnotideal,
			   int nhb_rep,int nhb,real cutoff,real rmsd)
{
  double ndtot=(natom*(natom-1.0))/2.0;
  
  fprintf(fp,"Determined %d distances out of a possible %.0f "
	  "(%.2f %%)\nby measuring from the structure within %f A\n",
	  nm,ndtot,(nm*100.0)/ndtot,cutoff);
  fprintf(fp,"%d distances from the database were replaced by measured ones\n",
	  nover);
  fprintf(fp,"Of these, %d were not within the bounds of the database.\n"
	  "RMS deviation from bound limits for these was %e A.\n",
	  nnotideal,rmsd);
  fprintf(fp,"For %d hydrogen bonds out of %d in the index file the margins " 
	  "were reduced\n",nhb_rep/2,nhb/3);
  fflush(fp);
}

static void measure_dist(FILE *log,t_dist *d,t_atoms *atoms,rvec x[],
			 real cutoff,real margin,real hbmargin,
			 int nhb,atom_id hb[])
{
  int    i,j,natom,dd,aa,ai,aj,nm,nover,nnotideal,nhb_rep;
  real   msd,rmsd;
  rvec   dx;
  real   ideal,lbfac,ubfac,vdw,lb,ub,nmargin;
  
  nm        = 0;
  nover     = 0;
  nnotideal = 0;
  msd       = 0;
  natom     = atoms->nr;
  
  lbfac     = 1-margin;
  ubfac     = 1+margin;

  for(ai=0; (ai<natom); ai++)
    for(aj=ai+1; (aj<natom); aj++) {
      rvec_sub(x[ai],x[aj],dx);
      ideal = 10*norm(dx);
      if (ideal == 0.0) {
	fprintf(stderr,"Warning distance between atoms %s and %s is zero\n",
		atomname(atoms,ai),atomname(atoms,aj));
      }
      else {
	if (!dist_set(d,natom,ai,aj)) {
	  /* This distance is not determined by the database. */
	  vdw = 0; /*vdwlen(atoms,ai,aj);*/
	  if ((ideal < cutoff)  || (cutoff == 0)) {
	    set_dist(d,natom,ai,aj,max(vdw,ideal*lbfac),ideal*ubfac,ideal);
	    nm++;
	  }
	}
	else {
	  /* These distances are already set by the database routines.
	   * However, we override the distances with the measured ones
	   * while keeping the original margins.
	   */
	  lb = d_lb(d,natom,ai,aj);
	  ub = d_ub(d,natom,ai,aj);
	  if ((ideal < lb) || (ideal > ub)) {
	    if (debug)
	      fprintf(debug,"Warning: d(%s,%s) = %8.4f. According to E&H"
		      " it should be %8.4f (dev. %.1f%%)\n",
		      atomname(atoms,ai),
		      atomname(atoms,aj),
		      ideal,(lb+ub)*0.5,50*fabs(ideal*2-lb-ub));
	    msd += (ideal < lb) ? sqr(ideal-lb) : sqr(ideal-ub);
	    nnotideal++;
	  }
	  nmargin = (ub-lb)/(ub+lb);
	  set_dist(d,natom,ai,aj,ideal*(1-nmargin),ideal*(1+nmargin),ideal);
	  nm++;
	  nover++;
	}
      }
    }
  /* Now we have set all the distances. We have to verify though
   * that h bonds are maintained, we therefore reduce their margins.
   */
  nhb_rep = 0;
  for(i=0; (i<nhb); i+= 3) {
    /* Acceptor atom */
    aa = hb[i+2];
    range_check(aa,0,natom);
    for(j=0; (j<2); j++) {
      /* Donor atom */
      dd      = hb[i+j];
      range_check(dd,0,natom);
      if (dist_set(d,natom,dd,aa)) {
	lb      = d_lb(d,natom,dd,aa);
	ub      = d_ub(d,natom,dd,aa);
	ideal   = d_len(d,natom,dd,aa);
	nmargin = (ub-lb)/(ub+lb);
	if (hbmargin < nmargin) {
	  set_dist(d,natom,dd,aa,ideal*(1-hbmargin),ideal*(1+hbmargin),ideal);
	  nhb_rep++;
	}
      }
      else
	gmx_fatal(FARGS,"Distance between %s and %s not set, while they do make"
		    " a hbond:\nincrease radius for measuring (now %f A)\n",
		    atomname(atoms,aa),atomname(atoms,dd),cutoff);
    }
  }
	
  if (nnotideal > 0)
    rmsd = sqrt(msd/nnotideal);
  else
    rmsd = 0;
  measure_report(log,nm,natom,nover,nnotideal,nhb_rep,nhb,cutoff,rmsd);
  measure_report(stderr,nm,natom,nover,nnotideal,nhb_rep,nhb,cutoff,rmsd);
}

static void dump_dist(char *fn,t_dist *d,int natom,bool bAddR)
{
  FILE *fp;
  int  i,j;
  real lb,ub,len;
  
  fp    = ffopen(fn,"w");
  fprintf(fp,"# distance file generated by cdist\n");
  
  for(i=0; (i<natom); i++) {
    for(j=i+1; (j<natom); j++) {
      if (dist_set(d,natom,i,j)) {
	lb  = d_lb(d,natom,i,j);
	ub  = d_ub(d,natom,i,j);
	len = d_len(d,natom,i,j);
	if (bAddR)
	  fprintf(fp,"%5d%5d%10.5f%10.5f\n",i+1,j+1,lb,ub);
	else
	  fprintf(fp,"%5d%5d%10.5f%10.5f%10.5f\n",i+1,j+1,len,lb,ub); 
	  /* Output for real-precision disco: */
	  /* fprintf(fp,"%5d%5d%15.10f%15.10f%15.10f\n",i+1,j+1,len,lb,ub); */ 
      }
    }
  }
  fclose(fp);
}

static real *read_weights(char *fn,int natom)
{
  FILE    *in;
  int     i,n;
  char    title[STRLEN];
  real    *w;
  matrix  box;
  t_atoms atoms;
  rvec    *x;
  
  /* Open file */
  in = ffopen(fn,"r");
  
  /* Check the number of atoms */
  get_pdb_coordnum(in,&n);
  if (n != natom)
    gmx_fatal(FARGS,"Number of atoms in pdb file (%d) does not match tpx (%d)",
		n,natom);
  
  /* Allocate space */
  snew(w,natom);
  snew(x,natom);
  init_t_atoms(&atoms,natom,TRUE);
  clear_mat(box);

  /* Now read it all */  
  rewind(in);
  read_pdbfile(in,title,NULL,&atoms,x,box,FALSE);
  fclose(in);
  fprintf(stderr,"%s\n",title);
  
  /* Now copy the weights */
  for(i=0; (i<natom); i++)
    w[i] = atoms.pdbinfo[i].occup;
    
  /* Now free temps */
  sfree(x);
  free_t_atoms(&atoms);
  
  /* And go back */
  return w;
}

void init_rot_matrix(real mat[3][3],real theta, real omega)
{
  mat[0][0] = -(cos(theta)*cos(omega));
  mat[1][0] = -(cos(theta)*sin(omega));
  mat[2][0] = sin(theta);

  mat[0][1] = sin(omega);
  mat[1][1] = -(cos(omega));
  mat[2][1] = 0.0;
  
  mat[0][2] = sin(theta)*cos(omega);
  mat[1][2] = sin(theta)*sin(omega);
  mat[2][2] = cos(theta);
}

void vect_matrix(real vect[3], real mat[3][3]) 
{

  real tmp[3];
  int i,j;

  tmp[0] = 0.0;
  tmp[1] = 0.0;
  tmp[2] = 0.0;

  for ( i=0 ; i <= 2 ; i++ ) {
    for ( j=0 ; j <=2 ; j++ ) {
      tmp[i]=tmp[i]+mat[i][j]*vect[j];
    }
  }
  for ( i=0 ; i <=2 ; i++ ) {
    vect[i]=tmp[i];
  }
}

void gauche_(int ai,int aj,int ak,int al,t_ilist ilist[],
	     t_iparams iparams[],real *lb,t_atoms *atoms,char *file,int line)
{

  /* Matrix based approach */
  int    i,j;
  real matrix1[3][3],matrix2[3][3];
  real vect1[3],vect2[3],vect3[3];
  real dist = 0.0;
  real pi      = M_PI;
  real half_pi = M_PI*0.5;
  real rij,rjk,rkl;
  real thijk,thjkl,theta1,theta2,omega1=pi,omega2 = pi+pi/3.0;

  rij    = lookup_bondlength_(ai,aj,ilist,iparams,TRUE,atoms,file,line);
  rjk    = lookup_bondlength_(aj,ak,ilist,iparams,TRUE,atoms,file,line);
  rkl    = lookup_bondlength_(ak,al,ilist,iparams,TRUE,atoms,file,line);
  thijk  = lookup_angle_(ai,aj,ak,ilist,iparams,atoms,file,line);
  thjkl  = lookup_angle_(aj,ak,al,ilist,iparams,atoms,file,line);
  theta1 = pi-thijk;
  theta2 = pi-thjkl;

  /*Initialise vectors*/
  vect1[0]=0.0;
  vect1[1]=0.0;
  vect1[2]=rij;
  vect2[0]=0.0;
  vect2[1]=0.0;
  vect2[2]=rjk;
  vect3[0]=0.0;
  vect3[1]=0.0;
  vect3[2]=rkl;

  /*Initialize matrices*/
  init_rot_matrix(matrix1,theta1,omega1);
  init_rot_matrix(matrix2,theta2,omega2);

  /* Express vect2 and vect3 in coord. system of vect1 */
  vect_matrix(vect2,matrix1);
  vect_matrix(vect3,matrix2);
  vect_matrix(vect3,matrix1);

  /* Add vectors */
  for ( i=0 ; i <=2 ; i++) {
    vect3[i] = vect1[i]+vect2[i]+vect3[i];
  }

  /* Calculate distance */
  *lb = sqrt(vect3[0]*vect3[0]+vect3[1]*vect3[1]+vect3[2]*vect3[2]);
}


void gauche15_(int ai,int aj,int ak,int al,int am,real omega1,real omega2,
	       real omega3,t_ilist ilist[],t_iparams iparams[],real *lb,
	       t_atoms *atoms,char *file,int line)
{

  /* Matrix based approach */
  int    i,j;
  real matrix1[3][3],matrix2[3][3],matrix3[3][3];
  real vect1[3],vect2[3],vect3[3],vect4[3];
  real pi      = M_PI;
  real half_pi = M_PI*0.5;
  real rij,rjk,rkl,rlm;
  real thijk,thjkl,thklm,theta1,theta2,theta3;

  rij    = lookup_bondlength_(ai,aj,ilist,iparams,TRUE,atoms,file,line);
  rjk    = lookup_bondlength_(aj,ak,ilist,iparams,TRUE,atoms,file,line);
  rkl    = lookup_bondlength_(ak,al,ilist,iparams,TRUE,atoms,file,line);
  rlm    = lookup_bondlength_(al,am,ilist,iparams,TRUE,atoms,file,line);
  thijk  = lookup_angle_(ai,aj,ak,ilist,iparams,atoms,file,line);
  thjkl  = lookup_angle_(aj,ak,al,ilist,iparams,atoms,file,line);
  thklm  = lookup_angle_(ak,al,am,ilist,iparams,atoms,file,line);
  theta1 = pi-thijk;
  theta2 = pi-thjkl;
  theta3 = pi-thklm;

  /*Initialise vectors*/
  vect1[0]=0.0;
  vect1[1]=0.0;
  vect1[2]=rij;
  vect2[0]=0.0;
  vect2[1]=0.0;
  vect2[2]=rjk;
  vect3[0]=0.0;
  vect3[1]=0.0;
  vect3[2]=rkl;
  vect4[0]=0.0;
  vect4[1]=0.0;
  vect4[2]=rlm;

  /*Initialize matrices*/
  init_rot_matrix(matrix1,theta1,omega1);
  init_rot_matrix(matrix2,theta2,omega2);
  init_rot_matrix(matrix3,theta3,omega3);

  /* Express vect2 and vect3 in coord. system of vect1 */
  vect_matrix(vect2,matrix1);
  vect_matrix(vect3,matrix2);
  vect_matrix(vect3,matrix1);
  vect_matrix(vect4,matrix3);
  vect_matrix(vect4,matrix2);
  vect_matrix(vect4,matrix1);

  /* Add vectors */
  for ( i=0 ; i <=2 ; i++) {
    vect4[i] = vect1[i]+vect2[i]+vect3[i]+vect4[i];
  }

  /* Calculate distance */
  *lb = sqrt(vect4[0]*vect4[0]+vect4[1]*vect4[1]+vect4[2]*vect4[2]);
}


static void dump_bonds(t_atoms *atoms,t_dist *d,
		       t_ilist ilist[],t_functype functype[],
		       t_iparams ip[],real bond_margin,real angle_margin,
		       real dih_margin,real idih_margin,
		       real weight[],bool bAddR)
{
  int  dodist[] = { F_BONDS, F_MORSE, F_SHAKE, F_G96BONDS,
		    F_G96ANGLES, F_ANGLES, F_IDIHS /*,F_RBDIHS, F_PDIHS*/ };
  int  i,j,atom1,atom2,ai,aj,ak,al,type,ftype,nra,ndist,odist,natoms;
  real blen,rij,rjk,c6,c12,lb=-1,ub,theta;
  
  natoms= atoms->nr;
  ndist = 0;
  odist = 0;
  
  for(j=0; (j<asize(dodist)); j++) {
    ftype = dodist[j];
    nra   = interaction_function[ftype].nratoms;
    for(i=0; (i<ilist[ftype].nr); i+=nra+1) {
      type  = ilist[ftype].iatoms[i];
      ai    = ilist[ftype].iatoms[i+1];
      aj    = ilist[ftype].iatoms[i+2];
      atom1 = ai;
      atom2 = aj;
      blen  = 0;
      switch (ftype) {
      case F_BONDS:
      case F_MORSE:
      case F_SHAKE:
	blen  = lookup_bondlength(ai,aj,ilist,ip,TRUE,atoms);
	break;
      case F_G96ANGLES:
      case F_ANGLES:
	ak    = ilist[ftype].iatoms[i+3];
	theta = lookup_angle(ai,aj,ak,ilist,ip,atoms);
	blen  = angle_length(ai,aj,ak,RAD2DEG*theta,ilist,ip,atoms);
	atom2 = ak;
	break;
      case F_PDIHS:
      case F_RBDIHS:
	ak    = ilist[ftype].iatoms[i+3];
	al    = ilist[ftype].iatoms[i+4];
	pdih_lengths(ai,aj,ak,al,ilist,ip,&blen,&ub,atoms);
	atom2 = al;
	break;
      case F_IDIHS:
	ak    = ilist[ftype].iatoms[i+3];
	al    = ilist[ftype].iatoms[i+4];
	blen  = idih_lengths(ai,aj,ak,al,ilist,ip,atoms);
	atom2 = al;
	break;
      default:
	break;
      }
      if ((blen != 0) && ((weight[atom1] != 0.0) || (weight[atom2] != 0.0))) {
	switch (ftype) {
	case F_BONDS:
	case F_MORSE:
	case F_SHAKE:
	case F_G96BONDS:
	  lb = (1.0-bond_margin)*blen;
	  ub = (1.0+bond_margin)*blen;
	  break;
	case F_ANGLES:
	case F_G96ANGLES:
	  lb = (1.0-angle_margin)*blen;
	  ub = (1.0+angle_margin)*blen;
	  break;
	case F_IDIHS:
	  lb = (1-idih_margin)*blen;
	  ub = (1+idih_margin)*blen;
	  break;
	case F_PDIHS:
	case F_RBDIHS:
	  lb = (1.0-dih_margin)*blen;
	  ub = (1.0+dih_margin)*ub;
	  break;
	}
	if (!dist_set(d,natoms,atom1,atom2)) {
	  set_dist(d,natoms,atom1,atom2,lb,ub,blen);
	  ndist ++;
	}
	else
	  odist ++;
      }
    }
  }  
  fprintf(stderr,"There are %d new bonded distances + %d old\n",ndist,odist);
}

static bool notexcl(t_block *excl,int ai,int aj)
{
  int  i;
  
  for(i=excl->index[ai]; (i<excl->index[ai+1]); i++) {
    if (aj == excl->a[i])
      return FALSE;
  }
  return TRUE;    
}

static void dump_nonbonds(t_dist *d,t_idef *idef,t_atoms *atoms,
			  real hblen,real weight[],
			  real nb_margin,bool bAddR,
			  real maxdist)
{
  int  i,j,k,natoms,ntype,tpi,tpj,ndist,odist;
  real **dmat,len,lb;
  bool *bDA;
  char element,ati,atj;
  
  ntype   = idef->atnr;
  natoms  = atoms->nr;
  
  /* Determine which atoms are capable of H bonding or not */
  snew(bDA,ntype);
  for(j=0; (j<ntype); j++)
    for(i=0; (i<natoms); i++) {
      if (atoms->atom[i].type == j) {
	element = (*atoms->atomname[i])[0];
	bDA[j] = ((element == 'N') || (element == 'O'));
	break;
      }
    }
  
  /* Make a distance matrix for VDWaals distances */  
  /*snew(dmat,ntype);
  for(i=0; (i<ntype); i++) {
    snew(dmat[i],ntype);
    for(j=0; (j<ntype); j++) {
      if (bDA[i] && bDA[j])
	dmat[i][j] = hblen;
      else 
        dmat[i][j] = vdwlen(idef,i,j,hblen); 
    }
  } */
  sfree(bDA);
  
  /* Now loop over atom pairs */
  ndist = 0;
  odist = 0;
  for(i=0; (i<natoms); i++) {
    tpi = atoms->atom[i].type;
    for(j=i+1; (j<natoms); j++) {
      if (((weight[i] != 0.0) || (weight[j] != 0.0))
	  /*&& (notexcl(&atoms->excl,i,j))*/) {
	if (!dist_set(d,natoms,i,j)) {
	  /* We don't care about virtual particles */
	  ati=(*atoms->atomname[i])[0];
	  atj=(*atoms->atomname[j])[0];
	  if ( !(ati == 'V' || atj == 'V') ) {
	    /* *** */
	    tpj  = atoms->atom[j].type;
	    len  = vdwlen(atoms,i,j);
	    lb   = (1-nb_margin)*len; 
	    if (len > 0) {
	      set_dist(d,natoms,i,j,lb,maxdist,0.0);
	      /* set_dist(d,natoms,i,j,lb,0.0,lb); */
	      ndist++;
	    }
	  }
	}
	else
	  odist++;
      }
    }
  }
  fprintf(stderr,"There are %d new non-bonded distances + %d old ones\n",
	  ndist,odist);
}

static void release_domains(FILE *fp,char *ndx,t_dist dist[],t_atoms *atoms,
			    real maxdist)
{
  t_block *doms;
  char    **grpname=NULL;
  int     i,j,k,ii,jj,ai,aj,natoms;
  real    lb;
  
  natoms  = atoms->nr;
  doms    = init_index(ndx,&grpname);
  if (doms->nr > 1) {
    fprintf(fp,"Found %d domains, named:\n",doms->nr);
    for(i=0; (i<doms->nr); i++)
      fprintf(fp,"%3d  %s\n",i,grpname[i]);
    fflush(fp);
    for(i=0; (i<doms->nr); i++) {
      for(ii=doms->index[i]; (ii<doms->index[i+1]); ii++) {
	ai = doms->a[ii];
	for(j=i+1; (j<doms->nr); j++) {
	  for(jj=doms->index[j]; (jj<doms->index[j+1]); jj++) {
	    aj    = doms->a[jj];
	    if (dist_set(dist,natoms,ai,aj)) {
	      lb    = vdwlen(atoms,ai,aj);
	      set_dist(dist,natoms,ai,aj,lb,maxdist,0.0);
	    }
	  }
	}
      }
    }
  }
  done_block(doms);
  sfree(doms);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "cdist read a [BB]tpx[bb] file and dumps an input file for disco.",
    "Bond lengths etc. are read from the topology. Pairs of atoms that can",
    "form hydrogen bonds are given a lowest possible distance of",
    "[BB]hblen[bb] (can be specified by the user). Other nonbonded pairs",
    "take their minimum distance from the Lennard Jones parameters",
    "(at the combined sigma).[PAR]",
    "The program uses proper dihedrals to give a distance too, as minimum",
    "respectively maximum the [IT]cis[it] and [IT]trans[it] configurations",
    "are taken. It is therefore beneficial to use the [BB]-alldih[bb] option",
    "of [TT]pdb2gmx[tt] to generate a topology with all dihedrals in there.",
    "If the optional pdb file is given, weights are read from the occupancy",
    "field, so that",
    "not all atoms are part of the disco run, only those of which one of the",
    "weights is non-zero.[PAR]",
    "If the option -engh is on (default) bond lengths and angles etc. are",
    "read from another database, which is basically the Engh-Huber data",
    "but refined to be completely self consistent. The database name is",
    "refi_aa.dat and it resides in the $GMXLIB directory, or in the current",
    "directory.[PAR]",
    "The program can read a file with distances from NMR distance restraints",
    "(-d option). Note that these distance are treated slightly different",
    "in the disco program, and therefore these distance should be NMR",
    "derived distance restraints only.[PAR]",
    "Furthermore, the program can read an index file with hydrogen bond",
    "information as generated by [TT]g_hbond[tt]. This is then used to set",
    "tighter restraints on the hydrogen bonded atoms than on the other",
    "non bonded atom pairs, in order to maintain secondary structure.",
    "This option is useful only in combination with the [TT]-measure[tt]",
    "option, when a sensible structure is known.[PAR]",
    "The option [TT]-dom[tt] can be used to release distances bounds between",
    "different domains to the lower bounds given by Van der Waals contacts.",
    "This way, different domains can move independently, but without",
    "overlapping. The index file should contain domains that do not overlap",
    "with each other."
  };
  t_filenm fnm[] = {
    { efTPS, "-s", NULL,    ffREAD  },
    { efLOG, "-g", "cdist", ffWRITE },
    { efPDB, "-q", NULL,    ffOPTRD },
    { efDAT, "-d", NULL,    ffOPTRD },
    { efDAT, "-o", "cdist", ffWRITE },
    { efNDX, "-n", "hbond", ffOPTRD },
    { efNDX, "-dom","domain",ffOPTRD }
  };
#define NFILE asize(fnm)

  FILE        *fp;
  t_topology  *top;
  t_dist      *dist;
  real        *weight;
  rvec        *x;
  matrix      box;
  char        *topfn,title[256];
  int         i,nhb;
  atom_id     *hb;
  char        *grpname;
  
  /* Tolerance used in smoothing functions (real precision)*/

  /* Hacked 990609 by Adam */

  /* Command line options */
  static real tol=1e-6;
  static real bond_margin  = 0.01;
  static real angle_margin = 0.01;
  /* static real pep_margin   = 0.005; */
  static real pep_margin   = 0.01;
  /* static real ring_margin  = 0.002; */
  static real ring_margin  = 0.01;
  /* static real arg_margin   = 0.002; */
  static real arg_margin   = 0.01;
  /* Use end_margin for asn and gln */
  /* static real end_margin   = 0.004; */
  static real end_margin   = 0.01;
  static real val_margin   = 0.01;
  static real leu_margin   = 0.01;
  static real ile_margin   = 0.03;
  static real dih_margin   = 0.01;
  static real idih_margin  = 0.01;
  static real nb_margin    = 0.05;
  static real hb_margin    = 0.02;
  static real hblen        = 2.3;
  static real rcut         = 0.0;
  static real maxdist      = 0.0;
  static bool bNB=TRUE,bBON=TRUE,bAddR=FALSE;
  static bool bVir=FALSE,bEnghHuber=TRUE;
  static char *smth_str[]  = { NULL, "none", "tri", "tetra", NULL };
  t_pargs pa[] = {
    { "-engh",FALSE,etBOOL,  {&bEnghHuber},
      "Use the Engh&Huber parameters for bond-lengths etc." },
    { "-tol", FALSE, etREAL, {&tol},
      "HIDDENTolerance for smoothing" },
    { "-bm", FALSE, etREAL,  {&bond_margin}, 
      "Relative margin for bond lengths" },
    { "-am", FALSE, etREAL,  {&angle_margin}, 
      "Relative margin for bond angle lengths" },
    { "-pm", FALSE, etREAL,  {&pep_margin},
      "Relative margin for peptidebond dihedrals" },
    { "-rr", FALSE, etREAL,  {&ring_margin},
      "Relative margin to keep rings flat (trp,tyr,phe,hisb)" },
    { "-ar", FALSE, etREAL,  {&arg_margin},
      "Relative margin for arginine" },
    { "-er", FALSE, etREAL,  {&end_margin},
      "Relative margin for asn and gln" },
    { "-vm", FALSE, etREAL,  {&val_margin},
      "Relative margin for valine (0 disables)" },
    { "-lm", FALSE, etREAL,  {&leu_margin},
      "Relative margin for leucine (0 disables)" },
    { "-il", FALSE, etREAL,  {&ile_margin},
      "Relative margin for isoleucine (0 disables)" },
    { "-dm", FALSE, etREAL,  {&dih_margin}, 
      "!inactive! Relative margin for dihedral lengths" },
    { "-im", FALSE, etREAL,  {&idih_margin}, 
      "Relative margin for improper dihedral lengths" },
    { "-nm", FALSE, etREAL,  {&nb_margin}, 
      "Relative margin for nonbonded lower bounds" },
    { "-hm", FALSE, etREAL,  {&hb_margin},
      "Relative margin for hydrogen bonded atoms, which must be specified in an index file, as generated by g_hbond" },
    { "-hb", FALSE, etREAL,  {&hblen},
      "Shortest possible distance for a hydrogen bond (in Angstrom!)" },
    { "-bon", FALSE, etBOOL, {&bBON},
      "Make bonded distance constraints" },
    { "-nb", FALSE, etBOOL,  {&bNB},
      "Make nonbonded distance constraints (lower bound only) " },
    { "-measure", FALSE, etREAL, {&rcut},
      "Add (nonbonded) distances by examining all atoms within the distance given (in Angstrom), and using the margin given by the -nm option." },
    { "-maxdist", FALSE, etREAL, {&maxdist},
      "Maximum distance between any pair of atoms" },
    { "-add",FALSE, etBOOL,  {&bAddR},
      "Write restraints in format of additional restraints for disco" },
    { "-vir",FALSE, etBOOL,  {&bVir},
      "Use virtual particles"},
    { "-sm", FALSE, etENUM,  {smth_str},
      "Smoothing: none, tri (Using triangle inequality), or tetra (Partial tetrangle inequaliy)" }
  };
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);

#ifndef DOUBLE
  fprintf(stderr,"WARNING: running %s in single precision may produce bad" 
	  " distances\n",argv[0]);
#endif
  
  fp=ffopen(ftp2fn(efLOG,NFILE,fnm),"w");
  
  if ((!bBON && !bNB) && (rcut==0.0))
    fprintf(stderr,"That was interesting... (nothing done)\n");
  else {
    snew(top,1);
    topfn = ftp2fn(efTPS,NFILE,fnm);
    if (!read_tps_conf(topfn,title,top,&x,NULL,box,FALSE))
      gmx_fatal(FARGS,"No topology in %s",topfn);
    fprintf(stderr,"Successfully read %s (%s)\n",topfn,title);
    
    if (!opt2parg_bSet("-maxdist",asize(pa),pa))
      maxdist = (top->atoms.nres+3)*3.5; /* Ca-Ca distance + some buffer */
	  
    if (opt2bSet("-q",NFILE,fnm)) {
      weight = read_weights(opt2fn("-q",NFILE,fnm),top->atoms.nr);
    }
    else {
      snew(weight,top->atoms.nr);
      for(i=0; (i<top->atoms.nr); i++)
	weight[i] = 1.0;
    }
    if (opt2bSet("-d",NFILE,fnm)) {
      dist = read_dist(fp,opt2fn("-d",NFILE,fnm),top->atoms.nr,weight);
    }
    else
      dist = new_dist(top->atoms.nr);

    if (bEnghHuber)
      read_O_dist();
      
    if (bBON) {
      simple_bonds_and_angles(fp,dist,&top->idef,&top->atoms,weight,
			      bond_margin,angle_margin);
      peptide_bonds(fp,dist,&top->idef,&top->atoms,weight,pep_margin,
		    top->idef.il,top->idef.iparams,bVir);
      arg(fp,dist,&top->idef,&top->atoms,weight,arg_margin,
	  top->idef.il,top->idef.iparams,bVir);
      asn(fp,dist,&top->idef,&top->atoms,weight,end_margin,
	  top->idef.il,top->idef.iparams,bVir);
      gln(fp,dist,&top->idef,&top->atoms,weight,end_margin,
	  top->idef.il,top->idef.iparams,bVir);
      phe(fp,dist,&top->idef,&top->atoms,weight,ring_margin,
	  top->idef.il,top->idef.iparams,bVir);
      tyr(fp,dist,&top->idef,&top->atoms,weight,ring_margin,
	  top->idef.il,top->idef.iparams,bVir);
      trp(fp,dist,&top->idef,&top->atoms,weight,ring_margin,
	  top->idef.il,top->idef.iparams,bVir);
      hisb(fp,dist,&top->idef,&top->atoms,weight,ring_margin,
	   top->idef.il,top->idef.iparams,bVir);
      if ( val_margin != 0 ) {
	val(fp,dist,&top->idef,&top->atoms,weight,val_margin,
	    top->idef.il,top->idef.iparams);
      }
      if ( leu_margin != 0 ) {
	leu(fp,dist,&top->idef,&top->atoms,weight,leu_margin,
	    top->idef.il,top->idef.iparams);
      }
      if ( ile_margin != 0 ) {
	ile(fp,dist,&top->idef,&top->atoms,weight,ile_margin,
	    top->idef.il,top->idef.iparams);
      }
      fflush(fp);
      fprintf(stderr,"Done residues...\n");
      dump_bonds(&top->atoms,dist,
		 top->idef.il,top->idef.functype,top->idef.iparams,
		 bond_margin,angle_margin,dih_margin,
		 idih_margin,weight,bAddR); 
    }
    fprintf(stderr,"Done bondeds\n");
    
    if (rcut > 0) {
      nhb = 0;
      hb  = NULL;
      grpname = NULL;
      if (ftp2bSet(efNDX,NFILE,fnm)) {
	rd_index(ftp2fn(efNDX,NFILE,fnm),1,&nhb,&hb,&grpname);
	fprintf(stderr,"There are %d hydrogen bonds\n",nhb/3);
      }
      measure_dist(fp,dist,&top->atoms,x,rcut,nb_margin,hb_margin,nhb,hb);
    }
    if (opt2bSet("-dom",NFILE,fnm)) {
      release_domains(fp,opt2fn("-dom",NFILE,fnm),dist,&top->atoms,maxdist);
    }
    if (bNB)
      dump_nonbonds(dist,&top->idef,&top->atoms,hblen,weight,
		    nb_margin,bAddR,maxdist);
    
    if (strcmp(smth_str[0],smth_str[1]) == 0)
      fprintf(stderr,"No smoothing\n");
    else if (strcmp(smth_str[0],smth_str[2]) == 0) {
      fprintf(stderr,"Triangle smoothing only\n");
      (void) do_triangle (dist,&top->atoms,tol);
    }
    else if (strcmp(smth_str[0],smth_str[3]) == 0) {
      fprintf(stderr,"Partial tetrangle + triangle smoothing\n");
      do_smooth(dist,&top->atoms,tol);
    }
    else
      gmx_fatal(FARGS,"Uh-oh, smth_str = %s, %s, line %d",
		  smth_str[0],__FILE__,__LINE__);
    
    dump_dist(opt2fn("-o",NFILE,fnm),dist,top->atoms.nr,bAddR);
  }
  
  ffclose(fp);
  
  thanx(stderr);

  return 0;
}
