/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_residues_c = "$Id$";

#include "cdist.h"

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


void arg (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	  real arg_margin,t_ilist ilist[],t_iparams iparams[],bool bVir)
{
  int natoms,ndist,i,j,q,logg[12],residnr,oldresidnr,n12dist,n13dist,n14dist,
    n15dist,nVdist;
  real blen,lb,ub,angle;

  natoms= atoms->nr;
  ndist   = 0;
  n12dist = 0;
  n13dist = 0;
  n14dist = 0;
  n15dist = 0;
  nVdist = 0;
  residnr = -1;
  oldresidnr=-1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),"ARG") == 0) {
      oldresidnr=residnr;
      residnr=atoms->atom[i].resnr;
      logg[11]=0;
      if ( oldresidnr == residnr ) {
	continue;
      }
      j=i;
      while ((atoms->atom[j].resnr) == residnr ) {
	if ( strcmp((*atoms->atomname[j]),"NE") == 0) {
	  logg[11]++;
	  logg[0]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"HE") == 0) {
	  logg[11]++;
	  logg[1]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CZ") == 0) {
	  logg[11]++;
	  logg[2]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"NH1") == 0) {
	  logg[11]++;
	  logg[3]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"HH11") == 0) {
	  logg[11]++;
	  logg[4]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HH12") == 0) {
	  logg[11]++;
	  logg[5]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"NH2") == 0) {
	  logg[11]++;
	  logg[6]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HH21") == 0) {
	  logg[11]++;
	  logg[7]=j; 
	}
        else if ( strcmp((*atoms->atomname[j]),"HH22") == 0) {
	  logg[11]++;
	  logg[8]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"CD") == 0) {
	  logg[11]++;
	  logg[9]=j;
	}
	else if ( (strcmp((*atoms->atomname[j]),"VF") == 0) && bVir) {
	  logg[11]++;
	  logg[10]=j;
	}
	if ( ((logg[11] == 10) && !bVir) || ((logg[11] == 11) && bVir) ) {
	  break;
	}
	j++;
      }
    if ( ((logg[11] == 10) && !bVir) || ((logg[11] == 11) && bVir) ) {

      if  ((logg[11] == 10) && !bVir) {
      fprintf(log,"logg (arg) = %d %d %d %d %d %d %d %d %d %d\n",logg[0],
	      logg[1],logg[2],logg[3],logg[4],logg[5],logg[6],logg[7],logg[8],
	      logg[9]);
      }
      else if ((logg[11] == 11) && bVir) {
      fprintf(log,"logg (arg) = %d %d %d %d %d %d %d %d %d %d %d\n",logg[0],
	      logg[1],logg[2],logg[3],logg[4],logg[5],logg[6],logg[7],logg[8],
	      logg[9],logg[10]);
      }

      /*SETDISTANCE for NE and HH12 (logg[0]&logg[5]) (transdihedral)*/
      pdih_lengths(logg[0],logg[2],logg[3],logg[5],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[5]))) {
	set_dist(d,natoms,logg[0],logg[5],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for NE and HH22 (logg[0]&logg[8]) (transdihedral)*/
      pdih_lengths(logg[0],logg[2],logg[6],logg[8],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[8]))) {
	set_dist(d,natoms,logg[0],logg[8],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HE and NH1 (logg[1]&logg[3]) (transdihedral)*/
      pdih_lengths(logg[1],logg[0],logg[2],logg[3],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[3]))) {
	set_dist(d,natoms,logg[1],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HH21 and NH1 (logg[7]&logg[3]) (transdihedral)*/
      pdih_lengths(logg[7],logg[6],logg[2],logg[3],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[3]))) {
	set_dist(d,natoms,logg[7],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HH11 and NH2 (logg[4]&logg[6]) (transdihedral)*/
      pdih_lengths(logg[4],logg[3],logg[2],logg[6],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[6]))) {
	set_dist(d,natoms,logg[4],logg[6],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CD and NH2 (logg[9]&logg[6]) (transdihedral)*/
      pdih_lengths(logg[9],logg[0],logg[2],logg[6],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[6]))) {
	set_dist(d,natoms,logg[9],logg[6],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CD and NH1 (logg[9]&logg[3]) (cisdihedral)*/
      pdih_lengths(logg[9],logg[0],logg[2],logg[3],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[3]))) {
	set_dist(d,natoms,logg[9],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HH12 and NH2 (logg[5]&logg[6]) (cisdihedral)*/
      pdih_lengths(logg[5],logg[3],logg[2],logg[6],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[6]))) {
	set_dist(d,natoms,logg[5],logg[6],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HH22 and NH1 (logg[8]&logg[3]) (cisdihedral)*/
      pdih_lengths(logg[8],logg[6],logg[2],logg[3],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[3]))) {
	set_dist(d,natoms,logg[8],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HE and NH2 (logg[1]&logg[6]) (cisdihedral)*/
      pdih_lengths(logg[1],logg[0],logg[2],logg[6],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[6]))) {
	set_dist(d,natoms,logg[1],logg[6],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for NE and HH21 (logg[0]&logg[7]) (cisdihedral)*/
      pdih_lengths(logg[0],logg[2],logg[6],logg[7],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[7]))) {
	set_dist(d,natoms,logg[0],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for NE and HH11 (logg[0]&logg[4]) (cisdihedral)*/
      pdih_lengths(logg[0],logg[2],logg[3],logg[4],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[4]))) {
	set_dist(d,natoms,logg[0],logg[4],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CD and NE (1-2) */
      blen=lookup_bondlength(logg[9],logg[0],ilist,iparams,TRUE,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[0]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[0]))) {
	set_dist(d,natoms,logg[9],logg[0],lb,ub,blen);
	ndist++;
	n12dist++;
      }       
      /*SETDISTANCE for NE and HE (1-2) */
      blen=lookup_bondlength(logg[0],logg[1],ilist,iparams,TRUE,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[1]))) {
	set_dist(d,natoms,logg[0],logg[1],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for NE and CZ (1-2) */
      blen=lookup_bondlength(logg[0],logg[2],ilist,iparams,TRUE,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[2]))) {
	set_dist(d,natoms,logg[0],logg[2],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CZ and NH1 (1-2) */
      blen=lookup_bondlength(logg[2],logg[3],ilist,iparams,TRUE,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[3]))) {
	set_dist(d,natoms,logg[2],logg[3],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CZ and NH2 (1-2) */
      blen=lookup_bondlength(logg[2],logg[6],ilist,iparams,TRUE,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[6]))) {
	set_dist(d,natoms,logg[2],logg[6],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for NH1 and HH11 (1-2) */
      blen=lookup_bondlength(logg[3],logg[4],ilist,iparams,TRUE,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[4]))) {
	set_dist(d,natoms,logg[3],logg[4],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for NH1 and HH12 (1-2) */
      blen=lookup_bondlength(logg[3],logg[5],ilist,iparams,TRUE,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[5]))) {
	set_dist(d,natoms,logg[3],logg[5],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for NH2 and HH21 (1-2) */
      blen=lookup_bondlength(logg[6],logg[7],ilist,iparams,TRUE,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[7]))) {
	set_dist(d,natoms,logg[6],logg[7],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for NH2 and HH22 (1-2) */
      blen=lookup_bondlength(logg[6],logg[8],ilist,iparams,TRUE,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[8]))) {
	set_dist(d,natoms,logg[6],logg[8],lb,ub,blen);
	ndist++;
	n12dist++;
      }

      /* new 1-3 distances added 981126 */
      /*SETDISTANCE for NH1 and NH2 (1-3) */
      angle=lookup_angle(logg[3],logg[2],logg[6],ilist,iparams,atoms);
      blen=angle_length(logg[3],logg[2],logg[6],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[6]))) {
	set_dist(d,natoms,logg[3],logg[6],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HH11 and HH12 (1-3) */
      angle=lookup_angle(logg[5],logg[3],logg[4],ilist,iparams,atoms);
      blen=angle_length(logg[5],logg[3],logg[4],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[4]))) {
	set_dist(d,natoms,logg[5],logg[4],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HH21 and HH22 (1-3) */
      angle=lookup_angle(logg[8],logg[6],logg[7],ilist,iparams,atoms);
      blen=angle_length(logg[8],logg[6],logg[7],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[7]))) {
	set_dist(d,natoms,logg[8],logg[7],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /* end of new 1-3 distances */

      /*SETDISTANCE for CD and HE (1-3) */
      angle=lookup_angle(logg[9],logg[0],logg[1],ilist,iparams,atoms);
      blen=angle_length(logg[9],logg[0],logg[1],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[1]))) {
	set_dist(d,natoms,logg[9],logg[1],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CD and CZ (1-3) */
      angle=lookup_angle(logg[9],logg[0],logg[2],ilist,iparams,atoms);
      blen=angle_length(logg[9],logg[0],logg[2],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[2]))) {
	set_dist(d,natoms,logg[9],logg[2],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HE and CZ (1-3) */
      angle=lookup_angle(logg[1],logg[0],logg[2],ilist,iparams,atoms);
      blen=angle_length(logg[1],logg[0],logg[2],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[2]))) {
	set_dist(d,natoms,logg[1],logg[2],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for NE and NH1 (1-3) */
      angle=lookup_angle(logg[0],logg[2],logg[3],ilist,iparams,atoms);
      blen=angle_length(logg[0],logg[2],logg[3],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[3]))) {
	set_dist(d,natoms,logg[0],logg[3],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for NE and NH2 (1-3) */
      angle=lookup_angle(logg[0],logg[2],logg[6],ilist,iparams,atoms);
      blen=angle_length(logg[0],logg[2],logg[6],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[6]))) {
	set_dist(d,natoms,logg[0],logg[6],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CZ and HH11 (1-3) */
      angle=lookup_angle(logg[2],logg[3],logg[4],ilist,iparams,atoms);
      blen=angle_length(logg[2],logg[3],logg[4],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[4]))) {
	set_dist(d,natoms,logg[2],logg[4],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CZ and HH12 (1-3) */
      angle=lookup_angle(logg[2],logg[3],logg[5],ilist,iparams,atoms);
      blen=angle_length(logg[2],logg[3],logg[5],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[5]))) {
	set_dist(d,natoms,logg[2],logg[5],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CZ and HH21 (1-3) */
      angle=lookup_angle(logg[2],logg[6],logg[7],ilist,iparams,atoms);
      blen=angle_length(logg[2],logg[6],logg[7],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[7]))) {
	set_dist(d,natoms,logg[2],logg[7],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CZ and HH22 (1-3) */
      angle=lookup_angle(logg[2],logg[6],logg[8],ilist,iparams,atoms);
      blen=angle_length(logg[2],logg[6],logg[8],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[8]))) {
	set_dist(d,natoms,logg[2],logg[8],lb,ub,blen);
	ndist++;
	n13dist++;
      }

      /* new 1-5 distances added 981126 */
      /*SETDISTANCE for HE and HH12 (1-5) (trans) */
      arg_15_CDHH2(logg[1],logg[0],logg[2],logg[3],logg[5],ilist,
		      iparams,&lb,&blen,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[5]))) {
	set_dist(d,natoms,logg[1],logg[5],lb,ub,blen);
	ndist++;
	n15dist++;
      }      
      /*SETDISTANCE for HE and HH11 (1-5) (trans) */
      arg_15_CDHH2(logg[1],logg[0],logg[2],logg[3],logg[4],ilist,
		      iparams,&blen,&ub,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[4]))) {
	set_dist(d,natoms,logg[1],logg[4],lb,ub,blen);
	ndist++;
	n15dist++;
      }      
      /*SETDISTANCE for HE and HH22 (1-5) (trans) */
      arg_15_CDHH1(logg[1],logg[0],logg[2],logg[6],logg[8],ilist,
		      iparams,&lb,&blen,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[8]))) {
	set_dist(d,natoms,logg[1],logg[8],lb,ub,blen);
	ndist++;
	n15dist++;
      }      
      /*SETDISTANCE for HE and HH21 (1-5) (trans) */
      arg_15_CDHH1(logg[1],logg[0],logg[2],logg[6],logg[7],ilist,
		      iparams,&blen,&ub,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[7]))) {
	set_dist(d,natoms,logg[1],logg[7],lb,ub,blen);
	ndist++;
	n15dist++;
      }      
      /*SETDISTANCE for HH11 and HH21 (1-5) (trans) */
      arg_15_CDHH2(logg[4],logg[3],logg[2],logg[6],logg[7],ilist,
		      iparams,&lb,&blen,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[7]))) {
	set_dist(d,natoms,logg[4],logg[7],lb,ub,blen);
	ndist++;
	n15dist++;
      }      
      /*SETDISTANCE for HH11 and HH22 (1-5) (trans) */
      arg_15_CDHH2(logg[4],logg[3],logg[2],logg[6],logg[8],ilist,
		      iparams,&blen,&ub,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[8]))) {
	set_dist(d,natoms,logg[4],logg[8],lb,ub,blen);
	ndist++;
	n15dist++;
      }      
      /*SETDISTANCE for HH12 and HH21 (1-5) (trans) */
      arg_15_CDHH1(logg[5],logg[3],logg[2],logg[6],logg[7],ilist,
		      iparams,&lb,&blen,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[7]))) {
	set_dist(d,natoms,logg[5],logg[7],lb,ub,blen);
	ndist++;
	n15dist++;
      }      
      /*SETDISTANCE for HH12 and HH22 (1-5) (trans) */
      arg_15_CDHH1(logg[5],logg[3],logg[2],logg[6],logg[8],ilist,
		      iparams,&blen,&ub,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[8]))) {
	set_dist(d,natoms,logg[5],logg[8],lb,ub,blen);
	ndist++;
	n15dist++;
      }  
      /* end of new 1-5 distances */

      /*SETDISTANCE for CD and HH12 (1-5) (trans) */
      arg_15_CDHH1(logg[9],logg[0],logg[2],logg[3],logg[5],ilist,
		      iparams,&lb,&blen,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[5]))) {
	set_dist(d,natoms,logg[9],logg[5],lb,ub,blen);
	ndist++;
	n15dist++;
      }      
      /*SETDISTANCE for CD and HH11 (1-5) (cis) */
      arg_15_CDHH1(logg[9],logg[0],logg[2],logg[3],logg[4],ilist,
		      iparams,&blen,&ub,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[4]))) {
	set_dist(d,natoms,logg[9],logg[4],lb,ub,blen);
	ndist++;
	n15dist++;
      }    
      /*SETDISTANCE for CD and HH21 (1-5) (cis) */
      arg_15_CDHH2(logg[9],logg[0],logg[2],logg[6],logg[7],ilist,
		      iparams,&blen,&ub,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[7]))) {
	set_dist(d,natoms,logg[9],logg[7],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for CD and HH22 (1-5) (trans) */
      arg_15_CDHH2(logg[9],logg[0],logg[2],logg[6],logg[8],ilist,
		      iparams,&lb,&blen,atoms);
      lb=(1.0-arg_margin)*blen;
      ub=(1.0+arg_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[8]))) {
	set_dist(d,natoms,logg[9],logg[8],lb,ub,blen);
	ndist++;
	n15dist++;
      } 

      if (bVir) {
	/* VIRTUAL DISTANCES */
	nVdist += set_virtual (logg,10,arg_margin,d,natoms);
	ndist += nVdist;
      }
    }
    logg[11]=0;
    }
  }

  fprintf(log,"There are %d new arginine distances\n",ndist);
  if (ndist > 0 ) {
    fprintf(log,"(%d 1-2, %d 1-3, %d 1-4, %d 1-5, %d virtual)\n",
	    n12dist,n13dist,n14dist,n15dist,nVdist);
  }
}

/**********************************************************
 *
 *     A S P A R A G I N E
 *
 **********************************************************/
void asn (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	  real end_margin,t_ilist ilist[],t_iparams iparams[],bool bVir)
{
  int natoms,ndist,i,j,q,logg[8],residnr,oldresidnr,n12dist,n13dist,n14dist,nVdist;
  real blen,lb,ub,angle;

  natoms= atoms->nr;
  ndist   = 0;
  n12dist = 0;
  n13dist = 0;
  n14dist = 0;
  nVdist = 0;
  residnr = -1;
  oldresidnr=-1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),"ASN") == 0) {
      oldresidnr=residnr;
      residnr=atoms->atom[i].resnr;
      logg[7]=0;
      if ( oldresidnr == residnr ) {
	continue;
      }
      j=i;
      while ((atoms->atom[j].resnr) == residnr ) {
	if ( strcmp((*atoms->atomname[j]),"CB") == 0) {
	  logg[7]++;
	  logg[0]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CG") == 0) {
	  logg[7]++;
	  logg[1]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"OD1") == 0) {
	  logg[7]++;
	  logg[2]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"ND2") == 0) {
	  logg[7]++;
	  logg[3]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"HD21") == 0) {
	  logg[7]++;
	  logg[4]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HD22") == 0) {
	  logg[7]++;
	  logg[5]=j;
	}
        else if ( (strcmp((*atoms->atomname[j]),"VF") == 0) && bVir) {
	  logg[7]++;
	  logg[6]=j;
	}
	if ( ((logg[7] == 6) && !bVir) || ((logg[7] == 7) && bVir) ) {
	  break;
	}
	j++;
      }
    if ( ((logg[7] == 6) && !bVir) || ((logg[7] == 7) && bVir) ) {

      if ((logg[7] == 6) && !bVir) {
      fprintf(log,"logg (asn) = %d %d %d %d %d %d\n",logg[0],
	      logg[1],logg[2],logg[3],logg[4],logg[5]);
      }
      else if ((logg[7] == 7) && bVir) {
      fprintf(log,"logg (asn) = %d %d %d %d %d %d %d\n",logg[0],
	      logg[1],logg[2],logg[3],logg[4],logg[5],logg[6]);
      }

      /*SETDISTANCE for CB and CG (1-2) */
      blen=lookup_bondlength(logg[0],logg[1],ilist,iparams,TRUE,atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[1]))) {
	set_dist(d,natoms,logg[0],logg[1],lb,ub,blen);
	ndist++;
	n12dist++;
      }       
      /*SETDISTANCE for CG and OD1 (1-2) */
      blen=lookup_bondlength(logg[1],logg[2],ilist,iparams,TRUE,atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[2]))) {
	set_dist(d,natoms,logg[1],logg[2],lb,ub,blen);
	ndist++;
	n12dist++;
      }     
      /*SETDISTANCE for CG and ND2 (1-2) */
      blen=lookup_bondlength(logg[1],logg[3],ilist,iparams,TRUE,atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[3]))) {
	set_dist(d,natoms,logg[1],logg[3],lb,ub,blen);
	ndist++;
	n12dist++;
      }     
      /*SETDISTANCE for ND2 and HD21 (1-2) */
      blen=lookup_bondlength(logg[3],logg[4],ilist,iparams,TRUE,atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[4]))) {
	set_dist(d,natoms,logg[3],logg[4],lb,ub,blen);
	ndist++;
	n12dist++;
      }     
      /*SETDISTANCE for ND2 and HD22 (1-2) */
      blen=lookup_bondlength(logg[3],logg[5],ilist,iparams,TRUE,atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[5]))) {
	set_dist(d,natoms,logg[3],logg[5],lb,ub,blen);
	ndist++;
	n12dist++;
      }     

      /* new 1-3 distances added 981126 */
      /*SETDISTANCE for HD21 and HD22 (1-3) */
      angle=lookup_angle(logg[4],logg[3],logg[5],ilist,iparams,atoms);
      blen=angle_length(logg[4],logg[3],logg[5],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[5]))) {
	set_dist(d,natoms,logg[4],logg[5],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /* end of new 1-3 distances */

      /*SETDISTANCE for CB and OD1 (1-3) */
      angle=lookup_angle(logg[0],logg[1],logg[2],ilist,iparams,atoms);
      blen=angle_length(logg[0],logg[1],logg[2],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[2]))) {
	set_dist(d,natoms,logg[0],logg[2],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CB and ND2 (1-3) */
      angle=lookup_angle(logg[0],logg[1],logg[3],ilist,iparams,atoms);
      blen=angle_length(logg[0],logg[1],logg[3],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[3]))) {
	set_dist(d,natoms,logg[0],logg[3],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for OD1 and ND2 (1-3) */
      angle=lookup_angle(logg[2],logg[1],logg[3],ilist,iparams,atoms);
      blen=angle_length(logg[2],logg[1],logg[3],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[3]))) {
	set_dist(d,natoms,logg[2],logg[3],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CG and HD21 (1-3) */
      angle=lookup_angle(logg[1],logg[3],logg[4],ilist,iparams,atoms);
      blen=angle_length(logg[1],logg[3],logg[4],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[4]))) {
	set_dist(d,natoms,logg[1],logg[4],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CG and HD22 (1-3) */
      angle=lookup_angle(logg[1],logg[3],logg[5],ilist,iparams,atoms);
      blen=angle_length(logg[1],logg[3],logg[5],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[5]))) {
	set_dist(d,natoms,logg[1],logg[5],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CB and HD22 (logg[0]&logg[5]) (transdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[3],logg[5],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[5]))) {
	set_dist(d,natoms,logg[0],logg[5],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for OD1 and HD21 (logg[2]&logg[4]) (transdihedral)*/
      pdih_lengths(logg[2],logg[1],logg[3],logg[4],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[4]))) {
	set_dist(d,natoms,logg[2],logg[4],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CB and HD21 (logg[0]&logg[4]) (cisdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[3],logg[4],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[4]))) {
	set_dist(d,natoms,logg[0],logg[4],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for OD1 and HD22 (logg[2]&logg[5]) (cisdihedral)*/
      pdih_lengths(logg[2],logg[1],logg[3],logg[5],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[5]))) {
	set_dist(d,natoms,logg[2],logg[5],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      if (bVir) {
	/* VIRTUAL DISTANCES */
	nVdist += set_virtual (logg,6,end_margin,d,natoms);
	ndist += nVdist;
      }
    }
    logg[7]=0;
    }
  }
  fprintf(log,"There are %d new asparagine distances\n",ndist);
  if (ndist > 0) {
    fprintf(log,"(%d 1-2, %d 1-3, %d 1-4, %d virtual)\n",
	    n12dist,n13dist,n14dist,nVdist);
  }
}

/**********************************************************
 *
 *     G L U T A M I N E
 *
 **********************************************************/
void gln (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	  real end_margin,t_ilist ilist[],t_iparams iparams[],bool bVir)
{
  int natoms,ndist,i,j,q,logg[8],residnr,oldresidnr,n12dist,n13dist,n14dist,nVdist;
  real blen,lb,ub,angle;

  natoms= atoms->nr;
  ndist   = 0;
  n12dist = 0;
  n13dist = 0;
  n14dist = 0;
  nVdist  = 0;
  residnr = -1;
  oldresidnr=-1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),"GLN") == 0) {
      oldresidnr=residnr;
      residnr=atoms->atom[i].resnr;
      logg[7]=0;
      if ( oldresidnr == residnr ) {
	continue;
      }
      j=i;
      while ((atoms->atom[j].resnr) == residnr ) {
	if ( strcmp((*atoms->atomname[j]),"CG") == 0) {
	  logg[7]++;
	  logg[0]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CD") == 0) {
	  logg[7]++;
	  logg[1]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"OE1") == 0) {
	  logg[7]++;
	  logg[2]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"NE2") == 0) {
	  logg[7]++;
	  logg[3]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"HE21") == 0) {
	  logg[7]++;
	  logg[4]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HE22") == 0) {
	  logg[7]++;
	  logg[5]=j;
	}
        else if ( (strcmp((*atoms->atomname[j]),"VF") == 0) && bVir) {
	  logg[7]++;
	  logg[6]=j;
	}
	if ( ((logg[7] == 6) && !bVir) || ((logg[7] == 7) && bVir) ) {
	  break;
	}
	j++;
      }
    if ( ((logg[7] == 6) && !bVir) || ((logg[7] == 7) && bVir) ) {

      if ((logg[7] == 6) && !bVir) {
	fprintf(log,"logg (arg) = %d %d %d %d %d %d\n",logg[0],
		logg[1],logg[2],logg[3],logg[4],logg[5]);
      }
      else if ((logg[7] == 7) && bVir) {
	fprintf(log,"logg (arg) = %d %d %d %d %d %d %d\n",logg[0],
		logg[1],logg[2],logg[3],logg[4],logg[5],logg[6]);
      }
      /*SETDISTANCE for CG and CD (1-2) */
      blen=lookup_bondlength(logg[0],logg[1],ilist,iparams,TRUE,atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[1]))) {
	set_dist(d,natoms,logg[0],logg[1],lb,ub,blen);
	ndist++;
	n12dist++;
      }       
      /*SETDISTANCE for CD and OE1 (1-2) */
      blen=lookup_bondlength(logg[1],logg[2],ilist,iparams,TRUE,atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[2]))) {
	set_dist(d,natoms,logg[1],logg[2],lb,ub,blen);
	ndist++;
	n12dist++;
      }     
      /*SETDISTANCE for CD and NE2 (1-2) */
      blen=lookup_bondlength(logg[1],logg[3],ilist,iparams,TRUE,atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[3]))) {
	set_dist(d,natoms,logg[1],logg[3],lb,ub,blen);
	ndist++;
	n12dist++;
      }     
      /*SETDISTANCE for NE2 and HE21 (1-2) */
      blen=lookup_bondlength(logg[3],logg[4],ilist,iparams,TRUE,atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[4]))) {
	set_dist(d,natoms,logg[3],logg[4],lb,ub,blen);
	ndist++;
	n12dist++;
      }     
      /*SETDISTANCE for NE2 and HE22 (1-2) */
      blen=lookup_bondlength(logg[3],logg[5],ilist,iparams,TRUE,atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[5]))) {
	set_dist(d,natoms,logg[3],logg[5],lb,ub,blen);
	ndist++;
	n12dist++;
      }     

      /* new 1-3 distances added 981126 */
      /*SETDISTANCE for HE21 and HE22 (1-3) */
      angle=lookup_angle(logg[4],logg[3],logg[5],ilist,iparams,atoms);
      blen=angle_length(logg[4],logg[3],logg[5],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[5]))) {
	set_dist(d,natoms,logg[4],logg[5],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /* end of new 1-3 distances */

      /*SETDISTANCE for CG and OE1 (1-3) */
      angle=lookup_angle(logg[0],logg[1],logg[2],ilist,iparams,atoms);
      blen=angle_length(logg[0],logg[1],logg[2],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[2]))) {
	set_dist(d,natoms,logg[0],logg[2],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CG and NE2 (1-3) */
      angle=lookup_angle(logg[0],logg[1],logg[3],ilist,iparams,atoms);
      blen=angle_length(logg[0],logg[1],logg[3],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[3]))) {
	set_dist(d,natoms,logg[0],logg[3],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for OE1 and NE2 (1-3) */
      angle=lookup_angle(logg[2],logg[1],logg[3],ilist,iparams,atoms);
      blen=angle_length(logg[2],logg[1],logg[3],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[3]))) {
	set_dist(d,natoms,logg[2],logg[3],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CD and HE21 (1-3) */
      angle=lookup_angle(logg[1],logg[3],logg[4],ilist,iparams,atoms);
      blen=angle_length(logg[1],logg[3],logg[4],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[4]))) {
	set_dist(d,natoms,logg[1],logg[4],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CD and HE22 (1-3) */
      angle=lookup_angle(logg[1],logg[3],logg[5],ilist,iparams,atoms);
      blen=angle_length(logg[1],logg[3],logg[5],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[5]))) {
	set_dist(d,natoms,logg[1],logg[5],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CG and HE22 (logg[0]&logg[5]) (transdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[3],logg[5],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[5]))) {
	set_dist(d,natoms,logg[0],logg[5],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for OD1 and HE21 (logg[2]&logg[4]) (transdihedral)*/
      pdih_lengths(logg[2],logg[1],logg[3],logg[4],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[4]))) {
	set_dist(d,natoms,logg[2],logg[4],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CG and HE21 (logg[0]&logg[4]) (cisdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[3],logg[4],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[4]))) {
	set_dist(d,natoms,logg[0],logg[4],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for OD1 and HE22 (logg[2]&logg[5]) (cisdihedral)*/
      pdih_lengths(logg[2],logg[1],logg[3],logg[5],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-end_margin)*blen;
      ub=(1.0+end_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[5]))) {
	set_dist(d,natoms,logg[2],logg[5],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      if (bVir) {
	/* VIRTUAL DISTANCES */
	nVdist += set_virtual (logg,6,end_margin,d,natoms);
	ndist += nVdist;
      }
    }
    logg[7]=0;
    }
  }
  fprintf(log,"There are %d new glutamine distances\n",ndist);
  if (ndist > 0 ) {
    fprintf(log,"(%d 1-2, %d 1-3, %d 1-4)\n",n12dist,n13dist,n14dist);
  }
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
  /*fprintf(stderr,"leaving routine\n");*/
}


void hisb (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	   real ring_margin,t_ilist ilist[],t_iparams iparams[],bool bVir)
{
  int natoms,ndist,i,j,q,logg[11],residnr,oldresidnr,n12dist,n13dist,n14dist,
    n15dist,nVdist;
  real blen,lb,ub,angle;

  natoms= atoms->nr;
  ndist   = 0;
  n12dist = 0;
  n13dist = 0;
  n14dist = 0;
  n15dist = 0;
  nVdist  = 0;
  residnr = -1;
  oldresidnr=-1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),"HISB") == 0) {
      oldresidnr=residnr;
      residnr=atoms->atom[i].resnr;
      logg[10]=0;
      if ( oldresidnr == residnr ) {
	continue;
      }
      j=i;
      while ((atoms->atom[j].resnr) == residnr ) {
	if ( strcmp((*atoms->atomname[j]),"CB") == 0) {
	  logg[10]++;
	  logg[0]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CG") == 0) {
	  logg[10]++;
	  logg[1]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"ND1") == 0) {
	  logg[10]++;
	  logg[2]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CE1") == 0) {
	  logg[10]++;
	  logg[3]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"NE2") == 0) {
	  logg[10]++;
	  logg[4]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"CD2") == 0) {
	  logg[10]++;
	  logg[5]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HE1") == 0) {
	  logg[10]++;
	  logg[6]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HE2") == 0) {
	  logg[10]++;
	  logg[7]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HD2") == 0) {
	  logg[10]++;
	  logg[8]=j;
	}
        else if ( (strcmp((*atoms->atomname[j]),"VF") == 0) && bVir) {
	  logg[10]++;
	  logg[9]=j;
	}
	if ( ((logg[10] == 9) && !bVir) || ((logg[10] == 10) && bVir) ) {
	  break;
	}
	j++;
      }
    if ( ((logg[10] == 9) && !bVir) || ((logg[10] == 10) && bVir) ) {

      if ((logg[10] == 9) && !bVir) {
	fprintf(log,"logg (hisb) = %d %d %d %d %d %d %d %d %d\n",
		logg[0],logg[1],logg[2],logg[3],logg[4],logg[5],logg[6],
		logg[7],logg[8]);
      }
      else if ((logg[10] == 10) && bVir) {
	fprintf(log,"logg (hisb) = %d %d %d %d %d %d %d %d %d %d\n",
		logg[0],logg[1],logg[2],logg[3],logg[4],logg[5],logg[6],
		logg[7],logg[8],logg[9]);
      }

      /*SETDISTANCE for CB and CG (1-2) */
      blen=lookup_bondlength(logg[0],logg[1],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[1]))) {
	set_dist(d,natoms,logg[0],logg[1],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CG and ND1 (1-2) */
      blen=lookup_bondlength(logg[1],logg[2],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[2]))) {
	set_dist(d,natoms,logg[1],logg[2],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for ND1 and CE1 (1-2) */
      blen=lookup_bondlength(logg[2],logg[3],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[3]))) {
	set_dist(d,natoms,logg[2],logg[3],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CE1 and NE2 (1-2) */
      blen=lookup_bondlength(logg[3],logg[4],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[4]))) {
	set_dist(d,natoms,logg[3],logg[4],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for NE2 and CD2 (1-2) */
      blen=lookup_bondlength(logg[4],logg[5],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[5]))) {
	set_dist(d,natoms,logg[4],logg[5],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CD2 and CG (1-2) */
      blen=lookup_bondlength(logg[5],logg[1],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[1]))) {
	set_dist(d,natoms,logg[5],logg[1],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CE1 and HE1 (1-2) */
      blen=lookup_bondlength(logg[3],logg[6],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[6]))) {
	set_dist(d,natoms,logg[3],logg[6],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for NE2 and HE2 (1-2) */
      blen=lookup_bondlength(logg[4],logg[7],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[7]))) {
	set_dist(d,natoms,logg[4],logg[7],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CD2 and HD2 (1-2) */
      blen=lookup_bondlength(logg[5],logg[8],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[8]))) {
	set_dist(d,natoms,logg[5],logg[8],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CB and ND1 (1-3) */
      angle=lookup_angle(logg[0],logg[1],logg[2],ilist,iparams,atoms);
      blen=angle_length(logg[0],logg[1],logg[2],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[2]))) {
	set_dist(d,natoms,logg[0],logg[2],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CB and CD2 (1-3) */
      angle=lookup_angle(logg[0],logg[1],logg[5],ilist,iparams,atoms);
      blen=angle_length(logg[0],logg[1],logg[5],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[5]))) {
	set_dist(d,natoms,logg[0],logg[5],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CD2 and ND1 (1-3) */
      angle=lookup_angle(logg[5],logg[1],logg[2],ilist,iparams,atoms);
      blen=angle_length(logg[5],logg[1],logg[2],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[2]))) {
	set_dist(d,natoms,logg[5],logg[2],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CG  and CE1 (1-3) */
      angle=lookup_angle(logg[1],logg[2],logg[3],ilist,iparams,atoms);
      blen=angle_length(logg[1],logg[2],logg[3],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[3]))) {
	set_dist(d,natoms,logg[1],logg[3],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for ND1 and NE2 (1-3) */
      angle=lookup_angle(logg[2],logg[3],logg[4],ilist,iparams,atoms);
      blen=angle_length(logg[2],logg[3],logg[4],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[4]))) {
	set_dist(d,natoms,logg[2],logg[4],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CE1 and CD2 (1-3) */
      angle=lookup_angle(logg[3],logg[4],logg[5],ilist,iparams,atoms);
      blen=angle_length(logg[3],logg[4],logg[5],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[5]))) {
	set_dist(d,natoms,logg[3],logg[5],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for NE2 and CG  (1-3) */
      angle=lookup_angle(logg[4],logg[5],logg[1],ilist,iparams,atoms);
      blen=angle_length(logg[4],logg[5],logg[1],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[1]))) {
	set_dist(d,natoms,logg[4],logg[1],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HD2 and CG (1-3) */
      angle=lookup_angle(logg[8],logg[5],logg[1],ilist,iparams,atoms);
      blen=angle_length(logg[8],logg[5],logg[1],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[1]))) {
	set_dist(d,natoms,logg[8],logg[1],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HD2 and NE2 (1-3) */
      angle=lookup_angle(logg[8],logg[5],logg[4],ilist,iparams,atoms);
      blen=angle_length(logg[8],logg[5],logg[4],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[4]))) {
	set_dist(d,natoms,logg[8],logg[4],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HE2 and CD2 (1-3) */
      angle=lookup_angle(logg[7],logg[4],logg[5],ilist,iparams,atoms);
      blen=angle_length(logg[7],logg[4],logg[5],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[5]))) {
	set_dist(d,natoms,logg[7],logg[5],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HE2 and CE1 (1-3) */
      angle=lookup_angle(logg[7],logg[4],logg[3],ilist,iparams,atoms);
      blen=angle_length(logg[7],logg[4],logg[3],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[3]))) {
	set_dist(d,natoms,logg[7],logg[3],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HE1 and NE2 (1-3) */
      angle=lookup_angle(logg[6],logg[3],logg[4],ilist,iparams,atoms);
      blen=angle_length(logg[6],logg[3],logg[4],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[4]))) {
	set_dist(d,natoms,logg[6],logg[4],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HE1 and ND1 (1-3) */
      angle=lookup_angle(logg[6],logg[3],logg[2],ilist,iparams,atoms);
      blen=angle_length(logg[6],logg[3],logg[2],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[2]))) {
	set_dist(d,natoms,logg[6],logg[2],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CB and CE1 (logg[0]&logg[3]) (transdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[2],logg[3],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[3]))) {
	set_dist(d,natoms,logg[0],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CB and NE2 (logg[0]&logg[4]) (transdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[5],logg[4],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[4]))) {
	set_dist(d,natoms,logg[0],logg[4],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HE1 and CG  (logg[6]&logg[1]) (transdihedral)*/
      pdih_lengths(logg[6],logg[3],logg[2],logg[1],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[1]))) {
	set_dist(d,natoms,logg[6],logg[1],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HE1 and CD2 (logg[6]&logg[5]) (transdihedral)*/
      pdih_lengths(logg[6],logg[3],logg[4],logg[5],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[5]))) {
	set_dist(d,natoms,logg[6],logg[5],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HE2 and ND1 (logg[7]&logg[2]) (transdihedral)*/
      pdih_lengths(logg[7],logg[4],logg[3],logg[2],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[2]))) {
	set_dist(d,natoms,logg[7],logg[2],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HE2 and CG  (logg[7]&logg[1]) (transdihedral)*/
      pdih_lengths(logg[7],logg[4],logg[5],logg[1],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[1]))) {
	set_dist(d,natoms,logg[7],logg[1],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HE1 and HE2 (logg[6]&logg[7]) (cisdihedral)*/
      pdih_lengths(logg[6],logg[3],logg[4],logg[7],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[7]))) {
	set_dist(d,natoms,logg[6],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HE2 and HD2 (logg[7]&logg[5]) (cisdihedral)*/
      pdih_lengths(logg[7],logg[4],logg[5],logg[8],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[8]))) {
	set_dist(d,natoms,logg[7],logg[8],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CB and HD2 (logg[0]&logg[8]) (cisdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[5],logg[8],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[8]))) {
	set_dist(d,natoms,logg[0],logg[8],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for ND1 and HD2 (logg[2]&logg[8]) (transdihedral)*/
      pdih_lengths(logg[2],logg[1],logg[5],logg[8],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[8]))) {
	set_dist(d,natoms,logg[2],logg[8],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CE1 and HD2 (logg[3]&logg[8]) (transdihedral)*/
      pdih_lengths(logg[3],logg[4],logg[5],logg[8],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[8]))) {
	set_dist(d,natoms,logg[3],logg[8],lb,ub,blen);
	ndist++;
	n14dist++;
      }

      /* new 1-5 distance added 981126 */

      /*SETDISTANCE for HE1 and HD2 (1-5) */
      hisb_15_type2(logg[6],logg[3],logg[4],logg[5],logg[8],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[8]))) {
	set_dist(d,natoms,logg[6],logg[8],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /* end new distances */

      /*SETDISTANCE for CB and HE1 (1-5) */
      hisb_15_type2(logg[0],logg[1],logg[2],logg[3],logg[6],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[6]))) {
	set_dist(d,natoms,logg[0],logg[6],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for CB and HE2 (1-5) */
      hisb_15_type2(logg[0],logg[1],logg[5],logg[4],logg[7],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[7]))) {
	set_dist(d,natoms,logg[0],logg[7],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      if (bVir) {
	/* VIRTUAL DISTANCES */
	nVdist = set_virtual (logg,9,ring_margin,d,natoms);
	ndist += nVdist;
      }
     }
    logg[10]=0;
    }
  }

  fprintf(log,"There are %d new histidine distances\n",ndist);

  if (ndist > 0 ) {
  fprintf(log,"(%d 1-2, %d 1-3, %d 1-4, %d 1-5 \n",n12dist,n13dist,
	  n14dist,n15dist);
  }
}

/**********************************************************
 *
 *     I S O L E U C I N E
 *
 **********************************************************/
void ile (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	  real ile_margin,t_ilist ilist[],t_iparams iparams[])
{
  /* Directly based on val.c */
  int natoms,ndist,i,j,logg[15],residnr,oldresidnr,n14dist,n15dist;
  real blen,lb,ub,angle;
  real pi = M_PI;

  natoms= atoms->nr;
  ndist   = 0;
  n14dist = 0;
  n15dist = 0;
  residnr = -1;
  oldresidnr=-1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),"ILE") == 0) {
      oldresidnr=residnr;
      residnr=atoms->atom[i].resnr;
      logg[14]=0;
      if ( oldresidnr == residnr ) {
	continue;
      }
      j=i;
      while ((atoms->atom[j].resnr) == residnr ) {
	if ( strcmp((*atoms->atomname[j]),"CA") == 0) {
	  logg[14]++;
	  logg[0]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CB") == 0) {
	  logg[14]++;
	  logg[1]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"HB") == 0) {
	  logg[14]++;
	  logg[2]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CD") == 0) {
	  logg[14]++;
	  logg[3]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"HD1") == 0) {
	  logg[14]++;
	  logg[4]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HD2") == 0) {
	  logg[14]++;
	  logg[5]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HD3") == 0) {
	  logg[14]++;
	  logg[6]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"CG1") == 0) {
	  logg[14]++;
	  logg[7]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HG11") == 0) {
	  logg[14]++;
	  logg[8]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HG12") == 0) {
	  logg[14]++;
	  logg[9]=j;
	}	
        else if ( strcmp((*atoms->atomname[j]),"CG2") == 0) {
	  logg[14]++;
	  logg[10]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HG21") == 0) {
	  logg[14]++;
	  logg[11]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HG22") == 0) {
	  logg[14]++;
	  logg[12]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HG23") == 0) {
	  logg[14]++;
	  logg[13]=j;
	}	
	if ( logg[14] == 14) {
	  break;
	}
	j++;
      }
    if ( logg[14] == 14 ) {
      fprintf(log,"logg (ile) =%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n"
	      ,logg[0],logg[1],logg[2],logg[3],logg[4],logg[5],logg[6],logg[7],
	      logg[8],logg[9],logg[10],logg[11],logg[12],logg[13]);

      /*SETDISTANCE for HD1 and CB */
      gauche(logg[4],logg[3],logg[7],logg[1],ilist,iparams,&blen,atoms);
      lb=(1.0-ile_margin)*blen;
      ub=(1.0+ile_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[1]))) {
	set_dist(d,natoms,logg[4],logg[1],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD1 and HG12 */
      gauche(logg[4],logg[3],logg[7],logg[9],ilist,iparams,&blen,atoms);
      lb=(1.0-ile_margin)*blen;
      ub=(1.0+ile_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[9]))) {
	set_dist(d,natoms,logg[4],logg[9],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD1 and HG11 */
      pdih_lengths(logg[4],logg[3],logg[7],logg[8],ilist,iparams,
		   &lb,&blen,atoms);
      lb=(1.0-ile_margin)*blen;
      ub=(1.0+ile_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[8]))) {
	set_dist(d,natoms,logg[4],logg[8],lb,ub,blen);
	ndist++;
	n14dist++;
      }

      /*SETDISTANCE for HD2 and CB  */
      gauche(logg[5],logg[3],logg[7],logg[1],ilist,iparams,&blen,atoms);
      lb=(1.0-ile_margin)*blen;
      ub=(1.0+ile_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[1]))) {
	set_dist(d,natoms,logg[5],logg[1],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD2 and HG11 */
      gauche(logg[5],logg[3],logg[7],logg[8],ilist,iparams,&blen,atoms);
      lb=(1.0-ile_margin)*blen;
      ub=(1.0+ile_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[8]))) {
	set_dist(d,natoms,logg[5],logg[8],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD2 and HG12 */
      pdih_lengths(logg[5],logg[3],logg[7],logg[9],ilist,iparams,
		   &lb,&blen,atoms);
      lb=(1.0-ile_margin)*blen;
      ub=(1.0+ile_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[9]))) {
	set_dist(d,natoms,logg[5],logg[9],lb,ub,blen);
	ndist++;
	n14dist++;
      }

      /*SETDISTANCE for HD3 and HG12 */
      gauche(logg[6],logg[3],logg[7],logg[9],ilist,iparams,&blen,atoms);
      lb=(1.0-ile_margin)*blen;
      ub=(1.0+ile_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[9]))) {
	set_dist(d,natoms,logg[6],logg[9],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD3 and HG11 */
      gauche(logg[6],logg[3],logg[7],logg[8],ilist,iparams,&blen,atoms);
      lb=(1.0-ile_margin)*blen;
      ub=(1.0+ile_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[8]))) {
	set_dist(d,natoms,logg[6],logg[8],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD3 and CB  */
      pdih_lengths(logg[6],logg[3],logg[7],logg[1],ilist,iparams,
		   &lb,&blen,atoms);
      lb=(1.0-ile_margin)*blen;
      ub=(1.0+ile_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[1]))) {
	set_dist(d,natoms,logg[6],logg[1],lb,ub,blen);
	ndist++;
	n14dist++;
      }

      /*SETDISTANCE for HG21 and CG1 */
      gauche(logg[11],logg[10],logg[1],logg[7],ilist,iparams,&blen,atoms);
      lb=(1.0-ile_margin)*blen;
      ub=(1.0+ile_margin)*blen;
      if (((weight[logg[11]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[11],logg[7]))) {
	set_dist(d,natoms,logg[11],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG21 and HB   */
      gauche(logg[11],logg[10],logg[1],logg[2],ilist,iparams,&blen,atoms);
      lb=(1.0-ile_margin)*blen;
      ub=(1.0+ile_margin)*blen;
      if (((weight[logg[11]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[11],logg[2]))) {
	set_dist(d,natoms,logg[11],logg[2],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG21 and CA   */
      pdih_lengths(logg[11],logg[10],logg[1],logg[0],ilist,iparams,
		   &lb,&blen,atoms);
      lb=(1.0-ile_margin)*blen;
      ub=(1.0+ile_margin)*blen;
      if (((weight[logg[11]] != 0.0) || (weight[logg[0]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[11],logg[0]))) {
	set_dist(d,natoms,logg[11],logg[0],lb,ub,blen);
	ndist++;
	n14dist++;
      }

      /*SETDISTANCE for HG22 and CG1 */
      gauche(logg[12],logg[10],logg[1],logg[7],ilist,iparams,&blen,atoms);
      lb=(1.0-ile_margin)*blen;
      ub=(1.0+ile_margin)*blen;
      if (((weight[logg[12]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[12],logg[7]))) {
	set_dist(d,natoms,logg[12],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG22 and CA   */
      gauche(logg[12],logg[10],logg[1],logg[0],ilist,iparams,&blen,atoms);
      lb=(1.0-ile_margin)*blen;
      ub=(1.0+ile_margin)*blen;
      if (((weight[logg[12]] != 0.0) || (weight[logg[0]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[12],logg[0]))) {
	set_dist(d,natoms,logg[12],logg[0],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG22 and HB   */
      pdih_lengths(logg[12],logg[10],logg[1],logg[2],ilist,iparams,
		   &lb,&blen,atoms);
      lb=(1.0-ile_margin)*blen;
      ub=(1.0+ile_margin)*blen;
      if (((weight[logg[12]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[12],logg[2]))) {
	set_dist(d,natoms,logg[12],logg[2],lb,ub,blen);
	ndist++;
	n14dist++;
      }

      /*SETDISTANCE for HG23 and CA   */
      gauche(logg[13],logg[10],logg[1],logg[0],ilist,iparams,&blen,atoms);
      lb=(1.0-ile_margin)*blen;
      ub=(1.0+ile_margin)*blen;
      if (((weight[logg[13]] != 0.0) || (weight[logg[0]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[13],logg[0]))) {
	set_dist(d,natoms,logg[13],logg[0],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG23 and HB   */
      gauche(logg[13],logg[10],logg[1],logg[2],ilist,iparams,&blen,atoms);
      lb=(1.0-ile_margin)*blen;
      ub=(1.0+ile_margin)*blen;
      if (((weight[logg[13]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[13],logg[2]))) {
	set_dist(d,natoms,logg[13],logg[2],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG23 and CG1*/
      pdih_lengths(logg[13],logg[10],logg[1],logg[7],ilist,iparams,
		   &lb,&blen,atoms);
      lb=(1.0-ile_margin)*blen;
      ub=(1.0+ile_margin)*blen;
      if (((weight[logg[13]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[13],logg[7]))) {
	set_dist(d,natoms,logg[13],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }

    }
    logg[14]=0;
    }
  }
  fprintf(log,"There are %d distances to keep isoleucine gauche\n",ndist);
  if (ndist > 0 ) {
    fprintf(log,"(%d 1-4, %d 1-5)\n",n14dist,n15dist);
  }
}

/**********************************************************
 *
 *     L E U C I N E
 *
 **********************************************************/
void leu (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	  real leu_margin,t_ilist ilist[],t_iparams iparams[])
{
  /* Directly based on val.c */
  int natoms,ndist,i,j,logg[12],residnr,oldresidnr,n14dist,n15dist;
  real blen,lb,ub,angle;
  real pi = M_PI;

  natoms= atoms->nr;
  ndist   = 0;
  n14dist = 0;
  n15dist = 0;
  residnr = -1;
  oldresidnr=-1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),"LEU") == 0) {
      oldresidnr=residnr;
      residnr=atoms->atom[i].resnr;
      logg[11]=0;
      if ( oldresidnr == residnr ) {
	continue;
      }
      j=i;
      while ((atoms->atom[j].resnr) == residnr ) {
	if ( strcmp((*atoms->atomname[j]),"CB") == 0) {
	  logg[11]++;
	  logg[0]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CG") == 0) {
	  logg[11]++;
	  logg[1]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"HG") == 0) {
	  logg[11]++;
	  logg[2]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CD1") == 0) {
	  logg[11]++;
	  logg[3]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"HD11") == 0) {
	  logg[11]++;
	  logg[4]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HD12") == 0) {
	  logg[11]++;
	  logg[5]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HD13") == 0) {
	  logg[11]++;
	  logg[6]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"CD2") == 0) {
	  logg[11]++;
	  logg[7]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HD21") == 0) {
	  logg[11]++;
	  logg[8]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HD22") == 0) {
	  logg[11]++;
	  logg[9]=j;
	}	
        else if ( strcmp((*atoms->atomname[j]),"HD23") == 0) {
	  logg[11]++;
	  logg[10]=j;
	}	
	if ( logg[11] == 11) {
	  break;
	}
	j++;
      }
    if ( logg[11] == 11 ) {
      fprintf(log,"logg (leu) = %d %d %d %d %d %d %d %d %d %d %d\n",
	      logg[0],logg[1],logg[2],logg[3],logg[4],logg[5],logg[6],
	      logg[7],logg[8],logg[9],logg[10]);
      /*SETDISTANCE for HD11 and CB */
      gauche(logg[4],logg[3],logg[1],logg[0],ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[0]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[0]))) {
	set_dist(d,natoms,logg[4],logg[0],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD11 and CD2 */
      gauche(logg[4],logg[3],logg[1],logg[7],ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[7]))) {
	set_dist(d,natoms,logg[4],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD11 and HG */
      pdih_lengths(logg[4],logg[3],logg[1],logg[2],ilist,iparams,
		   &lb,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[2]))) {
	set_dist(d,natoms,logg[4],logg[2],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD12 and CD2 */
      gauche(logg[5],logg[3],logg[1],logg[7],ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[7]))) {
	set_dist(d,natoms,logg[5],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD12 and HG */
      gauche(logg[5],logg[3],logg[1],logg[2],ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[2]))) {
	set_dist(d,natoms,logg[5],logg[2],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD12 and CB */
      pdih_lengths(logg[5],logg[3],logg[1],logg[0],ilist,iparams,
		   &lb,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[0]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[0]))) {
	set_dist(d,natoms,logg[5],logg[0],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD13 and HG */
      gauche(logg[6],logg[3],logg[1],logg[2],ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[2]))) {
	set_dist(d,natoms,logg[6],logg[2],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD13 and CB */
      gauche(logg[6],logg[3],logg[1],logg[0],ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[0]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[0]))) {
	set_dist(d,natoms,logg[6],logg[0],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD13 and CD2 */
      pdih_lengths(logg[6],logg[3],logg[1],logg[7],ilist,iparams,
		   &lb,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[7]))) {
	set_dist(d,natoms,logg[6],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD21 and HG */
      gauche(logg[8],logg[7],logg[1],logg[2],ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[2]))) {
	set_dist(d,natoms,logg[8],logg[2],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD21 and CD1 */
      gauche(logg[8],logg[7],logg[1],logg[3],ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[3]))) {
	set_dist(d,natoms,logg[8],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD21 and CB */
      pdih_lengths(logg[8],logg[7],logg[1],logg[0],ilist,iparams,
		   &lb,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[0]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[0]))) {
	set_dist(d,natoms,logg[8],logg[0],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD22 and CD1 */
      gauche(logg[9],logg[7],logg[1],logg[3],ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[3]))) {
	set_dist(d,natoms,logg[9],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD22 and CB */
      gauche(logg[9],logg[7],logg[1],logg[0],ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[0]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[0]))) {
	set_dist(d,natoms,logg[9],logg[0],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD22 and HG */
      pdih_lengths(logg[9],logg[7],logg[1],logg[2],ilist,iparams,
		   &lb,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[2]))) {
	set_dist(d,natoms,logg[9],logg[2],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD23 and CB */
      gauche(logg[10],logg[7],logg[1],logg[0],ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[0]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[0]))) {
	set_dist(d,natoms,logg[10],logg[0],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD23 and HG */
      gauche(logg[10],logg[7],logg[1],logg[2],ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[2]))) {
	set_dist(d,natoms,logg[10],logg[2],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD23 and CD1 */
      pdih_lengths(logg[10],logg[7],logg[1],logg[3],ilist,iparams,
		   &lb,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[3]))) {
	set_dist(d,natoms,logg[10],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }

      /* 111111111-----------5555555555 */
      /*SETDISTANCE for HD11 and HD21 */
      gauche15(logg[4],logg[3],logg[1],logg[7],logg[8],pi,pi+(pi/3.0),
	       pi+(pi/3.0),ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[8]))) {
	set_dist(d,natoms,logg[4],logg[8],lb,ub,blen);
	ndist++;
	n15dist++;
      }
      /*SETDISTANCE for HD11 and HD22 */
      gauche15(logg[4],logg[3],logg[1],logg[7],logg[9],pi,pi+pi/3.0,
	       pi-(pi/3.0),ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[9]))) {
	set_dist(d,natoms,logg[4],logg[9],lb,ub,blen);
	ndist++;
	n15dist++;
      }
      /*SETDISTANCE for HD11 and HD23 */
      gauche15(logg[4],logg[3],logg[1],logg[7],logg[10],pi,pi+(pi/3.0),0,
	       ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[10]))) {
	set_dist(d,natoms,logg[4],logg[10],lb,ub,blen);
	ndist++;
	n15dist++;
      }


      /*SETDISTANCE for HD12 and HD21 */
      gauche15(logg[5],logg[3],logg[1],logg[7],logg[8],pi,pi-pi/3.0,
	       pi+pi/3.0,ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[8]))) {
	set_dist(d,natoms,logg[5],logg[8],lb,ub,blen);
	ndist++;
	n15dist++;
      }
      /*SETDISTANCE for HD12 and HD22 */
      gauche15(logg[5],logg[3],logg[1],logg[7],logg[9],pi,pi-pi/3.0,
	       pi-pi/3.0,ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[9]))) {
	set_dist(d,natoms,logg[5],logg[9],lb,ub,blen);
	ndist++;
	n15dist++;
      }
      /*SETDISTANCE for HD12 and HD23 */
      gauche15(logg[5],logg[3],logg[1],logg[7],logg[10],pi,pi-pi/3.0,0,ilist,
	       iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[10]))) {
	set_dist(d,natoms,logg[5],logg[10],lb,ub,blen);
	ndist++;
	n15dist++;
      }

      /*SETDISTANCE for HD13 and HD21 */
      gauche15(logg[6],logg[3],logg[1],logg[7],logg[8],pi,0,pi+pi/3.0,
	       ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[8]))) {
	set_dist(d,natoms,logg[6],logg[8],lb,ub,blen);
	ndist++;
	n15dist++;
      }
      /*SETDISTANCE for HD13 and HD22 */
      gauche15(logg[6],logg[3],logg[1],logg[7],logg[9],pi,0,pi-pi/3.0,
	       ilist,iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[9]))) {
	set_dist(d,natoms,logg[6],logg[9],lb,ub,blen);
	ndist++;
	n15dist++;
      }
      /*SETDISTANCE for HD13 and HD23 */
      gauche15(logg[6],logg[3],logg[1],logg[7],logg[10],pi,0,0,ilist,
	       iparams,&blen,atoms);
      lb=(1.0-leu_margin)*blen;
      ub=(1.0+leu_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[10]))) {
	set_dist(d,natoms,logg[6],logg[10],lb,ub,blen);
	ndist++;
	n15dist++;
      }


    }
    logg[11]=0;
    }
  }
  fprintf(log,"There are %d distances to keep leucine gauche\n",ndist);
  if (ndist > 0 ) {
    fprintf(log,"(%d 1-4, %d 1-5)\n",n14dist,n15dist);
  }
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
  /*fprintf(stderr,"leaving routine\n");*/
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

  /*fprintf(stderr,"entering the routine\n");*/

  rij   = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  rjk   = lookup_bondlength(aj,ak,ilist,iparams,TRUE,atoms);
  rkl   = lookup_bondlength(ak,al,ilist,iparams,TRUE,atoms);
  rlm   = lookup_bondlength(al,am,ilist,iparams,TRUE,atoms);
  thijk = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  thjkl = lookup_angle(aj,ak,al,ilist,iparams,atoms);
  thklm = lookup_angle(ak,al,am,ilist,iparams,atoms);
  rik   = sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(thijk));
  rkm   = sqrt(rkl*rkl+rlm*rlm-2.0*rkl*rlm*cos(thklm));
  /*fprintf(stderr,"Only calculations left\n");*/

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




void phe (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	  real ring_margin,t_ilist ilist[],t_iparams iparams[],bool bVir)
{
  int natoms,ndist,i,j,q,logg[14],residnr,oldresidnr,n12dist,n13dist,n14dist,
    n15dist,n16dist,nVdist;
  real blen,lb,ub,angle;

  natoms= atoms->nr;
  ndist   = 0;
  n12dist = 0;
  n13dist = 0;
  n14dist = 0;
  n15dist = 0;
  n16dist = 0;
  nVdist  = 0;
  residnr = -1;
  oldresidnr=-1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),"PHE") == 0) {
      oldresidnr=residnr;
      residnr=atoms->atom[i].resnr;
      logg[13]=0;
      if ( oldresidnr == residnr ) {
	continue;
      }
      j=i;
      while ((atoms->atom[j].resnr) == residnr ) {
	if ( strcmp((*atoms->atomname[j]),"CB") == 0) {
	  logg[13]++;
	  logg[0]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CG") == 0) {
	  logg[13]++;
	  logg[1]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CD2") == 0) {
	  logg[13]++;
	  logg[2]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"HD2") == 0) {
	  logg[13]++;
	  logg[3]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CE2") == 0) {
	  logg[13]++;
	  logg[4]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HE2") == 0) {
	  logg[13]++;
	  logg[5]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"CZ") == 0) {
	  logg[13]++;
	  logg[6]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HZ") == 0) {
	  logg[13]++;
	  logg[7]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"CE1") == 0) {
	  logg[13]++;
	  logg[8]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HE1") == 0) {
	  logg[13]++;
	  logg[9]=j;
	}	
        else if ( strcmp((*atoms->atomname[j]),"CD1") == 0) {
	  logg[13]++;
	  logg[10]=j;
	}	
        else if ( strcmp((*atoms->atomname[j]),"HD1") == 0) {
	  logg[13]++;
	  logg[11]=j;
	}	
        else if ( (strcmp((*atoms->atomname[j]),"VF") == 0) && bVir) {
	  logg[13]++;
	  logg[12]=j;
	}	
	if ( ((logg[13] == 12) && !bVir) || ((logg[13] == 13) && bVir) ) {
	  break;
	}
	j++;
      }
    if ( ((logg[13] == 12) && !bVir) || ((logg[13] == 13) && bVir) ) {

      if ((logg[13] == 12) && !bVir) {
      fprintf(log,"logg (phe) = %d %d %d %d %d %d %d %d %d %d %d %d\n",
	      logg[0],logg[1],logg[2],logg[3],logg[4],logg[5],logg[6],
	      logg[7],logg[8],logg[9],logg[10],logg[11]);
      }
      else if ((logg[13] == 13) && bVir) {
      fprintf(log,"logg (phe) = %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
	      logg[0],logg[1],logg[2],logg[3],logg[4],logg[5],logg[6],
	      logg[7],logg[8],logg[9],logg[10],logg[11],logg[12]);
      }

      /*SETDISTANCE for CB and CG (1-2) */
      blen=lookup_bondlength(logg[0],logg[1],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[1]))) {
	set_dist(d,natoms,logg[0],logg[1],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CG and CD2 (1-2) */
      blen=lookup_bondlength(logg[1],logg[2],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[2]))) {
	set_dist(d,natoms,logg[1],logg[2],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CD2 and HD2 (1-2) */
      blen=lookup_bondlength(logg[2],logg[3],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[3]))) {
	set_dist(d,natoms,logg[2],logg[3],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CD2 and CE2 (1-2) */
      blen=lookup_bondlength(logg[2],logg[4],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[4]))) {
	set_dist(d,natoms,logg[2],logg[4],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CE2 and HE2 (1-2) */
      blen=lookup_bondlength(logg[4],logg[5],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[5]))) {
	set_dist(d,natoms,logg[4],logg[5],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CE2 and CZ (1-2) */
      blen=lookup_bondlength(logg[4],logg[6],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[6]))) {
	set_dist(d,natoms,logg[4],logg[6],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CZ and HZ (1-2) */
      blen=lookup_bondlength(logg[6],logg[7],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[7]))) {
	set_dist(d,natoms,logg[6],logg[7],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CZ and CE1 (1-2) */
      blen=lookup_bondlength(logg[6],logg[8],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[8]))) {
	set_dist(d,natoms,logg[6],logg[8],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CE1 and HE1 (1-2) */
      blen=lookup_bondlength(logg[8],logg[9],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[9]))) {
	set_dist(d,natoms,logg[8],logg[9],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CE1 and CD1 (1-2) */
      blen=lookup_bondlength(logg[8],logg[10],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[10]))) {
	set_dist(d,natoms,logg[8],logg[10],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CD1 and HD1 (1-2) */
      blen=lookup_bondlength(logg[10],logg[11],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[11]))) {
	set_dist(d,natoms,logg[10],logg[11],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CD1 and CG (1-2) */
      blen=lookup_bondlength(logg[1],logg[10],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[10]))) {
	set_dist(d,natoms,logg[1],logg[10],lb,ub,blen);
	ndist++;
	n12dist++;
      }

      /* new 1-3 distances added 981126 */
      /*SETDISTANCE for CG and HD1 (1-3) */
      angle=lookup_angle(logg[1],logg[10],logg[11],ilist,iparams,atoms);
      blen=angle_length(logg[1],logg[10],logg[11],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[11]))) {
	set_dist(d,natoms,logg[1],logg[11],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /* end of new 1-3 distances */

      /*SETDISTANCE for CD1 and CD2 (1-3) */
      angle=lookup_angle(logg[10],logg[1],logg[2],ilist,iparams,atoms);
      blen=angle_length(logg[10],logg[1],logg[2],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[2]))) {
	set_dist(d,natoms,logg[10],logg[2],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CG and CE2 (1-3) */
      angle=lookup_angle(logg[1],logg[2],logg[4],ilist,iparams,atoms);
      blen=angle_length(logg[1],logg[2],logg[4],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[4]))) {
	set_dist(d,natoms,logg[1],logg[4],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CD2 and CZ (1-3) */
      angle=lookup_angle(logg[2],logg[4],logg[6],ilist,iparams,atoms);
      blen=angle_length(logg[2],logg[4],logg[6],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[6]))) {
	set_dist(d,natoms,logg[2],logg[6],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CE2 and CE1 (1-3) */
      angle=lookup_angle(logg[4],logg[6],logg[8],ilist,iparams,atoms);
      blen=angle_length(logg[4],logg[6],logg[8],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[8]))) {
	set_dist(d,natoms,logg[4],logg[8],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CZ and CE1 (1-3) */
      angle=lookup_angle(logg[6],logg[8],logg[10],ilist,iparams,atoms);
      blen=angle_length(logg[6],logg[8],logg[10],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[10]))) {
	set_dist(d,natoms,logg[6],logg[10],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CE1 and CG (1-3) */
      angle=lookup_angle(logg[8],logg[10],logg[1],ilist,iparams,atoms);
      blen=angle_length(logg[8],logg[10],logg[1],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[1]))) {
	set_dist(d,natoms,logg[8],logg[1],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CB and CD1 (1-3) */
      angle=lookup_angle(logg[0],logg[1],logg[10],ilist,iparams,atoms);
      blen=angle_length(logg[0],logg[1],logg[10],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[10]))) {
	set_dist(d,natoms,logg[0],logg[10],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CB and CD2 (1-3) */
      angle=lookup_angle(logg[0],logg[1],logg[2],ilist,iparams,atoms);
      blen=angle_length(logg[0],logg[1],logg[2],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[2]))) {
	set_dist(d,natoms,logg[0],logg[2],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HD2 and CG (1-3) */
      angle=lookup_angle(logg[3],logg[2],logg[1],ilist,iparams,atoms);
      blen=angle_length(logg[3],logg[2],logg[1],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[1]))) {
	set_dist(d,natoms,logg[3],logg[1],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HD2 and CE2 (1-3) */
      angle=lookup_angle(logg[3],logg[2],logg[4],ilist,iparams,atoms);
      blen=angle_length(logg[3],logg[2],logg[4],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[4]))) {
	set_dist(d,natoms,logg[3],logg[4],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HE2 and CD2 (1-3) */
      angle=lookup_angle(logg[5],logg[4],logg[2],ilist,iparams,atoms);
      blen=angle_length(logg[5],logg[4],logg[2],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[2]))) {
	set_dist(d,natoms,logg[5],logg[2],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HE2 and CZ (1-3) */
      angle=lookup_angle(logg[5],logg[4],logg[6],ilist,iparams,atoms);
      blen=angle_length(logg[5],logg[4],logg[6],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[6]))) {
	set_dist(d,natoms,logg[5],logg[6],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HZ and CE2 (1-3) */
      angle=lookup_angle(logg[7],logg[6],logg[4],ilist,iparams,atoms);
      blen=angle_length(logg[7],logg[6],logg[4],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[4]))) {
	set_dist(d,natoms,logg[7],logg[4],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HZ and CE1 (1-3) */
      angle=lookup_angle(logg[7],logg[6],logg[8],ilist,iparams,atoms);
      blen=angle_length(logg[7],logg[6],logg[8],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[8]))) {
	set_dist(d,natoms,logg[7],logg[8],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HE1 and CZ (1-3) */
      angle=lookup_angle(logg[9],logg[8],logg[6],ilist,iparams,atoms);
      blen=angle_length(logg[9],logg[8],logg[6],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[6]))) {
	set_dist(d,natoms,logg[9],logg[6],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HE1 and CD1 (1-3) */
      angle=lookup_angle(logg[9],logg[8],logg[10],ilist,iparams,atoms);
      blen=angle_length(logg[9],logg[8],logg[10],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[10]))) {
	set_dist(d,natoms,logg[9],logg[10],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HD1 and CE1 (1-3) */
      angle=lookup_angle(logg[11],logg[10],logg[8],ilist,iparams,atoms);
      blen=angle_length(logg[11],logg[10],logg[8],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[11]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[11],logg[8]))) {
	set_dist(d,natoms,logg[11],logg[8],lb,ub,blen);
	ndist++;
	n13dist++;
      }

      /* new 1-4 distances added 981126 */
      /*SETDISTANCE for CD2 and HD1 (logg[2]&logg[11]) (transdihedral)*/
      pdih_lengths(logg[2],logg[1],logg[10],logg[11],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[11]))) {
	set_dist(d,natoms,logg[2],logg[11],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CZ and HD1 (logg[6]&logg[11]) (transdihedral)*/
      pdih_lengths(logg[6],logg[8],logg[10],logg[11],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[11]))) {
	set_dist(d,natoms,logg[6],logg[11],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CD2 and CE1 (logg[2]&logg[8]) (cisdihedral)*/
      pdih_lengths(logg[2],logg[1],logg[10],logg[8],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[8]))) {
	set_dist(d,natoms,logg[2],logg[8],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CE2 and CD1 (logg[4]&logg[10]) (cisdihedral)*/
      pdih_lengths(logg[4],logg[6],logg[8],logg[10],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[10]))) {
	set_dist(d,natoms,logg[4],logg[10],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /* end of new 1-4 distances */

      /*SETDISTANCE for CB and CE1 (logg[0]&logg[8]) (transdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[10],logg[8],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[8]))) {
	set_dist(d,natoms,logg[0],logg[8],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CB and CE2 (logg[0]&logg[4]) (transdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[2],logg[4],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[4]))) {
	set_dist(d,natoms,logg[0],logg[4],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CD1 and HD2 (logg[10]&logg[3]) (transdihedral)*/
      pdih_lengths(logg[10],logg[1],logg[2],logg[3],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[3]))) {
	set_dist(d,natoms,logg[10],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CZ and HD2 (logg[6]&logg[3]) (transdihedral)*/
      pdih_lengths(logg[6],logg[4],logg[2],logg[3],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[3]))) {
	set_dist(d,natoms,logg[6],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CG and HE2 (logg[1]&logg[5]) (transdihedral)*/
      pdih_lengths(logg[1],logg[2],logg[4],logg[5],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[5]))) {
	set_dist(d,natoms,logg[1],logg[5],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CE1 and HE2 (logg[8]&logg[5]) (transdihedral)*/
      pdih_lengths(logg[8],logg[6],logg[4],logg[5],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[5]))) {
	set_dist(d,natoms,logg[8],logg[5],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CD1 and HZ (logg[10]&logg[7]) (transdihedral)*/
      pdih_lengths(logg[10],logg[8],logg[6],logg[7],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[7]))) {
	set_dist(d,natoms,logg[10],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CD2 and HZ (logg[2]&logg[7]) (transdihedral)*/
      pdih_lengths(logg[2],logg[4],logg[6],logg[7],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[7]))) {
	set_dist(d,natoms,logg[2],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CE2 and HE1 (logg[4]&logg[9]) (transdihedral)*/
      pdih_lengths(logg[4],logg[6],logg[8],logg[9],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[9]))) {
	set_dist(d,natoms,logg[4],logg[9],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CG and HE1 (logg[1]&logg[9]) (transdihedral)*/
      pdih_lengths(logg[1],logg[10],logg[8],logg[9],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[9]))) {
	set_dist(d,natoms,logg[1],logg[9],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CB and HD1 (logg[0]&logg[11]) (cisdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[10],logg[11],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[11]))) {
	set_dist(d,natoms,logg[0],logg[11],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CB and HD2 (logg[0]&logg[3]) (cisdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[2],logg[3],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[3]))) {
	set_dist(d,natoms,logg[0],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CG and CZ (logg[1]&logg[6]) (cisdihedral)*/
      pdih_lengths(logg[1],logg[2],logg[4],logg[6],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[6]))) {
	set_dist(d,natoms,logg[1],logg[6],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD2 and HE2 (logg[3]&logg[5]) (cisdihedral)*/
      pdih_lengths(logg[3],logg[2],logg[4],logg[5],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[5]))) {
	set_dist(d,natoms,logg[3],logg[5],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HE2 and HZ (logg[5]&logg[7]) (cisdihedral)*/
      pdih_lengths(logg[5],logg[4],logg[6],logg[7],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[7]))) {
	set_dist(d,natoms,logg[5],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HZ and HE1 (logg[7]&logg[9]) (cisdihedral)*/
      pdih_lengths(logg[7],logg[6],logg[8],logg[9],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[9]))) {
	set_dist(d,natoms,logg[7],logg[9],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HE1 and HD1 (logg[9]&logg[11]) (cisdihedral)*/
      pdih_lengths(logg[9],logg[8],logg[10],logg[11],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[11]))) {
	set_dist(d,natoms,logg[9],logg[11],lb,ub,blen);
	ndist++;
	n14dist++;
      }

      /* new 1-5 distances added 981126 */
      /*SETDISTANCE for HD2 and HZ (1-5) */
      phe_15_CBHE(logg[3],logg[2],logg[4],logg[6],logg[7],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[7]))) {
	set_dist(d,natoms,logg[3],logg[7],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HD2 and HD1 (1-5) */
      phe_15_CBHE(logg[3],logg[2],logg[1],logg[10],logg[11],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[11]))) {
	set_dist(d,natoms,logg[3],logg[11],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HE2 and HE1 (1-5) */
      phe_15_CBHE(logg[5],logg[4],logg[6],logg[8],logg[9],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[9]))) {
	set_dist(d,natoms,logg[5],logg[9],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HZ and HD1 (1-5) */
      phe_15_CBHE(logg[7],logg[6],logg[8],logg[10],logg[11],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[11]))) {
	set_dist(d,natoms,logg[7],logg[11],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HE1 and CD2 (1-5) */
      phe_15_CBCZ(logg[9],logg[8],logg[6],logg[4],logg[2],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[2]))) {
	set_dist(d,natoms,logg[9],logg[2],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HD2 and CE1 (1-5) */
      phe_15_CBCZ(logg[3],logg[2],logg[4],logg[6],logg[8],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[8]))) {
	set_dist(d,natoms,logg[3],logg[8],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HD1 and CE2 (1-5) */
      phe_15_CBCZ(logg[11],logg[10],logg[1],logg[2],logg[4],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[11]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[11],logg[4]))) {
	set_dist(d,natoms,logg[11],logg[4],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HE2 and CD1 (1-5) */
      phe_15_CBCZ(logg[5],logg[4],logg[6],logg[8],logg[10],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[10]))) {
	set_dist(d,natoms,logg[5],logg[10],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /* end new 1-5 distances */

      /*fprintf(stderr,"Calling phe_15_CBHE\n");*/
      /*SETDISTANCE for CB and HE2 (1-5) */
      phe_15_CBHE(logg[0],logg[1],logg[2],logg[4],logg[5],ilist,
		      iparams,&blen,atoms);
      /*fprintf(stderr,"Exited phe_15_CBHE\n");*/
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[5]))) {
	set_dist(d,natoms,logg[0],logg[5],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*fprintf(stderr,"Calling phe_15_CBHE\n");*/     
      /*SETDISTANCE for CB and HE1 (1-5) */
      phe_15_CBHE(logg[0],logg[1],logg[10],logg[8],logg[9],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[9]))) {
	set_dist(d,natoms,logg[0],logg[9],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for CB and CZ (1-5) */
      phe_15_CBCZ(logg[0],logg[1],logg[2],logg[4],logg[6],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[6]))) {
	set_dist(d,natoms,logg[0],logg[6],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for CG and HZ (1-5) */
      phe_15_CGHZ(logg[1],logg[2],logg[4],logg[6],logg[7],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[7]))) {
	set_dist(d,natoms,logg[1],logg[7],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /* 1111111111-------------666666666666 */
      /*SETDISTANCE for CB and HZ (1-6) */
      phe_16_type2(logg[0],logg[1],logg[2],logg[4],logg[6],logg[7],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[7]))) {
	set_dist(d,natoms,logg[0],logg[7],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      /*SETDISTANCE for HD1 and HE2 (1-6) */
      phe_16_type2(logg[11],logg[10],logg[1],logg[2],logg[4],logg[5],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[11]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[11],logg[5]))) {
	set_dist(d,natoms,logg[11],logg[5],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      /*SETDISTANCE for HD2 and HE1 (1-6) */
      phe_16_type2(logg[3],logg[2],logg[4],logg[6],logg[8],logg[9],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[9]))) {
	set_dist(d,natoms,logg[3],logg[9],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      if (bVir) {
	/* VIRTUAL DISTANCES */
	nVdist += set_virtual (logg,12,ring_margin,d,natoms);
	ndist += nVdist;
      }
    }
    logg[13]=0;
    }
  }

  fprintf(log,"There are %d new phenylalanine distances\n",ndist);
  if (ndist > 0 ) {
  fprintf(log,"(%d 1-2, %d 1-3, %d 1-4, %d 1-5, %d 1-6, %d virtual)\n",n12dist,n13dist,
	  n14dist,n15dist,n16dist,nVdist);
  }
}

/**********************************************************
 *
 *      T Y R O S I N E
 *
 **********************************************************/
static void tyr_15_CBHE(int ai,int aj,int ak,int al,int am,
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
  /*fprintf(stderr,"leaving routine\n");*/
}

static void tyr_15_CBCZ(int ai,int aj,int ak,int al,int am,
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

static void tyr_15_CGOH(int ai,int aj,int ak,int al,int am,
			 t_ilist ilist[],t_iparams iparams[],
			 real *lb,t_atoms *atoms)
{
  /* Returns the length corresponding to a 1-5 distance */
  real rij,rjk,rkl,rlm,rik,rkm;
  real thijk,thjkl,thklm,thikj,thmkl;
  real real_pi = M_PI*2.0;

  /*fprintf(stderr,"entering the routine\n");*/

  rij   = lookup_bondlength(ai,aj,ilist,iparams,TRUE,atoms);
  rjk   = lookup_bondlength(aj,ak,ilist,iparams,TRUE,atoms);
  rkl   = lookup_bondlength(ak,al,ilist,iparams,TRUE,atoms);
  rlm   = lookup_bondlength(al,am,ilist,iparams,TRUE,atoms);
  thijk = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  thjkl = lookup_angle(aj,ak,al,ilist,iparams,atoms);
  thklm = lookup_angle(ak,al,am,ilist,iparams,atoms);
  rik   = sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(thijk));
  rkm   = sqrt(rkl*rkl+rlm*rlm-2.0*rkl*rlm*cos(thklm));
  /*fprintf(stderr,"Only calculations left\n");*/

  /* Compute angle thikj using law of sines */
  thikj = asin(rij*sin(thijk)/rik);
  
  /* Compute thmkl using law of sines */
  thmkl = asin(rlm*sin(thklm)/rkm);

  /* Compute cis length using law of cosines */
  *lb = sqrt(rik*rik+rkm*rkm-2.0*rik*rkm*cos(thjkl-thikj+thmkl));
}

static void tyr_16_type2(int ai,int aj,int ak,int al,int am,int an,
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

void tyr (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	  real ring_margin,t_ilist ilist[],t_iparams iparams[],bool bVir)
{
  int natoms,ndist,i,j,q,logg[14],residnr,oldresidnr,n12dist,n13dist,n14dist,
    n15dist,n16dist,nVdist;
  real blen,lb,ub,angle;

  natoms= atoms->nr;
  ndist   = 0;
  n12dist = 0;
  n13dist = 0;
  n14dist = 0;
  n15dist = 0;
  n16dist = 0;
  nVdist  = 0;
  residnr = -1;
  oldresidnr=-1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),"TYR") == 0) {
      oldresidnr=residnr;
      residnr=atoms->atom[i].resnr;
      logg[13]=0;
      if ( oldresidnr == residnr ) {
	continue;
      }
      j=i;
      while ((atoms->atom[j].resnr) == residnr ) {
	if ( strcmp((*atoms->atomname[j]),"CB") == 0) {
	  logg[13]++;
	  logg[0]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CG") == 0) {
	  logg[13]++;
	  logg[1]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CD2") == 0) {
	  logg[13]++;
	  logg[2]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"HD2") == 0) {
	  logg[13]++;
	  logg[3]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CE2") == 0) {
	  logg[13]++;
	  logg[4]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HE2") == 0) {
	  logg[13]++;
	  logg[5]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"CZ") == 0) {
	  logg[13]++;
	  logg[6]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"OH") == 0) {
	  logg[13]++;
	  logg[7]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"CE1") == 0) {
	  logg[13]++;
	  logg[8]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HE1") == 0) {
	  logg[13]++;
	  logg[9]=j;
	}	
        else if ( strcmp((*atoms->atomname[j]),"CD1") == 0) {
	  logg[13]++;
	  logg[10]=j;
	}	
        else if ( strcmp((*atoms->atomname[j]),"HD1") == 0) {
	  logg[13]++;
	  logg[11]=j;
	}	
        else if ( (strcmp((*atoms->atomname[j]),"VF") == 0) && bVir) {
	  logg[13]++;
	  logg[12]=j;
	}	
	if ( ((logg[13] == 12) && !bVir) || ((logg[13] == 13) && bVir) ) {
	  break;
	}
	j++;
      }
    if ( ((logg[13] == 12) && !bVir) || ((logg[13] == 13) && bVir) ) {

      if ((logg[13] == 12) && !bVir) {
	fprintf(log,"logg (tyr) = %d %d %d %d %d %d %d %d %d %d %d %d\n",
		logg[0],logg[1],logg[2],logg[3],logg[4],logg[5],logg[6],
		logg[7],logg[8],logg[9],logg[10],logg[11]);
      }
      else if ((logg[13] == 13) && bVir) {
	fprintf(log,"logg (tyr) = %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
		logg[0],logg[1],logg[2],logg[3],logg[4],logg[5],logg[6],
		logg[7],logg[8],logg[9],logg[10],logg[11],logg[12]);
      }

      /*SETDISTANCE for CB and CG (1-2) */
      blen=lookup_bondlength(logg[0],logg[1],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[1]))) {
	set_dist(d,natoms,logg[0],logg[1],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CG and CD2 (1-2) */
      blen=lookup_bondlength(logg[1],logg[2],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[2]))) {
	set_dist(d,natoms,logg[1],logg[2],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CD2 and HD2 (1-2) */
      blen=lookup_bondlength(logg[2],logg[3],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[3]))) {
	set_dist(d,natoms,logg[2],logg[3],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CD2 and CE2 (1-2) */
      blen=lookup_bondlength(logg[2],logg[4],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[4]))) {
	set_dist(d,natoms,logg[2],logg[4],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CE2 and HE2 (1-2) */
      blen=lookup_bondlength(logg[4],logg[5],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[5]))) {
	set_dist(d,natoms,logg[4],logg[5],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CE2 and CZ (1-2) */
      blen=lookup_bondlength(logg[4],logg[6],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[6]))) {
	set_dist(d,natoms,logg[4],logg[6],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CZ and OH (1-2) */
      blen=lookup_bondlength(logg[6],logg[7],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[7]))) {
	set_dist(d,natoms,logg[6],logg[7],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CZ and CE1 (1-2) */
      blen=lookup_bondlength(logg[6],logg[8],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[8]))) {
	set_dist(d,natoms,logg[6],logg[8],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CE1 and HE1 (1-2) */
      blen=lookup_bondlength(logg[8],logg[9],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[9]))) {
	set_dist(d,natoms,logg[8],logg[9],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CE1 and CD1 (1-2) */
      blen=lookup_bondlength(logg[8],logg[10],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[10]))) {
	set_dist(d,natoms,logg[8],logg[10],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CD1 and HD1 (1-2) */
      blen=lookup_bondlength(logg[10],logg[11],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[11]))) {
	set_dist(d,natoms,logg[10],logg[11],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CD1 and CG (1-2) */
      blen=lookup_bondlength(logg[1],logg[10],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[10]))) {
	set_dist(d,natoms,logg[1],logg[10],lb,ub,blen);
	ndist++;
	n12dist++;
      }

      /* new 1-3 distances added 981126 */
      /*SETDISTANCE for CG and HD1 (1-3) */
      angle=lookup_angle(logg[1],logg[10],logg[11],ilist,iparams,atoms);
      blen=angle_length(logg[1],logg[10],logg[11],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[11]))) {
	set_dist(d,natoms,logg[1],logg[11],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /* end of new 1-3 distances */

      /*SETDISTANCE for CD1 and CD2 (1-3) */
      angle=lookup_angle(logg[10],logg[1],logg[2],ilist,iparams,atoms);
      blen=angle_length(logg[10],logg[1],logg[2],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[2]))) {
	set_dist(d,natoms,logg[10],logg[2],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CG and CE2 (1-3) */
      angle=lookup_angle(logg[1],logg[2],logg[4],ilist,iparams,atoms);
      blen=angle_length(logg[1],logg[2],logg[4],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[4]))) {
	set_dist(d,natoms,logg[1],logg[4],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CD2 and CZ (1-3) */
      angle=lookup_angle(logg[2],logg[4],logg[6],ilist,iparams,atoms);
      blen=angle_length(logg[2],logg[4],logg[6],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[6]))) {
	set_dist(d,natoms,logg[2],logg[6],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CE2 and CE1 (1-3) */
      angle=lookup_angle(logg[4],logg[6],logg[8],ilist,iparams,atoms);
      blen=angle_length(logg[4],logg[6],logg[8],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[8]))) {
	set_dist(d,natoms,logg[4],logg[8],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CZ and CE1 (1-3) */
      angle=lookup_angle(logg[6],logg[8],logg[10],ilist,iparams,atoms);
      blen=angle_length(logg[6],logg[8],logg[10],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[10]))) {
	set_dist(d,natoms,logg[6],logg[10],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CE1 and CG (1-3) */
      angle=lookup_angle(logg[8],logg[10],logg[1],ilist,iparams,atoms);
      blen=angle_length(logg[8],logg[10],logg[1],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[1]))) {
	set_dist(d,natoms,logg[8],logg[1],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CB and CD1 (1-3) */
      angle=lookup_angle(logg[0],logg[1],logg[10],ilist,iparams,atoms);
      blen=angle_length(logg[0],logg[1],logg[10],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[10]))) {
	set_dist(d,natoms,logg[0],logg[10],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CB and CD2 (1-3) */
      angle=lookup_angle(logg[0],logg[1],logg[2],ilist,iparams,atoms);
      blen=angle_length(logg[0],logg[1],logg[2],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[2]))) {
	set_dist(d,natoms,logg[0],logg[2],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HD2 and CG (1-3) */
      angle=lookup_angle(logg[3],logg[2],logg[1],ilist,iparams,atoms);
      blen=angle_length(logg[3],logg[2],logg[1],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[1]))) {
	set_dist(d,natoms,logg[3],logg[1],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HD2 and CE2 (1-3) */
      angle=lookup_angle(logg[3],logg[2],logg[4],ilist,iparams,atoms);
      blen=angle_length(logg[3],logg[2],logg[4],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[4]))) {
	set_dist(d,natoms,logg[3],logg[4],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HE2 and CD2 (1-3) */
      angle=lookup_angle(logg[5],logg[4],logg[2],ilist,iparams,atoms);
      blen=angle_length(logg[5],logg[4],logg[2],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[2]))) {
	set_dist(d,natoms,logg[5],logg[2],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HE2 and CZ (1-3) */
      angle=lookup_angle(logg[5],logg[4],logg[6],ilist,iparams,atoms);
      blen=angle_length(logg[5],logg[4],logg[6],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[6]))) {
	set_dist(d,natoms,logg[5],logg[6],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for OH and CE2 (1-3) */
      angle=lookup_angle(logg[7],logg[6],logg[4],ilist,iparams,atoms);
      blen=angle_length(logg[7],logg[6],logg[4],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[4]))) {
	set_dist(d,natoms,logg[7],logg[4],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for OH and CE1 (1-3) */
      angle=lookup_angle(logg[7],logg[6],logg[8],ilist,iparams,atoms);
      blen=angle_length(logg[7],logg[6],logg[8],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[8]))) {
	set_dist(d,natoms,logg[7],logg[8],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HE1 and CZ (1-3) */
      angle=lookup_angle(logg[9],logg[8],logg[6],ilist,iparams,atoms);
      blen=angle_length(logg[9],logg[8],logg[6],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[6]))) {
	set_dist(d,natoms,logg[9],logg[6],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HE1 and CD1 (1-3) */
      angle=lookup_angle(logg[9],logg[8],logg[10],ilist,iparams,atoms);
      blen=angle_length(logg[9],logg[8],logg[10],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[10]))) {
	set_dist(d,natoms,logg[9],logg[10],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HD1 and CE1 (1-3) */
      angle=lookup_angle(logg[11],logg[10],logg[8],ilist,iparams,atoms);
      blen=angle_length(logg[11],logg[10],logg[8],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[11]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[11],logg[8]))) {
	set_dist(d,natoms,logg[11],logg[8],lb,ub,blen);
	ndist++;
	n13dist++;
      }

      /* new 1-4 distances added 981126 */
      /*SETDISTANCE for CD2 and HD1 (logg[2]&logg[11]) (transdihedral)*/
      pdih_lengths(logg[2],logg[1],logg[10],logg[11],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[11]))) {
	set_dist(d,natoms,logg[2],logg[11],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CZ and HD1 (logg[6]&logg[11]) (transdihedral)*/
      pdih_lengths(logg[6],logg[8],logg[10],logg[11],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[11]))) {
	set_dist(d,natoms,logg[6],logg[11],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CD2 and CE1 (logg[2]&logg[8]) (cisdihedral)*/
      pdih_lengths(logg[2],logg[1],logg[10],logg[8],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[8]))) {
	set_dist(d,natoms,logg[2],logg[8],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CE2 and CD1 (logg[4]&logg[10]) (cisdihedral)*/
      pdih_lengths(logg[4],logg[6],logg[8],logg[10],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[10]))) {
	set_dist(d,natoms,logg[4],logg[10],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /* end of new 1-4 distances */

      /*SETDISTANCE for CB and CE1 (logg[0]&logg[8]) (transdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[10],logg[8],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[8]))) {
	set_dist(d,natoms,logg[0],logg[8],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CB and CE2 (logg[0]&logg[4]) (transdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[2],logg[4],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[4]))) {
	set_dist(d,natoms,logg[0],logg[4],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CD1 and HD2 (logg[10]&logg[3]) (transdihedral)*/
      pdih_lengths(logg[10],logg[1],logg[2],logg[3],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[3]))) {
	set_dist(d,natoms,logg[10],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CZ and HD2 (logg[6]&logg[3]) (transdihedral)*/
      pdih_lengths(logg[6],logg[4],logg[2],logg[3],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[3]))) {
	set_dist(d,natoms,logg[6],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CG and HE2 (logg[1]&logg[5]) (transdihedral)*/
      pdih_lengths(logg[1],logg[2],logg[4],logg[5],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[5]))) {
	set_dist(d,natoms,logg[1],logg[5],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CE1 and HE2 (logg[8]&logg[5]) (transdihedral)*/
      pdih_lengths(logg[8],logg[6],logg[4],logg[5],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[5]))) {
	set_dist(d,natoms,logg[8],logg[5],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CD1 and OH (logg[10]&logg[7]) (transdihedral)*/
      pdih_lengths(logg[10],logg[8],logg[6],logg[7],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[7]))) {
	set_dist(d,natoms,logg[10],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CD2 and OH (logg[2]&logg[7]) (transdihedral)*/
      pdih_lengths(logg[2],logg[4],logg[6],logg[7],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[7]))) {
	set_dist(d,natoms,logg[2],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CE2 and HE1 (logg[4]&logg[9]) (transdihedral)*/
      pdih_lengths(logg[4],logg[6],logg[8],logg[9],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[9]))) {
	set_dist(d,natoms,logg[4],logg[9],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CG and HE1 (logg[1]&logg[9]) (transdihedral)*/
      pdih_lengths(logg[1],logg[10],logg[8],logg[9],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[9]))) {
	set_dist(d,natoms,logg[1],logg[9],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CB and HD1 (logg[0]&logg[11]) (cisdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[10],logg[11],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[11]))) {
	set_dist(d,natoms,logg[0],logg[11],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CB and HD2 (logg[0]&logg[3]) (cisdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[2],logg[3],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[3]))) {
	set_dist(d,natoms,logg[0],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CG and CZ (logg[1]&logg[6]) (cisdihedral)*/
      pdih_lengths(logg[1],logg[2],logg[4],logg[6],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[6]))) {
	set_dist(d,natoms,logg[1],logg[6],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HD2 and HE2 (logg[3]&logg[5]) (cisdihedral)*/
      pdih_lengths(logg[3],logg[2],logg[4],logg[5],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[5]))) {
	set_dist(d,natoms,logg[3],logg[5],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HE2 and OH (logg[5]&logg[7]) (cisdihedral)*/
      pdih_lengths(logg[5],logg[4],logg[6],logg[7],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[7]))) {
	set_dist(d,natoms,logg[5],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for OH and HE1 (logg[7]&logg[9]) (cisdihedral)*/
      pdih_lengths(logg[7],logg[6],logg[8],logg[9],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[9]))) {
	set_dist(d,natoms,logg[7],logg[9],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HE1 and HD1 (logg[9]&logg[11]) (cisdihedral)*/
      pdih_lengths(logg[9],logg[8],logg[10],logg[11],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[11]))) {
	set_dist(d,natoms,logg[9],logg[11],lb,ub,blen);
	ndist++;
	n14dist++;
      }

      /* new 1-5 distances added 981126 */
      /*SETDISTANCE for HD2 and OH (1-5) */
      tyr_15_CBHE(logg[3],logg[2],logg[4],logg[6],logg[7],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[7]))) {
	set_dist(d,natoms,logg[3],logg[7],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HD2 and HD1 (1-5) */
      tyr_15_CBHE(logg[3],logg[2],logg[1],logg[10],logg[11],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[11]))) {
	set_dist(d,natoms,logg[3],logg[11],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HE2 and HE1 (1-5) */
      tyr_15_CBHE(logg[5],logg[4],logg[6],logg[8],logg[9],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[9]))) {
	set_dist(d,natoms,logg[5],logg[9],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for OH and HD1 (1-5) */
      tyr_15_CBHE(logg[7],logg[6],logg[8],logg[10],logg[11],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[11]))) {
	set_dist(d,natoms,logg[7],logg[11],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HE1 and CD2 (1-5) */
      tyr_15_CBCZ(logg[9],logg[8],logg[6],logg[4],logg[2],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[2]))) {
	set_dist(d,natoms,logg[9],logg[2],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HD2 and CE1 (1-5) */
      tyr_15_CBCZ(logg[3],logg[2],logg[4],logg[6],logg[8],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[8]))) {
	set_dist(d,natoms,logg[3],logg[8],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HD1 and CE2 (1-5) */
      tyr_15_CBCZ(logg[11],logg[10],logg[1],logg[2],logg[4],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[11]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[11],logg[4]))) {
	set_dist(d,natoms,logg[11],logg[4],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HE2 and CD1 (1-5) */
      tyr_15_CBCZ(logg[5],logg[4],logg[6],logg[8],logg[10],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[10]))) {
	set_dist(d,natoms,logg[5],logg[10],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /* end new 1-5 distances */

      /*fprintf(stderr,"Calling tyr_15_CBHE\n");*/
      /*SETDISTANCE for CB and HE2 (1-5) */
      tyr_15_CBHE(logg[0],logg[1],logg[2],logg[4],logg[5],ilist,
		      iparams,&blen,atoms);
      /*fprintf(stderr,"Exited tyr_15_CBHE\n");*/
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[5]))) {
	set_dist(d,natoms,logg[0],logg[5],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*fprintf(stderr,"Calling tyr_15_CBHE\n");*/     
      /*SETDISTANCE for CB and HE1 (1-5) */
      tyr_15_CBHE(logg[0],logg[1],logg[10],logg[8],logg[9],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[9]))) {
	set_dist(d,natoms,logg[0],logg[9],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for CB and CZ (1-5) */
      tyr_15_CBCZ(logg[0],logg[1],logg[2],logg[4],logg[6],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[6]))) {
	set_dist(d,natoms,logg[0],logg[6],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*fprintf(stderr,"Calling tyr_15_CGOH\n");*/
      /*SETDISTANCE for CG and OH (1-5) */
      tyr_15_CGOH(logg[1],logg[2],logg[4],logg[6],logg[7],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[7]))) {
	set_dist(d,natoms,logg[1],logg[7],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /* 1111111111------------6666666666  */
      /*SETDISTANCE for CB and OH (1-6) */
      tyr_16_type2(logg[0],logg[1],logg[2],logg[4],logg[6],logg[7],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[7]))) {
	set_dist(d,natoms,logg[0],logg[7],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      /*SETDISTANCE for HD1 and HE2 (1-6) */
      tyr_16_type2(logg[11],logg[10],logg[1],logg[2],logg[4],logg[5],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[11]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[11],logg[5]))) {
	set_dist(d,natoms,logg[11],logg[5],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      /*SETDISTANCE for HD2 and HE1 (1-6) */
      tyr_16_type2(logg[3],logg[2],logg[4],logg[6],logg[8],logg[9],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[9]))) {
	set_dist(d,natoms,logg[3],logg[9],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      if (bVir) {
	/* VIRTUAL DISTANCES */
	nVdist += set_virtual (logg,12,ring_margin,d,natoms);
	ndist += nVdist;
	if (FALSE) {
	  for (q=0 ; q<12 ; q++ ) {
	    if (q == 10) {
	      blen = 8.0;
	      lb   = (1.0-ring_margin)*blen;
	      ub   = (1.0+ring_margin)*blen;
	    }
	    else {
	      blen = d_len(d,natoms,logg[q],logg[10]);
	      blen = sqrt(blen*blen+64.0);
	      lb   = (1.0-ring_margin)*blen;
	      ub   = (1.0+ring_margin)*blen;
	    }
        /* set distance to virtual atom */
	    set_dist(d,natoms,logg[q],logg[12],lb,ub,blen);
	    ndist++;
	    nVdist++;
	  }
	}

      }
      
    }
    logg[13]=0;
    }
  }
  
  fprintf(log,"There are %d new tyrosine distances\n",ndist);
  if (ndist > 0 ) {
    fprintf(log,"(%d 1-2, %d 1-3, %d 1-4, %d 1-5, %d 1-6, %d virtual)\n",n12dist,n13dist,
	  n14dist,n15dist,n16dist,nVdist);
  }
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
  /*fprintf(stderr,"leaving routine\n");*/
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

void trp (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	  real ring_margin,t_ilist ilist[],t_iparams iparams[],bool bVir)
{
  int natoms,ndist,i,j,q,logg[18],residnr,oldresidnr,n12dist,n13dist,n14dist,
    n15dist,n16dist,n17dist,nVdist;
  real blen,lb,ub,angle;

  natoms= atoms->nr;
  ndist   = 0;
  n12dist = 0;
  n13dist = 0;
  n14dist = 0;
  n15dist = 0;
  n16dist = 0;
  n17dist = 0;
  nVdist  = 0;
  residnr = -1;
  oldresidnr=-1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),"TRP") == 0) {
      oldresidnr=residnr;
      residnr=atoms->atom[i].resnr;
      logg[17]=0;
      if ( oldresidnr == residnr ) {
	continue;
      }
      j=i;
      while ((atoms->atom[j].resnr) == residnr ) {
	if ( strcmp((*atoms->atomname[j]),"CB") == 0) {
	  logg[17]++;
	  logg[0]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CG") == 0) {
	  logg[17]++;
	  logg[1]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CD2") == 0) {
	  logg[17]++;
	  logg[2]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CE3") == 0) {
	  logg[17]++;
	  logg[3]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CZ3") == 0) {
	  logg[17]++;
	  logg[4]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"CH2") == 0) {
	  logg[17]++;
	  logg[5]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"CZ2") == 0) {
	  logg[17]++;
	  logg[6]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"CE2") == 0) {
	  logg[17]++;
	  logg[7]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"NE1") == 0) {
	  logg[17]++;
	  logg[8]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"CD1") == 0) {
	  logg[17]++;
	  logg[9]=j;
	}	
        else if ( strcmp((*atoms->atomname[j]),"HE3") == 0) {
	  logg[17]++;
	  logg[10]=j;
	}	
        else if ( strcmp((*atoms->atomname[j]),"HZ3") == 0) {
	  logg[17]++;
	  logg[11]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HH2") == 0) {
	  logg[17]++;
	  logg[12]=j;
	}	
        else if ( strcmp((*atoms->atomname[j]),"HZ2") == 0) {
	  logg[17]++;
	  logg[13]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HE1") == 0) {
	  logg[17]++;
	  logg[14]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HD1") == 0) {
	  logg[17]++;
	  logg[15]=j;
	}
        else if ( (strcmp((*atoms->atomname[j]),"VF") == 0) && bVir) {
	  logg[17]++;
	  logg[16]=j;
	}
	if ( ((logg[17] == 16) && !bVir) || ((logg[17] == 17) && bVir) ) {
	  break;
	}
	j++;
      }
    if ( ((logg[17] == 16) && !bVir) || ((logg[17] == 17) && bVir) ) {

      if ((logg[17] == 16) && !bVir) {
      fprintf(log,"logg (trp) = %d %d %d %d %d %d %d %d %d %d %d %d",
	      logg[0],logg[1],logg[2],logg[3],logg[4],logg[5],logg[6],
	      logg[7],logg[8],logg[9],logg[10],logg[11]);
      fprintf(log," %d %d %d %d\n",logg[12],logg[13],logg[14],logg[15]);
      }
      else if  ((logg[17] == 17) && bVir) {
      fprintf(log,"logg (trp) = %d %d %d %d %d %d %d %d %d %d %d %d",
	      logg[0],logg[1],logg[2],logg[3],logg[4],logg[5],logg[6],
	      logg[7],logg[8],logg[9],logg[10],logg[11]);
      fprintf(log," %d %d %d %d %d\n",logg[12],logg[13],logg[14],logg[15],
	      logg[16]);
      }

      /*SETDISTANCE for CB and CG (1-2) */
      blen=lookup_bondlength(logg[0],logg[1],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[1]))) {
	set_dist(d,natoms,logg[0],logg[1],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CG and CD2 (1-2) */
      blen=lookup_bondlength(logg[1],logg[2],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[2]))) {
	set_dist(d,natoms,logg[1],logg[2],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CD2 and CE3 (1-2) */
      blen=lookup_bondlength(logg[2],logg[3],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[3]))) {
	set_dist(d,natoms,logg[2],logg[3],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CE3 and CZ3 (1-2) */
      blen=lookup_bondlength(logg[3],logg[4],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[4]))) {
	set_dist(d,natoms,logg[3],logg[4],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CZ3 and CH2 (1-2) */
      blen=lookup_bondlength(logg[4],logg[5],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[5]))) {
	set_dist(d,natoms,logg[4],logg[5],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CH2 and CZ2 (1-2) */
      blen=lookup_bondlength(logg[5],logg[6],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[6]))) {
	set_dist(d,natoms,logg[5],logg[6],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CZ2 and CE2 (1-2) */
      blen=lookup_bondlength(logg[6],logg[7],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[7]))) {
	set_dist(d,natoms,logg[6],logg[7],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CE2 and CD2 (1-2) */
      blen=lookup_bondlength(logg[7],logg[2],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[2]))) {
	set_dist(d,natoms,logg[7],logg[2],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CE2 and NE1 (1-2) */
      blen=lookup_bondlength(logg[7],logg[8],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[8]))) {
	set_dist(d,natoms,logg[7],logg[8],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for NE1 and CD1 (1-2) */
      blen=lookup_bondlength(logg[8],logg[9],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[9]))) {
	set_dist(d,natoms,logg[8],logg[9],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CD1 and CG (1-2) */
      blen=lookup_bondlength(logg[9],logg[1],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[1]))) {
	set_dist(d,natoms,logg[9],logg[1],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CE3 and HE3 (1-2) */
      blen=lookup_bondlength(logg[3],logg[10],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[10]))) {
	set_dist(d,natoms,logg[3],logg[10],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CZ3 and HZ3 (1-2) */
      blen=lookup_bondlength(logg[4],logg[11],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[11]))) {
	set_dist(d,natoms,logg[4],logg[11],lb,ub,blen);
	ndist++;
	n12dist++;
      }      
      /*SETDISTANCE for CH2 and HH2 (1-2) */
      blen=lookup_bondlength(logg[5],logg[12],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[12]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[12]))) {
	set_dist(d,natoms,logg[5],logg[12],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CZ2 and HZ2 (1-2) */
      blen=lookup_bondlength(logg[6],logg[13],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[13]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[13]))) {
	set_dist(d,natoms,logg[6],logg[13],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for NE1 and HE1 (1-2) */
      blen=lookup_bondlength(logg[8],logg[14],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[14]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[14]))) {
	set_dist(d,natoms,logg[8],logg[14],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CD1 and HD1 (1-2) */
      blen=lookup_bondlength(logg[9],logg[15],ilist,iparams,TRUE,atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[15]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[15]))) {
	set_dist(d,natoms,logg[9],logg[15],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for CB and CD1 (1-3) */
      angle=lookup_angle(logg[0],logg[1],logg[9],ilist,iparams,atoms);
      blen=angle_length(logg[0],logg[1],logg[9],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[9]))) {
	set_dist(d,natoms,logg[0],logg[9],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CB and CD2 (1-3) */
      angle=lookup_angle(logg[0],logg[1],logg[2],ilist,iparams,atoms);
      blen=angle_length(logg[0],logg[1],logg[2],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[2]))) {
	set_dist(d,natoms,logg[0],logg[2],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CG and CE3 (1-3) */
      angle=lookup_angle(logg[1],logg[2],logg[3],ilist,iparams,atoms);
      blen=angle_length(logg[1],logg[2],logg[3],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[3]))) {
	set_dist(d,natoms,logg[1],logg[3],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HE2 and CD2 (1-3) */
      angle=lookup_angle(logg[10],logg[3],logg[2],ilist,iparams,atoms);
      blen=angle_length(logg[10],logg[3],logg[2],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[2]))) {
	set_dist(d,natoms,logg[10],logg[2],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HE3 and CZ3 (1-3) */
      angle=lookup_angle(logg[10],logg[3],logg[4],ilist,iparams,atoms);
      blen=angle_length(logg[10],logg[3],logg[4],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[4]))) {
	set_dist(d,natoms,logg[10],logg[4],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HZ3 and CE3 (1-3) */
      angle=lookup_angle(logg[11],logg[4],logg[3],ilist,iparams,atoms);
      blen=angle_length(logg[11],logg[4],logg[3],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[11]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[11],logg[3]))) {
	set_dist(d,natoms,logg[11],logg[3],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HZ3 and CH2 (1-3) */
      angle=lookup_angle(logg[11],logg[4],logg[5],ilist,iparams,atoms);
      blen=angle_length(logg[11],logg[4],logg[5],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[11]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[11],logg[5]))) {
	set_dist(d,natoms,logg[11],logg[5],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HH2 and CZ3 (1-3) */
      angle=lookup_angle(logg[12],logg[5],logg[4],ilist,iparams,atoms);
      blen=angle_length(logg[12],logg[5],logg[4],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[12]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[12],logg[4]))) {
	set_dist(d,natoms,logg[12],logg[4],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HH2 and CZ2 (1-3) */
      angle=lookup_angle(logg[12],logg[5],logg[6],ilist,iparams,atoms);
      blen=angle_length(logg[12],logg[5],logg[6],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[12]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[12],logg[6]))) {
	set_dist(d,natoms,logg[12],logg[6],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HZ2 and CH2 (1-3) */
      angle=lookup_angle(logg[13],logg[6],logg[5],ilist,iparams,atoms);
      blen=angle_length(logg[13],logg[6],logg[5],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[13]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[13],logg[5]))) {
	set_dist(d,natoms,logg[13],logg[5],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HZ2 and CE2 (1-3) */
      angle=lookup_angle(logg[13],logg[6],logg[7],ilist,iparams,atoms);
      blen=angle_length(logg[13],logg[6],logg[7],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[13]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[13],logg[7]))) {
	set_dist(d,natoms,logg[13],logg[7],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CZ2 and NE1 (1-3) */
      angle=lookup_angle(logg[6],logg[7],logg[8],ilist,iparams,atoms);
      blen=angle_length(logg[6],logg[7],logg[8],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[8]))) {
	set_dist(d,natoms,logg[6],logg[8],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HE1 and CE2 (1-3) */
      angle=lookup_angle(logg[14],logg[8],logg[7],ilist,iparams,atoms);
      blen=angle_length(logg[14],logg[8],logg[7],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[14]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[14],logg[7]))) {
	set_dist(d,natoms,logg[14],logg[7],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HE1 and CD1 (1-3) */
      angle=lookup_angle(logg[14],logg[8],logg[9],ilist,iparams,atoms);
      blen=angle_length(logg[14],logg[8],logg[9],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[14]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[14],logg[9]))) {
	set_dist(d,natoms,logg[14],logg[9],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HD1 and NE1 (1-3) */
      angle=lookup_angle(logg[15],logg[9],logg[8],ilist,iparams,atoms);
      blen=angle_length(logg[15],logg[9],logg[8],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[15]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[15],logg[8]))) {
	set_dist(d,natoms,logg[15],logg[8],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for HD1 and CG (1-3) */
      angle=lookup_angle(logg[15],logg[9],logg[1],ilist,iparams,atoms);
      blen=angle_length(logg[15],logg[9],logg[1],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[15]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[15],logg[1]))) {
	set_dist(d,natoms,logg[15],logg[1],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CD1 and CD2 (1-3) */
      angle=lookup_angle(logg[9],logg[1],logg[2],ilist,iparams,atoms);
      blen=angle_length(logg[9],logg[1],logg[2],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[2]))) {
	set_dist(d,natoms,logg[9],logg[2],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CG and CE2 (1-3) */
      angle=lookup_angle(logg[1],logg[2],logg[7],ilist,iparams,atoms);
      blen=angle_length(logg[1],logg[2],logg[7],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[7]))) {
	set_dist(d,natoms,logg[1],logg[7],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CD2 and NE1 (1-3) */
      angle=lookup_angle(logg[2],logg[7],logg[8],ilist,iparams,atoms);
      blen=angle_length(logg[2],logg[7],logg[8],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[8]))) {
	set_dist(d,natoms,logg[2],logg[8],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CE2 and CD1 (1-3) */
      angle=lookup_angle(logg[7],logg[8],logg[9],ilist,iparams,atoms);
      blen=angle_length(logg[7],logg[8],logg[9],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[9]))) {
	set_dist(d,natoms,logg[7],logg[9],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for NE1 and CG (1-3) */
      angle=lookup_angle(logg[8],logg[9],logg[1],ilist,iparams,atoms);
      blen=angle_length(logg[8],logg[9],logg[1],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[1]))) {
	set_dist(d,natoms,logg[8],logg[1],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CE2 and CE3 (1-3) */
      angle=lookup_angle(logg[7],logg[2],logg[3],ilist,iparams,atoms);
      blen=angle_length(logg[7],logg[2],logg[3],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[3]))) {
	set_dist(d,natoms,logg[7],logg[3],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CD2 and CZ3 (1-3) */
      angle=lookup_angle(logg[2],logg[3],logg[4],ilist,iparams,atoms);
      blen=angle_length(logg[2],logg[3],logg[4],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[4]))) {
	set_dist(d,natoms,logg[2],logg[4],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CE3 and CH2 (1-3) */
      angle=lookup_angle(logg[3],logg[4],logg[5],ilist,iparams,atoms);
      blen=angle_length(logg[3],logg[4],logg[5],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[5]))) {
	set_dist(d,natoms,logg[3],logg[5],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CZ3 and CZ2 (1-3) */
      angle=lookup_angle(logg[4],logg[5],logg[6],ilist,iparams,atoms);
      blen=angle_length(logg[4],logg[5],logg[6],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[6]))) {
	set_dist(d,natoms,logg[4],logg[6],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CH2 and CE2 (1-3) */
      angle=lookup_angle(logg[5],logg[6],logg[7],ilist,iparams,atoms);
      blen=angle_length(logg[5],logg[6],logg[7],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[7]))) {
	set_dist(d,natoms,logg[5],logg[7],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CZ2 and CD2 (1-3) */
      angle=lookup_angle(logg[6],logg[7],logg[2],ilist,iparams,atoms);
      blen=angle_length(logg[6],logg[7],logg[2],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[2]))) {
	set_dist(d,natoms,logg[6],logg[2],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*if (1>5) {*/
      /*SETDISTANCE for CB and CE3 (logg[0]&logg[3]) (cisdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[2],logg[3],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[3]))) {
	set_dist(d,natoms,logg[0],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CB and HD1 (logg[0]&logg[15]) (cisdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[9],logg[15],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[15]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[15]))) {
	set_dist(d,natoms,logg[0],logg[15],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CG and HE3 (logg[1]&logg[10]) (cisdihedral)*/
      pdih_lengths(logg[1],logg[2],logg[3],logg[10],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[10]))) {
	set_dist(d,natoms,logg[1],logg[10],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HE3 and HZ3 (logg[10]&logg[11]) (cisdihedral)*/
      pdih_lengths(logg[10],logg[3],logg[4],logg[11],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[11]))) {
	set_dist(d,natoms,logg[10],logg[11],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HZ3 and HH2 (logg[11]&logg[12]) (cisdihedral)*/
      pdih_lengths(logg[11],logg[4],logg[5],logg[12],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[11]] != 0.0) || (weight[logg[12]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[11],logg[12]))) {
	set_dist(d,natoms,logg[11],logg[12],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HH2 and HZ2 (logg[12]&logg[13]) (cisdihedral)*/
      pdih_lengths(logg[12],logg[5],logg[6],logg[13],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[12]] != 0.0) || (weight[logg[13]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[12],logg[13]))) {
	set_dist(d,natoms,logg[12],logg[13],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HZ2 and NE1 (logg[13]&logg[8]) (cisdihedral)*/
      pdih_lengths(logg[13],logg[6],logg[7],logg[8],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[13]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[13],logg[8]))) {
	set_dist(d,natoms,logg[13],logg[8],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CZ2 and HE1 (logg[6]&logg[14]) (cisdihedral)*/
      pdih_lengths(logg[6],logg[7],logg[8],logg[14],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[14]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[14]))) {
	set_dist(d,natoms,logg[6],logg[14],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HE1 and HD1 (logg[14]&logg[15]) (cisdihedral)*/
      pdih_lengths(logg[14],logg[8],logg[9],logg[15],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[14]] != 0.0) || (weight[logg[15]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[14],logg[15]))) {
	set_dist(d,natoms,logg[14],logg[15],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CD2 and CH2 (logg[2]&logg[5]) (cisdihedral)*/
      pdih_lengths(logg[2],logg[3],logg[4],logg[5],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[5]))) {
	set_dist(d,natoms,logg[2],logg[5],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CE3 and CZ2 (logg[0]&logg[11]) (cisdihedral)*/
      pdih_lengths(logg[3],logg[4],logg[5],logg[6],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[6]))) {
	set_dist(d,natoms,logg[3],logg[6],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CZ3 and CE2 (logg[4]&logg[7]) (cisdihedral)*/
      pdih_lengths(logg[4],logg[5],logg[6],logg[7],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[7]))) {
	set_dist(d,natoms,logg[4],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /* trans 1-4's added 981124 */

      /*SETDISTANCE for CB and CE2 (logg[0]&logg[7]) (transdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[2],logg[7],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[7]))) {
	set_dist(d,natoms,logg[0],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CB and NE1 (logg[0]&logg[8]) (transdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[9],logg[8],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[8]))) {
	set_dist(d,natoms,logg[0],logg[8],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CG and HE1 (logg[1]&logg[14]) (transdihedral)*/
      pdih_lengths(logg[1],logg[9],logg[8],logg[14],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[14]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[14]))) {
	set_dist(d,natoms,logg[1],logg[14],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CD2 and HZ2 (logg[2]&logg[13]) (transdihedral)*/
      pdih_lengths(logg[2],logg[7],logg[6],logg[13],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[13]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[13]))) {
	set_dist(d,natoms,logg[2],logg[13],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CD2 and HE1 (logg[2]&logg[14]) (transdihedral)*/
      pdih_lengths(logg[2],logg[7],logg[8],logg[14],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[14]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[14]))) {
	set_dist(d,natoms,logg[2],logg[14],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CE3 and CD1 (logg[3]&logg[9]) (transdihedral)*/
      pdih_lengths(logg[3],logg[2],logg[1],logg[9],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[9]))) {
	set_dist(d,natoms,logg[3],logg[9],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CH2 and HE3 (logg[5]&logg[10]) (transdihedral)*/
      pdih_lengths(logg[5],logg[4],logg[3],logg[10],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[10]))) {
	set_dist(d,natoms,logg[5],logg[10],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CZ2 and HZ3 (logg[6]&logg[11]) (transdihedral)*/
      pdih_lengths(logg[6],logg[5],logg[4],logg[11],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[11]))) {
	set_dist(d,natoms,logg[6],logg[11],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CE2 and HE3 (logg[7]&logg[10]) (transdihedral)*/
      pdih_lengths(logg[7],logg[2],logg[3],logg[10],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[10]))) {
	set_dist(d,natoms,logg[7],logg[10],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CE2 and HH2 (logg[7]&logg[12]) (transdihedral)*/
      pdih_lengths(logg[7],logg[6],logg[5],logg[12],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[12]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[12]))) {
	set_dist(d,natoms,logg[7],logg[12],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /* end of 1-4's added 981124 */



      /*SETDISTANCE for CG and CZ3 (logg[1]&logg[4]) (transdihedral)*/
      pdih_lengths(logg[1],logg[2],logg[3],logg[4],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[4]))) {
	set_dist(d,natoms,logg[1],logg[4],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CD1 and CZ2 (logg[9]&logg[6]) (transdihedral)*/
      pdih_lengths(logg[9],logg[8],logg[7],logg[6],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[6]))) {
	set_dist(d,natoms,logg[9],logg[6],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for NE1 and CH2 (logg[8]&logg[5]) (transdihedral)*/
      pdih_lengths(logg[8],logg[7],logg[6],logg[5],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[5]))) {
	set_dist(d,natoms,logg[8],logg[5],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for NE1 and CE3 (logg[8]&logg[3]) (transdihedral)*/
      pdih_lengths(logg[8],logg[7],logg[2],logg[3],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[3]))) {
	set_dist(d,natoms,logg[8],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CD2 and HZ3 (logg[2]&logg[11]) (transdihedral)*/
      pdih_lengths(logg[2],logg[3],logg[4],logg[11],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[11]))) {
	set_dist(d,natoms,logg[2],logg[11],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CE3 and HH2 (logg[3]&logg[12]) (transdihedral)*/
      pdih_lengths(logg[3],logg[4],logg[5],logg[12],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[12]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[12]))) {
	set_dist(d,natoms,logg[3],logg[12],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CZ3 and HZ2 (logg[4]&logg[13]) (transdihedral)*/
      pdih_lengths(logg[4],logg[5],logg[6],logg[13],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[13]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[13]))) {
	set_dist(d,natoms,logg[4],logg[13],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CG and CZ2 (logg[1]&logg[6]) (transdihedral)*/
      pdih_lengths(logg[1],logg[2],logg[7],logg[6],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[6]))) {
	set_dist(d,natoms,logg[1],logg[6],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CE2 and HD1 (logg[7]&logg[15]) (transdihedral)*/
      pdih_lengths(logg[7],logg[8],logg[9],logg[15],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[7]] != 0.0) || (weight[logg[15]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[7],logg[15]))) {
	set_dist(d,natoms,logg[7],logg[15],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CD2 and HD1 (logg[2]&logg[15]) (transdihedral)*/
      pdih_lengths(logg[2],logg[1],logg[9],logg[15],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-ring_margin)*blen;
      ub=(1.0+ring_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[15]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[15]))) {
	set_dist(d,natoms,logg[2],logg[15],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /* OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO*/

      /* new 1-5 distances added 981124 */

      /*SETDISTANCE for CZ3 and NE1 (1-5) */
      trp_15_type1(logg[8],logg[7],logg[2],logg[3],logg[4],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[4]))) {
	set_dist(d,natoms,logg[8],logg[4],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for CD1 and HE3 (1-5) */
      trp_15_type1(logg[9],logg[1],logg[2],logg[3],logg[10],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[10]))) {
	set_dist(d,natoms,logg[9],logg[10],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for CD1 and HZ2 (1-5) */
      trp_15_type1(logg[9],logg[8],logg[7],logg[6],logg[13],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[13]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[13]))) {
	set_dist(d,natoms,logg[9],logg[13],lb,ub,blen);
	ndist++;
	n15dist++;
      } 

      /*SETDISTANCE for CG and HZ3 (1-5) */
      trp_15_type2(logg[1],logg[2],logg[3],logg[4],logg[11],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[11]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[11]))) {
	set_dist(d,natoms,logg[1],logg[11],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for CH2 and CD1 (1-5) */
      trp_15_type2(logg[5],logg[6],logg[7],logg[8],logg[9],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[9]))) {
	set_dist(d,natoms,logg[5],logg[9],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for NE1 and HE3 (1-5) */
      trp_15_type2(logg[8],logg[7],logg[2],logg[3],logg[10],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[10]))) {
	set_dist(d,natoms,logg[8],logg[10],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for NE1 and HH2 (1-5) */
      trp_15_type2(logg[8],logg[7],logg[6],logg[5],logg[12],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[12]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[12]))) {
	set_dist(d,natoms,logg[8],logg[12],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HE3 and HH2 (1-5) */
      trp_15_type2(logg[10],logg[3],logg[4],logg[5],logg[12],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[12]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[12]))) {
	set_dist(d,natoms,logg[10],logg[12],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HZ3 and HZ2 (1-5) */
      trp_15_type2(logg[11],logg[4],logg[5],logg[6],logg[13],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[11]] != 0.0) || (weight[logg[13]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[11],logg[13]))) {
	set_dist(d,natoms,logg[11],logg[13],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*end of 1-5 distances added */


      /*SETDISTANCE for CB and CZ3 (1-5) */
      trp_15_type1(logg[4],logg[3],logg[2],logg[1],logg[0],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[0]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[0]))) {
	set_dist(d,natoms,logg[4],logg[0],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for CG and CH2 (1-5) */
      trp_15_type1(logg[1],logg[2],logg[7],logg[6],logg[5],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[5]))) {
	set_dist(d,natoms,logg[1],logg[5],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for CD and HH2 (1-5) */
      trp_15_type1(logg[12],logg[5],logg[6],logg[7],logg[2],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[12]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[12],logg[2]))) {
	set_dist(d,natoms,logg[12],logg[2],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HE3 and CZ2 (1-5) */
      trp_15_type1(logg[10],logg[3],logg[2],logg[7],logg[6],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[6]))) {
	set_dist(d,natoms,logg[10],logg[6],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HZ2 and CE3 (1-5) */
      trp_15_type1(logg[13],logg[6],logg[7],logg[2],logg[3],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[13]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[13],logg[3]))) {
	set_dist(d,natoms,logg[13],logg[3],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HZ3 and CE2 (1-5) */
      trp_15_type1(logg[11],logg[4],logg[3],logg[2],logg[7],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[11]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[11],logg[7]))) {
	set_dist(d,natoms,logg[11],logg[7],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for HE1 and CB (1-5) */
      trp_15_type2(logg[14],logg[8],logg[9],logg[1],logg[0],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[14]] != 0.0) || (weight[logg[0]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[14],logg[0]))) {
	set_dist(d,natoms,logg[14],logg[0],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for CB and CZ2 (1-5) */
      trp_15_type2(logg[0],logg[1],logg[2],logg[7],logg[6],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[6]))) {
	set_dist(d,natoms,logg[0],logg[6],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for CZ3 and CD1 (1-5) */
      trp_15_type2(logg[9],logg[1],logg[2],logg[3],logg[4],ilist,
		      iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[4]))) {
	set_dist(d,natoms,logg[9],logg[4],lb,ub,blen);
	ndist++;
	n15dist++;
      } 
      /*SETDISTANCE for CB and HE3 (1-5) */
      trp_15_type3(logg[0],logg[1],logg[2],logg[3],logg[10],ilist,
		   iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[10]))) {
	set_dist(d,natoms,logg[0],logg[10],lb,ub,blen);
	ndist++;
	n15dist++;
      }      




     /*SETDISTANCE for HD1 and CE3 (1-5) */
      trp_15_type2(logg[15],logg[9],logg[1],logg[2],logg[3],ilist,
		   iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[15]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[15],logg[3]))) {
	set_dist(d,natoms,logg[15],logg[3],lb,ub,blen);
	ndist++;
	n15dist++;
      }

      /*SETDISTANCE for HD1 and CZ2 (1-5) */
      trp_15_type2(logg[15],logg[9],logg[8],logg[7],logg[6],ilist,
		   iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[15]] != 0.0) || (weight[logg[6]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[15],logg[6]))) {
	set_dist(d,natoms,logg[15],logg[6],lb,ub,blen);
	ndist++;
	n15dist++;
      }
      /*SETDISTANCE for HE1 and HZ2 (1-5) */
      trp_15_type3(logg[14],logg[8],logg[7],logg[6],logg[13],ilist,
		   iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[14]] != 0.0) || (weight[logg[13]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[14],logg[13]))) {
	set_dist(d,natoms,logg[14],logg[13],lb,ub,blen);
	ndist++;
	n15dist++;
      }      
      /*SETDISTANCE for HE1 and CH2 (1-5) */
      trp_15_type1(logg[5],logg[6],logg[7],logg[8],logg[14],ilist,
		   iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[14]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[14]))) {
	set_dist(d,natoms,logg[5],logg[14],lb,ub,blen);
	ndist++;
	n15dist++;
      }      
      /*SETDISTANCE for HE1 and CE3 (1-5) */
      trp_15_type2(logg[14],logg[8],logg[7],logg[2],logg[3],ilist,
		   iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[14]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[14],logg[3]))) {
	set_dist(d,natoms,logg[14],logg[3],lb,ub,blen);
	ndist++;
	n15dist++;
      }      
      /*SETDISTANCE for HZ2 and CG (1-5) */
      trp_15_type2(logg[13],logg[6],logg[7],logg[2],logg[1],ilist,
		   iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if (((weight[logg[13]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[13],logg[1]))) {
	set_dist(d,natoms,logg[13],logg[1],lb,ub,blen);
	ndist++;
	n15dist++;
      }      
      /*}*/




      /*1111111111111111---------------66666666666666*/

      /* new 1-6 distances added 981124 */

      /*SETDISTANCE for HE3 and HE1 (1-6) */
      trp_16_type4(logg[10],logg[3],logg[2],logg[7],logg[8],logg[14],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[10]] != 0.0 || (weight[logg[14]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[10],logg[14]))) {
	set_dist(d,natoms,logg[10],logg[14],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      /* end new distances */

      /*SETDISTANCE for CB and CH2 (1-6) */
      trp_16_type1(logg[0],logg[1],logg[2],logg[3],logg[4],logg[5],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[0]] != 0.0 || (weight[logg[5]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[0],logg[5]))) {
	set_dist(d,natoms,logg[0],logg[5],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      /*SETDISTANCE for CZ3 and HE1 (1-6) */
      trp_16_type1(logg[4],logg[5],logg[6],logg[7],logg[8],logg[14],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[4]] != 0.0 || (weight[logg[14]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[4],logg[14]))) {
	set_dist(d,natoms,logg[4],logg[14],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      /*SETDISTANCE for CG  and HH2 (1-6) */
      trp_16_type2(logg[12],logg[5],logg[6],logg[7],logg[2],logg[1],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[12]] != 0.0 || (weight[logg[1]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[12],logg[1]))) {
	set_dist(d,natoms,logg[12],logg[1],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      /*SETDISTANCE for HZ3 and NE1 (1-6) */
      trp_16_type2(logg[11],logg[4],logg[3],logg[2],logg[7],logg[8],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[11]] != 0.0 || (weight[logg[8]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[11],logg[8]))) {
	set_dist(d,natoms,logg[11],logg[8],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      /*SETDISTANCE for HZ2 and HE3 (1-6) */
      trp_16_type2(logg[13],logg[6],logg[5],logg[4],logg[3],logg[10],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[13]] != 0.0 || (weight[logg[10]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[13],logg[10]))) {
	set_dist(d,natoms,logg[13],logg[10],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      /*SETDISTANCE for CB  and HZ3 (1-6) */
      trp_16_type3(logg[0],logg[1],logg[2],logg[3],logg[4],logg[11],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[0]] != 0.0 || (weight[logg[11]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[0],logg[11]))) {
	set_dist(d,natoms,logg[0],logg[11],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      /*SETDISTANCE for HE1 and HH2 (1-6) */
      trp_16_type3(logg[14],logg[8],logg[7],logg[6],logg[5],logg[12],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[14]] != 0.0 || (weight[logg[12]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[14],logg[12]))) {
	set_dist(d,natoms,logg[14],logg[12],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      /*SETDISTANCE for HZ2 and HD1 (1-6) */
      trp_16_type3(logg[13],logg[6],logg[7],logg[8],logg[9],logg[15],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[13]] != 0.0 || (weight[logg[15]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[13],logg[15]))) {
	set_dist(d,natoms,logg[13],logg[15],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      /*SETDISTANCE for HE3 and HD1 (1-6) */
      trp_16_type3(logg[10],logg[3],logg[2],logg[1],logg[9],logg[15],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[10]] != 0.0 || (weight[logg[15]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[10],logg[15]))) {
	set_dist(d,natoms,logg[10],logg[15],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      /*SETDISTANCE for CB and HZ2 (1-6) */
      trp_16_type4(logg[0],logg[1],logg[2],logg[7],logg[6],logg[13],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[0]] != 0.0 || (weight[logg[13]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[0],logg[13]))) {
	set_dist(d,natoms,logg[0],logg[13],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      /*SETDISTANCE for CZ3 and HD1 (1-6) */
      trp_16_type4(logg[4],logg[3],logg[2],logg[1],logg[9],logg[15],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[4]] != 0.0 || (weight[logg[15]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[4],logg[15]))) {
	set_dist(d,natoms,logg[4],logg[15],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      /*SETDISTANCE for CH2 and HD1 (1-6) */
      trp_16_type4(logg[5],logg[6],logg[7],logg[8],logg[9],logg[15],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[5]] != 0.0 || (weight[logg[15]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[5],logg[15]))) {
	set_dist(d,natoms,logg[5],logg[15],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      /*SETDISTANCE for HZ3 and CD1 (1-6) */
      trp_16_type4(logg[11],logg[4],logg[3],logg[2],logg[1],logg[9],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[11]] != 0.0 || (weight[logg[9]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[11],logg[9]))) {
	set_dist(d,natoms,logg[11],logg[9],lb,ub,blen);
	ndist++;
	n16dist++;
      }
      /*SETDISTANCE for HH2 and CD1 (1-6) */
      trp_16_type4(logg[12],logg[5],logg[6],logg[7],logg[8],logg[9],
		   ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[12]] != 0.0 || (weight[logg[9]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[12],logg[9]))) {
	set_dist(d,natoms,logg[12],logg[9],lb,ub,blen);
	ndist++;
	n16dist++;
      }

      /* 111111111111---------7777777777 */


      /* new 1-7 distances added 981124 */

      /*SETDISTANCE for HH2 and HD1 (1-6) */
      trp_17_type2(logg[12],logg[5],logg[6],logg[7],logg[8],logg[9],
		   logg[15],ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[12]] != 0.0 || (weight[logg[15]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[12],logg[15]))) {
	set_dist(d,natoms,logg[12],logg[15],lb,ub,blen);
	ndist++;
	n17dist++;
      }
      /* end of new distances */

      /*SETDISTANCE for HZ3 and HE1 (1-6) */
      trp_17_type1(logg[11],logg[4],logg[3],logg[2],logg[7],logg[8],
		   logg[14],ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[11]] != 0.0 || (weight[logg[14]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[11],logg[14]))) {
	set_dist(d,natoms,logg[11],logg[14],lb,ub,blen);
	ndist++;
	n17dist++;
      }
      /*SETDISTANCE for HH2 and CB (1-6) */
      trp_17_type1(logg[12],logg[5],logg[6],logg[7],logg[2],logg[1],
		   logg[0],ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[12]] != 0.0 || (weight[logg[0]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[12],logg[0]))) {
	set_dist(d,natoms,logg[12],logg[0],lb,ub,blen);
	ndist++;
	n17dist++;
      }
      /*SETDISTANCE for HZ3 and HD1 (1-6) */
      trp_17_type2(logg[11],logg[4],logg[3],logg[2],logg[1],logg[9],
		   logg[15],ilist,iparams,&blen,atoms);
      lb=(1.0-(ring_margin*0.5))*blen;
      ub=(1.0+(ring_margin*0.5))*blen;
      if ((weight[logg[11]] != 0.0 || (weight[logg[15]] != 0.0)) && 
	  (!dist_set(d,natoms,logg[11],logg[15]))) {
	set_dist(d,natoms,logg[11],logg[15],lb,ub,blen);
	ndist++;
	n17dist++;
      }

      if (bVir) {
      /* VIRTUAL DISTANCES */
	nVdist += set_virtual (logg,16,ring_margin,d,natoms);
	ndist += nVdist;
	/* Try adding a virtual atom above CD1 (logg[9]) */
	if (FALSE) {
	for (q=0 ; q<16 ; q++ ) {
	  if (q == 9) {
	    blen = 8.0;
	    lb   = (1.0-ring_margin)*blen;
	    ub   = (1.0+ring_margin)*blen;
	  }
	  else {
	    blen = d_len(d,natoms,logg[q],logg[9]);
	    blen = sqrt(blen*blen+64.0);
	    lb   = (1.0-ring_margin)*blen;
          ub   = (1.0+ring_margin)*blen;
	  }
        /* set distance to virtual atom */
	  set_dist(d,natoms,logg[q],logg[16],lb,ub,blen);
	  ndist++;
	  nVdist++;       
	}
	}
      }
    }
    logg[17]=0;
    }
  }

  fprintf(log,"There are %d new tryptophan distances\n",ndist);
  if (ndist > 0 ) {
  fprintf(log,"(%d 1-2, %d 1-3, %d 1-4, %d 1-5, %d 1-6 %d 1-7, %d virtual)\n",n12dist,
	  n13dist,n14dist,n15dist,n16dist,n17dist,nVdist);
  }
} 
/**********************************************************
 *
 *     V A L I N E
 *
 **********************************************************/

void val (FILE *log,t_dist *d,t_idef *idef,t_atoms *atoms,real weight[],
	  real val_margin,t_ilist ilist[],t_iparams iparams[])
{
  int natoms,ndist,i,j,logg[12],residnr,oldresidnr,n14dist,n15dist;
  real blen,lb,ub,angle;
  real pi = M_PI;

  natoms= atoms->nr;
  ndist   = 0;
  n14dist = 0;
  n15dist = 0;
  residnr = -1;
  oldresidnr=-1;

  for (i=0; (i<natoms); i++) {
    if ( strcmp((*atoms->resname[atoms->atom[i].resnr]),"VAL") == 0) {
      oldresidnr=residnr;
      residnr=atoms->atom[i].resnr;
      logg[11]=0;
      if ( oldresidnr == residnr ) {
	continue;
      }
      j=i;
      while ((atoms->atom[j].resnr) == residnr ) {
	if ( strcmp((*atoms->atomname[j]),"CA") == 0) {
	  logg[11]++;
	  logg[0]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CB") == 0) {
	  logg[11]++;
	  logg[1]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"HB") == 0) {
	  logg[11]++;
	  logg[2]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"CG1") == 0) {
	  logg[11]++;
	  logg[3]=j;
	}
	else if ( strcmp((*atoms->atomname[j]),"HG11") == 0) {
	  logg[11]++;
	  logg[4]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HG12") == 0) {
	  logg[11]++;
	  logg[5]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HG13") == 0) {
	  logg[11]++;
	  logg[6]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"CG2") == 0) {
	  logg[11]++;
	  logg[7]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HG21") == 0) {
	  logg[11]++;
	  logg[8]=j;
	}
        else if ( strcmp((*atoms->atomname[j]),"HG22") == 0) {
	  logg[11]++;
	  logg[9]=j;
	}	
        else if ( strcmp((*atoms->atomname[j]),"HG23") == 0) {
	  logg[11]++;
	  logg[10]=j;
	}	
	if ( logg[11] == 11) {
	  break;
	}
	j++;
      }
    if ( logg[11] == 11 ) {
      fprintf(log,"logg (val) = %d %d %d %d %d %d %d %d %d %d %d\n",
	      logg[0],logg[1],logg[2],logg[3],logg[4],logg[5],logg[6],
	      logg[7],logg[8],logg[9],logg[10]);
      /*SETDISTANCE for HG11 and CA */
      gauche(logg[4],logg[3],logg[1],logg[0],ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[0]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[0]))) {
	set_dist(d,natoms,logg[4],logg[0],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG11 and CG2 */
      gauche(logg[4],logg[3],logg[1],logg[7],ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[7]))) {
	set_dist(d,natoms,logg[4],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG11 and HB */
      pdih_lengths(logg[4],logg[3],logg[1],logg[2],ilist,iparams,
		   &lb,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[2]))) {
	set_dist(d,natoms,logg[4],logg[2],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG12 and CG2 */
      gauche(logg[5],logg[3],logg[1],logg[7],ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[7]))) {
	set_dist(d,natoms,logg[5],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG12 and HB */
      gauche(logg[5],logg[3],logg[1],logg[2],ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[2]))) {
	set_dist(d,natoms,logg[5],logg[2],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG12 and CA */
      pdih_lengths(logg[5],logg[3],logg[1],logg[0],ilist,iparams,
		   &lb,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[0]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[0]))) {
	set_dist(d,natoms,logg[5],logg[0],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG13 and HB */
      gauche(logg[6],logg[3],logg[1],logg[2],ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[2]))) {
	set_dist(d,natoms,logg[6],logg[2],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG13 and CA */
      gauche(logg[6],logg[3],logg[1],logg[0],ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[0]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[0]))) {
	set_dist(d,natoms,logg[6],logg[0],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG13 and CG2 */
      pdih_lengths(logg[6],logg[3],logg[1],logg[7],ilist,iparams,
		   &lb,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[7]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[7]))) {
	set_dist(d,natoms,logg[6],logg[7],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG21 and HB */
      gauche(logg[8],logg[7],logg[1],logg[2],ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[2]))) {
	set_dist(d,natoms,logg[8],logg[2],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG21 and CG1 */
      gauche(logg[8],logg[7],logg[1],logg[3],ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[3]))) {
	set_dist(d,natoms,logg[8],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG21 and CA */
      pdih_lengths(logg[8],logg[7],logg[1],logg[0],ilist,iparams,
		   &lb,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[8]] != 0.0) || (weight[logg[0]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[8],logg[0]))) {
	set_dist(d,natoms,logg[8],logg[0],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG22 and CG1 */
      gauche(logg[9],logg[7],logg[1],logg[3],ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[3]))) {
	set_dist(d,natoms,logg[9],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG22 and CA */
      gauche(logg[9],logg[7],logg[1],logg[0],ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[0]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[0]))) {
	set_dist(d,natoms,logg[9],logg[0],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG22 and HB */
      pdih_lengths(logg[9],logg[7],logg[1],logg[2],ilist,iparams,
		   &lb,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[9]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[9],logg[2]))) {
	set_dist(d,natoms,logg[9],logg[2],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG23 and CA */
      gauche(logg[10],logg[7],logg[1],logg[0],ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[0]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[0]))) {
	set_dist(d,natoms,logg[10],logg[0],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG23 and HB */
      gauche(logg[10],logg[7],logg[1],logg[2],ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[2]))) {
	set_dist(d,natoms,logg[10],logg[2],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for HG23 and CG1 */
      pdih_lengths(logg[10],logg[7],logg[1],logg[3],ilist,iparams,
		   &lb,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[10]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[10],logg[3]))) {
	set_dist(d,natoms,logg[10],logg[3],lb,ub,blen);
	ndist++;
	n14dist++;
      }

      /* 111111111-----------5555555555 */
      /*SETDISTANCE for HG11 and HG21 */
      gauche15(logg[4],logg[3],logg[1],logg[7],logg[8],pi,pi+(pi/3.0),
	       pi+(pi/3.0),ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[8]))) {
	set_dist(d,natoms,logg[4],logg[8],lb,ub,blen);
	ndist++;
	n15dist++;
      }
      /*SETDISTANCE for HG11 and HG22 */
      gauche15(logg[4],logg[3],logg[1],logg[7],logg[9],pi,pi+pi/3.0,
	       pi-(pi/3.0),ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[9]))) {
	set_dist(d,natoms,logg[4],logg[9],lb,ub,blen);
	ndist++;
	n15dist++;
      }
      /*SETDISTANCE for HG11 and HG23 */
      gauche15(logg[4],logg[3],logg[1],logg[7],logg[10],pi,pi+(pi/3.0),0,
	       ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[10]))) {
	set_dist(d,natoms,logg[4],logg[10],lb,ub,blen);
	ndist++;
	n15dist++;
      }


      /*SETDISTANCE for HG12 and HG21 */
      gauche15(logg[5],logg[3],logg[1],logg[7],logg[8],pi,pi-pi/3.0,
	       pi+pi/3.0,ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[8]))) {
	set_dist(d,natoms,logg[5],logg[8],lb,ub,blen);
	ndist++;
	n15dist++;
      }
      /*SETDISTANCE for HG12 and HG22 */
      gauche15(logg[5],logg[3],logg[1],logg[7],logg[9],pi,pi-pi/3.0,
	       pi-pi/3.0,ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[9]))) {
	set_dist(d,natoms,logg[5],logg[9],lb,ub,blen);
	ndist++;
	n15dist++;
      }
      /*SETDISTANCE for HG12 and HG23 */
      gauche15(logg[5],logg[3],logg[1],logg[7],logg[10],pi,pi-pi/3.0,0,ilist,
	       iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[5]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[5],logg[10]))) {
	set_dist(d,natoms,logg[5],logg[10],lb,ub,blen);
	ndist++;
	n15dist++;
      }

      /*SETDISTANCE for HG13 and HG21 */
      gauche15(logg[6],logg[3],logg[1],logg[7],logg[8],pi,0,pi+pi/3.0,
	       ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[8]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[8]))) {
	set_dist(d,natoms,logg[6],logg[8],lb,ub,blen);
	ndist++;
	n15dist++;
      }
      /*SETDISTANCE for HG13 and HG22 */
      gauche15(logg[6],logg[3],logg[1],logg[7],logg[9],pi,0,pi-pi/3.0,
	       ilist,iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[9]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[9]))) {
	set_dist(d,natoms,logg[6],logg[9],lb,ub,blen);
	ndist++;
	n15dist++;
      }
      /*SETDISTANCE for HG13 and HG23 */
      gauche15(logg[6],logg[3],logg[1],logg[7],logg[10],pi,0,0,ilist,
	       iparams,&blen,atoms);
      lb=(1.0-val_margin)*blen;
      ub=(1.0+val_margin)*blen;
      if (((weight[logg[6]] != 0.0) || (weight[logg[10]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[6],logg[10]))) {
	set_dist(d,natoms,logg[6],logg[10],lb,ub,blen);
	ndist++;
	n15dist++;
      }


    }
    logg[11]=0;
    }
  }
  fprintf(log,"There are %d distances to keep valine gauche\n",ndist);
  if (ndist > 0 ) {
    fprintf(log,"(%d 1-4, %d 1-5)\n",n14dist,n15dist);
  }
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
   
  int natoms,ndist,odist,i,j,q,logg[8],n12dist,n13dist,n14dist,nVdist,
    pronr,oldpronr;
  real blen,lb,ub,angle;
  
  natoms= atoms->nr;
  ndist = 0;
  odist = 0;
  n12dist = 0;
  n13dist = 0;
  n14dist = 0;
  nVdist = 0;
  pronr   =-1;
  oldpronr=-1;
 
  for (i=0; (i<natoms); i++) {
    logg[7]=0;
    oldpronr=pronr;
    if ( strcmp((*atoms->atomname[i]),"CA") == 0) {
      logg[7]=1;
      logg[0]=i;
      /* Here comes an ugly initialisation. Put it in a separate function.*/
      logg[1]=-1;
      logg[2]=-1;
      logg[3]=-1;
      logg[4]=-1;
      logg[5]=-1;
      logg[6]=-1;
      /* end of ugly initialisation */
      /*for (j=(i+1); (j<natoms); j++) {*/
      j=i+1;
      while ( ((atoms->atom[j].resnr) <= ((atoms->atom[i].resnr)+1)) &&
	      (j<natoms) ) {
	if ( strcmp((*atoms->atomname[j]),"C") == 0) {
	  if (logg[1] == -1) {
	    logg[7]++;
	    logg[1]=j;
	  }
	}
	if ( strcmp((*atoms->atomname[j]),"O") == 0 ) {
	  if (logg[2] == -1) {
	    logg[7]++;
	    logg[2]=j;
	  }
	}
	if ( strcmp((*atoms->atomname[j]),"N") == 0 ) {
	  if (logg[3] == -1) {
	    logg[7]++;
	    logg[3]=j;
	  }
	}
	if (( strcmp((*atoms->resname[atoms->atom[j].resnr]),"PRO") == 0) &&
	    ( strcmp((*atoms->atomname[j]),"CD") == 0)) {
	  pronr=atoms->atom[j].resnr;
	  if (oldpronr != pronr) {
	      if (logg[4] == -1) {
		logg[7]++;
		logg[4]=j;
	      }
	  }
	}
	else if ( strcmp((*atoms->atomname[j]),"H") == 0 ) {
	  if (logg[4] == -1) {
	    logg[7]++;
	    logg[4]=j;
	  }
	}
	if ( strcmp((*atoms->atomname[j]),"CA") == 0) {
	  if (logg[5] == -1) {
	    logg[7]++;
	    logg[5]=j;
	  }
	}
	if ( (strcmp((*atoms->atomname[j]),"VP") == 0) && bVir) {
	  if (logg[6] == -1) {
	    logg[7]++;
	    logg[6]=j;
	  }
	}
	if ( ((logg[7] == 6) && !bVir) || ((logg[7] == 7) && bVir) ) {
	  break;
	}
	j++;
      } /* end of while-loop */

    }
    if ( ((logg[7] == 6) && !bVir) || ((logg[7] == 7) && bVir) ) {
      if ((logg[7] == 6) && !bVir) {
	
	fprintf(log,"logg (pep) = %d %d %d %d %d %d\n",logg[0],logg[1],
		logg[2],logg[3],logg[4],logg[5]);
      }
      else if ((logg[7] == 7) && bVir) {
	fprintf(log,"logg (pep) = %d %d %d %d %d %d %d\n",logg[0],logg[1],
		logg[2],logg[3],logg[4],logg[5],logg[6]);
      }
      /*SETDISTANCE for O and H (logg[2]&logg[4]) (transdihedral)*/
      pdih_lengths(logg[2],logg[1],logg[3],logg[4],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-pep_margin)*blen;
      ub=(1.0+pep_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  
	  (!dist_set(d,natoms,logg[2],logg[4]))) {
	set_dist(d,natoms,logg[2],logg[4],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for O and CA+ (logg[2]&logg[5]) (cisdihedral)*/
      pdih_lengths(logg[2],logg[1],logg[3],logg[5],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-pep_margin)*blen;
      ub=(1.0+pep_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[5]))) {
	set_dist(d,natoms,logg[2],logg[5],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CA- and CA+ (logg[0]&logg[5]) (transdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[3],logg[5],ilist,iparams,&lb,&blen,
		   atoms);
      lb=(1.0-pep_margin)*blen;
      ub=(1.0+pep_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[5]))) {
	set_dist(d,natoms,logg[0],logg[5],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CA- and H (logg[0]&logg[4]) (cisdihedral)*/
      pdih_lengths(logg[0],logg[1],logg[3],logg[4],ilist,iparams,&blen,&ub,
		   atoms);
      lb=(1.0-pep_margin)*blen;
      ub=(1.0+pep_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[4]))) {
	set_dist(d,natoms,logg[0],logg[4],lb,ub,blen);
	ndist++;
	n14dist++;
      }
      /*SETDISTANCE for CA- and C (1-2) */
      blen=lookup_bondlength(logg[0],logg[1],ilist,iparams,TRUE,atoms);
      lb=(1.0-pep_margin)*blen;
      ub=(1.0+pep_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[1]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[1]))) {
	set_dist(d,natoms,logg[0],logg[1],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for N and CA+ (1-2) */
      blen=lookup_bondlength(logg[3],logg[5],ilist,iparams,TRUE,atoms);
      lb=(1.0-pep_margin)*blen;
      ub=(1.0+pep_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[5]))) {
	set_dist(d,natoms,logg[3],logg[5],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for C and O (1-2) */
      blen=lookup_bondlength(logg[1],logg[2],ilist,iparams,TRUE,atoms);
      lb=(1.0-pep_margin)*blen;
      ub=(1.0+pep_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[2]))) {
	set_dist(d,natoms,logg[1],logg[2],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for C and N (1-2) */
      blen=lookup_bondlength(logg[1],logg[3],ilist,iparams,TRUE,atoms);
      lb=(1.0-pep_margin)*blen;
      ub=(1.0+pep_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[3]))) {
	set_dist(d,natoms,logg[1],logg[3],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for N and H (1-2) */
      blen=lookup_bondlength(logg[3],logg[4],ilist,iparams,TRUE,atoms);
      lb=(1.0-pep_margin)*blen;
      ub=(1.0+pep_margin)*blen;
      if (((weight[logg[3]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[3],logg[4]))) {
	set_dist(d,natoms,logg[3],logg[4],lb,ub,blen);
	ndist++;
	n12dist++;
      }
      /*SETDISTANCE for O and N (1-3) */
      angle=lookup_angle(logg[2],logg[1],logg[3],ilist,iparams,atoms);
      blen=angle_length(logg[2],logg[1],logg[3],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-pep_margin)*blen;
      ub=(1.0+pep_margin)*blen;
      if (((weight[logg[2]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[2],logg[3]))) {
	set_dist(d,natoms,logg[2],logg[3],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for C and H (1-3) */
      angle=lookup_angle(logg[1],logg[3],logg[4],ilist,iparams,atoms);
      /*      fprintf(stderr,"Hej!pept.c(C-H)1-3 has angle:%g\n",RAD2DEG*angle);*/
      blen=angle_length(logg[1],logg[3],logg[4],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-pep_margin)*blen;
      ub=(1.0+pep_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[4]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[4]))) {
	set_dist(d,natoms,logg[1],logg[4],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CA- and O (1-3) */
      angle=lookup_angle(logg[0],logg[1],logg[2],ilist,iparams,atoms);
      blen=angle_length(logg[0],logg[1],logg[2],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-pep_margin)*blen;
      ub=(1.0+pep_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[2]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[2]))) {
	set_dist(d,natoms,logg[0],logg[2],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for H and CA+ (1-3) */
      angle=lookup_angle(logg[4],logg[3],logg[5],ilist,iparams,atoms);
      /*     fprintf(stderr,"Hej!pept.c(H-CA+)1-3 has angle:%g\n",RAD2DEG*angle);*/
      blen=angle_length(logg[4],logg[3],logg[5],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-pep_margin)*blen;
      ub=(1.0+pep_margin)*blen;
      if (((weight[logg[4]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[4],logg[5]))) {
	set_dist(d,natoms,logg[4],logg[5],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for CA- and N (1-3) */
      angle=lookup_angle(logg[0],logg[1],logg[3],ilist,iparams,atoms);
      /*fprintf(stderr,"CA- N (1-3) angle:%g\n",RAD2DEG*angle);*/
      blen=angle_length(logg[0],logg[1],logg[3],RAD2DEG*angle,ilist,iparams,
			atoms);
      /*fprintf(stderr,"CA- N (1-3) dist:%g\n",blen);*/
      lb=(1.0-pep_margin)*blen;
      ub=(1.0+pep_margin)*blen;
      if (((weight[logg[0]] != 0.0) || (weight[logg[3]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[0],logg[3]))) {
	set_dist(d,natoms,logg[0],logg[3],lb,ub,blen);
	ndist++;
	n13dist++;
      }
      /*SETDISTANCE for C and CA+ (1-3) */
      angle=lookup_angle(logg[1],logg[3],logg[5],ilist,iparams,atoms);
      /*fprintf(stderr,"Hej!pept.c(C-CA+)1-3 has angle:%g\n",RAD2DEG*angle); */
      blen=angle_length(logg[1],logg[3],logg[5],RAD2DEG*angle,ilist,iparams,
			atoms);
      lb=(1.0-pep_margin)*blen;
      ub=(1.0+pep_margin)*blen;
      if (((weight[logg[1]] != 0.0) || (weight[logg[5]] != 0.0)) &&
	  (!dist_set(d,natoms,logg[1],logg[5]))) {
	set_dist(d,natoms,logg[1],logg[5],lb,ub,blen);
	ndist++;
	n13dist++;
      }

      if (bVir) {
      /* VIRTUAL DISTANCES */
	nVdist += set_virtual (logg,6,pep_margin,d,natoms);
	ndist += nVdist;
	if (FALSE) {
	  for (q=0 ; q<6 ; q++ ) {
	    if (q == 3) {
	      blen = 8.0;
	      lb   = (1.0-pep_margin)*blen;
	      ub   = (1.0+pep_margin)*blen;
	    }
	    else {
	      blen = d_len(d,natoms,logg[q],logg[3]);
	      blen = sqrt(blen*blen+64.0);
	      lb   = (1.0-pep_margin)*blen;
	      ub   = (1.0+pep_margin)*blen;
	    }
	    set_dist(d,natoms,logg[q],logg[6],lb,ub,blen);
	    ndist++;
	    nVdist++;
	  }
	}
      }
    }
  }
  fprintf(stderr,"There are %d new peptide distances\n",ndist);
  if (ndist > 0) {
    fprintf(stderr,"(%d 1-2, %d 1-3, %d 1-4, %d virtual)\n",
	    n12dist,n13dist,n14dist,nVdist);
  }
}
