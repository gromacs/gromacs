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
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_vcm_c = "$Id$";

#include "macros.h"
#include "vcm.h"
#include "vec.h"
#include "smalloc.h"
#include "do_fit.h"
 
void calc_vcm(FILE *log,int homenr,int start,real mass[],rvec v[],rvec vcm)
{
  int    i;
  real   m0;
  real   x,y,z;
  
  /* Calculate */
  x=y=z=0;
  for(i=start; (i<start+homenr); i++) {
    m0=mass[i];
    x+=m0*v[i][XX];
    y+=m0*v[i][YY];
    z+=m0*v[i][ZZ];
  }
  vcm[XX]=x;
  vcm[YY]=y;
  vcm[ZZ]=z;
}

void do_stopcm(FILE *log,int homenr,int start,rvec v[],rvec mvcm,
	       real tm,real invmass[])
{
  int  i;
  rvec vcm;
  
  vcm[XX]=mvcm[XX]/tm;
  vcm[YY]=mvcm[YY]/tm;
  vcm[ZZ]=mvcm[ZZ]/tm;
  
  for(i=start; (i<start+homenr); i++) 
    if (invmass[i] != 0)
      rvec_dec(v[i],vcm);
}

void check_cm(FILE *log,rvec mvcm,real tm)
{
  int    m;
  rvec   vcm;
  real   ekcm=0,max_vcm=0;

  for(m=0; (m<DIM); m++) {
    vcm[m]=mvcm[m]/tm;
    max_vcm=max(max_vcm,fabs(vcm[m]));
    ekcm+=vcm[m]*vcm[m];
  }
  if (max_vcm > 0.1) {
    ekcm*=0.5*tm;
    fprintf(log,"Large VCM: (%12.5f,  %12.5f,  %12.5f), ekin-cm: %12.5f\n",
	    vcm[XX],vcm[YY],vcm[ZZ],ekcm);
  }
}

void do_stoprot(FILE *log, int natoms, rvec box, rvec x[], real mass[])
{
  static atom_id *index=NULL;
  static rvec *old_x=NULL;
  rvec   half_box;
  int    i;
  
  for (i=0; (i<DIM); i++)
    half_box[i]=box[i]*0.5;
  if (!index) {
    snew(index,natoms);
    for (i=0; (i<natoms); i++)
      index[i]=i;
  }
  
  /* first center on origin (do_fit does not do this!) */
  reset_x(natoms,index,natoms,index,x,mass);
  
  if (old_x) {
    /* do least squares fit to previous structure */
    do_fit(natoms,mass,old_x,x);
  } else {
    /* no previous structure, make room to copy this one */
    snew(old_x, natoms);
  }
  /* save structure for next fit and move structure back to center of box */
  for (i=0; (i<natoms); i++) {
    copy_rvec(x[i],old_x[i]);
    rvec_inc(x[i],half_box);
  }
}

t_vcm *init_vcm(FILE *fp,t_topology *top,t_mdatoms *md,int nstcomm)
{
  t_vcm *vcm;
  int i,g;
  
  snew(vcm,1);
  
  vcm->nr = top->atoms.grps[egcVCM].nr;
  snew(vcm->group_mvcm,vcm->nr);
  snew(vcm->group_mass,vcm->nr);
  snew(vcm->group_name,vcm->nr);
  vcm->group_id = md->cVCM;
  
  /* Not parallel... */
  if (nstcomm != 0)
    for(i=0; (i<md->nr); i++) {
      g = vcm->group_id[i];
      vcm->group_mass[g] += md->massT[i];
    }
  
  /* Copy pointer to group names and print it. */
  fprintf(fp,"We have the following groups for center of mass motion removal:\n");
  for(g=0; (g<vcm->nr); g++) {
    vcm->group_name[g] = *top->atoms.grpname[top->atoms.grps[egcVCM].nm_ind[g]];
    fprintf(fp,"%3d:  %s, initial mass: %g\n",
	    g,vcm->group_name[g],vcm->group_mass[g]);
  }
  
  return vcm;
}

/* Center of mass code for groups */
void calc_vcm_grp(FILE *log,int homenr,int start,real mass[],rvec v[],
		  t_vcm *vcm)
{
  int    i,g;
  real   m0;
  
  /* Reset */
  for(g=0; (g<vcm->nr); g++) {
    clear_rvec(vcm->group_mvcm[g]);
    vcm->group_mass[g] = 0;
  }
    
  /* Calculate */
  for(i=start; (i<start+homenr); i++) {
    m0 = mass[i];
    g  = vcm->group_id[i];
    
    vcm->group_mass[g]     += m0;
    vcm->group_mvcm[g][XX] += m0*v[i][XX];
    vcm->group_mvcm[g][YY] += m0*v[i][YY];
    vcm->group_mvcm[g][ZZ] += m0*v[i][ZZ];
  }
}

void do_stopcm_grp(FILE *log,int homenr,int start,rvec v[],
		   t_vcm *vcm,real invmass[])
{
  int  i,g;
  rvec *my_vcm;
  real tm;
  
  snew(my_vcm,vcm->nr);
  for(g=0; (g<vcm->nr); g++) {
    tm = vcm->group_mass[g];
    if (tm != 0) {
      my_vcm[g][XX] = vcm->group_mvcm[g][XX]/tm;
      my_vcm[g][YY] = vcm->group_mvcm[g][YY]/tm;
      my_vcm[g][ZZ] = vcm->group_mvcm[g][ZZ]/tm;
    }
    /* Else it's zero anyway */
  }
  for(i=start; (i<start+homenr); i++) {
    g = vcm->group_id[i];
    if (invmass[i] != 0)
      rvec_dec(v[i],my_vcm[g]);
  }
  
  sfree(my_vcm);
}

void check_cm_grp(FILE *log,t_vcm *vcm)
{
  int  m,g;
  rvec my_vcm;
  real ekcm,max_vcm;

  for(g=0; (g<vcm->nr); g++) {
    ekcm    = 0;
    max_vcm = 0;
    if (vcm->group_mass[g] != 0) {
      for(m=0; (m<DIM); m++) {
	my_vcm[m] = vcm->group_mvcm[g][m]/vcm->group_mass[g];
	max_vcm = max(max_vcm,fabs(my_vcm[m]));
	ekcm += my_vcm[m]*my_vcm[m];
      }
      if (max_vcm > 0.1) {
	ekcm*=0.5*vcm->group_mass[g];
	fprintf(log,"Large VCM(group %s): (%12.5f,  %12.5f,  %12.5f), ekin-cm: %12.5f\n",
		vcm->group_name[g],my_vcm[XX],my_vcm[YY],my_vcm[ZZ],ekcm);
      }
    }
  }
}


