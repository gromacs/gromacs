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
 * GROwing Monsters And Cloning Shrimps
 */
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "macros.h"
#include "vcm.h"
#include "vec.h"
#include "smalloc.h"
#include "names.h"
#include "txtdump.h"
#include "network.h"
 
t_vcm *init_vcm(FILE *fp,t_topology *top,t_commrec *cr,t_mdatoms *md,
		int start,int homenr,int nstcomm,int comm_mode)
{
  t_vcm *vcm;
  real  *mass;
  int   i,g;
  
  snew(vcm,1);
  
  vcm->mode = (nstcomm > 0) ? comm_mode : ecmNO;
  
  if (vcm->mode != ecmNO) {
    vcm->nr = top->atoms.grps[egcVCM].nr;
    if (vcm->mode == ecmANGULAR) {
      snew(vcm->group_j,vcm->nr);
      snew(vcm->group_x,vcm->nr);
      snew(vcm->group_i,vcm->nr);
      snew(vcm->group_w,vcm->nr);
    }
    snew(vcm->group_p,vcm->nr);
    snew(vcm->group_v,vcm->nr);
    snew(vcm->group_mass,vcm->nr);
    snew(vcm->group_name,vcm->nr);
    vcm->group_id = md->cVCM;
    for(i=start; (i<start+homenr); i++) {
      g = vcm->group_id[i];
      vcm->group_mass[g] += md->massT[i];
    }
    snew(mass,vcm->nr);
    for(g=0; (g<vcm->nr); g++)
      mass[g] = vcm->group_mass[g];
    if(PAR(cr))
      gmx_sum(vcm->nr,mass,cr);

    /* Copy pointer to group names and print it. */
    fprintf(fp,"Center of mass motion removal mode is %s\n",ECOM(vcm->mode));
    fprintf(fp,"We have the following groups for center of"
	    " mass motion removal:\n");
    for(g=0; (g<vcm->nr); g++) {
      vcm->group_name[g] = 
	*top->atoms.grpname[top->atoms.grps[egcVCM].nm_ind[g]];
      fprintf(fp,"%3d:  %s, initial mass: %g\n",
	      g,vcm->group_name[g],mass[g]);
    }
    sfree(mass);
  }
  return vcm;
}

static void update_tensor(rvec x,real m0,tensor I)
{
  real xy,xz,yz;
  
  /* Compute inertia tensor contribution due to this atom */
  xy         = x[XX]*x[YY]*m0;
  xz         = x[XX]*x[ZZ]*m0;
  yz         = x[YY]*x[ZZ]*m0;
  I[XX][XX] += x[XX]*x[XX]*m0;
  I[YY][YY] += x[YY]*x[YY]*m0;
  I[ZZ][ZZ] += x[ZZ]*x[ZZ]*m0;
  I[XX][YY] += xy;
  I[YY][XX] += xy;
  I[XX][ZZ] += xz;
  I[ZZ][XX] += xz;
  I[YY][ZZ] += yz;
  I[ZZ][YY] += yz;
}

/* Center of mass code for groups */
void calc_vcm_grp(FILE *fp,int start,int homenr,real mass[],
		  rvec x[],rvec v[],t_vcm *vcm)
{
  int    i,g,m;
  real   m0,xx,xy,xz,yy,yz,zz;
  rvec   j0;
  
  if (vcm->mode != ecmNO) {
    for(g=0; (g<vcm->nr); g++) {
      /* Reset linear momentum */
      vcm->group_mass[g] = 0;
      clear_rvec(vcm->group_p[g]);
      
      if (vcm->mode == ecmANGULAR) {
	/* Reset anular momentum */
	clear_rvec(vcm->group_j[g]);
	clear_rvec(vcm->group_x[g]);
	clear_rvec(vcm->group_w[g]);
	clear_mat(vcm->group_i[g]);
      }
    }
    
    for(i=start; (i<start+homenr); i++) {
      m0 = mass[i];
      g  = vcm->group_id[i];
      
      /* Calculate linear momentum */
      vcm->group_mass[g]  += m0;
      for(m=0; (m<DIM);m++)
	vcm->group_p[g][m] += m0*v[i][m];

      if (vcm->mode == ecmANGULAR) {
	/* Calculate angular momentum */
	oprod(x[i],v[i],j0);
	
	for(m=0; (m<DIM); m++) {
	  vcm->group_j[g][m] += m0*j0[m];
	  vcm->group_x[g][m] += m0*x[i][m];
	}
	/* Update inertia tensor */
	update_tensor(x[i],m0,vcm->group_i[g]);
      }
    }
  }
}

void do_stopcm_grp(FILE *fp,int start,int homenr,rvec x[],rvec v[],
		   t_vcm *vcm)
{
  int  i,g,m;
  real tm,tm_1;
  rvec dv,dx;
  
  if (vcm->mode != ecmNO) {
    /* Subtract linear momentum */
    for(i=start; (i<start+homenr); i++) {
      g = vcm->group_id[i];
      rvec_dec(v[i],vcm->group_v[g]);
    }
    if (vcm->mode == ecmANGULAR) {
      /* Subtract angular momentum */
      for(i=start; (i<start+homenr); i++) {
	g = vcm->group_id[i];
	/* Compute the correction to the velocity for each atom */
	rvec_sub(x[i],vcm->group_x[g],dx);
	oprod(vcm->group_w[g],dx,dv);
	rvec_dec(v[i],dv);
      }
    }
  }
}

static void get_minv(tensor A,tensor B)
{
  int    m,n;
  double fac,rfac;
  tensor tmp;

  tmp[XX][XX] =  A[YY][YY] + A[ZZ][ZZ];
  tmp[YY][XX] = -A[XX][YY];
  tmp[ZZ][XX] = -A[XX][ZZ];
  tmp[XX][YY] = -A[XX][YY];
  tmp[YY][YY] =  A[XX][XX] + A[ZZ][ZZ];
  tmp[ZZ][YY] = -A[YY][ZZ];
  tmp[XX][ZZ] = -A[XX][ZZ];
  tmp[YY][ZZ] = -A[YY][ZZ];
  tmp[ZZ][ZZ] =  A[XX][XX] + A[YY][YY];
  
  /* This is a hack to prevent very large determinants */
  rfac  = (tmp[XX][XX]+tmp[YY][YY]+tmp[ZZ][ZZ])/3;
  if (rfac == 0.0) 
    gmx_fatal(FARGS,"Can not stop center of mass: maybe 2dimensional system");
  fac = 1.0/rfac;
  for(m=0; (m<DIM); m++)
    for(n=0; (n<DIM); n++)
      tmp[m][n] *= fac;
  m_inv(tmp,B);
  for(m=0; (m<DIM); m++)
    for(n=0; (n<DIM); n++)
      B[m][n] *= fac;
}

void check_cm_grp(FILE *fp,t_vcm *vcm)
{
  int    m,g;
  real   ekcm,ekrot,tm,tm_1;
  rvec   jcm;
  tensor Icm,Tcm;
    
  /* First analyse the total results */
  if (vcm->mode != ecmNO) {
    for(g=0; (g<vcm->nr); g++) {
      tm = vcm->group_mass[g];
      if (tm != 0) {
	tm_1 = 1.0/tm;
	svmul(tm_1,vcm->group_p[g],vcm->group_v[g]);
      }
      /* Else it's zero anyway! */
    }
    if (vcm->mode == ecmANGULAR) {
      for(g=0; (g<vcm->nr); g++) {
	tm = vcm->group_mass[g];
	if (tm != 0) {
	  tm_1 = 1.0/tm;
	  
	  /* Compute center of mass for this group */
	  for(m=0; (m<DIM); m++)
	    vcm->group_x[g][m] *= tm_1;
	  
	  /* Subtract the center of mass contribution to the 
	   * angular momentum 
	   */
 	  oprod(vcm->group_x[g],vcm->group_v[g],jcm);
	  for(m=0; (m<DIM); m++)
	    vcm->group_j[g][m] -= tm*jcm[m];
	  
	  /* Subtract the center of mass contribution from the inertia 
	   * tensor (this is not as trivial as it seems, but due to
	   * some cancellation we can still do it, even in parallel).
	   */
	  clear_mat(Icm);
	  update_tensor(vcm->group_x[g],tm,Icm);
	  m_sub(vcm->group_i[g],Icm,vcm->group_i[g]);
	  
	  /* Compute angular velocity, using matrix operation 
	   * Since J = I w
	   * we have
	   * w = I^-1 J
	   */
	  get_minv(vcm->group_i[g],Icm);
	  mvmul(Icm,vcm->group_j[g],vcm->group_w[g]);
	}
	/* Else it's zero anyway! */
      }
    }
  }
  for(g=0; (g<vcm->nr); g++) {
    ekcm    = 0;
    if (vcm->group_mass[g] != 0) {
      for(m=0; (m<DIM); m++) 
	ekcm += sqr(vcm->group_v[g][m]);
      ekcm *= 0.5*vcm->group_mass[g];
      
      if ((ekcm > 1) || debug)
	fprintf(fp,"Large VCM(group %s): %12.5f, %12.5f, %12.5f, ekin-cm: %12.5e\n",
		vcm->group_name[g],vcm->group_v[g][XX],
		vcm->group_v[g][YY],vcm->group_v[g][ZZ],ekcm);
      
      if (vcm->mode == ecmANGULAR) {
	ekrot = 0.5*iprod(vcm->group_j[g],vcm->group_w[g]);
	if ((ekrot > 1) || debug) {
	  tm    = vcm->group_mass[g];
	  fprintf(fp,"Group %s with mass %12.5e, Ekrot %12.5e Det(I) = %12.5e\n",
		  vcm->group_name[g],tm,ekrot,det(vcm->group_i[g]));
	  fprintf(fp,"  COM: %12.5f  %12.5f  %12.5f\n",
		  vcm->group_x[g][XX],vcm->group_x[g][YY],vcm->group_x[g][ZZ]);
	  fprintf(fp,"  P:   %12.5f  %12.5f  %12.5f\n",
		  vcm->group_p[g][XX],vcm->group_p[g][YY],vcm->group_p[g][ZZ]);
	  fprintf(fp,"  V:   %12.5f  %12.5f  %12.5f\n",
		  vcm->group_v[g][XX],vcm->group_v[g][YY],vcm->group_v[g][ZZ]);
	  fprintf(fp,"  J:   %12.5f  %12.5f  %12.5f\n",
		  vcm->group_j[g][XX],vcm->group_j[g][YY],vcm->group_j[g][ZZ]);
	  fprintf(fp,"  w:   %12.5f  %12.5f  %12.5f\n",
		vcm->group_w[g][XX],vcm->group_w[g][YY],vcm->group_w[g][ZZ]);
	  pr_rvecs(fp,0,"Inertia tensor",vcm->group_i[g],DIM);
	}
      }
    }
  }
}


