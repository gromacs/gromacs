/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "typedefs.h"
#include "gromacs/utility/smalloc.h"
#include "macros.h"
#include "copyrite.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/utility/fatalerror.h"
#include "random.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/utility/futil.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "names.h"
#include "ehdata.h"

typedef struct {
  int  maxparticle;
  int  maxstep;
  int  nsim;
  int  nsave;
  int  nana;
  int  seed;
  int  nevent;
  gmx_bool bForce;
  gmx_bool bScatter;
  gmx_bool bHole;
  real dt;
  real deltax;
  real epsr;
  real Alj;
  real Eauger;
  real Efermi;
  real Eband;
  real rho;
  real matom;
  real evdist;
  real size;
} t_eh_params;

#define ELECTRONMASS 5.447e-4
/* Resting mass of electron in a.m.u. */
#define HOLEMASS (0.8*ELECTRONMASS)
/* Effective mass of a hole! */
#define HBAR (PLANCK/2*M_PI)

static void calc_forces(int n,rvec x[],rvec f[],real q[],real ener[],real Alj)
{
  const real facel = FACEL;
  int  i,j,m;
  rvec dx;
  real qi,r2,r_1,r_2,fscal,df,vc,vctot,vlj,vljtot;
  
  for(i=0; (i<n); i++) 
    clear_rvec(f[i]);
  
  vctot = vljtot = 0;
  for(i=0; (i<n-1); i++) {
    qi = q[i]*facel;
    for(j=i+1; (j<n); j++) {
      rvec_sub(x[i],x[j],dx);
      r2      = iprod(dx,dx);
      r_1     = 1.0/sqrt(r2);
      r_2     = r_1*r_1;
      vc      = qi*q[j]*r_1;
      vctot  += vc;
      vlj     = Alj*(r_2*r_2*r_2);
      vljtot += vlj;
      fscal   = (6*vlj+vc)*r_2;
      for(m=0; (m<DIM); m++) {
	df = fscal*dx[m];
	f[i][m] += df;
	f[j][m] -= df;
      }
    }
  }
  ener[eCOUL]   = vctot;
  ener[eREPULS] = vljtot;
  ener[ePOT]    = vctot+vljtot;
}

static void calc_ekin(int nparticle,rvec v[],rvec vold[],
		      real q[],real m[],real ener[],real eparticle[])
{
  rvec vt;
  real ekh=0,eke=0,ee;
  int  i;
  
  for(i=0; (i<nparticle); i++) {
    rvec_add(v[i],vold[i],vt);
    ee = 0.125*m[i]*iprod(vt,vt);
    eparticle[i] = ee/ELECTRONVOLT;
    if (q[i] > 0)
      ekh += ee;
    else 
      eke += ee;
  }
  ener[eHOLE]     = ekh;
  ener[eELECTRON] = eke;
  ener[eKIN]      = ekh+eke+ener[eLATTICE];
}

static void polar2cart(real amp,real phi,real theta,rvec v)
{
  real ss = sin(theta);
  
  v[XX] = amp*cos(phi)*ss;
  v[YY] = amp*sin(phi)*ss;
  v[ZZ] = amp*cos(theta);
}

static void rand_vector(real amp,rvec v,int *seed)
{
  real theta,phi;

  theta = M_PI*rando(seed);
  phi   = 2*M_PI*rando(seed);
  polar2cart(amp,phi,theta,v);
}

static void rotate_theta(rvec v,real nv,real dth,int *seed,FILE *fp)
{
  real   dphi,theta0,phi0,cc,ss;
  matrix mphi,mtheta,mphi_1,mtheta_1; 
  rvec   vp,vq,vold;
  
  copy_rvec(v,vold);
  theta0 = acos(v[ZZ]/nv);
  phi0   = atan2(v[YY],v[XX]);
  if (fp)
    fprintf(fp,"Theta = %g  Phi = %g\n",theta0,phi0);
    
  clear_mat(mphi);
  cc = cos(-phi0);
  ss = sin(-phi0);
  mphi[XX][XX] = mphi[YY][YY] = cc;
  mphi[XX][YY] = -ss;
  mphi[YY][XX] = ss;
  mphi[ZZ][ZZ] = 1;
  m_inv(mphi,mphi_1);

  clear_mat(mtheta);
  cc = cos(-theta0);
  ss = sin(-theta0);
  mtheta[XX][XX] = mtheta[ZZ][ZZ] = cc;
  mtheta[XX][ZZ] = ss;
  mtheta[ZZ][XX] = -ss;
  mtheta[YY][YY] = 1;
  m_inv(mtheta,mtheta_1);
  
  dphi   = 2*M_PI*rando(seed);
  
  /* Random rotation */
  polar2cart(nv,dphi,dth,vp);
  
  mvmul(mtheta_1,vp,vq);
  mvmul(mphi_1,vq,v);
  
  if (fp) {
    real cold = cos_angle(vold,v);
    real cnew = cos(dth);
    if (fabs(cold-cnew) > 1e-4)
      fprintf(fp,"cos(theta) = %8.4f  should be %8.4f  dth = %8.4f  dphi = %8.4f\n",
	      cold,cnew,dth,dphi);
  }
}

static int create_electron(int index,rvec x[],rvec v[],rvec vold[],rvec vv,
			   real m[],real q[],
			   rvec center,real e0,int *seed)
{
  m[index] = ELECTRONMASS;
  q[index] = -1;

  clear_rvec(v[index]);
  svmul(sqrt(2*e0/m[index]),vv,v[index]);
  copy_rvec(v[index],vold[index]);
  copy_rvec(center,x[index]);
  
  return index+1;
}

static int create_pair(int index,rvec x[],rvec v[],rvec vold[],
		       real m[],real q[],
		       rvec center,real e0,t_eh_params *ehp,rvec dq)
{
  static real massfactor = 2*HOLEMASS/(ELECTRONMASS*(ELECTRONMASS+HOLEMASS));
  rvec x0;
  real ve,e1;
  
  m[index]        = ELECTRONMASS;
  m[index+1]      = HOLEMASS;
  q[index]        = -1;
  q[index+1]      = 1;
  
  rand_vector(0.5*ehp->deltax,x0,&ehp->seed);
  rvec_sub(center,x0,x[index]);
  rvec_add(center,x0,x[index+1]);

  ve = sqrt(massfactor*e0)/(0.5*ehp->deltax);
  svmul(-ve,x0,v[index]);
  svmul(ELECTRONMASS*ve/HOLEMASS,x0,v[index+1]);
  copy_rvec(v[index],vold[index]);
  copy_rvec(v[index+1],vold[index+1]);
  e1 = 0.5*(m[index]*iprod(v[index],v[index])+
	    m[index+1]*iprod(v[index+1],v[index+1]));
  if (fabs(e0-e1)/e0 > 1e-6)
    gmx_fatal(FARGS,"Error in create pair: e0 = %f, e1 = %f\n",e0,e1);
  
  return index+2;
}

static int scatter_all(FILE *fp,int nparticle,int nstep,
		       rvec x[],rvec v[],rvec vold[],
		       real mass[],real charge[],real ener[],real eparticle[],
		       t_eh_params *ehp,int *nelec,int *nhole,t_ana_scat s[])
{
  int  i,m,np;
  real p_el,p_inel,ptot,nv,ekin,omega,theta,costheta,Q,e0,ekprime,size2,fac;
  rvec dq,center,vv;
  
  size2 = sqr(ehp->size);
  np    = nparticle;  
  for(i=0; (i<nparticle); i++) {
    /* Check cross sections, assume same cross sections for holes
     * as for electrons, for elastic scattering
     */
    if ((size2 == 0) || (iprod(x[i],x[i]) < size2)) {
      nv   = norm(v[i]);
      ekin = eparticle[i];
      p_el = cross_el(ekin,ehp->rho,NULL)*nv*ehp->dt;
      
      /* Only electrons can scatter inelasticlly */
      if (charge[i] < 0)
	p_inel = cross_inel(ekin,ehp->rho,NULL)*nv*ehp->dt;
      else
	p_inel = 0;
      
      /* Test whether we have to scatter at all */
      ptot = (1 - (1-p_el)*(1-p_inel));
      if (debug && 0)
	fprintf(debug,"p_el = %10.3e  p_inel = %10.3e ptot = %10.3e\n",
		p_el,p_inel,ptot);
      if (rando(&ehp->seed) < ptot) {
	/* Test whether we have to scatter inelastic */
	ptot = p_inel/(p_el+p_inel);
	if (rando(&ehp->seed) < ptot) {
	  add_scatter_event(&(s[i]),x[i],TRUE,ehp->dt*nstep,ekin);
	  /* Energy loss in inelastic collision is omega */
	  if ((omega = get_omega(ekin,&ehp->seed,debug,NULL)) >= ekin)
	    gmx_fatal(FARGS,"Energy transfer error: omega = %f, ekin = %f",
			omega,ekin);
	  else {
	    /* Scattering angle depends on energy and energy loss */
	    Q = get_q_inel(ekin,omega,&ehp->seed,debug,NULL);
	    costheta = -0.5*(Q+omega-2*ekin)/sqrt(ekin*(ekin-omega));
	    
	    /* See whether we have gained enough energy to liberate another 
	     * hole-electron pair
	     */
	    e0      = band_ener(&ehp->seed,debug,NULL);
	    ekprime = e0 + omega - (ehp->Efermi+0.5*ehp->Eband);
	    /* Ouput */
	    fprintf(fp,"Inelastic %d: Ekin=%.2f Omega=%.2f Q=%.2f Eband=%.2f costheta=%.3f\n",
		    i+1,ekin,omega,Q,e0,costheta);
	    if ((costheta < -1) || (costheta > 1)) {
	      fprintf(fp,"Electron/hole creation not possible due to momentum constraints\n");
	      /* Scale the velocity according to the energy loss */
	      svmul(sqrt(1-omega/ekin),v[i],v[i]);
	      ener[eLATTICE] += omega*ELECTRONVOLT;
	    }
	    else {
	      theta = acos(costheta);
	      
	      copy_rvec(v[i],dq);
	      /* Rotate around theta with random delta phi */
	      rotate_theta(v[i],nv,theta,&ehp->seed,debug);
	      /* Scale the velocity according to the energy loss */
	      svmul(sqrt(1-omega/ekin),v[i],v[i]);
	      rvec_dec(dq,v[i]);
	      
	      if (ekprime > 0) {
		if (np >= ehp->maxparticle-2)
		  gmx_fatal(FARGS,"Increase -maxparticle flag to more than %d",
			      ehp->maxparticle);
		if (ehp->bHole) {
		  np = create_pair(np,x,v,vold,mass,charge,x[i],
				   ekprime*ELECTRONVOLT,ehp,dq);
		  (*nhole)++;
		}
		else {
		  copy_rvec(x[i],center);
		  center[ZZ] += ehp->deltax;
		  rand_vector(1,vv,&ehp->seed);
		  np = create_electron(np,x,v,vold,vv,mass,charge,
				       x[i],ekprime*ELECTRONVOLT,&ehp->seed);
		}
		ener[eLATTICE] += (omega-ekprime)*ELECTRONVOLT;
		(*nelec)++;
	      }
	      else
		ener[eLATTICE] += omega*ELECTRONVOLT;
	    }
	  }
	}
	else {
	  add_scatter_event(&(s[i]),x[i],FALSE,ehp->dt*nstep,ekin);
	  if (debug)
	    fprintf(debug,"Elastic scattering event\n");
	  
	  /* Scattering angle depends on energy only */
	  theta = get_theta_el(ekin,&ehp->seed,debug,NULL);
	  /* Rotate around theta with random delta phi */
	  rotate_theta(v[i],nv,theta,&ehp->seed,debug);
	}
      }
    }
  }
  return np;
}

static void integrate_velocities(int nparticle,rvec vcur[],rvec vnext[],
				 rvec f[],real m[],real dt)
{
  int i,k;
    
  for(i=0; (i<nparticle); i++) 
    for(k=0; (k<DIM); k++) 
      vnext[i][k] = vcur[i][k] + f[i][k]*dt/m[i];
}

static void integrate_positions(int nparticle,rvec x[],rvec v[],real dt)
{
  int i,k;
  
  for(i=0; (i<nparticle); i++) 
    for(k=0; (k<DIM); k++) 
      x[i][k] += v[i][k]*dt;
}

static void print_header(FILE *fp,t_eh_params *ehp)
{
  fprintf(fp,"Welcome to the electron-hole simulation!\n");
  fprintf(fp,"The energies printed in this file are in eV\n");
  fprintf(fp,"Coordinates are in nm because of fixed width format\n");
  fprintf(fp,"Atomtypes are used for coloring in rasmol\n");
  fprintf(fp,"O: electrons (red), N: holes (blue)\n");
  fprintf(fp,"Parametes for this simulation\n");
  fprintf(fp,"seed = %d maxstep = %d dt = %g\n",
	  ehp->seed,ehp->maxstep,ehp->dt);
  fprintf(fp,"nsave = %d nana = %d Force = %s Scatter = %s Hole = %s\n",
	  ehp->nsave,ehp->nana,gmx_bool_names[ehp->bForce],
	  gmx_bool_names[ehp->bScatter],gmx_bool_names[ehp->bHole]);
  if (ehp->bForce)
    fprintf(fp,"Force constant for repulsion Alj = %g\n",ehp->Alj);
}

static void do_sim(FILE *fp,char *pdbfn,t_eh_params *ehp,
		   int *nelec,int *nhole,t_ana_struct *total,
		   t_histo *hmfp,t_ana_ener *ae,int serial)
{
  FILE         *efp;
  int          nparticle[2];
  rvec         *x,*v[2],*f,center,vv;
  real         *charge,*mass,*ener,*eparticle;
  t_ana_struct *ana_struct;
  t_ana_scat   *ana_scat;
  int          step,i,cur = 0;
#define next (1-cur)

  /* Open output file */
  fprintf(fp,"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  fprintf(fp,"Simulation %d/%d\n",serial+1,ehp->nsim);
  
  ana_struct = init_ana_struct(ehp->maxstep,ehp->nana,ehp->dt,
			       ehp->maxparticle);
  /* Initiate arrays. The charge array determines whether a particle is 
   * a hole (+1) or an electron (-1)
   */
  snew(x,ehp->maxparticle);          /* Position  */
  snew(v[0],ehp->maxparticle);       /* Velocity  */
  snew(v[1],ehp->maxparticle);       /* Velocity  */
  snew(f,ehp->maxparticle);          /* Force     */
  snew(charge,ehp->maxparticle);     /* Charge    */
  snew(mass,ehp->maxparticle);       /* Mass      */
  snew(eparticle,ehp->maxparticle);  /* Energy per particle */
  snew(ana_scat,ehp->maxparticle);   /* Scattering event statistics */
  snew(ener,eNR);                    /* Eenergies */
  
  clear_rvec(center);
  /* Use first atom as center, it has coordinate 0,0,0 */
  if (ehp->bScatter) {
    /* Start with an Auger electron */
    nparticle[cur]=0;
    for(i=0; (i<ehp->nevent); i++) {
      if (ehp->nevent == 1) {
	clear_rvec(vv);
	vv[ZZ] = 1;
      }
      else
	rand_vector(1,vv,&ehp->seed);
      nparticle[cur]  = create_electron(nparticle[cur],x,v[cur],v[next],
					vv,mass,charge,center,
					ehp->Eauger*ELECTRONVOLT,&ehp->seed);
      rand_vector(ehp->evdist*0.1,vv,&ehp->seed);
      rvec_inc(center,vv);
    }
  }
  else if (ehp->bForce) {
    /* Start with two electron and hole pairs */
    nparticle[cur]  = create_pair(0,x,v[cur],v[next],mass,charge,center,
				  0.2*ehp->Eauger*ELECTRONVOLT,ehp,center);
    center[ZZ] = 0.5; /* nm */
    (*nelec)++;
    (*nhole)++;
  }
  else {
    fprintf(fp,"Nothing to do. Doei.\n");
    return;
  }
  nparticle[next] = nparticle[cur];
  for(step=0; (step<=ehp->maxstep); step++) {
    if (ehp->bScatter)
      nparticle[next] = scatter_all(fp,nparticle[cur],step,x,v[cur],v[next],
				    mass,charge,ener,eparticle,ehp,
				    nelec,nhole,ana_scat);
    
    if (ehp->bForce)
      calc_forces(nparticle[cur],x,f,charge,ener,ehp->Alj);
    
    integrate_velocities(nparticle[next],v[cur],v[next],f,mass,ehp->dt);
    
    calc_ekin(nparticle[next],v[cur],v[next],charge,mass,ener,eparticle);
    ener[eTOT] = ener[eKIN] + ener[ePOT];
    
    /* Produce output whenever the user says so, or when new
     * particles have been created.
     */
    if ((step == ehp->maxstep) ||
	((ehp->nana != 0) && ((step % ehp->nana) == 0))) {
      analyse_structure(ana_struct,(step*ehp->dt),center,x,
			nparticle[next],charge);
      add_ana_ener(ae,(step/ehp->nana),ener);
    }
    cur = next;
        
    integrate_positions(nparticle[cur],x,v[cur],ehp->dt);
  }
  for(i=0; (i<nparticle[cur]); i++) {
    analyse_scatter(&(ana_scat[i]),hmfp);
    done_scatter(&(ana_scat[i]));
  }
  sfree(ener);
  sfree(ana_scat); 
  sfree(eparticle); 
  sfree(mass);    
  sfree(charge); 
  sfree(f);
  sfree(v[1]);      
  sfree(v[0]); 
  sfree(x);
  dump_as_pdb(pdbfn,ana_struct);
  add_ana_struct(total,ana_struct);
  done_ana_struct(ana_struct);
  sfree(ana_struct);
}

void do_sims(int NFILE,t_filenm fnm[],t_eh_params *ehp)
{
  t_ana_struct *total;
  t_ana_ener   *ae;
  t_histo      *helec,*hmfp;
  int          *nelectron;
  int          i,imax,ne,nh;
  real         aver;
  FILE         *fp,*logfp;
  char         *pdbbuf,*ptr,*rptr;

  ptr  = ftp2fn(efPDB,NFILE,fnm);
  rptr = strdup(ptr);
  if ((ptr = strstr(rptr,".pdb")) != NULL)
    *ptr = '\0';
  snew(pdbbuf,strlen(rptr)+10);

  total = init_ana_struct(ehp->maxstep,ehp->nana,ehp->dt,1);
  hmfp  = init_histo((int)ehp->Eauger,0,(int)ehp->Eauger);
  helec = init_histo(500,0,500);
  snew(ae,1);

  logfp = gmx_ffopen(ftp2fn(efLOG,NFILE,fnm),"w");
  print_header(logfp,ehp);
    
  for(i=0; (i<ehp->nsim); i++) {
    nh = ne = 0;
    sprintf(pdbbuf,"%s-%d.pdb",rptr,i+1);
    do_sim(logfp,pdbbuf,ehp,&ne,&nh,total,hmfp,ae,i);
    add_histo(helec,ne,1);
    fprintf(stderr,"\rSim: %d/%d",i+1,ehp->nsim);
  }
  fprintf(stderr,"\n");
  gmx_ffclose(logfp);
  
  sfree(rptr);
  sfree(pdbbuf);
  dump_ana_struct(opt2fn("-maxdist",NFILE,fnm),opt2fn("-nion",NFILE,fnm),
		  opt2fn("-gyr_com",NFILE,fnm),opt2fn("-gyr_origin",NFILE,fnm),
		  total,ehp->nsim);
  dump_ana_ener(ae,ehp->nsim,ehp->dt*ehp->nana,
		opt2fn("-ener",NFILE,fnm),total);
  done_ana_struct(total);

  dump_histo(helec,opt2fn("-histo",NFILE,fnm),
	     "Number of cascade electrons","N","",enormFAC,1.0/ehp->nsim);
  dump_histo(hmfp,opt2fn("-mfp",NFILE,fnm),
	     "Mean Free Path","Ekin (eV)","MFP (nm)",enormNP,1.0);
}

int main(int argc,char *argv[])
{
  const char *desc[] = {
    "[TT]ehole[tt] performs a molecular dynamics simulation of electrons and holes",
    "in an implicit lattice. The lattice is modeled through scattering cross",
    "sections, for elastic and inelastic scattering.",
    "A detailed description of the scatterning processes simulated in ehole",
    "can be found in Timneanu et al. Chemical Physics 299 (2004) 277-283",
    "The paper also includes a description how to calculate the input files.[PAR]",
    "Description of the input files for [TT]ehole[tt]:[BR]",
    "[TT]-sigel.dat[tt]: elastic cross section (per atom). Two columns: Impact electron energy (eV) vs Elastic cross section (A2).[BR]",
    "[TT]-siginel.dat[tt]: inelastic cross section (per atom). Two columns: Impact electron energy (eV) vs Inelastic cross section (A2).[BR]",
    "[TT]-band-ener.dat[tt]: Probability of finding an electron in the valence band.",
    "Two columns: Impact electron energy (eV) vs Probability[BR]",
    "[TT]-eloss.dat[tt]: Probability of energy loss due to inelastic scattering. Three columns: Impact electron energy (eV) vs  Integrated probability vs Energy loss in inelastic scattering (eV).[BR]",
    "[TT]-theta-el.dat[tt]: Probability of elastic scattering angle. Three columns: Impact electron energy (eV) vs Integrated probability vs Scattering angle (rad).[BR]",
    "[TT]-qtrans.dat[tt]: Four columns: Impact electron energy (eV) vs Inelastic energy loss (eV) vs Integrated probability vs Scattering angle (rad).[PAR]",
    "The program produces a number of output files. It is important that",
    "the actual content is well-defined, sucht that no misunderstanding can",
    "occur (famous last words...). Anyway, the program does a number of",
    "simulations, and averages results over these. Here is a list of each of",
    "the results and how they are computed:[BR]",
    "[TT]-histo[tt] Distribution of nuber of liberated secondary electrons per simulation.[BR]",
    "[TT]-maxdist[tt] The maximum distance from the origin that any electron in any simulation reaches.[BR]",
    "[TT]-gyr_com[tt] The radius of gyration of the electron cloud with respect to its center of mass (contains 4 columns).[BR]",
    "[TT]-gyr_origin[tt] The radius of gyration of the electron cloud with respect to the origin (contains 4 columns).[BR]",
    "[TT]-mfp[tt] The mean free path of the electrons as a function of energy. If this is not a smooth curve you need to increase the number of simulations.[BR]",
    "[TT]-nion[tt] The number of ions as a function of time, averaged over simulations.[BR]",
    "[TT]-ener[tt] The energy terms in the simulation (note that there are multiple columns, so use [TT]xmgrace -nxy[tt]). This shows important information about the stability of the simulation, that is the total energy should be conserved. In this figure you can also inspect the kinetic energy per electron in order to check whether the electrons have thermalized.[BR]"
  };
  static t_eh_params ehp = {
    100,    /* Max number of particles. Is a parameter but should be dynamic */
    100000, /* Number of integration steps */
    1,      /* nsave */
    1,      /* nana */
    1,      /* nsim */
    1993,   /* Random seed */
    1,      /* Number of events */
    FALSE,  /* Use forces */
    TRUE,   /* Use scattering */
    FALSE,  /* Creat holes */
    1e-5,   /* Time step */
    0.05,   /* Distance (nm) between electron and hole when creating them */
    1.0,    /* Dielectric constant */
    0.1,    /* Force constant for repulsion function */
    250,    /* Starting energy for the first Auger electron */
    28.7,   /* Fermi level (eV) of diamond. */
    5.46,   /* Band gap energy (eV) of diamond */
    3.51,   /* Density of the solid */
    12.011, /* (Average) mass of the atom */
    10000.0,/* Distance between events */
    0.0     /* Size of the system */
  };
  static gmx_bool bTest    = FALSE;
  t_pargs pa[] = {
    { "-maxparticle", FALSE, etINT,  {&ehp.maxparticle},
      "Maximum number of particles" },
    { "-maxstep",     FALSE, etINT,  {&ehp.maxstep}, 
      "Number of integration steps" },
    { "-nsim",        FALSE, etINT,  {&ehp.nsim},
      "Number of independent simulations writing to different output files" },
    { "-nsave",       FALSE, etINT,  {&ehp.nsave}, 
      "Number of steps after which to save output. 0 means only when particles created. Final step is always written." },
    { "-nana",        FALSE, etINT,  {&ehp.nana}, 
      "Number of steps after which to do analysis." },
    { "-seed",        FALSE, etINT,  {&ehp.seed}, 
      "Random seed" },
    { "-dt",          FALSE, etREAL, {&ehp.dt}, 
      "Integration time step (ps)" },
    { "-rho",         FALSE, etREAL, {&ehp.rho}, 
      "Density of the sample (kg/l). Default is for Diamond" }, 
    { "-matom",       FALSE, etREAL, {&ehp.matom}, 
      "Mass (a.m.u.) of the atom in the solid. Default is C" },
    { "-fermi",       FALSE, etREAL, {&ehp.Efermi}, 
      "Fermi energy (eV) of the sample. Default is for Diamond" },
    { "-gap",         FALSE, etREAL, {&ehp.Eband}, 
      "Band gap energy (eV) of the sample. Default is for Diamond" },
    { "-auger",       FALSE, etREAL, {&ehp.Eauger}, 
      "Impact energy (eV) of first electron" },
    { "-dx",          FALSE, etREAL, {&ehp.deltax},
      "Distance between electron and hole when creating a pair" },
    { "-test",        FALSE, etBOOL, {&bTest},
      "Test table aspects of the program rather than running it for real" },
    { "-force",       FALSE, etBOOL, {&ehp.bForce},
      "Apply Coulomb/Repulsion forces" },
    { "-hole",        FALSE, etBOOL, {&ehp.bHole},
      "Create electron-hole pairs rather than electrons only" },
    { "-scatter",     FALSE, etBOOL, {&ehp.bScatter},
      "Do the scattering events" },
    { "-nevent",      FALSE, etINT,  {&ehp.nevent},
      "Number of initial Auger electrons" },
    { "-evdist",      FALSE, etREAL, {&ehp.evdist},
      "Average distance (A) between initial electronss" },
    { "-size",        FALSE, etREAL, {&ehp.size},
      "Size of the spherical system. If 0, then it is infinite" }
  };
#define NPA asize(pa)
  t_filenm fnm[] = {
    { efLOG, "-g",          "ehole",      ffWRITE },
    { efDAT, "-sigel",      "sigel",      ffREAD },
    { efDAT, "-sigin",      "siginel",    ffREAD },
    { efDAT, "-eloss",      "eloss",      ffREAD },
    { efDAT, "-qtrans",     "qtrans",     ffREAD },
    { efDAT, "-band",       "band-ener",  ffREAD },
    { efDAT, "-thetael",    "theta-el",   ffREAD },
    { efPDB, "-o",          "ehole",      ffWRITE },
    { efXVG, "-histo",      "histo",      ffWRITE },
    { efXVG, "-maxdist",    "maxdist",    ffWRITE },
    { efXVG, "-gyr_com",    "gyr_com",    ffWRITE },
    { efXVG, "-gyr_origin", "gyr_origin", ffWRITE },
    { efXVG, "-mfp",        "mfp",        ffWRITE },
    { efXVG, "-nion",       "nion",       ffWRITE },
    { efXVG, "-ener",       "ener",       ffWRITE }
  };
#define NFILE asize(fnm)
  int seed;
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,
		    NPA,pa,asize(desc),desc,0,NULL);
  please_cite(stdout,"Timneanu2004a");
  
  if (ehp.deltax <= 0)
    gmx_fatal(FARGS,"Delta X should be > 0");
  ehp.Alj = FACEL*pow(ehp.deltax,5);
  ehp.rho = (ehp.rho/ehp.matom)*AVOGADRO*1e-21;

  init_tables(NFILE,fnm,ehp.rho);

  if (bTest) 
    test_tables(&ehp.seed,ftp2fn(efPDB,NFILE,fnm),ehp.rho);  
  else 
    do_sims(NFILE,fnm,&ehp);
  
  gmx_thanx(stdout);
  
  return 0;
}
