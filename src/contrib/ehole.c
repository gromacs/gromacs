#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "copyrite.h"
#include "statutil.h"
#include "fatal.h"
#include "random.h"
#include "pdbio.h"
#include "futil.h"
#include "physics.h"
#include "xvgr.h"
#include "vec.h"
#include "names.h"
#include "ehdata.h"


#define ELECTRONMASS 5.447e-4
/* Resting mass of electron in a.m.u. */
#define HOLEMASS (0.8*ELECTRONMASS)
/* Effective mass of a hole! */
#define PLANCK 6.62606891e-34; 
/* Planck's constant */
#define HBAR (PLANCK/2*M_PI)

static real DELTAX = 0.05;
/* Distance (nm) between electron and hole when creating them */
static real Alj    = 0.1;
/* Force constant for repulsion function */
static real AUGER_ENERGY = 250;
/* Starting energy for the first Auger electron */
static real EFERMI = 28.7;
/* Fermi level (eV) of diamond. */
static real EBAND  = 5.46;
/* Band gap energy (eV) of diamond */
static int MAXPARTICLE=100;
/* Max number of particles. Is a parameter but should be dynamic */

enum { eCOUL, eREPULS, ePOT, eHOLE, eELECTRON, eLATTICE, eKIN, eTOT, eNR };

static void calc_forces(int n,rvec x[],rvec f[],real q[],real ener[])
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

static void write_data(FILE *fp,int step,real t,
		       int nparticle,rvec x[],rvec v[],rvec f[],
		       real q[],real ener[],real eparticle[])
{
  static char *enms[eNR] = {
    "Coulomb", "Repulsion", "Potential",
    "EkHole",  "EkElectron", "EkLattice", "Kinetic",
    "Total"
  };
  int i,j;
  
  fprintf(fp,"MODEL %d\n",step+1);
  fprintf(fp,"REMARK time = %10.3f fs, nparticle = %4d\n",1000*t,nparticle);
  for(i=0; (i<eNR); i+=2) {
    fprintf(fp,"REMARK");
    for(j=i; (j<eNR) && (j<i+2); j++)
      fprintf(fp," %12s = %12.5e",enms[j],ener[j]/ELECTRONVOLT);
    fprintf(fp,"\n");
  }
#ifdef DEBUG  
  fprintf(fp,"REMARK Coordinates\n");
#endif
  for(i=0; (i<nparticle); i++) {
    fprintf(fp,pdbformat,"ATOM",i+1,(q[i] < 0) ? "O" : "N","PLS",' ',1,
	    x[i][XX],x[i][YY],x[i][ZZ]);
    fprintf(fp,"  %8.3f\n",eparticle[i]);
  }
#ifdef DEBUG1
  fprintf(fp,"REMARK Velocities\n");
  for(i=0; (i<nparticle); i++) {
    fprintf(fp,pdbformat,"ATOM",i+1,(q[i] < 0) ? "O" : "N","PLS",' ',1,
	    v[i][XX],v[i][YY],v[i][ZZ]);
    fprintf(fp,"\n");
  }
  fprintf(fp,"REMARK Forces\n");
  for(i=0; (i<nparticle); i++) {
    fprintf(fp,pdbformat,"ATOM",i+1,(q[i] < 0) ? "O" : "N","PLS",' ',1,
	    f[i][XX],f[i][YY],f[i][ZZ]);
    fprintf(fp,"\n");
  }
#endif
  fprintf(fp,"ENDMDL\n");
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
      eke += ee;
    else 
      ekh += ee;
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

static int create_auger(int index,rvec x[],rvec v[],rvec vold[],
			real m[],real q[],
			rvec center,real e0,int *seed)
{
  m[index] = ELECTRONMASS;
  q[index] = -1;

  clear_rvec(v[index]);
  v[index][ZZ] = sqrt(2*e0/m[index]);
  copy_rvec(v[index],vold[index]);
  copy_rvec(center,x[index]);
  
  return index+1;
}

static int create_pair(int index,rvec x[],rvec v[],rvec vold[],
		       real m[],real q[],
		       rvec center,real e0,int *seed)
{
  static real massfactor = 2*HOLEMASS/(ELECTRONMASS*(ELECTRONMASS+HOLEMASS));
  rvec x0;
  real ve,e1;
  
  m[index]        = ELECTRONMASS;
  m[index+1]      = HOLEMASS;
  q[index]        = -1;
  q[index+1]      = 1;
  
  rand_vector(0.5*DELTAX,x0,seed);
  rvec_sub(center,x0,x[index]);
  rvec_add(center,x0,x[index+1]);

  ve = sqrt(massfactor*e0)/(0.5*DELTAX);
  svmul(-ve,x0,v[index]);
  svmul(ELECTRONMASS*ve/HOLEMASS,x0,v[index+1]);
  copy_rvec(v[index],vold[index]);
  copy_rvec(v[index+1],vold[index+1]);
  e1 = 0.5*(m[index]*iprod(v[index],v[index])+
	    m[index+1]*iprod(v[index+1],v[index+1]));
  if (fabs(e0-e1)/e0 > 1e-6)
    fatal_error(0,"Error in create pair: e0 = %f, e1 = %f\n",e0,e1);
  
  return index+2;
}

static int scatter_all(FILE *fp,int nparticle,
		       rvec x[],rvec v[],rvec vold[],
		       real mass[],real charge[],real ener[],real eparticle[],
		       int *seed,real dt,int nstep,int *nelec,int *nhole,t_ana_scat s[])
{
  int  i,m,np;
  real p_el,p_inel,ptot,nv,ekin,omega,theta,costheta,q,e0,ekprime;
  
  np = nparticle;  
  for(i=0; (i<nparticle); i++) {
    /* Check cross sections, assume same cross sections for holes
     * as for electrons, for elastic scattering
     */
    nv   = norm(v[i]);
    ekin = eparticle[i];
    p_el = cross_el(ekin)*nv*dt;
   
    /* Only electrons can scatter inelasticlly */
    if (charge[i] < 0)
      p_inel = cross_inel(ekin)*nv*dt;
    else
      p_inel = 0;
      
    /* Test whether we have to scatter at all */
    ptot = (1 - (1-p_el)*(1-p_inel));
    if (debug)
      fprintf(debug,"p_el = %10.3e  p_inel = %10.3e ptot = %10.3e\n",
	      p_el,p_inel,ptot);
    if (rando(seed) < ptot) {
      /* Test whether we have to scatter inelastic */
      ptot = p_inel/(p_el+p_inel);
      if (rando(seed) < ptot) {
	add_scatter_event(&(s[i]),x[i],TRUE,dt*nstep,ekin);
	/* Energy loss in inelastic collision is omega */
	omega = get_omega(ekin,seed,debug);
	/* Scattering angle depends on energy and energy loss */
	q = get_q_inel(ekin,omega,seed,debug);
	costheta = 0.5*(q+omega-2*ekin)/sqrt(ekin*(ekin-omega));

	/* See whether we have gained enough energy to liberate another 
	 * hole-electron pair
	 */
	e0 = band_ener(seed,debug);
	ekprime = e0 + omega - (EFERMI+0.5*EBAND);
	/* Ouput */
	fprintf(fp,"REMARK Inelastic %d: omega = %.2f q = %.2f costheta = %.3f E0 = %.1f\n",
		i+1,omega,q,costheta,e0);
	if ((costheta < -1) || (costheta > 1)) {
	  fprintf(fp,"REMARK Electron/hole creation not possible due to momentum constraints\n");
	  ener[eLATTICE] += omega*ELECTRONVOLT;
	}
	else {
	  theta = acos(costheta);
	  
	  /* Rotate around theta with random delta phi */
	  rotate_theta(v[i],nv,theta,seed,debug);
	  /* Scale the velocity according to the energy loss */
	  svmul(sqrt(1-omega/ekin),v[i],v[i]);
	  
	  if (ekprime > 0) {
	    if (np >= MAXPARTICLE-2)
	      fatal_error(0,"Increase -maxparticle flag to more than %d",
			  MAXPARTICLE);
	    np = create_pair(np,x,v,vold,mass,charge,x[i],
			     ekprime*ELECTRONVOLT,seed);
	    ener[eLATTICE] += (omega-ekprime)*ELECTRONVOLT;
	    (*nhole)++;
	    (*nelec)++;
	  }
	  else
	    ener[eLATTICE] += omega*ELECTRONVOLT;
	}
      }
      else {
	add_scatter_event(&(s[i]),x[i],FALSE,dt*nstep,ekin);
	if (debug)
	  fprintf(debug,"REMARK Elastic scattering event\n");
	
	/* Scattering angle depends on energy only */
	theta = get_theta_el(ekin,seed,debug);
	/* Rotate around theta with random delta phi */
	rotate_theta(v[i],nv,theta,seed,debug);
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

static void do_sim(char *fn,int *seed,int maxstep,real dt,
		   int nsave,bool bForce,bool bScatter,
		   int *nelec,int *nhole,int nana,t_ana_struct *total,
		   t_histo *hmfp)
{
  FILE         *fp;
  int          nparticle[2];
  rvec         *x,*v[2],*f,center;
  real         *charge,*mass,*ener,*eparticle;
  t_ana_struct *ana_struct;
  t_ana_scat   *ana_scat;
  int          step,i,cur = 0;
#define next (1-cur)

  /* Open output file */
  fp = ffopen(fn,"w");
  fprintf(fp,"REMARK Welcome to the electron-hole simulation!\n");
  fprintf(fp,"REMARK The energies printed in this file are in eV\n");
  fprintf(fp,"REMARK Coordinates are in nm because of fixed width format\n");
  fprintf(fp,"REMARK Atomtypes are used for coloring in rasmol\n");
  fprintf(fp,"REMARK O: electrons (red), N: holes (blue)\n");
  fprintf(fp,"REMARK Parametes for this simulation\n");
  fprintf(fp,"REMARK seed = %d maxstep = %d dt = %g\n",*seed,maxstep,dt);
  fprintf(fp,"REMARK nsave = %d Force = %s Scatter = %s\n",
	  nsave,bool_names[bForce],bool_names[bScatter]);
  if (bForce)
    fprintf(fp,"REMARK Force constant for repulsion Alj = %g\n",Alj);

  ana_struct = init_ana_struct(maxstep,nana,dt);
  /* Initiate arrays. The charge array determines whether a particle is 
   * a hole (+1) or an electron (-1)
   */
  snew(x,MAXPARTICLE);          /* Position  */
  snew(v[0],MAXPARTICLE);       /* Velocity  */
  snew(v[1],MAXPARTICLE);       /* Velocity  */
  snew(f,MAXPARTICLE);          /* Force     */
  snew(charge,MAXPARTICLE);     /* Charge    */
  snew(mass,MAXPARTICLE);       /* Mass      */
  snew(eparticle,MAXPARTICLE);  /* Energy per particle */
  snew(ana_scat,MAXPARTICLE);   /* Scattering event statistics */
  snew(ener,eNR);               /* Eenergies */
  
  clear_rvec(center);
  /* Use first atom as center, it has coordinate 0,0,0 */
  if (bScatter) {
    /* Start with a single Auger electron */
    nparticle[cur]  = create_auger(0,x,v[cur],v[next],mass,charge,center,
				   AUGER_ENERGY*ELECTRONVOLT,seed);
  }
  else if (bForce) {
    /* Start with two electron and hole pairs */
    nparticle[cur]  = create_pair(0,x,v[cur],v[next],mass,charge,center,
				  0.2*AUGER_ENERGY*ELECTRONVOLT,seed);
    center[ZZ] = 0.5; /* nm */
    (*nelec)++;
    (*nhole)++;
  }
  else {
    fprintf(fp,"Nothing to do. Doei.\n");
    return;
  }
  nparticle[next] = nparticle[cur];
  for(step=0; (step<=maxstep); step++) {
    if (bScatter)
      nparticle[next] = scatter_all(fp,nparticle[cur],x,v[cur],v[next],
				    mass,charge,ener,eparticle,seed,dt,step,
				    nelec,nhole,ana_scat);
    
    if (bForce)
      calc_forces(nparticle[cur],x,f,charge,ener);
    
    integrate_velocities(nparticle[next],v[cur],v[next],f,mass,dt);
    
    calc_ekin(nparticle[next],v[cur],v[next],charge,mass,ener,eparticle);
    ener[eTOT] = ener[eKIN] + ener[ePOT];
    
    /* Produce output whenever the user says so, or when new
     * particles have been created.
     */
    if ((nparticle[next] != nparticle[cur]) || 
	(step == maxstep) ||
	((nsave != 0) && ((step % nsave) == 0))) {
      write_data(fp,step,(step*dt),nparticle[next],
		 x,v[next],f,charge,ener,eparticle);
    }
    if ((nana != 0) && ((step % nana) == 0))
      analyse_structure(ana_struct,(step*dt),center,x,nparticle[next],charge);
    cur = next;
        
    integrate_positions(nparticle[cur],x,v[cur],dt);
  }
  for(i=0; (i<nparticle[cur]); i++) {
    analyse_scatter(&(ana_scat[i]),hmfp);
    done_scatter(&(ana_scat[i]));
  }
  sfree(ana_scat); 
  add_ana_struct(total,ana_struct);
  done_ana_struct(ana_struct);
  sfree(ana_struct);
  sfree(x); 
  sfree(v[0]); 
  sfree(v[1]);      
  sfree(f);
  sfree(charge); 
  sfree(mass);    
  sfree(eparticle); 
  sfree(ener);
  fclose(fp);
}

void do_sims(char *ptr,char *ndist,char *rmax,char *nion,char *gyr,char *mfp,
	     int *seed,int maxstep,real dt,int nsave,bool bForce,
	     bool bScatter,int nsim)
{
  t_ana_struct *total;
  t_histo      *helec,*hmfp;
  int          *nelectron;
  int          i,imax,ne,nh,nana;
  real         aver;
  FILE         *fp;
  char         *buf;

  nana  = 1;
  total = init_ana_struct(maxstep,nana,dt);
  hmfp  = init_histo(250,0,250);
  helec = init_histo(50,0,50);
  
  snew(buf,strlen(ptr)+10);
  for(i=0; (i<nsim); i++) {
    nh = ne = 0;
    sprintf(buf,"%s-%d.pdb",ptr,i+1);
    do_sim(buf,seed,maxstep,dt,nsave,bForce,bScatter,&ne,&nh,nana,total,hmfp);
    add_histo(helec,ne,1);
  }
  sfree(buf);
  dump_ana_struct(rmax,nion,gyr,total,nsim);
  done_ana_struct(total);

  dump_histo(helec,ndist,"Number of cascade electrons","N","",enormFAC,1.0/nsim);
  dump_histo(hmfp,mfp,"Mean Free Path","Ekin (eV)","MFP (nm)",enormNP,1.0);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "ehole performs a molecular dynamics simulation of electrons and holes"
  };
  static int  seed     = 1993;
  static int  maxstep  = 10000;
  static int  nsave    = 1;
  static int  nsim     = 1;
  static real dt       = 1e-5; /* ps */
  static bool bTest    = FALSE;
  static bool bForce   = TRUE;
  static bool bScatter = TRUE;
  static bool bSMP     = FALSE;
  t_pargs pa[] = {
    { "-seed",      FALSE, etINT,  {&seed}, 
      "Random seed" },
    { "-maxstep",   FALSE, etINT,  {&maxstep}, 
      "Number of integration steps" },
    { "-nsave",     FALSE, etINT,  {&nsave}, 
      "Number of steps after which to save output. 0 means only when particles created. Final step is always written." },
    { "-dt",        FALSE, etREAL, {&dt}, 
      "Integration time step (ps)" },
    { "-fermi",     FALSE, etREAL, {&EFERMI}, 
      "Fermi energy (eV) of the sample. Default is for Diamond" },
    { "-gap",       FALSE, etREAL, {&EBAND}, 
      "Band gap energy (eV) of the sample. Default is for Diamond" },
    { "-auger",     FALSE, etREAL, {&AUGER_ENERGY}, 
      "Impact energy (eV) of first electron" },
    { "-dx",        FALSE, etREAL, {&DELTAX},
      "Distance between electron and hole when creating a pair" },
    { "-fermi",     FALSE, etREAL, {&EFERMI}, 
      "Fermi energy of the sample. Default is for Diamond" },
    { "-test",      FALSE, etBOOL, {&bTest},
      "Test table aspects of the program rather than running it for real" },
    { "-force",     FALSE, etBOOL, {&bForce},
      "Apply Coulomb/Repulsion forces" },
    { "-scatter",   FALSE, etBOOL, {&bScatter},
      "Do the scattering events" },
    { "-maxparticle", FALSE, etINT, {&MAXPARTICLE},
      "Maximum number of particles" },
    { "-nsim",      FALSE, etINT,  {&nsim},
      "Number of independent simulations writing to different output files" }
  };
#define NPA asize(pa)
  t_filenm fnm[] = {
    { efPDB, "-o", "ehole",  ffWRITE },
    { efXVG, "-h", "histo",  ffWRITE },
    { efXVG, "-r", "radius", ffWRITE },
    { efXVG, "-g", "gyrate", ffWRITE },
    { efXVG, "-m", "mfp",    ffWRITE },
    { efXVG, "-n", "nion",   ffWRITE }
  };
#define NFILE asize(fnm)
  char *buf,*ptr,*rptr;
  int  i;
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);
  if (DELTAX <= 0)
    fatal_error(0,"DELTAX should be > 0");
  Alj = FACEL*pow(DELTAX,5);
  
  if (bTest) 
    test_tables(&seed,ftp2fn(efPDB,NFILE,fnm));  
  else {
    ptr  = ftp2fn(efPDB,NFILE,fnm);
    snew(buf,strlen(ptr)+10);
    rptr = strdup(ptr);
    if ((ptr  = strstr(rptr,".pdb")) != NULL)
      *ptr = '\0';
    do_sims(rptr,opt2fn("-h",NFILE,fnm),opt2fn("-r",NFILE,fnm),
	    opt2fn("-n",NFILE,fnm),opt2fn("-g",NFILE,fnm),opt2fn("-m",NFILE,fnm),
	    &seed,maxstep,dt,nsave,bForce,bScatter,nsim);
  }
  thanx(stdout);
  
  return 0;
}
