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
static char *SRCID_g_dipoles_c = "$Id$";

#include <string.h>
#include <math.h>
#include "macros.h"
#include "statutil.h"
#include "sysstuff.h"
#include "smalloc.h"
#include "vec.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "xvgr.h"
#include "txtdump.h"
#include "gstat.h"
#include "rdgroup.h"
#include "random.h"
#include "names.h"
#include "physics.h"
#include "calcmu.h"
#include "enxio.h"
#include "nrjac.h"

#define e2d(x) ENM2DEBYE*(x)
#define EANG2CM  E_CHARGE*1.0e-10       /* e Angstrom to Coulomb meter */
#define CM2D  SPEED_OF_LIGHT*1.0e+24    /* Coulomb meter to Debye */

typedef struct {
  int  nelem;
  real spacing,radius;
  real *elem;
  int  *count;
} t_gkrbin;

static t_gkrbin *mk_gkrbin(real radius)
{
  t_gkrbin *gb;
  char *ptr;

  snew(gb,1);
  
  if ((ptr = getenv("GKRWIDTH")) != NULL) {
    double bw;

    sscanf(ptr,"%lf",&bw);
    gb->spacing = bw; 
  }
  else
    gb->spacing = 0.01; /* nm */
  gb->nelem   = 1 + radius/gb->spacing;
  gb->radius  = radius;
  snew(gb->elem,gb->nelem);
  snew(gb->count,gb->nelem);
  
  return gb;
}

static void done_gkrbin(t_gkrbin **gb)
{
  sfree((*gb)->elem);
  sfree((*gb)->count);
  sfree((*gb));
  *gb = NULL;
}

void do_gkr(t_gkrbin *gb,int ngrp,atom_id grpindex[],
	    atom_id mindex[],atom_id ma[],rvec x[],rvec mu[],
	    matrix box,t_atom *atom,int nAtom)
{
  static rvec *xcm=NULL;
  int  gi,aj,j0,j1,i,j,k,index;
  real qtot,q,r2;
  rvec dx;
  
  if (!xcm)
    snew(xcm,ngrp);
    
  for(i=0; (i<ngrp); i++) {
    /* Calculate center of mass of molecule */
    gi = grpindex ? grpindex[i] : i;
    j0 = mindex[gi];
    
    if (nAtom > 0)
      copy_rvec(x[ma[j0+nAtom-1]],xcm[i]);
    else {
      j1 = mindex[gi+1];
      clear_rvec(xcm[i]);
      qtot = 0;
      for(j=j0; j<j1; j++) {
	aj = ma[j];
	q = fabs(atom[aj].q);
	qtot += q;
	for(k=0; k<DIM; k++)
	  xcm[i][k] += q*x[aj][k];
      }
      svmul(1/qtot,xcm[i],xcm[i]);
    }
  }
  
  for(i=0; i<ngrp; i++) {
    for(j=i+1; j<ngrp; j++) {
      /* Compute distance between molecules including PBC */
      pbc_dx(xcm[i],xcm[j],dx);
      index = (int)(norm(dx)/gb->spacing);
      gb->elem[index]  += cos_angle(mu[i],mu[j]);
      gb->count[index] ++;
    }
  }
}

void print_gkrbin(char *fn,t_gkrbin *gb,
		  int ngrp,int nframes,real volume)
{
  /* We compute Gk(r), gOO and hOO according to
   * Nymand & Linse, JCP 112 (2000) pp 6386-6395.
   * In this implementation the angle between dipoles is stored
   * rather than their inner product. This allows to take polarizible
   * models into account. The RDF is calculated as well, almost for free!
   */
  FILE   *fp;
  char   *leg[] = { "G\\sk\\N(r)", "< cos >", "h\\sOO\\N", "g\\sOO\\N" };
  int    i,last;
  real   x0,x1,ggg,Gkr,vol_s,rho,gOO,hOO,cosav;
  double fac;
    
  fp=xvgropen(fn,"Distance dependent Gk","r (nm)","G\\sk\\N(r)");
  xvgr_legend(fp,asize(leg),leg);
  
  Gkr = 1;  /* Self-dipole inproduct = 1 */
  rho = ngrp/volume;
  
  if (debug) {
    fprintf(debug,"Number density is %g molecules / nm^3\n",rho);
    fprintf(debug,"ngrp = %d, nframes = %d\n",ngrp,nframes);
  }
  
  last = gb->nelem-1;
  while(last>1 && gb->elem[last-1]==0)
    last--;

  /* Divide by dipole squared, by number of frames, by number of origins.
   * Multiply by 2 because we only take half the matrix of interactions
   * into account.
   */
  fac  = 2.0/((double) ngrp * (double) nframes);

  x0 = 0;
  for(i=0; i<last; i++) {
    /* Centre of the coordinate in the spherical layer */
    x1    = x0+gb->spacing;
    
    /* Volume of the layer */
    vol_s = (4.0/3.0)*M_PI*(x1*x1*x1-x0*x0*x0);
    
    /* gOO */
    gOO   = gb->count[i]*fac/(rho*vol_s);
    
    /* Dipole correlation hOO, normalized by the relative number density, like
     * in a Radial distribution function.
     */
    ggg  = gb->elem[i]*fac;
    hOO  = 3.0*ggg/(rho*vol_s);
    Gkr += ggg;
    if (gb->count[i])
      cosav = gb->elem[i]/gb->count[i];
    else
      cosav = 0;
    
    fprintf(fp,"%10.5e %12.5e %12.5e %12.5e %12.5e\n",x1,Gkr,cosav,hOO,gOO);
    
    /* Swap x0 and x1 */
    x0 = x1;
  }
  ffclose(fp);
}

bool read_mu_from_enx(int fmu,int Vol,ivec iMu,rvec mu,real *vol,
		      real *t,int step,int nre)
{
  int      i,ndr;
  bool     eof;
  t_energy *ee;  

  snew(ee,nre);

  eof = do_enx(fmu,t,&step,&nre,ee,&ndr,NULL);

  if (eof) {
    if (Vol != -1)          /* we've got Volume in the energy file */
      *vol = ee[Vol].e;
    for (i=0;(i<DIM);i++)
      mu[i]=ee[iMu[i]].e;
  }
 
  free(ee);
 
  return eof;
}


void mol_dip(int k0,int k1,atom_id ma[],rvec x[],t_atom atom[],rvec mu)
{
  int  k,kk,m;
  real q;
  
  clear_rvec(mu);
  for(k=k0; (k<k1); k++) {
    kk = ma[k];
    q  = e2d(atom[kk].q);
    for(m=0; (m<DIM); m++) 
      mu[m] += q*x[kk][m];
  }
}

#define NDIM 3          /* We will be using a numerical recipes routine */

void mol_quad(int k0,int k1,atom_id ma[],rvec x[],t_atom atom[],rvec quad)
{
  int  i,k,kk,m,n,niter;
  real q,r2,mass,masstot;
  rvec com;          /* center of mass */
  rvec r;            /* distance of atoms to center of mass */
  rvec xmass;
  real rcom_m,rcom_n;
  tensor quadrupole;
  double **inten;
  double dd[NDIM],**ev,tmp;

  snew(inten,NDIM);
  snew(ev,NDIM);
  for(i=0; (i<NDIM); i++) {
    snew(inten[i],NDIM);
    snew(ev[i],NDIM);
    dd[i]=0.0;
  }

  for(i=0; (i<NDIM); i++)
    for(m=0; (m<NDIM); m++)
      inten[i][m]=0;
  
  clear_rvec(quad);
  clear_rvec(com);
  clear_mat(quadrupole);

  clear_rvec(xmass);

  masstot = 0.0;

  for(k=k0; (k<k1); k++) {
    kk = ma[k];
    mass = atom[kk].m;
    masstot += mass;
    /*    svmul(mass,x[k],xmass[k]); */
    for(i=0; (i<DIM); i++) {
      xmass[i] = mass*x[k][i];
    }
    rvec_inc(com,xmass);
  }
  svmul((1.0/masstot),com,com);

  /* We want traceless quadrupole moments, so let us calculate the complete
   * quadrupole moment tensor and diagonalize this tensor to get
   * the individual components on the diagonal.
   */

#define delta(a,b) (( a == b ) ? 1.0 : 0.0)

  for(k=k0; (k<k1); k++) {       /* loop over atoms in a molecule */
    kk = ma[k];
    q  = (atom[kk].q)*100.0;
    rvec_sub(x[k],com,r);
    r2 = iprod(r,r);
    for(m=0; (m<DIM); m++) {
      for(n=0; (n<DIM); n++) {  
	quadrupole[m][n]+=0.5*q*(3.0*r[m]*r[n] - r2*delta(m,n))*EANG2CM*CM2D;
	inten[m][n]=quadrupole[m][n];
      }
    }
  }
#ifdef DEBUG
  pr_rvecs(stdout,0,"Quadrupole",quadrupole,DIM);
#endif

  /* We've got the quadrupole tensor, now diagonalize the sucker */
  
  jacobi(inten,3,dd,ev,&niter);

  /* Sort the eigenvalues, for water we know that the order is as follows:
   *
   * Q_yy, Q_zz, Q_xx
   *
   * At the moment I have no idea how this will work out for other molecules...
   */

#define SWAP(i) 			\
  if (dd[i+1] > dd[i]) {	        \
    tmp=dd[i];		                \
    dd[i]=dd[i+1];			\
    dd[i+1]=tmp;			\
  }
  SWAP(0);
  SWAP(1);
  SWAP(0);

  quad[0]=dd[2];  /* yy */
  quad[1]=dd[0];  /* zz */
  quad[2]=dd[1];  /* xx */

#ifdef DEBUG
  pr_rvec(stdout,0,"Quadrupole",quad,DIM);
#endif

  /* clean-up */
  for(i=0; (i<NDIM); i++) {
    sfree(inten[i]);
    sfree(ev[i]);
  }
  sfree(inten);
  sfree(ev);

  /*  sfree(xmass); */
}

/*
 * Calculates epsilon according to M. Neumann, Mol. Phys. 50, 841 (1983)
 */ 
real calc_eps(real M_diff,real volume,real epsRF,real temp)
{
  double eps,A,teller,noemer;
  double eps_0=8.854187817e-12;     /* epsilon_0 in C^2 J^-1 m^-1 */
  double fac=1.112650021e-59;       /* converts Debye^2 to C^2 m^2 */

  A = M_diff*fac/(3*eps_0*volume*NANO*NANO*NANO*BOLTZMANN*temp);
 
  if (epsRF == 0.0) {
    teller = 1 + A;
    noemer = 1; 
  } else { 
    teller = 1 + (A*2*epsRF/(2*epsRF+1));
    noemer = 1 - (A/(2*epsRF+1));
  }
  eps = teller / noemer;

  return eps;
}


static void do_dip(char *fn,char *topf,
		   char *out_mtot,char *out_eps,char *out_aver, 
		   char *dipdist,bool bAverCorr,
		   bool bCorr,   char *corf,
		   bool bGkr,    char *gkrfn,
		   bool bQuad,   char *quadfn,
		   bool bMU,     char *mufn,
		   int gnx,atom_id grpindex[],
		   real mu_max,real mu_aver,
		   real epsilonRF,real temp,
		   int gkatom,int skip)
{
  static char *leg_mtot[] = { 
    "< M\\sx \\N>", 
    "< M\\sy \\N>",
    "< M\\sz \\N>",
    "< |M\\stot \\N| >"
  };
#define NLEGMTOT asize(leg_mtot)
  static char *leg_eps[] = { 
    "epsilon",
    "G\\sk",
    "g\\sk"
  };
#define NLEGEPS asize(leg_eps)
  static char *leg_aver[] = { 
    "< |M|\\S2\\N >", 
    "< |M| >\\S2\\N",
    "< |M|\\S2\\N > - < |M| >\\S2\\N",
    "< |M| >\\S2\\N / < |M|\\S2\\N >"
  };
#define NLEGAVER asize(leg_aver)

  FILE  *outdd,*outmtot,*outaver,*outeps;
  rvec       *x,*dipole=NULL,mu_t,M_av,M_av2,Q_av,Q_av2,*quadrupole=NULL;
  t_gkrbin   *gkrbin;
  int        nframes=1000,fmu=0,nre,timecheck=0;
  char       **enm=NULL;
  real       rcut=0;
  matrix     box;
  bool       bCont;
  real       t,t0,t1,dt;
  double     M_diff=0,epsilon;
  double     mu_ave,mu_mol,M2_ave=0,M_ave2=0;
  rvec       quad_ave,quad_mol;
  ivec       iMu;
  int        iVol;
  real       M_XX,M_YY,M_ZZ,M_XX2,M_YY2,M_ZZ2,Gk=0,g_k=0;
  real       **muall=NULL,**quadall=NULL;
  t_topology *top;
  t_atom     *atom=NULL;
  t_block    *mols=NULL;
  int        i,j,k,m,natom=0,nmol,status,teller,tel3;
  int        *dipole_bin,ndipbin,ibin;
  real       volume,vol_aver;
  unsigned long mode;
  /* PvM, I need these to be able to call read_tpx */
  int        natoms,step;
  real       lambda;

  snew(top,1);
  read_tpx(topf,&step,&t,&lambda,NULL,box,
	   &natoms,NULL,NULL,NULL,top);
  volume   = det(box);
  vol_aver = 0.0;
  
  if (!grpindex)
    gnx = top->blocks[ebMOLS].nr;
    
  iVol=-1;
  if (bMU) {
    fmu = open_enx(mufn,"r");
    do_enxnms(fmu,&nre,&enm);

    /* Determine the indexes of the energy grps we need */
    for (i=0; (i<nre); i++) {
      if (strstr(enm[i],"Volume"))
	iVol=i;
      else if (strstr(enm[i],"Mu-X"))
	iMu[XX]=i;
      else if (strstr(enm[i],"Mu-Y"))
	iMu[YY]=i;
      else if (strstr(enm[i],"Mu-Z"))
	iMu[ZZ]=i;
    }
  }
  else {
    atom = top->atoms.atom;
    mols = &(top->blocks[ebMOLS]);
  }
  
  if (iVol == -1)
    printf("Using Volume from topology: %g nm^3\n",volume);

  /* Correlation stuff */ 
  if (bCorr) {
    if (bAverCorr) {
      snew(muall,1);
      snew(muall[0],nframes*DIM);
      if (bQuad) {
	snew(quadall,1);
	snew(quadall[0],nframes*DIM);
      }
    }
    else {
      snew(muall,gnx);
      for(i=0; (i<gnx); i++)
	snew(muall[i],nframes*DIM);
      if (bQuad) {
	snew(quadall,gnx);
	for(i=0; (i<gnx); i++)
	  snew(quadall[i],nframes*DIM);
      }
    }
  }

  /* Allocate array which contains for every molecule in a frame the
   * dipole moment.
   */
  if (!bMU)
    snew(dipole,gnx);
  if (bQuad)
    snew(quadrupole,gnx);
    
  /* Open all the files */
  outmtot = xvgropen(out_mtot,
		     "Total dipole moment of the simulation box vs. time",
		     "Time (ps)","Total Dipole Moment (Debye)");
  outeps  = xvgropen(out_eps,"Epsilon and Kirkwood factors",
		     "Time (ps)","");
  outaver = xvgropen(out_aver,"Total dipole moment",
		     "Time (ps)","D");
		     
  /* Write legends to all the files */
  xvgr_legend(outmtot,NLEGMTOT,leg_mtot);
  xvgr_legend(outaver,NLEGAVER,leg_aver);
  
  if (bMU && (mu_aver == -1))
    xvgr_legend(outeps,NLEGEPS-2,leg_eps);
  else
    xvgr_legend(outeps,NLEGEPS,leg_eps);
    
  teller = 0;
  /* Read the first frame from energy or traj file */
  if (bMU)
    do {
      bCont = read_mu_from_enx(fmu,iVol,iMu,mu_t,&volume,&t,teller,nre); 
      if (bCont) {  
	timecheck=check_times(t,t);
	if (timecheck < 0)
	  teller++;
	if ((teller % 10) == 0)
	  fprintf(stderr,"\r Skipping Frame %6d, time: %8.3f", teller, t);
      }
      else {
	printf("End of %s reached\n",mufn);
	break;
      }
    } while (bCont && (timecheck < 0));
  else
    natom  = read_first_x(&status,fn,&t,&x,box);
  
  /* Calculate spacing for dipole bin (simple histogram) */
  ndipbin = 1+(mu_max/0.01);
  snew(dipole_bin, ndipbin);
  epsilon    = 0.0;
  M_XX=M_XX2 = 0.0;
  M_YY=M_YY2 = 0.0;
  M_ZZ=M_ZZ2 = 0.0;
  mu_ave     = 0.0;

  if (bQuad)
    clear_rvec(quad_ave);

  if (bGkr) {
    /* Use 0.7 iso 0.5 to account for pressure scaling */
    rcut   = 0.7*sqrt(sqr(box[XX][XX])+sqr(box[YY][YY])+sqr(box[ZZ][ZZ]));
    gkrbin = mk_gkrbin(rcut); 
  }

  /* Start while loop over frames */
  t1 = t0 = t;
  teller = 0;
  do {
    if (bCorr && (teller >= nframes)) {
      nframes += 1000;
      if (bAverCorr) {
	srenew(muall[0],nframes*DIM);
	if (bQuad) 
	  srenew(quadall[0],nframes*DIM);
      }
      else {
	for(i=0; (i<gnx); i++)
	  srenew(muall[i],nframes*DIM);
	if (bQuad) {
	  for(i=0; (i<gnx); i++)
	    srenew(quadall[i],nframes*DIM);
	}
      }
    }
    t1 = t;

    /* Initialise */
    clear_rvec(M_av);
    clear_rvec(M_av2);

    if (bMU) {
      if (timecheck == 0) {
	for(m=0; (m<DIM); m++) {
	  M_av[m]  += mu_t[m];          /* M per frame */
	  M_av2[m] += mu_t[m]*mu_t[m];  /* M^2 per frame */
	}
      }
    } 
    else {
      /* Begin loop of all molecules in frame */
      for(i=0; (i<gnx); i++) {
	int gi = grpindex ? grpindex[i] : i;
	mol_dip(mols->index[gi],mols->index[gi+1],mols->a,x,atom,dipole[i]);
	if (bQuad)
	  mol_quad(mols->index[gi],mols->index[gi+1],
		   mols->a,x,atom,quadrupole[i]);
	
	if (bCorr && !bAverCorr) {
	  tel3=DIM*teller;
	  muall[i][tel3+XX] = dipole[i][XX];
	  muall[i][tel3+YY] = dipole[i][YY];
	  muall[i][tel3+ZZ] = dipole[i][ZZ];
	  if (bQuad) {
	    quadall[i][tel3+XX] = quadrupole[i][XX];
	    quadall[i][tel3+YY] = quadrupole[i][YY];
	    quadall[i][tel3+ZZ] = quadrupole[i][ZZ];
	  }
	}
	mu_mol = 0.0;
	for(m=0; (m<DIM); m++) {
	  M_av[m]  += dipole[i][m];               /* M per frame */
	  mu_mol   += dipole[i][m]*dipole[i][m];  /* calc. mu for distribution */
	}

	mu_ave += sqrt(mu_mol);                   /* calc. the average mu */
	if (bQuad) {
	  clear_rvec(quad_mol);
	  for(m=0; (m<DIM); m++) {
	    Q_av[m]  += quadrupole[i][m];                    /* Q per frame */
	    quad_mol[m]  += quadrupole[i][m]*quadrupole[i][m];  
	    quad_ave[m] += quadrupole[i][m];
	  }
	}
	
	/* Update the dipole distribution */
	ibin = (ndipbin*sqrt(mu_mol)/mu_max);
	if (ibin < ndipbin)
	  dipole_bin[ibin]++;
      } /* End loop of all molecules in frame */
      
      /* Compute square of total dipole an quadrupole */
      for(m=0; (m<DIM); m++) {
	M_av2[m] = sqr(M_av[m]);
	if (bQuad)
	  Q_av2[m] = sqr(Q_av[m]);
      }
    }    
    if (bGkr) {
      init_pbc(box,FALSE);
      do_gkr(gkrbin,gnx,grpindex,mols->index,mols->a,x,dipole,box,
	     atom,gkatom);
    }
    
    if (bAverCorr) {
      tel3=DIM*teller;
      muall[0][tel3+XX] = M_av[XX];
      muall[0][tel3+YY] = M_av[YY];
      muall[0][tel3+ZZ] = M_av[ZZ];
      if (bQuad) {
	quadall[0][tel3+XX] = Q_av[XX];
	quadall[0][tel3+YY] = Q_av[YY];
	quadall[0][tel3+ZZ] = Q_av[ZZ];
      }
    }

    /* Write to file the total dipole moment of the box, and its components 
     * for this frame.
     */
    if ((skip == 0) || ((teller % skip) == 0))
      fprintf(outmtot,"%10g  %12.8e %12.8e %12.8e %12.8e\n",
	      t,M_av[XX],M_av[YY],M_av[ZZ],norm(M_av));

    M_XX      += M_av[XX];
    M_XX2     += M_av2[XX];
    M_YY      += M_av[YY];
    M_YY2     += M_av2[YY];
    M_ZZ      += M_av[ZZ];
    M_ZZ2     += M_av2[ZZ];

    /* Increment loop counter */
    teller++;
    
    /* Calculate for output the running averages */
    M2_ave  = (M_XX2+M_YY2+M_ZZ2)/teller;
    
    /* Strange construction because teller*teller may go beyond INT_MAX */
    M_ave2  = (sqr(M_XX/teller)+sqr(M_YY/teller)+sqr(M_ZZ/teller));
    M_diff  = M2_ave - M_ave2;

    /* Compute volume from box in traj, else we use the one from above */
    if (!bMU)
      volume  = det(box);
    vol_aver += volume;
    
    epsilon = calc_eps(M_diff,volume,epsilonRF,temp);

    /* Calculate running average for dipole */
    if (mu_ave != 0) 
      mu_aver = (mu_ave/gnx)/teller;
    
    if ((skip == 0) || ((teller % skip) == 0)) {
      /* Write to file < |M|^2 >, < |M| >^2. And the difference between 
       * the two. Here M is sum mu_i. Further write the finite system
       * Kirkwood G factor and epsilon.
       */
      fprintf(outaver,"%10g  %10.3e %10.3e %10.3e %10.3e\n",
	      t,M2_ave,M_ave2,M_diff,M_ave2/M2_ave);
	      
      if (!bMU || (mu_aver != -1)) {
	/* Finite system Kirkwood G-factor */
	Gk = M_diff/(gnx*mu_aver*mu_aver);
	/* Infinite system Kirkwood G-factor */
	if (epsilonRF == 0.0) 
	  g_k = ((2*epsilon+1)*Gk/(3*epsilon));
	else 
	  g_k = ((2*epsilonRF+epsilon)*(2*epsilon+1)*
		 Gk/(3*epsilon*(2*epsilonRF+1)));
	
	fprintf(outeps,"%10g  %10.3e %10.3e %10.3e\n",t,epsilon,Gk,g_k);
      }
      else 
	fprintf(outeps,"%10g  %12.8e\n",t,epsilon);
    }
    
    if (bMU)
      bCont = read_mu_from_enx(fmu,iVol,iMu,mu_t,&volume,&t,teller,nre); 
    else
      bCont = read_next_x(status,&t,natom,x,box);
  } while(bCont);
  
  if (!bMU)
    close_trj(status);
    
  fclose(outmtot);
  fclose(outaver);
  fclose(outeps);

  vol_aver /= teller;
  printf("Average volume over run is %g\n",vol_aver);
  if (bGkr) 
    print_gkrbin(gkrfn,gkrbin,gnx,teller,vol_aver);

  /* Autocorrelation function */  
  if (bCorr) {
    if (teller < 2) {
      printf("Not enough frames for autocorrelation\n");
    }
    else {
      dt=(t1 - t0)/(teller-1);
      printf("t0 %g, t %g, teller %d\n", t0,t,teller);
      
      mode = eacVector;

      if (bAverCorr)
	do_autocorr(corf,"Autocorrelation Function of Total Dipole",
		    teller,1,muall,dt,mode,TRUE);
      else
	do_autocorr(corf,"Dipole Autocorrelation Function",
		    teller,gnx,muall,dt,mode,TRUE);
    }
  }
  if (!bMU) {
    printf("\n\nAverage dipole moment (Debye)\n");
    printf(" Tot= %g\n",  (mu_ave/gnx)/teller);
    if (bQuad) {
      printf("Average quadrupole moment (Debye-Ang)\n");
      printf(" XX=  %g  YY=  %g ZZ=  %g norm= %g asymm= %g\n\n",  
	      quad_ave[XX]/(gnx*teller),
	      quad_ave[YY]/(gnx*teller),
	      quad_ave[ZZ]/(gnx*teller),
	      norm(quad_ave)/(gnx*teller),
	      (quad_ave[ZZ] - quad_ave[XX])/ quad_ave[YY]);
    }
  }
  printf("The following averages for the complete trajectory have been calculated:\n\n");
  printf(" Total < M_x > = %g Debye\n", M_XX/teller);
  printf(" Total < M_y > = %g Debye\n", M_YY/teller);
  printf(" Total < M_z > = %g Debye\n\n", M_ZZ/teller);

  printf(" Total < M_x^2 > = %g Debye^2\n", M_XX2/teller);
  printf(" Total < M_y^2 > = %g Debye^2\n", M_YY2/teller);
  printf(" Total < M_z^2 > = %g Debye^2\n\n", M_ZZ2/teller);

  printf(" Total < |M|^2 > = %g Debye^2\n", M2_ave);
  printf(" Total < |M| >^2 = %g Debye^2\n\n", M_ave2);

  printf(" < |M|^2 > - < |M| >^2 = %g Debye^2\n\n", M_diff);
  if (!bMU || (mu_aver != -1)) {
    printf("Finite system Kirkwood g factor G_k = %g\n", Gk);
    printf("Infinite system Kirkwood g factor g_k = %g\n\n", g_k);
  }
  printf("Epsilon = %g\n", epsilon);

  if (!bMU) {
    /* Write to file the dipole moment distibution during the simulation.
     */
    outdd=xvgropen(dipdist,"Dipole Moment Distribution","mu (Debye)","");
    for(i=0; (i<ndipbin); i++)
      fprintf(outdd,"%10g  %d\n",(i*mu_max)/ndipbin,dipole_bin[i]);
    fclose(outdd);
    sfree(dipole_bin);
  }
  if (bGkr) 
    done_gkrbin(&gkrbin);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_dipoles computes the total dipole plus fluctuations of a simulation",
    "system. From this you can compute e.g. the dielectric constant for",
    "low dielectric media[PAR]",
    "The file dip.xvg contains the total dipole moment of a frame, the",
    "components as well as the norm of the vector.",
    "The file aver.xvg contains < |Mu|^2 > and < |Mu| >^2 during the",
    "simulation.",
    "The file dip.xvg contains the distribution of dipole moments during",
    "the simulation",
    "The mu_max is used as the highest value in the distribution graph.[PAR]",
    "Furthermore the dipole autocorrelation function will be computed, when",
    "option -c is used. It can be averaged over all molecules, ",
    "or (with option -avercorr) it can be computed as the autocorrelation",
    "of the total dipole moment of the simulation box.[PAR]",
    "At the moment the dielectric constant is calculated only correct if",
    "a rectangular or cubic simulation box is used.[PAR]",
    "Option [TT]-g[tt] produces a plot of the distance dependent Kirkwood",
    "G-factor, as well as the average cosine of the angle between the dipoles",
    "as a function of the distance. The plot also includes gOO and hOO",
    "according to Nymand & Linse, JCP 112 (2000) pp 6386-6395.[PAR]",
    "[PAR]",
    "EXAMPLES[PAR]",
    "g_dipoles -P1 -n mols -o dip_sqr -mu 2.273 -mumax 5.0",
    "-nofft[PAR]",
    "This will calculate the autocorrelation function of the molecular",
    "dipoles using a first order Legendre polynomial of the angle of the",
    "dipole vector and itself a time t later. For this calculation 1001",
    "frames will be used. Further the dielectric constant will be calculated",
    "using an epsilonRF of infinity (default), temperature of 300 K (default) and",
    "an average dipole moment of the molecule of 2.273 (SPC). For the",
    "distribution function a maximum of 5.0 will be used."
  };
  static real mu_max=5, mu_aver=-1;
  static real epsilonRF=0.0, temp=300;
  static bool bAverCorr=FALSE;
  static int  skip=0,nFA=0;
  t_pargs pa[] = {
    { "-mu",       FALSE, etREAL, {&mu_aver},
      "dipole of a single molecule (in Debye)" },
    { "-mumax",    FALSE, etREAL, {&mu_max},
      "max dipole in Debye (for histrogram)" },
    { "-epsilonRF",    FALSE, etREAL, {&epsilonRF},
      "epsilon of the reaction field used during the simulation, needed for dieclectric constant calculation. WARNING: 0.0 means infinity (default)" },
    { "-skip",    FALSE, etINT, {&skip},
      "Skip steps in the output (but not in the computations)" },
    { "-temp",    FALSE, etREAL, {&temp},
      "average temperature of the simulation (needed for dielectric constant calculation)" },
    { "-avercorr", FALSE, etBOOL, {&bAverCorr},
      "calculate AC function of average dipole moment of the simulation box rather than average of AC function per molecule" },
    { "-gkratom", FALSE, etINT, {&nFA},
      "Use the n-th atom of a molecule (starting from 1) to calculate the distance between molecules rather than the center of charge (when 0) in the calculation of distance dependent Kirkwood factors" }
  };
  int          gnx;
  atom_id      *grpindex;
  char         *grpname;
  bool         bCorr,bQuad,bGkr,bMU;  
  t_filenm fnm[] = {
    { efENX, "-enx", NULL,    ffOPTRD },
    { efTRX, "-f", NULL,      ffREAD },
    { efTPX, NULL, NULL,      ffREAD },
    { efNDX, NULL, NULL,      ffOPTRD },
    { efXVG, "-o", "Mtot",    ffWRITE },
    { efXVG, "-e", "epsilon", ffWRITE },
    { efXVG, "-a", "aver",    ffWRITE },
    { efXVG, "-d", "dipdist", ffWRITE },
    { efXVG, "-c", "dipcorr", ffOPTWR },
    { efXVG, "-g", "gkr",     ffOPTWR },
    { efXVG, "-q", "quadrupole", ffOPTWR },
  };
#define NFILE asize(fnm)
  int     npargs;
  t_pargs *ppa;
  
  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  ppa    = add_acf_pargs(&npargs,pa);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,TRUE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL);

  printf("Using %g as mu_max and %g as the dipole moment.\n", 
	  mu_max,mu_aver);
  if (epsilonRF == 0.0)
    printf("WARNING: EpsilonRF = 0.0, this really means EpsilonRF = infinity\n");

  bMU   = opt2bSet("-enx",NFILE,fnm);
  bQuad = opt2bSet("-q",NFILE,fnm);
  bGkr  = opt2bSet("-g",NFILE,fnm);
  if (bMU) {
    bAverCorr = TRUE;
    if (bQuad) {
      printf("WARNING: Can not determine quadrupoles from energy file\n");
      bQuad = FALSE;
    }
    if (bGkr) {
      printf("WARNING: Can not determine Gk(r) from energy file\n");
      bGkr  = FALSE;
    }
    if (mu_aver == -1) 
      printf("WARNING: Can not calculate Gk and gk, since you did\n"
	     "         not enter a valid dipole for the molecules\n");
  }
  
  if (ftp2bSet(efNDX,NFILE,fnm))
    rd_index(ftp2fn(efNDX,NFILE,fnm),1,&gnx,&grpindex,&grpname);
  else {
    gnx=1;
    grpindex=NULL;
  }
  bCorr   = (bAverCorr || opt2bSet("-c",NFILE,fnm));
  do_dip(ftp2fn(efTRX,NFILE,fnm),ftp2fn(efTPX,NFILE,fnm),
	 opt2fn("-o",NFILE,fnm),opt2fn("-e",NFILE,fnm),
	 opt2fn("-a",NFILE,fnm),opt2fn("-d",NFILE,fnm),
	 bAverCorr,bCorr,
	 opt2fn("-c",NFILE,fnm),
	 bGkr,    opt2fn("-g",NFILE,fnm),
	 bQuad,   opt2fn("-q",NFILE,fnm),
	 bMU,     opt2fn("-enx",NFILE,fnm),
	 gnx,grpindex,mu_max,mu_aver,epsilonRF,temp,nFA,skip);
  
  do_view(opt2fn("-o",NFILE,fnm),"-autoscale xy -nxy");
  do_view(opt2fn("-e",NFILE,fnm),"-autoscale xy -nxy");
  do_view(opt2fn("-a",NFILE,fnm),"-autoscale xy -nxy");
  do_view(opt2fn("-d",NFILE,fnm),"-autoscale xy");
  do_view(opt2fn("-c",NFILE,fnm),"-autoscale xy");

  thanx(stderr);
  
  return 0;
}
