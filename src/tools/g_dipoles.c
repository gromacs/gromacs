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

#define MAXBIN 200 
#define e2d(x) ENM2DEBYE*(x)
#define EANG2CM  E_CHARGE*1.0e-10       /* e Angstrom to Coulomb meter */
#define CM2D  SPEED_OF_LIGHT*1.0e+24    /* Coulomb meter to Debye */

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

void do_gkr(int ngrp,atom_id grpindex[],
	    atom_id mindex[],atom_id ma[],rvec x[],rvec mu[],
	    matrix box,int ngraph,real *graph,int *count,real rcut,
	    int nAtom)
{
  static rvec *xcm=NULL;
  int  gi,aj,j0,j1,i,j,k,index;
  real fac,r2,rc2,rc_2,mu2;
  rvec dx;
  
  if (!xcm)
    snew(xcm,ngrp);
    
  rc2  = rcut*rcut;
  rc_2 = 1.0/rc2;
  for(i=0; (i<ngrp); i++) {
    /* Calculate center of mass of molecule */
    gi = grpindex ? grpindex[i] : i;
    j0 = mindex[gi];
    
    if (nAtom >= 0)
      copy_rvec(x[ma[j0+nAtom]],xcm[i]);
    else {
      j1 = mindex[gi+1];
      clear_rvec(xcm[i]);
      fac=1.0/(j1-j0);
      for(j=j0; (j<j1); j++) {
      aj = ma[j];
      for(k=0; (k<DIM); k++)
	xcm[i][k] += fac*x[aj][k];
      }
    }
  }
  
  for(i=0; (i<ngrp); i++) {
    for(j=i+1; (j<ngrp); j++) {
      /* Compute distance between molecules including PBC */
      pbc_dx(xcm[i],xcm[j],dx);
      r2 = iprod(dx,dx);
      if (r2 < rc2) {
	index = ngraph*sqrt(r2*rc_2);
	if (index < ngraph) {
	  mu2           = iprod(mu[i],mu[j]);
	  graph[index] += mu2;
	  count[index] += 1;
	}
      }
    }
  }
}

void print_gkr(char *fn,int ngraph,real rcut,real graph[],int count[],
	       real mu)
{
  FILE *fp;
  int  i;
  real x,y,ytot;
  
  fp=xvgropen(fn,"Distance dependent Gk","r (nm)","G\\sk\\N(r)");
  ytot = 0;
  for(i=0; (i<ngraph); i++) {
    x = (i*rcut)/ngraph;
    if (count[i] > 0) {
      y = graph[i]/(mu*mu*count[i]);
    }
    else
      y = 0.0;
      
    ytot += y;
    
    fprintf(fp,"%10.5e  %12.7e  %12.7e  %12.7e\n",x,ytot,y,graph[i]);
  }
  ffclose(fp);
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


static void do_dip(char *fn,char *topf,char *outf,char *outfa, 
		   char *dipdist,bool bAverCorr,
		   bool bCorr,   char *corf,
		   bool bGkr,    char *gkrfn,
		   bool bQuad,   char *quadfn,
		   bool bMU,     char *mufn,
		   int gnx,atom_id grpindex[],
		   real mu_max,real mu,real epsilonRF,real temp,
		   bool bFA,int skip)
{
  FILE       *out,*outaver;
  static char *legoutf[] = { 
    "< M\\sx \\N>", 
    "< M\\sy \\N>",
    "< M\\sz \\N>",
    "< |M\\stot \\N| >"
  };
  static char *legoutfa[] = { 
    "< |M|\\S2\\N >", 
    "< |M| >\\S2\\N",
    "< |M|\\S2\\N > - < |M| >\\S2\\N",
    "G\\sk",
    "epsilon"
  };
  rvec       *x,*dipole=NULL,mu_t,M_av,M_av2,Q_av,Q_av2,*quadrupole=NULL;
  real       *gkr=NULL;
  int        *gkcount=NULL,nframes=1000;
  int        fmu=0,nre,timecheck=0;
  char       **enm=NULL;
  real       rcut=0;
  matrix     box;
  bool       bCont;
  real       t,t0,t1,dt;
  double     M_diff=0,epsilon;
  double     mu_ave,mu_mol,M2_ave=0,M_ave2=0;
  rvec       quad_ave,quad_mol;
  ivec       iMu;
  int        Vol=0;
  real       M_XX,M_YY,M_ZZ,M_XX2,M_YY2,M_ZZ2,Gk=0,g_k=0;
  real       **muall=NULL,**quadall=NULL;
  t_topology *top;
  t_atom     *atom=NULL;
  t_block    *mols=NULL;
  int        i,j,k,m,natom=0,nmol,status,teller,tel3;
  int        *bin,ibin;
  real       volume;
  unsigned long mode;
  /* PvM, I need these to be able to call read_tpx */
  int        natoms,step;
  real       lambda;

  snew(top,1);
  read_tpx(topf,&step,&t,&lambda,NULL,box,
	   &natoms,NULL,NULL,NULL,top);
  volume = det(box);
  if (!grpindex)
    gnx = top->blocks[ebMOLS].nr;
    
  if (bMU) {
    fmu = open_enx(mufn,"r");
    do_enxnms(fmu,&nre,&enm);

    Vol=-1;
    /* Determine the indexes of the energy grps we need */
    for (i=0; (i<nre); i++) {
      if (strstr(enm[i],"Volume"))
	Vol=i;
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
  
  if (Vol == -1)
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
  out=xvgropen(outf,"Total dipole moment of the simulation box vs. time",
	       "Time (ps)","Total Dipole Moment (Debye)");
  xvgr_legend(out,asize(legoutf),legoutf);
  xvgr_line_props(out, 0, elDashed, ecFrank);
  xvgr_line_props(out, 1, elDotDashed, ecFrank);
  xvgr_line_props(out, 2, elLongDashed, ecFrank);
  xvgr_line_props(out, 3, elSolid, ecFrank);
  outaver=xvgropen(outfa,"Averages of the total dipole moment vs. time",
		   "Time (ps)","<|M|\\S2\\N>, <|M|>\\S2\\N (Debye\\S2)");
  xvgr_legend(outaver,asize(legoutfa),legoutfa);

  teller = 0;
  /* Read the first frame from energy or traj file */
  if (bMU)
    do {
      bCont = read_mu_from_enx(fmu,Vol,iMu,mu_t,&volume,&t,teller,nre); 
      if (bCont) {  
	timecheck=check_times(t,t);
	if (timecheck < 0)
	  teller++;
	if ((teller % 10) == 0)
	  fprintf(stderr,"\r Skipping Frame %6d, time: %8.3f", teller, t);
      }
      else {
	fprintf(stderr,"End of %s reached\n",mufn);
	break;
      }
    } while (bCont && (timecheck < 0));
  else
    natom  = read_first_x(&status,fn,&t,&x,box);
  
  snew(bin, MAXBIN);
  epsilon    = 0.0;
  M_XX=M_XX2 = 0.0;
  M_YY=M_YY2 = 0.0;
  M_ZZ=M_ZZ2 = 0.0;
  mu_ave     = 0.0;

  if (bQuad)
    clear_rvec(quad_ave);

  /* Start while loop over frames */
  t1 = t0 = t;
  teller = 0;

  if (bGkr) {
    rcut       = 0.475*min(box[XX][XX],min(box[YY][YY],box[ZZ][ZZ]));
    snew(gkr,MAXBIN);
    snew(gkcount,MAXBIN);
  }
  teller=0;
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
	  M_av2[m] += dipole[i][m]*dipole[i][m];  /* M^2 per frame */
	  mu_mol   += dipole[i][m]*dipole[i][m];  /* calc. mu for distribution */
	}

	mu_ave += sqrt(mu_mol);                   /* calc. the average mu */
	if (bQuad) {
	  clear_rvec(quad_mol);
	  for(m=0; (m<DIM); m++) {
	    Q_av[m]  += quadrupole[i][m];                    /* Q per frame */
	    Q_av2[m] += quadrupole[i][m]*quadrupole[i][m];   /* Q^2 per frame */
	    quad_mol[m]  += quadrupole[i][m]*quadrupole[i][m];  
	    quad_ave[m] += quadrupole[i][m];
	  }
	}
	
	/* Update the dipole distribution */
	ibin = (MAXBIN*sqrt(mu_mol)/mu_max);
	if (ibin < MAXBIN)
	  bin[ibin]++;
      } /* End loop of all molecules in frame */
    } 
    
    if (bGkr) {
      init_pbc(box,FALSE);
      do_gkr(gnx,grpindex,mols->index,mols->a,x,dipole,
	     box,MAXBIN,gkr,gkcount,rcut,bFA);
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
      fprintf(out,"%10g  %12.8e %12.8e %12.8e %12.8e\n",
	      t, M_av[XX], M_av[YY], M_av[ZZ], norm(M_av));

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

    epsilon = calc_eps(M_diff,volume,epsilonRF,temp);

    /* Finite system Kirkwood G-factor */
    Gk = M_diff/(gnx*mu*mu);
    /* Infinite system Kirkwood G-factor */
    if (epsilonRF == 0.0) {
      g_k = ((2*epsilon+1)*Gk/(3*epsilon));
    } else {
       g_k = ((2*epsilonRF+epsilon)*(2*epsilon+1)*
	       Gk/(3*epsilon*(2*epsilonRF+1)));
    }

    /* Write to file < |M|^2 >, < |M| >^2. And the difference between 
     * the two. Here M is sum mu_i. Further write the finite system
     * Kirkwood G factor and epsilon.
     */
    if ((skip == 0) || ((teller % skip) == 0))
      fprintf(outaver,"%10g  %12.8e %12.8e %12.8e %12.8e %12.8e\n",
	      t, M2_ave, M_ave2, M_diff, Gk, epsilon);

    if (bMU)
      bCont = read_mu_from_enx(fmu,Vol,iMu,mu_t,&volume,&t,teller,nre); 
    else
      bCont = read_next_x(status,&t,natom,x,box);
  } while(bCont);
  
  if (!bMU)
    close_trj(status);
    
  fclose(out);
  fclose(outaver);

  if (bGkr) {
    print_gkr(gkrfn,MAXBIN,rcut,gkr,gkcount,mu);
    sfree(gkr);
    sfree(gkcount);
  }

  /* Autocorrelation function */  
  if (bCorr) {
    if (teller < 2) {
      fprintf(stderr,"Not enough frames for autocorrelation\n");
    }
    else {
      dt=(t1 - t0)/(teller-1);
      fprintf(stderr,"to %g, t %g, teller %d\n", t0,t,teller);
      
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
    fprintf(stderr,"\n\nAverage dipole moment (Debye)\n");
    fprintf(stderr," Tot= %g\n",  mu_ave/(gnx*teller));
    if (bQuad) {
      fprintf(stderr,"Average quadrupole moment (Debye-Ang)\n");
      fprintf(stderr," XX=  %g  YY=  %g ZZ=  %g norm= %g asymm= %g\n\n",  
	      quad_ave[XX]/(gnx*teller),
	      quad_ave[YY]/(gnx*teller),
	      quad_ave[ZZ]/(gnx*teller),
	      norm(quad_ave)/(gnx*teller),
	      (quad_ave[ZZ] - quad_ave[XX])/ quad_ave[YY]);
    }
  }
  fprintf(stderr,"The following averages for the complete trajectory have been calculated:\n\n");
  fprintf(stderr," Total < M_x > = %g Debye\n", M_XX/teller);
  fprintf(stderr," Total < M_y > = %g Debye\n", M_YY/teller);
  fprintf(stderr," Total < M_z > = %g Debye\n\n", M_ZZ/teller);

  fprintf(stderr," Total < M_x^2 > = %g Debye^2\n", M_XX2/teller);
  fprintf(stderr," Total < M_y^2 > = %g Debye^2\n", M_YY2/teller);
  fprintf(stderr," Total < M_z^2 > = %g Debye^2\n\n", M_ZZ2/teller);

  fprintf(stderr," Total < |M|^2 > = %g Debye^2\n", M2_ave);
  fprintf(stderr," Total < |M| >^2 = %g Debye^2\n\n", M_ave2);

  fprintf(stderr," < |M|^2 > - < |M| >^2 = %g Debye^2\n\n", M_diff);
  fprintf(stderr,"Finite system Kirkwood g factor G_k = %g\n", Gk);
  fprintf(stderr,"Infinite system Kirkwood g factor g_k = %g\n\n", g_k);
  fprintf(stderr,"Epsilon = %g\n", epsilon);

  if (!bMU) {
    /* Write to file the dipole moment distibution during the simulation.
     */
    out=xvgropen(dipdist,"Dipole Moment Distribution","mu (Debye)","");
    for(i=0; (i<MAXBIN); i++)
      fprintf(out,"%10g  %d\n",(i*mu_max)/MAXBIN,bin[i]);
    fclose(out);
    sfree(bin);
  }
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
  static real mu_max=5, mu=2.5;
  static real epsilonRF=0.0, temp=300;
  static bool bAverCorr=FALSE;
  static int  skip=0,nFA=0;
  t_pargs pa[] = {
    { "-mu",       FALSE, etREAL, {&mu},
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
      "Use the n-th atom of a molecule (water ?) to calculate the distance between molecules rather than the center of geometry (when -1) in the calculation of distance dependent Kirkwood factors" }
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
		    
  fprintf(stderr,"Using %g as mu_max and %g as the dipole moment.\n", 
	  mu_max, mu);
  if (epsilonRF == 0.0)
    fprintf(stderr,"WARNING: EpsilonRF = 0.0, this really means EpsilonRF = infinity\n");
  
  bMU = (opt2bSet("-enx",NFILE,fnm));
  if (bMU) {
    bAverCorr = TRUE;
    bQuad = FALSE;
    bGkr  = FALSE;
  }
  else {
    bQuad   = opt2bSet("-q",NFILE,fnm);
    bGkr    = opt2bSet("-g",NFILE,fnm);
  }
  if (ftp2bSet(efNDX,NFILE,fnm))
    rd_index(ftp2fn(efNDX,NFILE,fnm),1,&gnx,&grpindex,&grpname);
  else {
    gnx=1;
    grpindex=NULL;
  }
  bCorr   = (bAverCorr || opt2bSet("-c",NFILE,fnm));
  do_dip(ftp2fn(efTRX,NFILE,fnm),ftp2fn(efTPX,NFILE,fnm),
	 ftp2fn(efXVG,NFILE,fnm),
	 opt2fn("-a",NFILE,fnm),opt2fn("-d",NFILE,fnm),
	 bAverCorr,bCorr,
	 opt2fn("-c",NFILE,fnm),
	 bGkr,    opt2fn("-g",NFILE,fnm),
	 bQuad,   opt2fn("-q",NFILE,fnm),
	 bMU,     opt2fn("-enx",NFILE,fnm),
	 gnx,grpindex,mu_max,mu,epsilonRF,temp,nFA,skip);
  
  do_view(opt2fn("-o",NFILE,fnm),"-autoscale xy -nxy");
  do_view(opt2fn("-a",NFILE,fnm),"-autoscale xy -nxy");
  do_view(opt2fn("-d",NFILE,fnm),"-autoscale xy");
  do_view(opt2fn("-c",NFILE,fnm),"-autoscale xy");

  thanx(stderr);
  
  return 0;
}
