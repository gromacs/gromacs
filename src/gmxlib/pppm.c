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
 * Good ROcking Metal Altar for Chronical Sinners
 */
static char *SRCID_pppm_c = "$Id$";

#include <stdio.h>
#include <math.h>
#include "physics.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "xvgr.h"
#include "fatal.h"
#include "txtdump.h"
#include "network.h"
#include "nr.h"
#include "assert.h"
#include "nrnb.h"
#include "ghat.h"
#include "pppm.h"
#include "lrutil.h"
#include "mdrun.h"
#include "fftgrid.h"

#define llim  (-1)
#define ulim   (1)
#define llim2 (-3)
#define ulim2  (3)

static void calc_invh(rvec box,int nx,int ny,int nz,rvec invh)
{
  invh[XX] = nx/box[XX];
  invh[YY] = ny/box[YY];
  invh[ZZ] = nz/box[ZZ];
}

static void calc_weights(int nx,int ny,int nz,
			 rvec x,rvec box,rvec invh,ivec ixyz,real WXYZ[])
{
  const  real half=0.5;
  tensor wxyz;
  real   abc,ttt,bhh,fact;
#ifdef DEBUG
  real   wtot;
#endif
  ivec   nxyz;
  int    it,j,k,m,nm;
  real   Wx,Wy,Wzx,Wzy,Wzz;
  
  fact = 0.125;
  
  nxyz[XX] = nx;
  nxyz[YY] = ny;
  nxyz[ZZ] = nz;
  for(m=0; (m<DIM); m++) {
    /* Put particle in the box... */  
    ttt = x[m]*invh[m];
    
    /* Round to nearest grid point. Do the math in integer! */
    it  = (ttt+0.5);
    nm  = nxyz[m];
    if (it < 0) {
      ttt+=nm;
      it +=nm;
    }
    else if (it >= nm) {
      ttt-=nm;
      it -=nm;
    }
    if ((it < 0) && (it >= nxyz[m]))
      fatal_error(0,"it = %d, x=%f, ttt=%f",it,x[m],ttt);
    ixyz[m]    = it;
    
    /* Fraction (offset) from grid point */
    abc        = ttt - (real)ixyz[m];
    
    wxyz[m][0] = sqr((real)(half  - abc));
    wxyz[m][1] = 1.5 - 2.0*sqr(abc);
    wxyz[m][2] = sqr((real)(half  + abc));
  }
  Wzx=wxyz[ZZ][XX];
  Wzy=wxyz[ZZ][YY];
  Wzz=wxyz[ZZ][ZZ];
  for(j=m=0; (j<DIM); j++) {
    Wx = wxyz[XX][j]*fact;
    for(k=0; (k<DIM); k++,m+=3) {
      Wy        = Wx*wxyz[YY][k];
      WXYZ[m]   = Wy*Wzx;
      WXYZ[m+1] = Wy*Wzy;
      WXYZ[m+2] = Wy*Wzz;
    }
  }
#ifdef DEBUG
  wtot = 0;
  for(j=0; (j<27); j++)
    wtot+=WXYZ[j];
  fprintf(stderr,"wtot = %g\n",wtot);
#endif
#ifdef HACK
  for(j=0; (j<27); j++)
    WXYZ[j] = 0;
  WXYZ[13] = 1.0;
#endif
}
	
static void calc_nxyz(int nx,int ny,int nz,
		      int **nnx,int **nny,int **nnz)
{
  int i;
  
  snew(*nnx,3*nx);
  snew(*nny,3*ny);
  snew(*nnz,3*nz);
  for(i=0; (i<3*nx); i++)
    (*nnx)[i] = i % nx;
  for(i=0; (i<3*ny); i++)
    (*nny)[i] = i % ny;
  for(i=0; (i<3*nz); i++)
    (*nnz)[i] = i % nz;
}
	
static void spread_q(FILE *log,bool bVerbose,
		     int natoms,rvec x[],real charge[],rvec box,
		     t_fftgrid *grid,t_nrnb *nrnb)
{
  static bool bFirst = TRUE;
  static int  *nnx,*nny,*nnz;
  rvec   invh;
  real   qi,qwt;
#ifdef DEBUG
  real   qt;
#endif
  real   WXYZ[27];
  ivec   ixyz;
  int    i,iX,iY,iZ,index;
  int    jx,jy,jz,jcx,jcy,jcz;
  int    nxyz;
  int    nx,ny,nz,la1,la2,la12;
  t_fft_tp *ptr;
  
  unpack_fftgrid(grid,&nx,&ny,&nz,&la1,&la2,&la12,&ptr);
  
  calc_invh(box,nx,ny,nz,invh);

  if (bFirst) {
    fprintf(log,"Spreading Charges using Triangle Shaped on %dx%dx%d grid\n",
	    nx,ny,nz);
    fprintf(log,"invh = %10g,%10g,%10g\n",invh[XX],invh[YY],invh[ZZ]);
  
    calc_nxyz(nx,ny,nz,&nnx,&nny,&nnz);
    
    bFirst = FALSE;
  }

  for(i=0; (i<natoms); i++) {
    qi=charge[i];

    /* Each charge is spread over the nearest 27 grid cells,
     * thus we loop over -1..1 in X,Y and Z direction
     * We apply the TSC (triangle shaped charge)
     * see Luty et. al, JCP 103 (1995) 3014
     */
    
    if (qi != 0.0) {
      calc_weights(nx,ny,nz,x[i],box,invh,ixyz,WXYZ);
      iX  = ixyz[XX] + nx;
      iY  = ixyz[YY] + ny;
      iZ  = ixyz[ZZ] + nz;

#ifdef DEBUG
      qt=0;
#endif
      nxyz = 0;
      for(jx=-1; (jx<=1); jx++) {
	jcx = nnx[iX + jx];
	for(jy=-1; (jy<=1); jy++) {
	  jcy = nny[iY + jy];
	  for(jz=-1; (jz<=1); jz++,nxyz++) {
	    jcz   = nnz[iZ + jz];
	    index = INDEX(jcx,jcy,jcz);
	    qwt   = qi*WXYZ[nxyz];
	    GR_INC(ptr[index],qwt);
#ifdef DEBUG
	    qt   += qwt;
	    if (bVerbose)
	      fprintf(log,"spread %4d %2d %2d %2d  %10.3e, weight=%10.3e\n",
		      index,jcx,jcy,jcz,GR_VALUE(ptr[index]),WXYZ[nxyz]);
#endif
	  }
	}
      }
#ifdef DEBUG
      fprintf(log,"q[%3d] = %6.3f, qwt = %6.3f\n",i,qi,qt);
#endif
    }
  }
  inc_nrnb(nrnb,eNR_SPREADQ,9*natoms);
  inc_nrnb(nrnb,eNR_WEIGHTS,3*natoms);
}

real gather_inner(int JCXYZ[],real WXYZ[],int ixw[],int iyw[],int izw[],
		  int la1,int la2,int la12,
		  real c1x,real c1y,real c1z,real c2x,real c2y,real c2z,
		  real qi,rvec f,t_fft_tp ptr[])
{
  real pi,fX,fY,fZ,weight;
  int  jxyz,m,jcx,jcy,jcz;
  int  jcx0,jcy0,jcz0;
  
  pi = 0.0;
  fX = 0.0;
  fY = 0.0;
  fZ = 0.0;
  
  /* Now loop over 27 surrounding vectors */      
  for(jxyz=m=0; (jxyz < 27); jxyz++,m+=3) {
    jcx    = JCXYZ[m];
    jcy    = JCXYZ[m+1];
    jcz    = JCXYZ[m+2];
    weight = WXYZ[jxyz];
    
    jcx0   = ixw[jcx];
    jcy0   = iyw[jcy];
    jcz0   = izw[jcz];

    /* Electrostatic Potential ! */
    pi += weight * GR_VALUE(ptr[INDEX(jcx0,jcy0,jcz0)]);

    /* Forces */
    fX += weight * ((c1x*(GR_VALUE(ptr[INDEX(ixw[jcx-1],jcy0,jcz0)]) - 
			  GR_VALUE(ptr[INDEX(ixw[jcx+1],jcy0,jcz0)]) )) +
		    (c2x*(GR_VALUE(ptr[INDEX(ixw[jcx-2],jcy0,jcz0)]) - 
			  GR_VALUE(ptr[INDEX(ixw[jcx+2],jcy0,jcz0)]) )));
    fY += weight * ((c1y*(GR_VALUE(ptr[INDEX(jcx0,iyw[jcy-1],jcz0)]) -
			  GR_VALUE(ptr[INDEX(jcx0,iyw[jcy+1],jcz0)]) ))  +
		    (c2y*(GR_VALUE(ptr[INDEX(jcx0,iyw[jcy-2],jcz0)]) -
			  GR_VALUE(ptr[INDEX(jcx0,iyw[jcy+2],jcz0)]) )));
    fZ += weight * ((c1z*(GR_VALUE(ptr[INDEX(jcx0,jcy0,izw[jcz-1])]) -
			  GR_VALUE(ptr[INDEX(jcx0,jcy0,izw[jcz+1])]) ))  +
		    (c2z*(GR_VALUE(ptr[INDEX(jcx0,jcy0,izw[jcz-2])]) -
			  GR_VALUE(ptr[INDEX(jcx0,jcy0,izw[jcz+2])]) )));
  }
  f[XX] += qi*fX;
  f[YY] += qi*fY;
  f[ZZ] += qi*fZ;
  
  return pi;
}

static real gather_f(FILE *log,bool bVerbose,
		     int natoms,rvec x[],rvec f[],real charge[],rvec box,
		     real pot[],t_fftgrid *grid,rvec beta,t_nrnb *nrnb)
{
  static bool bFirst=TRUE;
  static int  *nnx,*nny,*nnz;
  static int  JCXYZ[81];
  int    i,m;
  real   energy;
  real   qi,pi;
  ivec   ixyz;
  rvec   invh,c1,c2;
  real   WXYZ[27];
  real   c1x,c1y,c1z,c2x,c2y,c2z;
  int    ixw[7],iyw[7],izw[7];
  int    ll;
  int    nx,ny,nz,la1,la2,la12;
  t_fft_tp *ptr;
  
  unpack_fftgrid(grid,&nx,&ny,&nz,&la1,&la2,&la12,&ptr);
  
  calc_invh(box,nx,ny,nz,invh);
  
  for(m=0; (m<DIM); m++) {
    c1[m] = (beta[m]/2.0)*invh[m];
    c2[m] = ((1.0-beta[m])/4.0)*invh[m];
  }
  c1x = c1[XX];
  c1y = c1[YY];
  c1z = c1[ZZ];
  c2x = c2[XX];
  c2y = c2[YY];
  c2z = c2[ZZ];

  if (bFirst) {
    fprintf(log,"Gathering Forces using Triangle Shaped on %dx%dx%d grid\n",
	    nx,ny,nz);
    fprintf(log,"beta = %10g,%10g,%10g\n",beta[XX],beta[YY],beta[ZZ]);
    fprintf(log,"c1   = %10g,%10g,%10g\n",c1[XX],c1[YY],c1[ZZ]);
    fprintf(log,"c2   = %10g,%10g,%10g\n",c2[XX],c2[YY],c2[ZZ]);
    fprintf(log,"invh = %10g,%10g,%10g\n",invh[XX],invh[YY],invh[ZZ]);

    calc_nxyz(nx,ny,nz,&nnx,&nny,&nnz);

    for(i=0; (i<27); i++) {
      JCXYZ[3*i]   = 2 + (i/9);
      JCXYZ[3*i+1] = 2 + (i/3) % 3;
      JCXYZ[3*i+2] = 2 + (i % 3); 
    }
    
    bFirst = FALSE;
  }

  energy=0.0;  	  
  for(i=0; (i<natoms); i++) {
    /* Each charge is spread over the nearest 27 grid cells,
     * thus we loop over -1..1 in X,Y and Z direction
     * We apply the TSC (triangle shaped charge)
     * see Luty et. al, JCP 103 (1995) 3014
     */
     
    calc_weights(nx,ny,nz,x[i],box,invh,ixyz,WXYZ);

    for(ll=llim2; (ll<=ulim2); ll++) {
      ixw[ll-llim2] = nnx[ixyz[XX]+ll+nx];
      iyw[ll-llim2] = nny[ixyz[YY]+ll+ny];
      izw[ll-llim2] = nnz[ixyz[ZZ]+ll+nz];
    }
    
    qi      = charge[i];
    pi      = gather_inner(JCXYZ,WXYZ,ixw,iyw,izw,la1,la2,la12,
			   c1x,c1y,c1z,c2x,c2y,c2z,
			   qi,f[i],ptr);
    
    energy += pi*qi;
    pot[i]  = pi;
  }
  
  inc_nrnb(nrnb,eNR_GATHERF,27*natoms);
  inc_nrnb(nrnb,eNR_WEIGHTS,3*natoms);
  
  return energy*0.5;
}

void convolution(FILE *fp,bool bVerbose,t_fftgrid *grid,real ***ghat)
{
  int      i,j,k,index;
  real     gk;
  int      nx,ny,nz,la1,la2,la12;
  t_fft_tp *ptr;
  int      *nTest;
  
  unpack_fftgrid(grid,&nx,&ny,&nz,&la1,&la2,&la12,&ptr);
  
  snew(nTest,grid->nptr);
  /* CHECK THIS.... */
#ifdef USE_SGI_FFT
  for(i=0; (i<nx/2); i++) {
#else
  for(i=0; (i<nx); i++) {
#endif
    for(j=0; (j<ny); j++) {
      for(k=0; (k<nz); k++) {
#ifdef USE_SGI_FFT
	gk    = ghat[k][j][i];
	index = INDEX(2*i,j,k);
	ptr[index] *= gk;
	nTest[index]++;
	/* fprintf(fp,"SGI: index = %d\n",index); */
	index = INDEX(2*i+1,j,k);
	ptr[index] *= gk;
	nTest[index]++;
	/* fprintf(fp,"SGI: index = %d\n",index); */
#else
	gk    = ghat[i][j][k];
	index = INDEX(i,j,k);
	ptr[index].re *= gk;
	ptr[index].im *= gk;
	nTest[index]++;
#endif
      }
    }
  }
  for(i=0; (i<nx); i++) {
    for(j=0; (j<ny); j++) {
      for(k=0; (k<nz); k++) {
	index = INDEX(i,j,k);
	if (nTest[index] != 1)
	  fprintf(fp,"Index %d sucks, set %d times\n",index,nTest[index]);
      }
    }
  }
  sfree(nTest);
}

void solve_pppm(FILE *fp,t_commrec *cr,
		t_fftgrid *grid,real ***ghat,rvec box,
		bool bVerbose,t_nrnb *nrnb)
{
  int  ntot,npppm;
  
  if (bVerbose) 
    print_fftgrid(fp,"Q-Real",grid,grid->nxyz,"qreal.pdb",box,TRUE);
  
  gmxfft3D(fp,bVerbose,grid,FFTW_FORWARD);
  
  if (bVerbose) {
    print_fftgrid(fp,"Q-k",grid,1.0,"qk-re.pdb",box,TRUE);
    print_fftgrid(fp,"Q-k",grid,1.0,"qk-im.pdb",box,FALSE);
    fprintf(stderr,"Doing convolution\n");
  }
  
  convolution(fp,bVerbose,grid,ghat); 
  
  if (bVerbose) 
    print_fftgrid(fp,"Convolution",grid,1.0,
		  "convolute.pdb",box,TRUE);
  
  gmxfft3D(fp,bVerbose,grid,FFTW_BACKWARD);
  
  if (bVerbose) 
    print_fftgrid(fp,"Potential",grid,1.0,"pot.pdb",box,TRUE);
  
  ntot  = grid->nxyz;  
  npppm = ntot*log((real)ntot)/log(2.0);
  inc_nrnb(nrnb,eNR_FFT,2*npppm);
  inc_nrnb(nrnb,eNR_CONV,ntot);
}

real do_pppm(FILE *log,       bool bVerbose,
	     bool bGenerGhat, char *ghatfn,
	     t_inputrec *ir,  int natoms,
	     rvec x[],        rvec f[],
	     real charge[],   rvec box,
	     real phi[],      t_commrec *cr,
	     t_nrnb *nrnb)
{
  static  bool bFirst = TRUE;
  static  real      ***ghat;
  static  t_fftgrid *grid;
  static  rvec      beta;
  static  int       porder;
  
  const     real tol = 1e-5;
  int       m;
  real      ener,r1,rc;
  ivec      grids;
  rvec      spacing;
  
  ener = 0.0;
  
  if (bFirst) {
    int nx,ny,nz;
    
    if (cr != NULL) {
      if (cr->nprocs > 1)
	fatal_error(0,"No parallel PPPM yet...");
    }
    if (bGenerGhat) {    
      fprintf(log,"Generating Ghat function\n");
      beta[XX]=beta[YY]=beta[ZZ]= 4.0/3.0;
      porder = 2;
      nx     = ir->nkx;
      ny     = ir->nky;
      nz     = ir->nkz;
      ghat   = mk_rgrid(nx,ny,nz);
      mk_ghat(NULL,nx,ny,nz,ghat,box,ir->rshort,ir->rlong,TRUE);
    }
    else {
      fprintf(stderr,"Reading Ghat function from %s\n",ghatfn);
      ghat = rd_ghat(log,ghatfn,grids,spacing,beta,&porder,&r1,&rc);
      
      /* Check whether cut-offs correspond */
      if ((fabs(r1-ir->rshort) > tol) || (fabs(rc-ir->rlong) > tol)) {
	fprintf(log,"rshort = %10.3e  rlong = %10.3e"
		"  r1 = %10.3e  rc = %10.3e\n",ir->rshort,ir->rlong,r1,rc);
	fflush(log);
	fatal_error(0,"Cut-off lengths in tpb file and Ghat file %s "
		    "do not match\nCheck your log file!",ghatfn);
      }
      
      /* Check whether boxes correspond */
      for(m=0; (m<DIM); m++)
	if (fabs(box[m]-grids[m]*spacing[m]) > tol) {
	  pr_rvec(log,0,"box",box,DIM);
	  pr_rvec(log,0,"grid-spacing",spacing,DIM);
	  pr_ivec(log,0,"grid size",grids,DIM);
	  fflush(log);
	  fatal_error(0,"Box sizes in tpb file and Ghat file %s do not match\n"
		      "Check your log file!",ghatfn);
	}
      nx = grids[XX];
      ny = grids[YY];
      nz = grids[ZZ];
    }

    fprintf(log,"Will use the PPPM algorithm for long-range electrostatics\n");
    fprintf(log,"Grid size is %d x %d x %d\n",nx,ny,nz);

    if ((nx < 4) || (ny < 4) || (nz < 4)) 
      fatal_error(0,"Grid must be at least 4 points in all directions");
    
    if (porder != 2)
      fatal_error(0,"porder = %d, should be 2 in %s",porder,ghatfn);
    
    grid = mk_fftgrid(nx,ny,nz);

    bFirst = FALSE;
  }
  else {
    /* We don't do PPPM the first time around, that is only
     * for setting things up.
     */

    /* Now start the actual PPPM procedure.
     * First step: spreading the charges over the grid.
     */
    /* Make the grid empty */
    clear_fftgrid(grid);

    spread_q(log,bVerbose,natoms,x,charge,box,grid,nrnb);
    
    /* Second step: solving the poisson equation in Fourier space */
    solve_pppm(log,NULL,grid,ghat,box,bVerbose,nrnb);
    
    /* Third and last step: gather the forces, energies and potential
     * from the grid.
     */
    ener=gather_f(log,bVerbose,natoms,x,f,charge,box,phi,grid,beta,nrnb);
  }
  
  return ener;
}

