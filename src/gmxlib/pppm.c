/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * Great Red Oystrich Makes All Chemists Sane
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
#include "grids.h"
#include "nrnb.h"
#include "ghat.h"
#include "pppm.h"
#include "lrutil.h"
#include "mdrun.h"

#ifdef USEF77
DECLAREF77(doconv) (int *n,real grid[],real gk[]);
#endif

#define INDEX(i,j,k,la1,la2)  ((k*la1*la2)+(j*la1)+i)

#define llim  (-1)
#define ulim   (1)
#define llim2 (-3)
#define ulim2  (3)

void rlft3f_()
{
  fprintf(stderr,"Just called rlft3f_!\n");
  exit(1);
}

static void calc_invh(rvec box,int nx,int ny,int nz,rvec invh)
{
  invh[XX] = nx/box[XX];
  invh[YY] = ny/box[YY];
  invh[ZZ] = nz/box[ZZ];
}

static void calc_weights(rvec x,rvec box,rvec invh,ivec ixyz,real WXYZ[])
{
  const  real half=0.5;
  tensor wxyz;
  real   abc,ttt,bhh,fact;
  int    i,j,k,l,m;
  real   Wx,Wy,Wzx,Wzy,Wzz;
  
  fact = 0.125;
  
  for(m=0; (m<DIM); m++) {
    /* Put particle in the box... */  
    ttt = x[m]*invh[m];
    bhh = box[m]*invh[m];
    if (ttt < 0)
      ttt += bhh;
    else if (ttt >= bhh)
      ttt -= bhh;
    
    /* Calculate nearest grid point, Round */
    ixyz[m]    = ttt+half;
    
    /* Fraction (offset) from grid point */
    abc        = ttt - (real)ixyz[m];
    
    wxyz[m][0] = sqr(half  - abc);
    wxyz[m][1] = 1.5 - 2.0*sqr(abc);
    wxyz[m][2] = sqr(half  + abc);
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
		     int nx,int ny,int nz,int la1,int la2,
		     real ***grid,t_nrnb *nrnb)
{
  static bool bFirst = TRUE;
  static int  *nnx,*nny,*nnz;
  int    i,j,k,m,indX,indY;
  rvec   invh;
  real   qi,rix,riy,riz;
  real   WXYZ[27];
  real   *gptr;
  ivec   ixyz;
  int    iX,iY,iZ,index;
  int    jx,jy,jz,jcx,jcy,jcz;
  int    nxyz;
  
  calc_invh(box,nx,ny,nz,invh);

  if (bFirst) {
    fprintf(log,"Spreading Charges using Triangle Shaped on %dx%dx%d grid\n",
	    nx,ny,nz);
    fprintf(log,"invh = %10g,%10g,%10g\n",invh[XX],invh[YY],invh[ZZ]);
  
    calc_nxyz(nx,ny,nz,&nnx,&nny,&nnz);
    
    bFirst = FALSE;
  }

  gptr=grid[0][0];
  for(i=0; (i<natoms); i++) {
    qi=charge[i];

    /* Each charge is spread over the nearest 27 grid cells,
     * thus we loop over -1..1 in X,Y and Z direction
     * We apply the TSC (triangle shaped charge)
     * see Luty et. al, JCP 103 (1995) 3014
     */
    
    if (qi != 0.0) {
      calc_weights(x[i],box,invh,ixyz,WXYZ);
      iX  = ixyz[XX] + nx;
      iY  = ixyz[YY] + ny;
      iZ  = ixyz[ZZ] + nz;
      
      nxyz = 0;
      for(jx=-1; (jx<=1); jx++) {
	jcx = nnx[iX + jx];
	for(jy=-1; (jy<=1); jy++) {
	  jcy = nny[iY + jy];
	  for(jz=-1; (jz<=1); jz++,nxyz++) {
	    jcz = nnz[iZ + jz];
#ifdef USE_SGI_FFT
	    index=INDEX(jcx,jcy,jcz,la1,la2);
	    gptr[index] += qi*WXYZ[nxyz];
#ifdef DEBUG
	    if (bVerbose)
	      fprintf(log,"spread %4d %2d %2d %2d  %10.3e\n",
		      index,jcx,jcy,jcz,grid[0][0][index]);
#endif
#else
	    grid[jcx][jcy][jcz] += qi*WXYZ[nxyz];
#endif
	  }
	}
      }
    }
  }
  inc_nrnb(nrnb,eNR_SPREADQ,9*natoms);
  inc_nrnb(nrnb,eNR_WEIGHTS,3*natoms);
}

real gather_inner(int JCXYZ[],real WXYZ[],int ixw[],int iyw[],int izw[],
		  int la1,int la2,
		  real c1x,real c1y,real c1z,real c2x,real c2y,real c2z,
		  real qi,rvec f,real ***grid)
{
  real *gptr;
  real pi,fX,fY,fZ,weight;
  int  jxyz,m,jcx,jcy,jcz;
  int  jcx0,jcy0,jcz0;
  
  pi  = 0.0;
  fX  = 0.0;
  fY  = 0.0;
  fZ  = 0.0;
  gptr = grid[0][0];
  
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
#ifdef USE_SGI_FFT
    pi += weight * gptr[INDEX(jcx0,jcy0,jcz0,la1,la2)];
#else
    pi += weight * grid[jcx0][jcy0][jcz0];
#endif

    /* Forces */
#ifdef USE_SGI_FFT
    fX += weight * ((c1x*( gptr[INDEX(ixw[jcx-1],jcy0,jcz0,la1,la2)] -
			   gptr[INDEX(ixw[jcx+1],jcy0,jcz0,la1,la2)] ))  +
		    (c2x*( gptr[INDEX(ixw[jcx-2],jcy0,jcz0,la1,la2)] -
			   gptr[INDEX(ixw[jcx+2],jcy0,jcz0,la1,la2)] )));
    fY += weight * ((c1y*( gptr[INDEX(jcx0,iyw[jcy-1],jcz0,la1,la2)] -
			   gptr[INDEX(jcx0,iyw[jcy+1],jcz0,la1,la2)] ))  +
		    (c2y*( gptr[INDEX(jcx0,iyw[jcy-2],jcz0,la1,la2)] -
			   gptr[INDEX(jcx0,iyw[jcy+2],jcz0,la1,la2)] )));
    fZ += weight * ((c1z*( gptr[INDEX(jcx0,jcy0,izw[jcz-1],la1,la2)] -
			   gptr[INDEX(jcx0,jcy0,izw[jcz+1],la1,la2)] ))  +
		    (c2z*( gptr[INDEX(jcx0,jcy0,izw[jcz-2],la1,la2)] -
			   gptr[INDEX(jcx0,jcy0,izw[jcz+2],la1,la2)] )));
#else
    fX += weight * ((c1x*( grid[ixw[jcx-1]][jcy0][jcz0] -
			   grid[ixw[jcx+1]][jcy0][jcz0] ))  +
		    (c2x*( grid[ixw[jcx-2]][jcy0][jcz0] -
			   grid[ixw[jcx+2]][jcy0][jcz0] )));
    fY += weight * ((c1y*( grid[jcx0][iyw[jcy-1]][jcz0] -
			   grid[jcx0][iyw[jcy+1]][jcz0] ))  +
		    (c2y*( grid[jcx0][iyw[jcy-2]][jcz0] -
			   grid[jcx0][iyw[jcy+2]][jcz0] )));
    fZ += weight * ((c1z*( grid[jcx0][jcy0][izw[jcz-1]] -
			   grid[jcx0][jcy0][izw[jcz+1]] ))  +
		    (c2z*( grid[jcx0][jcy0][izw[jcz-2]] -
			   grid[jcx0][jcy0][izw[jcz+2]] )));
#endif
  }
  f[XX] += qi*fX;
  f[YY] += qi*fY;
  f[ZZ] += qi*fZ;
  
  return pi;
}

static real gather_f(FILE *log,bool bVerbose,
		     int natoms,rvec x[],rvec f[],real charge[],rvec box,
		     int nx,int ny,int nz,int la1,int la2,real pot[],
		     real ***grid,real beta,
		     t_nrnb *nrnb)
{
  static bool bFirst=TRUE;
  static int  *nnx,*nny,*nnz;
  static int  JCXYZ[81];
  int    i,j,k,l,m;
  real   energy;
  real   w1,weight,fact,rhs,rix,riy,riz,qi,pi;
  ivec   ixyz;
  rvec   invh,c1,c2;
  real   WXYZ[27];
  real   c1x,c1y,c1z,c2x,c2y,c2z,fX,fY,fZ;
  int    ixw[7],iyw[7],izw[7];
  int    ll,jx,jy,jz,jcx,jcy,jcz;
  int    jcx_2,jcx_1,jcx0,jcx1,jcx2;
  int    jcy_2,jcy_1,jcy0,jcy1,jcy2;
  int    jcz_2,jcz_1,jcz0,jcz1,jcz2;
  int    jxyz;
    
  calc_invh(box,nx,ny,nz,invh);
  
  fact = 0.125;
  
  for(m=0; (m<DIM); m++) {
    c1[m] = (beta/2.0)*invh[m];
    c2[m] = ((1.0-beta)/4.0)*invh[m];
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
    fprintf(log,"beta = %10g\n",beta);
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
     
    calc_weights(x[i],box,invh,ixyz,WXYZ);

    for(ll=llim2; (ll<=ulim2); ll++) {
      ixw[ll-llim2] = nnx[ixyz[XX]+ll+nx];
      iyw[ll-llim2] = nny[ixyz[YY]+ll+ny];
      izw[ll-llim2] = nnz[ixyz[ZZ]+ll+nz];
    }
    
    qi      = charge[i];
    pi      = gather_inner(JCXYZ,WXYZ,ixw,iyw,izw,la1,la2,
			   c1x,c1y,c1z,c2x,c2y,c2z,
			   qi,f[i],grid);
    
    energy += pi*qi;
    pot[i]  = pi;
  }
  
  inc_nrnb(nrnb,eNR_GATHERF,27*natoms);
  inc_nrnb(nrnb,eNR_WEIGHTS,3*natoms);
  
  return energy*0.5;
}

void convolution(FILE *fp,bool bVerbose,
		 int nx,int ny,int nz,int la1,int la2,
		 real ***grid,real ***ghat)
{
  int  i,i0,i1,j,k,index;
  real gk,fac;
  
  fac=1.0;
  for(k=0; (k<nz); k++) {
    for(j=0; (j<ny); j++) {
      for(i=0; (i<=nx/2); i++) {
	gk    = fac*ghat[i][j][k];
	index = INDEX(2*i,j,k,la1,la2);
	grid[0][0][index] *= gk;
#ifdef DEBUG
	if (bVerbose)
	  fprintf(fp,"convol %4d %2d %2d %2d  %10.3e  %10.3e\n",
		  index,i,j,k,gk,grid[0][0][index]);
#endif
	index = index+1;
	grid[0][0][index] *= gk;
#ifdef DEBUG
	if (bVerbose)
	  fprintf(fp,"convol %4d %2d %2d %2d  %10.3e  %10.3e\n",
		  index,i,j,k,gk,grid[0][0][index]);
#endif
      }
    }
  }
}

void printit(FILE *fp,
	     char *s,int nx,int ny,int nz,int la1,int la2,real arr[],
	     real factor)
{
  int i,j,k;

  for(k=0; (k<nz); k++)
    for(j=0; (j<ny); j++) {
      fprintf(fp,"%7s[%2d][%2d](index=%4d) = ",s,k,j,INDEX(0,j,k,la1,la2));
      for(i=0; (i<nx); i++)
	fprintf(fp,"  %8.1e",factor*arr[INDEX(i,j,k,la1,la2)]);
      fprintf(fp,"\n");
    }
}

void solve_pppm(FILE *fp,t_commrec *cr,
		real ***grid,real ***ghat,rvec box,
		int nx,int ny,int nz,int la1,int la2,
		bool bVerbose,t_nrnb *nrnb)
{
  static    bool bFirst=TRUE;
#ifdef USE_SGI_FFT
#include <fft.h>
  static    real *coeff;
#else
  static    real *speq;
#endif
  real      *cptr,*gptr,*fqqq,fg,fac;
  int       ntot,i,j,k,m,n,sign;
  int       npppm;
  
  ntot    = nx*ny*nz;
  cptr    = grid[0][0];
  gptr    = ghat[0][0];
  
  if (bVerbose) {
    printit(fp,"b4",nx,ny,1,la1,la2,grid[0][0],1.0);
    print_rgrid_pdb("q-spread.pdb",nx,ny,nz,grid);
    fprintf(stderr,"Doing forward 3D-FFT\n");
  }
  if (bFirst) {
#ifdef USE_SGI_FFT
    fprintf(fp,"Going to use SGI optimized FFT routines.\n");
#ifdef DOUBLE
    coeff  = dzfft3dui(nx,ny,nz,NULL);
#else
    coeff  = scfft3dui(nx,ny,nz,NULL);
#endif
#else
    snew(speq,ny*nz);
#endif
    bFirst = FALSE;
  }
#ifdef USE_SGI_FFT
#ifdef DOUBLE
  dzfft3du(1,nx,ny,nz,cptr,la1,la2,coeff);
#else
  scfft3du(1,nx,ny,nz,cptr,la1,la2,coeff);
#endif
#else
  sign = 1;
  rlft3f_(cptr,speq,&nx,&ny,&nz,&sign);
#endif
  if (bVerbose) {
    (void) print_rgrid(fp,"Qk",nx,ny,nz,grid);
    printit(fp,"fft",la1,ny,1,la1,la2,grid[0][0],1.0);
    fprintf(stderr,"Doing convolution\n");
  }
  
  convolution(fp,bVerbose,nx,ny,nz,la1,la2,grid,ghat);

  if (bVerbose) 
    fprintf(stderr,"Doing backward 3D-FFT\n");
  
#ifdef USE_SGI_FFT
#ifdef DOUBLE
  zdfft3du(-1,nx,ny,nz,cptr,la1,la2,coeff);
#else
  csfft3du(-1,nx,ny,nz,cptr,la1,la2,coeff);
#endif
#else
  sign = -1;
  rlft3f_(cptr,speq,&nx,&ny,&nz,&sign);
#endif
  
  if (bVerbose) {
    printit(fp,"after",nx,ny,1,la1,la2,grid[0][0],1.0/(nx*ny*nz));
    print_rgrid_pdb("phi.pdb",nx,ny,nz,grid);
  }
  
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
  static  real      ***grid;
  static  real      beta;
  static  int       porder,nx,ny,nz,la1,la2;
  
  static  int       niter;
  
  const     real tol = 1e-6;
  int       i,m;
  real      ctot;
  real      aver,tot,ener,r1,rc;
  ivec      grids;
  rvec      spacing;
  
  ener = 0.0;
  
  fatal_error(0,"PPPM not debugged yet.");

  if (bFirst) {
    if (cr != NULL) {
      if (cr->nprocs > 1)
	fatal_error(0,"No parallel PPPM yet...");
    }
    niter = ir->niter;
    if (bGenerGhat) {    
      fprintf(log,"Generating Ghat function\n");
      beta   = 4.0/3.0;
      porder = 2;
      nx     = ir->nkx;
      ny     = ir->nky;
      nz     = ir->nkz;
      ghat   = mk_rgrid(nx,ny,nz);
      mk_ghat(NULL,nx,ny,nz,ghat,box,ir->rshort,ir->rlong,TRUE);
    }
    else {
      fprintf(stderr,"Reading Ghat function from %s\n",ghatfn);
      ghat = rd_ghat(log,ghatfn,grids,spacing,&beta,&porder,&r1,&rc);
      
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
    
    if (bVerbose) 
      print_rgrid(log,"ghat",nx,ny,nz,ghat);
#ifdef USE_SGI_FFT
    if (nx != nz) {
      fprintf(stderr,"You can't use SGI optimized FFT routines when the number of grid points is not the same in X and Z directions. Sorry\n");
      exit(1);
    }
    grid = mk_rgrid(nz,ny,nx+2);
    la1  = nx+2;
    la2  = ny;
#else
    grid = mk_rgrid(nx,ny,nz);
#endif  
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
    clear_rgrid(nx,ny,nz,grid);

    spread_q(log,bVerbose,natoms,x,charge,box,nx,ny,nz,la1,la2,grid,nrnb);
    
    if (bVerbose) {
      ctot = print_rgrid(log,"qqq-spread(P3M)",nx,ny,nz,grid);
      aver = ctot/(nx*ny*nz);
      fprintf(log,"qqqtot = %10.3e, qqqav = %10.3e\n",ctot,aver);
    }
    
    /* Second step: solving the poisson equation in Fourier space */
    solve_pppm(log,NULL,grid,ghat,box,nx,ny,nz,la1,la2,bVerbose,nrnb);
    
    if (bVerbose) {
      ctot = print_rgrid(log,"phi",nx,ny,nz,grid);
      aver = ctot/(nx*ny*nz);
      fprintf(log,"phitot = %10.3e, phiav = %10.3e\n",ctot,aver);
    }
    
    /* Third and last step: gather the forces, energies and potential
     * from the grid.
     */
    ener=gather_f(log,bVerbose,natoms,x,f,charge,box,nx,ny,nz,la1,la2,phi,
		  grid,beta,nrnb);
  }
  
  return ener;
}

void transpose(FILE *log,t_commrec *cr,real *fqqq,int nx,int ny,int nz,
	       bool bForward)
{
  real ***tmp;
  int  i,j,k,n,nf,nt,ntot;
  
  tmp=mk_rgrid(nx,ny,nz*2);
    
  n=1;
  for(i=1; (i<=nx); i++)
    for(j=1; (j<=ny); j++)
      for(k=1; (k<=2*nz); k++,n++)
	tmp[i][j][k] = fqqq[n];
	
  n=1;
  for(i=1; (i<=nx); i++)
    for(j=1; (j<=ny); j++) {
      if (bForward) 
	for(k=1; (k<=nz); k++) {
	  fqqq[n++] = tmp[nz-k+1][j][2*i-1];
	  fqqq[n++] = tmp[nz-k+1][j][2*i];
	}
      else
	for(k=1; (k<=nz); k++) {
	  fqqq[n++] = tmp[k][j][2*(nx-i)+1];
	  fqqq[n++] = tmp[k][j][2*(nx-i)+2];
	}
    }
  free_rgrid(tmp,nx,ny);
}

void do_par_solve(FILE *log,t_commrec *cr,
		  real *fqqq,real ***ghat,rvec box,
		  int nx,int ny,int nz,bool bVerbose)
{
  real fac;
  int  ntot,i,j,k,m,n,ndim[4];
  
  ndim[0]=0;
  ndim[1]=ny;
  ndim[2]=nz;
  ntot=ny*nz*2;
  
  fprintf(stderr,"Doing %d forward 2D-FFTs\n",nx);
  
  for(i=0; (i<nx); i++) 
    fourn(fqqq+ntot*i,ndim,2,1);
  
  if (1) {
  }
  else {
    fprintf(stderr,"Doing forward transpose\n");
    transpose(log,cr,fqqq,nx,ny,nz,TRUE);
    
    /*
      fprintf(stderr,"Doing %d forward 1D-FFTs\n",ntot);
    for(i=1; (i<nx); i++) {
      for(j=1; (j<ny); j++) {
	four1(fqqq+((i-1)*nx+j-1)*ny,nz,1);
      }
    }
    
    fprintf(stderr,"Doing convolution on transposed\n");
    
    n=1;
    for(i=0; (i<nx); i++) {
      for(j=0; (j<ny); j++) {
	for(k=0; (k<nz); k++) {
	  fac        = ghat[nz-k-1][j][i];
	  fqqq[n]   *= fac;
	  fqqq[n+1] *= fac;
	  n+=2;
	}
      }
    }
    fprintf(stderr,"Doing %d backward 1D-FFTs\n",ntot);
    for(i=1; (i<nx); i++) {
      for(j=1; (j<ny); j++) {
	four1(fqqq+((i-1)*nx+j-1)*ny,nz,-11);
      }
    }
    */
    fprintf(stderr,"Doing backward transpose\n");
    transpose(log,cr,fqqq,nx,ny,nz,FALSE);
  }
  
  fprintf(stderr,"Doing %d backward 2D-FFTs\n",nx);
  for(i=0; (i<nx); i++) 
    fourn(fqqq+ntot*i,ndim,2,-1);
  ntot=nx*ny*nz*2;
  fac=2.0/(ntot);
  for(i=0; (i<=ntot); i++)
    fqqq[i]*=fac;
}

