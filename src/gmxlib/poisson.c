#include <stdio.h>
#include <math.h>
#include "smalloc.h"
#include "typedefs.h"
#include "lrutil.h"
#include "macros.h"
#include "fftgrid.h"
#include "vec.h"
#include "pppm.h"

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
	
static void calc_invh(rvec box,int nx,int ny,int nz,rvec invh)
{
  invh[XX] = nx/box[XX];
  invh[YY] = ny/box[YY];
  invh[ZZ] = nz/box[ZZ];
}

static void spread_q_poisson(FILE *log,bool bVerbose,
			     int natoms,rvec x[],real charge[],rvec box,
			     real r1,real rc,
			     t_fftgrid *grid,t_nrnb *nrnb)
{
  static bool bFirst = TRUE;
  static int  *nnx,*nny,*nnz;
  rvec   invh;
  real   qi,qt,qwt;
  rvec   gridpoint,dx;
  real   WXYZ[27],bhh,r,half=0.5;
  ivec   ixyz;
  int    i,iX,iY,iZ,index,ttt,m;
  int    jx,jy,jz,jcx,jcy,jcz;
  int    nxyz,ncellsx,ncellsy,ncellsz;
  int    nx,ny,nz,la1,la2,la12;
  int    xmin,xmax,ymin,ymax,zmin,zmax;
  t_fft_tp *ptr;
  
  unpack_fftgrid(grid,&nx,&ny,&nz,&la1,&la2,&la12,&ptr);
  
  calc_invh(box,nx,ny,nz,invh);

  if (bFirst) {
    fprintf(log,"Spreading Charges using spread function on %dx%dx%d grid\n",
	    nx,ny,nz);
    fprintf(log,"invh = %10g,%10g,%10g\n",invh[XX],invh[YY],invh[ZZ]);
  
    calc_nxyz(nx,ny,nz,&nnx,&nny,&nnz);
    
    bFirst = FALSE;
  }
  ncellsx=(rc/box[XX])*nx;
  ncellsy=(rc/box[YY])*ny;
  ncellsz=(rc/box[ZZ])*nz;
  
  for(i=0; (i<natoms); i++) {
    qi=charge[i];
    
    if (qi != 0.0) {
      /* Determine position of particle in box */
      for(m=0; (m<DIM); m++) {
	/* Put particle in the box... */  
	ttt = x[i][m]*invh[m];
	bhh = box[m]*invh[m];
	if (ttt < 0)
	  ttt += bhh;
	else if (ttt >= bhh)
	  ttt -= bhh;
      
	/* Calculate nearest grid point, Round */
	ixyz[m]    = ttt+half;
      }
      xmin = ixyz[XX] + nx - ncellsx;
      ymin = ixyz[YY] + ny - ncellsy;
      zmin = ixyz[ZZ] + nz - ncellsz;
      xmax = xmin+min(nx,2*ncellsx+1);
      ymax = ymin+min(ny,2*ncellsy+1);
      zmax = zmin+min(nz,2*ncellsz+1);
      
      for(jx=xmin; (jx<=xmax); jx++) {
	jcx = nnx[jx];
	gridpoint[XX] = (jcx*box[XX])/nx;
	for(jy=ymin; (jy<=ymax); jy++) {
	  jcy = nny[jy];
	  gridpoint[YY] = (jcy*box[YY])/ny;
	  for(jz=zmin; (jz<=zmax); jz++) {
	    jcz = nnz[jz]; 
	    gridpoint[ZZ] = (jcz*box[ZZ])/nz;
	    
	    rvec_sub(gridpoint,x[i],dx);
	    r = norm(dx);
	    GR_INC(ptr[INDEX(jcx,jcy,jcz)],qi*spreadfunction(r1,rc,r));
	  }
	}
      }
    }
  }
}

void solve_poisson(FILE *log,t_fftgrid *grid,bool bVerbose,t_nrnb *nrnb)
{
  static bool bFirst = TRUE;
  static int  *nnx,*nny,*nnz;
  int    i,j,k;
  int    nx,ny,nz,la1,la2,la12;
  t_fft_tp *ptr;
  
  unpack_fftgrid(grid,&nx,&ny,&nz,&la1,&la2,&la12,&ptr);

  if (bFirst) {
    fprintf(log,"Solving Poisson Equation on %dx%dx%d grid\n",
	    nx,ny,nz);
  
    calc_nxyz(nx,ny,nz,&nnx,&nny,&nnz);
    
    bFirst = FALSE;
  }

  for(i=0; (i<nx); i++)
    ;
}

real do_poisson(FILE *log,       bool bVerbose,
		t_inputrec *ir,  int natoms,
		rvec x[],        rvec f[],
		real charge[],   rvec box,
		real phi[],      t_commrec *cr,
		t_nrnb *nrnb)
{
  static  bool bFirst = TRUE;
  static  t_fftgrid *grid;
  static  int       niter;
  static  real      r1,rc;
  static  rvec      beta;
  
  const     real tol = 1e-5;
  int       i,m;
  real      ctot;
  real      aver,tot,ener;
  ivec      grids;
  rvec      spacing;
  
  ener = 0.0;
  
  if (bFirst) {
    niter = ir->niter;

    fprintf(log,"Will use Poisson Solver for long-range electrostatics\n");
    fprintf(log,"Grid size is %d x %d x %d\n",ir->nkx,ir->nky,ir->nkz);

    if ((ir->nkx < 4) || (ir->nky < 4) || (ir->nkz < 4)) 
      fatal_error(0,"Grid must be at least 4 points in all directions");
    
    grid = mk_fftgrid(ir->nkx,ir->nky,ir->nkz);
    
    r1 = ir->rshort;
    rc = ir->rlong;
    for(m=0; (m<DIM); m++)
      beta[m] = 4.0/3.0;
      
    bFirst = FALSE;
  }
  else {
    /* Make the grid empty */
    clear_fftgrid(grid);
    
    spread_q_poisson(log,bVerbose,natoms,x,charge,box,r1,rc,grid,nrnb);
    
    /* Second step: solving the poisson equation in Fourier space */
    solve_poisson(log,grid,bVerbose,nrnb);
    
    /* Third and last step: gather the forces, energies and potential
     * from the grid.
     */
    ener=gather_f(log,bVerbose,natoms,x,f,charge,box,phi,grid,beta,nrnb);
  }
  
  return ener;
}

