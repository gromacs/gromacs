#include <stdio.h>
#include <math.h>
#include "smalloc.h"
#include "vec.h"
#include "nrnb.h"
#include "lrutil.h"
#include "poisson.h"
#include "physics.h"

/* a dirty fix for a problem with the variable hz on an SP2 */
#undef hz

void spread_q_poisson(FILE *log,bool bVerbose,bool bCoulomb,
		      int natoms,rvec x[],real prop[],rvec box,
		      real rc,t_PSgrid *grid,t_nrnb *nrnb,
		      bool bOld,real r1)
{
  /* Spreading charges or C6 (or any other property)
   * by convolution in real space of the spread function with   
   * the neighbouring grid points, (that is within the cut-off rc).
   * bOld and r1 are for backwards compatibility and testing.
   */
  static bool bFirst = TRUE;
  static int  *nnx,*nny,*nnz,NCELLS,MAXCELLS;
  static ivec *cells=NULL;
  rvec   invh,h;
  real   qi,dx2,dy2,dz2,r2,xi,yi,zi,sf,hx,hy,hz;
  real   A,B;
  real   bhh,r,half=0.5,rc2,inveps0;
  ivec   ixyz;
  int    i,j,k,iX,iY,iZ,ttt,m,n;
  int    jx,jy,jz,jcx,jcy,jcz;
  int    ncellsx,ncellsy,ncellsz;
  int    nx,ny,nz;
  real   ***rho;
  
  unpack_PSgrid(grid,&nx,&ny,&nz,&rho);
  
  calc_invh_h(box,nx,ny,nz,invh,h);
  ncellsx=(rc/box[XX])*nx;
  ncellsy=(rc/box[YY])*ny;
  ncellsz=(rc/box[ZZ])*nz;

  if (bFirst) {
    fprintf(log,"Spreading %s using spread function on %dx%dx%d grid\n",
	    bCoulomb ? "charges" : "C6",nx,ny,nz);
    fprintf(log,"invh = %10g,%10g,%10g\n",invh[XX],invh[YY],invh[ZZ]);
    fprintf(log,"ncells = %d,%d,%d\n",ncellsx,ncellsy,ncellsz);
    calc_nxyz(nx,ny,nz,&nnx,&nny,&nnz);

    MAXCELLS = (2*ncellsx+2)*(2*ncellsy+2)*(2*ncellsz+2);
    snew(cells,MAXCELLS);
    rc2    = sqr(rc+0.5*norm(h));
    NCELLS = 0;
    for(i=-ncellsx-1; i<=ncellsx+1; i++) {
      dx2 = sqr(i*h[XX]);
      for(j=-ncellsy-1; j<=ncellsy+1; j++) {
	dy2 = sqr(j*h[YY]);
	if (dx2 + dy2 < rc2) {
	  for(k=-ncellsz-1; k<=ncellsz+1; k++) {
	    dz2 = sqr(k*h[ZZ]);
	    if (dx2+dy2+dz2 < rc2) {
	      cells[NCELLS][XX] = i;
	      cells[NCELLS][YY] = j;
	      cells[NCELLS][ZZ] = k;
	      NCELLS++;
	    }
	  }
	}
      }
    }
    fprintf(log,"There are %d cells (maximum was %d)\n",
	    NCELLS,MAXCELLS);
    
    bFirst = FALSE;
  }

  rc2      = rc*rc; 
  inveps0  = 1.0/EPSILON0;
  if (bCoulomb) {
    A = -7.5*ONE_4PI_EPS0*pow(rc,-5.0);
    B =  7.5*ONE_4PI_EPS0*pow(rc,-3.0);
  }
  else {
    A = 120*pow(rc,-10.0);
    B = -90*pow(rc,-8.0);
  }
  hx = h[XX];
  hy = h[YY];
  hz = h[ZZ];
  
  /* Has to be parallellized too! */
  for(i=0; (i<natoms); i++) {
    qi=prop[i];
    
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
      xi = x[i][XX];
      yi = x[i][YY];
      zi = x[i][ZZ];
      iX = ixyz[XX];
      iY = ixyz[YY];
      iZ = ixyz[ZZ];
      
      for(n=0; (n<NCELLS); n++) {
	/* Compute cell number */
	jx  = iX + cells[n][XX];
	jy  = iY + cells[n][YY];
	jz  = iZ + cells[n][ZZ];
	
	/* Compute distance from atom to grid point */
	dx2 = sqr(xi - jx*hx);
	dy2 = sqr(yi - jy*hy);
	dz2 = sqr(zi - jz*hz);
	r2  = dx2+dy2+dz2;
	
	if (r2 < rc2) {
	  if (bOld) {
	    r  = sqrt(r2);
	    sf = spreadfunction(r1,rc,r)*inveps0;
	  } 
	  else
	    sf  = A*r2+B;
	  
	  /* Do modulo to compute real grid number */
	  jcx = nnx[jx+nx];
	  jcy = nny[jy+ny];
	  jcz = nnz[jz+nz]; 
	  rho[jcx][jcy][jcz] += qi*sf;
	}
      }
    }
  }
}

