#include <stdio.h>
#include <math.h>
#include "smalloc.h"
#include "vec.h"
#include "nrnb.h"
#include "lrutil.h"
#include "poisson.h"

void spread_q_poisson(FILE *log,bool bVerbose,
		      int natoms,rvec x[],real charge[],rvec box,
		      real r1,real rc,t_PSgrid *grid,t_nrnb *nrnb,
		      int ntab,real sftab[],real tabspace)
{
  /* Spreading charges by convolution in real space of the spread function with   * the neighbouring grid points, (that is within the cut-off)
   */
  static bool bFirst = TRUE;
  static int  *nnx,*nny,*nnz,NCELLS,MAXCELLS;
  static ivec *cells=NULL;
  rvec   invh,h;
  real   qi,qt,qwt,dx2,dy2,dz2,r2,xi,yi,zi,sf,dr,invspace;
  rvec   gridpoint,dx;
  real   WXYZ[27],bhh,r,half=0.5,rc2;
  ivec   ixyz;
  int    i,j,k,iX,iY,iZ,index,ttt,m,n,nr;
  int    jx,jy,jz,jcx,jcy,jcz;
  int    nxyz,ncellsx,ncellsy,ncellsz;
  int    nx,ny,nz,la1,la2,la12;
  int    xmin,xmax,ymin,ymax,zmin,zmax;
  real   ***rho;
  
  unpack_PSgrid(grid,&nx,&ny,&nz,&rho);
  
  calc_invh_h(box,nx,ny,nz,invh,h);
  ncellsx=(rc/box[XX])*nx;
  ncellsy=(rc/box[YY])*ny;
  ncellsz=(rc/box[ZZ])*nz;

  if (bFirst) {
    fprintf(log,"Spreading Charges using spread function on %dx%dx%d grid\n",
	    nx,ny,nz);
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
  invspace = 1.0/tabspace;
  
  /* Has to be parallellized too! */
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
      xi = x[i][XX];
      yi = x[i][YY];
      zi = x[i][ZZ];
      
      for(n=0; (n<NCELLS); n++) {
	/* Compute cell number */
	jx  = ixyz[XX] + cells[n][XX];
	jy  = ixyz[YY] + cells[n][YY];
	jz  = ixyz[ZZ] + cells[n][ZZ];
	
	/* Compute distance from atom to grid point */
	dx2 = sqr(xi - jx*h[XX]);
	dy2 = sqr(yi - jy*h[YY]);
	dz2 = sqr(zi - jz*h[ZZ]);
	r2  = dx2+dy2+dz2;
	
	if (r2 < rc2) {
	  r  = sqrt(r2);
	  nr = r*invspace;
	  dr = r-nr;
	  sf = (tabspace-dr)*sftab[nr]+dr*sftab[nr+1];
	  
	  /* Do modulo to compute real grid number */
	  jcx = nnx[jx+nx];
	  jcy = nny[jy+ny];
	  jcz = nnz[jz+nz]; 
	  rho[jcx][jcy][jcz] += qi*sf/*spreadfunction(r1,rc,r)*/ ;
	}
      }
    }
  }
}

