#include <stdio.h>
#include <math.h>
#include "poisson.h"
#include "nrnb.h"
#include "shift_util.h"
	
real ps_gather_inner(int JCXYZ[],real WXYZ[],int ixw[],int iyw[],int izw[],
		     real c1x,real c1y,real c1z,real c2x,real c2y,real c2z,
		     real qi,rvec f,real ***ptr)
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
    pi += weight * ptr[jcx0][jcy0][jcz0];
  
    /* Forces */
    fX += weight * ((c1x*(ptr[ixw[jcx-1]] [jcy0]       [jcz0] - 
			  ptr[ixw[jcx+1]] [jcy0]       [jcz0] )) +
		    (c2x*(ptr[ixw[jcx-2]] [jcy0]       [jcz0] - 
			  ptr[ixw[jcx+2]] [jcy0]       [jcz0] )));
    fY += weight * ((c1y*(ptr[jcx0]       [iyw[jcy-1]] [jcz0] -
			  ptr[jcx0]       [iyw[jcy+1]] [jcz0] ))  +
		    (c2y*(ptr[jcx0]       [iyw[jcy-2]] [jcz0] -
			  ptr[jcx0]       [iyw[jcy+2]] [jcz0] )));
    fZ += weight * ((c1z*(ptr[jcx0]       [jcy0]       [izw[jcz-1]] -
			  ptr[jcx0]       [jcy0]       [izw[jcz+1]] ))  +
		    (c2z*(ptr[jcx0]       [jcy0]       [izw[jcz-2]] -
			  ptr[jcx0]       [jcy0]       [izw[jcz+2]] )));
  }
  f[XX] += qi*fX;
  f[YY] += qi*fY;
  f[ZZ] += qi*fZ;
  
  return pi;
}

real ps_gather_f(FILE *log,bool bVerbose,
		 int natoms,rvec x[],rvec f[],real charge[],rvec box,
		 real pot[],t_PSgrid *grid,rvec beta,t_nrnb *nrnb)
{
  static bool bFirst=TRUE;
  static int  *nnx,*nny,*nnz;
  static int  JCXYZ[81];
  int    i,m;
  real   energy;
  real   qi,pi;
  ivec   ixyz;
  rvec   invh,h,c1,c2;
  real   WXYZ[27];
  real   c1x,c1y,c1z,c2x,c2y,c2z;
  int    ixw[7],iyw[7],izw[7];
  int    ll;
  int    nx,ny,nz;
  real   ***ptr;
  
  unpack_PSgrid(grid,&nx,&ny,&nz,&ptr);
  
  calc_invh_h(box,nx,ny,nz,invh,h);
  
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
     
    calc_weights(i,nx,ny,nz,x[i],box,invh,ixyz,WXYZ);

    for(ll=llim2; (ll<=ulim2); ll++) {
      ixw[ll-llim2] = nnx[ixyz[XX]+ll+nx];
      iyw[ll-llim2] = nny[ixyz[YY]+ll+ny];
      izw[ll-llim2] = nnz[ixyz[ZZ]+ll+nz];
    }
    
    qi      = charge[i];
    pi      = ps_gather_inner(JCXYZ,WXYZ,ixw,iyw,izw,
			      c1x,c1y,c1z,c2x,c2y,c2z,
			      qi,f[i],ptr);
    
    energy += pi*qi;
    pot[i]  = pi;
  }
  
  inc_nrnb(nrnb,eNR_GATHERF,27*natoms);
  inc_nrnb(nrnb,eNR_WEIGHTS,3*natoms);
  
  return energy*0.5;
}

