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
static char *SRCID_pme_c = "$Id$";

#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "txtdump.h"
#include "vec.h"
#include "complex.h"
#include "smalloc.h"
#include "futil.h"
#include "shift_util.h"
#include "ewald_util.h"
#include "fftgrid.h"
#include "fatal.h"
#include "ewald.h"
#include "pme.h"
#include "pppm.h"
#include "network.h"
#include "physics.h"
#include "nrnb.h"

#define DFT_TOL 1e-7

typedef real *splinevec[DIM];

static void calc_idx(int natoms,rvec invh,
		     rvec x[],ivec idx[],int nx,int ny,int nz,
		     int nnx[],int nny[],int nnz[])
{
  int  i;
  int  *idxptr;
  real *xptr,ix,iy,iz;
 
  ix = invh[XX];
  iy = invh[YY];
  iz = invh[ZZ];
  for(i=0; (i<natoms); i++) {
    xptr   = x[i];
    idxptr = idx[i];
    idxptr[XX] = nnx[(int)(nx+xptr[XX]*ix)];
    idxptr[YY] = nny[(int)(ny+xptr[YY]*iy)];
    idxptr[ZZ] = nnz[(int)(nz+xptr[ZZ]*iz)];
  }
}

void sum_qgrid(t_commrec *cr,t_nsborder *nsb,t_fftgrid *grid,bool bForward)
{
  static bool bFirst=TRUE;
  static t_fft_r *tmp;
  int i;
  static int localsize;
  static int maxproc;

#ifdef USE_MPI
  if(bFirst) {
    localsize=grid->la12r*grid->pfft.local_nx;
    if(!grid->workspace)
      snew(tmp,localsize);
    maxproc=grid->nx/grid->pfft.local_nx;
  }
  /* NOTE: FFTW doesnt necessarily use all processors for the fft;
     * above I assume that the ones that do have equal amounts of data.
     * this is bad since its not guaranteed by fftw, but works for now...
     * This will be fixed in the next release.
     */
  bFirst=FALSE;
  if(grid->workspace)
    tmp=grid->workspace;
  if(bForward) { /* sum contributions to local grid */
    for(i=0;i<maxproc;i++) {
      MPI_Reduce(grid->ptr+i*localsize, /*ptr arithm.     */
		 tmp,localsize,      
		 GMX_MPI_REAL,MPI_SUM,i,MPI_COMM_WORLD);
    }
    if(cr->pid<maxproc)
      memcpy(grid->ptr+cr->pid*localsize,tmp,localsize*sizeof(t_fft_r));
  }
  else { /* distribute local grid to all processors */
    for(i=0;i<maxproc;i++)
      MPI_Bcast(grid->ptr+i*localsize, /* ptr arithm     */
		localsize,       
		GMX_MPI_REAL,i,MPI_COMM_WORLD);
  }
#else
  fatal_error(0,"Parallel grid summation requires MPI.\n");    
#endif
}

void spread_q_bsplines(t_fftgrid *grid,rvec invh,t_nsborder *nsb,
		       ivec idx[],real *charge,splinevec theta,int order,
		       int nnx[],int nny[],int nnz[])
{
  /* spread charges from home atoms to local grid */
  t_fft_r *ptr;
  int      i,j,k,n,*i0,*j0,*k0,*ii0,*jj0,*kk0,ithx,ithy,ithz;
  int      nx,ny,nz,la2,la12,xidx,yidx,zidx;
  int      start,nr,norder,norder1,*idxptr,ind0;
  real     valx,valxy,qn;
  real     *thx,*thy,*thz;
  
  clear_fftgrid(grid);
  unpack_fftgrid(grid,&nx,&ny,&nz,&la2,&la12,TRUE,&ptr);
  start = START(nsb);
  nr    = HOMENR(nsb);
  ii0   = nnx+nx+1-order;
  jj0   = nny+ny+1-order;
  kk0   = nnz+nz+1-order;
  thx   = theta[XX];
  thy   = theta[YY];
  thz   = theta[ZZ];
  
  for(n=0; (n<nr); n++) {
    qn     = charge[start+n];
    idxptr = idx[n];
    
    if (qn != 0) {
      xidx    = idxptr[XX];
      yidx    = idxptr[YY];
      zidx    = idxptr[ZZ];
      i0      = ii0+xidx; /* Pointer arithmetic */
      norder  = n*order;
      norder1 = norder+order;
      
      for(ithx=norder; (ithx<norder1); ithx++,i0++) {
	i    = *i0;
	j0   = jj0+yidx; /* Pointer arithmetic */
	valx = qn*thx[ithx];
	
	for(ithy=norder; (ithy<norder1); ithy++,j0++) {
	  j     = *j0;
	  k0    = kk0+zidx; /* Pointer arithmetic */
	  valxy = valx*thy[ithy];
	  ind0  = INDEX(i,j,0);
	  
	  for(ithz=norder; (ithz<norder1); ithz++,k0++) {
	    k = *k0;
	    ptr[ind0+k] += valxy*thz[ithz];
	  }
	}
      }
    }
  }
}

real solve_pme(t_fftgrid *grid,real ewaldcoeff,real vol,
	       splinevec bsp_mod,rvec invh,matrix vir,t_commrec *cr)
{
  /* do recip sum over local cells in grid */
  t_fft_c *ptr,*p0;
  int     nx,ny,nz,la2,la12;
  int     kx,ky,kz,idx,idx0,maxkx,maxky,maxkz,kystart,kyend;
  real    m2,mx,my,mz;
  real    factor=M_PI*M_PI/(ewaldcoeff*ewaldcoeff);
  real    ets2,struct2,vfactor,ets2vf;
  real    eterm,d1,d2,energy=0;
  real    denom;
  real    bx,by;
  real    invx,invy,invz;
  real    virxx=0,virxy=0,virxz=0,viryy=0,viryz=0,virzz=0;
  bool    bPar = PAR(cr);
  
  unpack_fftgrid(grid,&nx,&ny,&nz,&la2,&la12,FALSE,(t_fft_r **)&ptr);
  clear_mat(vir);
    
  invx  = invh[XX]/nx;
  invy  = invh[YY]/ny;
  invz  = invh[ZZ]/nz;
    
  maxkx = (nx+1)/2;
  maxky = (ny+1)/2;
  maxkz = nz/2+1;
    
  if (bPar) { /* transpose X & Y and only sum local cells */
#ifdef USE_MPI
    kystart = grid->pfft.local_y_start_after_transpose;
    kyend   = kystart+grid->pfft.local_ny_after_transpose;
    if (debug)
      fprintf(debug,"solve_pme: kystart = %d, kyend=%d\n",kystart,kyend);
#else
    fatal_error(0,"Parallel PME attempted without MPI");
#endif /* end of parallel case loop */
  }
  else {
    kystart = 0;
    kyend   = ny;
  }
  for(ky=kystart; (ky<kyend); ky++) {  /* our local cells */
    if(ky<maxky)
      my = ky*invy;
    else
      my = (ky-ny)*invy;
    by = M_PI*vol*bsp_mod[YY][ky];
    
    for(kx=0; (kx<nx); kx++) {
      if(kx < maxkx)
	mx = kx*invx;
      else
	mx = (kx-nx)*invx;
      bx = bsp_mod[XX][kx];
      if (bPar)
	p0 = ptr + INDEX(ky,kx,0); /* Pointer Arithmetic */
      else
	p0 = ptr + INDEX(kx,ky,0); /* Pointer Arithmetic */
      for(kz=0; (kz<maxkz); kz++,p0++)  {
	if ((kx==0) && (ky==0) && (kz==0))
	  continue;
	d1      = p0->re;
	d2      = p0->im;
	mz      = kz*invz;
	m2      = mx*mx+my*my+mz*mz;
	denom   = m2*bx*by*bsp_mod[ZZ][kz];
	eterm   = ONE_4PI_EPS0*exp(-factor*m2)/denom;
	p0->re  = d1*eterm;
	p0->im  = d2*eterm;
	
	struct2 = d1*d1+d2*d2;
	if ((kz > 0) && (kz < (nz+1)/2))
	  struct2*=2;
	ets2     = eterm*struct2;
	vfactor  = (factor*m2+1)*2.0/m2;
	energy  += ets2;
	
	ets2vf   = ets2*vfactor;
	virxx   += ets2vf*mx*mx-ets2;
	virxy   += ets2vf*mx*my;   
	virxz   += ets2vf*mx*mz;  
	viryy   += ets2vf*my*my-ets2;
	viryz   += ets2vf*my*mz;
	virzz   += ets2vf*mz*mz-ets2;
      }
    }
  }
    
  /* Update virial with local values */  
  vir[XX][XX] = virxx;
  vir[XX][YY] = virxy;
  vir[XX][ZZ] = virxz;
  vir[YY][XX] = virxy;
  vir[YY][YY] = viryy;
  vir[YY][ZZ] = viryz;
  vir[ZZ][XX] = virxz; 
  vir[ZZ][YY] = viryz;
  vir[ZZ][ZZ] = virzz;
  
  for(nx=0;nx<DIM;nx++)
    for(ny=0;ny<DIM;ny++)
      vir[nx][ny]*=0.25;
  /* this virial seems ok for isotropic scaling, but I'm
   * experiencing problems on semiisotropic membranes */
   
  /* this energy should be corrected for a charged system */
  return(0.5*energy);
}

void gather_f_bsplines(t_fftgrid *grid,rvec invh,t_nsborder *nsb,
		       ivec idx[],rvec f[],real *charge,splinevec theta,
		       splinevec dtheta,int order,
		       int nnx[],int nny[],int nnz[])
{
  /* sum forces for local particles */  
  int     i,j,k,n,*i0,*j0,*k0,*ii0,*jj0,*kk0,ithx,ithy,ithz;
  int     nx,ny,nz,la2,la12;
  t_fft_r *ptr;
  int     xidx,yidx,zidx;
  real    tx,ty,dx,dy,qn;
  real    fx,fy,fz,gval,tgz;
  real    *thx,*thy,*thz,*dthx,*dthy,*dthz;
  int     start,nr,sn,norder,norder1,*idxptr,ind0;
    
  start = START(nsb);
  nr    = HOMENR(nsb);
  unpack_fftgrid(grid,&nx,&ny,&nz,&la2,&la12,TRUE,&ptr);
 
  thx  = theta[XX];
  thy  = theta[YY];
  thz  = theta[ZZ];
  dthx = dtheta[XX];
  dthy = dtheta[YY];
  dthz = dtheta[ZZ];
  ii0  = nnx+nx+1-order;
  jj0  = nny+ny+1-order;
  kk0  = nnz+nz+1-order;
  
  for(n=0; (n<nr); n++) {
    sn     = start+n;
    qn     = charge[sn];
    idxptr = idx[n];
    
    if (qn != 0) {
      xidx = idxptr[XX];
      yidx = idxptr[YY];
      zidx = idxptr[ZZ];
      
      fx      = 0.0;
      fy      = 0.0;
      fz      = 0.0;
      i0      = ii0+xidx; /* Pointer arithemtic */
      norder  = n*order;
      norder1 = norder+order;
      for(ithx=norder; (ithx<norder1); ithx++,i0++) {
	i     = *i0;
	tx    = thx[ithx];
	dx    = dthx[ithx];
	j0    = jj0+yidx; /* Pointer arithemtic */
	
	for(ithy=norder; (ithy<norder1); ithy++,j0++) {
	  j     = *j0;
	  ty    = thy[ithy];
	  dy    = dthy[ithy];
	  k0    = kk0+zidx; /* Pointer arithemtic */
	  ind0  = INDEX(i,j,0);
	  
	  for(ithz=norder; (ithz<norder1); ithz++,k0++) {
	    k     = *k0;
	    gval  = ptr[ind0+k];
	    tgz   = thz[ithz]*gval;
	    fx   -= dx*ty*tgz;
	    fy   -= tx*dy*tgz;
	    fz   -= tx*ty*dthz[ithz]*gval;
	  }
	}
      }
      f[sn][XX]+=invh[XX]*fx*qn;
      f[sn][YY]+=invh[YY]*fy*qn;
      f[sn][ZZ]+=invh[ZZ]*fz*qn;
    }
  }
  /* Since the energy and not forces are interpolated
   * the net force might not be exactly zero.
   * This can be solved by also interpolating F, but
   * that comes at a cost.
   * A better hack is to remove the net force every
   * step, but that must be done at a higher level
   * since this routine doesn't see all atoms if running
   * in parallel. Don't know how important it is?  EL 990726
   */
}


void make_bsplines(splinevec theta,splinevec dtheta,int order,int nx,
		   int ny,int nz,rvec x[],real *charge,
		   t_nsborder *nsb,rvec invh)
{
  /* construct splines for local atoms */
  int  i,j,k,l;
  real dr,div,rcons;
  real *data,*ddata,*xptr;
  int  start,nr,sn;
  rvec nk;
  
  start  = START(nsb);
  nr     = HOMENR(nsb);
  nk[XX] = nx;
  nk[YY] = ny;
  nk[ZZ] = nz;
  rcons  = max(nx,max(ny,nz))+1;
  
  for(i=0; (i<nr); i++) {
    sn = start+i;
    if (charge[sn] != 0.0) {
      xptr = x[sn];
      for(j=0; (j<DIM); j++) {
	dr  = rcons+xptr[j]*invh[j];
	dr -= (int)dr;
	
	/* dr is relative offset from lower cell limit */
	data=&(theta[j][i*order]);
	data[order-1]=0;
	data[1]=dr;
	data[0]=1-dr;
		
	for(k=3; (k<order); k++) {
	  div=1.0/(k-1.0);    
	  data[k-1]=div*dr*data[k-2];
	  for(l=1; (l<(k-1)); l++)
	    data[k-l-1]=div*((dr+l)*data[k-l-2]+(k-l-dr)*
			     data[k-l-1]);
	  data[0]=div*(1-dr)*data[0];
	}
	/* differentiate */
	ddata    = &(dtheta[j][i*order]);
	ddata[0] = -data[0];
	for(k=1; (k<order); k++)
	  ddata[k]=data[k-1]-data[k];
		
	div=1.0/(order-1);
	data[order-1]=div*dr*data[order-2];
	for(l=1; (l<(order-1)); l++)
	  data[order-l-1]=div*((dr+l)*data[order-l-2]+
			       (order-l-dr)*data[order-l-1]);
	data[0]=div*(1-dr)*data[0]; 
      }
    }
  }
}

    
void make_dft_mod(real *mod,real *data,int ndata)
{
  int i,j;
  real sc,ss,arg;
    
  for(i=0;i<ndata;i++) {
    sc=ss=0;
    for(j=0;j<ndata;j++) {
      arg=(2.0*M_PI*i*j)/ndata;
      sc+=data[j]*cos(arg);
      ss+=data[j]*sin(arg);
    }
    mod[i]=sc*sc+ss*ss;
  }
  for(i=0;i<ndata;i++)
    if(mod[i]<1e-7)
      mod[i]=(mod[i-1]+mod[i+1])*0.5;
}



void make_bspline_moduli(splinevec bsp_mod,int nx,int ny,int nz,int order)
{
  int nmax=max(nx,max(ny,nz));
  real *data,*ddata,*bsp_data;
  int i,k,l;
  real div;
    
  snew(data,order);
  snew(ddata,order);
  snew(bsp_data,nmax);

  data[order-1]=0;
  data[1]=0;
  data[0]=1;
	    
  for(k=3;k<order;k++) {
    div=1.0/(k-1.0);
    data[k-1]=0;
    for(l=1;l<(k-1);l++)
      data[k-l-1]=div*(l*data[k-l-2]+(k-l)*data[k-l-1]);
    data[0]=div*data[0];
  }
  /* differentiate */
  ddata[0]=-data[0];
  for(k=1;k<order;k++)
    ddata[k]=data[k-1]-data[k];
  div=1.0/(order-1);
  data[order-1]=0;
  for(l=1;l<(order-1);l++)
    data[order-l-1]=div*(l*data[order-l-2]+(order-l)*data[order-l-1]);
  data[0]=div*data[0]; 

  for(i=0;i<nmax;i++)
    bsp_data[i]=0;
  for(i=1;i<=order;i++)
    bsp_data[i]=data[i-1];
    
  make_dft_mod(bsp_mod[XX],bsp_data,nx);
  make_dft_mod(bsp_mod[YY],bsp_data,ny);
  make_dft_mod(bsp_mod[ZZ],bsp_data,nz);

  sfree(data);
  sfree(ddata);
  sfree(bsp_data);
}

static t_fftgrid *grid=NULL;
static int nx,ny,nz;
static    splinevec theta;
static    splinevec dtheta;
static    splinevec bsp_mod;


void init_pme(FILE *log,t_commrec *cr,t_nsborder *nsb,t_inputrec *ir)
{
  int i;
    
  fprintf(log,"Will do PME sum in reciprocal space.\n");
  nx = ir->nkx;
  ny = ir->nky;
  nz = ir->nkz;
    
  if (PAR(cr) && cr->nprocs>1) {
    fprintf(log,"Parallelized PME sum used.\n");
    if(nx%(cr->nprocs)!=0)
      fprintf(log,"Warning: For load balance, "
	      "fourier_nx should be divisible by NPROCS\n");
  } 
 
  /* allocate space for things */
  snew(bsp_mod[XX],nx);
  snew(bsp_mod[YY],ny);
  snew(bsp_mod[ZZ],nz);
  for(i=0;i<DIM;i++) {
    snew(theta[i],ir->pme_order*HOMENR(nsb)); 
    snew(dtheta[i],ir->pme_order*HOMENR(nsb));
  }

  grid=mk_fftgrid(log,PAR(cr),nx,ny,nz,ir->bOptFFT);
  make_bspline_moduli(bsp_mod,nx,ny,nz,ir->pme_order);   
}

real do_pme(FILE *logfile,   bool bVerbose,
	    t_inputrec *ir,  rvec x[],
	    rvec f[],        real charge[],
	    rvec box,	     t_commrec *cr,
	    t_nsborder *nsb, t_nrnb *nrnb,    
	    matrix vir,      real ewaldcoeff)
{ 
  static ivec *idx=NULL;
  static int *nnx,*nny,*nnz;
  int  i,ntot,npme;
  rvec invh;
  real energy,vol;
   
  energy=0;
  calc_invh(box,nx,ny,nz,invh);
  
  vol=box[XX]*box[YY]*box[ZZ];

  /* Compute fftgrid index for all atoms, with help of some extra variables */
  if (!idx) {
    snew(idx,HOMENR(nsb));
    snew(nnx,3*nx);
    snew(nny,3*ny);
    snew(nnz,3*nz);
    for(i=0; (i<3*nx); i++)
      nnx[i] = i % nx;
    for(i=0; (i<3*ny); i++)
      nny[i] = i % ny;
    for(i=0; (i<3*nz); i++)
      nnz[i] = i % nz;
  }    
  calc_idx(HOMENR(nsb),invh,x+START(nsb),idx,nx,ny,nz,nnx,nny,nnz);
    
  /* make local bsplines  */
  make_bsplines(theta,dtheta,ir->pme_order,nx,ny,nz,x,charge,nsb,invh);
  
  /* put local atoms on grid. */
  spread_q_bsplines(grid,invh,nsb,idx,charge,theta,ir->pme_order,nnx,nny,nnz);
		    
  inc_nrnb(nrnb,eNR_SPREADQBSP,
	   ir->pme_order*ir->pme_order*ir->pme_order*HOMENR(nsb));
	   
  /* sum contributions to local grid from other processors */
  if (PAR(cr))
    sum_qgrid(cr,nsb,grid,TRUE);
  
  /* do 3d-fft */ 
  gmxfft3D(grid,FFTW_FORWARD,cr);
  
  /* solve in k-space for our local cells */
  energy=solve_pme(grid,ewaldcoeff,vol,bsp_mod,invh,vir,cr);
  inc_nrnb(nrnb,eNR_SOLVEPME,nx*ny*nz*0.5);

  /* do 3d-invfft */
  gmxfft3D(grid,FFTW_BACKWARD,cr);
   
  /* distribute local grid to all processors */
  if (PAR(cr))
    sum_qgrid(cr,nsb,grid,FALSE);
    
  /* interpolate forces for our local atoms */
  gather_f_bsplines(grid,invh,nsb,idx,f,charge,theta,dtheta,ir->pme_order,
		    nnx,nny,nnz);
  /* gather_f_bsplines(grid,invh,nsb,x,f,charge,theta,dtheta,ir->pme_order);*/

  inc_nrnb(nrnb,eNR_GATHERFBSP,
	   ir->pme_order*ir->pme_order*ir->pme_order*HOMENR(nsb));

  ntot  = grid->nxyz;  
  npme  = ntot*log((real)ntot)/(cr->nprocs*log(2.0));
  inc_nrnb(nrnb,eNR_FFT,2*npme);

  return energy;  
}
