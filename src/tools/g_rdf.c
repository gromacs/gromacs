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
static char *SRCID_g_rdf_c = "$Id$";

#include <math.h>
#include <ctype.h>
#include "string2.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "rmpbc.h"
#include "xvgr.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "tpxio.h"
#include "rdgroup.h"
#include "smalloc.h"
#include "fftgrid.h"
#include "calcgrid.h"
#include "nrnb.h"
#include "shift_util.h"
#include "pme.h"
#include "gstat.h"
#include "matio.h"

static void do_rdf(char *fnNDX,char *fnTPS,char *fnTRX,
		   char *fnRDF,char *fnCNRDF, char *fnHQ,
		   bool bCM,real cutoff,real binwidth,real fade)
{
  FILE       *fp;
  int        status;
  char       outf1[STRLEN],outf2[STRLEN];
  char       title[STRLEN];
  int        natoms,i,j,k,nbin,j0,j1,n,nframes;
  int        *count;
  char       **grpname;
  int        isize[3],nrdf;
  atom_id    *index[3];
  unsigned long int sum;
  real       t,boxmin,hbox,hbox2,cut2,r,r2,invbinw,normfac;
  real       segvol,spherevol,prev_spherevol,*rdf;
  rvec       *x,xcom,dx,*x_i1,xi;
  real       *inv_segvol,vol,vol_sum,rho;
  bool       *bExcl,bTop,bNonSelfExcl;
  matrix     box;
  int        *npairs;
  atom_id    ix,jx,**pairs;
  t_topology top;
  t_block    *excl;
  excl=NULL;
  
  if (fnTPS) {
    bTop=read_tps_conf(fnTPS,title,&top,&x,NULL,box,TRUE);
    mk_single_top(&top);
    if (bTop && !bCM)
      /* get exclusions from topology */
      excl=&(top.atoms.excl);
  }
  snew(grpname,2);
  fprintf(stderr,"\nSelect groups for RDF computation:\n");
  if (fnTPS)
    get_index(&top.atoms,fnNDX,2,isize,index,grpname);
  else
    rd_index(fnNDX,2,isize,index,grpname);
  
  natoms=read_first_x(&status,fnTRX,&t,&x,box);
  if ( !natoms )
    fatal_error(0,"Could not read coordinates from statusfile\n");
  if (fnTPS)
    /* check with topology */
    if ( natoms > top.atoms.nr ) 
      fatal_error(0,"Trajectory (%d atoms) does not match topology (%d atoms)",
		  natoms,top.atoms.nr);
  /* check with index groups */
  for (i=0; i<2; i++)
    for (j=0; j<isize[i]; j++)
      if ( index[i][j] >= natoms )
	fatal_error(0,"Atom index (%d) in index group %s (%d atoms) larger "
		    "than number of atoms in trajectory (%d atoms)",
		    index[i][j],grpname[i],isize[i],natoms);
  
  if (bCM) {
    /* move index[0] to index[2] and make 'dummy' index[0] */
    isize[2]=isize[0];
    snew(index[2],isize[2]);
    for(i=0; i<isize[0]; i++)
      index[2][i]=index[0][i];
    isize[0]=1;
    index[0][0]=natoms;
    srenew(index[0],isize[0]);
    /* make space for center of mass */
    srenew(x,natoms+1);
  }
  
  /* initialize some handy things */
  boxmin = min( norm(box[XX]), min( norm(box[YY]), norm(box[ZZ]) ) );
  hbox   = boxmin / 2.0;
  nbin   = (int)(hbox / binwidth) + 1;
  invbinw = 1.0 / binwidth;
  hbox2  = sqr(hbox);
  cut2   = sqr(cutoff);
  
  /* this is THE array */
  snew(count,nbin+1);
  
  /* make pairlist array for groups and exclusions */
  snew(pairs,isize[0]);
  snew(npairs,isize[0]);
  snew(bExcl,natoms);
  for( i = 0; i < isize[0]; i++) {
    ix = index[0][i];
    for( j = 0; j < natoms; j++)
      bExcl[j] = FALSE;
    /* exclusions? */
    if (excl)
      for( j = excl->index[ix]; j < excl->index[ix+1]; j++)
	bExcl[excl->a[j]]=TRUE;
    k = 0;
    snew(pairs[i], isize[1]);
    bNonSelfExcl = FALSE;
    for( j = 0; j < isize[1]; j++) {
      jx = index[1][j];
      if (!bExcl[jx])
	pairs[i][k++]=jx;
      else
	/* Check if we have exclusions other than self exclusions */
	bNonSelfExcl = bNonSelfExcl || (ix != jx);
    }
    if (bNonSelfExcl) {
      npairs[i]=k;
      srenew(pairs[i],npairs[i]);
    } else {
      /* Save a LOT of memory and some cpu cycles */
      npairs[i]=-1;
      sfree(pairs[i]);
    }
  }
  sfree(bExcl);

  snew(x_i1,isize[1]);
  nframes = 0;
  vol_sum = 0;
  do {
    /* Must init pbc every step because of pressure coupling */
    init_pbc(box,FALSE);
    rm_pbc(&top.idef,natoms,box,x,x);
    
    vol = det(box);
    vol_sum += vol;
    
    if (bCM) {
      /* calculate centre of mass */
      clear_rvec(xcom);
      for(i=0; (i < isize[2]); i++) {
	ix = index[2][i];
	rvec_inc(xcom,x[ix]);
      }
      /* store it in the first 'group' */
      for(j=0; (j<DIM); j++)
	x[index[0][0]][j] = xcom[j] / isize[2];
    }

    /* Copy the indexed coordinates to a continuous array */
    for(i=0; i<isize[1]; i++)
      copy_rvec(x[index[1][i]],x_i1[i]);
    
    for(i=0; i<isize[0]; i++) {
      copy_rvec(x[index[0][i]],xi);
      if (npairs[i] >= 0)
	/* Expensive loop, because of indexing */
	for(j=0; j<npairs[i]; j++) {
	  jx=pairs[i][j];
	  pbc_dx(xi,x[jx],dx);
	  r2=iprod(dx,dx);
	  if (r2>cut2 && r2<=hbox2)
	    count[(int)(sqrt(r2)*invbinw)]++;
	}
      else
	/* Cheaper loop, no exclusions */
	for(j=0; j<isize[1]; j++) {
	  pbc_dx(xi,x_i1[j],dx);
	  r2=iprod(dx,dx);
	  if (r2>cut2 && r2<=hbox2)
	    count[(int)(sqrt(r2)*invbinw)]++;
	}
    }
    nframes++;
  } while (read_next_x(status,&t,natoms,x,box));
  fprintf(stderr,"\n");
  
  close_trj(status);
  
  sfree(x);
  
  /* Average volume */
  vol = vol_sum/nframes;
  
  /* Calculate volume of sphere segments */
  snew(inv_segvol,nbin);
  prev_spherevol=0;
  for(i=0; (i<nbin); i++) {
    r = (i+1)*binwidth;
    spherevol=(4.0/3.0)*M_PI*r*r*r;
    segvol=spherevol-prev_spherevol;
    inv_segvol[i]=1.0/segvol;
    prev_spherevol=spherevol;
  }
  
  /* We have to normalize by dividing by the number of frames */
  rho     = isize[1]/vol;
  normfac = 1.0/((rho*nframes)*isize[0]);
    
  /* Do the normalization */
  nrdf = max(nbin-1,1+(2*fade/binwidth));
  snew(rdf,nrdf);
  for(i=0; (i<nbin-1); i++) {
    r = (i+0.5)*binwidth;
    if ((fade > 0) && (r >= fade))
      rdf[i] = 1+(count[i]*inv_segvol[i]*normfac-1)*exp(-16*sqr(r/fade-1));
    else
      rdf[i] = count[i]*inv_segvol[i]*normfac;
  }
  for( ; (i<nrdf); i++)
    rdf[i] = 1.0;
    
  fp=xvgropen(fnRDF,"Radial Distribution","r","");
  fprintf(fp,"@ subtitle \"%s-%s\"\n",grpname[0],grpname[1]);
  for(i=0; (i<nrdf); i++)
    fprintf(fp,"%10g %10g\n", (i+0.5)*binwidth,rdf[i]);
  ffclose(fp);
  
  do_view(fnRDF,NULL);

  /* h(Q) function: fourier transform of rdf */  
  if (fnHQ) {
    int nhq = 401;
    real *hq,*integrand,Q;
    
    /* Get a better number density later! */
    rho = isize[1]/vol;
    snew(hq,nhq);
    snew(integrand,nrdf);
    for(i=0; (i<nhq); i++) {
      Q = i*0.5;
      integrand[0] = 0;
      for(j=1; (j<nrdf); j++) {
	r = (j+0.5)*binwidth;
	integrand[j]  = (Q == 0) ? 1.0 : sin(Q*r)/(Q*r);
	integrand[j] *= 4.0*M_PI*rho*r*r*(rdf[j]-1.0);
      }
      hq[i] = print_and_integrate(debug,nrdf,binwidth,integrand,NULL,0);
    }
    fp=xvgropen(fnHQ,"h(Q)","Q(/nm)","h(Q)");
    for(i=0; (i<nhq); i++) 
      fprintf(fp,"%10g %10g\n",i*0.5,hq[i]);
    ffclose(fp);
    do_view(fnHQ,NULL);
    sfree(hq);
    sfree(integrand);
  }
    
  if (fnCNRDF) {  
    normfac = 1.0/(isize[0]*nframes);
    fp=xvgropen(fnCNRDF,"Cumulative Number RDF","r","number");
    fprintf(fp,"@ subtitle \"%s-%s\"\n",grpname[0],grpname[1]);
    sum = 0;
    for(i=0; (i<nbin-1); i++) {
      fprintf(fp,"%10g %10g\n",i*binwidth,(real)((double)sum*normfac));
      sum += count[i];
    }
    ffclose(fp);
    
    do_view(fnCNRDF,NULL);
  }
  sfree(rdf);
}

typedef struct {
  int  ndata;
  real kkk;
  real intensity;
} t_xdata;

int comp_xdata(const void *a,const void *b)
{
  t_xdata *xa,*xb;
  real tmp;
  
  xa = (t_xdata *)a;
  xb = (t_xdata *)b;
  
  if (xa->ndata == 0)
    return 1;
  else if (xb->ndata == 0)
    return -1;
  else {
    tmp = xa->kkk - xb->kkk;
    if (tmp < 0)
      return -1;
    else if (tmp > 0)
      return 1;
    else return 0;
  }
}

static t_xdata *init_xdata(int nx,int ny)
{
  int     ix,iy,i,j,maxkx,maxky;
  t_xdata *data;
  
  maxkx = (nx+1)/2;
  maxky = (ny+1)/2;
  snew(data,nx*ny);
  for(i=0; (i<nx); i++) {
    for(j=0; (j<ny); j++) {
      ix = abs((i < maxkx) ? i : (i - nx)); 
      iy = abs((j < maxky) ? j : (j - ny)); 
      data[ix*ny+iy].kkk = sqrt(ix*ix+iy*iy);
    }
  }
  return data;
}

static void extract_sq(t_fftgrid *fftgrid,int nbin,real k_max,real lambda,
		       real count[],rvec box,int npixel,real *map[],
		       t_xdata data[])
{
  int nx,ny,nz,la2,la12;
  t_fft_c *ptr,*p0;
  int  i,j,k,maxkx,maxky,maxkz,n,ind,ix,iy;
  real k1,kxy2,kz2,k2,z,kxy,kxy_max,cos_theta2,ttt,factor;
  rvec lll,kk;
  
  /*calc_lll(box,lll);
    k_max   = nbin/factor;
    kxy_max = k_max/sqrt(3);*/
  unpack_fftgrid(fftgrid,&nx,&ny,&nz,&la2,&la12,FALSE,(t_fft_r **)&ptr);
  /* This bit copied from pme.c */
  maxkx = (nx+1)/2;
  maxky = (ny+1)/2;
  maxkz = nz/2+1;
  factor = nbin/k_max;
  for(i=0; (i<nx); i++) {
#define IDX(i,n)  (i<=n/2) ? (i) : (i-n)
    kk[XX] = IDX(i,nx);
    for(j=0; (j<ny); j++) {
      kk[YY] = IDX(j,ny);
      ind = INDEX(i,j,0);
      p0  = ptr + ind;
      for(k=0; (k<maxkz); k++,p0++) {
	if ((i==0) && (j==0) && (k==0))
	  continue;
	kk[ZZ] = IDX(k,nz);
	z   = sqrt(sqr(p0->re)+sqr(p0->im));
	kxy2 = sqr(kk[XX]) + sqr(kk[YY]);
	k2   = kxy2+sqr(kk[ZZ]);
	k1   = sqrt(k2);
	ind  = k1*factor;
	if (ind < nbin) {
	  /* According to:
	   * R. W. James (1962), 
	   * The Optical Principles of the Diffraction of X-Rays,
	   * Oxbow press, Woodbridge Connecticut
	   * the intensity is proportional to (eq. 9.10):
	   * I = C (1+cos^2 [2 theta])/2 FFT
	   * And since
	   * cos[2 theta] = cos^2[theta] - sin^2[theta] = 2 cos^2[theta] - 1 
	   * we can compute the prefactor straight from cos[theta]
	   */
	  cos_theta2  = kxy2/k2;
	  /*ttt         = z*0.5*(1+sqr(2*cos_theta2-1));*/
	  ttt         = z*0.5*(1+cos_theta2);
	  count[ind] += ttt;
	  ix = ((i < maxkx) ? i : (i - nx)); 
	  iy = ((j < maxky) ? j : (j - ny));
	  map[npixel/2+ix][npixel/2+iy] += ttt; 
	  data[abs(ix)*ny+abs(iy)].ndata     += 1;
	  data[abs(ix)*ny+abs(iy)].intensity += ttt;
	}
	else
	  fprintf(stderr,"k (%g) > k_max (%g)\n",k1,k_max);
      }
    }
  }
}

typedef struct {
  char *name;
  int  nelec;
} t_element;

static void do_sq(char *fnNDX,char *fnTPS,char *fnTRX,char *fnSQ,
		  char *fnXPM,real grid,real lambda,real distance,
		  int npixel,int nlevel)
{
  FILE       *fp;
  t_element  elem[] = { { "H", 1 }, { "C", 6 }, { "N", 7 }, { "O", 8 }, { "F", 9 }, { "S", 16 } };
#define NELEM asize(elem)
  int        status;
  char       title[STRLEN],*aname;
  int        natoms,i,j,k,nbin,j0,j1,n,nframes,pme_order;
  real       *count,**map;
  char       *grpname;
  int        isize;
  atom_id    *index;
  real       I0,C,t,k_max,factor,yfactor,segvol;
  rvec       *x,*xndx,box_size,kk,lll;
  real       fj0,*fj,max_spacing,r,lambda_1;
  bool       *bExcl,bTop;
  matrix     box;
  int        nx,ny,nz,nelectron;
  atom_id    ix,jx,**pairs;
  splinevec  *theta;
  t_topology top;
  t_fftgrid  *fftgrid;
  t_nrnb     nrnb;
  t_xdata    *data;
    
  bTop=read_tps_conf(fnTPS,title,&top,&x,NULL,box,TRUE);

  fprintf(stderr,"\nSelect group for structure factor computation:\n");
  get_index(&top.atoms,fnNDX,1,&isize,&index,&grpname);
  if (isize < top.atoms.nr)
    snew(xndx,isize);
  else
    xndx = x;
  natoms=read_first_x(&status,fnTRX,&t,&x,box);
  fprintf(stderr,"\n");
  
  init_nrnb(&nrnb);
    
  if ( !natoms )
    fatal_error(0,"Could not read coordinates from statusfile\n");
  /* check with topology */
  if ( natoms > top.atoms.nr ) 
    fatal_error(0,"Trajectory (%d atoms) does not match topology (%d atoms)",
		natoms,top.atoms.nr);
	
  /* Atomic scattering factors */
  snew(fj,isize);
  I0 = 0;
  nelectron = 0;
  for(i=0; (i<isize); i++) {
    aname = *(top.atoms.atomname[index[i]]);
    fj0 = 1;
    if (top.atoms.atom[i].ptype == eptAtom) {
      for(j=0; (j<NELEM); j++)
	if (aname[0] == elem[j].name[0]) {
	  fj0 = elem[j].nelec;
	  break;
	}
      if (j == NELEM)
	fprintf(stderr,"Warning: don't know number of electrons for atom %s\n",aname);
    }
    /* Correct for partial charge */
    fj[i] = fj0 - top.atoms.atom[index[i]].q;
    
    nelectron += fj[i];
    
    I0 += sqr(fj[i]);
  }
  if (debug) {
    /* Dump scattering factors */
    for(i=0; (i<isize); i++)
      fprintf(debug,"Atom %3s-%5d q = %10.4f  f = %10.4f\n",
	      *(top.atoms.atomname[index[i]]),index[i],
	      top.atoms.atom[index[i]].q,fj[i]);
  }

  /* Constant for scattering */
  C = sqr(1.0/(ELECTRONMASS_keV*KILO*ELECTRONVOLT*1e7*distance));
  fprintf(stderr,"C is %g\n",C);
  
  /* This bit is dimensionless */
  nx = ny = nz = 0;
  max_spacing = calc_grid(box,grid,&nx,&ny,&nz,1);	
  pme_order   = max(4,1+(0.2/grid));
  npixel      = max(nx,ny);
  data        = init_xdata(nx,ny);
  
  fprintf(stderr,"Largest grid spacing: %g nm, pme_order %d, %dx%d pixel on image\n",
	  max_spacing,pme_order,npixel,npixel);
  init_pme(stdout,NULL,nx,ny,nz,pme_order,isize,FALSE);
    
  /* Determine largest k vector length. */
  k_max = 1+sqrt(sqr(1+nx/2)+sqr(1+ny/2)+sqr(1+nz/2));

  /* this is the S(q) array */
  nbin = npixel;
  snew(count,nbin+1);
  snew(map,npixel);
  for(i=0; (i<npixel); i++)
    snew(map[i],npixel);
  
  nframes = 0;
  do {
    /* Put the atoms with scattering factor on a grid. Misuses
     * an old routine from the PPPM code.
     */
    for(j=0; (j<DIM); j++)
      box_size[j] = box[j][j];
    
    /* Scale coordinates to the wavelength */
    for(i=0; (i<isize); i++)
      copy_rvec(x[index[i]],xndx[i]);
      
    /* put local atoms on grid. */
    fftgrid = spread_on_grid(stdout,isize,pme_order,xndx,fj,box,FALSE);

    /* FFT the density */
    gmxfft3D(fftgrid,FFTW_FORWARD,NULL);  
    
    /* Extract the Sq function and sum it into the average array */
    extract_sq(fftgrid,nbin,k_max,lambda,count,box_size,npixel,map,data);
    
    nframes++;
  } while (read_next_x(status,&t,natoms,x,box));
  fprintf(stderr,"\n");
  
  close_trj(status);
  
  sfree(x);

  /* Normalize it ?? */  
  factor  = k_max/(nbin);
  yfactor = (1.0/nframes)/*(1.0/fftgrid->nxyz)*/;
  fp=xvgropen(fnSQ,"Structure Factor","q (1/nm)","S(q)");
  fprintf(fp,"@ subtitle \"Lambda = %g nm. Grid spacing = %g nm\"\n",
	  lambda,grid);
  factor *= lambda;
  for(i=0; i<nbin; i++) {
    r      = (i+0.5)*factor*2*M_PI;
    segvol = 4*M_PI*sqr(r)*factor;
    fprintf(fp,"%10g %10g\n",r,count[i]*yfactor/segvol);
  }
  ffclose(fp);
  
  do_view(fnSQ,NULL);

  if (fnXPM) {
    t_rgb rhi = { 0,0,0 }, rlo = { 1,1,1 };
    real *tx,*ty,hi,inv_nframes;
    int  maxkx,maxky;
    
    hi = 0;
    inv_nframes = 1.0/nframes;
    maxkx = (nx+1)/2;
    maxky = (ny+1)/2;
    snew(tx,npixel);
    snew(ty,npixel);
    for(i=0; (i<npixel); i++) {
      tx[i] = i-npixel/2;
      ty[i] = i-npixel/2;

      for(j=0; (j<npixel); j++) { 
	map[i][j] *= inv_nframes;
	hi         = max(hi,map[i][j]);
      }
    }
      
    fp = ffopen(fnXPM,"w");
    write_xpm(fp,"Diffraction Image","Intensity","kx","ky",
	      nbin,nbin,tx,ty,map,0,hi,rlo,rhi,&nlevel);
    fclose(fp);
    sfree(tx);
    sfree(ty);

    /* qsort(data,nx*ny,sizeof(data[0]),comp_xdata);    
       fp = ffopen("test.xvg","w");
       for(i=0; (i<nx*ny); i++) {
       if (data[i].ndata != 0) {
       fprintf(fp,"%10.3f  %10.3f\n",data[i].kkk,data[i].intensity/data[i].ndata);
       }
       }
       fclose(fp);
    */
  }
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "The structure of liquids can be studied by either neutron or X-ray",
    "scattering. The most common way to describe liquid structure is by a",
    "radial distribution function. However, this is not easy to obtain from",
    "a scattering experiment.[PAR]",
    "g_rdf calculates radial distribution functions in different ways.",
    "The normal method is around a (set of) particle(s), the other method",
    "is around the center of mass of a set of particles.[PAR]",
    "If a run input file is supplied ([TT]-s[tt]), exclusions defined",
    "in that file are taken into account when calculating the rdf.",
    "The option [TT]-cut[tt] is meant as an alternative way to avoid",
    "intramolecular peaks in the rdf plot.",
    "It is however better to supply a run input file with a higher number of",
    "exclusions. For eg. benzene a topology with nrexcl set to 5",
    "would eliminate all intramolecular contributions to the rdf.",
    "Note that all atoms in the selected groups are used, also the ones",
    "that don't have Lennard-Jones interactions.[PAR]",
    "Option [TT]-cn[tt] produces the cumulative number rdf.[PAR]"
    "To bridge the gap between theory and experiment structure factors can",
    "be computed (option [TT]-sq[tt]). The algorithm uses FFT, the grid"
    "spacing of which is determined by option [TT]-grid[tt]."
  };
  static bool bCM=FALSE;
  static real cutoff=0,binwidth=0.001,grid=0.05,fade=0.0,lambda=0.1,distance=10;
  static int  npixel=256,nlevel=20;
  t_pargs pa[] = {
    { "-bin",      FALSE, etREAL, {&binwidth},
      "Binwidth (nm)" },
    { "-com",      FALSE, etBOOL, {&bCM},
      "RDF with respect to the center of mass of first group" },
    { "-cut",      FALSE, etREAL, {&cutoff},
      "Shortest distance (nm) to be considered"},
    { "-fade",     FALSE, etREAL, {&fade},
      "From this distance onwards the RDF is tranformed by g'(r) = 1 + [g(r)-1] exp(-(r/fade-1)^2 to make it go to 1 smoothly. If fade is 0.0 nothing is done." },
    { "-grid",     FALSE, etREAL, {&grid},
      "Grid spacing (in nm) for FFTs when computing structure factors" },
    { "-npixel",   FALSE, etINT,  {&npixel},
      "[HIDDEN]# pixels per edge of the square detector plate" },
    { "-nlevel",   FALSE, etINT,  {&nlevel},
      "Number of different colors in the diffraction image" },
    { "-distance", FALSE, etREAL, {&distance},
      "[HIDDEN]Distance (in cm) from the sample to the detector" },
    { "-wave",     FALSE, etREAL, {&lambda},
      "Wavelength for X-rays/Neutrons for scattering. 0.1 nm corresponds to roughly 12 keV" }
  };
#define NPA asize(pa)
  char       *fnTPS,*fnNDX;
  bool       bSQ,bRDF;
  
  t_filenm   fnm[] = {
    { efTRX, "-f",  NULL,     ffREAD },
    { efTPS, NULL,  NULL,     ffOPTRD },
    { efNDX, NULL,  NULL,     ffOPTRD },
    { efXVG, "-o",  "rdf",    ffOPTWR },
    { efXVG, "-sq", "sq",     ffOPTWR },
    { efXVG, "-cn", "rdf_cn", ffOPTWR },
    { efXVG, "-hq", "hq",     ffOPTWR },
    { efXPM, "-image", "sq",  ffOPTWR }
  };
#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);

  fnTPS = ftp2fn_null(efTPS,NFILE,fnm);
  fnNDX = ftp2fn_null(efNDX,NFILE,fnm);
  bSQ   = opt2bSet("-sq",NFILE,fnm) || opt2parg_bSet("-grid",NPA,pa);
  bRDF  = opt2bSet("-o",NFILE,fnm) || !bSQ;
  
  if (bSQ) {
    if (!fnTPS)
      fatal_error(0,"Need a tps file for calculating structure factors\n");
  }
  else {
    if (!fnTPS && !fnNDX)
      fatal_error(0,"Neither index file nor topology file specified\n"
		  "             Nothing to do!");
  }
 
  if  (bSQ) 
    do_sq(fnNDX,fnTPS,ftp2fn(efTRX,NFILE,fnm),opt2fn("-sq",NFILE,fnm),
	  ftp2fn(efXPM,NFILE,fnm),grid,lambda,distance,npixel,nlevel);
  
  if (bRDF) 
    do_rdf(fnNDX,fnTPS,ftp2fn(efTRX,NFILE,fnm),
	   opt2fn("-o",NFILE,fnm),opt2fn_null("-cn",NFILE,fnm),
	   opt2fn_null("-hq",NFILE,fnm),
	   bCM,cutoff,binwidth,fade);

  thanx(stderr);
  
  return 0;
}
