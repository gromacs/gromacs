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
#include "pppm.h"
#include "gstat.h"

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
  real       *inv_segvol,vol,vol_sum;
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
  
  /* Calculate normalization factor and totals */
  sum = 0;
  for(i=0; (i<nbin-1); i++)
    sum += count[i];
  r = nbin*binwidth;
  normfac = (4.0/3.0)*M_PI * r*r*r / sum;

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
    real *hq,*integrand,Q,rho;
    
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
      hq[i] = print_and_integrate(debug,nrdf,binwidth,integrand,0);
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

static void extract_sq(t_fftgrid *fftgrid,int nbin,real factor,
		       real count[],rvec box)
{
  int nx,ny,nz,la2,la12;
  t_fft_c *ptr,*p0;
  int i,j,k,maxkx,maxky,maxkz,n,ind;
  real k1,z,k_max;
  rvec lll,kk;
  
  calc_lll(box,lll);
  k_max = nbin/factor;
  unpack_fftgrid(fftgrid,&nx,&ny,&nz,&la2,&la12,FALSE,(t_fft_r **)&ptr);
  /* This bit copied from pme.c */
  maxkx = (nx+1)/2;
  maxky = (ny+1)/2;
  maxkz = nz/2+1;
  for(i=0; (i<nx); i++)
    for(j=0; (j<ny); j++) {
      ind = INDEX(i,j,0);
      p0  = ptr + ind;
      for(k=0; (k<maxkz); k++,p0++) {
	if ((i==0) && (j==0) && (k==0))
	  continue;
	z   = sqrt(sqr(p0->re)+sqr(p0->im));
	calc_k(lll,i,j,k,nx,ny,nz,kk);
	k1  = norm(kk);
	ind = k1*factor;
	if (ind < nbin)
	  count[ind] += z;
	else
	  fprintf(stderr,"k (%g) > k_max (%g)\n",k1,k_max);
      }
    }
}

static void do_sq(char *fnNDX,char *fnTPS,char *fnTRX,char *fnSQ,
		  real grid)
{
  FILE       *fp;
  int        status;
  char       title[STRLEN],*aname;
  int        natoms,i,j,k,nbin,j0,j1,n,nframes;
  real       *count;
  char       *grpname;
  int        isize;
  atom_id    *index;
  real       t,k_max,factor,yfactor;
  rvec       *x,*xndx,box_size,kk,lll;
  real       *inv_segvol,*fj,max_spacing;
  bool       *bExcl,bTop;
  matrix     box;
  int        nx,ny,nz,nelectron;
  atom_id    ix,jx,**pairs;
  t_topology top;
  t_fftgrid  *fftgrid;
  t_nrnb     nrnb;
  
  bTop=read_tps_conf(fnTPS,title,&top,&x,NULL,box,TRUE);

  fprintf(stderr,"\nSelect group for structure factor computation:\n");
  get_index(&top.atoms,fnNDX,1,&isize,&index,&grpname);
  if (isize < top.atoms.nr)
    snew(xndx,isize);
  else
    xndx = x;
  
  natoms=read_first_x(&status,fnTRX,&t,&x,box);

  init_nrnb(&nrnb);
    
  if ( !natoms )
    fatal_error(0,"Could not read coordinates from statusfile\n");
  /* check with topology */
  if ( natoms > top.atoms.nr ) 
    fatal_error(0,"Trajectory (%d atoms) does not match topology (%d atoms)",
		natoms,top.atoms.nr);
	
  /* Atomic scattering factors */
  snew(fj,isize);
  nelectron = 0;
  for(i=0; (i<isize); i++) {
    aname = *(top.atoms.atomname[index[i]]);
    switch (aname[0]) {
    case 'H':
      fj[i] = 1;
      break;
    case 'C':
      fj[i] = 6;
      break;
    case 'N':
      fj[i] = 7;
      break;
    case 'O':
      fj[i] = 8;
      break;
    case 'F':
      fj[i] = 9;
      break;
    case 'S':
      fj[i] = 16;
      break;
    default:
      fprintf(stderr,"Warning: don't know number of electrons for atom %s\n",
	      aname);
      fj[i] = 1;
    }
    nelectron += fj[i];
  }

  nx = ny = nz = 0;
  max_spacing = calc_grid(box,grid,&nx,&ny,&nz,1);	

  fftgrid = mk_fftgrid(stdout,FALSE,nx,ny,nz,FALSE);
  
  /* Determine largest k vector length */
  for(i=0; (i<DIM); i++)
    box_size[i] = box[i][i];
  calc_lll(box_size,lll);
  calc_k(lll,nx/2,ny/2,nz/2,nx,ny,nz,kk);
  k_max = 1+norm(kk);
  
  /* this is the S(q) array */
  nbin = 200;
  snew(count,nbin+1);
  
  nframes = 0;
  do {
    /* Put the atoms with scattering factor on a grid. Misuses
     * an old routine from the PPPM code.
     */
    for(i=0; (i<DIM); i++)
      box_size[i] = box[i][i];

    /* Make the grid empty */
    clear_fftgrid(fftgrid);
  
    /* Put particles on a grid */
    if (xndx != x) 
      for(i=0; (i<isize); i++)
	copy_rvec(x[index[i]],xndx[i]);
    spread_q(stdout,TRUE,0,isize,xndx,fj,box_size,fftgrid,&nrnb);
 
    /* FFT the density */
    gmxfft3D(fftgrid,FFTW_FORWARD,NULL);  
    
    /* Extract the Sq function and sum it into the average array */
    factor = nbin/k_max;
    extract_sq(fftgrid,nbin,factor,count,box_size);
    
    nframes++;
  } while (read_next_x(status,&t,natoms,x,box));
  fprintf(stderr,"\n");
  
  close_trj(status);
  
  sfree(x);

  /* Normalize it ?? */  
  factor  = k_max*grid/(nbin);
  yfactor = nelectron*(1.0/nframes)*(1.0/fftgrid->nxyz);
  fp=xvgropen(fnSQ,"Structure Factor","q (/nm)","S(q)");
  for(i=0; i<nbin; i++)
    fprintf(fp,"%10g %10g\n", (i+0.5)*factor,count[i]*yfactor);
  ffclose(fp);
  
  do_view(fnSQ,NULL);
  
  done_fftgrid(fftgrid);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
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
  static real cutoff=0, binwidth=0.005, grid = 0.1, fade=0.0;
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
      "Grid spacing (in nm) for FFTs when computing structure factors" }
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
	  grid);
  
  if (bRDF) 
    do_rdf(fnNDX,fnTPS,ftp2fn(efTRX,NFILE,fnm),
	   opt2fn("-o",NFILE,fnm),opt2fn_null("-cn",NFILE,fnm),
	   opt2fn_null("-hq",NFILE,fnm),
	   bCM,cutoff,binwidth,fade);

  thanx(stderr);
  
  return 0;
}
