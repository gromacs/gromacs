/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include <string.h>
#include "statutil.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "fatal.h"
#include "vec.h"
#include "copyrite.h"
#include "futil.h"
#include "readinp.h"
#include "statutil.h"
#include "txtdump.h"
#include "gstat.h"
#include "xvgr.h"
#include "physics.h"
#include "pdbio.h"
#include "matio.h"

typedef struct{
  int Nx; /* x grid points in unit cell */
  int Ny; /* y grid points in unit cell */
  int Nz; /* z grid points in unit cell */
  int dmin[3]; /* starting point x,y,z*/
  int dmax[3]; /* ending point x,y,z */
  real cell[6]; /* usual cell parameters */
  real * ed; /* data */
} XplorMap;

static int lo_write_xplor(XplorMap * map,char * file)
{
  FILE * fp;
  int z,i,j,n;
  fp = fopen(file,"w");
  if(!fp)
    gmx_file(file);
  /* The REMARKS part is the worst part of the XPLOR format
     and may cause problems with some programs */
  fprintf(fp,"\n       2 !NTITLE\n") ;
  fprintf(fp," REMARKS Energy Landscape from GROMACS\n") ;
  fprintf(fp," REMARKS DATE: 2004-12-21 \n") ;
  fprintf(fp," %7d %7d %7d %7d %7d %7d %7d %7d %7d\n",
	  map->Nx, map->dmin[0], map->dmax[0],
	  map->Ny, map->dmin[1], map->dmax[1],
	  map->Nz, map->dmin[2], map->dmax[2]);
  fprintf(fp,"%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n",
	  map->cell[0],map->cell[1],map->cell[2],map->cell[3],map->cell[4],map->cell[5]);
  fprintf(fp, "ZYX\n") ;
  
  z = map->dmin[2];
  for (n = 0; n < map->Nz; n++, z++) {
    fprintf(fp, "%8d\n", z) ;
    for (i = 0; i < map->Nx*map->Ny; i += 6) {
      for (j = 0; j < 6; j++)
	if (i+j < map->Nx*map->Ny)
	  fprintf(fp, "%12.5E", map->ed[n*map->Nx*map->Ny+i+j]);
      fprintf(fp, "\n") ;
    }
  }
  fprintf(fp, "   -9999\n") ;
  fclose(fp) ;
  return 0;
}

static int write_xplor(char *file,real *data,int nx,int ny,int nz,
		       real dmin[],real dmax[])
{
  XplorMap *xm;
  
  snew(xm,1);
  xm->Nx = nx;
  xm->Ny = ny;
  xm->Nz = nz;
  xm->ed = data;
  xm->cell[0] = dmax[XX]-dmin[XX];
  xm->cell[1] = dmax[YY]-dmin[YY];
  xm->cell[2] = dmax[ZZ]-dmin[ZZ];
  xm->cell[3] = xm->cell[4] = xm->cell[5] = 90;
  clear_ivec(xm->dmin);
  xm->dmax[XX] = nx-1;
  xm->dmax[YY] = ny-1;

  clear_ivec(xm->dmin);
  xm->dmax[XX] = nx-1;
  xm->dmax[YY] = ny-1;

  clear_ivec(xm->dmin);
  xm->dmax[XX] = nx-1;
  xm->dmax[YY] = ny-1;

  clear_ivec(xm->dmin);
  xm->dmax[XX] = nx-1;
  xm->dmax[YY] = ny-1;
  xm->dmax[ZZ] = nz-1;
  return lo_write_xplor(xm,file);
}

static void normalize_p_e(int len,double *P,int *nbin,real *E,real pmin)
{
  int i;
  double Ptot=0;
  
  for(i=0; (i<len); i++) {
    Ptot += P[i];
    if (nbin[i] > 0)
      E[i] = E[i]/nbin[i];
  }
  printf("Ptot = %g\n",Ptot);
  for(i=0; (i<len); i++) {
    P[i] = P[i]/Ptot;
    /* Have to check for pmin after normalizing to prevent "stretching"
     * the energies.
     */
    if (P[i] < pmin)
      P[i] = 0;
  }
}

typedef struct {
  int index;
  real ener;
} t_minimum;

static int comp_minima(const void *a,const void *b)
{
  t_minimum *ma = (t_minimum *) a;
  t_minimum *mb = (t_minimum *) b;
  
  if (ma->ener < mb->ener)
    return -1;
  else if (ma->ener > mb->ener)
    return 1;
  else
    return 0;
}

static void pick_minima(char *logfile,int nbins,int ndim,real W[])
{
  FILE *fp;
  int  i,j,k,ijk,nmin;
  bool bMin;
  t_minimum *mm;
  
  snew(mm,gmx_nint(pow(nbins,ndim))+1);
  nmin = 0;
  fp = fopen(logfile,"w");
#define index2(x,y) (nbins*x+y)
#define index3(x,y,z) (nbins*(nbins*x+y)+z)
  for(i=0; (i<nbins); i++) {
    for(j=0; (j<nbins); j++) {
      if (ndim == 3) {
	for(k=0; (k<nbins); k++) {
	  ijk    = index3(i,j,k);
	  bMin   = (((i == 0)       || ((i > 0)       && 
					(W[ijk] < W[index3(i-1,j,k)]))) &&
		    ((i == nbins-1) || ((i < nbins-1) && 
					(W[ijk] < W[index3(i+1,j,k)]))) &&
		    ((j == 0)       || ((j > 0)       && 
					(W[ijk] < W[index3(i,j-1,k)]))) &&
		    ((j == nbins-1) || ((j < nbins-1) && 
					(W[ijk] < W[index3(i,j+1,k)]))) &&
		    ((k == 0)       || ((k > 0)       && 
					(W[ijk] < W[index3(i,j,k-1)]))) &&
		    ((k == nbins-1) || ((k < nbins-1) && 
					(W[ijk] < W[index3(i,j,k+1)]))));
	  if (bMin) {
	    fprintf(fp,"Minimum %d at index %6d energy %10.3f\n",
		    nmin,ijk,W[ijk]);
	    mm[nmin].index = ijk;
	    mm[nmin].ener  = W[ijk];
	    nmin++;
	  }
	}
      }
      else {
	ijk    = index2(i,j);
	bMin   = (((i == 0)       || ((i > 0)       && 
				      (W[ijk] < W[index2(i-1,j)]))) &&
		  ((i == nbins-1) || ((i < nbins-1) && 
				      (W[ijk] < W[index2(i+1,j)]))) &&
		  ((j == 0)       || ((j > 0)       && 
				      (W[ijk] < W[index2(i,j-1)]))) &&
		  ((j == nbins-1) || ((j < nbins-1) && 
				      (W[ijk] < W[index2(i,j+1)]))));
	if (bMin) {
	  fprintf(fp,"Minimum %d at index %6d energy %10.3f\n",
		  nmin,ijk,W[ijk]);
	  mm[nmin].index = ijk;
	  mm[nmin].ener  = W[ijk];
	  nmin++;
	}
      }
    }
  }
  qsort(mm,nmin,sizeof(mm[0]),comp_minima);
  fprintf(fp,"Minima sorted after energy\n");
  for(i=0; (i<nmin); i++) {
    fprintf(fp,"Minimum %d at index %6d energy %10.3f\n",
	    i,mm[i].index,mm[i].ener);
  }
  fclose(fp);
}

static void do_sham(char *fn,char *ndx,char *xpm,char *xpm2,
		    char *xpm3,char *xpm4,char *pdb,char *logf,
		    int n,int neig,real **eig,real **enerT,
		    int nmap,real *mapindex,real **map,
		    real Tref,int nbins,int nlevels,real pmin,
		    char *mname,bool bSham,
		    bool bXmin,real xmin,bool bXmax,real xmax,
		    bool bYmin,real ymin,bool bYmax,real ymax)
{
  FILE    *fp;
  real    *min_eig,*max_eig;
  real    *axis_x,*axis_y,*axis_z;
  real    *W,*E,**WW,**EE,*S,**SS,*M,**MM,*bE;
  rvec    xxx;
  char    *buf;
  double  *P,*bfac,efac,bref,Wmin,Wmax,Winf,Emin,Emax,Einf,Smin,Smax,Sinf,Mmin,Mmax,Minf;
  real    *delta;
  int     i,j,k,imin,len,index,d,*nbin,*bindex,bi;
  int     *nxyz;
  t_block *b;
  t_rgb   rlo  = { 0, 0, 0 };
  t_rgb   rhi  = { 1, 1, 1 };
  
  /* Determine extremes for the eigenvectors */
  snew(min_eig,neig);
  snew(max_eig,neig);
  snew(nxyz,neig);
  snew(bfac,neig);
  snew(delta,neig);
  for(i=0; (i<neig); i++) {
    for(j=0; (j<n); j++) {
      min_eig[i] = min(min_eig[i],eig[i][j]);
      max_eig[i] = max(max_eig[i],eig[i][j]);
      delta[i]  = (max_eig[i]-min_eig[i])/(2.0*nbins);
    }
  }
  /* Check for input constraints */
  if (neig == 2) {
    min_eig[0] = bXmin ? xmin : min_eig[0] - delta[0];
    max_eig[0] = bXmax ? xmax : max_eig[0] + delta[0];
    min_eig[1] = bYmin ? ymin : min_eig[1] - delta[1];
    max_eig[1] = bYmax ? ymax : max_eig[1] + delta[1];
    bfac[0]    = nbins/(max_eig[0]-min_eig[0]);
    bfac[1]    = nbins/(max_eig[1]-min_eig[1]);
  }
  else for(i=0; (i<neig); i++) {
    /* Add some extra space, half a bin on each side */
    max_eig[i] += delta[i];
    min_eig[i] -= delta[i];
    bfac[i]     = nbins/(max_eig[i]-min_eig[i]);
  }
  /* Do the binning */ 
  bref = 1/(BOLTZ*Tref);
  snew(bE,n);
  if (enerT) {
    Emin = 1e8;
    for(j=0; (j<n); j++) {
      bE[j] = (bref - 1/(BOLTZ*enerT[1][j]))*enerT[0][j];
      Emin  = min(Emin,bE[j]);
    }
  }
  else
    Emin = 0;
  len=1;
  for(i=0; (i<neig); i++) 
    len=len*nbins;
  printf("There are %d bins in the %d-dimensional histogram. Beta-Emin = %g\n",
	 len,neig,Emin);
  snew(P,len);
  snew(W,len);
  snew(E,len);
  snew(S,len);
  snew(nbin,len);
  snew(bindex,n);
  
  /* Loop over projections */
  for(j=0; (j<n); j++) {
    index = 0;
    /* Loop over dimensions */
    for(i=0; (i<neig); i++) {
      nxyz[i] = bfac[i]*(eig[i][j]-min_eig[i]);
      /* Compute index in 1-D array */
      d = 1;
      for(k=i; (k<neig-1); k++)
	d = d*nbins;
      index += nxyz[i]*d;
    }
    range_check(index,0,len);
    /* Compute the exponential factor */
    if (enerT)
      efac = exp(-bE[j]+Emin);
    else
      efac = 1;
    /* Update the probability */
    P[index] += efac;
    /* Update the energy */
    if (enerT)
      E[index] += enerT[0][j];
    /* Statistics: which "structure" in which bin */
    nbin[index]++;
    bindex[j]=index;
  }
  /* Normalize probability */
  normalize_p_e(len,P,nbin,E,pmin);
  /* Compute boundaries for the Free energy */
  Wmin = 1e8;
  imin = -1;
  Wmax = -1e8;
  /* Recompute Emin: it may have changed due to averaging */
  Emin = 1e8;
  Emax = -1e8;
  for(i=0; (i<len); i++)
    if (P[i] != 0) {
      W[i] = -BOLTZ*Tref*log(P[i]);
      if (W[i] < Wmin) {
	Wmin = W[i];
	imin = i;
      }
      Emin = min(E[i],Emin);
      Emax = max(E[i],Emax);
      Wmax = max(W[i],Wmax);
    } 
  Wmax -= Wmin;
  Winf = Wmax+1;
  Einf = Emax+1;
  Smin = Emin-Wmax;
  Smax = Emax-Smin;
  Sinf = Smax+1;
  /* Write out the free energy as a function of bin index */
  fp = ffopen(fn,"w");
  for(i=0; (i<len); i++)
    if (P[i] != 0) {
      W[i] -= Wmin;
      S[i] = E[i]-W[i]-Smin;
      fprintf(fp,"%5d  %10.5e  %10.5e  %10.5e\n",i,W[i],E[i],S[i]);
    }
    else {
      W[i] = Winf;
      E[i] = Einf;
      S[i] = Sinf;
    }
  fclose(fp);
  /* Organize the structures in the bins */
  snew(b,1);
  snew(b->index,len+1);
  snew(b->a,n);
  b->index[0] = 0;
  for(i=0; (i<len); i++) {
    b->index[i+1] = b->index[i]+nbin[i]; 
    nbin[i] = 0;
  }
  for(i=0; (i<n); i++) {
    bi = bindex[i];
    b->a[b->index[bi]+nbin[bi]] = i;
    nbin[bi]++;
  }
  /* Consistency check */
  for(i=0; (i<len); i++) {
    if (nbin[i] != (b->index[i+1] - b->index[i]))
      gmx_fatal(FARGS,"nbin[%d] = %d, should be %d",i,nbin[i],
		b->index[i+1] - b->index[i]);
  }
  /* Write the index file */
  fp = ffopen(ndx,"w");
  for(i=0; (i<len); i++) {
    if (nbin[i] > 0) {
      fprintf(fp,"[ %d ]\n",i);
      for(j=b->index[i]; (j<b->index[i+1]); j++)
	fprintf(fp,"%d\n",b->a[j]+1);
    }
  }  
  fclose(fp);
  snew(axis_x,nbins);
  snew(axis_y,nbins);
  snew(axis_z,nbins);
  snew(WW,nbins);
  snew(EE,nbins);
  snew(SS,nbins);
  for(i=0; (i<nbins); i++) {
    axis_x[i] = min_eig[XX]+i/bfac[XX];
    axis_y[i] = min_eig[YY]+i/bfac[YY];
    axis_z[i] = min_eig[ZZ]+i/bfac[ZZ];
  }
  if (neig == 2) {
    pick_minima(logf,nbins,2,W);
    /* Dump to XPM file */
    for(i=0; (i<nbins); i++) {
      WW[i] = &(W[i*nbins]);
      EE[i] = &(E[i*nbins]);
      SS[i] = &(S[i*nbins]);
    }
    fp = fopen(xpm,"w");
    write_xpm(fp,0,"Gibbs Energy Landscape","G (kJ/mol)","PC1","PC2",nbins,nbins,
	      axis_x,axis_y,WW,0,Winf,rlo,rhi,&nlevels);
    fclose(fp);
    fp = fopen(xpm2,"w");
    write_xpm(fp,0,"Enthalpy Landscape","H (kJ/mol)","PC1","PC2",nbins,nbins,
	      axis_x,axis_y,EE,Emin,Einf,rlo,rhi,&nlevels);
    fclose(fp);
    fp = fopen(xpm3,"w");
    write_xpm(fp,0,"Entropy Landscape","TDS (kJ/mol)","PC1","PC2",nbins,nbins,
	      axis_x,axis_y,SS,0,Sinf,rlo,rhi,&nlevels);
    fclose(fp);
    if (map) {
      snew(M,len);
      snew(MM,nbins);
      for(i=0; (i<nbins); i++) 
	MM[i] = &(M[i*nbins]);
      Mmin = 1e8;
      Mmax = -1e8;
      for(i=0; (i<nmap); i++) {
	Mmin = min(Mmin,map[0][i]);
	Mmax = max(Mmax,map[0][i]);
      }
      Minf = Mmax*1.05;
      for(i=0; (i<len); i++) 
	M[i] = Minf;
      for(i=0; (i<nmap); i++) {
	index = gmx_nint(mapindex[i]);
	if (P[index] != 0)
	  M[index] = map[0][i];
      }
      fp = fopen(xpm4,"w");
      write_xpm(fp,0,"Custom Landscape",mname,"PC1","PC2",nbins,nbins,
		axis_x,axis_y,MM,0,Minf,rlo,rhi,&nlevels);
      fclose(fp);
      sfree(MM);
      sfree(M);
    }
  }
  else if (neig == 3) {
    pick_minima(logf,nbins,3,W);
    /* Dump to PDB file */
    fp = fopen(pdb,"w");
    index = 0;
    for(i=0; (i<nbins); i++) {
      xxx[XX] = 3*(i+0.5-nbins/2);
      for(j=0; (j<nbins); j++) {
	xxx[YY] = 3*(j+0.5-nbins/2);
	for(k=0; (k<nbins); k++,index++) {
	  xxx[ZZ] = 3*(k+0.5-nbins/2);
	  if (P[index] > 0)
	    fprintf(fp,"%-6s%5u  %-4.4s%3.3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
		    "ATOM",(index+1) %10000,"H","H",(index+1)%10000,
		    xxx[XX],xxx[YY],xxx[ZZ],1.0,W[index]);
	}
      }
    }
    fclose(fp);
    write_xplor("out.xplor",W,nbins,nbins,nbins,min_eig,max_eig);
#define index3(x,y,z) (nbins*(nbins*x+y)+z)
    nxyz[XX] = imin/(nbins*nbins);
    nxyz[YY] = (imin-nxyz[XX]*nbins*nbins)/nbins;
    nxyz[ZZ] = imin % nbins;
    for(i=0; (i<nbins); i++) {
      snew(WW[i],nbins);
      for(j=0; (j<nbins); j++)
	WW[i][j] = W[index3(i,j,nxyz[ZZ])];
    }
    snew(buf,strlen(xpm)+4);
    sprintf(buf,"%s",xpm);
    sprintf(&buf[strlen(xpm)-4],"12.xpm");
    fp = fopen(buf,"w");
    write_xpm(fp,0,"Gibbs Energy Landscape","W (kJ/mol)","PC1","PC2",nbins,nbins,
	      axis_x,axis_y,WW,0,Winf,rlo,rhi,&nlevels);
    fclose(fp);
    for(i=0; (i<nbins); i++) {
      for(j=0; (j<nbins); j++)
	WW[i][j] = W[index3(i,nxyz[YY],j)];
    }
    sprintf(&buf[strlen(xpm)-4],"13.xpm");
    fp = fopen(buf,"w");
    write_xpm(fp,0,"SHAM Energy Landscape","kJ/mol","PC1","PC3",nbins,nbins,
	      axis_x,axis_z,WW,0,Winf,rlo,rhi,&nlevels);
    fclose(fp);
    for(i=0; (i<nbins); i++) {
      for(j=0; (j<nbins); j++)
	WW[i][j] = W[index3(nxyz[XX],i,j)];
    }
    sprintf(&buf[strlen(xpm)-4],"23.xpm");
    fp = fopen(buf,"w");
    write_xpm(fp,0,"SHAM Energy Landscape","kJ/mol","PC2","PC3",nbins,nbins,
	      axis_y,axis_z,WW,0,Winf,rlo,rhi,&nlevels);
    fclose(fp);
    sfree(buf);
  }
  sfree(WW);
}

static void ehisto(char *fh,int n,real **enerT)
{
  FILE *fp;
  int  i,j,k,nbin,blength;
  int  *bindex;
  real *T,bmin,bmax,bwidth;
  int  **histo;
  
  bmin =  1e8;
  bmax = -1e8;
  snew(bindex,n);
  snew(T,n);
  nbin = 0;
  for(j=1; (j<n); j++) {
    for(k=0; (k<nbin); k++) {
      if (T[k] == enerT[1][j]) {
	bindex[j] = k;
	break;
      }
    }
    if (k == nbin) {
      bindex[j] = nbin;
      T[nbin]   = enerT[1][j];
      nbin++;
    }
    bmin = min(enerT[0][j],bmin);
    bmax = max(enerT[0][j],bmax);
  }
  bwidth  = 1.0;
  blength = (bmax - bmin)/bwidth + 2;
  snew(histo,nbin);
  for(i=0; (i<nbin); i++) {
    snew(histo[i],blength);
  }
  for(j=0; (j<n); j++) {
    k = (enerT[0][j]-bmin)/bwidth;
    histo[bindex[j]][k]++;
  }
  fp = xvgropen(fh,"Energy distribution","E (kJ/mol)","");
  for(j=0; (j<blength); j++) {
    fprintf(fp,"%8.3f",bmin+j*bwidth);
    for(k=0; (k<nbin); k++) {
      fprintf(fp,"  %6d",histo[k][j]);
    }
    fprintf(fp,"\n");
  }  
  fclose(fp);
}

int gmx_sham(int argc,char *argv[])
{
  static char *desc[] = {
    "g_sham reads a number of xvg files and analyzes data sets.",
    "A line in the input file may start with a time",
    "(see option [TT]-time[tt]) and any number of y values may follow.",
    "Multiple sets can also be",
    "read when they are seperated by & (option [TT]-n[tt]),",
    "in this case only one y value is read from each line.",
    "All lines starting with # and @ are skipped."
  };
  static real tb=-1,te=-1,frac=0.5,filtlen=0,binwidth=0.1;
  static bool bHaveT=TRUE,bDer=FALSE,bSubAv=TRUE,bAverCorr=FALSE,bXYdy=FALSE;
  static bool bEESEF=FALSE,bEENLC=FALSE,bEeFitAc=FALSE,bPower=FALSE;
  static bool bShamEner=TRUE,bSham=TRUE; 
  static real Tref=298.15,pmin=0;
  static int  linelen=4096,nsets_in=1,nb_min=4,resol=10,nbin=32,nlevels=25;
  static real xmin=0,ymin=0,xmax=1,ymax=1,ttol=0;
  static char *mname="";
  t_pargs pa[] = {
    { "-linelen", FALSE, etINT, {&linelen},
      "HIDDENMaximum input line length" },
    { "-time",    FALSE, etBOOL, {&bHaveT},
      "Expect a time in the input" },
    { "-b",       FALSE, etREAL, {&tb},
      "First time to read from set" },
    { "-e",       FALSE, etREAL, {&te},
      "Last time to read from set" },
    { "-ttol",     FALSE, etREAL, {&ttol},
      "Tolerance on time in appropriate units (usually ps)" },
    { "-n",       FALSE, etINT, {&nsets_in},
      "Read # sets seperated by &" },
    { "-d",       FALSE, etBOOL, {&bDer},
	"Use the derivative" },
    { "-bw",      FALSE, etREAL, {&binwidth},
      "Binwidth for the distribution" },
    { "-sham",    FALSE, etBOOL, {&bSham},
      "Turn off energy weighting even if energies are given" },
    { "-tsham",   FALSE, etREAL, {&Tref},
      "Temperature for single histogram analysis" },
    { "-pmin",    FALSE, etREAL, {&pmin},
      "Minimum probability. Anything lower than this will be set to zero" },
    { "-nbins",   FALSE, etINT,  {&nbin},
      "Number of bins for single histogram analysis" },
    { "-nlevels", FALSE, etINT,  {&nlevels},
      "Number of levels for energy landscape from single histogram analysis" },
    { "-mname",   FALSE, etSTR,  {&mname},
      "Legend label for the custom landscape" },
    { "-xmin",    FALSE, etREAL, {&xmin},
      "Minimum for the X-axis in 2D landscape" },
    { "-xmax",    FALSE, etREAL, {&xmax},
      "Maximum for the X-axis in 2D landscape" },
    { "-ymin",    FALSE, etREAL, {&ymin},
      "Minimum for the Y-axis in 2D landscape" },
    { "-ymax",    FALSE, etREAL, {&ymax},
      "Maximum for the Y-axis in 2D landscape" }
  };
#define NPA asize(pa)

  FILE     *out;
  int      n,e_n,d_n,nlast,s,nset,e_nset,d_nset,i,j=0;
  real     **val,**et_val,**dt_val,*t,*e_t,e_dt,d_dt,*d_t,dt,tot,error;
  double   *av,*sig,cum1,cum2,cum3,cum4,db;
    
  t_filenm fnm[] = { 
    { efXVG, "-f",    "graph",    ffREAD   },
    { efXVG, "-ene",  "esham",    ffOPTRD  },
    { efXVG, "-dist", "ener",     ffOPTWR  },
    { efXVG, "-histo","edist",    ffOPTWR  },
    { efNDX, "-bin",  "bindex",   ffOPTWR  },
    { efXPM, "-ls",   "gibbs",    ffOPTWR  },
    { efXPM, "-lsh",  "enthalpy", ffOPTWR  },
    { efXPM, "-lss",  "entropy",  ffOPTWR  },
    { efXPM, "-map",  "map",      ffOPTWR  },
    { efPDB, "-ls3",  "gibbs3",   ffOPTWR  },
    { efXVG, "-mdata","mapdata",  ffOPTWR  },
    { efLOG, "-g",    "shamlog",  ffOPTWR  }
  }; 
#define NFILE asize(fnm) 

  int     npargs;

  npargs = asize(pa); 
  
  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_BE_NICE ,
		    NFILE,fnm,npargs,pa,asize(desc),desc,0,NULL); 

  val=read_val(opt2fn("-f",NFILE,fnm),bHaveT,
	       opt2parg_bSet("-b",npargs,pa),tb-ttol,
	       opt2parg_bSet("-e",npargs,pa),te+ttol,
	       nsets_in,&nset,&n,&dt,&t,linelen);
  printf("Read %d sets of %d points, dt = %g\n\n",nset,n,dt);
  
  if (opt2fn_null("-ene",NFILE,fnm) != NULL) {
    et_val=read_val(opt2fn("-ene",NFILE,fnm),bHaveT,
		    opt2parg_bSet("-b",npargs,pa),tb-ttol,
		    opt2parg_bSet("-e",npargs,pa),te+ttol,
		    nsets_in,&e_nset,&e_n,&e_dt,&e_t,linelen);
    if (e_nset != 2) 
      gmx_fatal(FARGS,"Can only handle one energy component and one T in %s",
		opt2fn("-ene",NFILE,fnm));
    if (e_n != n)
      gmx_fatal(FARGS,"Number of energies (%d) does not match number of entries (%d) in %s",e_n,n,opt2fn("-f",NFILE,fnm));
  }
  else 
    et_val = NULL;
    
  if (opt2fn_null("-mdata",NFILE,fnm) != NULL) {
    dt_val=read_val(opt2fn("-mdata",NFILE,fnm),bHaveT,
		    FALSE,tb,FALSE,te,
		    nsets_in,&d_nset,&d_n,&d_dt,&d_t,linelen);
    if (d_nset != 1)
      gmx_fatal(FARGS,"Can only handle one mapping data column in %s",
		opt2fn("-mdata",NFILE,fnm));
  }
  else
    dt_val = NULL;

  if (et_val)
    ehisto(opt2fn("-histo",NFILE,fnm),e_n,et_val);
  
  do_sham(opt2fn("-dist",NFILE,fnm),opt2fn("-bin",NFILE,fnm),
	  opt2fn("-ls",NFILE,fnm),opt2fn("-lsh",NFILE,fnm),
	  opt2fn("-lss",NFILE,fnm),opt2fn("-map",NFILE,fnm),
	  opt2fn("-ls3",NFILE,fnm),opt2fn("-g",NFILE,fnm),
	  n,nset,val,et_val,d_n,d_t,dt_val,Tref,nbin,nlevels,pmin,
	  mname,bSham,
	  opt2parg_bSet("-xmin",NPA,pa),xmin,
	  opt2parg_bSet("-xmax",NPA,pa),xmax,
	  opt2parg_bSet("-ymin",NPA,pa),ymin,
	  opt2parg_bSet("-ymax",NPA,pa),ymax);
  
  thanx(stderr);

  return 0;
}
  
