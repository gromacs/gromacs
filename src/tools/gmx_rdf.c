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
#include "physics.h"
#include "index.h"
#include "smalloc.h"
#include "fftgrid.h"
#include "calcgrid.h"
#include "nrnb.h"
#include "coulomb.h"
#include "gstat.h"
#include "matio.h"

typedef struct
{
  char *label;
  int  elem,mass;
  real a[4], b[4], c;
} t_CM_table;

/*
 * 
 * f0[k] = c + [SUM a_i*EXP(-b_i*(k^2)) ]
 *             i=1,4
 */

const t_CM_table CM_t[] =
{

  { "H", 1,  1,   { 0.489918, 0.262003, 0.196767, 0.049879 },
    { 20.6593, 7.74039, 49.5519, 2.20159 },
    0.001305 },
  { "C", 6,  12, { 2.26069, 1.56165, 1.05075, 0.839259 },
    { 22.6907, 0.656665, 9.75618, 55.5949 },
    0.286977 },
  { "N", 7,  14,   { 12.2126, 3.13220, 2.01250, 1.16630 },     
    { 0.005700, 9.89330, 28.9975, 0.582600 },
    -11.529 },
  { "O", 8,  16, { 3.04850, 2.28680, 1.54630, 0.867000 },  
    { 13.2771, 5.70110, 0.323900, 32.9089 },
    0.250800 },
  { "Na", 11, 23,  { 3.25650, 3.93620, 1.39980, 1.00320 },       /*  Na 1+ */
    { 2.66710, 6.11530, 0.200100, 14.0390 }, 
    0.404000 }
};

#define NCMT asize(CM_t)

typedef struct
{
  int     n_angles;
  int     n_groups;
  double  lambda;
  double  energy;
  double  momentum;
  double  ref_k;
  double  **F;
  int     nSteps;
  int     total_n_atoms;
} structure_factor;

typedef struct
{
  rvec x;
  int  t;
} reduced_atom;

static void check_box_c(matrix box)
{
  if (fabs(box[ZZ][XX]) > GMX_REAL_EPS*box[ZZ][ZZ] ||
      fabs(box[ZZ][YY]) > GMX_REAL_EPS*box[ZZ][ZZ])
    gmx_fatal(FARGS,
	      "The last box vector is not parallel to the z-axis: %f %f %f",
	      box[ZZ][XX],box[ZZ][YY],box[ZZ][ZZ]);
}

static void do_rdf(char *fnNDX,char *fnTPS,char *fnTRX,
		   char *fnRDF,char *fnCNRDF, char *fnHQ,
		   bool bCM,bool bXY,bool bPBC,
		   real cutoff,real binwidth,real fade,int ng)
{
  FILE       *fp;
  int        status;
  char       outf1[STRLEN],outf2[STRLEN];
  char       title[STRLEN];
  int        g,natoms,i,j,k,nbin,j0,j1,n,nframes;
  int        **count;
  char       **grpname;
  int        *isize,isize_cm=0,nrdf=0,max_i;
  atom_id    **index,*index_cm=NULL;
#if (defined SIZEOF_LONG_LONG_INT) && (SIZEOF_LONG_LONG_INT >= 8)    
  long long int *sum;
#else
  double     *sum;
#endif
  real       t,rmax2,cut2,r,r2,invbinw,normfac;
  real       segvol,spherevol,prev_spherevol,**rdf;
  rvec       *x,xcom,dx,*x_i1,xi;
  real       *inv_segvol,invvol,invvol_sum,rho;
  bool       *bExcl,bTop,bNonSelfExcl;
  matrix     box,box_pbc;
  int        **npairs,*nself;
  atom_id    ix,jx,***pairs;
  t_topology top;
  t_block    *excl;
  t_pbc      pbc;

  excl=NULL;
  
  if (fnTPS) {
    bTop=read_tps_conf(fnTPS,title,&top,&x,NULL,box,TRUE);
    mk_single_top(&top);
    if (bTop && !bCM)
      /* get exclusions from topology */
      excl=&(top.blocks[ebEXCLS]);
  }
  snew(grpname,ng+1);
  snew(isize,ng+1);
  snew(index,ng+1);
  fprintf(stderr,"\nSelect a reference group and %d group%s\n",
	  ng,ng==1?"":"s");
  if (fnTPS)
    get_index(&top.atoms,fnNDX,ng+1,isize,index,grpname);
  else
    rd_index(fnNDX,ng+1,isize,index,grpname);
  
  natoms=read_first_x(&status,fnTRX,&t,&x,box);
  if ( !natoms )
    gmx_fatal(FARGS,"Could not read coordinates from statusfile\n");
  if (fnTPS)
    /* check with topology */
    if ( natoms > top.atoms.nr ) 
      gmx_fatal(FARGS,"Trajectory (%d atoms) does not match topology (%d atoms)",
		  natoms,top.atoms.nr);
  /* check with index groups */
  for (i=0; i<=ng; i++)
    for (j=0; j<isize[i]; j++)
      if ( index[i][j] >= natoms )
	gmx_fatal(FARGS,"Atom index (%d) in index group %s (%d atoms) larger "
		    "than number of atoms in trajectory (%d atoms)",
		    index[i][j],grpname[i],isize[i],natoms);
  
  if (bCM) {
    /* move index[0] to index_cm and make 'dummy' index[0] */
    isize_cm=isize[0];
    snew(index_cm,isize_cm);
    for(i=0; i<isize[0]; i++)
      index_cm[i]=index[0][i];
    isize[0]=1;
    index[0][0]=natoms;
    srenew(index[0],isize[0]);
    /* make space for center of mass */
    srenew(x,natoms+1);
  }
  
  /* initialize some handy things */
  copy_mat(box,box_pbc);
  if (bXY) {
    check_box_c(box);
    /* Make sure the z-height does not influence the cut-off */
    box_pbc[ZZ][ZZ] = 2*max(box[XX][XX],box[YY][YY]);
  }
  rmax2   = 0.99*0.99*max_cutoff2(box_pbc);
  nbin    = (int)(sqrt(rmax2) / binwidth) + 1;
  invbinw = 1.0 / binwidth;
  cut2   = sqr(cutoff);

  snew(count,ng);
  snew(pairs,ng);
  snew(npairs,ng);
  snew(nself,ng);

  snew(bExcl,natoms);
  max_i = 0;
  for(g=0; g<ng; g++) {
    if (isize[g+1] > max_i)
      max_i = isize[g+1];

    /* this is THE array */
    snew(count[g],nbin+1);
  
    /* make pairlist array for groups and exclusions */
    snew(pairs[g],isize[0]);
    snew(npairs[g],isize[0]);
    for(i=0; i<isize[0]; i++) {
      ix = index[0][i];
      for(j=0; j < natoms; j++)
	bExcl[j] = FALSE;
      /* exclusions? */
      if (excl)
	for( j = excl->index[ix]; j < excl->index[ix+1]; j++)
	  bExcl[excl->a[j]]=TRUE;
      k = 0;
      snew(pairs[g][i], isize[g+1]);
      bNonSelfExcl = FALSE;
      for(j=0; j<isize[g+1]; j++) {
	jx = index[g+1][j];
	if (!bExcl[jx])
	  pairs[g][i][k++]=jx;
	else if (ix != jx)
	  /* Check if we have exclusions other than self exclusions */
	  bNonSelfExcl = TRUE;
	if (ix == jx)
	  nself[g]++;
      }
      if (bNonSelfExcl) {
	npairs[g][i]=k;
	srenew(pairs[g][i],npairs[g][i]);
      } else {
	/* Save a LOT of memory and some cpu cycles */
	npairs[g][i]=-1;
	sfree(pairs[g][i]);
      }
    }
  }
  sfree(bExcl);

  snew(x_i1,max_i);
  nframes = 0;
  invvol_sum = 0;
  do {
    /* Must init pbc every step because of pressure coupling */
    copy_mat(box,box_pbc);
    if (bPBC) {
      rm_pbc(&top.idef,natoms,box,x,x);
      if (bXY) {
	check_box_c(box);
	box_pbc[ZZ][ZZ] = 2*max(box[XX][XX],box[YY][YY]);
      }
      set_pbc(&pbc,box_pbc);

      if (bXY)
	/* Set z-size to 1 so we get the surface iso the volume */
	box_pbc[ZZ][ZZ] = 1;
    }
    invvol = 1/det(box_pbc);
    invvol_sum += invvol;

    if (bCM) {
      /* calculate centre of mass */
      clear_rvec(xcom);
      for(i=0; (i < isize_cm); i++) {
	ix = index_cm[i];
	rvec_inc(xcom,x[ix]);
      }
      /* store it in the first 'group' */
      for(j=0; (j<DIM); j++)
	x[index[0][0]][j] = xcom[j] / isize_cm;
    }

    for(g=0; g<ng; g++) {
      /* Copy the indexed coordinates to a continuous array */
      for(i=0; i<isize[g+1]; i++)
	copy_rvec(x[index[g+1][i]],x_i1[i]);
    
      for(i=0; i<isize[0]; i++) {
	copy_rvec(x[index[0][i]],xi);
	if (npairs[g][i] >= 0)
	  /* Expensive loop, because of indexing */
	  for(j=0; j<npairs[g][i]; j++) {
	    jx=pairs[g][i][j];
	    if (bPBC)
	      pbc_dx(&pbc,xi,x[jx],dx);
	    else
	      rvec_sub(xi,x[jx],dx);
	      
	    if (bXY)
	      r2 = dx[XX]*dx[XX] + dx[YY]*dx[YY];
	    else 
	      r2=iprod(dx,dx);
	    if (r2>cut2 && r2<=rmax2)
	      count[g][(int)(sqrt(r2)*invbinw)]++;
	  }
	else
	  /* Cheaper loop, no exclusions */
	  for(j=0; j<isize[g+1]; j++) {
	    if (bPBC)
	      pbc_dx(&pbc,xi,x_i1[j],dx);
	    else
	      rvec_sub(xi,x_i1[j],dx);
	    if (bXY)
	      r2 = dx[XX]*dx[XX] + dx[YY]*dx[YY];
	    else 
	      r2=iprod(dx,dx);
	    if (r2>cut2 && r2<=rmax2)
	      count[g][(int)(sqrt(r2)*invbinw)]++;
	  }
      }
    }
    nframes++;
  } while (read_next_x(status,&t,natoms,x,box));
  fprintf(stderr,"\n");
  
  close_trj(status);
  
  sfree(x);
  
  /* Average volume */
  invvol = invvol_sum/nframes;

  /* Calculate volume of sphere segments or length of circle segments */
  snew(inv_segvol,nbin);
  prev_spherevol=0;
  for(i=0; (i<nbin); i++) {
    r = (i+1)*binwidth;
    if (bXY) {
      spherevol=M_PI*r*r;
    } else {
      spherevol=(4.0/3.0)*M_PI*r*r*r;
    }
    segvol=spherevol-prev_spherevol;
    inv_segvol[i]=1.0/segvol;
    prev_spherevol=spherevol;
  }
  
  snew(rdf,ng);
  for(g=0; g<ng; g++) {
    /* We have to normalize by dividing by the number of frames */
    normfac = 1.0/(nframes*invvol*(isize[0]*isize[g+1] - nself[g]));
      
    /* Do the normalization */
    nrdf = max(nbin-1,1+(2*fade/binwidth));
    snew(rdf[g],nrdf);
    for(i=0; (i<nbin-1); i++) {
      r = (i+0.5)*binwidth;
      if ((fade > 0) && (r >= fade))
	rdf[g][i] = 1+(count[g][i]*inv_segvol[i]*normfac-1)*exp(-16*sqr(r/fade-1));
      else
	rdf[g][i] = count[g][i]*inv_segvol[i]*normfac;
    }
    for( ; (i<nrdf); i++)
      rdf[g][i] = 1.0;
  }
    
  fp=xvgropen(fnRDF,"Radial Distribution","r","");
  if (ng==1)
    fprintf(fp,"@ subtitle \"%s-%s\"\n",grpname[0],grpname[1]);
  else {
    fprintf(fp,"@ subtitle \"reference %s\"\n",grpname[0]);
    xvgr_legend(fp,ng,grpname+1);
  }
  for(i=0; (i<nrdf); i++) {
    fprintf(fp,"%10g", (i+0.5)*binwidth);
    for(g=0; g<ng; g++)
      fprintf(fp," %10g",rdf[g][i]);
    fprintf(fp,"\n");
  }
  ffclose(fp);
  
  do_view(fnRDF,NULL);

  /* h(Q) function: fourier transform of rdf */  
  if (fnHQ) {
    int nhq = 401;
    real *hq,*integrand,Q;
    
    /* Get a better number density later! */
    rho = isize[1]*invvol;
    snew(hq,nhq);
    snew(integrand,nrdf);
    for(i=0; (i<nhq); i++) {
      Q = i*0.5;
      integrand[0] = 0;
      for(j=1; (j<nrdf); j++) {
	r = (j+0.5)*binwidth;
	integrand[j]  = (Q == 0) ? 1.0 : sin(Q*r)/(Q*r);
	integrand[j] *= 4.0*M_PI*rho*r*r*(rdf[0][j]-1.0);
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
    if (ng==1)
      fprintf(fp,"@ subtitle \"%s-%s\"\n",grpname[0],grpname[1]);
    else {
      fprintf(fp,"@ subtitle \"reference %s\"\n",grpname[0]);
      xvgr_legend(fp,ng,grpname+1);
    }
    snew(sum,ng);
    for(i=0; (i<nbin-1); i++) {
      fprintf(fp,"%10g",i*binwidth);
      for(g=0; g<ng; g++) {
	fprintf(fp," %10g",(real)((double)sum[g]*normfac));
	sum[g] += count[g][i];
      }
      fprintf(fp,"\n");
    }
    ffclose(fp);
    sfree(sum);
    
    do_view(fnCNRDF,NULL);
  }

  for(g=0; g<ng; g++)
    sfree(rdf[g]);
  sfree(rdf);
}

t_complex *** rc_tensor_allocation(int x, int y, int z)
{
  t_complex ***t;
  int i,j;
  
  snew(t,x);
  t = (t_complex ***)calloc(x,sizeof(t_complex**));
  if(!t) exit(fprintf(stderr,"\nallocation error"));
  t[0] = (t_complex **)calloc(x*y,sizeof(t_complex*));
  if(!t[0]) exit(fprintf(stderr,"\nallocation error"));
  t[0][0] = (t_complex *)calloc(x*y*z,sizeof(t_complex));
  if(!t[0][0]) exit(fprintf(stderr,"\nallocation error"));
  
  for( j = 1 ; j < y ; j++) 
    t[0][j] = t[0][j-1] + z;
  for( i = 1 ; i < x ; i++) {
    t[i] = t[i-1] + y;
    t[i][0] = t[i-1][0] + y*z;
    for( j = 1 ; j < y ; j++) 
      t[i][j] = t[i][j-1] + z;
  }
  return t;
}
    
int return_atom_type (char *name)
{
  typedef struct {
    char *name;
    int  nh;
  } t_united_h;
  t_united_h uh[] = {
    { "CH1", 1 }, { "CH2", 2 }, { "CH3", 3 }, 
    { "CS1", 1 }, { "CS2", 2 }, { "CS3", 3 }, 
    { "CP1", 1 }, { "CP2", 2 }, { "CP3", 3 }
  };
  int i;

  for(i=0; (i<asize(uh)); i++) 
    if (strcmp(name,uh[i].name) == 0)
      return NCMT-1+uh[i].nh;
      
  for(i=0; (i<NCMT); i++)
    if (strncmp (name, CM_t[i].label,strlen(CM_t[i].label)) == 0)
      return i;
  gmx_fatal(FARGS,"\nError: atom (%s) not in list (%d types checked)!\n", 
	    name,i);
  
  return 0;
}

double CMSF (int type,int nh,double lambda, double sin_theta)
/* 
 * return Cromer-Mann fit for the atomic scattering factor:
 * sin_theta is the sine of half the angle between incoming and scattered
 * vectors. See g_sq.h for a short description of CM fit.
 */
{
  int i;
  double tmp = 0.0, k2;
  
  /*
   *  united atoms case
   *  CH2 / CH3 groups  
   */
  if (nh > 0) {
    tmp = (CMSF (return_atom_type ("C"),0,lambda, sin_theta) +
	   nh*CMSF (return_atom_type ("H"),0,lambda, sin_theta));
  }
  /* all atom case */
  else {
    k2 = (sqr (sin_theta) / sqr (10.0 * lambda));
    tmp = CM_t[type].c;
    for (i = 0; (i < 4); i++)
      tmp += CM_t[type].a[i] * exp (-CM_t[type].b[i] * k2);
  }
  return tmp;
}

real **compute_scattering_factor_table (structure_factor * sf,int *nsftable)
{
/*
 *  this function build up a table of scattering factors for every atom
 *  type and for every scattering angle.
 */
    int i, j;
    double sin_theta,q,hc=1239.842;
    real ** sf_table;

    /* \hbar \omega \lambda = hc = 1239.842 eV * nm */
    sf->momentum = ((double) (2. * 1000.0 * M_PI * sf->energy) / hc);
    sf->lambda = hc / (1000.0 * sf->energy);
    fprintf (stderr, "\nwavelenght = %f nm\n", sf->lambda);
    *nsftable = NCMT+3;
    snew (sf_table,*nsftable);
    for (i = 0; (i < *nsftable); i++) {
	snew (sf_table[i], sf->n_angles);
	for (j = 0; j < sf->n_angles; j++) {
	    q = ((double) j * sf->ref_k);
	    /* theta is half the angle between incoming 
	       and scattered wavevectors. */
	    sin_theta = q / (2.0 * sf->momentum);
	    if (i < NCMT)
	      sf_table[i][j] = CMSF (i,0,sf->lambda, sin_theta);
	    else
	      sf_table[i][j] = CMSF (i,i-NCMT+1,sf->lambda, sin_theta);
	}
    }
    return sf_table;
}

int * create_indexed_atom_type (reduced_atom * atm, int size)
{
/* 
 * create an index of the atom types found in a  group
 * i.e.: for water index_atp[0]=type_number_of_O and 
 *                 index_atp[1]=type_number_of_H
 * 
 * the last element is set to 0 
 */
    int *index_atp, i, i_tmp, j;

    snew (index_atp, 1);
    i_tmp = 1;
    index_atp[0] = atm[0].t;
    for (i = 1; i < size; i++) {
	for (j = 0; j < i_tmp; j++)
	    if (atm[i].t == index_atp[j])
		break;
	if (j == i_tmp) {	/* i.e. no indexed atom type is  == to atm[i].t */
	    i_tmp++;
	    srenew (index_atp, i_tmp * sizeof (int));
	    index_atp[i_tmp - 1] = atm[i].t;
	}
    }
    i_tmp++;
    srenew (index_atp, i_tmp * sizeof (int));
    index_atp[i_tmp - 1] = 0;
    return index_atp;
}

void rearrange_atoms (reduced_atom * positions, t_trxframe *fr, atom_id * index,
		      int isize, t_topology * top, bool flag)
/* given the group's index, return the (continuous) array of atoms */
{
  int i;
  
  if (flag)
    for (i = 0; i < isize; i++)
      positions[i].t =
	return_atom_type (*(top->atoms.atomname[index[i]]));
  for (i = 0; i < isize; i++)
    copy_rvec (fr->x[index[i]], positions[i].x);
}


int atp_size (int *index_atp)
{
    int i = 0;

    while (index_atp[i])
	i++;
    return i;
}

void compute_structure_factor (structure_factor * sf, matrix box,
			       reduced_atom * red, int isize, real start_q,
			       real end_q, int group,real **sf_table)
{
    t_complex ***tmpSF;
    rvec k_factor;
    real kdotx, asf, kx, ky, kz, krr;
    int kr, maxkx, maxky, maxkz, i, j, k, p, *counter;


    k_factor[XX] = 2 * M_PI / box[XX][XX];
    k_factor[YY] = 2 * M_PI / box[YY][YY];
    k_factor[ZZ] = 2 * M_PI / box[ZZ][ZZ];

    maxkx = (int) rint (end_q / k_factor[XX]);
    maxky = (int) rint (end_q / k_factor[YY]);
    maxkz = (int) rint (end_q / k_factor[ZZ]);

    snew (counter, sf->n_angles);

    tmpSF = rc_tensor_allocation(maxkx,maxky,maxkz);
/*
 * The big loop...
 * compute real and imaginary part of the structure factor for every
 * (kx,ky,kz)) 
 */
    fprintf(stderr,"\n");
    for (i = 0; i < maxkx; i++) {
	fprintf (stderr,"\rdone %3.1f%%     ", (double)(100.0*(i+1))/maxkx);
	kx = i * k_factor[XX];
	for (j = 0; j < maxky; j++) {
	    ky = j * k_factor[YY];
	    for (k = 0; k < maxkz; k++)
		if (i != 0 || j != 0 || k != 0) {
		    kz = k * k_factor[ZZ];
		    krr = sqrt (sqr (kx) + sqr (ky) + sqr (kz));
		    if (krr >= start_q && krr <= end_q) {
			kr = (int) rint (krr/sf->ref_k);
			if (kr < sf->n_angles) {
			    counter[kr]++;  /* will be used for the copmutation 
					       of the average*/
			    for (p = 0; p < isize; p++) {
				    
				asf = sf_table[red[p].t][kr];

				kdotx = kx * red[p].x[XX] +
				    ky * red[p].x[YY] + kz * red[p].x[ZZ];
				
				tmpSF[i][j][k].re += cos (kdotx) * asf;
				tmpSF[i][j][k].im += sin (kdotx) * asf;
			    }
			}
		    }
		}
	}
    }				/* end loop on i */
/*
 *  compute the square modulus of the structure factor, averaging on the surface
 *  kx*kx + ky*ky + kz*kz = krr*krr 
 *  note that this is correct only for a (on the macroscopic scale)
 *  isotropic system. 
 */
    for (i = 0; i < maxkx; i++) {
	kx = i * k_factor[XX]; for (j = 0; j < maxky; j++) {
	    ky = j * k_factor[YY]; for (k = 0; k < maxkz; k++) {
		kz = k * k_factor[ZZ]; krr = sqrt (sqr (kx) + sqr (ky)
		+ sqr (kz)); if (krr >= start_q && krr <= end_q) {
		    kr = (int) rint (krr / sf->ref_k); if (kr <
		    sf->n_angles && counter[kr] != 0)
			sf->F[group][kr] +=
			    (sqr (tmpSF[i][j][k].re) +
			     sqr (tmpSF[i][j][k].im))/ counter[kr];
		}
	    }
	}
    } sfree (counter); free(tmpSF[0][0]); free(tmpSF[0]); free(tmpSF);
}

void save_data (structure_factor * sf, char *file, int ngrps, real start_q,
	   real end_q)
{

    FILE *fp;
    int i, g = 0;
    double *tmp, polarization_factor, A;

    fp = xvgropen (file, "Scattering Intensity", "q (1/nm)",
		   "Intensity (a.u.)");

    snew (tmp, ngrps);

    for (g = 0; g < ngrps; g++)
	for (i = 0; i < sf->n_angles; i++) {
/*
 *          theta is half the angle between incoming and scattered vectors.
 *          
 *          polar. fact. = 0.5*(1+cos^2(2*theta)) = 1 - 0.5 * sin^2(2*theta)
 *          
 *          sin(theta) = q/(2k) := A  ->  sin^2(theta) = 4*A^2 (1-A^2) ->
 *          -> 0.5*(1+cos^2(2*theta)) = 1 - 2 A^2 (1-A^2)
 */
	    A = (double) (i * sf->ref_k) / (2.0 * sf->momentum);
	    polarization_factor = 1 - 2.0 * sqr (A) * (1 - sqr (A));
	    sf->F[g][i] *= polarization_factor;
	}
    for (i = 0; i < sf->n_angles; i++) {
	if (i * sf->ref_k >= start_q && i * sf->ref_k <= end_q) {
	    fprintf (fp, "%10.5f  ", i * sf->ref_k);
	    for (g = 0; g < ngrps; g++)
               fprintf (fp, "  %10.5f ", (sf->F[g][i]) /( sf->total_n_atoms*
				                          sf->nSteps));   
	    fprintf (fp, "\n");
	}
    }
    ffclose (fp);
}

int
do_scattering_intensity (char* fnTPS, char* fnNDX, char* fnXVG, char *fnTRX,
		         real start_q,real end_q, real energy,int ng)
{
    int i,*isize,status,flags = TRX_READ_X,**index_atp;
    char **grpname,title[STRLEN];
    atom_id **index;
    t_topology top;
    t_trxframe fr;
    reduced_atom **red;
    structure_factor *sf;
    rvec *xtop;
    real **sf_table;
    int nsftable;
    matrix box;
    double r_tmp;

    snew (sf, 1);
    sf->energy = energy;
    /* Read the topology informations */
    
    read_tps_conf (fnTPS, title, &top, &xtop, NULL, box, TRUE);
    sfree (xtop);

    /* groups stuff... */
    snew (isize, ng);
    snew (index, ng);
    snew (grpname, ng);

    fprintf (stderr, "\nSelect %d group%s\n", ng,
	     ng == 1 ? "" : "s");
    if (fnTPS)
	get_index (&top.atoms, fnNDX, ng, isize, index, grpname);
    else
	rd_index (fnNDX, ng, isize, index, grpname);

    /* The first time we read data is a little special */
    read_first_frame (&status, fnTRX, &fr, flags);

    sf->total_n_atoms = fr.natoms;
    
    snew (red, ng);
    snew (index_atp, ng);

    r_tmp = max (box[XX][XX], box[YY][YY]);
    r_tmp = (double) max (box[ZZ][ZZ], r_tmp);

    sf->ref_k = (2.0 * M_PI) / (r_tmp);
    /* ref_k will be the reference momentum unit */
    sf->n_angles = (int) rint (end_q / sf->ref_k);

    snew (sf->F, ng);
    for (i = 0; i < ng; i++)
	snew (sf->F[i], sf->n_angles);
    for (i = 0; i < ng; i++) {
	snew (red[i], isize[i]);
	rearrange_atoms (red[i], &fr, index[i], isize[i], &top, TRUE);
	index_atp[i] = create_indexed_atom_type (red[i], isize[i]);
    }
    sf_table = compute_scattering_factor_table (sf,&nsftable);
    /* This is the main loop over frames */

    do {
	sf->nSteps++;
	for (i = 0; i < ng; i++) {
	    rearrange_atoms (red[i], &fr, index[i], isize[i], &top,FALSE);

	    compute_structure_factor (sf, box, red[i], isize[i],
				      start_q, end_q, i, sf_table);
	}
    }
    while (read_next_frame (status, &fr));

    save_data (sf, fnXVG, ng, start_q, end_q);

    return 0;
}

int gmx_rdf(int argc,char *argv[])
{
  static char *desc[] = {
    "The structure of liquids can be studied by either neutron or X-ray",
    "scattering. The most common way to describe liquid structure is by a",
    "radial distribution function. However, this is not easy to obtain from",
    "a scattering experiment.[PAR]",
    "g_rdf calculates radial distribution functions in different ways.",
    "The normal method is around a (set of) particle(s), the other method",
    "is around the center of mass of a set of particles.",
    "With both methods rdf's can also be calculated around axes parallel",
    "to the z-axis with option [TT]-xy[tt].[PAR]",
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
  static bool bCM=FALSE,bXY=FALSE,bPBC=TRUE;
  static real cutoff=0,binwidth=0.002,grid=0.05,fade=0.0,lambda=0.1,distance=10;
  static int  npixel=256,nlevel=20,ngroups=1;
  static real start_q=0.0, end_q=60.0, energy=12.0;
  t_pargs pa[] = {
    { "-bin",      FALSE, etREAL, {&binwidth},
      "Binwidth (nm)" },
    { "-com",      FALSE, etBOOL, {&bCM},
      "RDF with respect to the center of mass of first group" },
    { "-pbc",      FALSE, etBOOL, {&bPBC},
      "Use periodic boundary conditions for computing distances" },
    { "-xy",       FALSE, etBOOL, {&bXY},
      "Use only the x and y components of the distance" },
    { "-cut",      FALSE, etREAL, {&cutoff},
      "Shortest distance (nm) to be considered"},
    { "-ng",       FALSE, etINT, {&ngroups},
      "Number of secondary groups to compute RDFs around a central group" },
    { "-fade",     FALSE, etREAL, {&fade},
      "From this distance onwards the RDF is tranformed by g'(r) = 1 + [g(r)-1] exp(-(r/fade-1)^2 to make it go to 1 smoothly. If fade is 0.0 nothing is done." },
    { "-grid",     FALSE, etREAL, {&grid},
      "[HIDDEN]Grid spacing (in nm) for FFTs when computing structure factors" },
    { "-npixel",   FALSE, etINT,  {&npixel},
      "[HIDDEN]# pixels per edge of the square detector plate" },
    { "-nlevel",   FALSE, etINT,  {&nlevel},
      "Number of different colors in the diffraction image" },
    { "-distance", FALSE, etREAL, {&distance},
      "[HIDDEN]Distance (in cm) from the sample to the detector" },
    { "-wave",     FALSE, etREAL, {&lambda},
      "[HIDDEN]Wavelength for X-rays/Neutrons for scattering. 0.1 nm corresponds to roughly 12 keV" },
    
    {"-startq", FALSE, etREAL, {&start_q},
     "Starting q (1/nm) "},
    {"-endq", FALSE, etREAL, {&end_q},
     "Ending q (1/nm)"},
    {"-energy", FALSE, etREAL, {&energy},
     "Energy of the incoming X-ray (keV) "}
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
/*    { efXPM, "-image", "sq",  ffOPTWR }*/
  };
#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);

  fnTPS = ftp2fn_null(efTPS,NFILE,fnm);
  fnNDX = ftp2fn_null(efNDX,NFILE,fnm);
  /*  bSQ   = opt2bSet("-sq",NFILE,fnm) || opt2parg_bSet("-grid",NPA,pa); */
  bSQ   = opt2bSet("-sq",NFILE,fnm);
  bRDF  = opt2bSet("-o",NFILE,fnm) || !bSQ;
  
  if (bSQ) {
    if (!fnTPS)
      gmx_fatal(FARGS,"Need a tps file for calculating structure factors\n");
  }
  else {
    if (!fnTPS && !fnNDX)
      gmx_fatal(FARGS,"Neither index file nor topology file specified\n"
		  "             Nothing to do!");
  }
 
  if  (bSQ) 
   do_scattering_intensity(fnTPS,fnNDX,opt2fn("-sq",NFILE,fnm),ftp2fn(efTRX,NFILE,fnm),
		           start_q, end_q, energy, ngroups  );

  if (bRDF) 
    do_rdf(fnNDX,fnTPS,ftp2fn(efTRX,NFILE,fnm),
	   opt2fn("-o",NFILE,fnm),opt2fn_null("-cn",NFILE,fnm),
	   opt2fn_null("-hq",NFILE,fnm),
	   bCM,bXY,bPBC,cutoff,binwidth,fade,ngroups);

  thanx(stderr);
  
  return 0;
}
