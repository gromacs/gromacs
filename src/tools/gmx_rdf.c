/*
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

#include <histogram.h>
#include <trajana.h>

typedef struct {
  int        nr;
  int       *index;
} t_pairlist;

typedef struct {
  gmx_ana_selection_t  *refsel;
  gmx_ana_indexmap_t   *cmap;
  gmx_histogram_t     **count;
  double    *pcount;
  double     isize0;
  real       invvol_sum;
  real       rmax2,cut2;
  t_pairlist **pairs;
  bool                 bXY, bMassDen;
} t_rdfupdatedata;

static void check_box_c(matrix box)
{
  if (fabs(box[ZZ][XX]) > GMX_REAL_EPS*box[ZZ][ZZ] ||
      fabs(box[ZZ][YY]) > GMX_REAL_EPS*box[ZZ][ZZ])
    gmx_fatal(FARGS,
	      "The last box vector is not parallel to the z-axis: %f %f %f",
	      box[ZZ][XX],box[ZZ][YY],box[ZZ][ZZ]);
}

static t_pairlist **
create_excl_pairlists(gmx_ana_selection_t *refsel,
                      int nr, gmx_ana_selection_t *sel[],
                      t_blocka *excl, int natoms)
{
  t_pairlist **p;
  int          g, i, j, k;
  bool        *bExcl,bNonSelfExcl;
  atom_id      ix,jx;

  if (!excl)
    return 0;
  snew(p, nr-1);
  snew(bExcl, natoms);
  for (g = 0; g < nr; g++) {
    /* make pairlist array for groups and exclusions */
    snew(p[g],refsel->g->isize);
    for (i = 0; i < refsel->g->isize; i++) {
      ix = refsel->g->index[i];
      for (j = 0; j < natoms; j++)
        bExcl[j] = FALSE;
      for (j = excl->index[ix]; j < excl->index[ix+1]; j++)
        bExcl[excl->a[j]] = TRUE;
      k = 0;
      snew(p[g][i].index, sel[g]->g->isize);
      bNonSelfExcl = FALSE;
      for (j = 0; j < sel[g]->g->isize; j++) {
        jx = sel[g]->g->index[j];
        if (!bExcl[jx])
          p[g][i].index[k++] = jx;
        else if (ix != jx)
          /* Check if we have exclusions other than self exclusions */
          bNonSelfExcl = TRUE;
      }
      if (bNonSelfExcl) {
        p[g][i].nr = k;
        srenew(p[g][i].index, p[g][i].nr);
      } else {
        /* Save a LOT of memory and some cpu cycles */
        p[g][i].nr = -1;
        sfree(p[g][i].index);
      }
    }
  }
  sfree(bExcl);
  return p;
}

static int
update_rdfs(t_topology *top, t_trxframe *fr, t_pbc *pbc,
            int nr, gmx_ana_selection_t *sel[], void *data)
{
  t_rdfupdatedata *d = (t_rdfupdatedata *)data;
  matrix   box_pbc;
  rvec    *x;
  int      g, i, j, ii, isize0, isize_g;
  atom_id  jx;
  rvec     xi, dx;
  real     r2, r2ii;

  x = fr->x;
  copy_mat(fr->box,box_pbc);
  if (pbc) {
    if (d->bXY) {
      check_box_c(fr->box);
      clear_rvec(box_pbc[ZZ]);
    }
    set_pbc(pbc,fr->ePBC,box_pbc);
  }
  if (d->bXY)
    /* Set z-size to 1 so we get the surface iso the volume */
    box_pbc[ZZ][ZZ] = 1;
  d->invvol_sum += 1/det(box_pbc);


  if (d->cmap) {
    gmx_ana_indexmap_update(d->cmap, d->refsel->g, FALSE);
    isize0 = d->cmap->nr;
  } else
    isize0 = d->refsel->p.nr;

  d->isize0 += isize0;

  for(g=0; g<nr; g++) {
    for(i=0; i<isize0; i++) {
      if (d->cmap) {
        /* Special loop, since we need to determine the minimum distance
         * over all selected atoms in the reference molecule/residue.
         */
        isize_g = sel[g]->p.nr;
        for(j=0; j<isize_g; j++) {
          r2 = 1e30;
          /* Loop over the positions within the reference. */
          for (ii = d->cmap->mapb.index[i]; ii < d->cmap->mapb.index[i+1]; ++ii) {
            if (pbc)
              pbc_dx(pbc,d->refsel->p.x[d->refsel->g->index[ii]],sel[g]->p.x[j],dx);
            else
              rvec_sub(d->refsel->p.x[d->refsel->g->index[ii]],sel[g]->p.x[j],dx);
            if (d->bXY)
              r2ii = dx[XX]*dx[XX] + dx[YY]*dx[YY];
            else
              r2ii = iprod(dx,dx);
            if (r2ii < r2)
              r2 = r2ii;
          }
          if (r2>d->cut2 && r2<=d->rmax2) {
            gmx_histogram_increment(d->count[g], sqrt(r2));
          }
        }
        d->pcount[g] += isize_g;
      } else {
        copy_rvec(d->refsel->p.x[i],xi);
        if (d->pairs && d->pairs[g][i].nr >= 0) {
          /* Expensive loop, because of indexing */
          for(j=0; j<d->pairs[g][i].nr; j++) {
            jx=d->pairs[g][i].index[j];
            if (pbc)
              pbc_dx(pbc,xi,x[jx],dx);
            else
              rvec_sub(xi,x[jx],dx);

            if (d->bXY)
              r2 = dx[XX]*dx[XX] + dx[YY]*dx[YY];
            else
              r2=iprod(dx,dx);
            if (r2>d->cut2 && r2<=d->rmax2) {
              if (d->bMassDen)
                gmx_histogram_add(d->count[g], sqrt(r2), top->atoms.atom[jx].m);
              else
                gmx_histogram_increment(d->count[g], sqrt(r2));
            }
          }
          d->pcount[g] += d->pairs[g][i].nr;
        } else {
          /* Cheaper loop, no exclusions */
          isize_g = sel[g]->p.nr;
          for(j=0; j<isize_g; j++) {
            if (pbc)
              pbc_dx(pbc,xi,sel[g]->p.x[j],dx);
            else
              rvec_sub(xi,sel[g]->p.x[j],dx);
            if (d->bXY)
              r2 = dx[XX]*dx[XX] + dx[YY]*dx[YY];
            else
              r2=iprod(dx,dx);
            if (r2>d->cut2 && r2<=d->rmax2) {
              if (d->bMassDen)
                gmx_histogram_add(d->count[g], sqrt(r2), sel[g]->m[j]);
              else
                gmx_histogram_increment(d->count[g], sqrt(r2));
            }
          }
          d->pcount[g] += isize_g;
        }
      }
    }
    gmx_histogram_finish_frame(d->count[g]);
  }
  return 0;
}

static void do_rdf(gmx_ana_traj_t *trj,
		   char *fnRDF,char *fnCNRDF, char *fnHQ,
                   const char **closet,const char **normt,bool bXY,
                   real cutoff,real binwidth,real fade)
{
  FILE                 *fp;
  char                  gtitle[STRLEN];
  int                   g,i,j,k,nbin;
  bool                  bClose;
  bool                  bCM;
  bool                  bDyn;
  bool                  bAtom;
  real                  r, normfac;
  real                  segvol,spherevol,prev_spherevol;
  gmx_histogram_t     **rdf, **cum_rdf;
  real                  m;
  real                 *inv_segvol,invvol,rho;
  matrix                box_pbc;
  t_topology           *top;
  bool                  bTop;
  int                   ngrps;
  gmx_ana_selection_t **sel;
  char                **grpnames;
  t_trxframe           *fr;
  int                   nframes;
  t_rdfupdatedata       d;
  e_index_t             itype;

  bClose = (closet[0][0] != 'n');
  d.bXY = bXY;
  d.bMassDen = (normt[0][0] == 'd' && normt[0][4] == 'm');
  gmx_ana_get_refsel(trj, 0, &d.refsel);
  gmx_ana_get_nanagrps(trj, &ngrps);
  gmx_ana_get_anagrps(trj, &sel);
  gmx_ana_get_grpnames(trj, &grpnames);

  gmx_ana_get_topology(trj, d.bMassDen, &top, &bTop);
  
  if (!top && d.bMassDen)
    gmx_fatal(FARGS, "No masses available while mass weighting was requested");

  bCM   = (d.refsel->p.m.type == INDEX_ALL);
  bDyn  = d.refsel->bDynamic;
  bAtom = (d.refsel->p.m.type == INDEX_ATOM);
  for (g = 0; g < ngrps; ++g) {
    if (sel[g]->bDynamic)
      bDyn = TRUE;
    if (sel[g]->p.m.type != INDEX_ATOM)
      bAtom = FALSE;
  }

  if (bCM && normt[0][0] != 'n') {
    for (g = 0; g < ngrps; ++g)
      if (sel[g]->bDynamic) {
        if (!gmx_ana_selection_init_coverfrac(sel[g], CFRAC_SOLIDANGLE))
          fprintf(stderr, "warning: selection '%s' will not be normalized by the angular fraction\n",
                  sel[g]->name);
      }
  }

  if (bClose) {
    e_index_t ctype;
    if (d.refsel->p.m.type != INDEX_ATOM)
      gmx_fatal(FARGS, "The first selection should give atom positions with -surf");
    snew(d.cmap, 1);
    if (closet[0][0] == 'r')
      ctype = INDEX_RES;
    else if (closet[0][0] == 'm')
      ctype = INDEX_MOL;
    else
      ctype = INDEX_ALL;
    gmx_ana_indexmap_init(d.cmap, d.refsel->g, top, ctype);
  } else
    d.cmap = NULL;
  
  gmx_ana_get_first_frame(trj, &fr);
  
  /* initialize some handy things */
  copy_mat(fr->box,box_pbc);
  if (bXY) {
    check_box_c(fr->box);
    /* Make sure the z-height does not influence the cut-off */
    box_pbc[ZZ][ZZ] = 2*max(fr->box[XX][XX], fr->box[YY][YY]);
  }
  if (gmx_ana_has_pbc(trj))
    d.rmax2   = 0.99*0.99*max_cutoff2(bXY ? epbcXY : epbcXYZ,box_pbc);
  else
    d.rmax2   = sqr(3*max(fr->box[XX][XX],max(fr->box[YY][YY],fr->box[ZZ][ZZ])));
  if (debug)
    fprintf(debug,"rmax2 = %g\n",d.rmax2);

  /* We use the double amount of bins, so we can correctly
   * write the rdf and rdf_cn output at i*binwidth values.
   */
  nbin     = (int)(sqrt(d.rmax2) / binwidth);
  d.cut2   = sqr(cutoff);

  snew(d.count, ngrps);
  for (g = 0; g < ngrps; ++g) {
    gmx_histogram_create(&d.count[g], d.bMassDen ? HIST_WEIGHT : HIST_SIMPLE,
                         nbin*2);
    gmx_histogram_set_binwidth(d.count[g], 0.0, binwidth/2);
  }
  snew(d.pcount, ngrps);
  /* We can only have exclusions with atomic rdfs */
  if (!bAtom || !top || !bTop || bClose || bDyn)
    d.pairs = NULL;
  else
    d.pairs = create_excl_pairlists(d.refsel, ngrps, sel, &(top->excls),
                                    fr->natoms);

  d.isize0 = 0;
  d.invvol_sum = 0;

  /* Do the analysis */
  gmx_ana_do(trj, 0, &update_rdfs, &d);
  gmx_ana_get_nframes(trj, &nframes);
  
  /* Average volume */
  invvol = d.invvol_sum/nframes;
  d.isize0 /= nframes;
  for (g = 0; g < ngrps; ++g)
    d.pcount[g] /= nframes;
  for (g = 0; g < ngrps; ++g)
    gmx_histogram_finish(d.count[g]);

  for (g = 0; g < ngrps; ++g)
    if (sel[g]->cfractype != CFRAC_NONE)
      fprintf(stderr, "Average angular fraction covered by selection '%s': %.1f%%\n",
              sel[g]->name, sel[g]->avecfrac * 100);

  /* Calculate volume of sphere segments or length of circle segments */
  snew(inv_segvol,nbin);
  prev_spherevol=0;
  for (i = 0; i < nbin; ++i) {
    r = (i + 0.5)*binwidth;
    if (bXY) {
      spherevol=M_PI*r*r;
    } else {
      spherevol=(4.0/3.0)*M_PI*r*r*r;
    }
    segvol=spherevol-prev_spherevol;
    inv_segvol[i]=1.0/segvol;
    prev_spherevol=spherevol;
  }
  
  snew(rdf, ngrps);
  for (g = 0; g < ngrps; ++g) {
    gmx_histogram_resample_dblbw(&rdf[g], d.count[g], TRUE);
    /* Normalize */
    if (normt[0][0] == 'n') {
      gmx_histogram_scale(rdf[g], 1.0/(binwidth*d.isize0));
    } else {
      gmx_histogram_scale_vec(rdf[g], inv_segvol);
      if (normt[0][0] == 'd')
        normfac = 1.0/(d.isize0) * (d.bMassDen ? AMU/(NANO*NANO*NANO) : 1);
      else
        normfac = 1.0/(invvol*d.pcount[g]);
      if (sel[g]->cfractype != CFRAC_NONE)
        normfac /= sel[g]->avecfrac;
      gmx_histogram_scale(rdf[g], normfac);
    }
    /* Do the fade */
    if (fade > 0 && normt[0][0] == 'r') {
      double *hist = gmx_histogram_get_values(rdf[g]);
      for (i = 0; i < nbin; ++i) {
        r = i*binwidth;
        if (r >= fade)
          hist[i] = 1 + (hist[i]-1)*exp(-16*sqr(r/fade-1));
      }
    }
  }

  /* Check if the positions have a common type */
  itype = sel[g]->p.m.type;
  for (g = 1; g < ngrps; ++g) {
    if (itype != sel[g]->p.m.type) {
      itype = INDEX_UNKNOWN;
      break;
    }
  }
  if (itype == INDEX_RES || itype == INDEX_MOL) {
    sprintf(gtitle,"Radial distribution of %s %s",
	    itype == INDEX_MOL ? "molecule" : "residue",
            "centers");
/*	    rdft[0][6]=='m' ? "COM" : "COG");*/
  } else {
    sprintf(gtitle,d.bMassDen ? "Mass density" : "Radial distribution");
  }
  fp=xvgropen(fnRDF,gtitle,"r","");
  xvgr_selections(fp, trj);
  if (ngrps == 1) {
    if (bPrintXvgrCodes())
      fprintf(fp,"@ subtitle \"%s%s - %s\"\n",
	      d.refsel->name,bCM ? " COM" : "",sel[0]->name);
  }
  else {
    if (bPrintXvgrCodes())
      fprintf(fp,"@ subtitle \"reference %s%s\"\n",
	      d.refsel->name,bCM ? " COM" : "");
    xvgr_legend(fp, ngrps, grpnames);
  }
  gmx_histogram_write_array(fp, ngrps, rdf, TRUE, FALSE);
  ffclose(fp);
  
  do_view(fnRDF,NULL);

  /* h(Q) function: fourier transform of rdf */  
  if (fnHQ) {
    int nhq = 401;
    real *hq,*integrand,Q;
    double *hist = gmx_histogram_get_values(rdf[0]);
    
    /* Get a better number density later! */
    rho = sel[0]->p.nr*invvol;
    snew(hq,nhq);
    snew(integrand,nbin);
    for(i=0; (i<nhq); i++) {
      Q = i*0.5;
      integrand[0] = 0;
      for(j=1; (j<nbin); j++) {
	r = j*binwidth;
	integrand[j]  = (Q == 0) ? 1.0 : sin(Q*r)/(Q*r);
	integrand[j] *= 4.0*M_PI*rho*r*r*(hist[j]-1.0);
      }
      hq[i] = print_and_integrate(debug,nbin,binwidth,integrand,NULL,0);
    }
    fp=xvgropen(fnHQ,"h(Q)","Q(/nm)","h(Q)");
    xvgr_selections(fp, trj);
    for(i=0; (i<nhq); i++) 
      fprintf(fp,"%10g %10g\n",i*0.5,hq[i]);
    ffclose(fp);
    do_view(fnHQ,NULL);
    sfree(hq);
    sfree(integrand);
  }
  
  if (fnCNRDF) {
    snew(cum_rdf, ngrps);
    for (g = 0; g < ngrps; ++g) {
      gmx_histogram_resample_dblbw(&cum_rdf[g], d.count[g], FALSE);
      gmx_histogram_scale(cum_rdf[g], 1.0/d.isize0);
    }
    fp=xvgropen(fnCNRDF,"Cumulative Number RDF","r","number");
    xvgr_selections(fp, trj);
    if (ngrps == 1) {
      if (bPrintXvgrCodes())
	fprintf(fp,"@ subtitle \"%s-%s\"\n",
                d.refsel->name,sel[0]->name);
    }
    else {
      if (bPrintXvgrCodes())
	fprintf(fp,"@ subtitle \"reference %s\"\n",d.refsel->name);
      xvgr_legend(fp, ngrps, grpnames);
    }
    gmx_histogram_write_cum_array(fp, ngrps, cum_rdf, TRUE, FALSE);
    ffclose(fp);
    for (g = 0; g < ngrps; ++g)
      gmx_histogram_free(cum_rdf[g]);
    sfree(cum_rdf);
    
    do_view(fnCNRDF,NULL);
  }

  for (g = 0; g < ngrps; ++g)
    gmx_histogram_free(rdf[g]);
  sfree(rdf);
}


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
  rvec x;
  int  t;
} reduced_atom;

typedef struct
{
  real    start_q, end_q;
  int     n_angles;
  int     n_groups;
  matrix  box;
  double  lambda;
  double  energy;
  double  momentum;
  double  ref_k;
  double  **F;
  int     nSteps;
  int     total_n_atoms;
  reduced_atom **red;
  real         **table;
} structure_factor;

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
    
int return_atom_type (const char *name)
{
  typedef struct {
    const char *name;
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

static int
update_structure_factor(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                        int nr, gmx_ana_selection_t *sel[], void *data)
{
    structure_factor *sf = (structure_factor *)data;
    t_complex ***tmpSF;
    rvec k_factor;
    real kdotx, asf, kx, ky, kz, krr;
    int kr, maxkx, maxky, maxkz, i, j, k, p, *counter;
    int g;

    k_factor[XX] = 2 * M_PI / sf->box[XX][XX];
    k_factor[YY] = 2 * M_PI / sf->box[YY][YY];
    k_factor[ZZ] = 2 * M_PI / sf->box[ZZ][ZZ];

    maxkx = (int) (sf->end_q / k_factor[XX] + 0.5);
    maxky = (int) (sf->end_q / k_factor[YY] + 0.5);
    maxkz = (int) (sf->end_q / k_factor[ZZ] + 0.5);

    snew (counter, sf->n_angles);

    tmpSF = rc_tensor_allocation(maxkx,maxky,maxkz);
    for (g = 0; g < nr; g++)
    {
        rearrange_atoms(sf->red[g], fr, sel[g]->g->index, sel[g]->g->isize, top, FALSE);
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
		    if (krr >= sf->start_q && krr <= sf->end_q) {
			kr = (int) (krr/sf->ref_k + 0.5);
			if (kr < sf->n_angles) {
			    counter[kr]++;  /* will be used for the copmutation 
					       of the average*/
			    for (p = 0; p < sel[g]->g->isize; p++) {
				    
				asf = sf->table[sf->red[g][p].t][kr];

				kdotx = kx * sf->red[g][p].x[XX] +
                                    ky * sf->red[g][p].x[YY] +
                                    kz * sf->red[g][p].x[ZZ];
				
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
	kx = i * k_factor[XX];
        for (j = 0; j < maxky; j++) {
	    ky = j * k_factor[YY];
            for (k = 0; k < maxkz; k++) {
		kz = k * k_factor[ZZ];
                krr = sqrt (sqr (kx) + sqr (ky) + sqr (kz));
                if (krr >= sf->start_q && krr <= sf->end_q) {
		    kr = (int) (krr / sf->ref_k + 0.5);
                    if (kr < sf->n_angles && counter[kr] != 0)
			sf->F[g][kr] +=
			    (sqr (tmpSF[i][j][k].re) +
			     sqr (tmpSF[i][j][k].im))/ counter[kr];
		}
	    }
	}
    }
    }   /* end loop on g */
    sfree(counter);
    free(tmpSF[0][0]);
    free(tmpSF[0]);
    free(tmpSF);

    return 0;
}

void save_data (structure_factor * sf, char *file, int ngrps)
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
	if (i * sf->ref_k >= sf->start_q && i * sf->ref_k <= sf->end_q) {
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
do_scattering_intensity(gmx_ana_traj_t *trj, char* fnXVG,
		        real start_q, real end_q, real energy)
{
    int                   i;
    t_topology           *top;
    int                   ngrps;
    gmx_ana_selection_t **sel;
    t_trxframe           *fr;
    structure_factor *sf;
    int nsftable;
    double r_tmp;

    gmx_ana_get_ngrps(trj, &ngrps);
    gmx_ana_get_anagrps(trj, &sel);
    for (i = 0; i < ngrps; ++i) {
        if (sel[i]->bDynamic)
            gmx_fatal(FARGS,"Structure factor calculation with dynamic selection not supported");
        if (sel[i]->p.m.type != INDEX_ATOM)
            gmx_fatal(FARGS,"Structure factor calculation with non-atom positions not supported");
    }

    snew (sf, 1);
    sf->start_q = start_q;
    sf->end_q   = end_q;
    sf->energy  = energy;
    gmx_ana_get_topology(trj, TRUE, &top, NULL);
    gmx_ana_get_topconf(trj, NULL, sf->box, NULL);

    /* The first time we read data is a little special */
    gmx_ana_get_first_frame(trj, &fr);

    sf->total_n_atoms = fr->natoms;
    
    snew (sf->red, ngrps);

    r_tmp = max (sf->box[XX][XX], sf->box[YY][YY]);
    r_tmp = (double) max (sf->box[ZZ][ZZ], r_tmp);

    sf->ref_k = (2.0 * M_PI) / (r_tmp);
    /* ref_k will be the reference momentum unit */
    sf->n_angles = (int) (sf->end_q / sf->ref_k + 0.5);

    snew (sf->F, ngrps);
    for (i = 0; i < ngrps; i++)
       snew (sf->F[i], sf->n_angles);
    for (i = 0; i < ngrps; i++) {
	snew (sf->red[i], sel[i]->g->isize);
	rearrange_atoms (sf->red[i], fr, sel[i]->g->index, sel[i]->g->isize, top, TRUE);
    }
    sf->table = compute_scattering_factor_table (sf, &nsftable);

    /* This is the main loop over frames */
    gmx_ana_do(trj, 0, &update_structure_factor, sf);
    gmx_ana_get_nframes(trj, &sf->nSteps);

    save_data (sf, fnXVG, ngrps);

    return 0;
}

int gmx_rdf(int argc,char *argv[])
{
  const char *desc[] = {
    "The structure of liquids can be studied by either neutron or X-ray",
    "scattering. The most common way to describe liquid structure is by a",
    "radial distribution function. However, this is not easy to obtain from",
    "a scattering experiment.[PAR]",
    "g_rdf calculates radial distribution functions in different ways.",
    "The normal method is around a (set of) positions(s),",
    "an alternative method is with respect to the closest particle",
    "in a set ([TT]-surf[tt]).",
    "In all cases, rdf's can also be calculated around axes parallel",
    "to the z-axis with option [TT]-xy[tt].[PAR]",
    "Different kinds of normalization are possible with [TT]-norm[tt].",
    "Normalization cannot be used with the option [TT]-surf[tt].[PAR]",
    "Rdfs can be calculated for centers of mass/geometry of",
    "residues/molecules instead of atom positions by providing",
    "selections that evaluate to such positions.",
    "Other weighting than COM or COG can currently only be achieved",
    "by providing a run input file with different masses.",
    "If a run input file is supplied ([TT]-s[tt]) and static ",
    "atom selections are provided, exclusions defined",
    "in that file are taken into account when calculating the rdf.",
    "The option [TT]-cut[tt] is meant as an alternative way to avoid",
    "intramolecular peaks in the rdf plot.",
    "It is however better to supply a run input file with a higher number of",
    "exclusions. For eg. benzene a topology with nrexcl set to 5",
    "would eliminate all intramolecular contributions to the rdf.",
    "Note that all atoms in the selected groups are used, also the ones",
    "that don't have Lennard-Jones interactions.[PAR]",
    "Option [TT]-cn[tt] produces the cumulative number rdf,",
    "i.e. the average number of particles within a distance r.[PAR]",
    "To bridge the gap between theory and experiment structure factors can",
    "be computed (option [TT]-sq[tt]).",
  };
  static bool bXY=FALSE;
  static real cutoff=0,binwidth=0.002,grid=0.05,fade=0.0,lambda=0.1,distance=10;
  static int  npixel=256,nlevel=20;
  static real start_q=0.0, end_q=60.0, energy=12.0;

  static const char *closet[] = { NULL, "no", "res", "mol", "all", NULL };
  static const char *normt[]={ NULL, "rdf", "den_num", "den_mass", "none", NULL };

  t_pargs pa[] = {
    { "-bin",      FALSE, etREAL, {&binwidth},
      "Binwidth (nm)" },
    { "-surf",     FALSE, etENUM, {closet},
      "RDF with respect to the surface of the first group" },
    { "-norm",     FALSE, etENUM, {normt},
      "Normalization type" },
    { "-xy",       FALSE, etBOOL, {&bXY},
      "Use only the x and y components of the distance" },
    { "-cut",      FALSE, etREAL, {&cutoff},
      "Shortest distance (nm) to be considered"},
    { "-fade",     FALSE, etREAL, {&fade},
      "From this distance onwards the RDF is tranformed by g'(r) = 1 + [g(r)-1] exp(-(r/fade-1)^2 to make it go to 1 smoothly. If fade is 0.0 nothing is done." },
    { "-grid",     FALSE, etREAL, {&grid},
      "[HIDDEN]Grid spacing (in nm) for FFTs when computing structure factors" },
    { "-npixel",   FALSE, etINT,  {&npixel},
      "[HIDDEN]# pixels per edge of the square detector plate" },
    { "-nlevel",   FALSE, etINT,  {&nlevel},
      "[HIDDEN]Number of different colors in the diffraction image" },
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
  bool       bSQ,bRDF;
  gmx_ana_traj_t *trj;
  
  t_filenm   fnm[] = {
    { efXVG, "-o",  "rdf",    ffOPTWR },
    { efXVG, "-sq", "sq",     ffOPTWR },
    { efXVG, "-cn", "rdf_cn", ffOPTWR },
    { efXVG, "-hq", "hq",     ffOPTWR },
/*    { efXPM, "-image", "sq",  ffOPTWR }*/
  };
#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  gmx_ana_traj_create(&trj, 0);
  gmx_ana_set_nrefgrps(trj, 1);
  gmx_ana_set_nanagrps(trj, -1);
  parse_trjana_args(trj, &argc,argv,PCA_CAN_VIEW,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);

  bSQ   = opt2bSet("-sq",NFILE,fnm);
  bRDF  = opt2bSet("-o",NFILE,fnm) || !bSQ;

  if (closet[0][0] != 'n' && normt[0][0] != 'n')
  {
    fprintf(stderr, "Turning off normalization because of option -surf\n");
    normt[0] = normt[4];
  }
 
  if  (bSQ) 
   do_scattering_intensity(trj, opt2fn("-sq",NFILE,fnm), start_q, end_q, energy);

  if (bRDF) 
    do_rdf(trj,opt2fn("-o",NFILE,fnm),opt2fn_null("-cn",NFILE,fnm),
	   opt2fn_null("-hq",NFILE,fnm),
	   closet,normt,bXY,cutoff,binwidth,fade);

  thanx(stderr);
  
  return 0;
}
