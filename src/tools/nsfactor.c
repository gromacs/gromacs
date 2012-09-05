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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include <string.h>
#include "futil.h"
#include "gmx_random.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "strdb.h"
#include "vec.h"
#include "nsfactor.h"
#include "gmx_omp.h"

void check_binwidth(real binwidth) {
    real smallest_bin=0.1;
    if (binwidth<smallest_bin)
        gmx_fatal(FARGS,"Binwidth shouldnt be smaller then smallest bond length (H-H bond ~0.1nm) in a box");
}

void check_mcover(real mcover) {
    if (mcover>1.0) {
        gmx_fatal(FARGS,"mcover should be -1 or (0,1]");
    } else if ((mcover<0)&(mcover != -1)) {
        gmx_fatal(FARGS,"mcover should be -1 or (0,1]");
    } else {
        return;
    }
}

void normalize_probability(int n,double *a){
    int i;
    double norm=0.0;
    for (i=0;i<n;i++) norm +=a[i];
    for (i=0;i<n;i++) a[i]/=norm;
}

gmx_nentron_atomic_structurefactors_t *gmx_neutronstructurefactors_init(const char *datfn) {
    /* read nsfactor.dat */
    FILE    *fp;
    char    line[STRLEN];
    int     nralloc=10;
    int     n,p;
    int     i, line_no;
    char    atomnm[8];
    double  slength;
    gmx_nentron_atomic_structurefactors_t   *gnsf;

    fp=libopen(datfn);
    line_no = 0;
    /* allocate memory for structure */
    snew(gnsf,nralloc);
    snew(gnsf->atomnm,nralloc);
    snew(gnsf->p,nralloc);
    snew(gnsf->n,nralloc);
    snew(gnsf->slength,nralloc);

    gnsf->nratoms=line_no;

    while(get_a_line(fp,line,STRLEN)) {
        i=line_no;
        if (sscanf(line,"%s %d %d %lf",atomnm,&p,&n,&slength) == 4) {
            gnsf->atomnm[i]=strdup(atomnm);
            gnsf->n[i]=n;
            gnsf->p[i]=p;
            gnsf->slength[i]=slength;
            line_no++;
            gnsf->nratoms=line_no;
            if (line_no==nralloc){
                nralloc++;
                srenew(gnsf->atomnm,nralloc);
                srenew(gnsf->p,nralloc);
                srenew(gnsf->n,nralloc);
                srenew(gnsf->slength,nralloc);
            }
        } else
            fprintf(stderr,"WARNING: Error in file %s at line %d ignored\n",
                    datfn,line_no);
    }
    srenew(gnsf->atomnm,gnsf->nratoms);
    srenew(gnsf->p,gnsf->nratoms);
    srenew(gnsf->n,gnsf->nratoms);
    srenew(gnsf->slength,gnsf->nratoms);

    fclose(fp);

    return (gmx_nentron_atomic_structurefactors_t *) gnsf;
}

gmx_sans_t *gmx_sans_init (t_topology *top, gmx_nentron_atomic_structurefactors_t *gnsf) {
    gmx_sans_t    *gsans=NULL;
    int     i,j;
    /* Try to assing scattering length from nsfactor.dat */
    snew(gsans,1);
    snew(gsans->slength,top->atoms.nr);
    /* copy topology data */
    gsans->top = top;
    for(i=0;i<top->atoms.nr;i++) {
        for(j=0;j<gnsf->nratoms;j++) {
            if(top->atoms.atom[i].atomnumber == gnsf->p[j]) {
                /* we need special case for H and D */
                if(top->atoms.atom[i].atomnumber == 1) {
                    if(top->atoms.atom[i].m == 1.008000) {
                        gsans->slength[i] = gnsf->slength[0];
                    } else
                        gsans->slength[i] = gnsf->slength[1];
                } else
                    gsans->slength[i] = gnsf->slength[j];
            }
        }
    }

    return (gmx_sans_t *) gsans;
}

gmx_radial_distribution_histogram_t *calc_radial_distribution_histogram (
                            gmx_sans_t *gsans,
                            rvec *x,
                            matrix box,
                            atom_id *index,
                            int isize,
                            double binwidth,
                            gmx_bool bMC,
                            real mcover,
                            unsigned int seed) {
    gmx_radial_distribution_histogram_t    *pr=NULL;
    rvec            dist;
    double          rmax;
    int             i,j;
#ifdef GMX_OPENMP
    double          **tgr;
    int             tid;
    int             nthreads;
    gmx_rng_t       *trng=NULL;
#endif
    gmx_large_int_t mc=0,max;
    gmx_rng_t       rng=NULL;

    /* allocate memory for pr */
    snew(pr,1);
    /* set some fields */
    pr->binwidth=binwidth;

    /*
    * create max dist rvec
    * dist = box[xx] + box[yy] + box[zz]
    */
    rvec_add(box[XX],box[YY],dist);
    rvec_add(box[ZZ],dist,dist);

    rmax=norm(dist);

    pr->grn=(int)floor(rmax/pr->binwidth)+1;
    rmax=pr->grn*pr->binwidth;

    snew(pr->gr,pr->grn);

    if(bMC) {
        /* Special case for setting automaticaly number of mc iterations to 1% of total number of direct iterations */
        if (mcover==-1) {
            max=(gmx_large_int_t)floor(0.5*0.01*isize*(isize-1));
        } else {
            max=(gmx_large_int_t)floor(0.5*mcover*isize*(isize-1));
        }
        rng=gmx_rng_init(seed);
#ifdef GMX_OPENMP
        nthreads = gmx_omp_get_max_threads();
        snew(tgr,nthreads);
        snew(trng,nthreads);
        for(i=0;i<nthreads;i++){
            snew(tgr[i],pr->grn);
            trng[i]=gmx_rng_init(gmx_rng_uniform_uint32(rng));
        }
        #pragma omp parallel shared(tgr,trng,mc) private(tid,i,j)
        {
            tid = gmx_omp_get_thread_num();
            /* now starting parallel threads */
            #pragma omp for
            for(mc=0;mc<max;mc++) {
                i=(int)floor(gmx_rng_uniform_real(trng[tid])*isize);
                j=(int)floor(gmx_rng_uniform_real(trng[tid])*isize);
                if(i!=j) {
                    tgr[tid][(int)floor(sqrt(distance2(x[index[i]],x[index[j]]))/binwidth)]+=gsans->slength[index[i]]*gsans->slength[index[j]];
                }
            }
        }
        /* collecting data from threads */
        for(i=0;i<pr->grn;i++) {
            for(j=0;j<nthreads;j++) {
                pr->gr[i] += tgr[j][i];
            }
        }
        /* freeing memory for tgr and destroying trng */
        for(i=0;i<nthreads;i++) {
            sfree(tgr[i]);
            gmx_rng_destroy(trng[i]);
        }
        sfree(tgr);
        sfree(trng);
#else
        for(mc=0;mc<max;mc++) {
            i=(int)floor(gmx_rng_uniform_real(rng)*isize);
            j=(int)floor(gmx_rng_uniform_real(rng)*isize);
            if(i!=j)
                pr->gr[(int)floor(sqrt(distance2(x[index[i]],x[index[j]]))/binwidth)]+=gsans->slength[index[i]]*gsans->slength[index[j]];
        }
#endif
        gmx_rng_destroy(rng);
    } else {
#ifdef GMX_OPENMP
        nthreads = gmx_omp_get_max_threads();
        /* Allocating memory for tgr arrays */
        snew(tgr,nthreads);
        for(i=0;i<nthreads;i++) {
            snew(tgr[i],pr->grn);
        }
        #pragma omp parallel shared(tgr) private(tid,i,j)
        {
            tid = gmx_omp_get_thread_num();
            /* starting parallel threads */
            #pragma omp for
            for(i=0;i<isize;i++) {
                for(j=0;j<i;j++) {
                    tgr[tid][(int)floor(sqrt(distance2(x[index[i]],x[index[j]]))/binwidth)]+=gsans->slength[index[i]]*gsans->slength[index[j]];
                }
            }
        }
        /* collecating data for pr->gr */
        for(i=0;i<pr->grn;i++) {
            for(j=0;j<nthreads;j++) {
                pr->gr[i] += tgr[j][i];
            }
        }
        /* freeing memory for tgr */
        for(i=0;i<nthreads;i++) {
            sfree(tgr[i]);
        }
        sfree(tgr);
#else
        for(i=0;i<isize;i++) {
            for(j=0;j<i;j++) {
                pr->gr[(int)floor(sqrt(distance2(x[index[i]],x[index[j]]))/binwidth)]+=gsans->slength[index[i]]*gsans->slength[index[j]];
            }
        }
#endif
    }

    /* normalize */
    normalize_probability(pr->grn,pr->gr);
    snew(pr->r,pr->grn);
    for(i=0;i<pr->grn;i++)
        pr->r[i]=(pr->binwidth*i+pr->binwidth*0.5);

    return (gmx_radial_distribution_histogram_t *) pr;
}

gmx_static_structurefator_t *convert_histogram_to_intensity_curve (gmx_radial_distribution_histogram_t *pr, double start_q, double end_q, double q_step) {
    gmx_static_structurefator_t    *sq=NULL;
    int         i,j;
    /* init data */
    snew(sq,1);
    sq->qn=(int)floor((end_q-start_q)/q_step);
    snew(sq->q,sq->qn);
    snew(sq->s,sq->qn);
    for(i=0;i<sq->qn;i++)
        sq->q[i]=start_q+i*q_step;

    if(start_q==0.0) {
        sq->s[0]=1.0;
        for(i=1;i<sq->qn;i++) {
            for(j=0;j<pr->grn;j++)
                sq->s[i]+=(pr->gr[j]/pr->r[j])*sin(sq->q[i]*pr->r[j]);
            sq->s[i] /= sq->q[i];
        }
    } else {
        for(i=0;i<sq->qn;i++) {
            for(j=0;j<pr->grn;j++)
                sq->s[i]+=(pr->gr[j]/pr->r[j])*sin(sq->q[i]*pr->r[j]);
            sq->s[i] /= sq->q[i];
        }
    }

    return (gmx_static_structurefator_t *) sq;
}
