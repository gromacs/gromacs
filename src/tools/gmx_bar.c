/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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
#include <ctype.h>
#include <math.h>
#include <float.h>

#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "futil.h"
#include "statutil.h"
#include "copyrite.h"
#include "macros.h"
#include "enxio.h"
#include "physics.h"
#include "gmx_fatal.h"
#include "xvgr.h"
#include "gmx_ana.h"
#include "maths.h"


/* the dhdl.xvg data from a simulation (actually obsolete, but still
    here for reading the dhdl.xvg file*/
typedef struct xvg_t
{
    char   *filename;
    int    ftp;     /* file type */
    int    nset;    /* number of lambdas, including dhdl */
    int *np;        /* number of data points (du or hists) per lambda */
    int  np_alloc;  /* number of points (du or hists) allocated */
    double temp;    /* temperature */
    double *lambda; /* the lambdas (of first index for y). */
    double *t;      /* the times (of second index for y) */
    double **y;     /* the dU values. y[0] holds the derivative, while
                       further ones contain the energy differences between
                       the native lambda and the 'foreign' lambdas. */

    double native_lambda; /* the native lambda */

    struct xvg_t *next, *prev; /*location in the global linked list of xvg_ts*/
} xvg_t;



typedef struct hist_t
{
    unsigned int *bin[2];       /* the (forward + reverse) histogram values */
    double dx[2];               /* the histogram spacing. The reverse
                                   dx is the negative of the forward dx.*/
    gmx_large_int_t x0[2];      /* the (forward + reverse) histogram start 
                                   point(s) as int */

    int nbin[2];                /* the (forward+reverse) number of bins */
    gmx_large_int_t sum;        /* the total number of counts. Must be
                                   the same for forward + reverse.  */
    int nhist;                  /* number of hist datas (forward or reverse) */

    double start_time, delta_time;  /* start time, end time of histogram */
} hist_t;


/* an aggregate of samples for partial free energy calculation */
typedef struct samples_t 
{    
    double native_lambda; 
    double foreign_lambda;
    double temp; /* the temperature */
    gmx_bool derivative; /* whether this sample is a derivative */

    /* The samples come either as either delta U lists: */
    int ndu; /* the number of delta U samples */
    double *du; /* the delta u's */
    double *t; /* the times associated with those samples, or: */
    double start_time,delta_time;/*start time and delta time for linear time*/

    /* or as histograms: */
    hist_t *hist; /* a histogram */

    /* allocation data: (not NULL for data 'owned' by this struct) */
    double *du_alloc, *t_alloc; /* allocated delta u arrays  */
    size_t ndu_alloc, nt_alloc; /* pre-allocated sizes */
    hist_t *hist_alloc; /* allocated hist */

    gmx_large_int_t ntot; /* total number of samples */
    const char *filename; /* the file name this sample comes from */
} samples_t;

/* a sample range (start to end for du-style data, or boolean 
    for both du-style data and histograms */
typedef struct sample_range_t
{
    int start, end; /* start and end index for du style data */
    gmx_bool use; /* whether to use this sample */

    samples_t *s; /* the samples this range belongs to */
} sample_range_t;


/* a collection of samples for a partial free energy calculation 
    (i.e. the collection of samples from one native lambda to one 
    foreign lambda) */
typedef struct sample_coll_t
{
    double native_lambda;  /* these should be the same for all samples in the */
    double foreign_lambda; /* collection */
    double temp; /* the temperature */

    int nsamples; /* the number of samples */
    samples_t **s; /* the samples themselves */
    sample_range_t *r; /* the sample ranges */
    int nsamples_alloc; /* number of allocated samples */ 

    gmx_large_int_t ntot; /* total number of samples in the ranges of 
                             this collection */

    struct sample_coll_t *next, *prev; /* next and previous in the list */
} sample_coll_t;

/* all the samples associated with a lambda point */
typedef struct lambda_t
{
    double lambda; /* the native lambda (at start time if dynamic) */
    double temp; /* temperature */

    sample_coll_t *sc; /* the samples */

    sample_coll_t sc_head; /*the pre-allocated list head for the linked list.*/

    struct lambda_t *next, *prev; /* the next and prev in the list */
} lambda_t;


/* calculated values. */
typedef struct {
    sample_coll_t *a, *b; /* the simulation data */

    double dg; /* the free energy difference */
    double dg_err; /* the free energy difference */

    double dg_disc_err; /* discretization error */
    double dg_histrange_err; /* histogram range error */

    double sa; /* relative entropy of b in state a */
    double sa_err; /* error in sa */
    double sb; /* relative entropy of a in state b */
    double sb_err; /* error in sb */

    double dg_stddev; /* expected dg stddev per sample */
    double dg_stddev_err; /* error in dg_stddev */
} barres_t;




static void hist_init(hist_t *h, int nhist, int *nbin)
{
    int i;
    if (nhist>2)
    {
        gmx_fatal(FARGS, "histogram with more than two sets of data!");
    }
    for(i=0;i<nhist;i++)
    {
        snew(h->bin[i], nbin[i]);
        h->x0[i]=0;
        h->nbin[i]=nbin[i];
        h->start_time=h->delta_time=0;
        h->dx[i]=0;
    }
    h->sum=0;
    h->nhist=nhist;
}

static void hist_destroy(hist_t *h)
{
    sfree(h->bin);
}


static void xvg_init(xvg_t *ba)
{
    ba->filename=NULL;
    ba->nset=0;
    ba->np_alloc=0;
    ba->np=NULL;
    ba->y=NULL;
}

static void samples_init(samples_t *s, double native_lambda,
                         double foreign_lambda, double temp,
                         gmx_bool derivative, const char *filename)
{
    s->native_lambda=native_lambda;
    s->foreign_lambda=foreign_lambda;
    s->temp=temp;
    s->derivative=derivative;

    s->ndu=0;
    s->du=NULL;
    s->t=NULL;
    s->start_time = s->delta_time = 0;
    s->hist=NULL;
    s->du_alloc=NULL;
    s->t_alloc=NULL;
    s->hist_alloc=NULL;
    s->ndu_alloc=0;
    s->nt_alloc=0;

    s->ntot=0;
    s->filename=filename;
}

/* destroy the data structures directly associated with the structure, not
   the data it points to */
static void samples_destroy(samples_t *s)
{
    if (s->du_alloc)
    {
        sfree(s->du_alloc);
    }
    if (s->t_alloc)
    {
        sfree(s->t_alloc);
    }
    if (s->hist_alloc)
    {
        hist_destroy(s->hist_alloc);
        sfree(s->hist_alloc);
    }
}

static void sample_range_init(sample_range_t *r, samples_t *s)
{
    r->start=0;
    r->end=s->ndu;
    r->use=TRUE;
    r->s=NULL;
}

static void sample_coll_init(sample_coll_t *sc, double native_lambda,
                             double foreign_lambda, double temp)
{
    sc->native_lambda = native_lambda;
    sc->foreign_lambda = foreign_lambda;
    sc->temp = temp;

    sc->nsamples=0;
    sc->s=NULL;
    sc->r=NULL;
    sc->nsamples_alloc=0;

    sc->ntot=0;
    sc->next=sc->prev=NULL;
}

static void sample_coll_destroy(sample_coll_t *sc)
{
    /* don't free the samples themselves */
    sfree(sc->r);
    sfree(sc->s);
}


static void lambda_init(lambda_t *l, double native_lambda, double temp)
{
    l->lambda=native_lambda;
    l->temp=temp;

    l->next=NULL;
    l->prev=NULL;

    l->sc=&(l->sc_head);

    sample_coll_init(l->sc, native_lambda, 0., 0.);
    l->sc->next=l->sc;
    l->sc->prev=l->sc;
}

static void barres_init(barres_t *br)
{
    br->dg=0;
    br->dg_err=0;
    br->sa=0;
    br->sa_err=0;
    br->sb=0;
    br->sb_err=0;
    br->dg_stddev=0;
    br->dg_stddev_err=0;

    br->a=NULL;
    br->b=NULL;
}




static gmx_bool lambda_same(double lambda1, double lambda2)
{
    return gmx_within_tol(lambda1, lambda2, 10*GMX_REAL_EPS);
}

/* calculate the total number of samples in a sample collection */
static void sample_coll_calc_ntot(sample_coll_t *sc)
{
    int i;

    sc->ntot=0;
    for(i=0;i<sc->nsamples;i++)
    {
        if (sc->r[i].use)
        {
            if (sc->s[i]->hist)
            {
                sc->ntot += sc->s[i]->ntot;
            }
            else
            {
                sc->ntot += sc->r[i].end - sc->r[i].start;
            }
        }
    }
}


/* find the barsamples_t associated with a lambda that corresponds to
   a specific foreign lambda */
static sample_coll_t *lambda_find_sample_coll(lambda_t *l, 
                                              double foreign_lambda)
{
    sample_coll_t *sc=l->sc->next;

    while(sc != l->sc)
    {
        if (lambda_same(sc->foreign_lambda,foreign_lambda))
        {
            return sc;
        }
        sc=sc->next;
    }

    return NULL;
}

/* insert li into an ordered list of lambda_colls */
static void lambda_insert_sample_coll(lambda_t *l, sample_coll_t *sc)
{
    sample_coll_t *scn=l->sc->next;
    while ( (scn!=l->sc) )
    {
        if (scn->foreign_lambda > sc->foreign_lambda)
            break;
        scn=scn->next;
    }
    /* now insert it before the found scn */
    sc->next=scn;
    sc->prev=scn->prev;
    scn->prev->next=sc;
    scn->prev=sc;
}

/* insert li into an ordered list of lambdas */
static void lambda_insert_lambda(lambda_t *head, lambda_t *li)
{
    lambda_t *lc=head->next;
    while (lc!=head) 
    {
        if (lc->lambda > li->lambda)
            break;
        lc=lc->next;
    }
    /* now insert ourselves before the found lc */
    li->next=lc;
    li->prev=lc->prev;
    lc->prev->next=li;
    lc->prev=li;
}

/* insert a sample and a sample_range into a sample_coll. The
    samples are stored as a pointer, the range is copied. */
static void sample_coll_insert_sample(sample_coll_t *sc, samples_t *s,
                                      sample_range_t *r)
{
    /* first check if it belongs here */
    if (sc->temp != s->temp)
    {
        gmx_fatal(FARGS, "Temperatures in files %s and %s are not the same!",
                   s->filename, sc->next->s[0]->filename);
    }
    if (sc->native_lambda != s->native_lambda)
    {
        gmx_fatal(FARGS, "Native lambda in files %s and %s are not the same (and they should be)!",
                   s->filename, sc->next->s[0]->filename);
    }
    if (sc->foreign_lambda != s->foreign_lambda)
    {
        gmx_fatal(FARGS, "Foreign lambda in files %s and %s are not the same (and they should be)!",
                   s->filename, sc->next->s[0]->filename);
    }

    /* check if there's room */
    if ( (sc->nsamples + 1) > sc->nsamples_alloc)
    {
        sc->nsamples_alloc = max(2*sc->nsamples_alloc, 2);
        srenew(sc->s, sc->nsamples_alloc);
        srenew(sc->r, sc->nsamples_alloc);
    }
    sc->s[sc->nsamples]=s;
    sc->r[sc->nsamples]=*r;
    sc->nsamples++;

    sample_coll_calc_ntot(sc);
}

/* insert a sample into a lambda_list, creating the right sample_coll if 
   neccesary */
static void lambda_list_insert_sample(lambda_t *head, samples_t *s)
{
    gmx_bool found=FALSE;
    sample_coll_t *sc;
    sample_range_t r;

    lambda_t *l=head->next;

    /* first search for the right lambda_t */
    while(l != head)
    {
        if (lambda_same(l->lambda, s->native_lambda) )
        {
            found=TRUE;
            break;
        }
        l=l->next;
    }

    if (!found)
    {
        snew(l, 1); /* allocate a new one */
        lambda_init(l, s->native_lambda, s->temp); /* initialize it */
        lambda_insert_lambda(head, l); /* add it to the list */
    }

    /* now look for a sample collection */
    sc=lambda_find_sample_coll(l, s->foreign_lambda);
    if (!sc)
    {
        snew(sc, 1); /* allocate a new one */
        sample_coll_init(sc, s->native_lambda, s->foreign_lambda, s->temp);
        lambda_insert_sample_coll(l, sc);
    }

    /* now insert the samples into the sample coll */
    sample_range_init(&r, s);
    sample_coll_insert_sample(sc, s, &r);
}


/* make a histogram out of a sample collection */
static void sample_coll_make_hist(sample_coll_t *sc, int **bin, 
                                  int *nbin_alloc, int *nbin, 
                                  double *dx, double *xmin, int nbin_default)
{
    int i,j,k;
    gmx_bool dx_set=FALSE;
    gmx_bool xmin_set=FALSE;

    gmx_bool xmax_set=FALSE;
    gmx_bool xmax_set_hard=FALSE; /* whether the xmax is bounded by the 
                                     limits of a histogram */
    double xmax=-1;

    /* first determine dx and xmin; try the histograms */
    for(i=0;i<sc->nsamples;i++)
    {
        if (sc->s[i]->hist)
        {
            hist_t *hist=sc->s[i]->hist;
            for(k=0;k<hist->nhist;k++)
            {
                double hdx=hist->dx[k];
                double xmax_now=(hist->x0[k]+hist->nbin[k])*hdx;

                /* we use the biggest dx*/
                if ( (!dx_set) || hist->dx[0] > *dx)
                {
                    dx_set=TRUE;
                    *dx = hist->dx[0];
                }
                if ( (!xmin_set) || (hist->x0[k]*hdx) < *xmin)
                {
                    xmin_set=TRUE;
                    *xmin = (hist->x0[k]*hdx);
                }

                if ( (!xmax_set) || (xmax_now>xmax && !xmax_set_hard) )
                {
                    xmax_set=TRUE;
                    xmax = xmax_now;
                    if (hist->bin[k][hist->nbin[k]-1] != 0)
                        xmax_set_hard=TRUE;
                }
                if ( hist->bin[k][hist->nbin[k]-1]!=0 && (xmax_now < xmax) )
                {
                    xmax_set_hard=TRUE;
                    xmax = xmax_now;
                }
            }
        }
    }
    /* and the delta us */
    for(i=0;i<sc->nsamples;i++)
    {
        if (sc->s[i]->ndu > 0)
        {
            /* determine min and max */
            int starti=sc->r[i].start;
            int endi=sc->r[i].end;
            double du_xmin=sc->s[i]->du[starti]; 
            double du_xmax=sc->s[i]->du[starti];
            for(j=starti+1;j<endi;j++)
            {
                if (sc->s[i]->du[j] < du_xmin)
                    du_xmin = sc->s[i]->du[j];
                if (sc->s[i]->du[j] > du_xmax)
                    du_xmax = sc->s[i]->du[j];
            }

            /* and now change the limits */
            if ( (!xmin_set) || (du_xmin < *xmin) )
            {
                xmin_set=TRUE;
                *xmin=du_xmin;
            }
            if ( (!xmax_set) || ((du_xmax > xmax) &&  !xmax_set_hard) )
            {
                xmax_set=TRUE;
                xmax=du_xmax;
            }
        }
    }

    if (!xmax_set || !xmin_set)
    {
        *nbin=0;
        return;
    }


    if (!dx_set)
    {
        *nbin=nbin_default;
        *dx=(xmax-(*xmin))/((*nbin)-2); /* -2 because we want the last bin to
                                           be 0, and we count from 0 */
    }
    else
    {
        *nbin=(xmax-(*xmin))/(*dx);
    }

    if (*nbin > *nbin_alloc)
    {
        *nbin_alloc=*nbin;
        srenew(*bin, *nbin_alloc);
    }

    /* reset the histogram */
    for(i=0;i<(*nbin);i++)
    {
        (*bin)[i] = 0;
    }

    /* now add the actual data */   
    for(i=0;i<sc->nsamples;i++)
    {
        if (sc->s[i]->hist)
        {
            hist_t *hist=sc->s[i]->hist;
            for(k=0;k<hist->nhist;k++)
            {
                double hdx = hist->dx[k];
                double xmin_hist=hist->x0[k]*hdx;
                for(j=0;j<hist->nbin[k];j++)
                {
                    /* calculate the bin corresponding to the middle of the 
                       original bin */
                    double x=hdx*(j+0.5) + xmin_hist;
                    int binnr=(int)((x-(*xmin))/(*dx));

                    if (binnr >= *nbin || binnr<0)
                        binnr = (*nbin)-1;

                    (*bin)[binnr] += hist->bin[k][j]; 
                }
            }
        }
        else
        {
            int starti=sc->r[i].start;
            int endi=sc->r[i].end;
            for(j=starti;j<endi;j++)
            {
                int binnr=(int)((sc->s[i]->du[j]-(*xmin))/(*dx));
                if (binnr >= *nbin || binnr<0)
                    binnr = (*nbin)-1;

                (*bin)[binnr] ++;
            }
        }
    }
}

/* write a collection of histograms to a file */
void lambdas_histogram(lambda_t *bl_head, const char *filename, 
                       int nbin_default, const output_env_t oenv)
{
    char label_x[STRLEN];
    const char *dhdl="dH/d\\lambda",*deltag="\\DeltaH",*lambda="\\lambda";
    const char *title="N(\\DeltaH)";
    const char *label_y="Samples";
    FILE *fp;
    lambda_t *bl;
    int nsets=0;
    char **setnames=NULL;
    gmx_bool first_set=FALSE;
    /* histogram data: */
    int *hist=NULL;
    int nbin=0;
    int nbin_alloc=0;
    double dx=0;
    double min=0;
    int i;

    printf("\nWriting histogram to %s\n", filename);
    sprintf(label_x, "\\DeltaH (%s)", unit_energy);

    fp=xvgropen_type(filename, title, label_x, label_y, exvggtXNY, oenv);

    /* first get all the set names */
    bl=bl_head->next;
    /* iterate over all lambdas */
    while(bl!=bl_head)
    {
        sample_coll_t *sc=bl->sc->next;

        /* iterate over all samples */
        while(sc!=bl->sc)
        {
            nsets++;
            srenew(setnames, nsets); 
            snew(setnames[nsets-1], STRLEN);
            if (!lambda_same(sc->foreign_lambda, sc->native_lambda))
            {
                sprintf(setnames[nsets-1], "N(%s(%s=%g) | %s=%g)",
                        deltag, lambda, sc->foreign_lambda, lambda,
                        sc->native_lambda);
            }
            else
            {
                sprintf(setnames[nsets-1], "N(%s | %s=%g)",
                        dhdl, lambda, sc->native_lambda);
            }
            sc=sc->next;
        }

        bl=bl->next;
    }
    xvgr_legend(fp, nsets, (const char**)setnames, oenv);


    /* now make the histograms */
    bl=bl_head->next;
    /* iterate over all lambdas */
    while(bl!=bl_head)
    {
        sample_coll_t *sc=bl->sc->next;

        /* iterate over all samples */
        while(sc!=bl->sc)
        {
            if (!first_set)
            {
                xvgr_new_dataset(fp, 0, 0, NULL, oenv);
            }
    
            sample_coll_make_hist(sc, &hist, &nbin_alloc, &nbin, &dx, &min,
                                  nbin_default);

            for(i=0;i<nbin;i++)
            {
                double xmin=i*dx + min;
                double xmax=(i+1)*dx + min;

                fprintf(fp, "%g %d\n%g %d\n", xmin, hist[i], xmax, hist[i]);
            }

            first_set=FALSE;
            sc=sc->next;
        }

        bl=bl->next;
    }

    if(hist)
        sfree(hist);    

    xvgrclose(fp); 
}

/* create a collection (array) of barres_t object given a ordered linked list 
   of barlamda_t sample collections */
static barres_t *barres_list_create(lambda_t *bl_head, int *nres)
{
    lambda_t *bl;
    int nlambda=0;
    barres_t *res;
    int i;
    gmx_bool dhdl=FALSE;
    gmx_bool first=TRUE;

    /* first count the lambdas */
    bl=bl_head->next;
    while(bl!=bl_head)
    {
        nlambda++;    
        bl=bl->next;
    }
    snew(res, nlambda-1);

    /* next put the right samples in the res */
    *nres=0;
    bl=bl_head->next->next; /* we start with the second one. */
    while(bl!=bl_head)
    {
        sample_coll_t *sc, *scprev;
        barres_t *br=&(res[*nres]);
        /* there is always a previous one. we search for that as a foreign 
           lambda: */
        scprev=lambda_find_sample_coll(bl->prev, bl->lambda);
        sc=lambda_find_sample_coll(bl, bl->prev->lambda);

        barres_init(br);

        if (!scprev && !sc)
        {
            /* we use dhdl */

            scprev=lambda_find_sample_coll(bl->prev, bl->prev->lambda);
            sc=lambda_find_sample_coll(bl, bl->lambda);

            if (first)
            {
                printf("\nWARNING: Using the derivative data (dH/dlambda) to extrapolate delta H values.\nThis will only work if the Hamiltonian is linear in lambda.\n");
                dhdl=TRUE;
            }
            if (!dhdl)
            {
                gmx_fatal(FARGS,"Some dhdl files contain only one value (dH/dl), while others \ncontain multiple values (dH/dl and/or Delta H), will not proceed \nbecause of possible inconsistencies.\n");
            }
        }
        
        /* normal delta H */
        if (!scprev)
        {
            gmx_fatal(FARGS,"Could not find a set for lambda = %g in the files for lambda = %g",bl->lambda,bl->prev->lambda);
        }
        if (!sc)
        {
            gmx_fatal(FARGS,"Could not find a set for lambda = %g in the files for lambda = %g",bl->prev->lambda,bl->lambda);
        }
        br->a = scprev;
        br->b = sc;

        first=FALSE;
        (*nres)++;
        bl=bl->next;
    }
    return res;
}

/* estimate the maximum discretization error */
static double barres_list_max_disc_err(barres_t *res, int nres)
{
    int i,j;
    double disc_err=0.;
    double delta_lambda;

    for(i=0;i<nres;i++)
    {
        barres_t *br=&(res[i]);

        delta_lambda=fabs(br->b->native_lambda-br->a->native_lambda);

        for(j=0;j<br->a->nsamples;j++)
        {
            if (br->a->s[j]->hist)
            {
                double Wfac=1.;
                if (br->a->s[j]->derivative)
                    Wfac =  delta_lambda;

                disc_err=max(disc_err, Wfac*br->a->s[j]->hist->dx[0]);
            }
        }
        for(j=0;j<br->b->nsamples;j++)
        {
            if (br->b->s[j]->hist)
            {
                double Wfac=1.;
                if (br->b->s[j]->derivative)
                    Wfac =  delta_lambda;
                disc_err=max(disc_err, Wfac*br->b->s[j]->hist->dx[0]);
            }
        }
    } 
    return disc_err;
}


/* impose start and end times on a sample collection, updating sample_ranges */
static void sample_coll_impose_times(sample_coll_t *sc, double begin_t, 
                                     double end_t)
{
    int i;
    for(i=0;i<sc->nsamples;i++)
    {
        samples_t *s=sc->s[i];
        sample_range_t *r=&(sc->r[i]);
        if (s->hist)
        {
            double end_time=s->hist->delta_time*s->hist->sum + 
                            s->hist->start_time;
            if (s->hist->start_time < begin_t || end_time > end_t)
            {
                r->use=FALSE;
            }
        }
        else
        {
            if (!s->t)
            {
                double end_time;
                if (s->start_time < begin_t)
                {
                    r->start=(int)((begin_t - s->start_time)/s->delta_time);
                }
                end_time=s->delta_time*s->ndu + s->start_time;
                if (end_time > end_t)
                {
                    r->end=(int)((end_t - s->start_time)/s->delta_time);
                }
            }
            else
            {
                int j;
                for(j=0;j<s->ndu;j++)
                {
                    if (s->t[j] <begin_t)
                    {
                        r->start = j;
                    }

                    if (s->t[j] >= end_t)
                    {
                        r->end = j;
                        break;
                    }
                }
            }
            if (r->start > r->end)
            {
                r->use=FALSE;
            }
        }
    }
    sample_coll_calc_ntot(sc);
}

static void lambdas_impose_times(lambda_t *head, double begin, double end)
{
    double first_t, last_t;
    double begin_t, end_t;
    lambda_t *lc;
    int j;

    if (begin<=0 && end<0) 
    {
        return;
    }

    /* first determine the global start and end times */
    first_t = -1;
    last_t = -1;
    lc=head->next;
    while(lc!=head)
    {
        sample_coll_t *sc=lc->sc->next;
        while(sc != lc->sc)
        {
            for(j=0;j<sc->nsamples;j++)
            {
                double start_t,end_t;

                start_t = sc->s[j]->start_time;
                end_t =   sc->s[j]->start_time;
                if (sc->s[j]->hist)
                {
                    end_t += sc->s[j]->delta_time*sc->s[j]->hist->sum;
                }
                else
                {
                    if (sc->s[j]->t)
                    {
                        end_t = sc->s[j]->t[sc->s[j]->ndu-1];
                    }
                    else
                    {
                        end_t += sc->s[j]->delta_time*sc->s[j]->ndu;
                    }
                }

                if (start_t < first_t || first_t<0)
                {
                    first_t=start_t;
                }
                if (end_t > last_t)
                {
                    last_t=end_t;
                }
            }
            sc=sc->next;
        }
        lc=lc->next;
    }

    /* calculate the actual times */
    if (begin > 0 )
    {
        begin_t = begin;
    }
    else
    {
        begin_t = first_t;
    }

    if (end >0 )
    {
        end_t = end;
    }
    else
    {
        end_t = last_t;
    }
    printf("\n   Samples in time interval: %.3f - %.3f\n", first_t, last_t);

    if (begin_t > end_t)
    {
        return;
    }
    printf("Removing samples outside of: %.3f - %.3f\n", begin_t, end_t);

    /* then impose them */
    lc=head->next;
    while(lc!=head)
    {
        sample_coll_t *sc=lc->sc->next;
        while(sc != lc->sc)
        {
            sample_coll_impose_times(sc, begin_t, end_t);
            sc=sc->next;
        }
        lc=lc->next;
    }
}


/* create subsample i out of ni from an existing sample_coll */
static gmx_bool sample_coll_create_subsample(sample_coll_t  *sc, 
                                             sample_coll_t *sc_orig, 
                                             int i, int ni)
{
    int j;
    int hist_start, hist_end;

    gmx_large_int_t ntot_start;
    gmx_large_int_t ntot_end;
    gmx_large_int_t ntot_so_far;

    *sc = *sc_orig; /* just copy all fields */

    /* allocate proprietary memory */
    snew(sc->s, sc_orig->nsamples);
    snew(sc->r, sc_orig->nsamples);

    /* copy the samples */
    for(j=0;j<sc_orig->nsamples;j++)
    {
        sc->s[j] = sc_orig->s[j];
        sc->r[j] = sc_orig->r[j]; /* copy the ranges too */
    }
    
    /* now fix start and end fields */
    /* the casts avoid possible overflows */
    ntot_start=(gmx_large_int_t)(sc_orig->ntot*(double)i/(double)ni);
    ntot_end  =(gmx_large_int_t)(sc_orig->ntot*(double)(i+1)/(double)ni);
    ntot_so_far = 0;
    for(j=0;j<sc->nsamples;j++)
    {
        gmx_large_int_t ntot_add;
        gmx_large_int_t new_start, new_end;

        if (sc->r[j].use)
        {
            if (sc->s[j]->hist)
            {
                ntot_add = sc->s[j]->hist->sum;
            }
            else 
            {
                ntot_add = sc->r[j].end - sc->r[j].start;
            }
        }
        else
        {
            ntot_add = 0;
        }

        if (!sc->s[j]->hist)
        { 
            if (ntot_so_far < ntot_start)
            {
                /* adjust starting point */
                new_start = sc->r[j].start + (ntot_start - ntot_so_far);
            }
            else
            {
                new_start = sc->r[j].start;
            }
            /* adjust end point */
            new_end = sc->r[j].start + (ntot_end - ntot_so_far);
            if (new_end > sc->r[j].end)
            {
                new_end=sc->r[j].end;
            }

            /* check if we're in range at all */
            if ( (new_end < new_start) || (new_start > sc->r[j].end) )
            {
                new_start=0;
                new_end=0;
            }
            /* and write the new range */
            sc->r[j].start=(int)new_start;
            sc->r[j].end=(int)new_end;
        }
        else
        {
            if (sc->r[j].use)
            {
                double overlap;
                double ntot_start_norm, ntot_end_norm;
                /* calculate the amount of overlap of the 
                   desired range (ntot_start -- ntot_end) onto
                   the histogram range (ntot_so_far -- ntot_so_far+ntot_add)*/

                /* first calculate normalized bounds 
                   (where 0 is the start of the hist range, and 1 the end) */
                ntot_start_norm = (ntot_start-ntot_so_far)/(double)ntot_add;
                ntot_end_norm = (ntot_end-ntot_so_far)/(double)ntot_add;

                /* now fix the boundaries */
                ntot_start_norm = min(1, max(0., ntot_start_norm));
                ntot_end_norm = max(0, min(1., ntot_end_norm));

                /* and calculate the overlap */
                overlap = ntot_end_norm - ntot_start_norm;

                if (overlap > 0.95) /* we allow for 5% slack */
                {
                    sc->r[j].use = TRUE;
                }
                else if (overlap < 0.05)
                {
                    sc->r[j].use = FALSE;
                }
                else
                {
                    return FALSE;
                }
            }
        }
        ntot_so_far += ntot_add;
    }
    sample_coll_calc_ntot(sc);

    return TRUE;
}

/* calculate minimum and maximum work values in sample collection */
static void sample_coll_min_max(sample_coll_t *sc, double Wfac,
                                double *Wmin, double *Wmax)
{
    int i,j;

    *Wmin=FLT_MAX;
    *Wmax=-FLT_MAX;

    for(i=0;i<sc->nsamples;i++)
    {
        samples_t *s=sc->s[i];
        sample_range_t *r=&(sc->r[i]);
        if (r->use)
        {
            if (!s->hist)
            {
                for(j=r->start; j<r->end; j++)
                {
                    *Wmin = min(*Wmin,s->du[j]*Wfac);
                    *Wmax = max(*Wmax,s->du[j]*Wfac);
                }
            }
            else
            {
                int hd=0; /* determine the histogram direction: */
                double dx;
                if ( (s->hist->nhist>1) && (Wfac<0) )
                {
                    hd=1;
                }
                dx=s->hist->dx[hd];

                for(j=s->hist->nbin[hd]-1;j>=0;j--)
                {
                    *Wmin=min(*Wmin,Wfac*(s->hist->x0[hd])*dx);
                    *Wmax=max(*Wmax,Wfac*(s->hist->x0[hd])*dx);
                    /* look for the highest value bin with values */
                    if (s->hist->bin[hd][j]>0)
                    {
                        *Wmin=min(*Wmin,Wfac*(j+s->hist->x0[hd]+1)*dx);
                        *Wmax=max(*Wmax,Wfac*(j+s->hist->x0[hd]+1)*dx);
                        break;
                    }
                }
            }
        }
    }
}


static double calc_bar_sum(int n,const double *W,double Wfac,double sbMmDG)
{
    int    i;
    double sum;
    
    sum = 0;
    
    for(i=0; i<n; i++)
    {
        sum += 1./(1. + exp(Wfac*W[i] + sbMmDG));
    }
    
    return sum;
}

/* calculate the BAR average given a histogram 

    if type== 0, calculate the best estimate for the average,
    if type==-1, calculate the minimum possible value given the histogram 
    if type== 1, calculate the maximum possible value given the histogram */
static double calc_bar_sum_hist(const hist_t *hist, double Wfac, double sbMmDG,
                                int type)
{
    double sum=0.;
    int i;
    int max; 
    /* normalization factor multiplied with bin width and
       number of samples (we normalize through M): */
    double normdx = 1.;
    int hd=0; /* determine the histogram direction: */
    double dx;

    if ( (hist->nhist>1) && (Wfac<0) )
    {
        hd=1;
    }
    dx=hist->dx[hd];
    max=hist->nbin[hd]-1;
    if (type==1) 
    {
        max=hist->nbin[hd]; /* we also add whatever was out of range */
    }

    for(i=0;i<max;i++)
    {
        double x=Wfac*((i+hist->x0[hd])+0.5)*dx; /* bin middle */
        double pxdx=hist->bin[0][i]*normdx; /* p(x)dx */
   
        sum += pxdx/(1. + exp(x + sbMmDG));
    }

    return sum;
}

static double calc_bar_lowlevel(sample_coll_t *ca, sample_coll_t *cb,
                                double temp, double tol, int type)
{
    double kT,beta,M;
    double DG;
    int i,j;
    double Wfac1,Wfac2,Wmin,Wmax;
    double DG0,DG1,DG2,dDG1;
    double sum1,sum2;
    double n1, n2; /* numbers of samples as doubles */
    
    kT   = BOLTZ*temp;
    beta = 1/kT;
  
    /* count the numbers of samples */ 
    n1 = ca->ntot;
    n2 = cb->ntot;

    M = log(n1/n2);

    if (!lambda_same(ca->native_lambda, ca->foreign_lambda))
    {
        /* this is the case when the delta U were calculated directly
           (i.e. we're not scaling dhdl) */
        Wfac1 = beta;
        Wfac2 = beta;
    }
    else
    {
        /* we're using dhdl, so delta_lambda needs to be a 
           multiplication factor.  */
        double delta_lambda=cb->native_lambda-ca->native_lambda;
        Wfac1 =  beta*delta_lambda;
        Wfac2 = -beta*delta_lambda;
    }

    if (beta < 1)
    {
        /* We print the output both in kT and kJ/mol.
         * Here we determine DG in kT, so when beta < 1
         * the precision has to be increased.
         */
        tol *= beta;
    }

    /* Calculate minimum and maximum work to give an initial estimate of 
     * delta G  as their average.
     */
    {
        double Wmin1, Wmin2, Wmax1, Wmax2;
        sample_coll_min_max(ca, Wfac1, &Wmin1, &Wmax1);
        sample_coll_min_max(cb, Wfac2, &Wmin2, &Wmax2);

        Wmin=min(Wmin1, Wmin2);
        Wmax=max(Wmax1, Wmax2);
    }

    DG0 = Wmin;
    DG2 = Wmax;
    
    if (debug)
    {
        fprintf(debug,"DG %9.5f %9.5f\n",DG0,DG2);
    }
    /* We approximate by bisection: given our initial estimates
       we keep checking whether the halfway point is greater or
       smaller than what we get out of the BAR averages.

       For the comparison we can use twice the tolerance. */
    while (DG2 - DG0 > 2*tol)
    {
        DG1 = 0.5*(DG0 + DG2);

        /*printf("Wfac1=%g, Wfac2=%g, beta=%g, DG1=%g\n",Wfac1,Wfac2,beta,
          DG1);*/

        /* calculate the BAR averages */
        dDG1=0.;

        for(i=0;i<ca->nsamples;i++)
        {
            samples_t *s=ca->s[i];
            sample_range_t *r=&(ca->r[i]);
            if (r->use)
            {
                if (s->hist)
                {
                    dDG1 += calc_bar_sum_hist(s->hist, Wfac1, (M-DG1), type);
                }
                else
                {
                    dDG1 += calc_bar_sum(r->end - r->start, s->du + r->start,
                                         Wfac1, (M-DG1));
                }
            }
        }
        for(i=0;i<cb->nsamples;i++)
        {
            samples_t *s=cb->s[i];
            sample_range_t *r=&(cb->r[i]);
            if (r->use)
            {
                if (s->hist)
                {
                    dDG1 -= calc_bar_sum_hist(s->hist, Wfac2, -(M-DG1), type);
                }
                else
                {
                    dDG1 -= calc_bar_sum(r->end - r->start, s->du + r->start,
                                         Wfac2, -(M-DG1));
                }
            }
        }
       
        if (dDG1 < 0)
        {
            DG0 = DG1;
        }
        else
        {
            DG2 = DG1;
        }
        if (debug)
        {
            fprintf(debug,"DG %9.5f %9.5f\n",DG0,DG2);
        }
    }
    
    return 0.5*(DG0 + DG2);
}

static void calc_rel_entropy(sample_coll_t *ca, sample_coll_t *cb,
                             double temp, double dg, double *sa, double *sb)
{
    int i,j;
    double W_ab=0.;
    double W_ba=0.;
    double kT, beta;
    double Wfac1, Wfac2;
    double n1,n2;

    kT   = BOLTZ*temp;
    beta = 1/kT;

    /* count the numbers of samples */ 
    n1 = ca->ntot;
    n2 = cb->ntot;

    /* to ensure the work values are the same as during the delta_G */
    if (!lambda_same(ca->native_lambda, ca->foreign_lambda))
    {
        /* this is the case when the delta U were calculated directly
           (i.e. we're not scaling dhdl) */
        Wfac1 = beta;
        Wfac2 = beta;
    }
    else
    {
        /* we're using dhdl, so delta_lambda needs to be a 
           multiplication factor.  */
        double delta_lambda=cb->native_lambda-ca->native_lambda;
        Wfac1 =  beta*delta_lambda;
        Wfac2 = -beta*delta_lambda;
    }

    /* first calculate the average work in both directions */
    for(i=0;i<ca->nsamples;i++)
    {
        samples_t *s=ca->s[i];
        sample_range_t *r=&(ca->r[i]);
        if (r->use)
        {
            if (!s->hist)
            {
                for(j=r->start;j<r->end;j++)
                    W_ab += Wfac1*s->du[j];
            }
            else
            {
                /* normalization factor multiplied with bin width and
                   number of samples (we normalize through M): */
                double normdx = 1.;
                double dx;
                int hd=0; /* histogram direction */
                if ( (s->hist->nhist>1) && (Wfac1<0) )
                {
                    hd=1;
                }
                dx=s->hist->dx[hd];

                for(j=0;j<s->hist->nbin[0];j++)
                {
                    double x=Wfac1*((j+s->hist->x0[0])+0.5)*dx; /*bin ctr*/
                    double pxdx=s->hist->bin[0][j]*normdx; /* p(x)dx */
                    W_ab += pxdx*x;
                }
            }
        }
    }
    W_ab/=n1;

    for(i=0;i<cb->nsamples;i++)
    {
        samples_t *s=cb->s[i];
        sample_range_t *r=&(cb->r[i]);
        if (r->use)
        {
            if (!s->hist)
            {
                for(j=r->start;j<r->end;j++)
                    W_ba += Wfac1*s->du[j];
            }
            else
            {
                /* normalization factor multiplied with bin width and
                   number of samples (we normalize through M): */
                double normdx = 1.;
                double dx;
                int hd=0; /* histogram direction */
                if ( (s->hist->nhist>1) && (Wfac2<0) )
                {
                    hd=1;
                }
                dx=s->hist->dx[hd];

                for(j=0;j<s->hist->nbin[0];j++)
                {
                    double x=Wfac1*((j+s->hist->x0[0])+0.5)*dx;/*bin ctr*/
                    double pxdx=s->hist->bin[0][j]*normdx; /* p(x)dx */
                    W_ba += pxdx*x;
                }
            }
        }
    }
    W_ba/=n2;
   
    /* then calculate the relative entropies */
    *sa = (W_ab - dg);
    *sb = (W_ba + dg);
}

static void calc_dg_stddev(sample_coll_t *ca, sample_coll_t *cb,
                           double temp, double dg, double *stddev)
{
    int i,j;
    double M;
    double sigmafact=0.;
    double kT, beta;
    double Wfac1, Wfac2;
    double n1, n2;

    kT   = BOLTZ*temp;
    beta = 1/kT;

    /* count the numbers of samples */ 
    n1 = ca->ntot;
    n2 = cb->ntot;

    /* to ensure the work values are the same as during the delta_G */
    if (!lambda_same(ca->native_lambda, ca->foreign_lambda))
    {
        /* this is the case when the delta U were calculated directly
           (i.e. we're not scaling dhdl) */
        Wfac1 = beta;
        Wfac2 = beta;
    }
    else
    {
        /* we're using dhdl, so delta_lambda needs to be a 
           multiplication factor.  */
        double delta_lambda=cb->native_lambda-ca->native_lambda;
        Wfac1 =  beta*delta_lambda;
        Wfac2 = -beta*delta_lambda;
    }

    M = log(n1/n2);


    /* calculate average in both directions */
    for(i=0;i<ca->nsamples;i++)
    {
        samples_t *s=ca->s[i];
        sample_range_t *r=&(ca->r[i]);
        if (r->use)
        {
            if (!s->hist)
            {
                for(j=r->start;j<r->end;j++)
                {
                    sigmafact += 1./(2. + 2.*cosh((M + Wfac1*s->du[j] - dg)));
                }
            }
            else
            {
                /* normalization factor multiplied with bin width and
                   number of samples (we normalize through M): */
                double normdx = 1.;
                double dx;
                int hd=0; /* histogram direction */
                if ( (s->hist->nhist>1) && (Wfac1<0) )
                {
                    hd=1;
                }
                dx=s->hist->dx[hd];

                for(j=0;j<s->hist->nbin[0];j++)
                {
                    double x=Wfac1*((j+s->hist->x0[0])+0.5)*dx; /*bin ctr*/
                    double pxdx=s->hist->bin[0][j]*normdx; /* p(x)dx */

                    sigmafact += pxdx/(2. + 2.*cosh((M + x - dg)));
                }
            }
        }
    }
    for(i=0;i<cb->nsamples;i++)
    {
        samples_t *s=cb->s[i];
        sample_range_t *r=&(cb->r[i]);
        if (r->use)
        {
            if (!s->hist)
            {
                for(j=r->start;j<r->end;j++)
                {
                    sigmafact += 1./(2. + 2.*cosh((M - Wfac2*s->du[j] - dg)));
                }
            }
            else
            {
                /* normalization factor multiplied with bin width and
                   number of samples (we normalize through M): */
                double normdx = 1.;
                double dx;
                int hd=0; /* histogram direction */
                if ( (s->hist->nhist>1) && (Wfac2<0) )
                {
                    hd=1;
                }
                dx=s->hist->dx[hd];

                for(j=0;j<s->hist->nbin[0];j++)
                {
                    double x=Wfac2*((j+s->hist->x0[0])+0.5)*dx;/*bin ctr*/
                    double pxdx=s->hist->bin[0][j]*normdx; /* p(x)dx */

                    sigmafact += pxdx/(2. + 2.*cosh((M - x - dg)));
                }
            }
        }
    }

    sigmafact /= (n1 + n2);
 
  
    /* Eq. 10 from 
       Shirts, Bair, Hooker & Pande, Phys. Rev. Lett 91, 140601 (2003): */
    *stddev = sqrt(((1./sigmafact) - ( (n1+n2)/n1 + (n1+n2)/n2 )));
}



static void calc_bar(barres_t *br, double tol, 
                     int npee_min, int npee_max, gmx_bool *bEE, 
                     double *partsum)
{
    int npee,p;
    double dg_sig2,sa_sig2,sb_sig2,stddev_sig2; /* intermediate variance values
                                                   for calculated quantities */
    int nsample1, nsample2;
    double temp=br->a->temp;
    int i,j;
    double dg_min, dg_max;
    gmx_bool have_hist=FALSE;

    br->dg=calc_bar_lowlevel(br->a, br->b, temp, tol, 0);

    br->dg_disc_err = 0.;
    br->dg_histrange_err = 0.;

    /* check if there are histograms */
    for(i=0;i<br->a->nsamples;i++)
    {
        if (br->a->r[i].use && br->a->s[i]->hist)
        {
            have_hist=TRUE;
            break;
        }
    }
    if (!have_hist)
    {
        for(i=0;i<br->b->nsamples;i++)
        {
            if (br->b->r[i].use && br->b->s[i]->hist)
            {
                have_hist=TRUE;
                break;
            }
        }
    }

    /* calculate histogram-specific errors */
    if (have_hist)
    {
        dg_min=calc_bar_lowlevel(br->a, br->b, temp, tol, -1);
        dg_max=calc_bar_lowlevel(br->a, br->b, temp, tol, 1);

        if (fabs(dg_max - dg_min) > GMX_REAL_EPS*10)
        {
            /* the histogram range  error is the biggest of the differences 
               between the best estimate and the extremes */
            br->dg_histrange_err = fabs(dg_max - dg_min);
        }
        br->dg_disc_err = 0.;
        for(i=0;i<br->a->nsamples;i++)
        {
            if (br->a->s[i]->hist)
                br->dg_disc_err=max(br->dg_disc_err, br->a->s[i]->hist->dx[0]);
        }
        for(i=0;i<br->b->nsamples;i++)
        {
            if (br->b->s[i]->hist)
                br->dg_disc_err=max(br->dg_disc_err, br->b->s[i]->hist->dx[0]);
        }
    }
    calc_rel_entropy(br->a, br->b, temp, br->dg, &(br->sa), &(br->sb));
                     
    calc_dg_stddev(br->a, br->b, temp, br->dg, &(br->dg_stddev) );

    dg_sig2 = 0;
    sa_sig2 = 0;
    sb_sig2 = 0;
    stddev_sig2 = 0;

    *bEE=TRUE;
    {
        sample_coll_t ca, cb;

        /* initialize the samples */
        sample_coll_init(&ca, br->a->native_lambda, br->a->foreign_lambda, 
                        br->a->temp);
        sample_coll_init(&cb, br->b->native_lambda, br->b->foreign_lambda, 
                        br->b->temp);

        for(npee=npee_min; npee<=npee_max; npee++)
        {
            double dgs      = 0;
            double dgs2     = 0;
            double dsa      = 0;
            double dsb      = 0;
            double dsa2     = 0;
            double dsb2     = 0;
            double dstddev  = 0;
            double dstddev2 = 0;
 
            
            for(p=0; p<npee; p++)
            {
                double dgp;
                double stddevc;
                double sac, sbc;
                gmx_bool cac, cbc;

                cac=sample_coll_create_subsample(&ca, br->a, p, npee);
                cbc=sample_coll_create_subsample(&cb, br->b, p, npee);

                if (!cac || !cbc)
                {
                    printf("WARNING: histogram number incompatible with block number for averaging: can't do error estimate\n");
                    *bEE=FALSE;
                    if (cac)
                        sample_coll_destroy(&ca);
                    if (cbc)
                        sample_coll_destroy(&cb);
                    return;
                }

                dgp=calc_bar_lowlevel(&ca, &cb, temp, tol, 0);
                dgs  += dgp;
                dgs2 += dgp*dgp;

                partsum[npee*(npee_max+1)+p] += dgp;

                calc_rel_entropy(&ca, &cb, temp, dgp, &sac, &sbc); 
                dsa  += sac;
                dsa2 += sac*sac;
                dsb  += sbc;
                dsb2 += sbc*sbc;
                calc_dg_stddev(&ca, &cb, temp, dgp, &stddevc );

                dstddev  += stddevc;
                dstddev2 += stddevc*stddevc;

                sample_coll_destroy(&ca);
                sample_coll_destroy(&cb);
            }
            dgs  /= npee;
            dgs2 /= npee;
            dg_sig2 += (dgs2-dgs*dgs)/(npee-1);

            dsa  /= npee;
            dsa2 /= npee;
            dsb  /= npee;
            dsb2 /= npee;
            sa_sig2 += (dsa2-dsa*dsa)/(npee-1);
            sb_sig2 += (dsb2-dsb*dsb)/(npee-1);

            dstddev  /= npee;
            dstddev2 /= npee;
            stddev_sig2 += (dstddev2-dstddev*dstddev)/(npee-1);
        }
        br->dg_err = sqrt(dg_sig2/(npee_max - npee_min + 1));
        br->sa_err = sqrt(sa_sig2/(npee_max - npee_min + 1));
        br->sb_err = sqrt(sb_sig2/(npee_max - npee_min + 1));
        br->dg_stddev_err = sqrt(stddev_sig2/(npee_max - npee_min + 1));
    }
}


static double bar_err(int nbmin, int nbmax, const double *partsum)
{
    int nb,b;
    double svar,s,s2,dg;

    svar = 0;
    for(nb=nbmin; nb<=nbmax; nb++)
    {
        s  = 0;
        s2 = 0;
        for(b=0; b<nb; b++)
        {
            dg  = partsum[nb*(nbmax+1)+b];
            s  += dg;
            s2 += dg*dg;
        }
        s  /= nb;
        s2 /= nb;
        svar += (s2 - s*s)/(nb - 1);
    }

    return sqrt(svar/(nbmax + 1 - nbmin));
}

/* deduce lambda value from legend. 
input:
    bdhdl = if true, value may be a derivative. 
output:
    bdhdl = whether the legend was for a derivative.
    */
static double legend2lambda(char *fn,const char *legend,gmx_bool *bdhdl)
{
    double lambda=0;
    const char   *ptr;
    gmx_bool ok=FALSE;

    if (legend == NULL)
    {
        gmx_fatal(FARGS,"There is no legend in file '%s', can not deduce lambda",fn);
    }
    ptr = strrchr(legend,' ');

    if (strstr(legend,"dH"))
    {
        if (! (*bdhdl))
        {
            ok=FALSE;
        }
        else
        {
            ok=TRUE;
        }
    }
    else
    {
        if (strchr(legend,'D') != NULL && strchr(legend,'H') != NULL)
        {
            ok=TRUE;
            *bdhdl=FALSE;
        }
    }
    if (!ptr)
    {
        ok=FALSE;
    }

    if (!ok)
    {
        printf("%s\n", legend);
        gmx_fatal(FARGS,"There is no proper lambda legend in file '%s', can not deduce lambda",fn);
    }
    if (sscanf(ptr,"%lf",&lambda) != 1)
    {
        gmx_fatal(FARGS,"There is no proper lambda legend in file '%s', can not deduce lambda",fn);
    }

    return lambda;
}

static gmx_bool subtitle2lambda(const char *subtitle,double *lambda)
{
    gmx_bool bFound;
    char *ptr;

    bFound = FALSE;

    /* plain text lambda string */
    ptr = strstr(subtitle,"lambda");
    if (ptr == NULL)
    {
        /* xmgrace formatted lambda string */
        ptr = strstr(subtitle,"\\xl\\f{}");
    }
    if (ptr == NULL)
    {
        /* xmgr formatted lambda string */
        ptr = strstr(subtitle,"\\8l\\4");
    }
    if (ptr != NULL)
    {
        ptr = strstr(ptr,"=");
    }
    if (ptr != NULL)
    {
        bFound = (sscanf(ptr+1,"%lf",lambda) == 1);
    }

    return bFound;
}

static double filename2lambda(char *fn)
{
    double lambda;
    char   *ptr,*endptr,*digitptr;
    int     dirsep;
    ptr = fn;
    /* go to the end of the path string and search backward to find the last 
       directory in the path which has to contain the value of lambda 
     */
    while (ptr[1] != '\0')
    {
        ptr++;
    }
    /* searching backward to find the second directory separator */
    dirsep = 0;
    digitptr = NULL;
    while (ptr >= fn)
    {
        if (ptr[0] != DIR_SEPARATOR && ptr[1] == DIR_SEPARATOR)
        {            
            if (dirsep == 1) break;
            dirsep++;
        }
        /* save the last position of a digit between the last two 
           separators = in the last dirname */
        if (dirsep > 0 && isdigit(*ptr))
        {
            digitptr = ptr;
        }
        ptr--;
    }
    if (!digitptr)
    {
        gmx_fatal(FARGS,"While trying to read the lambda value from the file path:"
                    " last directory in the path '%s' does not contain a number",fn);
    }
    if (digitptr[-1] == '-')
    {
        digitptr--;
    }
    lambda = strtod(digitptr,&endptr);
    if (endptr == digitptr)
    {
        gmx_fatal(FARGS,"Malformed number in file path '%s'",fn);
    }

    return lambda;
}

static void read_bar_xvg_lowlevel(char *fn, real *temp, xvg_t *ba)
{
    int  i;
    char *subtitle,**legend,*ptr;
    int np;
    gmx_bool native_lambda_read=FALSE;

    xvg_init(ba);

    ba->filename = fn;

    np = read_xvg_legend(fn,&ba->y,&ba->nset,&subtitle,&legend);
    if (!ba->y)
    {
        gmx_fatal(FARGS,"File %s contains no usable data.",fn);
    }
    ba->t  = ba->y[0];

    snew(ba->np,ba->nset);
    for(i=0;i<ba->nset;i++)
        ba->np[i]=np;


    ba->temp = -1;
    if (subtitle != NULL)
    {
        ptr = strstr(subtitle,"T =");
        if (ptr != NULL)
        {
            ptr += 3;
            if (sscanf(ptr,"%lf",&ba->temp) == 1)
            {
                if (ba->temp <= 0)
                {
                    gmx_fatal(FARGS,"Found temperature of %g in file '%s'",
                              ba->temp,fn);
                }
            }
        }
    }
    if (ba->temp < 0)
    {
        if (*temp <= 0)
        {
            gmx_fatal(FARGS,"Did not find a temperature in the subtitle in file '%s', use the -temp option of g_bar",fn);
        }
        ba->temp = *temp;
    }

    /* Try to deduce lambda from the subtitle */
    if (subtitle)
    {
        if (subtitle2lambda(subtitle,&(ba->native_lambda)))
        {
            native_lambda_read=TRUE;
        }
    }
    snew(ba->lambda,ba->nset-1);
    if (legend == NULL)
    {
        /* Check if we have a single set, no legend, nset=2 means t and dH/dl */
        if (ba->nset == 2)
        {
            if (!native_lambda_read)
            {
                /* Deduce lambda from the file name */
                ba->native_lambda = filename2lambda(fn);
                native_lambda_read=TRUE;
            }
            ba->lambda[0] = ba->native_lambda;
        }
        else
        {
            gmx_fatal(FARGS,"File %s contains multiple sets but no legends, can not determine the lambda values",fn);
        }
    }
    else
    {
        for(i=0; i<ba->nset-1; i++)
        {
            gmx_bool is_dhdl=(i==0);
            /* Read lambda from the legend */
            ba->lambda[i] = legend2lambda(fn,legend[i], &is_dhdl);

            if (is_dhdl && !native_lambda_read)
            {
                ba->native_lambda = ba->lambda[i];
                native_lambda_read=TRUE;
            }
        }
    }

    if (!native_lambda_read)
    {
        gmx_fatal(FARGS,"File %s contains multiple sets but no indication of the native lambda",fn);
    }
    
    /* Reorder the data */
    for(i=1; i<ba->nset; i++)
    {
        ba->y[i-1] = ba->y[i];
    }
    if (legend != NULL)
    {
        for(i=0; i<ba->nset-1; i++)
        {
            sfree(legend[i]);
        }
        sfree(legend);
    }
    ba->nset--;
}

static void read_bar_xvg(char *fn, real *temp, lambda_t *lambda_head)
{
    xvg_t *barsim;
    samples_t *s;
    int i;
    double *lambda;

    snew(barsim,1);

    read_bar_xvg_lowlevel(fn, temp, barsim);

    if (barsim->nset <1 )
    {
        gmx_fatal(FARGS,"File '%s' contains fewer than two columns", fn);
    }

    if ( ( *temp != barsim->temp) && (*temp > 0) )
    {
        gmx_fatal(FARGS,"Temperature in file %s different from earlier files or setting\n", fn);
    }
    *temp=barsim->temp;

    /* now create a series of samples_t */
    snew(s, barsim->nset);
    for(i=0;i<barsim->nset;i++)
    {
        samples_init(s+i, barsim->native_lambda, barsim->lambda[i], 
                     barsim->temp, lambda_same(barsim->native_lambda,
                                               barsim->lambda[i]), 
                     fn);
        s[i].du=barsim->y[i];
        s[i].ndu=barsim->np[i];
        s[i].t=barsim->t;

        lambda_list_insert_sample(lambda_head, s+i);
    }
    printf("%s: %.1f - %.1f; lambda = %.3f\n    foreign lambdas:", 
           fn, s[0].t[0], s[0].t[s[0].ndu-1], s[0].native_lambda);
    for(i=0;i<barsim->nset;i++)
    {
        printf(" %.3f (%d pts)", s[i].foreign_lambda, s[i].ndu);
    }
    printf("\n\n");
}

static void read_edr_rawdh_block(samples_t **smp, int *ndu, t_enxblock *blk, 
                                 double start_time, double delta_time,
                                 double native_lambda, double temp,
                                 double *last_t, const char *filename)
{
    int j;
    gmx_bool allocated;
    double foreign_lambda;
    int derivative;
    samples_t *s; /* convenience pointer */
    int startj;

    /* check the block types etc. */
    if ( (blk->nsub < 3) ||
         (blk->sub[0].type != xdr_datatype_int) ||
         (blk->sub[1].type != xdr_datatype_double) ||
         (
          (blk->sub[2].type != xdr_datatype_float) &&
          (blk->sub[2].type != xdr_datatype_double) 
         ) ||
         (blk->sub[0].nr < 1) ||
         (blk->sub[1].nr < 1) )
    {
        gmx_fatal(FARGS, 
                  "Unexpected/corrupted block data in file %s around time %g.", 
                  filename, start_time);
    }
   
    derivative = blk->sub[0].ival[0]; 
    foreign_lambda = blk->sub[1].dval[0];

    if (! *smp)
    {
        /* initialize the samples structure if it's empty. */
        snew(*smp, 1);
        samples_init(*smp, native_lambda, foreign_lambda, temp,
                     derivative!=0, filename);
        (*smp)->start_time=start_time;
        (*smp)->delta_time=delta_time;
    }

    /* set convenience pointer */
    s=*smp;

    /* now double check */
    if ( ! lambda_same(s->foreign_lambda, foreign_lambda) ||
         (  (derivative!=0) != (s->derivative!=0) ) )
    {
        fprintf(stderr, "Got foreign lambda=%g, expected: %g\n", 
                foreign_lambda, s->foreign_lambda);
        fprintf(stderr, "Got derivative=%d, expected: %d\n", 
                derivative, s->derivative);
        gmx_fatal(FARGS, "Corrupted data in file %s around t=%g.", 
                  filename, start_time);
    }

    /* make room for the data */
    if (s->ndu_alloc < (s->ndu + blk->sub[2].nr) )
    {
        s->ndu_alloc += (s->ndu_alloc < blk->sub[2].nr) ?  
                            blk->sub[2].nr*2 : s->ndu_alloc;
        srenew(s->du_alloc, s->ndu_alloc);
        s->du=s->du_alloc;
    }
    startj = s->ndu;
    s->ndu += blk->sub[2].nr;
    s->ntot += blk->sub[2].nr;
    *ndu = blk->sub[2].nr;

    /* and copy the data*/
    for(j=0;j<blk->sub[2].nr;j++)
    {
        if (blk->sub[2].type == xdr_datatype_float)
        {
            s->du[startj+j] = blk->sub[2].fval[j];
        }
        else
        {
            s->du[startj+j] = blk->sub[2].dval[j];
        }
    }
    if (start_time + blk->sub[2].nr*delta_time > *last_t)
    {
        *last_t = start_time + blk->sub[2].nr*delta_time;
    }
}

static samples_t *read_edr_hist_block(int *nsamples, t_enxblock *blk,
                                      double start_time, double delta_time,
                                      double native_lambda, double temp,
                                      double *last_t, const char *filename)
{
    int i,j;
    samples_t *s;
    int nhist;
    double foreign_lambda;
    int derivative;
    int nbins[2];

    /* check the block types etc. */
    if ( (blk->nsub < 2) ||
         (blk->sub[0].type != xdr_datatype_double) ||
         (blk->sub[1].type != xdr_datatype_large_int) ||
         (blk->sub[0].nr < 2)  ||
         (blk->sub[1].nr < 2) )
    {
        gmx_fatal(FARGS, 
                  "Unexpected/corrupted block data in file %s around time %g", 
                  filename, start_time);
    }

    nhist=blk->nsub-2;
    if (nhist == 0)
    {
        return NULL;
    }
    if (nhist > 2)
    {
        gmx_fatal(FARGS, 
                  "Unexpected/corrupted block data in file %s around time %g", 
                  filename, start_time);
    }

    snew(s, 1);
    *nsamples=1;

    foreign_lambda=blk->sub[0].dval[0];
    derivative=(int)(blk->sub[1].lval[1]);
    if (derivative)
        foreign_lambda=native_lambda;

    samples_init(s, native_lambda, foreign_lambda, temp,
                 derivative!=0, filename);
    snew(s->hist, 1);

    for(i=0;i<nhist;i++)
    {
        nbins[i] = blk->sub[i+2].nr;
    }

    hist_init(s->hist, nhist, nbins);

    for(i=0;i<nhist;i++)
    {
        s->hist->x0[i]=blk->sub[1].lval[2+i];
        s->hist->dx[i] = blk->sub[0].dval[1];
        if (i==1)
            s->hist->dx[i] = - s->hist->dx[i];
    }

    s->hist->start_time = start_time;
    s->hist->delta_time = delta_time;
    s->start_time = start_time;
    s->delta_time = delta_time;

    for(i=0;i<nhist;i++)
    {
        int nbin;
        gmx_large_int_t sum=0;

        for(j=0;j<s->hist->nbin[i];j++)
        { 
            int binv=(int)(blk->sub[i+2].ival[j]);

            s->hist->bin[i][j] = binv;
            sum += binv;

        }
        if (i==0)
        {
            s->ntot = sum;
            s->hist->sum = sum;
        }
        else
        {
            if (s->ntot != sum)
            {
                gmx_fatal(FARGS, "Histogram counts don't match in %s", 
                          filename);
            }
        }
    }

    if (start_time + s->hist->sum*delta_time > *last_t)
    {
        *last_t = start_time + s->hist->sum*delta_time;
    }
    return s;
}


static void read_barsim_edr(char *fn, real *temp, lambda_t *lambda_head)
{
    int i;
    ener_file_t fp;
    t_enxframe *fr; 
    int nre;
    gmx_enxnm_t *enm=NULL;
    double first_t=-1;
    double last_t=-1;
    samples_t **samples_rawdh=NULL; /* contains samples for raw delta_h  */
    int *nhists=NULL;       /* array to keep count & print at end */
    int *npts=NULL;         /* array to keep count & print at end */
    double *lambdas=NULL;   /* array to keep count & print at end */
    double native_lambda=-1;
    double end_time;        /* the end time of the last batch of samples */
    int nsamples=0;

    fp = open_enx(fn,"r");
    do_enxnms(fp,&nre,&enm);
    snew(fr, 1);

    while(do_enx(fp, fr))
    {
        /* count the data blocks */
        int nblocks_raw=0;
        int nblocks_hist=0;
        int nlam=0;
        int k;
        /* DHCOLL block information: */
        double start_time=0, delta_time=0, start_lambda=0, delta_lambda=0;
        double rtemp=0;

        /* count the blocks and handle collection information: */
        for(i=0;i<fr->nblock;i++)
        {
            if (fr->block[i].id == enxDHHIST )
                nblocks_hist++;
            if ( fr->block[i].id == enxDH )
                nblocks_raw++;
            if (fr->block[i].id == enxDHCOLL )
            {
                nlam++;
                if ( (fr->block[i].nsub < 1) || 
                     (fr->block[i].sub[0].type != xdr_datatype_double) ||
                     (fr->block[i].sub[0].nr < 5))
                {
                    gmx_fatal(FARGS, "Unexpected block data in file %s", fn);
                }

                /* read the data from the DHCOLL block */
                rtemp =        fr->block[i].sub[0].dval[0];
                start_time =   fr->block[i].sub[0].dval[1];
                delta_time =   fr->block[i].sub[0].dval[2];
                start_lambda = fr->block[i].sub[0].dval[3];
                delta_lambda = fr->block[i].sub[0].dval[4];

                if (delta_lambda>0)
                {
                    gmx_fatal(FARGS, "Lambda values not constant in %s: can't apply BAR method", fn);
                }
                if ( ( *temp != rtemp) && (*temp > 0) )
                {
                    gmx_fatal(FARGS,"Temperature in file %s different from earlier files or setting\n", fn);
                }
                *temp=rtemp;

                if (first_t < 0)
                    first_t=start_time;
            }
        }

        if (nlam != 1)
        {
            gmx_fatal(FARGS, "Did not find a delta h information in file %s" , fn);
        }
        if (nblocks_raw > 0 && nblocks_hist > 0 )
        {
            gmx_fatal(FARGS, "Can't handle both raw delta U data and histograms in the same file %s", fn);
        }

        if (nsamples > 0)
        {
            /* check the native lambda */
            if (!lambda_same(start_lambda, native_lambda) )
            {
                gmx_fatal(FARGS, "Native lambda not constant in file %s: started at %g, and becomes %g at time %g", 
                          fn, native_lambda, start_lambda, start_time);
            }
            /* check the number of samples against the previous number */
            if ( ((nblocks_raw+nblocks_hist)!=nsamples) || (nlam!=1 ) )
            {
                gmx_fatal(FARGS, "Unexpected block count in %s: was %d, now %d\n",
                          fn, nsamples+1, nblocks_raw+nblocks_hist+nlam);
            }
            /* check whether last iterations's end time matches with 
               the currrent start time */
            if ( (fabs(last_t - start_time) > 2*delta_time)  && last_t>=0)
            {
                /* it didn't. We need to store our samples and reallocate */
                for(i=0;i<nsamples;i++)
                {
                    if (samples_rawdh[i])
                    {
                        /* insert it into the existing list */
                        lambda_list_insert_sample(lambda_head, 
                                                  samples_rawdh[i]);
                        /* and make sure we'll allocate a new one this time
                           around */
                        samples_rawdh[i]=NULL;
                    }
                }
            }
        }
        else
        {
            /* this is the first round; allocate the associated data 
               structures */
            native_lambda=start_lambda;
            nsamples=nblocks_raw+nblocks_hist;
            snew(nhists, nsamples);
            snew(npts, nsamples);
            snew(lambdas, nsamples);
            snew(samples_rawdh, nsamples);
            for(i=0;i<nsamples;i++)
            {
                nhists[i]=0;
                npts[i]=0;
                lambdas[i]=-1;
                samples_rawdh[i]=NULL; /* init to NULL so we know which
                                          ones contain values */
            }
        }

        /* and read them */
        k=0; /* counter for the lambdas, etc. arrays */
        for(i=0;i<fr->nblock;i++)
        {
            if (fr->block[i].id == enxDH)
            {
                int ndu;
                read_edr_rawdh_block(&(samples_rawdh[k]),
                                     &ndu,
                                     &(fr->block[i]), 
                                     start_time, delta_time, 
                                     start_lambda, rtemp, 
                                     &last_t, fn);
                npts[k] += ndu;
                if (samples_rawdh[k])
                {
                    lambdas[k]=samples_rawdh[k]->foreign_lambda;
                }
                k++;
            }
            else if (fr->block[i].id == enxDHHIST)
            {
                int j;
                int nb=0;
                samples_t *s; /* this is where the data will go */
                s=read_edr_hist_block(&nb, &(fr->block[i]), 
                                      start_time, delta_time, 
                                      start_lambda, rtemp, 
                                      &last_t, fn);
                nhists[k] += nb;
                if (nb>0)
                {
                    lambdas[k]= s->foreign_lambda;
                }
                k++;
                /* and insert the new sample immediately */
                for(j=0;j<nb;j++)
                {
                    lambda_list_insert_sample(lambda_head, s+j);
                }
            }
        }
    }
    /* Now store all our extant sample collections */
    for(i=0;i<nsamples;i++)
    {
        if (samples_rawdh[i])
        {
            /* insert it into the existing list */
            lambda_list_insert_sample(lambda_head, samples_rawdh[i]);
        }
    }


    fprintf(stderr, "\n");
    printf("%s: %.1f - %.1f; lambda = %.3f\n    foreign lambdas:", 
           fn, first_t, last_t, native_lambda);
    for(i=0;i<nsamples;i++)
    {
        if (nhists[i] > 0)
        {
            printf(" %.3f (%d hists)", lambdas[i], nhists[i]);
        }
        else
        {
            printf(" %.3f (%d pts)", lambdas[i], npts[i]);
        }
    }
    printf("\n\n");
    sfree(npts);
    sfree(nhists);
    sfree(lambdas);
}


int gmx_bar(int argc,char *argv[])
{
    static const char *desc[] = {
        "g_bar calculates free energy difference estimates through ",
        "Bennett's acceptance ratio method (BAR). It also automatically",
        "adds series of individual free energies obtained with BAR into",
        "a combined free energy estimate.[PAR]",

        "Every individual BAR free energy difference relies on two ",
        "simulations at different states: say state A and state B, as",
        "controlled by a parameter 'lambda' (see the mdp parameter",
        "'init_lambda'). The BAR method calculates a ratio of weighted",
        "average of the Hamiltonian difference of state B given state A and",
        "vice versa. If the Hamiltonian does not linearly depend on lambda",
        "(in which case we can extrapolate the derivative of the Hamiltonian",
        "w.r.t. lambda, as is the default when 'free_energy' is on), the",
        "energy differences to the other state need to be calculated",
        "explicitly during the simulation. This can be controlled with",
        "the mdp option 'foreign_lambda'.[PAR]",

        "Input option [TT]-f[tt] expects multiple dhdl files. ",
        "Two types of input files are supported:[BR]",
        "* Files with only one y-value, for such files it is assumed ",
        "that the y-value is dH/dlambda and that the Hamiltonian depends ",
        "linearly on lambda. The lambda value of the simulation is inferred ",
        "from the subtitle if present, otherwise from a number in the",
        "subdirectory in the file name.",
        "[BR]",
        "* Files with more than one y-value. The files should have columns ",
        "with dH/dlambda and Delta lambda. The lambda values are inferred ",
        "from the legends: ",
        "lambda of the simulation from the legend of dH/dlambda ",
        "and the foreign lambda's from the legends of Delta H.[BR]",
        "The lambda of the simulation is parsed from dhdl.xvg file's legend ",
        "containing the string 'dH', the foreign lambda's from the legend ",
        "containing the capitalized letters 'D' and 'H'. The temperature ",
        "is parsed from the legend line containing 'T ='.[PAR]",

        "The input option [TT]-g[tt] expects multiple .edr files. ",
        "These can contain either lists of energy differences (see the",
        "mdp option separate_dhdl_file), or a series of histograms",
        "(see the mdp options dh_hist_size and dh_hist_spacing).",
        "The temperature and lambda values are automatically deduced from",
        "the ener.edr file.[PAR]"

        "The free energy estimates are determined using BAR with bisection, ",
        "the precision of the output is set with [TT]-prec[tt]. ",
        "An error estimate taking into account time correlations ",
        "is made by splitting the data into blocks and determining ",
        "the free energy differences over those blocks and assuming ",
        "the blocks are independent. ",
        "The final error estimate is determined from the average variance ",
        "over 5 blocks. A range of blocks numbers for error estimation can ",
        "be provided with the options [TT]-nbmin[tt] and [TT]-nbmax[tt].[PAR]",

        "g_bar tries to aggregate samples with the same 'native' and 'foreign'",
        "lambda values, but always assumes independent samples: note that",
        "when aggregating energy differences/derivatives with different",
        "sampling intervals, this is almost certainly not correct: usually",
        "subsequent energies are correlated and different time intervals mean",
        "different degrees of correlation between samples.[PAR]",

        "The results are split in two parts: the last part contains the final ",
        "results in kJ/mol, together with the error estimate for each part ",
        "and the total. The first part contains detailed free energy ",
        "difference estimates and phase space overlap measures in units of ",
        "kT (together with their computed error estimate). The printed ",
        "values are:[BR]",
        "*  lam_A: the lambda values for point A.[BR]",
        "*  lam_B: the lambda values for point B.[BR]",
        "*     DG: the free energy estimate.[BR]",
        "*    s_A: an estimate of the relative entropy of B in A.[BR]",
        "*    s_A: an estimate of the relative entropy of A in B.[BR]",
        "*  stdev: an estimate expected per-sample standard deviation.[PAR]",
        
        "The relative entropy of both states in each other's ensemble can be ",
        "interpreted as a measure of phase space overlap: ", 
        "the relative entropy s_A of the work samples of lambda_B in the ",
        "ensemble of lambda_A (and vice versa for s_B), is a ", 
        "measure of the 'distance' between Boltzmann distributions of ",
        "the two states, that goes to zero for identical distributions. See ",
        "Wu & Kofke, J. Chem. Phys. 123 084109 (2009) for more information.",
        "[PAR]",
        "The estimate of the expected per-sample standard deviation, as given ",
        "in Bennett's original BAR paper: ",
        "Bennett, J. Comp. Phys. 22, p 245 (1976), Eq. 10 gives an estimate ",
        "of the quality of sampling (not directly of the actual statistical ", 
        "error, because it assumes independent samples).[PAR]",

        "To get a visual estimate of the phase space overlap, use the ",
        "-oh option to write series of histograms, together with the ",
        "-nbin option.[PAR]"
    };
    static real begin=0,end=-1,temp=-1;
    int nd=2,nbmin=5,nbmax=5;
    int nbin=100;
    gmx_bool calc_s,calc_v;
    t_pargs pa[] = {
        { "-b",    FALSE, etREAL, {&begin},  "Begin time for BAR" },
        { "-e",    FALSE, etREAL, {&end},    "End time for BAR" },
        { "-temp", FALSE, etREAL, {&temp},   "Temperature (K)" },
        { "-prec", FALSE, etINT,  {&nd},     "The number of digits after the decimal point" },
        { "-nbmin",  FALSE, etINT,  {&nbmin}, "Minimum number of blocks for error estimation" },
        { "-nbmax",  FALSE, etINT,  {&nbmax}, "Maximum number of blocks for error estimation" },
        { "-nbin",  FALSE, etINT, {&nbin}, "Number of bins for histogram output"}
    };
    
    t_filenm   fnm[] = {
        { efXVG, "-f",  "dhdl",   ffOPTRDMULT },
        { efEDR, "-g",  "ener",   ffOPTRDMULT },
        { efXVG, "-o",  "bar",    ffOPTWR },
        { efXVG, "-oi", "barint", ffOPTWR }, 
        { efXVG, "-oh", "histogram", ffOPTWR }
    };
#define NFILE asize(fnm)
    
    int      f,i,j;
    int      nf=0; /* file counter */
    int      nbs;
    int      nfile_tot; /* total number of input files */
    int      nxvgfile=0;
    int      nedrfile=0;
    char     **fxvgnms;
    char     **fedrnms;
    lambda_t *lb;    /* the pre-processed lambda data (linked list head) */
    lambda_t lambda_head; /* the head element */
    barres_t *results;  /* the results */
    int    nresults;  /* number of results in results array */

    double   *partsum;
    double   prec,dg_tot,dg,sig, dg_tot_max, dg_tot_min;
    FILE     *fpb,*fpi;
    char     lamformat[20];
    char     dgformat[20],xvg2format[STRLEN],xvg3format[STRLEN],buf[STRLEN];
    char     ktformat[STRLEN], sktformat[STRLEN];
    char     kteformat[STRLEN], skteformat[STRLEN];
    output_env_t oenv;
    double   kT, beta;
    gmx_bool     result_OK=TRUE,bEE=TRUE;

    gmx_bool     disc_err=FALSE;
    double   sum_disc_err=0.; /* discretization error */
    gmx_bool     histrange_err=FALSE;
    double   sum_histrange_err=0.; /* histogram range error */
    double   stat_err=0.; /* statistical error */
    
    CopyRight(stderr,argv[0]);
    parse_common_args(&argc,argv,
                      PCA_CAN_VIEW,
                      NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv);

    if (opt2bSet("-f", NFILE, fnm))
    {
        nxvgfile = opt2fns(&fxvgnms,"-f",NFILE,fnm);
    }
    if (opt2bSet("-g", NFILE, fnm))
    {
        nedrfile = opt2fns(&fedrnms,"-g",NFILE,fnm);
    }

    /* make linked list */
    lb=&lambda_head;
    lambda_init(lb, 0, 0);
    lb->next=lb;
    lb->prev=lb;


    nfile_tot = nxvgfile + nedrfile;

    if (nfile_tot == 0)
    {
        gmx_fatal(FARGS,"No input files!");
    }

    if (nd < 0)
    {
        gmx_fatal(FARGS,"Can not have negative number of digits");
    }
    prec = pow(10,-nd);

    snew(partsum,(nbmax+1)*(nbmax+1));
    nf = 0;

    /* read in all files. First xvg files */
    for(f=0; f<nxvgfile; f++)
    {
        read_bar_xvg(fxvgnms[f],&temp,lb);
        nf++;
    }
    /* then .edr files */
    for(f=0; f<nedrfile; f++)
    {
        read_barsim_edr(fedrnms[f],&temp,lb);;
        nf++;
    }

    /* fix the times to allow for equilibration */
    lambdas_impose_times(lb, begin, end);

    if (opt2bSet("-oh",NFILE,fnm))
    {
        lambdas_histogram(lb, opt2fn("-oh",NFILE,fnm), nbin, oenv);
    }
   
    /* assemble the output structures from the lambdas */
    results=barres_list_create(lb, &nresults);

    sum_disc_err=barres_list_max_disc_err(results, nresults);

    if (nresults == 0)
    {
        printf("\nNo results to calculate.\n");
        return 0;
    }

#if 1
    if (sum_disc_err > prec)
    {
        prec=sum_disc_err;
        nd = ceil(-log10(prec));
        printf("WARNING: setting the precision to %g because that is the minimum\n         reasonable number, given the expected discretization error.\n", prec);
    }
#endif

    sprintf(lamformat,"%%6.3f");
    sprintf( dgformat,"%%%d.%df",3+nd,nd);
    /* the format strings of the results in kT */
    sprintf( ktformat,"%%%d.%df",5+nd,nd);
    sprintf( sktformat,"%%%ds",6+nd);
    /* the format strings of the errors in kT */
    sprintf( kteformat,"%%%d.%df",3+nd,nd);
    sprintf( skteformat,"%%%ds",4+nd);
    sprintf(xvg2format,"%s %s\n","%g",dgformat);
    sprintf(xvg3format,"%s %s %s\n","%g",dgformat,dgformat);



    fpb = NULL;
    if (opt2bSet("-o",NFILE,fnm))
    {
        sprintf(buf,"%s (%s)","\\DeltaG",unit_energy);
        fpb = xvgropen_type(opt2fn("-o",NFILE,fnm),"Free energy differences",
                            "\\lambda",buf,exvggtXYDY,oenv);
    }
    
    fpi = NULL;
    if (opt2bSet("-oi",NFILE,fnm))
    {
        sprintf(buf,"%s (%s)","\\DeltaG",unit_energy);
        fpi = xvgropen(opt2fn("-oi",NFILE,fnm),"Free energy integral",
                      "\\lambda",buf,oenv);
    }



    if (nbmin > nbmax)
        nbmin=nbmax;

    /* first calculate results */
    bEE = TRUE;
    disc_err = FALSE;
    for(f=0; f<nresults; f++)
    {
        /* Determine the free energy difference with a factor of 10
         * more accuracy than requested for printing.
         */
        calc_bar(&(results[f]), 0.1*prec, nbmin, nbmax,
                 &bEE, partsum);

        if (results[f].dg_disc_err > prec/10.)
            disc_err=TRUE;
        if (results[f].dg_histrange_err > prec/10.)
            histrange_err=TRUE;
    }

    /* print results in kT */
    kT   = BOLTZ*temp;
    beta = 1/kT;

    printf("\nTemperature: %g K\n", temp);

    printf("\nDetailed results in kT (see help for explanation):\n\n");
    printf("%6s ", " lam_A");
    printf("%6s ", " lam_B");
    printf(sktformat,  "DG ");
    if (bEE)
        printf(skteformat, "+/- ");
    if (disc_err)
        printf(skteformat, "disc ");
    if (histrange_err)
        printf(skteformat, "range ");
    printf(sktformat,  "s_A ");
    if (bEE)
        printf(skteformat, "+/- " );
    printf(sktformat,  "s_B ");
    if (bEE)
        printf(skteformat, "+/- " );
    printf(sktformat,  "stdev ");
    if (bEE)
        printf(skteformat, "+/- ");
    printf("\n");
    for(f=0; f<nresults; f++)
    {
        printf(lamformat, results[f].a->native_lambda);
        printf(" ");
        printf(lamformat, results[f].b->native_lambda);
        printf(" ");
        printf(ktformat,  results[f].dg);
        printf(" ");
        if (bEE)
        {
            printf(kteformat, results[f].dg_err);
            printf(" ");
        }
        if (disc_err)
        {
            printf(kteformat, results[f].dg_disc_err);
            printf(" ");
        }
        if (histrange_err)
        {
            printf(kteformat, results[f].dg_histrange_err);
            printf(" ");
        }
        printf(ktformat,  results[f].sa);
        printf(" ");
        if (bEE)
        {
            printf(kteformat, results[f].sa_err);
            printf(" ");
        }
        printf(ktformat,  results[f].sb);
        printf(" ");
        if (bEE)
        {
            printf(kteformat, results[f].sb_err);
            printf(" ");
        }
        printf(ktformat,  results[f].dg_stddev);
        printf(" ");
        if (bEE)
        {
            printf(kteformat, results[f].dg_stddev_err);
        }
        printf("\n");

        /* Check for negative relative entropy with a 95% certainty. */
        if (results[f].sa < -2*results[f].sa_err ||
            results[f].sb < -2*results[f].sb_err)
        {
            result_OK=FALSE;
        }
    }

    if (!result_OK)
    {
        printf("\nWARNING: Some of these results violate the Second Law of "
               "Thermodynamics: \n"
               "         This is can be the result of severe undersampling, or "
               "(more likely)\n" 
               "         there is something wrong with the simulations.\n");
    }
 

    /* final results in kJ/mol */
    printf("\n\nFinal results in kJ/mol:\n\n");
    dg_tot  = 0;
    for(f=0; f<nresults; f++)
    {
        
        if (fpi != NULL)
        {
            fprintf(fpi, xvg2format, results[f].a->native_lambda, dg_tot);
        }


        if (fpb != NULL)
        {
            fprintf(fpb, xvg3format,
                    0.5*(results[f].a->native_lambda + 
                         results[f].b->native_lambda),
                    results[f].dg,results[f].dg_err);
        }

        /*printf("lambda %4.2f - %4.2f, DG ", results[f].lambda_a,
                                              results[f].lambda_b);*/
        printf("lambda ");
        printf(lamformat, results[f].a->native_lambda);
        printf(" - ");
        printf(lamformat, results[f].b->native_lambda);
        printf(",   DG ");

        printf(dgformat,results[f].dg*kT);
        if (bEE)
        {
            printf(" +/- ");
            printf(dgformat,results[f].dg_err*kT);
        }
        if (histrange_err)
        {
            printf(" (max. range err. = ");
            printf(dgformat,results[f].dg_histrange_err*kT);
            printf(")");
            sum_histrange_err += results[f].dg_histrange_err*kT;
        }

        printf("\n");
        dg_tot += results[f].dg;
    }
    printf("\n");
    printf("total  ");
    printf(lamformat, results[0].a->native_lambda);
    printf(" - ");
    printf(lamformat, results[nresults-1].b->native_lambda);
    printf(",   DG ");

    printf(dgformat,dg_tot*kT);
    if (bEE)
    {
        stat_err=bar_err(nbmin,nbmax,partsum)*kT;
        printf(" +/- ");
        printf(dgformat, max(max(stat_err, sum_disc_err), sum_histrange_err));
    }
    printf("\n");
    if (disc_err)
    {
        printf("\nmaximum discretization error = ");
        printf(dgformat,sum_disc_err);
        if (bEE && stat_err < sum_disc_err)
        {
            printf("WARNING: discretization error (%g) is larger than statistical error.\n       Decrease histogram spacing for more accurate results\n", stat_err);
        }
    }
    if (histrange_err)
    {
        printf("\nmaximum histogram range error = ");
        printf(dgformat,sum_histrange_err);
        if (bEE && stat_err < sum_histrange_err)
        {
            printf("WARNING: histogram range error (%g) is larger than statistical error.\n       Increase histogram range for more accurate results\n", stat_err);
        }

    }
    printf("\n");


    if (fpi != NULL)
    {
        fprintf(fpi, xvg2format,
                results[nresults-1].b->native_lambda, dg_tot);
        ffclose(fpi);
    }

    do_view(oenv,opt2fn_null("-o",NFILE,fnm),"-xydy");
    do_view(oenv,opt2fn_null("-oi",NFILE,fnm),"-xydy");
    
    thanx(stderr);
    
    return 0;
}
