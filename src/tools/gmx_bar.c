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


typedef struct
{
    unsigned int *bin;          /* the histogram values */
    double dx;                  /* the histogram spacing */
    gmx_large_int_t x0;         /* the histogram start point as int */

    int nbin;                   /* the number of bins */
    gmx_large_int_t sum;        /* the total number of counts */

    double starttime, endtime;  /* start time, end time of histogram */
} barhist_t;


/* the raw data from a simulation */
typedef struct {
    char   *filename;
    int    ftp;     /* file type */
    int    nset;    /* number of lambdas, including dhdl */
    int *np;        /* number of data points (du or hists) per lambda */
    int  np_alloc;  /* number of points (du or hists) allocated */
    double temp;    /* temperature */
    double *lambda; /* the lambdas (of first index for y). The first one
                       is just the 'native' lambda */
    double *t;      /* the times (of second index for y) */
    double **y;     /* the dU values. y[0] holds the derivative, while
                       further ones contain the energy differences between
                       the native lambda and the 'foreign' lambdas. */
    barhist_t **hists; /* histograms. */
} barsim_t;

/* an aggregate of samples for partial free energy calculation */
typedef struct barsamples_t 
{
    double native_lambda;
    double foreign_lambda;
    double temp;

    int nndu; /* number of delta u sample collections */
    int *ndu; /* the number of delta U samples */
    double **du; /* the delta u's */
    double **t; /* the times associated with those samples */

    int nhist; /* the number of histograms */ 
    barhist_t *hists; /* the histograms */

    gmx_large_int_t ntot; /* total number of samples */

    struct barsamples_t  *next, *prev; /* the next and prev in the list */
} barsamples_t;

/* all the samples associated with a lambda point */
typedef struct barlambda_t
{
    double lambda; /* the native lambda */
    double temp; /* temperature */

    int nbs; /* number of barsamples_ts with their own lambda */
    barsamples_t *bs; /* the samples */

    barsamples_t bs_head; /* the pre-allocated list head for the linked list. 
                             Also used to contain the dhdl data. */

    struct barlambda_t  *next, *prev; /* the next and prev in the list */
} barlambda_t;


/* calculated values. */
typedef struct {
    /*barsim_t *a, *b; *//* the simulation data */
    barsamples_t *a, *b;

    /*double lambda_a, lambda_b; *//* the lambda values at a and b */

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



static void barsim_init(barsim_t *ba)
{
    ba->filename=NULL;
    ba->nset=0;
    ba->np_alloc=0;
    ba->np=NULL;
    ba->y=NULL;
    ba->hists=NULL;
}

static void barsamples_init(barsamples_t *bs, double native_lambda,
                            double foreign_lambda, double temp)
{
    bs->native_lambda=native_lambda;
    bs->foreign_lambda=foreign_lambda;
    bs->temp=temp;
    bs->nndu=0;
    bs->nhist=0;
    bs->ntot=0;
    bs->ndu=NULL;
    bs->du=NULL;
    bs->t=NULL;
    bs->hists=NULL;

    bs->next=NULL;
    bs->prev=NULL;
}

/* destroy the data structures directly associated with the structure, not
   the data it points to */
static void barsamples_destroy(barsamples_t *bs)
{
    sfree(bs->ndu);
    sfree(bs->du);
    sfree(bs->t);
    sfree(bs->hists);
}

static void barlambda_init(barlambda_t *bl, double native_lambda, double temp)
{
    bl->lambda=native_lambda;
    bl->temp=temp;
    bl->nbs=0;
    bl->next=NULL;
    bl->prev=NULL;
    bl->bs=&(bl->bs_head);
    barsamples_init(bl->bs, native_lambda, 0., 0.);
    bl->bs->next=bl->bs;
    bl->bs->prev=bl->bs;
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


static bool lambda_same(double lambda1, double lambda2)
{
    return gmx_within_tol(lambda1, lambda2, 10*GMX_REAL_EPS);
}

/* find the barsamples_t associated with a barlambda that corresponds to
   a specific foreign lambda */
barsamples_t *barlambda_find_barsample(barlambda_t *bl, double foreign_lambda)
{
    barsamples_t *bs=bl->bs->next;

    while(bs != bl->bs)
    {
        if (lambda_same(bs->foreign_lambda,foreign_lambda))
        {
            return bs;
        }
        bs=bs->next;
    }

    return NULL;
}

/* create subsample i out of ni from an existing barsample_t */
void barsamples_create_subsample(barsamples_t *bs, barsamples_t *bs_orig, 
                                int i, int ni)
{
    int j;
    int hist_start, hist_end;

    *bs = *bs_orig; /* just copy all fields */

    /* allocate proprietary memory */
    snew(bs->ndu, bs_orig->nndu);
    snew(bs->du, bs_orig->nndu);
    snew(bs->t, bs_orig->nndu);
    snew(bs->hists, bs_orig->nhist);

    /* fix all the du fields */
    for(j=0;j<bs_orig->nndu;j++)
    {
        /* these ugly casts avoid a nasty overflow if ndu is very big. */
        int start=(int)(bs_orig->ndu[j]*((double)(i)/ni));
        int end=(int)(bs_orig->ndu[j]*((double)(i+1.)/ni))-1;

        bs->ndu[j]=end-start+1;
        bs->du[j]=&(bs_orig->du[j][start]);
        bs->t[j]=&(bs_orig->t[j][start]);
    }
    /* and all histograms */
    hist_start = (int)(bs_orig->nhist*((double)(i)/ni));
    hist_end = (int)(bs_orig->nhist*((double)(i+1)/ni))-1;
    bs->nhist=(hist_end-hist_start)+1;
    for(j=0;j<bs->nhist;j++)
        bs->hists[j] = bs_orig->hists[j+hist_start];
 
    /* and count ntot */
    bs->ntot = 0;
    for(i=0;i<bs->nndu;i++)
        bs->ntot += bs->ndu[i];
    for(i=0;i<bs->nhist;i++)
        bs->ntot += bs->hists[i].sum;
}


/* add simulation data to a barlambda structure */
void barlambda_add_sim(barlambda_t *bl, barsim_t *ba, double begin, double end)
{
    int i,j;

    if (!lambda_same(bl->lambda, ba->lambda[0]))
        gmx_fatal(FARGS, "barlambda_add_sim lambdas inconsistent!");

    for(i=0;i<ba->nset;i++)
    {
        if (ba->np[i] > 0)
        {
            barsamples_t *bs=NULL;

            if (i>0) 
            {
                bs=barlambda_find_barsample(bl, ba->lambda[i]);
            }

            if (!bs)
            {
                /* we make a new barsamples_t */
                snew(bs, 1);
                barsamples_init(bs, bl->lambda, ba->lambda[i], bl->temp);
                /* insert it */
                bs->next=bl->bs;
                bs->prev=bl->bs->prev;
                bs->next->prev=bs;
                bs->prev->next=bs;
            }
            /* else 
               there already exists a barsamples_t with this foreign lambda 
               and we don't need to do anything */

            /* and add our samples */
            if (ba->y && ba->y[i])
            {
                int ndu=ba->np[i];
                double *y=ba->y[i];
                double *t=ba->t;

                /* get the right time */
                while( ndu>0 && *t < begin )
                {
                    ndu--;
                    y++;
                    t++;
                }
                if (end > begin)
                {
                    while(t[ndu-1] > end)
                    {
                        ndu--;
                    }
                }

                bs->nndu++;
                srenew(bs->ndu, bs->nndu);
                srenew(bs->du, bs->nndu);
                srenew(bs->t, bs->nndu);
                bs->ndu[bs->nndu-1]=ndu;
                bs->du[bs->nndu-1]=y;
                bs->t[bs->nndu-1]=t;
                bs->ntot += ndu;
            }
            if (ba->hists && ba->hists[i])
            {
                int nhist_prev = bs->nhist;
                int starti=0;
                int endi=ba->np[i];

                while( starti<endi && ba->hists[i][starti].endtime<begin)
                {
                    starti++;
                }
                if (end > begin)
                {
                    while((endi>starti) && (ba->hists[i][endi].starttime>end))
                    {
                        endi--;
                    }
                }
                bs->nhist += endi-starti;
                srenew(bs->hists, bs->nhist);
                for(j=starti;j<endi;j++)
                {
                    int hi=nhist_prev+(j-starti);
                    bs->hists[hi] = ba->hists[i][j];
                    bs->ntot += bs->hists[hi].sum;
                }
            }
        }
    }
}


/* assemble an ordered list of barlambda_ts */
barlambda_t *barlambdas_list_create(barsim_t *ba, int nfile, 
                                    int *nbl, double begin, double end)
{
    barsim_t ba_tmp;
    barlambda_t *bl_head; /* the head of the list */
    int i;

    snew(bl_head, 1); /* the first element is a dummy element */
    barlambda_init(bl_head, 0, 0);
    bl_head->next=bl_head;
    bl_head->prev=bl_head;
    *nbl=0;

    for(i=0; i<nfile; i++)
    {
        barlambda_t *bl=bl_head->next;
        bool inserted=FALSE;

        while(bl != bl_head)
        {
            if (lambda_same(ba[i].lambda[0], bl->lambda))
            {
                /* this lambda is the same as a previous lambda; add 
                   our samples */
                barlambda_add_sim(bl, &(ba[i]), begin, end);
                inserted=TRUE;
                break;
            }
            if ( bl->lambda > ba[i].lambda[0] )
            {
                break;
            }
            bl=bl->next;
        }

        if (!inserted)
        {
            barlambda_t *bl_new;
            
            snew(bl_new, 1); 
            (*nbl)++;
            barlambda_init(bl_new,ba[i].lambda[0], ba[i].temp);
            /* update linked list */
            bl_new->next=bl;
            bl_new->prev=bl->prev;
            bl_new->prev->next=bl_new;
            bl_new->next->prev=bl_new;
            barlambda_add_sim(bl_new, &(ba[i]), begin, end);
        }
    }
    return bl_head;
}

/* make a histogram out of a sample collection */
void barsamples_make_hist(barsamples_t *bs, int **bin, 
                          int *nbin_alloc, int *nbin, 
                          double *dx, double *xmin, int nbin_default)
{
    int i,j;
    bool dx_set=FALSE;
    bool xmin_set=FALSE;

    bool xmax_set=FALSE;
    bool xmax_set_hard=FALSE; /* whether the xmax is bounded by the limits of 
                                  a histogram */
    double xmax=-1;

    /* first determine dx and xmin; try the histograms */
    for(i=0;i<bs->nhist;i++)
    {
        double xmax_now=(bs->hists[i].x0+bs->hists[i].nbin)*bs->hists[i].dx;

        /* we use the biggest dx*/
        if ( (!dx_set) || bs->hists[i].dx > *dx)
        {
            dx_set=TRUE;
            *dx = bs->hists[i].dx;
        }
        if ( (!xmin_set) || (bs->hists[i].x0*bs->hists[i].dx) < *xmin)
        {
            xmin_set=TRUE;
            *xmin = (bs->hists[i].x0*bs->hists[i].dx);
        }
        
        if ( (!xmax_set) || (xmax_now>xmax && !xmax_set_hard) )
        {
            xmax_set=TRUE;
            xmax = xmax_now;
            if (bs->hists[i].bin[bs->hists[i].nbin-1] != 0)
                xmax_set_hard=TRUE;
        }
        if ( bs->hists[i].bin[bs->hists[i].nbin-1]!=0 && (xmax_now < xmax) )
        {
            xmax_set_hard=TRUE;
            xmax = xmax_now;
        }
    }
    /* and the delta us */
    for(i=0;i<bs->nndu;i++)
    {
        if (bs->ndu[i]>0)
        {
            /* determine min and max */
            double du_xmin=bs->du[i][0]; 
            double du_xmax=bs->du[i][0];
            for(j=1;j<bs->ndu[i];j++)
            {
                if (bs->du[i][j] < du_xmin)
                    du_xmin = bs->du[i][j];
                if (bs->du[i][j] > du_xmax)
                    du_xmax = bs->du[i][j];
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

    /* now add the acutal data */   
    for(i=0;i<bs->nhist;i++)
    {
        double xmin_hist=bs->hists[i].x0*bs->hists[i].dx;
        for(j=0;j<bs->hists[i].nbin;j++)
        {
            /* calculate the bin corresponding to the middle of the original
               bin */
            double x=bs->hists[i].dx*(j+0.5) + xmin_hist;
            int binnr=(int)((x-(*xmin))/(*dx));

            if (binnr >= *nbin || binnr<0)
                binnr = (*nbin)-1;

            (*bin)[binnr] += bs->hists[i].bin[j]; 
        }
    }
    for(i=0;i<bs->nndu;i++)
    {
        for(j=0;j<bs->ndu[i];j++)
        {
            int binnr=(int)((bs->du[i][j]-(*xmin))/(*dx));
            if (binnr >= *nbin || binnr<0)
                binnr = (*nbin)-1;

            (*bin)[binnr] ++;
        }
    }
}

/* write a collection of histograms to a file */
void barlambdas_histogram(barlambda_t *bl_head, const char *filename, 
                          int nbin_default, const output_env_t oenv)
{
    char label_x[STRLEN];
    const char *title="N(\\Delta H)";
    const char *label_y="Samples";
    FILE *fp;
    barlambda_t *bl;
    int nsets=0;
    char **setnames=NULL;
    bool first_set=FALSE;
    /* histogram data: */
    int *hist=NULL;
    int nbin=0;
    int nbin_alloc=0;
    double dx=0;
    double min;
    int i;

    sprintf(label_x, "\\Delta H (%s)", unit_energy);

    fp=xvgropen_type(filename, title, label_x, label_y, exvggtXNY, oenv);

    /* first get all the set names */
    bl=bl_head->next;
    /* iterate over all lambdas */
    while(bl!=bl_head)
    {
        barsamples_t *bs=bl->bs->next;

        /* iterate over all samples */
        while(bs!=bl->bs)
        {
            nsets++;
            srenew(setnames, nsets); 
            snew(setnames[nsets-1], STRLEN);
            if (!lambda_same(bs->foreign_lambda, bs->native_lambda))
            {
                sprintf(setnames[nsets-1], 
                        "N(H(\\lambda=%g) - H(\\lambda=%g) | \\lambda=%g)", 
                        bs->foreign_lambda, bs->native_lambda, 
                        bs->native_lambda);
            }
            else
            {
                sprintf(setnames[nsets-1], 
                        "N(dH/d\\lambda | \\lambda=%g)", 
                        bs->native_lambda);
            }
            bs=bs->next;
        }

        bl=bl->next;
    }
    xvgr_legend(fp, nsets, (const char**)setnames, oenv);


    /* now make the histograms */
    bl=bl_head->next;
    /* iterate over all lambdas */
    while(bl!=bl_head)
    {
        barsamples_t *bs=bl->bs->next;

        /* iterate over all samples */
        while(bs!=bl->bs)
        {
            if (!first_set)
                xvgr_new_dataset(fp, oenv);
    
            barsamples_make_hist(bs, &hist, &nbin_alloc, &nbin, &dx, &min,
                                 nbin_default);

            for(i=0;i<nbin;i++)
            {
                double xmin=i*dx;
                double xmax=(i+1)*dx;

                fprintf(fp, "%g %d\n%g %d\n", xmin, hist[i], xmax, hist[i]);
            }

            first_set=FALSE;
            bs=bs->next;
        }

        bl=bl->next;
    }

    if(hist)
        sfree(hist);    

    xvgrclose(fp); 
}

/* create a collection (array) of barres_t object given a ordered linked list 
   of barlamda_t sample collections */
static barres_t *barres_list_create(barlambda_t *bl_head, int *nres)
{
    barlambda_t *bl;
    int nlambda=0;
    barres_t *res;
    int i;
    bool dhdl=FALSE;
    bool first=TRUE;

    /* first count the barlambdas */
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
        barsamples_t *bs,*bsprev;
        barres_t *br=&(res[*nres]);
        /* there is always a previous one. we search for that as a foreign 
           lambda: */
        bsprev=barlambda_find_barsample(bl->prev, bl->lambda);
        bs=barlambda_find_barsample(bl, bl->prev->lambda);

        barres_init(br);

        if (!bsprev && !bs)
        {
            /* we use dhdl */

            bsprev=barlambda_find_barsample(bl->prev, bl->prev->lambda);
            bs=barlambda_find_barsample(bl, bl->lambda);

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
        if (!bsprev)
        {
            gmx_fatal(FARGS,"Could not find a set for lambda = %g in the files for lambda = %g",bl->lambda,bl->prev->lambda);
        }
        if (!bs)
        {
            gmx_fatal(FARGS,"Could not find a set for lambda = %g in the files for lambda = %g",bl->prev->lambda,bl->lambda);
        }
        br->a = bsprev;
        br->b = bs;

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

    for(i=0;i<nres;i++)
    {
        barres_t *br=&(res[i]);

        for(j=0;j<br->a->nhist;j++)
        {
            disc_err=max(disc_err, br->a->hists[j].dx);
        }
        for(j=0;j<br->b->nhist;j++)
        {
            disc_err=max(disc_err, br->b->hists[j].dx);
        }
    } 
    return disc_err;
}

static double calc_bar_sum(int n,double *W,double Wfac,double sbMmDG)
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
static double calc_bar_sum_hist(barhist_t *hist, double Wfac, double sbMmDG,
                                int type)
{
    double sum=0.;
    int i;
    int max=hist->nbin-1;
    /* normalization factor multiplied with bin width and
       number of samples (we normalize through M): */
    double normdx = 1.;

    if (type==1) 
    {
        max=hist->nbin; /* we also add whatever was out of range */
    }

    for(i=0;i<max;i++)
    {
        double x=Wfac*((i+hist->x0)+0.5)*hist->dx; /* bin middle */
        double pxdx=hist->bin[i]*normdx; /* p(x)dx */
    
        sum += pxdx/(1. + exp(x + sbMmDG));
    }

    return sum;
}

static double calc_bar_lowlevel(barsamples_t *ba, barsamples_t *bb,
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
    n1 = ba->ntot;
    n2 = bb->ntot;

    M = log(n1/n2);

    if (!lambda_same(ba->native_lambda, ba->foreign_lambda))
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
        double delta_lambda=bb->native_lambda-ba->native_lambda;
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
    Wmin=FLT_MAX;
    Wmax=-FLT_MAX;
    for(i=0;i<ba->nndu;i++)
    {
        for(j=1; j<ba->ndu[i]; j++)
        {
            Wmin = min(Wmin,ba->du[i][j]*Wfac1);
            Wmax = max(Wmax,ba->du[i][j]*Wfac1);
        }
    }
    for(i=0;i<ba->nhist;i++)
    {
        Wmin = min(Wmin,ba->hists[i].x0*ba->hists[i].dx);
        for(j=ba->hists[i].nbin-1;j>=0;j--)
        {
            /* look for the highest value bin with values */
            if (ba->hists[i].bin[j]>0)
            {
                Wmax=max(Wmax,(j+ba->hists[i].x0+1)*ba->hists[i].dx);
                break;
            }
        }
    }
    for(i=0;i<bb->nndu;i++)
    {
        for(j=1; j<bb->ndu[i]; j++)
        {
            Wmin = min(Wmin,bb->du[i][j]*Wfac1);
            Wmax = max(Wmax,bb->du[i][j]*Wfac1);
        }
    }
    for(i=0;i<bb->nhist;i++)
    {
        Wmin = min(Wmin,bb->hists[i].x0*bb->hists[i].dx);
        for(j=bb->hists[i].nbin-1;j>=0;j--)
        {
            /* look for the highest value bin with values */
            if (bb->hists[i].bin[j]>0)
            {
                Wmax=max(Wmax,(j+bb->hists[i].x0+1)*bb->hists[i].dx);
                break;
            }
        }
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
        for(i=0;i<ba->nndu;i++)
        { 
            /* first the du lists */
            dDG1 += calc_bar_sum(ba->ndu[i], ba->du[i], Wfac1, (M-DG1));
        }
        for(i=0;i<ba->nhist;i++)
        {
            /* then the histograms */
            dDG1 += calc_bar_sum_hist(&(ba->hists[i]), Wfac1, (M-DG1), type);
        }

        for(i=0;i<bb->nndu;i++)
        {
            dDG1 -= calc_bar_sum(bb->ndu[i], bb->du[i], Wfac2, -(M-DG1));
        }
        for(i=0;i<bb->nhist;i++)
        {
            dDG1 -= calc_bar_sum_hist(&(bb->hists[i]), Wfac2, -(M-DG1), type);
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

static void calc_rel_entropy(barsamples_t *ba, barsamples_t *bb,
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
    n1 = ba->ntot;
    n2 = bb->ntot;

    /* to ensure the work values are the same as during the delta_G */
    if (!lambda_same(ba->native_lambda, ba->foreign_lambda))
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
        double delta_lambda=bb->native_lambda-ba->native_lambda;
        Wfac1 =  beta*delta_lambda;
        Wfac2 = -beta*delta_lambda;
    }

    /* first calculate the average work in both directions */
    for(i=0;i<ba->nndu;i++)
    { 
        for(j=0;j<ba->ndu[i];j++)
        {
            W_ab += Wfac1*ba->du[i][j];
        }
    }
    for(i=0;i<ba->nhist;i++)
    {
        /* then the histograms */
        /* normalization factor multiplied with bin width and
           number of samples (we normalize through M): */
        double normdx = 1.;
        barhist_t *hist=&(ba->hists[i]);

        for(j=0;j<hist->nbin;j++)
        {
            double x=Wfac1*((j+hist->x0)+0.5)*hist->dx; /* bin middle */
            double pxdx=hist->bin[j]*normdx; /* p(x)dx */

            W_ab += pxdx*x;
        }
    }
    W_ab/=n1;
    for(i=0;i<bb->nndu;i++)
    { 
        for(j=0;j<bb->ndu[i];j++)
        {
            W_ba += Wfac2*bb->du[i][j];
        }
    }
    for(i=0;i<bb->nhist;i++)
    {
        /* then the histograms */
        /* normalization factor multiplied with bin width and
           number of samples (we normalize through M): */
        double normdx = 1.;
        barhist_t *hist=&(bb->hists[i]);

        for(j=0;j<hist->nbin;j++)
        {
            double x=Wfac2*((j+hist->x0)+0.5)*hist->dx; /* bin middle */
            double pxdx=hist->bin[j]*normdx; /* p(x)dx */

            W_ba += pxdx*x;
        }
    }
    W_ba/=n2;
   
    /* then calculate the relative entropies */
    *sa = (W_ab - dg);
    *sb = (W_ba + dg);
}

static void calc_dg_stddev(barsamples_t *ba, barsamples_t *bb,
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
    n1 = ba->ntot;
    n2 = bb->ntot;

    /* to ensure the work values are the same as during the delta_G */
    if (!lambda_same(ba->native_lambda, ba->foreign_lambda))
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
        double delta_lambda=bb->native_lambda-ba->native_lambda;
        Wfac1 =  beta*delta_lambda;
        Wfac2 = -beta*delta_lambda;
    }

    M = log(n1/n2);

    /* calculate average in both directions */
    for(i=0;i<ba->nndu;i++)
    { 
        for(j=0;j<ba->ndu[i];j++)
        {
            sigmafact += 1./(2. + 2.*cosh((M + Wfac1*ba->du[i][j] - dg)));
        }
    }
    for(i=0;i<ba->nhist;i++)
    {
        /* then the histograms */
        /* normalization factor multiplied with bin width and
           number of samples (we normalize through M): */
        double normdx = 1.;
        barhist_t *hist=&(ba->hists[i]);

        for(j=0;j<hist->nbin;j++)
        {
            double x=Wfac1*((j+hist->x0)+0.5)*hist->dx; /* bin middle */
            double pxdx=hist->bin[j]*normdx; /* p(x)dx */
            
            sigmafact += pxdx/(2. + 2.*cosh((M + x - dg)));
        }
    }
    for(i=0;i<bb->nndu;i++)
    { 
        for(j=0;j<bb->ndu[i];j++)
        {
            sigmafact += 1./(2. + 2.*cosh((M - Wfac2*bb->du[i][j] - dg)));
        }
    }
    for(i=0;i<bb->nhist;i++)
    {
        /* then the histograms */
        /* normalization factor multiplied with bin width and
           number of samples (we normalize through M): */
        double normdx = 1.;
        barhist_t *hist=&(bb->hists[i]);

        for(j=0;j<hist->nbin;j++)
        {
            double x=Wfac2*((j+hist->x0)+0.5)*hist->dx; /* bin middle */
            double pxdx=hist->bin[j]*normdx; /* p(x)dx */

            sigmafact += pxdx/(2. + 2.*cosh((M - x - dg)));
        }
    }
    sigmafact /= (n1 + n2);
 
  
    /* Eq. 10 from 
       Shirts, Bair, Hooker & Pande, Phys. Rev. Lett 91, 140601 (2003): */
    *stddev = sqrt(((1./sigmafact) - ( (n1+n2)/n1 + (n1+n2)/n2 )));
}



static void calc_bar(barres_t *br, double tol, 
                     int npee_min, int npee_max, bool *bEE, 
                     double *partsum)
{
    int npee,p;
    double dg_sig2,sa_sig2,sb_sig2,stddev_sig2; /* intermediate variance values
                                                   for calculated quantities */
    int nsample1, nsample2;
    double temp=br->a->temp;
    int i,j;
    double dg_min, dg_max;

    br->dg=calc_bar_lowlevel(br->a, br->b, temp, tol, 0);

    br->dg_disc_err = 0.;
    br->dg_histrange_err = 0.;
    if ((br->a->nhist > 0) || (br->b->nhist > 0) )
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
        for(i=0;i<br->a->nhist;i++)
        {
            br->dg_disc_err=max(br->dg_disc_err, br->a->hists[i].dx);
        }
        for(i=0;i<br->b->nhist;i++)
        {
            br->dg_disc_err=max(br->dg_disc_err, br->b->hists[i].dx);
        }
    }

    calc_rel_entropy(br->a, br->b, temp, br->dg, &(br->sa), &(br->sb));
                     
    calc_dg_stddev(br->a, br->b, temp, br->dg, &(br->dg_stddev) );

    dg_sig2 = 0;
    sa_sig2 = 0;
    sb_sig2 = 0;
    stddev_sig2 = 0;

    /* we look for the smallest sample size */
    nsample1=INT_MAX;
    for(i=0;i<br->a->nndu;i++)
        nsample1 = min( nsample1, br->a->ndu[i] );
    if (br->a->nhist > 0)
        nsample1 = min( nsample1, br->a->nhist );

    nsample2=INT_MAX;
    for(i=0;i<br->b->nndu;i++)
        nsample2 = min( nsample2, br->b->ndu[i] );
    if (br->b->nhist > 0)
        nsample2 = min( nsample2, br->b->nhist );


    printf("nsample1=%d, nsample2=%d\n", nsample1, nsample2);
    if (nsample1 >= npee_max && nsample2 >= npee_max)
    {
        barsamples_t ba, bb;
        if ( (2*npee_max > nsample1) || (2*npee_max > nsample2) )
        {
            if (npee_min != min(nsample1, nsample2) && 
                npee_max != min(nsample1, nsample2) )
            {

                npee_min = npee_max = min(nsample1, nsample2);
                printf("NOTE: redefining nbin and nbmax to %d at lambda=%g-%g\n     because of the small number of blocks available\n", 
                       npee_min, br->a->native_lambda, br->b->native_lambda);
            }
        }

        barsamples_init(&ba, br->a->native_lambda, br->a->foreign_lambda, 
                        br->a->temp);
        barsamples_init(&bb, br->b->native_lambda, br->b->foreign_lambda, 
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

                barsamples_create_subsample(&ba, br->a, p, npee);
                barsamples_create_subsample(&bb, br->b, p, npee);

                dgp=calc_bar_lowlevel(&ba, &bb, temp, tol, 0);
                dgs  += dgp;
                dgs2 += dgp*dgp;

                partsum[npee*(npee_max+1)+p] += dgp;

                calc_rel_entropy(&ba, &bb, temp, dgp, &sac, &sbc); 
                dsa  += sac;
                dsa2 += sac*sac;
                dsb  += sbc;
                dsb2 += sbc*sbc;
                calc_dg_stddev(&ba, &bb, temp, dgp, &stddevc );

                dstddev  += stddevc;
                dstddev2 += stddevc*stddevc;

                barsamples_destroy(&ba);
                barsamples_destroy(&bb);
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
    else
    {
        *bEE = FALSE;
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


static double legend2lambda(char *fn,const char *legend,bool bdhdl)
{
    double lambda=0;
    const char   *ptr;

    if (legend == NULL)
    {
        gmx_fatal(FARGS,"There is no legend in file '%s', can not deduce lambda",fn);
    }
    ptr = strrchr(legend,' ');
    if (( bdhdl &&  strstr(legend,"dH") == NULL) ||
        (!bdhdl && (strchr(legend,'D') == NULL ||
                    strchr(legend,'H') == NULL)) ||
        ptr == NULL)
    {
        gmx_fatal(FARGS,"There is no proper lambda legend in file '%s', can not deduce lambda",fn);
    }
    if (sscanf(ptr,"%lf",&lambda) != 1)
    {
        gmx_fatal(FARGS,"There is no proper lambda legend in file '%s', can not deduce lambda",fn);
    }

    return lambda;
}

static bool subtitle2lambda(const char *subtitle,double *lambda)
{
    bool bFound;
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

static void read_barsim_xvg(char *fn,double begin,double end,real temp,
                            barsim_t *ba)
{
    int  i;
    char *subtitle,**legend,*ptr;
    int np;

    barsim_init(ba);

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
        if (temp <= 0)
        {
            gmx_fatal(FARGS,"Did not find a temperature in the subtitle in file '%s', use the -temp option of g_bar",fn);
        }
        ba->temp = temp;
    }

    snew(ba->lambda,ba->nset-1);
    if (legend == NULL)
    {
        /* Check if we have a single set, nset=2 means t and dH/dl */
        if (ba->nset == 2)
        {
            /* Try to deduce lambda from the subtitle */
            if (subtitle != NULL &&
                !subtitle2lambda(subtitle,&ba->lambda[0]))
            {
                /* Deduce lambda from the file name */
                ba->lambda[0] = filename2lambda(fn);
            }
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
            /* Read lambda from the legend */
            ba->lambda[i] = legend2lambda(fn,legend[i],i==0);
        }
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

static void read_edr_rawdh(barsim_t *ba, t_enxblock *blk, int id, 
                           bool lambda_set, double starttime, 
                           double endtime)
{
    int j;
    bool allocated;

    if (starttime < 0 || endtime < 0)
    {
        gmx_file("start time or end time <0. Probably due to wrong block order\n");
    }


    /* check the block types etc. */
    if ( (blk->nsub < 2) ||
         (blk->sub[0].type != xdr_datatype_double) ||
         (
          (blk->sub[1].type != xdr_datatype_float) &&
          (blk->sub[1].type != xdr_datatype_double) 
         ) ||
         (blk->sub[0].nr < 1) )
    {
        gmx_fatal(FARGS, "Unexpected block data in file %s", ba->filename);
    }

    /* we assume the blocks come in the same order */
    if (!lambda_set)
    {
        ba->lambda[id] = blk->sub[0].dval[0];
    }
    else
    {
        if (fabs(ba->lambda[id] - blk->sub[0].dval[0])>1e-12)
        {
            gmx_fatal(FARGS, "lambdas change in %s", ba->filename);
        }
    }

    /* now make room for the data */
    if (ba->np[id] + blk->sub[1].nr > ba->np_alloc )
    {
        ba->np_alloc = ba->np[id] + blk->sub[1].nr + 2;
        srenew(ba->t, ba->np_alloc);
        for(j=0;j<ba->nset;j++)
        {
            srenew(ba->y[j], ba->np_alloc);
        }
        allocated=TRUE;
    }
    /* and copy the data */
    for(j=0;j<blk->sub[1].nr;j++)
    {
        unsigned int index=ba->np[id];
        ba->np[id]++;
        if (allocated)
            ba->t[index] = ((endtime-starttime)*j)/blk->sub[1].nr;
        if (blk->sub[1].type == xdr_datatype_float)
        {
            ba->y[id][index] = blk->sub[1].fval[j];
        }
        else
        {
            ba->y[id][index] = blk->sub[1].dval[j];
        }
    }
}

static void read_edr_hist(barsim_t *ba, t_enxblock *blk, int id,
                          bool lambda_set, 
                          double starttime, double endtime)
{
    int j;
    barhist_t *bh;

    if (starttime < 0 || endtime < 0)
    {
        gmx_file("start time or end time <0. Probably due to wrong block order\n");
    }

    /* check the block types etc. */
    if ( (blk->nsub < 3) ||
         (blk->sub[0].type != xdr_datatype_double) ||
         (blk->sub[1].type != xdr_datatype_large_int) ||
         (blk->sub[2].type != xdr_datatype_int) ||
         (blk->sub[0].nr < 2)  ||
         (blk->sub[1].nr < 1) )
    {
        gmx_fatal(FARGS, "Unexpected block data in file %s", ba->filename);
    }

    if (ba->y != NULL)
    {
        gmx_fatal(FARGS, "Can't have both histograms and raw delta U data in one file %s", ba->filename);
    }

    if (!lambda_set)
    {
        ba->lambda[id] = blk->sub[0].dval[0];
    }
    else
    {
        if (fabs(ba->lambda[id] - blk->sub[0].dval[0])>1e-12)
        {
            gmx_fatal(FARGS, "lambdas change in %s: %g, %g", ba->filename,
                      ba->lambda[id], blk->sub[0].dval[0]);
        }
    }
    if (blk->sub[2].nr > 0)
    {
        /* make room for the data */
        if (ba->np[id] + 1 > ba->np_alloc)
        {
            ba->np_alloc = ba->np[id] + 2;
            /*srenew(ba->t, ba->np_alloc);*/
            for(j=0;j<ba->nset;j++)
            {
                srenew(ba->hists[j], ba->np_alloc);
            }
        }

        bh=&(ba->hists[id][ba->np[id]]);

        bh->dx = blk->sub[0].dval[1];
        bh->starttime = starttime;
        bh->endtime = endtime;
        bh->x0 = blk->sub[1].lval[0];

        bh->nbin = blk->sub[2].nr;
        bh->sum = 0;
        snew(bh->bin, bh->nbin);

        for(j=0;j<bh->nbin;j++)
        {
            bh->bin[j] = (int)blk->sub[2].ival[j];
            bh->sum += bh->bin[j];
        }

        ba->np[id]++;
    }
}



static void read_barsim_edr(char *fn,double begin,double end,real temp,
                            barsim_t *ba)
{
    int i;
    ener_file_t fp;
    t_enxframe *fr; 
    int nblocks=0;
    bool lambda_set=FALSE;
    int nre;
    gmx_enxnm_t *enm=NULL;


    barsim_init(ba);

    fp = open_enx(fn,"r");
    do_enxnms(fp,&nre,&enm);
    snew(fr, 1);

    ba->filename = fn;

    while(do_enx(fp, fr))
    {
        /* count the data blocks */
        int nblocks_raw=0;
        int nblocks_hist=0;
        int nlam=0;
        int nb=0;
        double starttime=-1;
        double endtime=-1;

        for(i=0;i<fr->nblock;i++)
        {
            if (fr->block[i].id == enxDHHIST )
                nblocks_hist++;
            if ( fr->block[i].id == enxDH )
                nblocks_raw++;
            if (fr->block[i].id == enxDHCOLL )
                nlam++;
        }

        if (nblocks_raw > 0 && nblocks_hist > 0 )
        {
            gmx_fatal(FARGS, "Can't handle both raw delta U data and histograms in the same file %s", fn);
        }
        if ((nblocks > 0 && (nblocks_raw+nblocks_hist)!=nblocks) || (nlam!=1 ))
        {
            gmx_fatal(FARGS, "Unexpected block count in %s: was %d, now %d\n",
                      fn, nblocks, nblocks_raw+nblocks_hist);
        }

        nblocks=nblocks_raw + nblocks_hist;
        ba->nset=nblocks+1;

        if (!ba->lambda) 
            snew(ba->lambda, ba->nset);
        if (!ba->np) 
        {
            snew(ba->np, ba->nset);
            for(i=0;i<ba->nset;i++)
                ba->np[i]=0.;
        }
        if (!ba->y && nblocks_raw>0) 
        {
            snew(ba->y, ba->nset);
            for(i=0;i<ba->nset;i++)
                ba->y[i]=NULL;
        }
        if (!ba->hists && nblocks_hist>0) 
        {
            snew(ba->hists, ba->nset);
            for(i=0;i<ba->nset;i++)
                ba->hists[i]=NULL;
        }
 
        
        for(i=0;i<fr->nblock;i++)
        {
            /* try to find the enxDHCOLL block */
            if (fr->block[i].id == enxDHCOLL)
            {
                if ( (fr->block[i].nsub < 1) || 
                     (fr->block[i].sub[0].type != xdr_datatype_double) ||
                     (fr->block[i].sub[0].nr < 4))
                {
                    gmx_fatal(FARGS, "Unexpected block data in file %s", fn);
                }

                ba->temp =      fr->block[i].sub[0].dval[0];
                ba->lambda[0] = fr->block[i].sub[0].dval[1];
                starttime =     fr->block[i].sub[0].dval[2];
                endtime =       fr->block[i].sub[0].dval[3];
            }

            if (fr->block[i].id == enxDH)
            {

                read_edr_rawdh(ba, &(fr->block[i]), nb+1, lambda_set, 
                               starttime, endtime);
                nb++;
            }

            if (fr->block[i].id == enxDHHIST)
            {
                read_edr_hist(ba, &(fr->block[i]), nb+1, lambda_set, 
                              starttime, endtime);
                nb++;
               
            }
        }
        lambda_set=TRUE;
    }
    fprintf(stderr, "\n");
}


int gmx_bar(int argc,char *argv[])
{
    static const char *desc[] = {
        "g_bar calculates free energy difference estimates through ",
        "Bennett's acceptance ratio method. ",
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
        "and the foreign lambda's from the legends of Delta H.[PAR]",

        "The lambda of the simulation is parsed from dhdl.xvg file's legend ",
        "containing the string 'dH', the foreign lambda's from the legend ",
        "containing the capitalized letters 'D' and 'H'. The temperature ",
        "is parsed from the legend line containing 'T ='.[PAR]",

        "The free energy estimates are determined using BAR with bisection, ",
        "the precision of the output is set with [TT]-prec[tt]. ",
        "An error estimate taking into account time correlations ",
        "is made by splitting the data into blocks and determining ",
        "the free energy differences over those blocks and assuming ",
        "the blocks are independent. ",
        "The final error estimate is determined from the average variance ",
        "over 5 blocks. A range of blocks numbers for error estimation can ",
        "be provided with the options [TT]-nbmin[tt] and [TT]-nbmax[tt].[PAR]",

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

    };
    static real begin=0,end=-1,temp=-1;
    int nd=2,nbmin=5,nbmax=5;
    int nbin=100;
    bool calc_s,calc_v;
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
        { efXVG, "-o",  "bar",    ffOPTWR },
        { efXVG, "-oi", "barint", ffOPTWR }, 
        { efXVG, "-oh", "histogram", ffOPTWR }, 
        { efEDR, "-g",  "energy", ffOPTRDMULT }
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
    barsim_t *ba;       /* the raw input data */
    barlambda_t *bl;    /* the pre-processed lambda data (linked list head) */
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
    bool     result_OK=TRUE,bEE=TRUE;

    bool     disc_err=FALSE;
    double   sum_disc_err=0.; /* discretization error */
    bool     histrange_err=FALSE;
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

    snew(ba,nfile_tot);
    snew(results,nfile_tot-1);
    snew(partsum,(nbmax+1)*(nbmax+1));
    nf = 0;

    /* read in all files. First xvg files */
    for(f=0; f<nxvgfile; f++)
    {
        read_barsim_xvg(fxvgnms[f],begin,end,temp,&ba[nf]);
        if (f > 0 && ba[nf].temp != ba[0].temp)
        {
            printf("\nWARNING: temperature for file '%s' (%g) is not equal to that of file '%s' (%g)\n\n",fxvgnms[nf],ba[nf].temp,fxvgnms[0],ba[0].temp);
        }

        if (ba[nf].nset == 0)
        {
            gmx_fatal(FARGS,"File '%s' contains fewer than two columns",
                      fxvgnms[nf]);
        }
        nf++;
    }
    /* then .edr files */
    for(f=0; f<nedrfile; f++)
    {
        read_barsim_edr(fedrnms[f],begin,end,temp,&ba[nf]);
        if (nf > 0 && ba[nf].temp != ba[0].temp)
        {
            printf("\nWARNING: temperature for file '%s' (%g) is not equal to that of file '%s' (%g)\n\n",fedrnms[nf],ba[nf].temp,fedrnms[0],ba[0].temp);
        }

        if (ba[nf].nset == 0)
        {
            gmx_fatal(FARGS,"File '%s' contains fewer than two columns",
                      fedrnms[nf]);
        }
        nf++;
    }

    /* print input file data summary */
    for(i=0;i<nf;i++)
    {
        int np;
        double begint, endt;
        if (ba[i].y)
        {
            if (ba[i].nset>0)
                np=ba[i].np[1];
            else
                np=ba[i].np[0];

            begint=ba[i].t[0];
            endt=ba[i].t[np-1];
        }
        else
        {
            np=ba[i].np[1];
            begint=ba[i].hists[1][0].starttime;
            endt=ba[i].hists[1][np-1].endtime;
        }
        printf("\n%s: %.1f - %.1f; lambdas:",ba[i].filename, begint, endt);

        for(j=0;j<ba[i].nset;j++)
        {
            printf(" %.3f", ba[i].lambda[j]);
            if (ba[i].np[j]>0)
            {
                if (ba[i].y)
                    printf(" (%d pts)",  ba[i].np[j]);
                if (ba[i].hists) 
                    printf(" (%d hists)", ba[i].np[j]);
            }
        }
    }
    printf("\n");

    /* Sort the data sets on lambda */
    bl=barlambdas_list_create(ba, nfile_tot, &nbs, begin, end);

    if (opt2bSet("-oh",NFILE,fnm))
    {
        barlambdas_histogram(bl, opt2fn("-oh",NFILE,fnm), nbin, oenv);
    }
   
    /* assemble the output structures from the barlambdas */
    results=barres_list_create(bl, &nresults);

    sum_disc_err=barres_list_max_disc_err(results, nresults);
    if (sum_disc_err > prec)
    {
        prec=sum_disc_err;
        nd = ceil(-log10(prec));
        printf("WARNING: setting the precision to %g because that is the minimum\n         reasonable number, given the expected discretization error.\n", prec);
    }
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
    kT   = BOLTZ*ba[0].temp;
    beta = 1/kT;

    printf("\nTemperature: %g K\n", ba[0].temp);

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
            fprintf(fpi, xvg2format, ba[f].lambda[0], dg_tot);
        }


        if (fpb != NULL)
        {
            fprintf(fpb, xvg3format,
                    0.5*(ba[f].lambda[0] + ba[f+1].lambda[0]),
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
    printf(lamformat, ba[0].lambda[0]);
    printf(" - ");
    printf(lamformat, ba[nfile_tot-1].lambda[0]);
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
                ba[nfile_tot-1].lambda[0], dg_tot);
        ffclose(fpi);
    }

    do_view(oenv,opt2fn_null("-o",NFILE,fnm),"-xydy");
    do_view(oenv,opt2fn_null("-oi",NFILE,fnm),"-xydy");
    
    thanx(stderr);
    
    return 0;
}
