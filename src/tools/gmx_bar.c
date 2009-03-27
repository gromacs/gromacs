/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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

#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "futil.h"
#include "statutil.h"
#include "copyrite.h"
#include "macros.h"
#include "physics.h"
#include "gmx_fatal.h"
#include "xvgr.h"

typedef struct {
    int    nset;
    int    np;
    double *lambda;
    double *t;
    double **y;
} barsim_t;

static double calc_bar_sum(int n,double *W,double MmDG,double sb)
{
    int    i;
    double sum;
    
    sum = 0;
    
    for(i=0; i<n; i++)
    {
        sum += 1/(1 + exp(sb*(W[i] + MmDG)));
    }
    
    return sum;
}

static double calc_bar_lowlevel(int n1,double *W1,int n2,double *W2,
                                real temp,double tol)
{
    double kT,beta,M;
    double DG;
    int    i;
    double Wmin,Wmax;
    double DG0,DG1,DG2,dDG1;
    double sum1,sum2;
    
    kT   = BOLTZ*temp;
    beta = 1/kT;
    
    M = kT*log((double)n1/(double)n2);
    
    Wmin = W1[0];
    Wmax = W1[0];
    for(i=0; i<n1; i++)
    {
        Wmin = min(Wmin,W1[i]);
        Wmax = max(Wmax,W1[i]);
    }
    for(i=0; i<n2; i++)
    {
        Wmin = min(Wmin,-W2[i]);
        Wmax = max(Wmax,-W2[i]);
    }
    DG0 = Wmin;
    DG2 = Wmax;
    
    /* For the comparison we can use twice the tolerance */
    if (debug)
    {
        fprintf(debug,"DG %9.5f %9.5f\n",DG0,DG2);
    }
    while (DG2 - DG0 > 2*tol*0.5*(fabs(DG0) + fabs(DG2)))
    {
        DG1 = 0.5*(DG0 + DG2);
        dDG1 =
            calc_bar_sum(n1,W1,  M-DG1 ,beta) -
            calc_bar_sum(n2,W2,-(M-DG1),beta);
        
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

static void get_begin_end(barsim_t *ba,real begin,real end,int *b,int *e)
{
    int i;
    
    i = 0;
    while (i + 1 < ba->np && ba->t[i] < begin)
    {
        i++;
    }
    if (i >= ba->np)
    {
        gmx_fatal(FARGS,"Some data send end before the start time %g",begin);
    }
    *b = i;

    i = ba->np;
    if (end >= begin)
    {
        while (i > *b && ba->t[i-1] > end)
        {
            i--;
        }
    }
    *e = i;
}

static int get_lam_set(barsim_t *ba,double lambda)
{
    int i;

    i = 1;
    while (i + 1 < ba->nset && ba->lambda[i] != lambda)
    {
        i++;
    }
    if (i == ba->nset)
    {
        gmx_fatal(FARGS,"Could not find a set for lambda = %g in the file of lambda = %g",lambda,ba->lambda[0]);
    }

    return i;
}

static void calc_bar(barsim_t *ba1,barsim_t *ba2,
                     real begin,real end,real temp,double tol,int npee,
                     double *dg,double *sig)
{
    int b1,b2,e1,e2,s1,s2,p;
    double dgp,dgs,dgs2;

    get_begin_end(ba1,begin,end,&b1,&e1);
    get_begin_end(ba2,begin,end,&b2,&e2);

    s1 = get_lam_set(ba1,ba2->lambda[0]);
    s2 = get_lam_set(ba2,ba1->lambda[0]);

    *dg = calc_bar_lowlevel(e1-b1,ba1->y[s1]+b1,
                            e2-b2,ba2->y[s2]+b2,
                            temp,tol);

    dgs  = 0;
    dgs2 = 0;
    if (e1 - b1 >= npee && e2 - b2 >= npee)
    {
        for(p=0; p<npee; p++)
        {
            dgp = calc_bar_lowlevel((e1-b1)/npee,ba1->y[s1]+b1+p*(e1-b1)/npee,
                                    (e2-b2)/npee,ba2->y[s2]+b2+p*(e2-b2)/npee,
                                    temp,tol);
            dgs  += dgp;
            dgs2 += dgp*dgp;
        }
        dgs  /= npee;
        dgs2 /= npee;
    }
    *sig = sqrt((dgs2-dgs*dgs)/(npee-1));
}

static void read_barsim(char *fn,barsim_t *ba)
{
    int  i;
    char **legend,*ptr;

    printf("Reading file '%s', lambda:",fn);
    ba->np = read_xvg_legend(fn,&ba->y,&ba->nset,&legend);
    snew(ba->lambda,ba->nset-1);
    for(i=0; i<ba->nset-1; i++)
    {
        if (legend[i] == 0)
        {
            gmx_fatal(FARGS,"There is no legend in file '%s', can not deduce lambda",fn);
        }
        ptr = strrchr(legend[i],' ');
        if ((i==0 && strstr(legend[i],"dG") == NULL) ||
            (i>0  && (strchr(legend[i],'D') == NULL ||
                      strchr(legend[i],'G') == NULL)) ||
            ptr == NULL)
        {
            gmx_fatal(FARGS,"There is no proper lambda legend in file '%s', can not deduce lambda",fn);
        }
        if (sscanf(ptr,"%lf",&ba->lambda[i]) != 1)
        {
            gmx_fatal(FARGS,"There is no proper lambda legend in file '%s', can not deduce lambda",fn);
        }
        printf(" %g",ba->lambda[i]);
    }
    printf("\n");
    
    /* Reorder the data */
    ba->t  = ba->y[0];
    for(i=1; i<ba->nset; i++)
    {
        ba->y[i-1] = ba->y[i];
    }
    for(i=0; i<ba->nset-1; i++)
    {
        sfree(legend[i]);
    }
    sfree(legend);
    ba->nset--;
}

int gmx_bar(int argc,char *argv[])
{
    static char *desc[] = {
        "g_bar does Bennett acceptance ratio free energy difference estimates.",
        "Multple dgdl files can be read at once, the (foreign) lambda values",
        "are determined from the xvg legends.",
        "A rough error estimate taking into account time correlations",
        "is done by splitting the data into n blocks and (re-)determining",
        "the free energy differences over those block and assuming",
        "those are independent."
    };
    static real begin=0,end=-1,temp=298;
    static int nblock=4;
    t_pargs pa[] = {
        { "-b",    FALSE, etREAL, {&begin},  "Begin time (ps)" },
        { "-e",    FALSE, etREAL, {&end},    "End time (ps)" },
        { "-temp", FALSE, etREAL, {&temp},   "Temperature (K)" },
        { "-nb",   FALSE, etINT,  {&nblock}, "The number of blocks for error estimation" }
    };
    
    t_filenm   fnm[] = {
        { efXVG, "-f", "dgdl",  ffRDMULT }
    };
#define NFILE asize(fnm)
    
    int      nfile,f,f2,fm;
    char     **fnms;
    barsim_t *ba,ba_tmp;
    double   dg_tot,var_tot,dg,sig;
    
    CopyRight(stderr,argv[0]);
    parse_common_args(&argc,argv,
                      PCA_CAN_VIEW,
                      NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
    
    nfile = opt2fns(&fnms,"-f",NFILE,fnm);
    if (nfile == 0)
    {
        gmx_fatal(FARGS,"No input files!");
    }

    snew(ba,nfile);
    for(f=0; f<nfile; f++)
    {
        read_barsim(fnms[f],&ba[f]);
    }
    printf("\n");

    /* Sort the data sets on lambda */
    for(f=0; f<nfile-1; f++)
    {
        fm = f;
        for(f2=f+1; f2<nfile; f2++)
        {
            if (ba[f2].lambda[0] == ba[fm].lambda[0])
            {
                gmx_fatal(FARGS,"There are multiple files with lambda = %g",
                          ba[fm].lambda[0]);
            }
            else if (ba[f2].lambda[0] < ba[fm].lambda[0])
            {
                fm = f2;
            }
        }
        ba_tmp = ba[f];
        ba[f]  = ba[fm];
        ba[fm] = ba_tmp;
    }
    
    dg_tot  = 0;
    var_tot = 0;
    for(f=0; f<nfile-1; f++)
    {
        calc_bar(&ba[f],&ba[f+1],begin,end,temp,0.0001,nblock,&dg,&sig);
        printf("lambda %4.2f - %4.2f, DG %9.5f (+-%.3f)\n",
               ba[f].lambda[0],ba[f+1].lambda[0],dg,sig);
        dg_tot  += dg;
        var_tot += sig*sig;
    }
    printf("\n");
    printf("total  %4.2f - %4.2f, DG %9.5f (+-%.3f)\n",
           ba[0].lambda[0],ba[nfile-1].lambda[0],dg_tot,sqrt(var_tot));

    /* do_view(opt2fn("-o",NFILE,fnm),"-nxy"); */
    
    thanx(stderr);
    
    return 0;
}
