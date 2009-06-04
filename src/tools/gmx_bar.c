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
    int    begin;
    int    end;
    double *lambda;
    double *t;
    double **y;
} barsim_t;

static double calc_bar_sum(int n,double *W,double beta_Wfac,double sbMmDG)
{
    int    i;
    double sum;
    
    sum = 0;
    
    for(i=0; i<n; i++)
    {
        sum += 1/(1 + exp(beta_Wfac*W[i] + sbMmDG));
    }
    
    return sum;
}

static double calc_bar_lowlevel(int n1,double *W1,int n2,double *W2,
                                double delta_lambda,real temp,double prec)
{
    double kT,beta,beta_dl,M;
    double DG;
    int    i;
    double Wfac1,Wfac2,Wmin,Wmax;
    double DG0,DG1,DG2,dDG1;
    double sum1,sum2;
    
    kT   = BOLTZ*temp;
    beta = 1/kT;
    
    M = kT*log((double)n1/(double)n2);

    if (delta_lambda == 0)
    {
        Wfac1 = 1;
        Wfac2 = 1;
    }
    else
    {
        Wfac1 =  delta_lambda;
        Wfac2 = -delta_lambda;
    }

    Wmin = W1[0];
    Wmax = W1[0];
    for(i=0; i<n1; i++)
    {
        Wmin = min(Wmin,W1[i]*Wfac1);
        Wmax = max(Wmax,W1[i]*Wfac1);
    }
    for(i=0; i<n2; i++)
    {
        Wmin = min(Wmin,-W2[i]*Wfac2);
        Wmax = max(Wmax,-W2[i]*Wfac2);
    }
    DG0 = Wmin;
    DG2 = Wmax;
    
    /* For the comparison we can use twice the tolerance */
    if (debug)
    {
        fprintf(debug,"DG %9.5f %9.5f\n",DG0,DG2);
    }
    while (DG2 - DG0 > 2*prec)
    {
        DG1 = 0.5*(DG0 + DG2);
        dDG1 =
            calc_bar_sum(n1,W1,beta*Wfac1, beta*(M-DG1)) -
            calc_bar_sum(n2,W2,beta*Wfac2,-beta*(M-DG1));
        
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
        gmx_fatal(FARGS,"Some data end before the start time %g",begin);
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

static void calc_bar(barsim_t *ba1,barsim_t *ba2,bool bUsedhdl,
                     real temp,double tol,
                     int npee0,int npee1,
                     double *dg,double *sig)
{
    int np1,np2,s1,s2,npee,p;
    double delta_lambda,dgp,dgs,dgs2;

    if (bUsedhdl)
    {
        s1 = 0;
        s2 = 0;

        delta_lambda = ba2->lambda[0] - ba1->lambda[0];
    }
    else
    {
        s1 = get_lam_set(ba1,ba2->lambda[0]);
        s2 = get_lam_set(ba2,ba1->lambda[0]);

        delta_lambda = 0;
    }

    np1 = ba1->end - ba1->begin;
    np2 = ba2->end - ba2->begin;

    *dg = calc_bar_lowlevel(np1,ba1->y[s1]+ba1->begin,
                            np2,ba2->y[s2]+ba2->begin,
                            delta_lambda,temp,tol);

    *sig = 0;
    if (np1 >= npee1 && np2 >= npee1)
    {
        for(npee=npee0; npee<=npee1; npee++)
        {
            dgs  = 0;
            dgs2 = 0;
            for(p=0; p<npee; p++)
            {
                dgp = calc_bar_lowlevel(np1/npee,
                                        ba1->y[s1]+ba1->begin+p*(np1/npee),
                                        np2/npee,
                                        ba2->y[s2]+ba2->begin+p*(np2/npee),
                                        delta_lambda,temp,tol);
                dgs  += dgp;
                dgs2 += dgp*dgp;
            }
            dgs  /= npee;
            dgs2 /= npee;
            *sig += sqrt((dgs2-dgs*dgs)/(npee-1));
        }
        *sig /= npee1 - npee0 + 1;
    }
}

static double legend2lambda(char *fn,const char *legend,bool bdhdl)
{
    double lambda=0;
    char   *ptr;

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

static double filename2lambda(char *fn)
{
    double lambda;
    char   *ptr,*endptr;
    
    ptr = fn;
    while (ptr[0] != '\0' && !isdigit(ptr[0]))
    {
        ptr++;
    }
    if (!isdigit(ptr[0]))
    {
        gmx_fatal(FARGS,"While trying to read the lambda value from the filename: filename '%s' does not contain a number",fn);
    }
    if (ptr > fn && fn[ptr-fn-1] == '-')
    {
        ptr--;
    }

    lambda = strtod(ptr,&endptr);
    if (endptr == ptr)
    {
        gmx_fatal(FARGS,"Malformed number in filename '%s'",fn);
    }

    return lambda;
}

static void read_barsim(char *fn,double begin,double end,barsim_t *ba)
{
    int  i;
    char **legend,*ptr;

    printf("'%s' ",fn);

    ba->np = read_xvg_legend(fn,&ba->y,&ba->nset,&legend);
    ba->t  = ba->y[0];

    get_begin_end(ba,begin,end,&ba->begin,&ba->end);
    printf("%.1f - %.1f, %6d points, lam:",
           ba->t[ba->begin],ba->t[ba->end-1],ba->end-ba->begin);

    snew(ba->lambda,ba->nset-1);
    if (legend == NULL)
    {
        /* Check if we have a single set, nset=2 means t and dH/dl */
        if (ba->nset == 2)
        {
            /* Deduce lambda from the file name */
            ba->lambda[0] = filename2lambda(fn);
            printf(" %g",ba->lambda[0]);
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
            printf(" %g",ba->lambda[i]);
        }
    }
    printf("\n");
    
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

int gmx_bar(int argc,char *argv[])
{
    static const char *desc[] = {
        "g_bar calculates free energy difference estimates through",
        "Bennett's acceptance ratio method.",
        "Input option [TT]-f[tt] expects multiple dhdl files.",
        "Two types of input files are supported:[BR]",
        "* Files with only one y-value, for such files it is assumed",
        "that the y-value is dH/dlambda and that the Hamiltonian depends",
        "linearly on lambda. The lambda value of the simulation is inferred",
        "from the legend if present, otherwise from a number in the file name.",
        "[BR]",
        "* Files with more than one y-value. The files should have columns",
        "with dH/dlambda and Delta lambda. The lambda values are inferred",
        "from the legends:",
        "lambda of the simulation from the legend of dH/dlambda",
        "and the foreign lambda's from the legends of Delta H.[BR]",
        "The lambda of the simulation is parsed from the legend containing",
        "the string dH, the foreign lambda's from the legend containing",
        "the capitalized letters D and H.[BR]",
        "The free energy estimates are determined using BAR with bisection,",
        "the precision of the output is set with [TT]-prec[tt].",
        "An error estimate taking into account time correlations",
        "is made by splitting the data into blocks and determining",
        "the free energy differences over those blocks and assuming",
        "the blocks are independent.",
        "The final error estimate is determined from the average variance",
        "over 4 and 5 blocks to lower the chance of an low error estimate",
        "when the differences between the block are coincidentally low",
        "for a certain number of blocks."
    };
    static real begin=0,end=-1,temp=298;
    static int nd=2,nb0=4,nb1=5;
    t_pargs pa[] = {
        { "-b",    FALSE, etREAL, {&begin},  "Begin time for BAR" },
        { "-e",    FALSE, etREAL, {&end},    "End time for BAR" },
        { "-temp", FALSE, etREAL, {&temp},   "Temperature (K)" },
        { "-prec", FALSE, etINT,  {&nd},     "The number of digits after the decimal point" },
        { "-nb0",  FALSE, etINT,  {&nb0}, "HIDDENMinimum number of blocks for error estimation" },
        { "-nb1",  FALSE, etINT,  {&nb1}, "HIDDENMaximum number of blocks for error estimation" }
    };
    
    t_filenm   fnm[] = {
        { efXVG, "-f", "dhdl",  ffRDMULT },
        { efXVG, "-o", "bar",   ffOPTWR }
    };
#define NFILE asize(fnm)
    
    int      nfile,f,f2,fm,n1,nm;
    char     **fnms;
    barsim_t *ba,ba_tmp;
    double   prec,dg_tot,var_tot,dg,sig;
    FILE     *fp;
    char     dgformat[20],xvgformat[STRLEN],buf[STRLEN];
    
    CopyRight(stderr,argv[0]);
    parse_common_args(&argc,argv,
                      PCA_CAN_VIEW,
                      NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
    
    nfile = opt2fns(&fnms,"-f",NFILE,fnm);
    if (nfile == 0)
    {
        gmx_fatal(FARGS,"No input files!");
    }

    if (nd < 0)
    {
        gmx_fatal(FARGS,"Can not have negative number of digits");
    }
    prec = pow(10,-nd);
    sprintf( dgformat,"%%%d.%df",3+nd,nd);
    sprintf(xvgformat,"%s %s %s\n","%g",dgformat,dgformat);


    snew(ba,nfile);
    n1 = 0;
    nm = 0;
    for(f=0; f<nfile; f++)
    {
        read_barsim(fnms[f],begin,end,&ba[f]);
        
        if (ba[f].nset == 0)
        {
            gmx_fatal(FARGS,"File '%s' contains less than two columns",fnms[f]);
        }
        else if (ba[f].nset == 1)
        {
            n1++;
        }
        else
        {
            nm++;
        }
    }
    printf("\n");

    if (n1 > 0 && nm > 0)
    {
        gmx_fatal(FARGS,"Some dhdl files contain only one value (assuming dH/dl), while others contain multiple values (assuming dH/dl and Delta H), will not proceed because of possible inconsistencies");
    }

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
    
    if (n1 > 0)
    {
        printf("Only one y value in all files,\n"
               "assuming the Hamiltonian depends linearly on lambda\n\n");
    }

    if (opt2bSet("-o",NFILE,fnm))
    {
        sprintf(buf,"%s (%s)","\\8D\\4G",unit_energy);
        fp = xvgropen(opt2fn("-o",NFILE,fnm),"Free energy differences",
                      "\\8l\\4",unit_energy);
        if (bPrintXvgrCodes())
        {
            fprintf(fp,"@TYPE xydy\n");
        }
    }
    else
    {
        fp = NULL;
    }

    dg_tot  = 0;
    var_tot = 0;
    for(f=0; f<nfile-1; f++)
    {
        if (fp != NULL)
        {
            fprintf(fp,xvgformat,
                    ba[f].lambda[0],dg_tot,sqrt(var_tot));
        }

        calc_bar(&ba[f],&ba[f+1],n1>0,temp,prec,nb0,nb1,&dg,&sig);

        printf("lambda %4.2f - %4.2f, DG ",ba[f].lambda[0],ba[f+1].lambda[0]);
        printf(dgformat,dg);
        printf(" err");
        printf(dgformat,sig);
        printf("\n");
        dg_tot  += dg;
        var_tot += sig*sig;
    }
    printf("\n");
    printf("total  %4.2f - %4.2f, DG ",ba[0].lambda[0],ba[nfile-1].lambda[0]);
    printf(dgformat,dg_tot);
    printf(" err");
    printf(dgformat,sqrt(var_tot));
    printf("\n");

    if (fp != NULL)
    {
        fprintf(fp,xvgformat,
                ba[nfile-1].lambda[0],dg_tot,sqrt(var_tot));
        fclose(fp);
    }

    do_view(opt2fn_null("-o",NFILE,fnm),"-xydy");
    
    thanx(stderr);
    
    return 0;
}
