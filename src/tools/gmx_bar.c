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
#include "gmx_ana.h"
#include "maths.h"

typedef struct {
    char   *filename;
    int    nset;
    int    np;
    int    begin;
    int    end;
    double temp;
    double *lambda;
    double *t;
    double **y;
} barsim_t;

/* calculated values */
typedef struct {
    barsim_t *a, *b; /* the simulation data */

    double lambda_a, lambda_b; /* the lambda values at a and b */

    double dg; /* the free energy difference */
    double dg_err; /* the free energy difference */

    double sa; /* relative entropy of b in state a */
    double sa_err; /* error in sa */
    double sb; /* relative entropy of a in state b */
    double sb_err; /* error in sb */

    double dg_stddev; /* expected dg stddev per sample */
    double dg_stddev_err; /* error in dg_stddev */
} barres_t;


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

static double calc_bar_lowlevel(int n1,double *W1,int n2,double *W2,
                                double delta_lambda,double temp,double tol)
{
    double kT,beta,beta_dl,M;
    double DG;
    int    i;
    double Wfac1,Wfac2,Wmin,Wmax;
    double DG0,DG1,DG2,dDG1;
    double sum1,sum2;
    
    kT   = BOLTZ*temp;
    beta = 1/kT;
    
    M = log((double)n1/(double)n2);

    if (delta_lambda == 0)
    {
        Wfac1 = beta;
        Wfac2 = beta;
    }
    else
    {
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
    
    /* For the comparison we can use twice the tolerance. */
    if (debug)
    {
        fprintf(debug,"DG %9.5f %9.5f\n",DG0,DG2);
    }
    while (DG2 - DG0 > 2*tol)
    {
        DG1 = 0.5*(DG0 + DG2);

        /*printf("Wfac1=%g, Wfac2=%g, beta=%g, DG1=%g\n",Wfac1,Wfac2,beta,
          DG1);*/
        dDG1 =
            calc_bar_sum(n1,W1,Wfac1, (M-DG1)) -
            calc_bar_sum(n2,W2,Wfac2,-(M-DG1));
        
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

static void calc_rel_entropy(int n1,double *W1,int n2,double *W2,
                             double delta_lambda, double temp, 
                             double dg, double *sa, double *sb)
{
    int i;
    double W_ab=0.;
    double W_ba=0.;
    double kT, beta;
    double Wfac1, Wfac2;

    kT   = BOLTZ*temp;
    beta = 1/kT;

    /* to ensure the work values are the same as during the delta_G */
    if (delta_lambda == 0)
    {
        Wfac1 = beta;
        Wfac2 = beta;
    }
    else
    {
        Wfac1 =  beta*delta_lambda;
        Wfac2 =  -beta*delta_lambda;
    }

    /* first calculate the average work in both directions */
    for(i=0;i<n1;i++)
    {
        W_ab += Wfac1*W1[i];
    }
    W_ab/=n1;
    for(i=0;i<n2;i++)
    {
        W_ba += Wfac2*W2[i]; 
    }
    W_ba/=n2;
   
    /* then calculate the relative entropies */
    *sa = (W_ab - dg);
    *sb = (W_ba + dg);
}

static void calc_dg_stddev(int n1, double *W1, int n2, double *W2, 
                             double delta_lambda, double temp, 
                             double dg, double *stddev)
{
    int i;
    double M;
    double sigmafact=0.;
    double kT, beta;
    double Wfac1, Wfac2;

    double nn1=n1; /* this makes the fraction in the *stddev eq below nicer */
    double nn2=n2;

    kT   = BOLTZ*temp;
    beta = 1/kT;

    /* to ensure the work values are the same as during the delta_G */
    if (delta_lambda == 0)
    {
        Wfac1 = beta;
        Wfac2 = beta;
    }
    else
    {
        Wfac1 =  beta*delta_lambda;
        Wfac2 =  -beta*delta_lambda;
    }

    M = log(nn1/nn2);

    /* calculate average in both directions */
    for(i=0;i<n1;i++)
    {
        sigmafact += 1./(2. + 2.*cosh((M + Wfac1*W1[i] - dg)));
    }
    for(i=0;i<n2;i++)
    {
        sigmafact += 1./(2. + 2.*cosh((M - Wfac2*W2[i] - dg)));
    }
    sigmafact /= (n1 + n2);
   
    /* Eq. 10 from 
       Shirts, Bair, Hooker & Pande, Phys. Rev. Lett 91, 140601 (2003): */
    *stddev = sqrt(((1./sigmafact) - ( (nn1+nn2)/nn1 + (nn1+nn2)/nn2 )));
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
    while (i < ba->nset &&
           !gmx_within_tol(ba->lambda[i],lambda,10*GMX_REAL_EPS))
    {
        i++;
    }
    if (i  == ba->nset)
    {
        gmx_fatal(FARGS,"Could not find a set for lambda = %g in the file '%s' of lambda = %g",lambda,ba->filename,ba->lambda[0]);
    }

    return i;
}

static void calc_bar(barsim_t *ba1,barsim_t *ba2,bool bUsedhdl,
                     double tol, int npee_min,int npee_max,
                     barres_t *br, bool *bEE, double *partsum)
{
    int np1,np2,s1,s2,npee,p;
    double delta_lambda; 
    double dg_sig2,sa_sig2,sb_sig2,stddev_sig2; /* intermediate variance values
                                                   for calculated quantities */
    br->a = ba1;
    br->b = ba2;
    br->lambda_a = ba1->lambda[0];
    br->lambda_b = ba2->lambda[0];

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

    br->dg = calc_bar_lowlevel(np1,ba1->y[s1]+ba1->begin,
                               np2,ba2->y[s2]+ba2->begin,
                               delta_lambda,ba1->temp,tol);


    calc_rel_entropy(np1, ba1->y[s1]+ba1->begin,
                     np2, ba2->y[s2]+ba2->begin,
                     delta_lambda, ba1->temp, br->dg, &(br->sa), &(br->sb));
    calc_dg_stddev(np1, ba1->y[s1]+ba1->begin,
                   np2, ba2->y[s2]+ba2->begin,
                   delta_lambda, ba1->temp, br->dg, &(br->dg_stddev) );


    dg_sig2 = 0;
    sa_sig2 = 0;
    sb_sig2 = 0;
    stddev_sig2 = 0;
    if (np1 >= npee_max && np2 >= npee_max)
    {
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
                dgp = calc_bar_lowlevel(np1/npee,
                                        ba1->y[s1]+ba1->begin+p*(np1/npee),
                                        np2/npee,
                                        ba2->y[s2]+ba2->begin+p*(np2/npee),
                                        delta_lambda,ba1->temp,tol);
                dgs  += dgp;
                dgs2 += dgp*dgp;

                partsum[npee*(npee_max+1)+p] += dgp;

                calc_rel_entropy(np1/npee, 
                                 ba1->y[s1]+ba1->begin+p*(np1/npee),
                                 np2/npee, 
                                 ba2->y[s2]+ba2->begin+p*(np2/npee),
                                 delta_lambda, ba1->temp, dgp, &sac, &sbc); 
                dsa  += sac;
                dsa2 += sac*sac;
                dsb  += sbc;
                dsb2 += sbc*sbc;
                calc_dg_stddev(np1/npee, 
                               ba1->y[s1]+ba1->begin+p*(np1/npee),
                               np2/npee, 
                               ba2->y[s2]+ba2->begin+p*(np2/npee),
                               delta_lambda, ba1->temp, dgp, &stddevc );

                dstddev  += stddevc;
                dstddev2 += stddevc*stddevc;
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

static void read_barsim(char *fn,double begin,double end,real temp,
                        barsim_t *ba)
{
    int  i;
    char *subtitle,**legend,*ptr;

    ba->filename = fn;

    printf("'%s' ",ba->filename);

    ba->np = read_xvg_legend(fn,&ba->y,&ba->nset,&subtitle,&legend);
    if (!ba->y)
    {
        gmx_fatal(FARGS,"File %s contains no usable data.",fn);
    }
    ba->t  = ba->y[0];

    get_begin_end(ba,begin,end,&ba->begin,&ba->end);
    printf("%.1f - %.1f, %6d points, lam:",
           ba->t[ba->begin],ba->t[ba->end-1],ba->end-ba->begin);

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
    static int nd=2,nbmin=5,nbmax=5;
    bool calc_s,calc_v;
    t_pargs pa[] = {
        { "-b",    FALSE, etREAL, {&begin},  "Begin time for BAR" },
        { "-e",    FALSE, etREAL, {&end},    "End time for BAR" },
        { "-temp", FALSE, etREAL, {&temp},   "Temperature (K)" },
        { "-prec", FALSE, etINT,  {&nd},     "The number of digits after the decimal point" },
        { "-nbmin",  FALSE, etINT,  {&nbmin}, "Minimum number of blocks for error estimation" },
        { "-nbmax",  FALSE, etINT,  {&nbmax}, "Maximum number of blocks for error estimation" }
    };
    
    t_filenm   fnm[] = {
        { efXVG, "-f",  "dhdl",   ffRDMULT },
        { efXVG, "-o",  "bar",    ffOPTWR },
        { efXVG, "-oi", "barint", ffOPTWR }
    };
#define NFILE asize(fnm)
    
    int      nfile,f,f2,fm,n1,nm;
    char     **fnms;
    barsim_t *ba,ba_tmp;
    barres_t *results;
    double   *partsum;
    double   prec,dg_tot,dg,sig;
    FILE     *fpb,*fpi;
    char     lamformat[20];
    char     dgformat[20],xvg2format[STRLEN],xvg3format[STRLEN],buf[STRLEN];
    char     ktformat[STRLEN], sktformat[STRLEN];
    char     kteformat[STRLEN], skteformat[STRLEN];
    output_env_t oenv;
    double   kT, beta;
    bool     result_OK=TRUE,bEE=TRUE;
    
    CopyRight(stderr,argv[0]);
    parse_common_args(&argc,argv,
                      PCA_CAN_VIEW,
                      NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv);
    
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


    snew(ba,nfile);
    snew(results,nfile-1);
    snew(partsum,(nbmax+1)*(nbmax+1));
    n1 = 0;
    nm = 0;
    for(f=0; f<nfile; f++)
    {
        read_barsim(fnms[f],begin,end,temp,&ba[f]);
        if (f > 0 && ba[f].temp != ba[0].temp)
        {
            printf("\nWARNING: temperature for file '%s' (%g) is not equal to that of file '%s' (%g)\n\n",fnms[f],ba[f].temp,fnms[0],ba[0].temp);
        }

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

    /* first calculate results */
    bEE = TRUE;
    for(f=0; f<nfile-1; f++)
    {
        /* Determine the free energy difference with a factor of 10
         * more accuracy than requested for printing.
         */
        calc_bar(&ba[f], &ba[f+1], n1>0, 0.1*prec, nbmin, nbmax,
                 &(results[f]), &bEE, partsum);
    }

    /* print results in kT */
    kT   = BOLTZ*ba[0].temp;
    beta = 1/kT;

    printf("\nTemperature: %g K\n", ba[0].temp);

    printf("\nDetailed results in kT (see help for explanation):\n\n");
    printf("%6s ", " lam_A");
    printf("%6s ", " lam_B");
    printf(sktformat,  "DG ");
    printf(skteformat, "+/- ");
    printf(sktformat,  "s_A ");
    printf(skteformat, "+/- " );
    printf(sktformat,  "s_B ");
    printf(skteformat, "+/- " );
    printf(sktformat,  "stdev ");
    printf(skteformat, "+/- ");
    printf("\n");
    for(f=0; f<nfile-1; f++)
    {
        printf(lamformat, results[f].lambda_a);
        printf(" ");
        printf(lamformat, results[f].lambda_b);
        printf(" ");
        printf(ktformat,  results[f].dg);
        printf(" ");
        printf(kteformat, results[f].dg_err);
        printf(" ");
        printf(ktformat,  results[f].sa);
        printf(" ");
        printf(kteformat, results[f].sa_err);
        printf(" ");
        printf(ktformat,  results[f].sb);
        printf(" ");
        printf(kteformat, results[f].sb_err);
        printf(" ");
        printf(ktformat,  results[f].dg_stddev);
        printf(" ");
        printf(kteformat, results[f].dg_stddev_err);
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
    for(f=0; f<nfile-1; f++)
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
        printf(lamformat, results[f].lambda_a);
        printf(" - ");
        printf(lamformat, results[f].lambda_b);
        printf(",   DG ");

        printf(dgformat,results[f].dg*kT);
        printf(" +/- ");
        printf(dgformat,results[f].dg_err*kT);

        printf("\n");
        dg_tot += results[f].dg;
    }
    printf("\n");
    printf("total  ");
    printf(lamformat, ba[0].lambda[0]);
    printf(" - ");
    printf(lamformat, ba[nfile-1].lambda[0]);
    printf(",   DG ");

    printf(dgformat,dg_tot*kT);
    if (bEE)
    {
        printf(" +/- ");
        printf(dgformat,bar_err(nbmin,nbmax,partsum)*kT);
    }
    printf("\n");

    if (fpi != NULL)
    {
        fprintf(fpi, xvg2format,
                ba[nfile-1].lambda[0], dg_tot);
        ffclose(fpi);
    }

    do_view(oenv,opt2fn_null("-o",NFILE,fnm),"-xydy");
    do_view(oenv,opt2fn_null("-oi",NFILE,fnm),"-xydy");
    
    thanx(stderr);
    
    return 0;
}
