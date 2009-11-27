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

    double dg_var; /* expected dg variance per sample */
    double dg_var_err; /* error in dg_var */
} barres_t;


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

    /* calculate minimum and maximum work to give an initial estimate of 
       delta G  as their average */
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

        /*printf("Wfac1=%g, Wfac2=%g, beta=%g, DG1=%g\n",Wfac1,Wfac2,beta,
          DG1);*/
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
        Wfac1 = 1;
        Wfac2 = 1;
    }
    else
    {
        Wfac1 =  delta_lambda;
        Wfac2 =  -delta_lambda;
    }

    /*printf("delta_lambda=%g, Wfac1=%g, Wfac2=%g, beta=%g\n", 
           delta_lambda, Wfac1, Wfac2, beta);*/
 
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
    *sa = beta*(W_ab - dg);
    *sb = beta*(W_ba + dg);
}

static void calc_dg_variance(int n1, double *W1, int n2, double *W2, 
                             double delta_lambda, double temp, 
                             double dg, double *var)
{
    int i;
    double M;
    double sigmafact=0.;
    double kT, beta;
    double Wfac1, Wfac2;

    kT   = BOLTZ*temp;
    beta = 1/kT;

    /* to ensure the work values are the same as during the delta_G */
    if (delta_lambda == 0)
    {
        Wfac1 = 1;
        Wfac2 = 1;
    }
    else
    {
        Wfac1 =  delta_lambda;
        Wfac2 =  -delta_lambda;
    }

    M=kT*log(((double)n1)/((double)n2));

    /* calculate average in both directions */
    for(i=0;i<n1;i++)
    {
        sigmafact += 1./(2. + 2.*cosh(beta*(M + Wfac1*W1[i] - dg)));
    }
    for(i=0;i<n2;i++)
    {
        sigmafact += 1./(2. + 2.*cosh(beta*(M + Wfac2*W2[i] - dg)));
    }
    sigmafact/=(n1+n2);
   
    /* Eq. 10 from 
       Shirts, Bair, Hooker & Pande, Phys. Rev. Lett 91, 140601 (2003): */
    *var = kT*kT*((1./sigmafact) - ( (n1+n2)/n1 + (n1+n2)/n2 ));
    /*printf("sigmafact=%g, var=%g\n",sigmafact, *var);*/
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
                     barres_t *br, bool calc_s, bool calc_var)
{
    int np1,np2,s1,s2,npee,p;
    double delta_lambda; 
    double dg_sig, sa_sig, sb_sig, var_sig; /* intermediate variance values for 
                                               calculated quantities */
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
                               delta_lambda,temp,tol);


    if (calc_s)
    {
        calc_rel_entropy(np1, ba1->y[s1]+ba1->begin,
                         np2, ba2->y[s2]+ba2->begin,
                         delta_lambda, temp, br->dg, &(br->sa), &(br->sb));
    }
    if (calc_var)
    {
        calc_dg_variance(np1, ba1->y[s1]+ba1->begin,
                         np2, ba2->y[s2]+ba2->begin,
                         delta_lambda, temp, br->dg, &(br->dg_var) );
    }


    dg_sig = 0;
    sa_sig = 0;
    sb_sig = 0;
    var_sig = 0;
    if (np1 >= npee1 && np2 >= npee1)
    {
        for(npee=npee0; npee<=npee1; npee++)
        {
            double dgs  = 0;
            double dgs2 = 0;
            double dsa  = 0;
            double dsb  = 0;
            double dsa2 = 0;
            double dsb2 = 0;
            double dvar = 0;
            double dvar2= 0;
 
            
            for(p=0; p<npee; p++)
            {
                double dgp;
                dgp = calc_bar_lowlevel(np1/npee,
                                        ba1->y[s1]+ba1->begin+p*(np1/npee),
                                        np2/npee,
                                        ba2->y[s2]+ba2->begin+p*(np2/npee),
                                        delta_lambda,temp,tol);
                dgs  += dgp;
                dgs2 += dgp*dgp;

                if (calc_s)
                {
                    double sac, sbc;
                    calc_rel_entropy(np1/npee, 
                                     ba1->y[s1]+ba1->begin+p*(np1/npee),
                                     np2/npee, 
                                     ba2->y[s2]+ba2->begin+p*(np2/npee),
                                     delta_lambda, temp, dgp, &sac, &sbc); 
                    dsa  += sac;
                    dsa2 += sac*sac;
                    dsb  += sbc;
                    dsb2 += sbc*sbc;
                }
                if (calc_var)
                {
                    double varc;
                    calc_dg_variance(np1/npee, 
                                     ba1->y[s1]+ba1->begin+p*(np1/npee),
                                     np2/npee, 
                                     ba2->y[s2]+ba2->begin+p*(np2/npee),
                                     delta_lambda, temp, dgp, &varc );

                    dvar  += varc;
                    dvar2 += varc*varc;
                }

            }
            dgs  /= npee;
            dgs2 /= npee;
            dg_sig += sqrt((dgs2-dgs*dgs)/(npee-1));

            if (calc_s)
            {
                dsa  /= npee;
                dsa2 /= npee;
                dsb  /= npee;
                dsb2 /= npee;
                sa_sig += sqrt((dsa2-dsa*dsa)/(npee-1));
                sb_sig += sqrt((dsb2-dsb*dsb)/(npee-1));
            }
            if (calc_var)
            {
                dvar  /= npee;
                dvar2 /= npee;
                var_sig += sqrt((dvar2-dvar*dvar)/(npee-1));
            }
        }
        br->dg_err = dg_sig/(npee1 - npee0 + 1);
        if (calc_s)
        {
            br->sa_err = sa_sig/(npee1 - npee0 + 1);
            br->sb_err = sb_sig/(npee1 - npee0 + 1);
        }
        if (calc_var)
        {
            br->dg_var_err = var_sig/(npee1 - npee0 + 1);
        }
 
    }
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
        "for a certain number of blocks.[BR]",
        "An estimate of the expected per-sample variance (as given in ",
        "Bennett's original BAR paper: Bennett, J. Comp. Phys. 22, p 245, ",    
        "(1976), Eq. 10) can be obtained with the [TT]-v[tt] option. ",
        "Note that this only gives an estimate of the quality of sampling, ",
        "not of the actual statistical error, because it assumes independent ",
        "samples.[BR]",
        "As a measure of phase space overlap, the relative entropy of ",
        "both states in each other's ensemble (i.e. the relative entropy s_ab ",
        "of the work samples of state b in the ensemble of state a, and vice ",
        "versa for s_ba) can be estimated using the [TT]-s[tt] option. The ",
        "shared entropy is a positive measure of the overlap of ",
        "the Boltzmann distributions of the two states that goes to zero ",
        "for identical distributions. See ",
        "Wu & Kofke, J. Chem. Phys. 123 084109 (2009) for more information."
    };
    static real begin=0,end=-1,temp=298;
    static int nd=2,nb0=4,nb1=5;
    bool calc_s,calc_v;
    t_pargs pa[] = {
        { "-b",    FALSE, etREAL, {&begin},  "Begin time for BAR" },
        { "-e",    FALSE, etREAL, {&end},    "End time for BAR" },
        { "-temp", FALSE, etREAL, {&temp},   "Temperature (K)" },
        { "-prec", FALSE, etINT,  {&nd},     "The number of digits after the decimal point" },
        { "-s",    FALSE, etBOOL, {&calc_s}, "Calculate relative entropy"},
        { "-v",    FALSE, etBOOL, {&calc_v}, "Calculate expected per-sample variance"},
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
    barres_t *results;
    double   prec,dg_tot,var_tot,dg,sig;
    FILE     *fp;
    char     dgformat[20],xvgformat[STRLEN],buf[STRLEN];
    output_env_t oenv;
    
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
    sprintf( dgformat,"%%%d.%df",3+nd,nd);
    sprintf(xvgformat,"%s %s %s\n","%g",dgformat,dgformat);


    snew(ba,nfile);
    snew(results,nfile-1);
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
                      "\\8l\\4",unit_energy,oenv);
        if (get_print_xvgr_codes(oenv))
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

        calc_bar(&ba[f], &ba[f+1], n1>0, temp, prec, nb0, nb1,
                 &(results[f]), calc_s, calc_v);

        printf("lambda %4.2f - %4.2f, DG ", results[f].lambda_a,
                                            results[f].lambda_b);
        printf(dgformat,results[f].dg);
        printf(" err");
        printf(dgformat,results[f].dg_err);
        if (calc_s)
        {
            printf("   s_ab "); 
            printf(dgformat, results[f].sa);
            printf(" err"); 
            printf(dgformat, results[f].sa_err);
            printf("  s_ba "); 
            printf(dgformat, results[f].sb);
            printf(" err"); 
            printf(dgformat, results[f].sb_err);
        }
        if (calc_v)
        {
            printf("   var est ");
            printf(dgformat, results[f].dg_var);
            printf(" err"); 
            printf(dgformat, results[f].dg_var_err);
        }
        printf("\n");
        dg_tot  += results[f].dg;
        var_tot += results[f].dg_err*results[f].dg_err;
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

    do_view(oenv,opt2fn_null("-o",NFILE,fnm),"-xydy");
    
    thanx(stderr);
    
    return 0;
}
