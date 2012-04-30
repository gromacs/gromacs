/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: tune_dip.c,v 1.39 2009/04/12 21:24:26 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREADS
#include "tmpi.h"
#endif
#include "maths.h"
#include "macros.h"
#include "copyrite.h"
#include "bondf.h"
#include "string2.h"
#include "smalloc.h"
#include "strdb.h"
#include "sysstuff.h"
#include "confio.h"
#include "physics.h"
#include "futil.h"
#include "statutil.h"
#include "vec.h"
#include "3dview.h"
#include "txtdump.h"
#include "readinp.h"
#include "names.h"
#include "vec.h"
#include "atomprop.h"
#include "xvgr.h"
#include "mdatoms.h"
#include "force.h"
#include "vsite.h"
#include "shellfc.h"
#include "network.h"
#include "viewit.h"
#include "gmx_random.h"
#include "gmx_wallcycle.h"
#include "gmx_statistics.h"
#include "convparm.h"
#include "gpp_atomtype.h"
#include "grompp.h"
#include "gen_ad.h"
#include "slater_integrals.h"
#include "gentop_qgen.h"
#include "gentop_core.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "molprop.h"
#include "molprop_xml.h"
#include "molprop_util.h"
#include "molselect.h"
#include "mtop_util.h"
#include "gentop_comm.h"
#include "nmsimplex.h"
#include "mymol.h"

static void print_stats(FILE *fp,const char *prop,gmx_stats_t lsq,gmx_bool bHeader,
                        char *xaxis,char *yaxis)
{
    real a,da,b,db,chi2,rmsd,Rfit;
    int  n;
    
    if (bHeader)
    {
        fprintf(fp,"Fitting data to y = ax+b, where x = %s and y = %s\n",
                xaxis,yaxis);
        fprintf(fp,"%-12s %5s %13s %13s %8s %8s\n",
                "Property","N","a","b","R","RMSD");
        fprintf(fp,"---------------------------------------------------------------\n");
    }
    gmx_stats_get_ab(lsq,elsqWEIGHT_NONE,&a,&b,&da,&db,&chi2,&Rfit);
    gmx_stats_get_rmsd(lsq,&rmsd);
    gmx_stats_get_npoints(lsq,&n);
    fprintf(fp,"%-12s %5d %6.3f(%5.3f) %6.3f(%5.3f) %7.2f%% %8.4f\n",
            prop,n,a,da,b,db,Rfit*100,rmsd);
}

static void print_lsq_set(FILE *fp,gmx_stats_t lsq)
{
    real   x,y;
    
    fprintf(fp,"@type xy\n");
    while (gmx_stats_get_point(lsq,&x,&y,NULL,NULL,0) == estatsOK)
    {
        fprintf(fp,"%10g  %10g\n",x,y);
    }
    fprintf(fp,"&\n");
}

static void xvgr_symbolize(FILE *xvgf,int nsym,const char *leg[],
                           const output_env_t oenv)
{
    int i;

    xvgr_legend(xvgf,nsym,leg,oenv);
    for(i=0; (i<nsym); i++)
    {
        xvgr_line_props(xvgf,i,elNone,ecBlack+i,oenv);
        fprintf(xvgf,"@ s%d symbol %d\n",i,i+1);
    }        
}

static void calc_moldip_deviation(t_moldip *md)
{
    int    i,j,count,atomnr;
    double qq,qtot,rr2,ener[ermsNR],etot[ermsNR];
    real   t = 0;
    rvec   mu_tot = {0,0,0};
    gmx_enerdata_t *epot;
    tensor force_vir={{0,0,0},{0,0,0},{0,0,0}};
    t_nrnb   my_nrnb;
    gmx_wallcycle_t wcycle;
    gmx_bool     bConverged;
    t_mymol *mymol;
    int     eQ;
    gmx_mtop_atomloop_all_t aloop;
    t_atom *atom; 
    int    at_global,resnr;
    
    if (PAR(md->cr)) 
    {
        gmx_bcast(sizeof(md->bDone),&md->bDone,md->cr);
        gmx_bcast(sizeof(md->bFinal),&md->bFinal,md->cr);
    }
    if (md->bDone)
        return;
    if (PAR(md->cr)) 
    {
        gmx_poldata_comm_eemprops(md->pd,md->cr);
    }
    init_nrnb(&my_nrnb);
    snew(epot,1);
  
    wcycle  = wallcycle_init(stdout,0,md->cr,1);
    for(j=0; (j<ermsNR); j++)
    {
        etot[j] = 0;
        md->ener[j] = 0;
    }
    for(i=0; (i<md->nmol); i++) {
        mymol = &(md->mymol[i]);
        if ((mymol->eSupport == eSupportLocal) ||
            (md->bFinal && (mymol->eSupport == eSupportRemote)))
        {
            /* Reset energies */
            for(j=0; (j<ermsNR); j++)
                ener[j] = 0;
            
            /* Now optimize the shell positions */
            if (mymol->shell) {
                atoms2md(&mymol->mtop,&(mymol->ir),0,NULL,0,
                         mymol->mtop.natoms,mymol->md);
                count = 
                    relax_shell_flexcon(debug,md->cr,FALSE,0,
                                        &(mymol->ir),TRUE,
                                        GMX_FORCE_ALLFORCES,FALSE,
                                        mymol->ltop,NULL,NULL,NULL,
                                        NULL,&(mymol->state),
                                        mymol->f,force_vir,mymol->md,
                                        &my_nrnb,wcycle,NULL,
                                        &(mymol->mtop.groups),
                                        mymol->shell,mymol->fr,FALSE,t,mu_tot,
                                        mymol->mtop.natoms,&bConverged,NULL,NULL);
            }

            for(j=0; (j<ermsNR); j++)
                etot[j] += ener[j];
        }
    }
    for(j=0; (j<ermsTOT); j++) {
        md->ener[j]       += md->fc[j]*etot[j]/md->nmol_support;
        md->ener[ermsTOT] += md->ener[j];
    }
    sfree(epot);
    if (debug)
    {
        fprintf(debug,"ENER:");
        for(j=0; (j<ermsNR); j++)
            fprintf(debug,"  %8.3f",etot[j]);
        fprintf(debug,"\n");
    }
    /* Global sum energies */
    if (PAR(md->cr)) 
    {
#ifdef GMX_DOUBLE
        gmx_sumd(ermsNR,md->ener,md->cr);
#else
        gmx_sumf(ermsNR,md->ener,md->cr);
#endif
    }
}

static double dipole_function(void *params,double v[])
{
    t_moldip *md = (t_moldip *) params;
    int      i,j,k,zz,nzeta;
    double   chi0,z,J0,bounds=0;
    char     *name,*zeta,*qstr,*rowstr;
    char     zstr[STRLEN],buf[STRLEN];
    
#define HARMONIC(x,xmin,xmax) (x < xmin) ? (sqr(x-xmin)) : ((x > xmax) ? (sqr(x-xmax)) : 0)
  
    /* Set parameters in eem record. There is a penalty if parameters
     * go out of bounds as well.
     */
    k=0;
    while ((name = opt_index_count(md->ic)) != NULL)
    {
        J0 = v[k++];
        bounds += HARMONIC(J0,md->J0_0,md->J0_1);
        if (strcasecmp(name,md->fixchi) != 0) 
        {
            chi0 = v[k++];
            bounds += HARMONIC(chi0,md->Chi0_0,md->Chi0_1);
        }
        else 
            chi0 = gmx_poldata_get_chi0(md->pd,md->iModel,name);
    
        qstr = gmx_poldata_get_qstr(md->pd,md->iModel,name);
        rowstr = gmx_poldata_get_rowstr(md->pd,md->iModel,name);
        nzeta = gmx_poldata_get_nzeta(md->pd,md->iModel,name);
        zstr[0] = '\0';
        for(zz=0; (zz<nzeta); zz++) 
        {
            z = gmx_poldata_get_zeta(md->pd,md->iModel,name,zz);
            if ((0 != z) && (md->bFitZeta))
            {
                z = v[k++];
                bounds += HARMONIC(z,md->w_0,md->w_1);
            }
            sprintf(buf,"  %g",z);
            strcat(zstr,buf);
        }
        gmx_poldata_set_eemprops(md->pd,md->iModel,name,J0,chi0,
                                 zstr,qstr,rowstr);
    }
    if (md->bOptHfac) 
    {
        md->hfac = v[k++];
        if (md->hfac >  md->hfac0) 
            bounds += 100*sqr(md->hfac - md->hfac0);
        else if (md->hfac < -md->hfac0)
            bounds += 100*sqr(md->hfac + md->hfac0);
    }
    for(j=0; (j<ermsNR); j++)
        md->ener[j] = 0;
    calc_moldip_deviation(md);
  
    /* This contribution is not scaled with force constant because
     * it are essential to charge generation convergence and can hence
     * not be left out of the optimization.
     */
    md->ener[ermsBOUNDS] += bounds;
    md->ener[ermsTOT]    += bounds;
    
    return md->ener[ermsTOT];
}

static real guess_new_param(real x,real step,real x0,real x1,gmx_rng_t rng,
                            gmx_bool bRandom)
{
    real r = gmx_rng_uniform_real(rng);
  
    if (bRandom) 
        x = x0+(x1-x0)*r;
    else
        x = x*(1-step+2*step*r);

    if (x < x0)
        return x0;
    else if (x > x1)
        return x1;
    else
        return x;
}

static int guess_all_param(FILE *fplog,t_moldip *md,int run,int iter,real stepsize,
                           gmx_bool bRandom,gmx_rng_t rng,
                           double orig_param[],double test_param[])
{
    double J00,xxx,chi0,zeta;
    char   *name,*qstr,*rowstr;
    char   zstr[STRLEN],buf[STRLEN];
    gmx_bool   bStart = (/*(run == 0) &&*/ (iter == 0));
    gmx_bool   bRand  = bRandom && (iter == 0);
    int    zz,nzeta,k = 0;
    
    fprintf(fplog,"%-5s %10s %10s %10s Run/Iter %d/%d - %s randomization\n","Name",
            "J00","chi0","zeta",run,iter,
            bRand ? "Complete" : (bStart ? "Initial" : "Limited"));
    while ((name = opt_index_count(md->ic)) != NULL) 
    {
        if (bStart)
        {
            J00 = gmx_poldata_get_j00(md->pd,md->iModel,name);
            xxx = guess_new_param(J00,stepsize,md->J0_0,md->J0_1,rng,bRand);
            if (bRand)
                orig_param[k] = xxx;
            else
                orig_param[k] = J00;
            J00 = xxx;
        }
        else 
            J00 = guess_new_param(orig_param[k],stepsize,md->J0_0,md->J0_1,rng,bRand);
        test_param[k++] = J00;
        
        chi0 = gmx_poldata_get_chi0(md->pd,md->iModel,name);
        if (strcasecmp(name,md->fixchi) != 0) 
        {
            if (bStart)
            {
                xxx = guess_new_param(chi0,stepsize,md->Chi0_0,md->Chi0_1,rng,bRand);
                if (bRand)
                    orig_param[k] = xxx;
                else
                    orig_param[k] = chi0;
                chi0 = xxx;
            }
            else 
                chi0 = guess_new_param(orig_param[k],stepsize,md->Chi0_0,md->Chi0_1,rng,bRand);
            test_param[k++] = chi0;
        }
        if ((qstr = gmx_poldata_get_qstr(md->pd,md->iModel,name)) == NULL)
            gmx_fatal(FARGS,"No qstr for atom %s model %d\n",name,md->iModel);
        if ((rowstr = gmx_poldata_get_rowstr(md->pd,md->iModel,name)) == NULL)
            gmx_fatal(FARGS,"No rowstr for atom %s model %d\n",name,md->iModel);
        nzeta   = gmx_poldata_get_nzeta(md->pd,md->iModel,name);
        zstr[0] = '\0';
        for(zz=0; (zz<nzeta); zz++) 
        {
            zeta = gmx_poldata_get_zeta(md->pd,md->iModel,name,zz);
            if ((md->bFitZeta) && (0 != zeta))
            {
                if (bStart)
                {
                    xxx = guess_new_param(zeta,stepsize,md->w_0,md->w_1,rng,bRand);
                    orig_param[k] = (bRand) ? xxx : zeta;
                    zeta = xxx;
                }
                else 
                    zeta = guess_new_param(orig_param[k],stepsize,md->w_0,md->w_1,rng,bRand);
                test_param[k++] = zeta;            
            }
            sprintf(buf,"  %10g",zeta);
            strcat(zstr,buf);
        }
        gmx_poldata_set_eemprops(md->pd,md->iModel,name,J00,chi0,
                                 zstr,qstr,rowstr);
        fprintf(fplog,"%-5s %10g %10g %10s\n",name,J00,chi0,zstr);
    }
    fprintf(fplog,"\n");
    fflush(fplog);
    if (md->bOptHfac) 
        test_param[k++] = md->hfac;
    return k;
}

static void optimize_moldip(FILE *fp,FILE *fplog,const char *convfn,
                            t_moldip *md,int maxiter,real tol,
                            int nrun,int reinit,real stepsize,int seed,
                            gmx_bool bRandom,real stol,output_env_t oenv)
{
    FILE   *cfp=NULL;
    double chi2,chi2_min,wj,rms_nw;
    int    nzeta,zz;
    int    status = 0;
    int    i,k,index,n,nparam;
    double *test_param,*orig_param,*best_param,*start;
    gmx_bool   bMinimum=FALSE;
    double J00,chi0,zeta;
    char   *name,*qstr,*rowstr;
    char   zstr[STRLEN],buf[STRLEN];
    gmx_rng_t rng;
  
    if (MASTER(md->cr)) 
    {
        rng = gmx_rng_init(seed);
  
        nparam = 0;
        while ((name = opt_index_count(md->ic)) != NULL) 
        {
            /* One parameter for J00 and one for chi0 */
            nparam++;
            if (strcasecmp(name,md->fixchi) != 0) 
                nparam++;
            if (md->bFitZeta) 
            {
                nzeta  = gmx_poldata_get_nzeta(md->pd,md->iModel,name);
                for(i=0; (i<nzeta); i++)
                {
                    zeta = gmx_poldata_get_zeta(md->pd,md->iModel,name,i);
                    if (zeta > 0)
                        nparam++;
                }
            }
        }
        /* Check whether we have to optimize the fudge factor for J00 H */
        if (md->hfac != 0)
            nparam++;
        snew(test_param,nparam+1);
        snew(orig_param,nparam+1);
        snew(best_param,nparam+1);
        
        /* Starting point */
        snew(start,nparam);
        
        /* Monitor convergence graphically */
        if (NULL != convfn) 
        {
            cfp = xvgropen(convfn,"Convergence","Value","Iter",oenv);
        }
        chi2_min = GMX_REAL_MAX; 
        for(n=0; (n<nrun); n++) {    
            if ((NULL != fp) && (0 == n)) {
                fprintf(fp,"\nStarting run %d out of %d\n",n+1,nrun);
                fprintf(fp,"%5s %8s %8s %8s %8s %8s %8s %8s\n",
                        "Run","d2","Total","Bounds",
                        "Dipole","Quad.","Charge","ESP");
            }
            if (NULL != cfp)
                fflush(cfp);
                    
            k = guess_all_param(fplog,md,n,0,stepsize,bRandom,rng,
                                orig_param,test_param);
            if (k != nparam)
                gmx_fatal(FARGS,"Inconsistency in guess_all_param: k = %d, should be %d",k,nparam);
            for(k=0; (k<nparam); k++) 
            {
                start[k] = test_param[k];
            }
                
            nmsimplex(cfp,(void *)md,dipole_function,start,nparam,
                      tol,1,maxiter,&chi2);
                
            if (chi2 < chi2_min) {
                bMinimum = TRUE;
                /* Print convergence if needed */
                if (NULL != cfp) 
                    fprintf(cfp,"%5d  ",n*maxiter);
                for(k=0; (k<nparam); k++)
                {
                    best_param[k] = start[k];
                    if (NULL != cfp)
                        fprintf(cfp," %10g",best_param[k]);
                }
                if (NULL != cfp)
                    fprintf(cfp,"\n");
                chi2_min   = chi2;
            }
            
            if (fp) 
                fprintf(fp,"%5d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                        n,chi2,
                        sqrt(md->ener[ermsTOT]),
                        sqrt(md->ener[ermsBOUNDS]),
                        sqrt(md->ener[ermsMU]),
                        sqrt(md->ener[ermsQUAD]),
                        sqrt(md->ener[ermsCHARGE]),
                        sqrt(md->ener[ermsESP]));
        }
        
        if (bMinimum) {
            for(k=0; (k<nparam); k++)
                start[k] = best_param[k];
            k = 0;
            while ((name = opt_index_count(md->ic)) != NULL) 
            {
                J00    = start[k++];
                chi0   = gmx_poldata_get_chi0(md->pd,md->iModel,name);
                if (strcasecmp(name,md->fixchi) != 0) 
                    chi0 = start[k++];
                qstr   = gmx_poldata_get_qstr(md->pd,md->iModel,name);
                rowstr = gmx_poldata_get_rowstr(md->pd,md->iModel,name);
                nzeta  = gmx_poldata_get_nzeta(md->pd,md->iModel,name);
                zstr[0] = '\0';
                for(zz=0; (zz<nzeta); zz++)
                {
                    zeta = gmx_poldata_get_zeta(md->pd,md->iModel,name,zz);
                    if ((0 != zeta) && md->bFitZeta)
                        zeta = start[k++];
                    sprintf(buf," %g",zeta);
                    strcat(zstr,buf);
                }
                gmx_poldata_set_eemprops(md->pd,md->iModel,name,J00,chi0,
                                         zstr,qstr,rowstr);
            }
            if (md->bOptHfac) 
                md->hfac = start[k++];
            gmx_assert(k,nparam);
            
            calc_moldip_deviation(md);
            chi2  = sqrt(md->ener[ermsTOT]);
            md->bFinal = TRUE;
            calc_moldip_deviation(md);
            if (fplog)
            {
                fprintf(fplog,"\nMinimum value for RMSD during optimization: %.3f.\n",chi2);
            }
        }
        md->bDone = TRUE;
        calc_moldip_deviation(md);
        
        if (NULL != cfp)
            fclose(cfp);
    }
    else 
    {
        /* Slave calculators */
        do {
            calc_moldip_deviation(md);
        } while (!md->bDone);
    }
}

static real quality_of_fit(real chi2,int N)
{
    return -1;
}

int main(int argc, char *argv[])
{
    static const char *desc[] = {
        "tune_dip read a series of molecules and corresponding experimental",
        "dipole moments from a file, and tunes parameters in an algorithm",
        "until the experimental dipole moments are reproduced by the",
        "charge generating algorithm AX as implemented in the gentop program.[PAR]",
        "Minima and maxima for the parameters can be set, these are however",
        "not strictly enforced, but rather they are penalized with a harmonic",
        "function, for which the force constant can be set explicitly.[PAR]",
        "At every reinit step parameters are changed by a random amount within",
        "the fraction set by step size, and within the boundaries given",
        "by the minima and maxima. If the [TT]-random[tt] flag is",
        "given a completely random set of parameters is generated at the start",
        "of each run. At reinit steps however, the parameters are only changed",
        "slightly, in order to speed-up local search but not global search."
        "In other words, complete random starts are done only at the beginning of each",
        "run, and only when explicitly requested.[PAR]",
        "The absolut dipole moment of a molecule remains unchanged if all the",
        "atoms swap the sign of the charge. To prevent this kind of mirror",
        "effects a penalty is added to the square deviation ",
        "if hydrogen atoms have a negative charge. Similarly a penalty is",
        "added if atoms from row VI or VII in the periodic table have a positive",
        "charge. The penalty is equal to the force constant given on the command line",
        "time the square of the charge.[PAR]",
        "One of the electronegativities (chi) is redundant in the optimization,",
        "only the relative values are meaningful.",
        "Therefore by default we fix the value for hydrogen to what is written",
        "in the eemprops.dat file (or whatever is given with the [tt]-d[TT] flag).",
        "A suitable value would be 2.3, the original, value due to Pauling,",
        "this can by overridden by setting the [tt]-fixchi[TT] flag to something else (e.g. a non-existing atom).[PAR]",
        "A selection of molecules into a training set and a test set (or ignore set)",
        "can be made using option [TT]-sel[tt]. The format of this file is:[BR]",
        "iupac|Train[BR]",
        "iupac|Test[BR]",
        "iupac|Ignore[BR]",
        "and you should ideally have a line for each molecule in the molecule database",
        "([TT]-f[tt] option). Missing molecules will be ignored."
    };
  
    t_filenm fnm[] = {
        { efDAT, "-f", "allmols",    ffREAD  },
        { efDAT, "-d", "gentop",     ffOPTRD },
        { efDAT, "-o", "tunedip",    ffWRITE },
        { efDAT, "-sel", "molselect",ffREAD },
        { efLOG, "-g", "charges",    ffWRITE },
        { efXVG, "-x", "dipcorr",    ffWRITE },
        { efXVG, "-qhisto", "q_histo",    ffWRITE },
        { efXVG, "-qdiff", "q_diff",    ffWRITE },
        { efXVG, "-mudiff", "mu_diff",    ffWRITE },
        { efXVG, "-thetadiff", "theta_diff",    ffWRITE },
        { efXVG, "-espdiff", "esp_diff",    ffWRITE },
        { efXVG, "-conv", "convergence", ffOPTWR }
    };
#define NFILE asize(fnm)
    static int  nrun=1,maxiter=100,reinit=0,seed=1993;
    static int  minimum_data=3,compress=1;
    static real tol=1e-3,stol=1e-6,watoms=1;
    static gmx_bool bRandom=FALSE,bZero=TRUE,bWeighted=TRUE,bOptHfac=FALSE,bQM=FALSE,bCharged=TRUE,bGaussianBug=TRUE,bPol=FALSE,bFitZeta=TRUE;
    static real J0_0=5,Chi0_0=1,w_0=5,step=0.01,hfac=0,rDecrZeta=-1;
    static real J0_1=30,Chi0_1=30,w_1=50,epsr=1;
    static real fc_mu=1,fc_bound=1,fc_quad=1,fc_charge=0,fc_esp=0;
    static real th_toler=170,ph_toler=5,dip_toler=0.5,quad_toler=5,q_toler=0.25;
    static char *opt_elem = NULL,*const_elem=NULL,*fixchi=(char *)"H";
    static char *lot = (char *)"B3LYP/aug-cc-pVTZ";
    static char *qgen[] = { NULL,(char *)"AXp", (char *)"AXs", (char *)"AXg", NULL };
    t_pargs pa[] = {
        { "-tol",   FALSE, etREAL, {&tol},
          "Tolerance for convergence in optimization" },
        { "-maxiter",FALSE, etINT, {&maxiter},
          "Max number of iterations for optimization" },
        { "-reinit", FALSE, etINT, {&reinit},
          "After this many iterations the search vectors are randomized again. A vlue of 0 means this is never done at all." },
        { "-stol",   FALSE, etREAL, {&stol},
          "If reinit is -1 then a reinit will be done as soon as the simplex size is below this treshold." },
        { "-nrun",   FALSE, etINT,  {&nrun},
          "This many runs will be done, before each run a complete randomization will be done" },
        { "-qm",     FALSE, etBOOL, {&bQM},
          "Use only quantum chemistry results (from the levels of theory below) in order to fit the parameters. If not set, experimental values will be used as reference with optional quantum chemistry results, in case no experimental results are available" },
        { "-lot",    FALSE, etSTR,  {&lot},
      "Use this method and level of theory when selecting coordinates and charges. Multiple levels can be specified which will be used in the order given, e.g.  B3LYP/aug-cc-pVTZ:HF/6-311G**" },
        { "-fc_bound",    FALSE, etREAL, {&fc_bound},
          "Force constant in the penalty function for going outside the borders given with the above six options." },
        { "-step",  FALSE, etREAL, {&step},
          "Step size in parameter optimization. Is used as a fraction of the starting value, should be less than 10%. At each reinit step the step size is updated." },
        { "-min_data",  FALSE, etINT, {&minimum_data},
          "Minimum number of data points in order to be able to optimize the parameters for a given atomtype" },
        { "-opt_elem",  FALSE, etSTR, {&opt_elem},
          "Space-separated list of elements to optimize, e.g. \"H C Br\". The other available elements in gentop.dat are left unmodified. If this variable is not set, all elements will be optimized." },
        { "-const_elem",  FALSE, etSTR, {&const_elem},
          "Space-separated list of elements to include but keep constant, e.g. \"O N\". These elements from gentop.dat are left unmodified" },
        { "-seed", FALSE, etINT, {&seed},
          "Random number seed for reinit" },
        { "-random", FALSE, etBOOL, {&bRandom},
          "Generate completely random starting parameters within the limits set by the options. This will be done at the very first step and before each subsequent run." },
        { "-weight", FALSE, etBOOL, {&bWeighted},
          "Perform a weighted fit, by using the errors in the dipoles presented in the input file. This may or may not improve convergence." },
        { "-th_toler", FALSE, etREAL, {&th_toler},
          "Minimum angle to be considered a linear A-B-C bond" },
        { "-ph_toler", FALSE, etREAL, {&ph_toler},
          "Maximum angle to be considered a planar A-B-C/B-C-D torsion" },
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML file" }
    };
    t_moldip  *md;
    FILE      *fp,*out;
    int       iModel;
    t_commrec *cr;
    output_env_t oenv;
    gmx_molselect_t gms;
    time_t    my_t;
    char      pukestr[STRLEN];

    cr = init_par(&argc,&argv);

    if (MASTER(cr))
        CopyRight(stdout,argv[0]);

    parse_common_args(&argc,argv,PCA_CAN_VIEW | (MASTER(cr) ? 0 : PCA_QUIET),
                      NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv);

    if (qgen[0]) 
        iModel = name2eemtype(qgen[0]);
    else
        iModel = eqgAXg;
    
    if (MASTER(cr)) {
        fp = ffopen(opt2fn("-g",NFILE,fnm),"w");

        time(&my_t);
        fprintf(fp,"# This file was created %s",ctime(&my_t));
        fprintf(fp,"# by the following command:\n# %s\n#\n",command_line());
        fprintf(fp,"# %s is part of G R O M A C S:\n#\n",ShortProgram());
        bromacs(pukestr,99);
        fprintf(fp,"# %s\n#\n",pukestr);
    }
    else
        fp = NULL;
    if (MASTER(cr)) 
        gms = gmx_molselect_init(opt2fn("-sel",NFILE,fnm));
    else
        gms = NULL;
    md = init_moldip(cr,bQM,bGaussianBug,iModel,rDecrZeta,epsr,
                     J0_0,Chi0_0,w_0,J0_1,Chi0_1,w_1,
                     fc_bound,fc_mu,fc_quad,fc_charge,
                     fc_esp,fixchi,bOptHfac,hfac,bPol,bFitZeta);
    read_moldip(md,fp ? fp : (debug ? debug : NULL),
                opt2fn("-f",NFILE,fnm),
                opt2fn_null("-d",NFILE,fnm),
                minimum_data,bZero,bWeighted,
                opt_elem,const_elem,
                lot,bCharged,oenv,gms,th_toler,ph_toler,dip_toler,
                watoms,FALSE);
    
    optimize_moldip(MASTER(cr) ? stderr : NULL,fp,opt2fn_null("-conv",NFILE,fnm),
                    md,maxiter,tol,nrun,reinit,step,seed,
                    bRandom,stol,oenv);
    if (MASTER(cr)) 
    {
        ffclose(fp);
  
        thanx(stdout);
    }
    
#ifdef GMX_MPI
    gmx_finalize();
#endif

    return 0;
}
