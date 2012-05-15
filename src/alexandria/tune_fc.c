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
#include <assert.h>
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

typedef struct {
    int natom;
    char **ai;
    int nparam;
    double *param;
    int noccur;
} opt_bad_t;

enum { ebadBOND, ebadANGLE, ebadDIH, ebadNR };

typedef struct {
    t_moldip *md;
    int nbad[ebadNR];
    opt_bad_t *bad[ebadNR];
    int nptot;
    double *ptot;
} opt_param_t;

static void add_obt(opt_bad_t *obt,int natom,char **ai,int nparam,double *param)
{
    int i;
    
    obt->natom = natom;
    snew(obt->ai,natom);
    for(i=0; (i<natom); i++) 
        obt->ai[i] = strdup(ai[i]);
    obt->nparam = nparam;
    snew(obt->param,nparam);
    for(i=0; (i<nparam); i++)
        obt->param[i] = param[i];
}

static void add_opt(opt_param_t *opt,int otype,int natom,char **ai,int nparam,double *param)
{
    assert(otype < ebadNR);
    assert(otype >= 0);
    srenew(opt->bad[otype],++opt->nbad[otype]);
    add_obt(&opt->bad[otype][opt->nbad[otype]-1],natom,ai,nparam,param);
}

static void opt2list(opt_param_t *opt)
{
    int i,j,p,n,nptot;
    
    for(i=n=0; (i<ebadNR); i++) {
        for(j=0; (j<opt->nbad[i]); j++) {
            nptot = n+opt->bad[i][j].nparam;
            if (nptot > opt->nptot) {
                srenew(opt->ptot,nptot);
                opt->nptot = nptot;
            }
            for(p=0; (p<opt->bad[i][j].nparam); p++) {
                opt->ptot[n++] = opt->bad[i][j].param[p];
            }
        }
    }
}

static void list2opt(opt_param_t *opt)
{
    gmx_poldata_t pd = opt->md->pd;
    int i,j,p,n,nptot;
    char buf[STRLEN];
    
    for(i=n=0; (i<ebadNR); i++) {
        for(j=0; (j<opt->nbad[i]); j++) {
            buf[0] = '\0';
            for(p=0; (p<opt->bad[i][j].nparam); p++) {
                opt->bad[i][j].param[p] = opt->ptot[n++];
                strcat(buf," ");
                strcat(buf,gmx_ftoa(opt->bad[i][j].param[p]));
            }
            switch (i) {
            case ebadBOND:
                gmx_poldata_set_gt_bond_params(pd,opt->bad[i][j].ai[0],
                                               opt->bad[i][j].ai[1],
                                               buf);
                break;
            case ebadANGLE:
                gmx_poldata_set_gt_angle_params(pd,opt->bad[i][j].ai[0],
                                                opt->bad[i][j].ai[1],
                                                opt->bad[i][j].ai[2],
                                                buf);
                break;
            case ebadDIH:
                gmx_poldata_set_gt_dihedral_params(pd,opt->bad[i][j].ai[0],
                                                   opt->bad[i][j].ai[1],
                                                   opt->bad[i][j].ai[2],
                                                   opt->bad[i][j].ai[3],
                                                   buf);
                break;
            }
        }
    }
}

static opt_param_t *init_opt(t_moldip *md,int *nparam)
{
    opt_param_t *opt;
    char   *ai[4];
    char   *params,**ptr;
    int    n,maxfc = 0;
    double *fc=NULL;
    
    snew(opt,1);
    opt->md = md;
    
    *nparam = 0;
    while(1 == gmx_poldata_get_gt_bond(md->pd,&(ai[0]),&(ai[1]),NULL,NULL,NULL,&params)) {
        ptr = split(' ',params);
        for(n = 0; (NULL != ptr[n]); n++) {
            if (n>=maxfc) {
                srenew(fc,++maxfc);
            }
            fc[n] = atof(ptr[n]);
        }
        add_opt(opt,ebadBOND,2,ai,n,fc);
        *nparam += n;
    }
    while(1 == gmx_poldata_get_gt_angle(md->pd,&(ai[0]),&(ai[1]),&(ai[2]),NULL,NULL,&params)) {
        ptr = split(' ',params);
        for(n = 0; (NULL != ptr[n]); n++) {
            if (n>maxfc) {
                srenew(fc,++maxfc);
            }
            fc[n] = atof(ptr[n]);
        }
        add_opt(opt,ebadANGLE,3,ai,n,fc);
        *nparam += n;
    }
    while (1 == gmx_poldata_get_gt_dihedral(md->pd,&(ai[0]),&(ai[1]),&(ai[2]),&(ai[3]),NULL,NULL,&params)) {
        ptr = split(' ',params);
        for(n = 0; (NULL != ptr[n]); n++) {
            if (n>maxfc) {
                srenew(fc,++maxfc);
            }
            fc[n] = atof(ptr[n]);
        }
        add_opt(opt,ebadDIH,4,ai,n,fc);
        *nparam += n;
    }
    sfree(fc);
    
    opt2list(opt);
    
    return opt;
}

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

static void analyze_idef(FILE *fp,int nmol,t_mymol mm[],gmx_poldata_t pd)
{
    int  gt,i,n,tp,ai,aj,ak,al;
    char *aai,*aaj,*aak,*aal,*params,**ptr;
    t_mymol *mymol;
    int  NGTB,NGTA,NGTD,nbtot,natot,ndtot;
    int  *ngtb=NULL,*ngta=NULL,*ngtd=NULL;
    char **cgtb=NULL,**cgta=NULL,**cgtd=NULL;
    
    NGTB = gmx_poldata_get_ngt_bond(pd);
    snew(ngtb,NGTB);
    snew(cgtb,NGTB);
    NGTA = gmx_poldata_get_ngt_angle(pd);
    snew(ngta,NGTA);
    snew(cgta,NGTA);
    NGTD = gmx_poldata_get_ngt_dihedral(pd);
    snew(ngtd,NGTD);
    snew(cgtd,NGTD);

    nbtot = natot = ndtot = 0;
    for(n=0; (n<nmol); n++) { 
        mymol = &(mm[n]);
        
        for(i=0; (i<mymol->ltop->idef.il[F_BONDS].nr); i+=interaction_function[F_BONDS].nratoms+1) {
            tp = mymol->ltop->idef.il[F_BONDS].iatoms[i];
            ai = mymol->ltop->idef.il[F_BONDS].iatoms[i+1];
            aj = mymol->ltop->idef.il[F_BONDS].iatoms[i+2];
            aai = *mymol->atoms->atomtype[ai];
            aaj = *mymol->atoms->atomtype[aj];
            if ((gt = gmx_poldata_search_gt_bond(pd,aai,aaj,NULL,NULL,&params)) != 0) {
                ngtb[gt-1]++;
                if (NULL == cgtb[gt-1]) {
                    snew(cgtb[gt-1],strlen(aai)+strlen(aaj)+2);
                    sprintf(cgtb[gt-1],"%s-%s",aai,aaj);
                }
                nbtot++;
            }
        }
        for(i=0; (i<mymol->ltop->idef.il[F_ANGLES].nr); i+=interaction_function[F_ANGLES].nratoms+1) {
            tp = mymol->ltop->idef.il[F_ANGLES].iatoms[i];
            ai = mymol->ltop->idef.il[F_ANGLES].iatoms[i+1];
            aj = mymol->ltop->idef.il[F_ANGLES].iatoms[i+2];
            ak = mymol->ltop->idef.il[F_ANGLES].iatoms[i+3];
            aai = *mymol->atoms->atomtype[ai];
            aaj = *mymol->atoms->atomtype[aj];
            aak = *mymol->atoms->atomtype[ak];
            if ((gt = gmx_poldata_search_gt_angle(pd,aai,aaj,aak,NULL,NULL,&params)) != 0) {
                ngta[gt-1]++;
                if (NULL == cgta[gt-1]) {
                    snew(cgta[gt-1],strlen(aai)+strlen(aaj)+strlen(aak)+3);
                    sprintf(cgta[gt-1],"%s-%s-%s",aai,aaj,aak);
                }
                natot++;
            }
        }
        for(i=0; (i<mymol->ltop->idef.il[F_PDIHS].nr); i+=interaction_function[F_PDIHS].nratoms+1) {
            tp = mymol->ltop->idef.il[F_PDIHS].iatoms[i];
            ai = mymol->ltop->idef.il[F_PDIHS].iatoms[i+1];
            aj = mymol->ltop->idef.il[F_PDIHS].iatoms[i+2];
            ak = mymol->ltop->idef.il[F_PDIHS].iatoms[i+3];
            al = mymol->ltop->idef.il[F_PDIHS].iatoms[i+4];
            aai = *mymol->atoms->atomtype[ai];
            aaj = *mymol->atoms->atomtype[aj];
            aak = *mymol->atoms->atomtype[ak];
            aal = *mymol->atoms->atomtype[al];
            if ((gt = gmx_poldata_search_gt_dihedral(pd,aai,aaj,aak,aal,
                                                     NULL,NULL,&params)) != 0) {
                ngtd[gt-1]++;
                if (NULL == cgtd[gt-1]) {
                    snew(cgtd[gt-1],strlen(aai)+strlen(aaj)+strlen(aak)+strlen(aal)+4);
                    sprintf(cgtd[gt-1],"%s-%s-%s-%s",aai,aaj,aak,aal);
                }
            }
            ndtot++;
        }
    }
    for(i=0; (i<NGTB); i++) {
        if (ngtb[i] > 0) {
            fprintf(fp,"BOND     %6d  %-20s  %d\n",i,cgtb[i],ngtb[i]);
            sfree(cgtb[i]);
        }
    }
    for(i=0; (i<NGTA); i++) {
        if (ngta[i] > 0) {
            fprintf(fp,"ANGLE    %6d  %-20s  %d\n",i,cgta[i],ngta[i]);
            sfree(cgta[i]);
        }
    }
    for(i=0; (i<NGTD); i++) {
        if (ngtd[i] > 0) {
            fprintf(fp,"DIHEDRAL %6d  %-20s  %d\n",i,cgtd[i],ngtd[i]);
            sfree(cgtd[i]);
        }
    }
    sfree(ngtd);
    sfree(ngta);
    sfree(ngtb);
    sfree(cgtd);
    sfree(cgta);
    sfree(cgtb);
    fprintf(fp,"In the total data set of %d molecules we have:\n",nmol);
    fprintf(fp,"%6d bonds     of %4d types\n",nbtot,NGTB);
    fprintf(fp,"%6d angles    of %4d types\n",natot,NGTA);
    fprintf(fp,"%6d dihedrals of %4d types.\n",ndtot,NGTD);
}

static void update_idef(t_mymol *mymol,gmx_poldata_t pd)
{
    int  gt,i,tp,ai,aj,ak,al;
    char *aai,*aaj,*aak,*aal,*params,**ptr;
    
    for(i=0; (i<mymol->ltop->idef.il[F_BONDS].nr); i+=interaction_function[F_BONDS].nratoms+1) {
        tp = mymol->ltop->idef.il[F_BONDS].iatoms[i];
        ai = mymol->ltop->idef.il[F_BONDS].iatoms[i+1];
        aj = mymol->ltop->idef.il[F_BONDS].iatoms[i+2];
        aai = *mymol->atoms->atomtype[ai];
        aaj = *mymol->atoms->atomtype[aj];
        if ((gt = gmx_poldata_search_gt_bond(pd,aai,aaj,NULL,NULL,&params)) != 0) {
            ptr = split(' ',params);
            if (NULL != ptr[0]) {
                mymol->mtop.ffparams.iparams[tp].harmonic.rA = atof(ptr[0]);
                sfree(ptr[0]);
            }
            if (NULL != ptr[1]) {
                mymol->mtop.ffparams.iparams[tp].harmonic.krA = atof(ptr[1]);
                sfree(ptr[1]);
            }
            sfree(ptr);
        }
        else {
            gmx_fatal(FARGS,"There are no parameters for bond %s-%s in the force field",aai,aaj);
        }
    }
    for(i=0; (i<mymol->ltop->idef.il[F_ANGLES].nr); i+=interaction_function[F_ANGLES].nratoms+1) {
        tp = mymol->ltop->idef.il[F_ANGLES].iatoms[i];
        ai = mymol->ltop->idef.il[F_ANGLES].iatoms[i+1];
        aj = mymol->ltop->idef.il[F_ANGLES].iatoms[i+2];
        ak = mymol->ltop->idef.il[F_ANGLES].iatoms[i+3];
        aai = *mymol->atoms->atomtype[ai];
        aaj = *mymol->atoms->atomtype[aj];
        aak = *mymol->atoms->atomtype[ak];
        if ((gt = gmx_poldata_search_gt_angle(pd,aai,aaj,aak,NULL,NULL,&params)) != 0) {
            ptr = split(' ',params);
            if (NULL != ptr[0]) {
                mymol->mtop.ffparams.iparams[tp].harmonic.rA = atof(ptr[0]);
                sfree(ptr[0]);
            }
            if (NULL != ptr[1]) {
                mymol->mtop.ffparams.iparams[tp].harmonic.krA = atof(ptr[1]);
                sfree(ptr[1]);
            }
            sfree(ptr);
        }
        else {
            gmx_fatal(FARGS,"There are no parameters for angle %s-%s-%s in the force field",aai,aaj,aak);
        }
    }
    for(i=0; (i<mymol->ltop->idef.il[F_PDIHS].nr); i+=interaction_function[F_PDIHS].nratoms+1) {
        tp = mymol->ltop->idef.il[F_PDIHS].iatoms[i];
        ai = mymol->ltop->idef.il[F_PDIHS].iatoms[i+1];
        aj = mymol->ltop->idef.il[F_PDIHS].iatoms[i+2];
        ak = mymol->ltop->idef.il[F_PDIHS].iatoms[i+3];
        al = mymol->ltop->idef.il[F_PDIHS].iatoms[i+4];
        aai = *mymol->atoms->atomtype[ai];
        aaj = *mymol->atoms->atomtype[aj];
        aak = *mymol->atoms->atomtype[ak];
        aal = *mymol->atoms->atomtype[al];
        if ((gt = gmx_poldata_search_gt_dihedral(pd,aai,aaj,aak,aal,
                                                 NULL,NULL,&params)) != 0) {
            ptr = split(' ',params);
            if (NULL != ptr[0]) {
                mymol->mtop.ffparams.iparams[tp].pdihs.phiA = atof(ptr[0]);
                sfree(ptr[0]);
            }
            if (NULL != ptr[1]) {
                mymol->mtop.ffparams.iparams[tp].pdihs.cpA = atof(ptr[1]);
                sfree(ptr[1]);
            }
            if (NULL != ptr[2]) {
                mymol->mtop.ffparams.iparams[tp].pdihs.mult = atof(ptr[2]);
                sfree(ptr[2]);
            }
            sfree(ptr);
        }
        else {
            gmx_fatal(FARGS,"There are no parameters for angle %s-%s-%s in the force field",aai,aaj,aak);
        }
    }
}

static double calc_opt_deviation(opt_param_t *opt)
{
    int    i,j,count,atomnr;
    double qq,qtot,rr2,ener[ermsNR],etot[ermsNR];
    real   t = 0;
    rvec   mu_tot = {0,0,0};
    tensor force_vir={{0,0,0},{0,0,0},{0,0,0}};
    t_nrnb   my_nrnb;
    gmx_wallcycle_t wcycle;
    gmx_bool     bConverged;
    t_mymol *mymol;
    int     eQ;
    gmx_mtop_atomloop_all_t aloop;
    gmx_enerdata_t gmd;
    t_atom *atom; 
    int    at_global,resnr;
    
    if (PAR(opt->md->cr)) 
    {
        gmx_bcast(sizeof(opt->md->bDone),&opt->md->bDone,opt->md->cr);
        gmx_bcast(sizeof(opt->md->bFinal),&opt->md->bFinal,opt->md->cr);
    }
    if (opt->md->bDone)
        return opt->md->ener[ermsTOT];
    if (PAR(opt->md->cr)) 
    {
        gmx_poldata_comm_eemprops(opt->md->pd,opt->md->cr);
    }
    init_nrnb(&my_nrnb);
  
    wcycle  = wallcycle_init(stdout,0,opt->md->cr,1);
    for(j=0; (j<ermsNR); j++)
    {
        etot[j] = 0;
        opt->md->ener[j] = 0;
    }
    for(i=0; (i<opt->md->nmol); i++) {
        mymol = &(opt->md->mymol[i]);
        if ((mymol->eSupport == eSupportLocal) ||
            (opt->md->bFinal && (mymol->eSupport == eSupportRemote)))
        {
            /* Update topology for this molecule */
            update_idef(mymol,opt->md->pd);
            
            /* Reset energies */
            for(j=0; (j<ermsNR); j++)
                ener[j] = 0;
            
            /* Now compute energy */
            atoms2md(&mymol->mtop,&(mymol->ir),0,NULL,0,
                     mymol->mtop.natoms,mymol->md);
            
            /* Now optimize the shell positions */
            if (mymol->shell) { 
                count = 
                    relax_shell_flexcon(debug,opt->md->cr,FALSE,0,
                                        &(mymol->ir),TRUE,
                                        GMX_FORCE_ALLFORCES,FALSE,
                                        mymol->ltop,NULL,NULL,&(mymol->enerd),
                                        NULL,&(mymol->state),
                                        mymol->f,force_vir,mymol->md,
                                        &my_nrnb,wcycle,NULL,
                                        &(mymol->mtop.groups),
                                        mymol->shell,mymol->fr,FALSE,t,mu_tot,
                                        mymol->mtop.natoms,&bConverged,NULL,NULL);
            }
            else {
                do_force(debug,opt->md->cr,&(mymol->ir),0,
                         &my_nrnb,wcycle,mymol->ltop,
                         &mymol->mtop,&(mymol->mtop.groups),
                         mymol->box,mymol->x,NULL,/*history_t *hist,*/
                         mymol->f,force_vir,mymol->md,
                         &mymol->enerd,NULL,/*mymol->fcd,*/
                         0.0,NULL,/*mymol->graph,*/
                         mymol->fr,
                         NULL,mu_tot,t,NULL,NULL,FALSE,GMX_FORCE_ALLFORCES | GMX_FORCE_STATECHANGED);

            }
            ener[ermsEPOT] = 1e-6*sqr(mymol->enerd.term[F_EPOT]-mymol->Hform);
            for(j=0; (j<ermsNR); j++)
                etot[j] += ener[j];
        }
    }
    for(j=0; (j<ermsTOT); j++) {
        opt->md->ener[j]       += opt->md->fc[j]*etot[j]/opt->md->nmol_support;
        opt->md->ener[ermsTOT] += opt->md->ener[j];
    }
    if (debug)
    {
        fprintf(debug,"ENER:");
        for(j=0; (j<ermsNR); j++)
            fprintf(debug,"  %8.3f",etot[j]);
        fprintf(debug,"\n");
    }
    /* Global sum energies */
    if (PAR(opt->md->cr)) 
    {
#ifdef GMX_DOUBLE
        gmx_sumd(ermsNR,opt->md->ener,opt->md->cr);
#else
        gmx_sumf(ermsNR,opt->md->ener,opt->md->cr);
#endif
    }
    return opt->md->ener[ermsTOT];
}

static double energy_function(void *params,double v[])
{
    opt_param_t *opt = (opt_param_t *)params;
    t_moldip *md = opt->md;
    int      i;
    
    /* Copy parameters to topologies */
    for(i=0; (i<opt->nptot); i++) {
        opt->ptot[i] = v[i];
    }
    list2opt(opt);
  
    return calc_opt_deviation(opt);
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

static void guess_all_param(FILE *fplog,opt_param_t *opt,int iter,real stepsize,
                            gmx_bool bRandom,gmx_rng_t rng,
                            double orig_param[],double test_param[])
{
    double   ppp,xxx;
    gmx_bool bStart = (iter == 0);
    gmx_bool bRand  = bRandom && (iter == 0);
    int      n;
    
    for(n=0; (n<opt->nptot); n++) 
    {
        if (bStart)
        {
            ppp = opt->ptot[n];
            xxx = guess_new_param(ppp,stepsize,0.8*ppp,1.25*ppp,rng,bRand);
            if (bRand)
                orig_param[n] = xxx;
            else
                orig_param[n] = ppp;
            ppp = xxx;
        }
        else 
            ppp = guess_new_param(orig_param[n],stepsize,0.8*ppp,1.25*ppp,rng,bRand);
        test_param[n] = ppp;
    }
}


static void optimize_moldip(FILE *fp,FILE *fplog,
                            t_moldip *md,int maxiter,real tol,
                            int nrun,int reinit,real stepsize,int seed,
                            gmx_bool bRandom,real stol,output_env_t oenv)
{
    double chi2,chi2_min,wj,rms_nw;
    int    status = 0;
    int    i,k,index,n,nparam;
    double *orig_param,*best_param,*start;
    gmx_bool bMinimum=FALSE;
    char   *name,*qstr,*rowstr;
    char   buf[STRLEN];
    gmx_rng_t rng;
    opt_param_t *opt;
    
    opt = init_opt(md,&nparam);
    if (MASTER(md->cr)) 
    {
        rng = gmx_rng_init(seed);
  
        snew(orig_param,nparam+1);
        snew(best_param,nparam+1);
        
        /* Starting point */
        snew(start,nparam+1);
        
        chi2_min = GMX_REAL_MAX; 
        for(n=0; (n<nrun); n++) {
            if ((NULL != fp) && (0 == n)) {
                fprintf(fp,"\nStarting run %d out of %d\n",n+1,nrun);
            }
                    
            guess_all_param(fplog,opt,0,stepsize,bRandom,rng,orig_param,start);
            
            nmsimplex(NULL,(void *)opt,energy_function,start,nparam,
                      tol,1.0,maxiter,&chi2);
                
            if (chi2 < chi2_min) {
                bMinimum = TRUE;
                /* Print convergence if needed */
                for(k=0; (k<opt->nptot); k++)
                {
                    best_param[k] = start[k];
                }
                chi2_min   = chi2;
            }
            
            if (fp) 
                fprintf(fp,"%5d %8.3f  %8.3f\n",n,1000*chi2,1000*sqrt(md->ener[ermsTOT]));
        }
        
        if (bMinimum) {
            for(k=0; (k<nparam); k++)
                start[k] = best_param[k];
            
            energy_function(opt,best_param);
            if (fplog)
            {
                fprintf(fplog,"\nMinimum chi^2 value during optimization: %.3f.\n",
                        1000*chi2_min);
                fprintf(fplog,"\nMinimum RMSD value during optimization: %.3f.\n",
                        1000*opt->md->ener[ermsTOT]);
                for(k=0; (k<nparam); k++)
                    fprintf(fplog,"Param %5d  = %g\n",k,best_param[k]);
            }
        }
        md->bDone = TRUE;
        calc_opt_deviation(opt);
        
    }
    else 
    {
        /* Slave calculators */
        do {
            calc_opt_deviation(opt);
        } while (!md->bDone);
    }
}

static real quality_of_fit(real chi2,int N)
{
    return -1;
}

static void print_moldip_specs(FILE *fp,t_moldip *md)
{
    int i;
    
    fprintf(fp,"Nr.   %-40s  %10s  %10s\n","Molecule","DHf@298K","Emol@0K");
    for(i=0; (i<md->nmol); i++) {
        fprintf(fp,"%-5d %-40s  %10g  %10g\n",i,md->mymol[i].molname,
                md->mymol[i].Hform,md->mymol[i].Emol);
    }
    fprintf(fp,"\n");
    
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
        { efDAT, "-o", "tune_fc",    ffWRITE },
        { efDAT, "-sel", "molselect",ffREAD },
        { efLOG, "-g", "tune_fc",    ffWRITE }
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
    static int  nthreads=0; /* set to determine # of threads automatically */
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
#ifdef GMX_THREAD_MPI
    { "-nt",      FALSE, etINT, {&nthreads},
      "Number of threads to start (0 is guess)" },
#endif
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

#ifndef GMX_THREAD_MPI
    nthreads=1;
#endif
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
        gms = gmx_molselect_init(opt2fn_null("-sel",NFILE,fnm));
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
                TRUE,TRUE,TRUE,2,watoms,FALSE);
    print_moldip_specs(fp,md);
    analyze_idef(fp,md->nmol,md->mymol,md->pd);
    exit(1);
    optimize_moldip(MASTER(cr) ? stderr : NULL,fp,
                    md,maxiter,tol,nrun,reinit,step,seed,
                    bRandom,stol,oenv);
    if (MASTER(cr)) 
    {
        ffclose(fp);
  
        thanx(stdout);
    }
    
#ifdef GMX_MPI
    if (gmx_parallel_env_initialized())
        gmx_finalize();
#endif

    return 0;
}
