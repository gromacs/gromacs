/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: molprop_tables.c,v 1.21 2009/05/29 15:01:18 spoel Exp $
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
#include <string.h>
#include <math.h>
#include "maths.h"
#include "futil.h"
#include "smalloc.h"
#include "string2.h"
#include "vec.h"
#include "xvgr.h"
#include "statutil.h"
#include "copyrite.h"
#include "gmx_statistics.h"
#include "gmx_fatal.h"
#include "gromacs/linearalgebra/matrix.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "molselect.h"
#include "molprop.hpp"
#include "molprop_util.hpp"
#include "molprop_tables.hpp"

static const char *SPOEL = "spoel";

typedef struct 
{
    char *cat;
    gmx_bool bPrint;
    int count;
    real rms;
} t_catli;

static int comp_catli(const void *a,const void *b)
{
    t_catli *ca = (t_catli *)a;
    t_catli *cb = (t_catli *)b;
    
    return strcmp(ca->cat,cb->cat);
}

typedef struct 
{
    int ncat;
    t_catli *catli;
} t_cats;

static void add_cat(t_cats *cats,const char *catname)
{
    int i;
  
    for(i=0; (i<cats->ncat); i++) 
        if (strcasecmp(catname,cats->catli[i].cat) == 0) 
        {
            cats->catli[i].count++;
            break;
        }
    if (i == cats->ncat) 
    {
        cats->ncat++;
        srenew(cats->catli,cats->ncat);
        cats->catli[cats->ncat-1].cat = strdup(catname);
        cats->catli[cats->ncat-1].count = 1;
        cats->catli[cats->ncat-1].rms = 0;
    }
}

static void decrease_table_number(FILE *fp)
{
    fprintf(fp,"\\addtocounter{table}{-1}\n");
}

static void stats_footer(FILE *fp,gmx_bool bSideways)
{
    fprintf(fp,"\\hline\n\\end{tabular}\n\\end{%stable}\n",
            (bSideways) ? "sideways" : "");
}

static void stats_header(FILE *fp,int emp,int caption,gmx_bool bSideways,
                         t_qmcount *qmc,int ims)
{
    char unit[32];
    int  i;
    
    if (!caption)
        decrease_table_number(fp);
    
    fprintf(fp,"\\begin{%stable}[H]\n\\centering\n",bSideways ? "sideways" : "");
    if (caption) 
    {
        switch (emp) 
        {
        case empPOLARIZABILITY:
            sprintf(unit,"\\AA$^3$");
            break;
        case empDIPOLE:
            sprintf(unit,"D");
            break;
        case empQUADRUPOLE:
            sprintf(unit,"B");
            break;
        case empENERGY:
            sprintf(unit,"kJ/mol");
            break;
        default:
            gmx_fatal(FARGS,"Unknown property %s",emp_name[emp]);
        }
        fprintf(fp,"\\caption{Performance of the different methods for predicting the molecular %s for molecules containing different chemical groups, given as the RMSD from experimental values (%s), and in brackets the number of molecules in this particular subset. {\\bf Data set: %s.} At the bottom the correlation coefficient R, the regression coefficient a and the intercept b are given as well as the normalized quality of the fit $\\chi^2$.}\n\\label{%s_rmsd}\n",
                emp_name[emp],unit,ims_names[ims],emp_name[emp]);
    }
    else 
    {
        fprintf(fp,"\\caption{Continued.}\n");
    }
    fprintf(fp,"\\begin{tabular}{l");
    for(i=0; (i<qmc->n); i++)
    {
        fprintf(fp,"c");
    }
    fprintf(fp,"}\n\\hline\n");
    
    fprintf(fp,"Method ");

    for(i=0; (i<qmc->n); i++) 
    {
        fprintf(fp," & %s",qmc->method[i]);
    }
    fprintf(fp,"\\\\\n");
  
    for(i=0; (i<qmc->n); i++) 
        fprintf(fp,"& %s ",qmc->basis[i]);
    fprintf(fp,"\\\\\n");
  
    for(i=0; (i<qmc->n); i++) 
        fprintf(fp,"& %s ",qmc->type[i]);
    fprintf(fp,"\\\\\n\\hline\n");
}

void gmx_molprop_stats_table(FILE *fp,int emp,
                             std::vector<alexandria::MolProp> mp,
                             int ntot,
                             t_qmcount *qmc,int iQM,char *lot,
                             double outlier,gmx_molselect_t gms,int ims)
{
    std::vector<alexandria::MolProp>::iterator mpi;
    std::vector<std::string>::iterator si;
    int    i,line,k,m,N,nprint,nqmres,nexpres;
    double exp_val,qm_val;
    real   rms,R,a,da,b,db,chi2;
    char   lbuf[256],tmp[32],catbuf[STRLEN];
    t_cats *cats;
    gmx_stats_t lsq,*lsqtot;
    gmx_bool   bSideways = (qmc->n > 3);

    nprint = 0;
    snew(cats,1);
    for(mpi=mp.begin(); (mpi<mp.end()); mpi++) 
    {
        if ((ims == gmx_molselect_status(gms,mpi->GetIupac().c_str())) &&
            mpi->HasComposition(SPOEL))
        {
            nprint++;
            for(si=mpi->BeginCategory(); (si<mpi->EndCategory()); si++) 
            {
                add_cat(cats,si->c_str());
            }
        }
    }
    if (nprint <= 0)
    {
        sfree(cats);
        return;
    }
    
    stats_header(fp,emp,1,bSideways,qmc,ims);
    
    qsort(cats->catli,cats->ncat,sizeof(cats->catli[0]),comp_catli);
    line = 0;
    for(i=0; (i<cats->ncat); i++) 
    {
        sprintf(catbuf," %s (%d) ",cats->catli[i].cat,cats->catli[i].count);
        nqmres = 0;
        nexpres = 0;
        for(k=0; (k<qmc->n); k++) 
        {
            sprintf(lbuf,"%s/%s",qmc->method[k],qmc->basis[k]);
            lsq = gmx_stats_init();
            for(mpi=mp.begin(); (mpi<mp.end()); mpi++) 
            {
                if ((gmx_molselect_status(gms,mpi->GetIupac().c_str()) == ims) &&
                    (mpi->HasComposition(SPOEL)))
                {
                    if (mpi->SearchCategory(cats->catli[i].cat) == 1)
                    {
                        /*for(m=0; (m<qmc->nconf); m++) 
                          {*/
                            if (mp_get_prop(*mpi,emp,iqmExp,NULL,
                                            NULL/*qmc->conf[m]*/,NULL,&exp_val) > 0)
                            {
                                nexpres = 1;
                                if (mp_get_prop(*mpi,emp,iqmQM,lbuf,NULL/*qmc->conf[m]*/,
                                                qmc->type[k],&qm_val) > 0)
                                {
                                    if (debug)
                                        fprintf(debug,"%s %s k = %d - TAB4\n",
                                                mpi->GetMolname().c_str(),cats->catli[i].cat,k);
                                    gmx_stats_add_point(lsq,exp_val,qm_val,0,0);
                                }
                            }
                            /*}*/
                    }
                }
            }
            if (outlier > 0)
                gmx_stats_remove_outliers(lsq,outlier);
            if (gmx_stats_get_rmsd(lsq,&rms) == estatsOK)
            { 
                sprintf(tmp,"& %8.2f",rms);
                strncat(catbuf,tmp,STRLEN-1);
                nqmres++;
            }
            else
            {
                strncat(catbuf,"& -",STRLEN-1);
            }
            
            if (gmx_stats_done(lsq) == estatsOK)
                sfree(lsq);
        }
        cats->catli[i].bPrint = ((nqmres > 0) && (nexpres > 0));
        if (cats->catli[i].bPrint)
        {
            strncat(catbuf,"\\\\\n",STRLEN-1);
            fprintf(fp,"%s",catbuf);
            line++;
        }
        if (line == 20)
        {
            stats_footer(fp,bSideways);
            stats_header(fp,emp,0,bSideways,qmc,ims);
            line = 0;
        }
    }
    fprintf(fp,"All");
    snew(lsqtot,qmc->n);
    for(k=0; (k<qmc->n); k++) 
    {
        sprintf(lbuf,"%s/%s",qmc->method[k],qmc->basis[k]);
        lsqtot[k] = gmx_stats_init();
        for(mpi=mp.begin(); (mpi<mp.end()); mpi++) 
        {
            if ((gmx_molselect_status(gms,mpi->GetIupac().c_str()) == ims) &&
                (mpi->HasComposition(SPOEL)))
            {
                for(m=0; (m<qmc->nconf); m++) {
                    if ((mp_get_prop(*mpi,emp,iqmExp,NULL,qmc->conf[m],NULL,&exp_val) > 0) &&
                        (mp_get_prop(*mpi,emp,iqmQM,lbuf,qmc->conf[m],qmc->type[k],&qm_val) > 0)) 
                    {
                        gmx_stats_add_point(lsqtot[k],exp_val,qm_val,0,0);
                    }
                }
            }
        }
        if (gmx_stats_get_rmsd(lsqtot[k],&rms) == estatsOK)
            fprintf(fp,"& %8.2f",rms);
        else
            fprintf(fp,"& -");
    }
    fprintf(fp,"\\\\\n");
    fprintf(fp,"\\hline\n");
    fprintf(fp,"N");
    for(k=0; (k<qmc->n); k++) 
    {
        /* gmx_stats_remove_outliers(lsqtot[k],3); */
        gmx_stats_get_npoints(lsqtot[k],&N);
        fprintf(fp,"& %d",N);
    }
    fprintf(fp,"\\\\\n");
    fprintf(fp,"a");
    for(k=0; (k<qmc->n); k++) 
    {
        if (gmx_stats_get_ab(lsqtot[k],elsqWEIGHT_NONE,&a,&b,&da,&db,&chi2,&R) == 
            estatsOK) 
            fprintf(fp,"& %8.2f(%4.2f)",a,da);
        
        else
            fprintf(fp,"& -");
    }
    fprintf(fp,"\\\\\n");
    fprintf(fp,"b");
    for(k=0; (k<qmc->n); k++) 
    {
        if (gmx_stats_get_ab(lsqtot[k],elsqWEIGHT_NONE,&a,&b,&da,&db,&chi2,&R) ==
            estatsOK)
            fprintf(fp,"& %8.2f(%4.2f)",b,db);
        else
            fprintf(fp,"& -");
    }
    fprintf(fp,"\\\\\n");
    fprintf(fp,"R (\\%%)");
    for(k=0; (k<qmc->n); k++) 
    {
        if (gmx_stats_get_corr_coeff(lsqtot[k],&R) == estatsOK)
            fprintf(fp,"& %8.2f",100*R);
        else
            fprintf(fp,"& -");
    }
    fprintf(fp,"\\\\\n");
    fprintf(fp,"$\\chi^2$");
    for(k=0; (k<qmc->n); k++) 
    {
        if (gmx_stats_get_ab(lsqtot[k],elsqWEIGHT_NONE,&a,&b,&da,&db,&chi2,&R) ==
            estatsOK)
            fprintf(fp,"& %8.2f",chi2);
        else
            fprintf(fp,"& -");
    }
    fprintf(fp,"\\\\\n");
    stats_footer(fp,bSideways);
    
    for(k=0; (k<qmc->n); k++) 
    {
        if (gmx_stats_done(lsqtot[k]) == estatsOK)
            sfree(lsqtot[k]);
    }
    sfree(lsqtot);
}
    
static void composition_header(FILE *fp,int caption,int ims)
{
    if (!caption)
        decrease_table_number(fp);
    fprintf(fp,"\\begin{sidewaystable}[H]\n\\centering\n");
    if (caption)
        fprintf(fp,"\\caption{Decomposition of molecules into Alexandria atom types. {\\bf Data set: %s.} Charge is given when not zero, multiplicity is given when not 1.}\n\\label{frag_defs}\n",ims_names[ims]);
    else 
        fprintf(fp,"\\caption{Continued}\n");
    fprintf(fp,"\\begin{tabular}{lll}\n\\hline\n");
    fprintf(fp,"Molecule & Formula  & Types \\\\\n");
    fprintf(fp,"\\hline\n");
}

void gmx_molprop_composition_table(FILE *fp,std::vector<alexandria::MolProp> mp,
                                   gmx_molselect_t gms,int ims)
{
    std::vector<alexandria::MolProp>::iterator mpi;
    alexandria::MolecularCompositionIterator mci;
    
    int k,q,m,iline,nprint;
    char qbuf[32];
    const char *comps[2] = { "spoel", "miller" };
    int  kmax = 1;
    
    nprint = 0;
    for(mpi=mp.begin(); (mpi<mp.end()); mpi++)
    {
        if ((ims == gmx_molselect_status(gms,mpi->GetIupac().c_str())) &&
            (mpi->HasComposition(SPOEL)))
            nprint++;
    }
    if (nprint <= 0)
        return;

    composition_header(fp,1,ims);  
    iline = 0;
    for(mpi=mp.begin(); (mpi<mp.end()); mpi++)
    {
        if ((ims == gmx_molselect_status(gms,mpi->GetIupac().c_str())) &&
            (mpi->HasComposition(SPOEL)))
        {
            q = mpi->GetCharge();
            m = mpi->GetMultiplicity();
            if ((q != 0) || (m != 1))
            {
                if ((q != 0) && (m != 1)) 
                    sprintf(qbuf," (q=%c%d, mult=%d)",(q < 0) ? '-' : '+',abs(q),m);
                else if (q != 0)
                    sprintf(qbuf," (q=%c%d)",(q < 0) ? '-' : '+',abs(q));
                else
                    sprintf(qbuf," (mult=%d)",m);
            }
            else
            {
                qbuf[0] = '\0';
            }
            fprintf(fp,"%3d. %s%s & %s &",++iline,
                    mpi->GetIupac().c_str(),qbuf,
                    mpi->GetFormula().c_str());
            for(k=0; (k<kmax); k++) 
            {
                mci = mpi->SearchMolecularComposition(comps[k]);
                if (mci != mpi->EndMolecularComposition())
                {
                    alexandria::AtomNumIterator ani;
                    
                    for(ani=mci->BeginAtomNum(); (ani<mci->EndAtomNum()); ani++)
                    {
                        fprintf(fp,"%d%s\t",ani->GetNumber(),ani->GetAtom().c_str());
                    }
                }
                if (k == kmax-1)
                    fprintf(fp," \\\\\n");
                else
                    fprintf(fp," &");
            }
            if (((iline % 28) == 0) && (mpi<mp.end()-1)) 
            {
                fprintf(fp,"\\hline\n\\end{tabular}\n\\end{sidewaystable}\n");
                composition_header(fp,0,ims);
            }
        }
    }
    fprintf(fp,"\\hline\n\\end{tabular}\n\\end{sidewaystable}\n");
}

typedef struct {
    char *cat;
    int nmol;
    char **molec;
} t_cat_stat;

static int comp_cats(const void *a,const void *b)
{
    t_cat_stat *ca = (t_cat_stat *)a;
    t_cat_stat *cb = (t_cat_stat *)b;
    
    return strcmp(ca->cat,cb->cat);
}

static void add_cats(int *ncs,t_cat_stat **cs,const char *iupac,const char *mycat)
{
    int j,nmol;

    for(j=0; (j<*ncs); j++)
    {
        if (strcmp((*cs)[j].cat,mycat) == 0)
            break;
    }
    if (j == *ncs) 
    {
        srenew(*cs,++(*ncs));
        (*cs)[j].cat = strdup(mycat);
        (*cs)[j].nmol = 0;
        (*cs)[j].molec = NULL;
    }
    nmol = ++(*cs)[j].nmol;
    srenew((*cs)[j].molec,nmol);
    (*cs)[j].molec[nmol-1] = strdup(iupac);
}

static void category_header(FILE *fp,int start,int ims)
{
    fprintf(fp,"\\newpage\n");
    if (start == 0)
        decrease_table_number(fp);
    fprintf(fp,"\\begin{table}[H]\n");
    fprintf(fp,"\\caption{Molecules that are part of each category used for statistics%s.}\n",
            (start == 1) ? "" : " (continued)");
    fprintf(fp,"\\begin{tabularx}{16cm}{lcX}\n");
    fprintf(fp,"\\hline\n");
    fprintf(fp,"Category & N & Molecule(s)\\\\\n");
    fprintf(fp,"\\hline\n");
}

static void category_footer(FILE *fp)
{
    fprintf(fp,"\\hline\n");
    fprintf(fp,"\\end{tabularx}\n");
    fprintf(fp,"\\end{table}\n");
}

void gmx_molprop_category_table(FILE *fp,
                                std::vector<alexandria::MolProp> mp,
                                gmx_molselect_t gms,int ims)
{
    std::vector<alexandria::MolProp>::iterator mpi;
    std::vector<std::string>::iterator si;
    int i,j,iline,ncs,prmax=50;
    t_cat_stat *cs=NULL;
    const char *iupac;
    
    ncs = 0;
    for(mpi=mp.begin(); (mpi<mp.end()); mpi++)
    {
        iupac = mpi->GetIupac().c_str();
        if (ims == gmx_molselect_status(gms,iupac))
        {
            for(si=mpi->BeginCategory(); (si<mpi->EndCategory()); si++)
            {
                add_cats(&ncs,&cs,iupac,si->c_str());
            }
        }
    }
    
    if (ncs > 0) 
    {
        qsort(cs,ncs,sizeof(cs[0]),comp_cats);
        category_header(fp,1,ims);  
        iline = 0;
        for(i=0; (i<ncs); i++) 
        {
            fprintf(fp,"%s & %d &",cs[i].cat,cs[i].nmol);
            for(j=0; (j<cs[i].nmol); j++) 
            {
                fprintf(fp," %s",cs[i].molec[j]);
                iline++;
                if (j<cs[i].nmol-1)
                {
                    if (iline >= prmax)
                    {
                        fprintf(fp,"\\\\\n");
                        category_footer(fp);
                        category_header(fp,0,ims);
                        fprintf(fp,"%s & (ctd) &",cs[i].cat);
                        iline = 0;
                    }
                    else 
                        fprintf(fp,",");
                }
            }
            if (iline != 0)
                fprintf(fp,"\\\\\n");
        }
        category_footer(fp);
    }
}

static void atomtype_tab_header(FILE *fp,int caption)
{
    if (!caption)
        decrease_table_number(fp);
    fprintf(fp,"\\begin{table}[H]\n\\centering\n");
    if (caption)
        fprintf(fp,"\\caption{Definition of atom types. The number of occurences  of atom N$_{Exp}$ and N$_{QM}$ indicate how many of the polarizability values were fitted to experimental data or quantum chemistry, respectively. The columns Ahc and Ahp contain group polarizabilites computed using Miller's equation~\\protect\\cite{Miller1979a} and parameters~\\protect\\cite{Miller1990a} and Kang and Jhon's method~\\cite{Kang1982a} with the parametrization of Miller~\\protect\\cite{Miller1990a} respectively. The column BS contains the equivalent using the polarizabilities of Bosque and Sales~\\protect\\cite{Bosque2002a}.}\n\\label{fragments}\n");
    else
        fprintf(fp,"\\caption{Continued}");
    fprintf(fp,"\\begin{tabular}{cccccccc}\n\\hline\n");
    fprintf(fp,"Name  & N$_{Exp}$ & N$_{QM}$ & \\multicolumn{4}{c}{Polarizability} & Ref.\\\\\n");
    fprintf(fp,"& & & Ahc & Ahp & BS & AX ($\\sigma_{AX}$) & \\\\\n");
    fprintf(fp,"\\hline\n");
}

static void atomtype_tab_end(FILE *fp)
{
    fprintf(fp,"\\hline\n");
    fprintf(fp,"\\end{tabular}\n\\end{table}\n\n");
}

typedef struct {
    char *gt_type,*elem,*miller_equiv,*spref;
    gmx_stats_t lsq,nexp,nqm;
} t_sm_lsq;
    
static void gmx_molprop_atomtype_polar_table(FILE *fp,int npd,gmx_poldata_t pd[],
                                             gmx_poldata_t pd_aver,
                                             std::vector<alexandria::MolProp> mp,
                                             int iQM,char *lot,output_env_t oenv,const char *histo)
{
    std::vector<alexandria::MolProp>::iterator mpi;
    FILE   *xvg;
    int    i,j,pp,ntab,nfitexp,nfitqm,atomnumber;
    int    nfirst = 20,nsecond = 24;
    double ahc,ahp,bos_pol,alexandria_pol,sig_pol;
    double exp_val,qm_val;
    real   alexandria_aver,alexandria_sigma,nnn;
    char   *gt_type[2] = { NULL, NULL };
    char   *elem,*miller_equiv,group[32],*spref;
    char   *desc;
    char   *charge;
    int    nsmlsq=0,estats,N,nset;
    t_sm_lsq *smlsq=NULL;
    int    cur = 0,emp=empPOLARIZABILITY;
#define prev (1-cur)

    /* First gather statistics from different input files. Note that the input files
     * do not need to have the same sets of types, as we check for the type name.
     */
    for(pp = 0; (pp<npd); pp++) 
    {
        while(1 == gmx_poldata_get_atype(pd[pp],&elem,&desc,&(gt_type[cur]),
                                         &miller_equiv,&charge,
                                         &alexandria_pol,&sig_pol,&spref))
        {
            if (((NULL == gt_type[prev]) || (strcmp(gt_type[cur],gt_type[prev]) != 0)) &&
                (alexandria_pol > 0))
            {
                for(j=0; (j<nsmlsq); j++)
                    if (strcasecmp(gt_type[cur],smlsq[j].gt_type) == 0)
                        break;
                if (j == nsmlsq) 
                {
                    srenew(smlsq,++nsmlsq);
                    smlsq[j].gt_type = strdup(gt_type[cur]);
                    smlsq[j].lsq = gmx_stats_init();
                    smlsq[j].nexp = gmx_stats_init();
                    smlsq[j].nqm = gmx_stats_init();
                    smlsq[j].elem = strdup(elem);
                    smlsq[j].miller_equiv = strdup(miller_equiv);
                    smlsq[j].spref = strdup(spref);
                }
                if ((estats = gmx_stats_get_npoints(smlsq[j].lsq,&N)) != estatsOK)
                    gmx_fatal(FARGS,"Statistics problems: %s",gmx_stats_message(estats));
                if ((estats = gmx_stats_add_point(smlsq[j].lsq,N,alexandria_pol,0,0)) != estatsOK)
                    gmx_fatal(FARGS,"Statistics problems: %s",gmx_stats_message(estats));
                nfitexp = nfitqm = 0;
                for(mpi=mp.begin(); (mpi<mp.end()); mpi++) 
                {
                    alexandria::MolecularCompositionIterator mci = mpi->SearchMolecularComposition(SPOEL);
                    
                    if (mci != mpi->EndMolecularComposition())
                    {
                        int ngt = mci->CountAtoms(gt_type[cur]);
                        if (mp_get_prop(*mpi,emp,iqmExp,lot,NULL,NULL,&exp_val) > 0)
                            nfitexp += ngt;
                        else if (mp_get_prop(*mpi,emp,iqmQM,lot,NULL,NULL,&qm_val) > 0)
                            nfitqm += ngt;
                    }
                }
                if ((estats = gmx_stats_get_npoints(smlsq[j].nexp,&N)) != estatsOK)
                    gmx_fatal(FARGS,"Statistics problems: %s",gmx_stats_message(estats));
                if ((estats = gmx_stats_add_point(smlsq[j].nexp,N,nfitexp,0,0)) != estatsOK)
                    gmx_fatal(FARGS,"Statistics problems: %s",gmx_stats_message(estats));
                if ((estats = gmx_stats_get_npoints(smlsq[j].nqm,&N)) != estatsOK)
                    gmx_fatal(FARGS,"Statistics problems: %s",gmx_stats_message(estats));
                if ((estats = gmx_stats_add_point(smlsq[j].nqm,N,nfitqm,0,0)) != estatsOK)
                    gmx_fatal(FARGS,"Statistics problems: %s",gmx_stats_message(estats));
                
                cur = prev;
            }
        }
    }
    if (NULL != histo)
    {
        nset = 0;
        xvg = xvgropen(histo,"Distribution of polarizabilities","Polarizability","a.u.",oenv);
        for(j = 0; (j<nsmlsq); j++) {
            if ((gmx_stats_get_npoints(smlsq[j].lsq,&N) == estatsOK) && (N > 0))
                fprintf(xvg,"@ s%d legend \"%s\"\n",nset++,smlsq[j].gt_type);
        }
    }
    else
        xvg = NULL;
        
    /* Now print it! */
    atomtype_tab_header(fp,1);
    nset = 0;
    ntab = 0;
    for(j = 0; (j<nsmlsq); j++) {
        /* Determine Miller and Bosque polarizabilities for this Spoel element */
        ahc = ahp = 0;
        if (NULL != gmx_poldata_get_miller(pd[0],smlsq[j].miller_equiv,
                                           &atomnumber,&ahc,&ahp,NULL)) 
        {
            ahc = (4.0/atomnumber)*sqr(ahc);
        }
        if (NULL == gmx_poldata_get_bosque(pd[0],smlsq[j].gt_type,smlsq[j].elem,&bos_pol))
            bos_pol = 0;
                
        /* Compute how many molecules contributed to the optimization 
         * of this polarizability.
         */
        /* Construct group name from element composition */
        strncpy(group,smlsq[j].elem,sizeof(group));
        if ((estats = gmx_stats_get_npoints(smlsq[j].lsq,&N)) != estatsOK)
            gmx_fatal(FARGS,"Statistics problems: %s",gmx_stats_message(estats));
        if (estatsOK != (estats = gmx_stats_get_average(smlsq[j].lsq,&alexandria_aver)))
            gmx_fatal(FARGS,"Statistics problem: %s. gt_type = %s. N = %d.",
                      gmx_stats_message(estats),smlsq[j].gt_type,N);
        if (estatsOK != (estats = gmx_stats_get_sigma(smlsq[j].lsq,&alexandria_sigma)))
            gmx_fatal(FARGS,"Statistics problem: %s. gt_type = %s. N = %d.",
                      gmx_stats_message(estats),smlsq[j].gt_type,N);
        /* Store average values */
        if (NULL != pd_aver) 
        {
            gmx_poldata_set_atype_polarizability(pd_aver,smlsq[j].gt_type,
                                                 alexandria_aver,alexandria_sigma);
        }
        if (estatsOK != (estats = gmx_stats_get_average(smlsq[j].nexp,&nnn)))
            gmx_fatal(FARGS,"Statistics problem: %s. gt_type = %s. N = %d.",
                      gmx_stats_message(estats),smlsq[j].gt_type,N);
        nfitexp = gmx_nint(nnn);
        if (estatsOK != (estats = gmx_stats_get_average(smlsq[j].nqm,&nnn)))
            gmx_fatal(FARGS,"Statistics problem: %s. gt_type = %s. N = %d.",
                      gmx_stats_message(estats),smlsq[j].gt_type,N);
        nfitqm = gmx_nint(nnn);

        fprintf(fp,"%s & %s & %s & %s & %s & %s & %s (%s)",
                smlsq[j].gt_type,
                (nfitexp > 0)     ? gmx_itoa(nfitexp)     : "",
                (nfitqm > 0)      ? gmx_itoa(nfitqm)      : "",
                (ahc > 0)         ? gmx_ftoa(ahc)         : "",
                (ahp > 0)         ? gmx_ftoa(ahp)         : "",
                (bos_pol > 0)     ? gmx_ftoa(bos_pol)     : "",
                (alexandria_aver > 0)  ? gmx_ftoa(alexandria_aver)  : "",
                (alexandria_sigma > 0) ? gmx_ftoa(alexandria_sigma) : "-");
        if (strcasecmp(smlsq[j].spref,"Maaren2009a") == 0)
            fprintf(fp,"& (*)\\\\\n");
        else
            fprintf(fp,"& \\cite{%s}\\\\\n",smlsq[j].spref);
        ntab++;
    
        if (((ntab == nfirst) || (((ntab - nfirst) % nsecond) == 0)) &&
            (j < nsmlsq-1))
        {
            atomtype_tab_end(fp);
            atomtype_tab_header(fp,0);
        }
        if ((NULL != xvg) && (gmx_stats_get_npoints(smlsq[j].lsq,&N) == estatsOK) && (N > 0))
        {
            real *x=NULL,*y=NULL;
            int nbins = 20;
            
            if (gmx_stats_make_histogram(smlsq[j].lsq,0,&nbins,ehistoY,1,&x,&y) == estatsOK)
            {
                fprintf(xvg,"@type xy\n");
                for(i=0; (i<nbins); i++)
                    fprintf(xvg,"%g  %g\n",x[i],y[i]);
                fprintf(xvg,"&\n");
            }
            if (NULL != x)
                sfree(x);
            if (NULL != y)
                sfree(y);
        }
    }
    if (NULL != xvg)
        fclose(xvg);
        
    atomtype_tab_end(fp);
    fflush(fp);
}

static void gmx_molprop_atomtype_dip_table(FILE *fp,gmx_poldata_t pd)
{
    int    i,k,m,cur=0;
    double alexandria_pol,sig_pol;
    char   *elem,*gt_type[2] = { NULL, NULL };
    char   *charge,*miller_equiv,*spref,*desc;
#define prev (1-cur)
#define NEQG 5
    int    eqgcol[NEQG] = { eqgAXp, eqgAXs, eqgAXg };
    int    npcol[NEQG]  = { 2, 3, 3 };
    const char   *clab[3] = { "$J_0$", "$\\chi_0$", "$\\zeta$" };
    int    ncol;
    
    ncol = 1;
    for(i=0; (i<NEQG); i++)
    {
        ncol += npcol[i];
    }
    fprintf(fp,"\\begin{sidewaystable}[H]\n\\centering\n");
    fprintf(fp,"\\caption{Electronegativity equalization parameters for Alexandria models. $J_0$ and $\\chi_0$ in eV, $\\zeta$ in 1/nm.}\n");
    fprintf(fp,"\\label{eemparams}\n");
    fprintf(fp,"\\begin{tabular}{l");
    for(i=1; (i<ncol); i++)
        fprintf(fp,"c");
    fprintf(fp,"}\n\\hline\n");
    for(i=0; (i<NEQG); i++) 
        fprintf(fp," & \\multicolumn{%d}{c}{%s}",npcol[i],get_eemtype_name(eqgcol[i]));
    fprintf(fp,"\\\\\n");
    for(i=0; (i<NEQG); i++) 
        for(m=0; (m<npcol[i]); m++)
            fprintf(fp," & %s",clab[m]);
    fprintf(fp,"\\\\\n\\hline\n");
    while(1 == gmx_poldata_get_atype(pd,&elem,&desc,&(gt_type[cur]),
                                     &miller_equiv,&charge,
                                     &alexandria_pol,&sig_pol,&spref))
    {
        if (((NULL == gt_type[prev]) || (strcmp(gt_type[cur],gt_type[prev]) != 0)))
        {
            fprintf(fp,"%s\n",gt_type[cur]);
            for(k=0; (k<NEQG); k++)
            {
                if (gmx_poldata_have_eem_support(pd,eqgcol[k],gt_type[cur],FALSE))
                {
                    fprintf(fp," & %.3f",gmx_poldata_get_j00(pd,eqgcol[k],gt_type[cur]));
                    fprintf(fp," & %.3f",gmx_poldata_get_chi0(pd,eqgcol[k],gt_type[cur]));
                    if (npcol[k] == 3)
                        fprintf(fp," & %.3f",gmx_poldata_get_zeta(pd,eqgcol[k],gt_type[cur],1));
                    if (npcol[k] == 4)
                        fprintf(fp," & %.3f",gmx_poldata_get_zeta(pd,eqgcol[k],gt_type[cur],2));
                }
                else 
                {
                    for(m=0; (m<npcol[k]); m++)
                        fprintf(fp," & ");
                }
            }
            fprintf(fp,"\\\\\n");
        }
        cur = prev;
    }
    fprintf(fp,"\\hline\n\\end{tabular}\\end{sidewaystable}\n\n");
}

void gmx_molprop_atomtype_table(FILE *fp,gmx_bool bPolar,
                                int npd,gmx_poldata_t pd[],
                                gmx_poldata_t pd_aver,
                                std::vector<alexandria::MolProp> mp,
                                int iQM,char *lot,output_env_t oenv,const char *histo)
{
    if (bPolar)
        gmx_molprop_atomtype_polar_table(fp,npd,pd,pd_aver,mp,iQM,lot,oenv,histo);
    else
        gmx_molprop_atomtype_dip_table(fp,pd[0]);
}
                               

static void prop_header(FILE *fp,int caption,const char *property,real rtoler,real atoler,
                        t_qmcount *qmc,gmx_bool bSideways,int ims,gmx_bool bPrintConf,
                        gmx_bool bPrintBasis,gmx_bool bPrintMultQ)
{
    int  i;
    
    if (!caption)
        decrease_table_number(fp);
    fprintf(fp,"\\begin{%stable}[H]\n",bSideways ? "sideways" : "");
    if (caption) 
    {
        fprintf(fp,"\\caption{Comparison of experimental %s to calculated values. {\\bf Data set: %s}. Calculated numbers that are more than %.0f%s off the experimental values are printed in bold, more than %.0f%s off in bold red.}\n\\label{%s}\n",
                property,ims_names[ims],
                (atoler > 0) ? atoler   : 100*rtoler,(atoler > 0) ? "" : "\\%",
                (atoler > 0) ? 2*atoler : 200*rtoler,(atoler > 0) ? "" : "\\%",
                ims_names[ims]);
    }
    else 
    {
        fprintf(fp,"\\caption{%s, continued}\n",property);
    }
    
    fprintf(fp,"\\begin{tabularx}{23cm}{Xcc");
    if (bPrintMultQ)
        fprintf(fp,"cc");
    if (bPrintConf)
        fprintf(fp,"c");
    for(i=0; (i<qmc->n); i++)
        fprintf(fp,"c");
    fprintf(fp,"}\n\\hline\n");
    fprintf(fp,"Molecule & Form. %s %s & Exper. ",
            bPrintMultQ ? "& q & mult" : "",
            bPrintConf  ? "& Conf." : "");
    for(i=0; (i<qmc->n); i++)
        fprintf(fp,"& %s",qmc->method[i]);
    fprintf(fp,"\\\\\n");
    if (bPrintBasis) 
    {
        fprintf(fp," & & %s %s",
                bPrintMultQ ? "& &" : "",
                bPrintConf ? "&" : "");
        for(i=0; (i<qmc->n); i++)
            fprintf(fp,"& %s",qmc->basis[i]);
        fprintf(fp,"\\\\\n");
    }
    fprintf(fp,"Type & &%s %s",
            bPrintMultQ ? "& &" : "",
            bPrintConf ? "&" : "");
    for(i=0; (i<qmc->n); i++)
        fprintf(fp,"& %s",qmc->type[i]);
    fprintf(fp,"\\\\\n\\hline\n");
}

static void prop_end(FILE *fp,gmx_bool bSideways)
{
    fprintf(fp,"\\hline\n\\end{tabularx}\n\\end{%stable}\n",
            bSideways ? "sideways" : "");
}

static int outside(real vexp,real vcalc,real rtoler,real atoler)
{
    real rdv,adv = fabs(vexp-vcalc);
    if (atoler > 0) 
    {
        if (adv > 2*atoler)
            return 2;
        else if (adv > atoler)
            return 1;
        return 0;
    }
    else 
    {
        if (vexp == 0)
            return 0;
        rdv = adv/vexp;
        
        if (rdv > 2*rtoler)
            return 2;
        else if (rdv > rtoler)
            return 1;
        return 0;
    }
}

void gmx_molprop_prop_table(FILE *fp,int emp,real rtoler,real atoler,
                            std::vector<alexandria::MolProp> mp,
                            int bDS,
                            t_qmcount *qmc,gmx_bool bPrintAll,
                            gmx_bool bPrintBasis,gmx_bool bPrintMultQ,
                            gmx_molselect_t gms,int ims)
{
    std::vector<alexandria::MolProp>::iterator mpi;
    alexandria::ExperimentIterator ei;
    alexandria::CalculationIterator ci;
    alexandria::MolecularDipPolarIterator mdi;
    alexandria::MolecularQuadrupoleIterator qi;
    alexandria::MolecularEnergyIterator mei;
    std::vector<double> val_exp,val_calc,err_exp,err_calc;
    std::vector<std::string> ref_exp,conf_exp,type_exp;
    std::vector<int> found_calc;
    std::string exp_conf;
    
    int    i,j,iprint=0,ne,iline,caption=1,maxline;
#define BLEN 1024
    char   lbuf[BLEN],myline[BLEN],mylbuf[BLEN],vbuf[BLEN];
    double calc_val,calc_err,vc;
    int    nprint,header;
    double dvec[DIM];
    tensor quadrupole;
    real   ds_fac;
    gmx_bool bSideways,bPrintConf,bOutlier;
  
    bSideways = (qmc->n > 1);
    
    if (bSideways)
        maxline = 20;
    else
        maxline = 40;
    /* Double spaced text ? */
    ds_fac = 1-0.4*bDS;

    nprint = 0;
    for(mpi=mp.begin(); (mpi<mp.end()); mpi++)
    {
        if ((ims == gmx_molselect_status(gms,mpi->GetIupac().c_str())) &&
            (mpi->HasComposition(SPOEL)))
            nprint++;
    }
    if (nprint <= 0)
        return;
    bPrintConf = (emp == empDIPOLE);
    prop_header(fp,caption,emp_name[emp],rtoler,atoler,qmc,bSideways,ims,bPrintConf,
                bPrintBasis,bPrintMultQ);  
    header  = 1;
    iline   = 0;
    for(mpi=mp.begin(); (mpi<mp.end()); mpi++)
    {
        if ((ims == gmx_molselect_status(gms,mpi->GetIupac().c_str())) && 
            (mpi->HasComposition(SPOEL)))
        {
            for(ei=mpi->BeginExperiment(); (ei<mpi->EndExperiment()); ei++)
            {
                switch (emp) 
                {
                case empDIPOLE:
                    for(mdi=ei->BeginDipole(); (mdi<ei->EndDipole()); mdi++)
                    {
                        val_exp.push_back(mdi->GetAver());
                        err_exp.push_back(mdi->GetError());
                        ref_exp.push_back(ei->GetReference());
                        conf_exp.push_back(ei->GetConformation());
                        type_exp.push_back(mdi->GetType());
                    }
                    break;
                case empPOLARIZABILITY:
                    for(mdi=ei->BeginPolar(); (mdi<ei->EndPolar()); mdi++)
                    {
                        val_exp.push_back(mdi->GetAver());
                        err_exp.push_back(mdi->GetError());
                        ref_exp.push_back(ei->GetReference());
                        conf_exp.push_back(ei->GetConformation());
                        type_exp.push_back(mdi->GetType());
                    }
                    break;
                case empENERGY:
                    for(mei=ei->BeginEnergy(); (mei<ei->EndEnergy()); mei++)
                    {
                        val_exp.push_back(mei->GetValue());
                        err_exp.push_back(mei->GetError());
                        ref_exp.push_back(ei->GetReference());
                        conf_exp.push_back(ei->GetConformation());
                        type_exp.push_back(mei->GetType());
                    }
                    break;
                default:
                    gmx_fatal(FARGS,"No support for for emp %d",emp);
                    break;
                }
            }
            for(ne=0; (ne < (val_exp.size() > 1 ? val_exp.size() : 1)); ne++)
            {
                if (val_exp.size() <= 1) 
                {
                    exp_conf.erase();
                }
                else
                {
                    exp_conf = conf_exp[ne];
                }
                for(j=0; (j<qmc->n); j++) 
                {
                    sprintf(lbuf,"%s/%s",qmc->method[j],qmc->basis[j]);
                    if (mp_get_prop_ref(*mpi,emp,iqmQM,lbuf,exp_conf.c_str(),
                                        qmc->type[j],&calc_val,&calc_err,
                                        NULL,NULL,dvec,quadrupole) == 1) 
                    {
                        found_calc.push_back(1);
                        val_calc.push_back(calc_val);
                        err_calc.push_back(calc_err);
                    }
                    else 
                    {
                        found_calc.push_back(0);
                        val_calc.push_back(0);
                        err_calc.push_back(0);
                    }
                }
                if ((bPrintAll || (val_exp.size() > 0))  && (found_calc.size() > 0))
                {
                    if (0 == ne)
                        iprint++;
                    myline[0] = '\0';
                    if ((ne == 0))
                    {
                        if (bPrintMultQ)
                            sprintf(myline,"%d. %-15s & %s & %d & %d",
                                    iprint,mpi->GetIupac().c_str(),
                                    mpi->GetFormula().c_str(),
                                    mpi->GetCharge(),
                                    mpi->GetMultiplicity());
                        else
                            sprintf(myline,"%d. %-15s & %s",
                                    iprint,mpi->GetIupac().c_str(),
                                    mpi->GetFormula().c_str());
                        header = 0;
                    }
                    else 
                    {
                        sprintf(myline," & ");
                    }
                    if (bPrintConf)
                    {
                        sprintf(mylbuf,"      & %s ",(exp_conf.size() > 0) ? exp_conf.c_str() : "-");
                        strncat(myline,mylbuf,BLEN-strlen(myline));
                    }
                    if (val_exp.size() > 0) 
                    {
                        sprintf(mylbuf,"& %8.3f",val_exp[ne]);
                        strncat(myline,mylbuf,BLEN-strlen(myline));
                        if (strcmp(ref_exp[ne].c_str(),"Maaren2013a") == 0)
                            sprintf(mylbuf," (*)");
                        else
                            sprintf(mylbuf,"~\\cite{%s} ",ref_exp[ne].c_str());
                        strncat(myline,mylbuf,BLEN-strlen(myline));
                    }
                    else 
                    {
                        sprintf(mylbuf,"& - ");
                        strncat(myline,mylbuf,BLEN-strlen(myline));
                    }
                    bOutlier = FALSE;
                    for(j=0; (j<qmc->n); j++) 
                    { 
                        if (found_calc[j] > 0) 
                        {
                            vc = val_calc[j];
                            if (err_calc[j] > 0) 
                            {
                                sprintf(vbuf,"%8.2f(%.2f)",vc,err_calc[j]);
                            }
                            else
                            {
                                sprintf(vbuf,"%8.2f",vc);
                            }
                            if (val_exp.size() > 0) 
                            {
                                int oo = outside(val_exp[ne],vc,rtoler,atoler);
                                switch(oo) {
                                case 2:
                                    sprintf(mylbuf,"& \\textcolor{Red}{\\bf %s} ",vbuf);
                                    bOutlier = TRUE;
                                    break;
                                case 1:
                                    sprintf(mylbuf,"& {\\bf %s} ",vbuf);
                                    bOutlier = TRUE;
                                    break;
                                default:
                                    sprintf(mylbuf,"& %s ",vbuf);
                                }
                            }
                            else
                                sprintf(mylbuf,"& %s ",vbuf);
                            strncat(myline,mylbuf,BLEN-strlen(myline));
                        }
                        else
                        {
                            sprintf(mylbuf,"& ");
                            strncat(myline,mylbuf,BLEN-strlen(myline));
                        }
                    }
                    sprintf(mylbuf,"\\\\\n");
                    strncat(myline,mylbuf,BLEN-strlen(myline));
                    /*if ((toler > 0) && bOutlier)
                      {*/
                        fprintf(fp,"%s",myline);
                        iline++;/*=(1+strlen(molname)/25);*/
                        /*}*/
                }
                if ((iline >= (maxline-7*caption)*ds_fac)) 
                {
                    caption = 0;
                    prop_end(fp,bSideways);
                    prop_header(fp,caption,emp_name[emp],rtoler,atoler,qmc,
                                bSideways,ims,bPrintConf,bPrintBasis,bPrintMultQ);
                    header  = 1;
                    iline   = 0;
                }
            }
        }
    }
    prop_end(fp,bSideways);
}

