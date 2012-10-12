/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: analyze_mp.c,v 1.9 2009/05/26 14:40:49 spoel Exp $
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
#include "maths.h"
#include "futil.h"
#include "smalloc.h"
#include "string2.h"
#include "vec.h"
#include "xvgr.h"
#include "statutil.h"
#include "copyrite.h"
#include "physics.h"
#include "gstat.h"
#include "gmx_fatal.h"
#include "molprop.h"
#include "molprop_util.h"
#include "gromacs/linearalgebra/matrix.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "molprop_xml.h"
#include "molprop_tables.h"
#include "molselect.h"

static void calc_frag_miller(int bTrain,gmx_poldata_t pd,
                             int np,gmx_molprop_t mp[],gmx_molselect_t gms)
{
    int    j,k,Nelec,natom,natom_tot,charge,hybridization,nhydrogen,atomnumber,calcref;
    double ahc,ahp,bos,spoel,sig_alexandria,bos0,polar,sig_pol,blength;
    double tau_ahc,alpha_ahp;
    char   *atomname,*gt_type,*elem,*miller_equiv,*ref,*neighbor1,*neighbor2;
    const char *comp,*iupac;
    int    alexandria_support,miller_support,bosque_support,ims;
    char *null = (char *)"0",*empirical = (char *)"empirical",*minimum = (char *)"minimum",*minus = (char *)"-";
    const char *program;
    const char *ang3;
    
    ang3    = unit2string(eg2cAngstrom3);
    program = ShortProgram();
    if (gmx_poldata_get_bosque(pd,null,null,&bos0) == NULL)
        gmx_fatal(FARGS,"Can not find Bosque polarizability for 0");
  
    for(j=0; (j<np); j++) 
    {
        iupac = gmx_molprop_get_iupac(mp[j]);
        if (strcasecmp(iupac,"bromomethane") == 0)
            fprintf(stderr,"KOKO: %s\n",iupac);
        ims = gmx_molselect_status(gms,iupac);
        if ((ims == imsTrain) || (ims == imsTest)) 
        {
            ahc = ahp = spoel = sig_alexandria = 0;
            bos = bos0;
            Nelec = 0;
            natom_tot = 0;
            alexandria_support  = -1;
            miller_support = -1;
            bosque_support = -1;
      
            while ((comp = gmx_molprop_get_composition(mp[j])) != NULL) {
                while (gmx_molprop_get_composition_atom(mp[j],comp,&atomname,&natom) == 1) 
                {
                    if (strcasecmp(comp,"spoel") == 0) 
                    {
                        if (gmx_poldata_type_polarizability(pd,atomname,&polar,&sig_pol) == 1) 
                        {
                            spoel     += polar*natom;
                            sig_alexandria += sqr(sig_pol)*natom;
                            natom_tot += natom;
                            if (alexandria_support == -1)
                                alexandria_support = 1;
                        }
                        else 
                            alexandria_support = 0;
                    }
                    else if (strcasecmp(comp,"miller") == 0) 
                    {
                        if (gmx_poldata_get_miller(pd,atomname,&atomnumber,
                                                   &tau_ahc,&alpha_ahp,NULL) != NULL) 
                        {
                            ahc   += tau_ahc*natom; 
                            ahp   += alpha_ahp*natom;
                            Nelec += atomnumber*natom;
                            if (miller_support == -1)
                                miller_support = 1;
                        }
                        else
                            miller_support = 0;
                    }
                    else if (strcasecmp(comp,"bosque") == 0) 
                    {
                        if (gmx_poldata_get_bosque(pd,NULL,atomname,&polar) != NULL) 
                        {
                            bos += polar*natom;
                            if (bosque_support == -1)
                                bosque_support = 1;
                        }
                        else
                            bosque_support = 0;
                    }
                    else
                        gmx_fatal(FARGS,"Unknown composition %s\n",comp);
                }
            }
            
            if (alexandria_support == 1) 
            {
                gmx_molprop_add_calculation(mp[j],program,(char *)"AX",minus,
                                            (char *)"Maaren2010a",minimum,&calcref);
                sig_alexandria = sqrt(sig_alexandria/natom_tot);
                gmx_molprop_add_polar(mp[j],calcref,empirical,ang3,
                                      0,0,0,spoel,sig_alexandria);
            }
            if (bosque_support == 1) 
            {
                gmx_molprop_add_calculation(mp[j],program,(char *)"BS",minus,
                                            (char *)"Bosque2002a",minimum,&calcref);
                gmx_molprop_add_polar(mp[j],calcref,empirical,ang3,
                                      0,0,0,bos,0);
            }
            if (miller_support == 1) 
            {
                ahc = 4*sqr(ahc)/Nelec;
                gmx_molprop_add_calculation(mp[j],program,(char *)"MK",(char *)"ahc",
                                            (char *)"Miller1990a,Khan",minimum,&calcref);
                gmx_molprop_add_polar(mp[j],calcref,empirical,ang3,
                                      0,0,0,ahc,0);
                gmx_molprop_add_calculation(mp[j],program,(char *)"MK",(char *)"ahp",
                                            (char *)"Miller1990a,Khan",minimum,&calcref);
                gmx_molprop_add_polar(mp[j],calcref,empirical,ang3,
                                      0,0,0,ahp,0);
            }
        }
    }
}

static void write_corr_xvg(const char *fn,int np,gmx_molprop_t pd[],
                           int emp,t_qmcount *qmc,
                           int iQM,char *lot,double outlier_frac,
                           output_env_t oenv,gmx_molselect_t gms)
{
    FILE   *fp;
    int    i,j,k,n,nout,ims;
    char   lbuf[256];
    double exp_val,qm_val; 
    
    fp  = xvgropen(fn,"","Exper.","Calc. - Exper.",oenv);
    for(i=0; (i<qmc->n); i++) 
    {
        fprintf(fp,"@s%d legend \"%s/%s-%s\"\n",i,qmc->method[i],qmc->basis[i],
                qmc->type[i]);
        fprintf(fp,"@s%d line linestyle 0\n",i);
        fprintf(fp,"@s%d symbol %d\n",i,i+1);
        fprintf(fp,"@s%d symbol size %g\n",i,0.5);
        fprintf(fp,"@s%d symbol fill color %d\n",i,i+1);
        fprintf(fp,"@s%d symbol fill pattern 1\n",i);
    }
    for(i=0; (i<qmc->n); i++)
    {
        fprintf(fp,"@type xy\n");
        sprintf(lbuf,"%s/%s",qmc->method[i],qmc->basis[i]);
        if (debug)
            fprintf(debug,"QM: %s\n",qmc->lot[i]);
        nout = 0;
        for(j=0; (j<np); j++) 
        {
            ims =gmx_molselect_status(gms,gmx_molprop_get_iupac(pd[j]));
            if ((ims == imsTrain) || (ims == imsTest))
            {
                for(k=0; (k<qmc->nconf); k++)
                {
                    if ((mp_get_prop(pd[j],emp,iqmExp,NULL,NULL,NULL,&exp_val) > 0)  &&
                        (mp_get_prop(pd[j],emp,iqmQM,lbuf,qmc->conf[k],qmc->type[i],&qm_val) > 0)) 
                    {
                        fprintf(fp,"%8.3f  %8.3f\n",exp_val,qm_val-exp_val);
                        if (debug && (fabs(qm_val-exp_val)/exp_val > outlier_frac))
                        {
                            fprintf(debug,"OUTLIER: %s Exp: %g, Calc: %g\n",
                                    gmx_molprop_get_iupac(pd[j]),exp_val,qm_val);
                            nout++;
                        }
                    }
                }
            }
        }
        fprintf(fp,"&\n");
        if (debug)
            fprintf(debug,"There were %3d outliers (more than %g%% off) for %s.\n",
                    nout,100*outlier_frac,qmc->lot[i]);
    }
    fclose(fp);
    do_view(oenv,fn,NULL);
}

typedef struct {
    int nref;
    int *rcount;
    char **refs;
} t_refcount;

static void add_refc(t_refcount *rc,char *ref)
{
    int i;
    
    for(i=0; (i<rc->nref); i++) 
    {
        if (strcmp(ref,rc->refs[i]) == 0) 
        {
            rc->rcount[i]++;
            break;
        }
    }
    if (i == rc->nref) 
    {
        rc->nref++;
        srenew(rc->rcount,rc->nref);
        srenew(rc->refs,rc->nref);
        rc->rcount[i] = 1;
        rc->refs[i] = strdup(ref);
    }
}

static void gmx_molprop_analyze(int np,gmx_molprop_t mp[],int npd,
                                gmx_poldata_t *pd,gmx_poldata_t pdref,
                                gmx_bool bCalcPol,int prop,
                                gmx_atomprop_t ap,int iQM,char *lot,real toler,
                                char *fc_str,gmx_bool bPrintAll,double outlier,
                                gmx_bool bStatsTable,
                                gmx_bool bPropTable,gmx_bool bCompositionTable,
                                const char *texfn,
                                const char *xvgfn,
                                const char *histo,
                                const char *atype,
                                output_env_t oenv,
                                gmx_molselect_t gms,
                                const char *selout)
{
    FILE *fp,*gp;
    int  i,ntot,cur=0;
    t_qmcount *qmc;
    t_refcount *rc;
    char *reference,*conformation,*mylot; 
    const char *molname[2];
    double value,error,vec[3];
    tensor quadrupole;
#define prev (1-cur)
    int  expref;
    const char   *iupac;
    char *ref;
  
    if (NULL != atype)
    {
        fp = fopen(atype,"w");
        gmx_molprop_atomtype_table(fp,(prop == empPOLARIZABILITY),
                                   npd,pd,pdref,np,mp,iQM,lot,oenv,histo);
        fclose(fp);
        do_view(oenv,histo,NULL);
    }
    
    if (bCalcPol)
        calc_frag_miller(FALSE,pdref ? pdref : pd[0],np,mp,gms);
        
    qmc = find_calculations(np,mp,prop,fc_str);
    
    snew(rc,1);
    for(i=0; (i<np); i++) 
    {
        molname[cur] = gmx_molprop_get_molname(mp[i]);
        if ((gmx_molprop_get_experiment(mp[i],&reference,&conformation,&expref) == 1) &&
            (get_val(mp[i],expref,NULL,prop,&value,&error,vec,quadrupole)  == 1)) 
        {
            add_refc(rc,reference);
            sfree(reference);
            sfree(conformation);
        }
        if (debug && ((i > 0) && (strcasecmp(molname[cur],molname[prev]) == 0)))
            fprintf(debug,"Double entry %s\n",molname[cur]);
        cur=prev;
    }
    
    printf("--------------------------------------------------\n");
    printf("      Some statistics for %s\n",emp_name[prop]);
    ntot = 0;
    for(i=0; (i<rc->nref); i++)
    {
        printf("There are %d experiments with %s as reference\n",
               rc->rcount[i],rc->refs[i]);
        ntot += rc->rcount[i];
    }
    printf("There are %d entries with experimental data\n",ntot);
    for(i=0; (i<qmc->n); i++) 
        printf("There are %d calculation results using %s/%s type %s\n",
               qmc->count[i],qmc->method[i],qmc->basis[i],qmc->type[i]);
    printf("--------------------------------------------------\n");
    
    fp = ffopen(texfn,"w");
    if (bStatsTable) 
    {
        gmx_molprop_stats_table(fp,prop,np,mp,ntot,qmc,iQM,lot,outlier,gms,imsTrain);
        gmx_molprop_stats_table(fp,prop,np,mp,ntot,qmc,iQM,lot,outlier,gms,imsTest);
    }
    if (bPropTable)
    {
        gmx_molprop_prop_table(fp,prop,toler,np,mp,0,qmc,bPrintAll,gms,imsTrain);
        gmx_molprop_prop_table(fp,prop,toler,np,mp,0,qmc,bPrintAll,gms,imsTest);
        if (NULL != selout) 
        {
            gp = fopen(selout,"w");
            for(i=0; (i<np); i++) 
            {
                iupac = gmx_molprop_get_iupac(mp[i]);
                if ((NULL != iupac) && (strlen(iupac) > 0))
                {
                    if (mp_get_prop_ref(mp[i],prop,iqmBoth,lot,NULL,NULL,
                                        &value,&error,&ref,&mylot,
                                        vec,quadrupole) == 1)
                    {
                        fprintf(gp,"%s|Train\n",iupac);
                    }
                }
            }
            fclose(gp);
        }
    }
    if (bCompositionTable)
    {
        gmx_molprop_composition_table(fp,np,mp,gms,imsTrain);
        gmx_molprop_composition_table(fp,np,mp,gms,imsTest);
    }
    fclose(fp);
    if (NULL != xvgfn)
        write_corr_xvg(xvgfn,np,mp,prop,qmc,FALSE,lot,outlier,oenv,gms);
}


int main(int argc,char *argv[])
{
    static const char *desc[] = {
        "analyze_mp reads force field information and a molecule database",
        "and produces tables and figures to describe the data.[PAR]",
        "A series of force field files can be read for generating statistics",
        "on the parameters, like average and standard deviation.",
        "The average values and errors can be written to a new force field file",
        "([TT]-pout[tt] flag) based on a reference file given by the",
        "[TT]-pref[tt] flag.[PAR]",
        "A selection of molecules into a training set and a test set (or ignore set)",
        "can be made using option [TT]-sel[tt]. The format of this file is:[BR]",
        "iupac|Train[BR]",
        "iupac|Test[BR]",
        "iupac|Ignore[BR]",
        "and you should ideally have a line for each molecule in the molecule database",
        "([TT]-m[tt] option). Missing molecules will be ignored. You can also write a",
        "selection file ([TT]-selout[tt]) that contains all the molecules in the",
        "output corresponding to the [TT]-prop[tt] flag. You can steer the content of",
        "this file e.g. by running the program with an empty [TT]-lot[tt] flag,",
        "yielding only those molecules for which experimental data is available."
    };
    t_filenm fnm[] = {
        { efDAT, "-p",    "poldata",   ffRDMULT },
        { efDAT, "-pref", "gentop",    ffOPTRD  },
        { efDAT, "-pout", "gentop",    ffOPTWR  },
        { efDAT, "-m",    "allmols",   ffRDMULT },
        { efTEX, "-t",    "table",     ffWRITE  },
        { efTEX, "-atype","atomtypes", ffOPTWR  },
        { efDAT, "-sel",  "molselect", ffREAD   },
        { efDAT, "-selout","selout",   ffOPTWR  },
        { efXVG, "-c",    "correl",    ffWRITE  },
        { efXVG, "-his",  "polhisto",  ffWRITE  }
    };
#define NFILE (sizeof(fnm)/sizeof(fnm[0]))
    static char *sort[] = { NULL, (char *)"molname", (char *)"formula", (char *)"composition", (char *)"selection", NULL };
    static char *prop[] = { NULL, (char *)"potential", (char *)"dipole", (char *)"quadrupole", (char *)"polarizability", (char *)"energy", (char *)"DHf(298.15K)", NULL };
    static char *fc_str = (char *)"";
    static char *lot = (char *)"B3LYP/aug-cc-pVTZ";
    static real toler = 0.15,outlier=0;
    static real th_toler=170,ph_toler=5;
    static gmx_bool bMerge = TRUE,bAll = FALSE,bCalcPol=TRUE;
    static gmx_bool bStatsTable = TRUE,bCompositionTable=FALSE,bPropTable=TRUE;
    t_pargs pa[] = {
        { "-sort",   FALSE, etENUM, {sort},
          "Key to sort the final data file on." },
        { "-tol",    FALSE, etREAL, {&toler},
          "Relative tolerance for printing in bold in tables" },
        { "-merge",  FALSE, etBOOL, {&bMerge},
          "Merge molecule records in the input file and generate atomtype compositions based on the calculated geometry." },
        { "-lot",    FALSE, etSTR,  {&lot},
          "Indicate the method and level of theory that were used together with experimental data in refining polarizabilities. If empty, is is assumed that only experimental data were used." },
        { "-prop",   FALSE, etENUM, {prop},
          "Property to print" },
        { "-calcpol", FALSE, etBOOL, {&bCalcPol},
          "Calculate polarizabilities based on empirical methods" },
        { "-composition", FALSE, etBOOL, {&bCompositionTable},
          "Print a table of composition of the molecules" },
        { "-proptable", FALSE, etBOOL, {&bPropTable},
          "Print a table of properties (slect which with the [TT]-prop[tt] flag)." },
        { "-statstable", FALSE, etBOOL, {&bStatsTable},
          "Print a table of statistics per category" },
        { "-th_toler",FALSE, etREAL, {&th_toler},
          "If bond angles are larger than this value the group will be treated as a linear one and a virtual site will be created to keep the group linear" },
        { "-ph_toler", FALSE, etREAL, {&ph_toler},
          "If dihedral angles are less than this (in absolute value) the atoms will be treated as a planar group with an improper dihedral being added to keep the group planar" },
        { "-outlier", FALSE, etREAL, {&outlier},
          "Outlier indicates a level (in relative units (i.e. (Calc-Exp)/Exp), one standard deviation). Calculations that deviate more than this level from the experiment are not taken into account when computing statistics. Moreover, the outliers are printed to the standard error. If outlier is 0, no action is take." },
        { "-all",    FALSE, etBOOL, {&bAll},
          "Print all values, also those for no data have been computed (e.g. strange metals)" },
        { "-fc_str", FALSE, etSTR, {&fc_str},
          "Selection of the stuff you want in the tables, given as a single string with spaces like: method1/basis1/type1:method2/basis2/type2 (you may have to put quotes around the whole thing in order to prevent the shell from interpreting it)." }
    };
    int        npa;
    FILE       *fp,*gp;
    int        i,alg,np,nspoel,nbosque,nhandbook,ntot,nqm,eMP,eprop;
    double     *fit2,*test2;
    int        cur = 0;
#define prev   (1-cur)

    gmx_molprop_t   *mp=NULL;
    gmx_atomprop_t  ap;
    gmx_poldata_t   *pd,pdref=NULL;
    output_env_t    oenv;
    gmx_molselect_t gms;
    char            **mpname=NULL,**fns=NULL;
    const char      *fn;
    int             npdfile,nmpfile;
    
    CopyRight(stdout,argv[0]);
    
    npa = sizeof(pa)/sizeof(pa[0]);
    parse_common_args(&argc,argv,PCA_CAN_VIEW,NFILE,fnm,
                      npa,pa,
                      sizeof(desc)/sizeof(desc[0]),desc,
                      0,NULL,&oenv);
		    
    ap = gmx_atomprop_init();
    
    nmpfile = opt2fns(&mpname,"-m",NFILE,fnm);
    npdfile = opt2fns(&fns,"-p",NFILE,fnm);
    gms = gmx_molselect_init(opt2fn("-sel",NFILE,fnm));
    
    alg = -1;
    if (opt2parg_bSet("-sort",npa,pa)) 
    {
        for(i=0; (i<empSORT_NR); i++)
            if (strcasecmp(sort[0],sort[i+1]) == 0) 
            {
                alg = i;
                break;
            }
    }
    eprop = -1;
    if (opt2parg_bSet("-prop",npa,pa)) 
    {
        for(i=0; (i<empNR); i++)
            if (strcasecmp(prop[0],prop[i+1]) == 0) 
            {
                eprop = i;
                break;
            }
    }
    if (eprop == -1)
    {
        eprop = empDIPOLE;
    }
    snew(pd,npdfile);
    for(i=0; (i<npdfile); i++) {
        if ((pd[i] = gmx_poldata_read(fns[i],ap)) == NULL)
            gmx_fatal(FARGS,"Can not read the force field information. File missing or incorrect.");
    }
    
    fn = opt2fn_null("-pref",NFILE,fnm);
    if (NULL != fn)
        pdref = gmx_poldata_read(fn,ap);
    if (bMerge)
        mp = merge_xml(nmpfile,mpname,NULL,NULL,NULL,&np,ap,pd[0],TRUE,TRUE,th_toler,ph_toler);
    else 
        mp = gmx_molprops_read(mpname[0],&np);
    if (alg != -1)
    {
        gmx_molprop_sort(np,mp,alg,ap,gms);
    }
    gmx_molprop_analyze(np,mp,npdfile,pd,pdref,bCalcPol,
                        eprop,ap,0,lot,toler,fc_str,bAll,outlier,
                        bStatsTable,bPropTable,bCompositionTable,
                        opt2fn("-t",NFILE,fnm),opt2fn("-c",NFILE,fnm),
                        opt2fn("-his",NFILE,fnm),opt2fn("-atype",NFILE,fnm),oenv,gms,
                        opt2fn_null("-selout",NFILE,fnm));

    if (NULL != fn)
        gmx_poldata_write(opt2fn("-pout",NFILE,fnm),pdref,ap,FALSE);
    
    thanx(stdout);
  
    return 0;
}
