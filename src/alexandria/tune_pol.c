/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: tune_pol.c,v 1.32 2009/05/20 06:25:20 spoel Exp $
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
#include "statutil.h"
#include "copyrite.h"
#include "gstat.h"
#include "gmx_fatal.h"
#include "names.h"
#include "molprop.h"
#include "molprop_util.h"
#include "gromacs/linearalgebra/matrix.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "molprop_xml.h"
#include "molprop_tables.h"
#include "molselect.h"

static int decompose_frag(FILE *fp,int bTrain,
                          gmx_poldata_t pd,int np,gmx_molprop_t mp[],
                          int iQM,char *lot,int mindata,gmx_molselect_t gms,
                          real sigma,gmx_bool bZero,gmx_bool bForceFit)
{
    FILE   *csv;
    double *x,*atx;
    double **a,**at,**ata,*fpp;
    double pol,poltot,a0,da0,ax,sig_pol,blength,pval,chi2;
    char   *elem,*miller_equiv,**atype=NULL,*spref,*gt_type,*atomname;
    const char *iupac,*molname;
    char   *neighbors,*geometry,*desc,*charge;
    int    i,j,niter=0,nusemol=0,nn,nnn,nexp,numbonds,ntp,natom,natom_tot,ca;
    double valence;
    int    *test=NULL,ntest=0,row,ims;
    gmx_bool   bPol,*bUseMol;

    snew(bUseMol,np+1);
    snew(x,np+1);
    poltot = 0;
    for(j=0; (j<np); j++) 
    {
        iupac = gmx_molprop_get_iupac(mp[j]);
        molname  = gmx_molprop_get_molname(mp[j]);
        pol  = 0;
        ims  = gmx_molselect_status(gms,iupac);
        bPol = (mp_get_prop(mp[j],empPOLARIZABILITY,iQM,lot,NULL,NULL,&pol) > 0);
        bUseMol[j] = ((imsTrain == ims) && bPol && (pol > 0));

        /* Now check for the right composition */        
        if (bUseMol[j])
        {
            natom_tot = 0;
            do {
                ca = gmx_molprop_get_composition_atom(mp[j],(char *)"spoel",&atomname,&natom);
                if (ca != 0) 
                {
                    bUseMol[j] = gmx_poldata_have_pol_support(pd,atomname);
                    sfree(atomname);
                    if (bUseMol[j])
                        natom_tot += natom;
                }
                else if (natom_tot == 0)
                    bUseMol[j] = FALSE;
            } while (bUseMol[j] && (ca != 0));
        }
        if (NULL != debug)
        {
            fprintf(debug,"Mol: %s, IUPAC: %s, ims: %s, bPol: %s, pol: %g - %s\n",
                    molname,iupac,ims_names[ims],bool_names[bPol],pol,
                    bUseMol[j] ? "Used" : "Not used");
        }
        
        if (bUseMol[j]) 
        {        
            x[nusemol] = pol;
            poltot += pol;
            nusemol++;
        }
    }        
    while (1 == gmx_poldata_get_atype(pd,&elem,&desc,&gt_type,&miller_equiv,
                                      &charge,&pol,&sig_pol,&spref)) 
    {
        if ((pol == 0) || bForceFit) {
            ntp = 0;
            for(j=0; (j<np); j++) 
            {
                if ((bUseMol[j]) &&
                    ((nnn = gmx_molprop_count_composition_atoms(mp[j],(char *)"spoel",gt_type)) > 0))
                {
                    ntp += nnn;
                }
            }
            if (ntp >= mindata) 
            {
                for(i=0; (i<ntest); i++) 
                {
                    if (strcmp(atype[i],gt_type) == 0) 
                        break;
                }
                if (i == ntest) 
                {
                    srenew(test,ntest+1);
                    srenew(atype,ntest+1);
                    test[ntest] = ntp;
                    atype[ntest] = strdup(gt_type);
                    ntest++;
                }
            }
            else
                fprintf(stderr,"Not enough polarizability data points (%d out of %d required) to optimize %s\n",ntp,mindata,gt_type);
        }
        else 
        {
            fprintf(stderr,"Polarizability for %s is %g. Not optimizing this one.\n",
                    gt_type,pol);
        }
    }
    if (ntest == 0) 
        gmx_fatal(FARGS,"Nothing to optimize. Check your input");

    a      = alloc_matrix(nusemol,ntest);
    at     = alloc_matrix(ntest,nusemol);
    ata    = alloc_matrix(ntest,ntest);
    
    printf("There are %d different atomtypes to optimize the polarizabilities\n",
           ntest);
    printf("There are %d (experimental) reference polarizabilities.\n",nusemol);
    csv = ffopen("out.csv","w");
    for(i=nn=0; (i<np); i++) 
    {
        if (bUseMol[i])
        {
            fprintf(csv,"\"%s\",",gmx_molprop_get_molname(mp[i]));
            for(j=0; (j<ntest); j++) {
                a[nn][j] = at[j][nn] = 
                    gmx_molprop_count_composition_atoms(mp[i],(char *)"spoel",atype[j]);
                fprintf(csv,"\"%g\",",a[nn][j]);
            }
            fprintf(csv,"\"%.3f\"\n",x[nn]);
            nn++;
        }
    }
    fclose(csv);
    if (nusemol != nn) 
        gmx_fatal(FARGS,"Consistency error: nusemol = %d, nn = %d",nusemol,nn);
    if (fp)
        for(i=0; (i<ntest); i++)
            fprintf(fp,"Optimizing polarizability for %s with %d copies\n",
                    atype[i],test[i]);

    matrix_multiply(fp,nusemol,ntest,a,at,ata);
    if ((row = matrix_invert(fp,ntest,ata)) != 0) {
        gmx_fatal(FARGS,"Matrix inversion failed. Incorrect row = %d.\nThis probably indicates that you do not have sufficient data points, or that some parameters are linearly dependent.",
                  row);
    }
    snew(atx,ntest);
    snew(fpp,ntest);
    a0 = 0;
    do {
        for(i=0; (i<ntest); i++)  
        {
            atx[i] = 0;
            for(j=0; (j<nusemol); j++)
                atx[i] += at[i][j]*(x[j]-a0);
        }
        for(i=0; (i<ntest); i++) 
        {
            fpp[i] = 0;
            for(j=0; (j<ntest); j++)
                fpp[i] += ata[i][j]*atx[j];
        }
        da0 = 0;
        chi2 = 0;
        if (bZero)
        {
            for(j=0; (j<nusemol); j++)
            {
                ax = a0;
                for(i=0; (i<ntest); i++)  
                    ax += fpp[i]*a[j][i];
                da0 += (x[j]-ax);
                chi2 += sqr(x[j]-ax);
            }
            da0 = da0 / nusemol;
            a0 += da0;
            niter++;
            printf("iter: %d <pol> = %g, a0 = %g, chi2 = %g\n",
                   niter,poltot/nusemol,a0,chi2/nusemol);
        }
    } while (bZero && (fabs(da0) > 1e-5) && (niter < 1000));
    
    for(i=0; (i<ntest); i++) 
    {
        gmx_poldata_set_atype_polarizability(pd,atype[i],fpp[i],sigma*sqrt(ata[i][i]));
    }
    if (bZero)
        gmx_poldata_add_atype(pd,(char *)"0",(char *)"NUL",(char *)"NUL",(char *)"NUL",
                              (char *)"",a0,0,(char *)"");

    sfree(bUseMol);
    sfree(fpp);
    
    return ntest;
}

int main(int argc,char *argv[])
{
    static const char *desc[] = 
    {
        "tune_pol optimes atomic polarizabilities that together build",
        "an additive model for polarization. The set of atomtypes used",
        "is determined by the input force field file ([TT]-di[tt] option). All",
        "atomtypes for which the polarizability is zero, and for which",
        "there is sufficient data in the input data file ([TT]-f[tt] option)",
        "will be used in the least-squares fit (done by matrix inversion).[PAR]"
        "A selection of molecules into a training set and a test set (or ignore set)",
        "can be made using option [TT]-sel[tt]. The format of this file is:[BR]",
        "iupac|Train[BR]",
        "iupac|Test[BR]",
        "iupac|Ignore[BR]",
        "and you should ideally have a line for each molecule in the molecule database",
        "([TT]-f[tt] option). Missing molecules will be ignored."
    };
    t_filenm fnm[] = 
    {
        { efDAT, "-f",  "data",      ffRDMULT },
        { efDAT, "-o",  "allmols",   ffWRITE },
        { efDAT, "-di", "gentop",    ffOPTRD },
        { efDAT, "-do", "tune_pol",  ffWRITE },
        { efDAT, "-sel", "molselect",ffREAD },
        { efXVG, "-x",  "pol_corr",  ffWRITE }
    };
    int NFILE = (sizeof(fnm)/sizeof(fnm[0]));
    static char *sort[] = { NULL, (char *)"molname", (char *)"formula", (char *)"composition", NULL };
    static int iQM = FALSE,mindata=1;
    static char *lot = "B3LYP/aug-cc-pVTZ";
    static real th_toler=170,ph_toler=5;
    static real sigma=0;
    static int seed=1993;
    static gmx_bool bZero=FALSE,bForceFit=FALSE,bCompress=TRUE;
    t_pargs pa[] = 
    {
        { "-sort",   FALSE, etENUM, {sort},
          "Key to sort the final data file on." },
        { "-qm",     FALSE, etBOOL, {&iQM},
          "Use QM data for optimizing the empirical polarizabilities as well." },
        { "-lot",    FALSE, etSTR,  {&lot},
          "Use this method and level of theory" },
        { "-mindata", FALSE, etINT, {&mindata},
          "Minimum number of data points to optimize a polarizability value" },
        { "-sigma",  FALSE, etREAL, {&sigma},
          "Assumed standard deviation (relative) in the experimental data points" },
        { "-zero",    FALSE, etBOOL, {&bZero},
          "Use a zero term (like in Bosque and Sales)" },
        { "-force",   FALSE, etBOOL, {&bForceFit},
          "Reset all polarizablities to zero in order to re-optimize based on a gentop.dat file with previous values" },
        { "-th_toler",FALSE, etREAL, {&th_toler},
          "If bond angles are larger than this value the group will be treated as a linear one and a virtual site will be created to keep the group linear" },
        { "-ph_toler", FALSE, etREAL, {&ph_toler},
          "If dihedral angles are less than this (in absolute value) the atoms will be treated as a planar group with an improper dihedral being added to keep the group planar" },
        { "-compress", FALSE, etBOOL, {&bCompress},
          "Compress output XML files" }
    };
    FILE  *fp,*gp;
    int    i,alg,np,nspoel,nbosque,nhandbook,ntot,nqm,nalexandria_atypes,eMP,*qm_count;
    char   **fns;
    int    nfiles;
    char   *ref,*molname[2],*prop_name;
    int    cur = 0;
#define prev (1-cur)
    gmx_molprop_t   *mp=NULL;
    gmx_atomprop_t  ap;
    gmx_poldata_t   pd;
    output_env_t    oenv;
    gmx_molselect_t gms;
    int npa = sizeof(pa)/sizeof(pa[0]);
    CopyRight(stdout,argv[0]);
    
    parse_common_args(&argc,argv,PCA_NOEXIT_ON_ARGS,NFILE,fnm,
                      npa,pa,sizeof(desc)/sizeof(desc[0]),desc,
                      0,NULL,&oenv);
    ap = gmx_atomprop_init();
    if ((pd = gmx_poldata_read(opt2fn_null("-di",NFILE,fnm),ap)) == NULL)
        gmx_fatal(FARGS,"Can not read the force field information. File missing or incorrect.");
    nfiles = opt2fns(&fns,"-f",NFILE,fnm);
    mp = merge_xml(nfiles,fns,NULL,
                   NULL,(char *)"double_dip.dat",&np,ap,pd,TRUE,TRUE,th_toler,ph_toler);
    for(i=0; (i<np); i++) 
        gmx_molprop_check_consistency(mp[i]);

    gms = gmx_molselect_init(opt2fn("-sel",NFILE,fnm));
    nalexandria_atypes = decompose_frag(debug,0,pd,np,mp,iQM,lot,mindata,gms,sigma,bZero,bForceFit);

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
    if (alg != -1)
    {
        gmx_molprop_sort(np,mp,alg,ap,NULL);
    }
    gmx_poldata_write(opt2fn("-do",NFILE,fnm),pd,ap,bCompress);
    gmx_molprops_write(opt2fn("-o",NFILE,fnm),np,mp,bCompress);
    thanx(stdout);
  
    return 0;
}
