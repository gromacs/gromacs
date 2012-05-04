/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: bastat.c,v 1.8 2009/04/12 21:24:26 spoel Exp $
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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "typedefs.h"
#include "physics.h"
#include "pdbio.h"
#include "smalloc.h"
#include "atomprop.h"
#include "macros.h"
#include "maths.h"
#include "vec.h"
#include "xvgr.h"
#include "futil.h"
#include "main.h"
#include "copyrite.h"
#include "statutil.h"
#include "mymol.h"
#include "gmx_statistics.h"

typedef struct {
    char *a1,*a2;
    int  nhisto;
    int  *histo;
    gmx_stats_t lsq;
} t_bond;

typedef struct {
    char *a1,*a2,*a3;
    int  nhisto;
    int  *histo;
    gmx_stats_t lsq;
} t_angle;

typedef struct {
    int nbond;
    t_bond *bond;
    int nangle;
    t_angle *angle;
} t_bonds;

int bondcmp(const void *a,const void *b)
{
    t_bond *ba = (t_bond *)a;
    t_bond *bb = (t_bond *)b;
    int d;
  
    if ((d = strcmp(ba->a1,bb->a1)) == 0)
        return strcmp(ba->a2,bb->a2);
    else
        return d;
}

int anglecmp(const void *a,const void *b)
{
    t_angle *ba = (t_angle *)a;
    t_angle *bb = (t_angle *)b;
    int d;
  
    if ((d = strcmp(ba->a1,bb->a1)) == 0)
        if ((d = strcmp(ba->a2,bb->a2)) == 0)
            return strcmp(ba->a3,bb->a3);
  
    return d;
}

void sort_bonds(t_bonds *b)
{
    qsort(b->bond,b->nbond,sizeof(b->bond[0]),bondcmp);
    qsort(b->angle,b->nangle,sizeof(b->angle[0]),anglecmp);
}

void add_bond(t_bonds *b,char *a1,char *a2,double blen,double spacing)
{
    int i,j,index;
  
    index = gmx_nint(blen/spacing);
    for(i=0; (i<b->nbond); i++) {
        if ((((strcmp(a1,b->bond[i].a1) == 0) && (strcmp(a2,b->bond[i].a2)) == 0)) ||
            (((strcmp(a2,b->bond[i].a1) == 0) && (strcmp(a1,b->bond[i].a2)) == 0))) {
            break;
        }
    }
    if (i == b->nbond) {
        b->nbond++;
        srenew(b->bond,b->nbond);
        if (strcmp(a1,a2) < 0) 
        {
            b->bond[i].a1     = strdup(a1);
            b->bond[i].a2     = strdup(a2);
        }
        else 
        {
            b->bond[i].a1     = strdup(a2);
            b->bond[i].a2     = strdup(a1);
        }
        b->bond[i].nhisto = 2*index+1;
        snew(b->bond[i].histo,b->bond[i].nhisto);
        b->bond[i].lsq = gmx_stats_init();
    }
    if (index >= b->bond[i].nhisto) {
        srenew(b->bond[i].histo,index+100);
        for(j=b->bond[i].nhisto; (j<index+100); j++)
            b->bond[i].histo[j] = 0;
        b->bond[i].nhisto = index+100;
    }
    gmx_stats_add_point(b->bond[i].lsq,0,blen,0,0);
    b->bond[i].histo[index]++;
}

void add_angle(t_bonds *b,char *a1,char *a2,char *a3,double angle,double spacing)
{
    int i,j,index;
  
    index = gmx_nint(angle/spacing);
    for(i=0; (i<b->nangle); i++) {
        if ((strcmp(a2,b->angle[i].a2) == 0) && 
            (((strcmp(a1,b->angle[i].a1) == 0) && (strcmp(a3,b->angle[i].a3) == 0)) ||
             ((strcmp(a3,b->angle[i].a1) == 0) && (strcmp(a1,b->angle[i].a3) == 0)))) {
            break;
        }
    }
    if (i == b->nangle) {
        b->nangle++;
        srenew(b->angle,b->nangle);
        b->angle[i].a2     = strdup(a2);
        if (strcmp(a1,a3) < 0) 
        {
            b->angle[i].a1     = strdup(a1);
            b->angle[i].a3     = strdup(a3);
        }
        else 
        {
            b->angle[i].a1     = strdup(a3);
            b->angle[i].a3     = strdup(a1);
        }
        b->angle[i].nhisto = (int) (180/spacing) + 1;
        snew(b->angle[i].histo,b->angle[i].nhisto);
        b->angle[i].lsq = gmx_stats_init();
    }
    gmx_stats_add_point(b->angle[i].lsq,0,angle,0,0);

    b->angle[i].histo[index]++;
}

void lo_dump_histo(FILE *fp,int n,int histo[],double spacing)
{
    int j,j0,j1;
    double sum;
  
    for(j0=0; (j0<n) && (histo[j0] == 0); j0++)
        ;
    j0 = max(j0-1,0);
    for(j1=n-1; (j1>0) && (histo[j1] == 0); j1--)
        ;
    j1 = max(j1+1,n-1);
    sum=0;
    for(j=j0; (j<=j1); j++) 
        sum+=histo[j];
    if (sum == 0) {
        fprintf(stderr,"Empty histogram\n");
        exit(1);
    }
    for(j=j0; (j<=j1); j++) 
        fprintf(fp,"%g  %g\n",spacing*j,histo[j]/sum);
}

void dump_histo(t_bonds *b,double bspacing,double aspacing,output_env_t oenv)
{
    int  i,j,j0,j1;
    FILE *fp;
    char buf[256];
  
    for(i=0; (i<b->nbond); i++) {
        sprintf(buf,"bond-%s-%s.xvg",b->bond[i].a1,b->bond[i].a2);
        fp = xvgropen(buf,buf,"Distance (pm)","P (a.u.)",oenv);
        lo_dump_histo(fp,b->bond[i].nhisto,b->bond[i].histo,bspacing);
        fclose(fp);
    }
    for(i=0; (i<b->nangle); i++) {
        sprintf(buf,"angle-%s-%s-%s.xvg",
                b->angle[i].a1,b->angle[i].a2,b->angle[i].a3);
        fp = xvgropen(buf,buf,"Angle (deg.)","P (a.u.)",oenv);
        lo_dump_histo(fp,b->angle[i].nhisto,b->angle[i].histo,aspacing);
        fclose(fp);
    }
}

void dump_table(t_bonds *b,double bspacing,double aspacing,
                gmx_atomprop_t aps)
{
    int i,j,k,N;
    real x,y,av,sig;
    FILE *fp;
  
    fp = fopen("gt_bonds.xml","w");
    fprintf(fp,"  <gt_bonds length_unit=\"pm\">\n");
    for(i=0; (i<b->nbond); i++) {
        for(j=0; (j<b->bond[i].nhisto); j++) {
            x = j*bspacing;
            y = b->bond[i].histo[j];
        }
        gmx_stats_get_average(b->bond[i].lsq,&av);
        gmx_stats_get_sigma(b->bond[i].lsq,&sig);
        gmx_stats_get_npoints(b->bond[i].lsq,&N);
        
        fprintf(fp,"    <gt_bond atom1=\"%s\" atom2=\"%s\" length=\"%.1f\" stddev=\%.1f\" npoints=\"%d\"/>\n",
                b->bond[i].a1,b->bond[i].a2,av,sig,N);
    }
    fprintf(fp,"  </gt_bonds>\n");
    fprintf(fp,"  <gt_angles angle_unit=\"degree\">\n");
    for(i=0; (i<b->nangle); i++) {
        for(j=0; (j<b->angle[i].nhisto); j++) {
            x = j*aspacing;
            y = b->angle[i].histo[j];
        }
        gmx_stats_get_average(b->angle[i].lsq,&av);
        gmx_stats_get_sigma(b->angle[i].lsq,&sig);
        gmx_stats_get_npoints(b->angle[i].lsq,&N);
        
        fprintf(fp,"    <gt_angle atom1=\"%s\" atom2=\"%s\" atom3=\"%s\" angle=\"%.1f\" stddev=\"%.1f\" npoints=\"%d\"/>\n",
                b->angle[i].a1,b->angle[i].a2,b->angle[i].a3,av,sig,N);
    }
    fprintf(fp,"  </gt_angles>\n");
    fclose(fp);
}

int main(int argc,char *argv[])
{
    static const char *desc[] = {
        "bastat read a series of molecules and corresponding experimental",
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
        { efDAT, "-o", "bastat",    ffWRITE },
        { efDAT, "-sel", "molselect",ffREAD },
        { efLOG, "-g", "bastat",    ffWRITE }
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
    gmx_molselect_t gms;
    time_t    my_t;
    char      pukestr[STRLEN];
    t_bonds   *b;
    FILE *in;
    int i,j,l,ai,aj,ak;
    char *cai,*caj,*cak;
    t_atoms atoms;
    char title[2560];
    rvec *x,dx,dx2;
    matrix box;
    gmx_conect conect;
    gmx_atomprop_t aps;
    double ang;
    double bspacing = 1; /* pm */
    double aspacing = 0.1; /* degree */
    output_env_t oenv = NULL;
    
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
                TRUE,TRUE,TRUE,2,watoms,FALSE);
                    
#define ATP(ii) get_atomtype_name(md->mymol[i].atoms->atom[ii].type,md->mymol[i].atype)
    snew(b,1);
    for(i=0; (i<md->nmol); i++) {
        for(j=0; (j<md->mymol[i].ltop->idef.il[F_BONDS].nr); j+=interaction_function[F_BONDS].nratoms+1) {
            ai = md->mymol[i].ltop->idef.il[F_BONDS].iatoms[j+1];
            aj = md->mymol[i].ltop->idef.il[F_BONDS].iatoms[j+2];
            rvec_sub(md->mymol[i].x[ai],md->mymol[i].x[aj],dx);
            cai = ATP(ai);
            caj = ATP(aj);
            add_bond(b,cai,caj,1000*norm(dx),bspacing);
        }
        for(j=0; (j<md->mymol[i].ltop->idef.il[F_ANGLES].nr); j+=interaction_function[F_ANGLES].nratoms+1) {
            ai = md->mymol[i].ltop->idef.il[F_ANGLES].iatoms[j+1];
            aj = md->mymol[i].ltop->idef.il[F_ANGLES].iatoms[j+2];
            ak = md->mymol[i].ltop->idef.il[F_ANGLES].iatoms[j+3];
            rvec_sub(md->mymol[i].x[ai],md->mymol[i].x[aj],dx);
            rvec_sub(md->mymol[i].x[ak],md->mymol[i].x[aj],dx2);
            ang = RAD2DEG*acos(cos_angle(dx,dx2));
            cai = ATP(ai);
            caj = ATP(aj);
            cak = ATP(ak);
            add_angle(b,cai,caj,cak,ang,aspacing);
        }
    }
    sort_bonds(b);
    dump_histo(b,bspacing,aspacing,oenv);
    dump_table(b,bspacing,aspacing,aps);
  
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
