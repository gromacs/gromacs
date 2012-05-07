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
#include "pbc.h"
#include "smalloc.h"
#include "bondf.h"
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
#include "poldata_xml.h"

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
    char *a1,*a2,*a3,*a4;
    int  nhisto;
    int  *histo;
    gmx_stats_t lsq;
} t_dih;

typedef struct {
    int nbond;
    t_bond *bond;
    int nangle;
    t_angle *angle;
    int ndih;
    t_dih *dih;
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

int dihcmp(const void *a,const void *b)
{
    t_dih *ba = (t_dih *)a;
    t_dih *bb = (t_dih *)b;
    int d;
  
    if ((d = strcmp(ba->a1,bb->a1)) == 0)
        if ((d = strcmp(ba->a2,bb->a2)) == 0)
            if ((d = strcmp(ba->a3,bb->a3)) == 0)
                return strcmp(ba->a4,bb->a4);
  
    return d;
}

void sort_bonds(t_bonds *b)
{
    qsort(b->bond,b->nbond,sizeof(b->bond[0]),bondcmp);
    qsort(b->angle,b->nangle,sizeof(b->angle[0]),anglecmp);
    qsort(b->dih,b->ndih,sizeof(b->dih[0]),dihcmp);
}

void add_bond(char *molname,t_bonds *b,char *a1,char *a2,double blen,double spacing)
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
    if (NULL != debug) {
        fprintf(debug,"%s bond %s-%s %g\n",molname,a1,a2,blen);
    }
}

void add_angle(char *molname,t_bonds *b,char *a1,char *a2,char *a3,double angle,double spacing)
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
    if (NULL != debug) {
        fprintf(debug,"%s angle %s-%s-%s %g\n",molname,a1,a2,a3,angle);
    }
}

void add_dih(char *molname,t_bonds *b,char *a1,char *a2,char *a3,char *a4,
             double angle,double spacing)
{
    int i,j,index;
    
    if (angle < 0)
        angle += 360;
    index = gmx_nint(angle/spacing);
    for(i=0; (i<b->ndih); i++) {
        if (((strcmp(a1,b->dih[i].a1) == 0) && (strcmp(a2,b->dih[i].a2) == 0) && 
             (strcmp(a3,b->dih[i].a3) == 0) && (strcmp(a4,b->dih[i].a4) == 0)) ||
            ((strcmp(a1,b->dih[i].a4) == 0) && (strcmp(a2,b->dih[i].a3) == 0) && 
             (strcmp(a3,b->dih[i].a2) == 0) && (strcmp(a4,b->dih[i].a1) == 0))) {
            break;
        }
    }
    if (i == b->ndih) {
        b->ndih++;
        srenew(b->dih,b->ndih);
        if (strcmp(a1,a4) < 0) 
        {
            b->dih[i].a1     = strdup(a1);
            b->dih[i].a2     = strdup(a2);
            b->dih[i].a3     = strdup(a3);
            b->dih[i].a4     = strdup(a4);
        }
        else 
        {
            b->dih[i].a4     = strdup(a1);
            b->dih[i].a3     = strdup(a2);
            b->dih[i].a2     = strdup(a3);
            b->dih[i].a1     = strdup(a4);
        }
        b->dih[i].nhisto = (int) (360/spacing) + 1;
        snew(b->dih[i].histo,b->dih[i].nhisto);
        b->dih[i].lsq = gmx_stats_init();
    }
    gmx_stats_add_point(b->dih[i].lsq,0,angle,0,0);
    b->dih[i].histo[index]++;
    if (NULL != debug) {
        fprintf(debug,"%s dihedrals %s-%s-%s-%s %g\n",molname,a1,a2,a3,a4,angle);
    }
}

void lo_dump_histo(char *fn,char *xaxis,output_env_t oenv,
                   int n,int histo[],double spacing)
{
    FILE *fp;
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
    if (sum > 0) {
        fp = xvgropen(fn,fn,xaxis,"P (a.u.)",oenv);
        for(j=j0; (j<=j1); j++) 
            fprintf(fp,"%g  %g\n",spacing*j,histo[j]/sum);
        fclose(fp);
    }
}

void dump_histo(t_bonds *b,double bspacing,double aspacing,double dspacing,
                output_env_t oenv)
{
    int  i,j,j0,j1;
    char buf[256];
  
    for(i=0; (i<b->nbond); i++) {
        sprintf(buf,"bond-%s-%s.xvg",b->bond[i].a1,b->bond[i].a2);
        lo_dump_histo(buf,"Distance (pm)",oenv,
                      b->bond[i].nhisto,b->bond[i].histo,bspacing);
    }
    for(i=0; (i<b->nangle); i++) {
        sprintf(buf,"angle-%s-%s-%s.xvg",
                b->angle[i].a1,b->angle[i].a2,b->angle[i].a3);
        lo_dump_histo(buf,"Angle (deg.)",oenv,
                      b->angle[i].nhisto,b->angle[i].histo,aspacing);
    }
    for(i=0; (i<b->ndih); i++) {
        sprintf(buf,"dih-%s-%s-%s-%s.xvg",
                b->dih[i].a1,b->dih[i].a2,b->dih[i].a3,b->dih[i].a4);
        lo_dump_histo(buf,"Dihedral angle (deg.)",oenv,
                      b->dih[i].nhisto,b->dih[i].histo,aspacing);
    }
}

void update_pd(FILE *fp,t_bonds *b,gmx_poldata_t pd,gmx_atomprop_t aps)
{
    int i,N;
    real av,sig;
    char pbuf[256];
    
    gmx_poldata_set_length_unit(pd,"pm");
    for(i=0; (i<b->nbond); i++) {
        gmx_stats_get_average(b->bond[i].lsq,&av);
        gmx_stats_get_sigma(b->bond[i].lsq,&sig);
        gmx_stats_get_npoints(b->bond[i].lsq,&N);
        sprintf(pbuf,"%g",av);
        gmx_poldata_add_gt_bond(pd,b->bond[i].a1,b->bond[i].a2,av,sig,1.0,pbuf);
        fprintf(fp,"bond %s-%s len %g sigma %g (pm)\n",
                b->bond[i].a1,b->bond[i].a2,av,sig);
    }
    for(i=0; (i<b->nangle); i++) {
        gmx_stats_get_average(b->angle[i].lsq,&av);
        gmx_stats_get_sigma(b->angle[i].lsq,&sig);
        gmx_stats_get_npoints(b->angle[i].lsq,&N);
        sprintf(pbuf,"%g",av);
        gmx_poldata_add_gt_angle(pd,b->angle[i].a1,b->angle[i].a2,
                                 b->angle[i].a3,av,sig,pbuf);
        fprintf(fp,"angle %s-%s-%s angle %g sigma %g (deg)\n",
                b->angle[i].a1,b->angle[i].a2,b->angle[i].a3,av,sig);
    }
    for(i=0; (i<b->ndih); i++) {
        gmx_stats_get_average(b->dih[i].lsq,&av);
        gmx_stats_get_sigma(b->dih[i].lsq,&sig);
        gmx_stats_get_npoints(b->dih[i].lsq,&N);
        sprintf(pbuf,"%g",av);
        gmx_poldata_add_gt_dihedral(pd,b->dih[i].a1,b->dih[i].a2,
                                    b->dih[i].a3,b->dih[i].a4,av,sig,pbuf);
        fprintf(fp,"dihedral %s-%s-%s-%s angle %g sigma %g (deg)\n",
                b->dih[i].a1,b->dih[i].a2,b->dih[i].a3,b->dih[i].a4,av,sig);
    }
}

int main(int argc,char *argv[])
{
    static const char *desc[] = {
        "bastat read a series of molecules and extracts average geometries from",
        "those. First atomtypes are determined and then bond-lengths, bond-angles",
        "and dihedral angles are extracted. The results are stored in a gentop.dat file."
    };
  
    t_filenm fnm[] = {
        { efDAT, "-f", "allmols",    ffREAD  },
        { efDAT, "-d", "gentop",     ffOPTRD },
        { efDAT, "-o", "bastat",    ffWRITE },
        { efDAT, "-sel", "molselect",ffREAD },
        { efLOG, "-g", "bastat",    ffWRITE }
    };
#define NFILE asize(fnm)
    static int  compress=0;
    static gmx_bool bHisto=FALSE,bZero=TRUE,bWeighted=TRUE,bOptHfac=FALSE,bQM=FALSE,bCharged=TRUE,bGaussianBug=TRUE,bPol=FALSE,bFitZeta=TRUE;
    int minimum_data = 3;
    real watoms = 1;
    static real J0_0=5,Chi0_0=1,w_0=5,step=0.01,hfac=0,rDecrZeta=-1;
    static real J0_1=30,Chi0_1=30,w_1=50,epsr=1;
    static real fc_mu=1,fc_bound=1,fc_quad=1,fc_charge=0,fc_esp=0;
    static real th_toler=170,ph_toler=5,dip_toler=0.5,quad_toler=5,q_toler=0.25;
    static char *opt_elem = NULL,*const_elem=NULL,*fixchi=(char *)"H";
    static char *lot = (char *)"B3LYP/aug-cc-pVTZ";
    static char *qgen[] = { NULL,(char *)"AXp", (char *)"AXs", (char *)"AXg", NULL };
    static int  nthreads=0; /* set to determine # of threads automatically */
    t_pargs pa[] = {
        { "-histo", FALSE, etBOOL, {&bHisto},
          "Print (hundreds of) xvg files containing histograms for bonds, angles and dihedrals" },
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
    int i,j,l,ai,aj,ak,al;
    char *cai,*caj,*cak,*cal;
    t_atoms atoms;
    char title[2560];
    rvec *x,dx,dx2,r_ij,r_kj,r_kl,mm,nn;
    t_pbc pbc;
    int t1,t2,t3;
    matrix box;
    real sign;
    gmx_conect conect;
    gmx_atomprop_t aps;
    double ang;
    double bspacing = 1; /* pm */
    double aspacing = 0.5; /* degree */
    double dspacing = 1; /* degree */
    output_env_t oenv = NULL;
    char *molname;
    
    cr = init_par(&argc,&argv);

    CopyRight(stdout,argv[0]);

    parse_common_args(&argc,argv,PCA_CAN_VIEW,
                      NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv);

#ifndef GMX_THREAD_MPI
    nthreads=1;
#endif
    if (qgen[0]) 
        iModel = name2eemtype(qgen[0]);
    else
        iModel = eqgAXg;
    
    fp = ffopen(opt2fn("-g",NFILE,fnm),"w");

    time(&my_t);
    fprintf(fp,"# This file was created %s",ctime(&my_t));
    fprintf(fp,"# by the following command:\n# %s\n#\n",command_line());
    fprintf(fp,"# %s is part of G R O M A C S:\n#\n",ShortProgram());
    bromacs(pukestr,99);
    fprintf(fp,"# %s\n#\n",pukestr);
        
    gms = gmx_molselect_init(opt2fn("-sel",NFILE,fnm));
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
        molname = md->mymol[i].molname;
        for(j=0; (j<md->mymol[i].ltop->idef.il[F_BONDS].nr); j+=interaction_function[F_BONDS].nratoms+1) {
            ai = md->mymol[i].ltop->idef.il[F_BONDS].iatoms[j+1];
            aj = md->mymol[i].ltop->idef.il[F_BONDS].iatoms[j+2];
            rvec_sub(md->mymol[i].x[ai],md->mymol[i].x[aj],dx);
            cai = ATP(ai);
            caj = ATP(aj);
            add_bond(molname,b,cai,caj,1000*norm(dx),bspacing);
        }
        for(j=0; (j<md->mymol[i].ltop->idef.il[F_ANGLES].nr); j+=interaction_function[F_ANGLES].nratoms+1) {
            ai = md->mymol[i].ltop->idef.il[F_ANGLES].iatoms[j+1];
            aj = md->mymol[i].ltop->idef.il[F_ANGLES].iatoms[j+2];
            ak = md->mymol[i].ltop->idef.il[F_ANGLES].iatoms[j+3];
            rvec_sub(md->mymol[i].x[ai],md->mymol[i].x[aj],dx);
            rvec_sub(md->mymol[i].x[ak],md->mymol[i].x[aj],dx2);
            ang = RAD2DEG*gmx_angle(dx,dx2);
            cai = ATP(ai);
            caj = ATP(aj);
            cak = ATP(ak);
            add_angle(molname,b,cai,caj,cak,ang,aspacing);
        }
        clear_mat(box);
        set_pbc(&pbc,epbcNONE,box);
        
        for(j=0; (j<md->mymol[i].ltop->idef.il[F_PDIHS].nr); j+=interaction_function[F_PDIHS].nratoms+1) {
            ai = md->mymol[i].ltop->idef.il[F_PDIHS].iatoms[j+1];
            aj = md->mymol[i].ltop->idef.il[F_PDIHS].iatoms[j+2];
            ak = md->mymol[i].ltop->idef.il[F_PDIHS].iatoms[j+3];
            al = md->mymol[i].ltop->idef.il[F_PDIHS].iatoms[j+4];
            ang = RAD2DEG*dih_angle(md->mymol[i].x[ai],md->mymol[i].x[aj],
                                    md->mymol[i].x[ak],md->mymol[i].x[al],
                                    &pbc,r_ij,r_kj,r_kl,mm,nn, /* out */
                                    &sign,&t1,&t2,&t3);
            cai = ATP(ai);
            caj = ATP(aj);
            cak = ATP(ak);
            cal = ATP(al);
            add_dih(molname,b,cai,caj,cak,cal,ang,dspacing);
        }
    }
    sort_bonds(b);
    if (bHisto)
        dump_histo(b,bspacing,aspacing,dspacing,oenv);
    update_pd(fp,b,md->pd,md->atomprop);
    
    gmx_poldata_write(opt2fn("-o",NFILE,fnm),md->pd,md->atomprop,compress);    
    
    printf("Extracted %d bondtypes, %d angletypes and %d dihedraltypes.\n",
           b->nbond,b->nangle,b->ndih);
    
    ffclose(fp);
    thanx(stdout);
    
#ifdef GMX_MPI
    if (gmx_parallel_env_initialized())
        gmx_finalize();
#endif

    return 0;
}
