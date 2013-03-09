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
    int ndih,nimp;
    t_dih *dih,*imp;
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
    qsort(b->imp,b->nimp,sizeof(b->imp[0]),dihcmp);
}

void add_bond(FILE *fplog,char *molname,t_bonds *b,char *a1,char *a2,
              double blen,double spacing)
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
    if (NULL != fplog) {
        fprintf(fplog,"%s bond-%s-%s %g\n",molname,a1,a2,blen);
    }
}

void add_angle(FILE *fplog,char *molname,t_bonds *b,char *a1,char *a2,char *a3,double angle,double spacing)
{
    int i,index;
  
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
    if (NULL != fplog) {
        fprintf(fplog,"%s angle-%s-%s-%s %g\n",molname,a1,a2,a3,angle);
    }
}

void add_dih(FILE *fplog,char *molname,t_bonds *b,char *a1,char *a2,char *a3,char *a4,
             double angle,double spacing,int egd)
{
    int i,*nd,index;
    t_dih **ddd;
    
    if (angle < 0)
        angle += 360;
    if (egd == egdPDIHS) {
        ddd = &b->dih;
        nd = &b->ndih;
    }
    else {
        ddd = &b->imp;
        nd = &b->nimp;
        while (angle > 176)
            angle -= 180;
    }
    index = gmx_nint(angle/spacing);
    if (index < 0)
        index = 0;
    for(i=0; (i<(*nd)); i++) {
        if (((strcmp(a1,(*ddd)[i].a1) == 0) && (strcmp(a2,(*ddd)[i].a2) == 0) && 
             (strcmp(a3,(*ddd)[i].a3) == 0) && (strcmp(a4,(*ddd)[i].a4) == 0)) ||
            ((strcmp(a1,(*ddd)[i].a4) == 0) && (strcmp(a2,(*ddd)[i].a3) == 0) && 
             (strcmp(a3,(*ddd)[i].a2) == 0) && (strcmp(a4,(*ddd)[i].a1) == 0))) {
            break;
        }
    }
    if (i == (*nd)) {
        (*nd)++;
        srenew((*ddd),(*nd));
        if (strcmp(a1,a4) < 0) 
        {
            (*ddd)[i].a1     = strdup(a1);
            (*ddd)[i].a2     = strdup(a2);
            (*ddd)[i].a3     = strdup(a3);
            (*ddd)[i].a4     = strdup(a4);
        }
        else 
        {
            (*ddd)[i].a4     = strdup(a1);
            (*ddd)[i].a3     = strdup(a2);
            (*ddd)[i].a2     = strdup(a3);
            (*ddd)[i].a1     = strdup(a4);
        }
        (*ddd)[i].nhisto = (int) (360/spacing) + 1;
        snew((*ddd)[i].histo,(*ddd)[i].nhisto);
        (*ddd)[i].lsq = gmx_stats_init();
    }
    gmx_stats_add_point((*ddd)[i].lsq,0,angle,0,0);
    (*ddd)[i].histo[index]++;
    if (NULL != fplog) {
        fprintf(fplog,"%s %s-%s-%s-%s-%s %g\n",molname,(egd==egdPDIHS) ? "dih" : "imp",
                a1,a2,a3,a4,angle);
    }
}

static int imax(int a,int b)
{
    if (a > b)
        return a;
    else
        return b;
}

void lo_dump_histo(char *fn,char *xaxis,output_env_t oenv,
                   int n,int histo[],double spacing)
{
    FILE *fp;
    int j,j0,j1;
    double sum;
  
    for(j0=0; (j0<n) && (histo[j0] == 0); j0++)
        ;
    j0 = imax(j0-1,0);
    for(j1=n-1; (j1>0) && (histo[j1] == 0); j1--)
        ;
    j1 = imax(j1+1,n-1);
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
    int  i;
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
    for(i=0; (i<b->nimp); i++) {
        sprintf(buf,"imp-%s-%s-%s-%s.xvg",
                b->imp[i].a1,b->imp[i].a2,b->imp[i].a3,b->imp[i].a4);
        lo_dump_histo(buf,"Improper angle (deg.)",oenv,
                      b->imp[i].nhisto,b->imp[i].histo,aspacing);
    }
}

void update_pd(FILE *fp,t_bonds *b,gmx_poldata_t pd,gmx_atomprop_t aps,
               real Dm,real beta,real kt,real kp)
{
    int i,N;
    real av,sig;
    char pbuf[256];
    double bondorder;
        
    gmx_poldata_set_length_unit(pd,"pm");
    for(i=0; (i<b->nbond); i++) {
        gmx_stats_get_average(b->bond[i].lsq,&av);
        gmx_stats_get_sigma(b->bond[i].lsq,&sig);
        gmx_stats_get_npoints(b->bond[i].lsq,&N);
        sprintf(pbuf,"%g  %g",Dm,beta);
        bondorder = 1;
        gmx_poldata_add_bond(pd,b->bond[i].a1,b->bond[i].a2,av,sig,N,bondorder,pbuf);
        fprintf(fp,"bond-%s-%s len %g sigma %g (pm) N = %d%s\n",
                b->bond[i].a1,b->bond[i].a2,av,sig,N,
                (sig > 1.5) ? " WARNING" : "");
                
    }
    for(i=0; (i<b->nangle); i++) {
        gmx_stats_get_average(b->angle[i].lsq,&av);
        gmx_stats_get_sigma(b->angle[i].lsq,&sig);
        gmx_stats_get_npoints(b->angle[i].lsq,&N);
        sprintf(pbuf,"%g",kt);
        gmx_poldata_add_angle(pd,b->angle[i].a1,b->angle[i].a2,
                              b->angle[i].a3,av,sig,N,pbuf);
        fprintf(fp,"angle-%s-%s-%s angle %g sigma %g (deg) N = %d%s\n",
                b->angle[i].a1,b->angle[i].a2,b->angle[i].a3,av,sig,N,
                (sig > 3) ? " WARNING" : "");
    }
    for(i=0; (i<b->ndih); i++) {
        gmx_stats_get_average(b->dih[i].lsq,&av);
        gmx_stats_get_sigma(b->dih[i].lsq,&sig);
        gmx_stats_get_npoints(b->dih[i].lsq,&N);
        sprintf(pbuf,"%g",kp);
        gmx_poldata_add_dihedral(pd,egdPDIHS,
                                 b->dih[i].a1,b->dih[i].a2,
                                 b->dih[i].a3,b->dih[i].a4,av,sig,N,pbuf);
        fprintf(fp,"dihedral-%s-%s-%s-%s angle %g sigma %g (deg)\n",
                b->dih[i].a1,b->dih[i].a2,b->dih[i].a3,b->dih[i].a4,av,sig);
    }
    for(i=0; (i<b->nimp); i++) {
        gmx_stats_get_average(b->imp[i].lsq,&av);
        gmx_stats_get_sigma(b->imp[i].lsq,&sig);
        gmx_stats_get_npoints(b->imp[i].lsq,&N);
        sprintf(pbuf,"%g",kp);
        gmx_poldata_add_dihedral(pd,egdIDIHS,
                                 b->imp[i].a1,b->imp[i].a2,
                                 b->imp[i].a3,b->imp[i].a4,av,sig,N,pbuf);
        fprintf(fp,"improper-%s-%s-%s-%s angle %g sigma %g (deg)\n",
                b->imp[i].a1,b->imp[i].a2,b->imp[i].a3,b->imp[i].a4,av,sig);
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
    static real Dm = 400,kt = 400,kp = 5, beta = 20;
    static real J0_0=5,Chi0_0=1,w_0=5,hfac=0,rDecrZeta=-1;
    static real J0_1=30,Chi0_1=30,w_1=50,epsr=1;
    static real fc_mu=1,fc_bound=1,fc_quad=1,fc_charge=0,fc_esp=0;
    static real th_toler=170,ph_toler=5,dip_toler=0.5;
    static char *opt_elem = NULL,*const_elem=NULL,*fixchi=(char *)"H";
    static char *lot = (char *)"B3LYP/aug-cc-pVTZ";
    static char *qgen[] = { NULL,(char *)"AXp", (char *)"AXs", (char *)"AXg", NULL };
    t_pargs pa[] = {
        { "-lot",    FALSE, etSTR,  {&lot},
          "Use this method and level of theory when selecting coordinates and charges" },
        { "-Dm",    FALSE, etREAL, {&Dm},
          "Dissociation energy (kJ/mol)" },
        { "-beta",    FALSE, etREAL, {&beta},
          "Steepness of the Morse potential (1/nm)" },
        { "-kt",    FALSE, etREAL, {&kt},
          "Angle force constant (kJ/mol/rad^2)" },
        { "-kp",    FALSE, etREAL, {&kp},
          "Dihedral angle force constant (kJ/mol/rad^2)" },
        { "-histo", FALSE, etBOOL, {&bHisto},
          "Print (hundreds of) xvg files containing histograms for bonds, angles and dihedrals" },
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML file" }
    };
    t_moldip  *md;
    FILE      *fp;
    int       iModel;
    t_commrec *cr;
    gmx_molselect_t gms;
    time_t    my_t;
    char      pukestr[STRLEN];
    t_bonds   *b;
    int i,j,ai,aj,ak,al;
    char *cai,*caj,*cak,*cal;
    rvec dx,dx2,r_ij,r_kj,r_kl,mm,nn;
    t_pbc pbc;
    int t1,t2,t3;
    int ftb,fta,ftd,fti;
    matrix box;
    real sign;
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
                TRUE,TRUE,TRUE,watoms,FALSE);
                        
#define ATP(ii) get_atomtype_name(md->mymol[i].topology->atoms.atom[ii].type,md->mymol[i].atype)
    ftb = gmx_poldata_get_bond_ftype(md->pd);
    fta = gmx_poldata_get_angle_ftype(md->pd);
    ftd = gmx_poldata_get_dihedral_ftype(md->pd,egdPDIHS);
    fti = gmx_poldata_get_dihedral_ftype(md->pd,egdIDIHS);
    snew(b,1);
    for(i=0; (i<md->nmol); i++) {
        molname = md->mymol[i].molname;
        for(j=0; (j<md->mymol[i].ltop->idef.il[ftb].nr); j+=interaction_function[ftb].nratoms+1) {
            ai = md->mymol[i].ltop->idef.il[ftb].iatoms[j+1];
            aj = md->mymol[i].ltop->idef.il[ftb].iatoms[j+2];
            rvec_sub(md->mymol[i].x[ai],md->mymol[i].x[aj],dx);
            cai = ATP(ai);
            caj = ATP(aj);
            add_bond(fp,molname,b,cai,caj,1000*norm(dx),bspacing);
        }
        for(j=0; (j<md->mymol[i].ltop->idef.il[fta].nr); j+=interaction_function[fta].nratoms+1) {
            ai = md->mymol[i].ltop->idef.il[fta].iatoms[j+1];
            aj = md->mymol[i].ltop->idef.il[fta].iatoms[j+2];
            ak = md->mymol[i].ltop->idef.il[fta].iatoms[j+3];
            rvec_sub(md->mymol[i].x[ai],md->mymol[i].x[aj],dx);
            rvec_sub(md->mymol[i].x[ak],md->mymol[i].x[aj],dx2);
            ang = RAD2DEG*gmx_angle(dx,dx2);
            cai = ATP(ai);
            caj = ATP(aj);
            cak = ATP(ak);
            add_angle(fp,molname,b,cai,caj,cak,ang,aspacing);
        }
        clear_mat(box);
        set_pbc(&pbc,epbcNONE,box);
        
        for(j=0; (j<md->mymol[i].ltop->idef.il[ftd].nr); j+=interaction_function[ftd].nratoms+1) {
            ai = md->mymol[i].ltop->idef.il[ftd].iatoms[j+1];
            aj = md->mymol[i].ltop->idef.il[ftd].iatoms[j+2];
            ak = md->mymol[i].ltop->idef.il[ftd].iatoms[j+3];
            al = md->mymol[i].ltop->idef.il[ftd].iatoms[j+4];
            ang = RAD2DEG*dih_angle(md->mymol[i].x[ai],md->mymol[i].x[aj],
                                    md->mymol[i].x[ak],md->mymol[i].x[al],
                                    &pbc,r_ij,r_kj,r_kl,mm,nn, /* out */
                                    &sign,&t1,&t2,&t3);
            cai = ATP(ai);
            caj = ATP(aj);
            cak = ATP(ak);
            cal = ATP(al);
            add_dih(fp,molname,b,cai,caj,cak,cal,ang,dspacing,egdPDIHS);
        }

        for(j=0; (j<md->mymol[i].ltop->idef.il[fti].nr); j+=interaction_function[fti].nratoms+1) {
            ai = md->mymol[i].ltop->idef.il[fti].iatoms[j+1];
            aj = md->mymol[i].ltop->idef.il[fti].iatoms[j+2];
            ak = md->mymol[i].ltop->idef.il[fti].iatoms[j+3];
            al = md->mymol[i].ltop->idef.il[fti].iatoms[j+4];
            ang = RAD2DEG*dih_angle(md->mymol[i].x[ai],md->mymol[i].x[aj],
                                    md->mymol[i].x[ak],md->mymol[i].x[al],
                                    &pbc,r_ij,r_kj,r_kl,mm,nn, /* out */
                                    &sign,&t1,&t2,&t3);
            cai = ATP(ai);
            caj = ATP(aj);
            cak = ATP(ak);
            cal = ATP(al);
            add_dih(fp,molname,b,cai,caj,cak,cal,ang,dspacing,egdIDIHS);
        }
    }
    sort_bonds(b);
    if (bHisto)
        dump_histo(b,bspacing,aspacing,dspacing,oenv);
    update_pd(fp,b,md->pd,md->atomprop,Dm,beta,kt,kp);
    
    gmx_poldata_write(opt2fn("-o",NFILE,fnm),md->pd,md->atomprop,compress);    
    
    printf("Extracted %d bondtypes, %d angletypes, %d dihedraltypes and %d impropertypes.\n",
           b->nbond,b->nangle,b->ndih,b->nimp);
    
    ffclose(fp);
    thanx(stdout);
    
#ifdef GMX_MPI
    if (gmx_mpi_initialized())
        gmx_finalize_par();
#endif

    return 0;
}
