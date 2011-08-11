/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2008
 * David van der Spoel, Erik Lindahl, Berk Hess, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */

/* The source code in this file should be thread-safe. 
 Please keep it that way. */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <time.h>

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif


#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
/* _chsize_s */
#include <io.h>
#include <sys/locking.h>
#endif


#include "filenm.h"
#include "names.h"
#include "typedefs.h"
#include "smalloc.h"
#include "gmxfio.h"
#include "xdrf.h"
#include "statutil.h"
#include "txtdump.h"
#include "vec.h"
#include "network.h"
#include "gmx_random.h"
#include "checkpoint.h"
#include "futil.h"
#include "string2.h"
#include <fcntl.h>


#ifdef GMX_FAHCORE
#include "corewrap.h"
#endif


/* Portable version of ctime_r implemented in src/gmxlib/string2.c, but we do not want it declared in public installed headers */
char *
gmx_ctime_r(const time_t *clock,char *buf, int n);


#define CPT_MAGIC1 171817
#define CPT_MAGIC2 171819

/* cpt_version should normally only be changed
 * when the header of footer format changes.
 * The state data format itself is backward and forward compatible.
 * But old code can not read a new entry that is present in the file
 * (but can read a new format when new entries are not present).
 */
static const int cpt_version = 12;


const char *est_names[estNR]=
{
    "FE-lambda",
    "box", "box-rel", "box-v", "pres_prev",
    "nosehoover-xi", "thermostat-integral",
    "x", "v", "SDx", "CGp", "LD-rng", "LD-rng-i",
    "disre_initf", "disre_rm3tav",
    "orire_initf", "orire_Dtav",
    "svir_prev", "nosehoover-vxi", "v_eta", "vol0", "nhpres_xi", "nhpres_vxi", "fvir_prev",
};

enum { eeksEKIN_N, eeksEKINH, eeksDEKINDL, eeksMVCOS, eeksEKINF, eeksEKINO, eeksEKINSCALEF, eeksEKINSCALEH, eeksVSCALE, eeksEKINTOTAL, eeksNR };

const char *eeks_names[eeksNR]=
{
    "Ekin_n", "Ekinh", "dEkindlambda", "mv_cos",
    "Ekinf", "Ekinh_old", "EkinScaleF_NHC", "EkinScaleH_NHC","Vscale_NHC","Ekin_Total"
};

enum { eenhENERGY_N, eenhENERGY_AVER, eenhENERGY_SUM, eenhENERGY_NSUM,
       eenhENERGY_SUM_SIM, eenhENERGY_NSUM_SIM,
       eenhENERGY_NSTEPS, eenhENERGY_NSTEPS_SIM, 
       eenhENERGY_DELTA_H_NN,
       eenhENERGY_DELTA_H_LIST, 
       eenhENERGY_DELTA_H_STARTTIME, 
       eenhENERGY_DELTA_H_STARTLAMBDA, 
       eenhNR };

const char *eenh_names[eenhNR]=
{
    "energy_n", "energy_aver", "energy_sum", "energy_nsum",
    "energy_sum_sim", "energy_nsum_sim",
    "energy_nsteps", "energy_nsteps_sim", 
    "energy_delta_h_nn",
    "energy_delta_h_list", 
    "energy_delta_h_start_time", 
    "energy_delta_h_start_lambda"
};



#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
static int
gmx_wintruncate(const char *filename, __int64 size)
{
#ifdef GMX_FAHCORE
    /*we do this elsewhere*/
    return 0;
#else
    FILE *fp;
    int   rc;
    
    fp=fopen(filename,"rb+");
    
    if(fp==NULL)
    {
        return -1;
    }
    
    return _chsize_s( fileno(fp), size);
#endif
}
#endif


enum { ecprREAL, ecprRVEC, ecprMATRIX };

static const char *st_names(int cptp,int ecpt)
{
    switch (cptp)
    {
    case 0: return est_names [ecpt]; break;
    case 1: return eeks_names[ecpt]; break;
    case 2: return eenh_names[ecpt]; break;
    }

    return NULL;
}

static void cp_warning(FILE *fp)
{
    fprintf(fp,"\nWARNING: Checkpoint file is corrupted or truncated\n\n");
}

static void cp_error()
{
    gmx_fatal(FARGS,"Checkpoint file corrupted/truncated, or maybe you are out of quota?");
}

static void do_cpt_string_err(XDR *xd,gmx_bool bRead,const char *desc,char **s,FILE *list)
{
#define CPTSTRLEN 1024
    bool_t res=0;
    
    if (bRead)
    {
        snew(*s,CPTSTRLEN);
    }
    res = xdr_string(xd,s,CPTSTRLEN);
    if (res == 0)
    {
        cp_error();
    }
    if (list)
    {
        fprintf(list,"%s = %s\n",desc,*s);
        sfree(*s);
    }
}

static int do_cpt_int(XDR *xd,const char *desc,int *i,FILE *list)
{
    bool_t res=0;
    
    res = xdr_int(xd,i);
    if (res == 0)
    {
        return -1;
    }
    if (list)
    {
        fprintf(list,"%s = %d\n",desc,*i);
    }
    return 0;
}

static int do_cpt_u_chars(XDR *xd,const char *desc,int n,unsigned char *i,FILE *list)
{
    bool_t res=1;
    int j;
    if (list)
    {
        fprintf(list,"%s = ",desc);
    }
    for (j=0; j<n && res; j++)
    {
        res &= xdr_u_char(xd,&i[j]);
        if (list)
        {
            fprintf(list,"%02x",i[j]);
        }
    }
    if (list)
    {
        fprintf(list,"\n");
    }
    if (res == 0)
    {
        return -1;
    }

    return 0;
}

static void do_cpt_int_err(XDR *xd,const char *desc,int *i,FILE *list)
{
    if (do_cpt_int(xd,desc,i,list) < 0)
    {
        cp_error();
    }
}

static void do_cpt_step_err(XDR *xd,const char *desc,gmx_large_int_t *i,FILE *list)
{
    bool_t res=0;
    char   buf[STEPSTRSIZE];

    res = xdr_gmx_large_int(xd,i,"reading checkpoint file");
    if (res == 0)
    {
        cp_error();
    }
    if (list)
    {
        fprintf(list,"%s = %s\n",desc,gmx_step_str(*i,buf));
    }
}

static void do_cpt_double_err(XDR *xd,const char *desc,double *f,FILE *list)
{
    bool_t res=0;
    
    res = xdr_double(xd,f);
    if (res == 0)
    {
        cp_error();
    }
    if (list)
    {
        fprintf(list,"%s = %f\n",desc,*f);
    }
}

/* If nval >= 0, nval is used; on read this should match the passed value.
 * If nval n<0, *nptr is used; on read the value is stored in nptr
 */
static int do_cpte_reals_low(XDR *xd,int cptp,int ecpt,int sflags,
                             int nval,int *nptr,real **v,
                             FILE *list,int erealtype)
{
    bool_t res=0;
#ifndef GMX_DOUBLE
    int  dtc=xdr_datatype_float; 
#else
    int  dtc=xdr_datatype_double;
#endif
    real *vp,*va=NULL;
    float  *vf;
    double *vd;
    int  nf,dt,i;
    
    if (list == NULL)
    {
        if (nval >= 0)
        {
            nf = nval;
        }
        else
        {
        if (nptr == NULL)
        {
            gmx_incons("*ntpr=NULL in do_cpte_reals_low");
        }
        nf = *nptr;
        }
    }
    res = xdr_int(xd,&nf);
    if (res == 0)
    {
        return -1;
    }
    if (list == NULL)
    {
        if (nval >= 0)
        {
            if (nf != nval)
            {
                gmx_fatal(FARGS,"Count mismatch for state entry %s, code count is %d, file count is %d\n",st_names(cptp,ecpt),nval,nf);
            }
        }
        else
        {
            *nptr = nf;
        }
    }
    dt = dtc;
    res = xdr_int(xd,&dt);
    if (res == 0)
    {
        return -1;
    }
    if (dt != dtc)
    {
        fprintf(stderr,"Precision mismatch for state entry %s, code precision is %s, file precision is %s\n",
                st_names(cptp,ecpt),xdr_datatype_names[dtc],
                xdr_datatype_names[dt]);
    }
    if (list || !(sflags & (1<<ecpt)))
    {
        snew(va,nf);
        vp = va;
    }
    else
    {
        if (*v == NULL)
        {
            snew(*v,nf);
        }
        vp = *v;
    }
    if (dt == xdr_datatype_float)
    {
        if (dtc == xdr_datatype_float)
        {
            vf = (float *)vp;
        }
        else
        {
            snew(vf,nf);
        }
        res = xdr_vector(xd,(char *)vf,nf,
                         (unsigned int)sizeof(float),(xdrproc_t)xdr_float);
        if (res == 0)
        {
            return -1;
        }
        if (dtc != xdr_datatype_float)
        {
            for(i=0; i<nf; i++)
            {
                vp[i] = vf[i];
            }
            sfree(vf);
        }
    }
    else
    {
        if (dtc == xdr_datatype_double)
        {
            vd = (double *)vp;
        }
        else
        {
            snew(vd,nf);
        }
        res = xdr_vector(xd,(char *)vd,nf,
                         (unsigned int)sizeof(double),(xdrproc_t)xdr_double);
        if (res == 0)
        {
            return -1;
        }
        if (dtc != xdr_datatype_double)
        {
            for(i=0; i<nf; i++)
            {
                vp[i] = vd[i];
            }
            sfree(vd);
        }
    }
    
    if (list)
    {
        switch (erealtype)
        {
        case ecprREAL:
            pr_reals(list,0,st_names(cptp,ecpt),vp,nf);
            break;
        case ecprRVEC:
            pr_rvecs(list,0,st_names(cptp,ecpt),(rvec *)vp,nf/3);
            break;
        default:
            gmx_incons("Unknown checkpoint real type");
        }
    }
    if (va)
    {
        sfree(va);
    }

    return 0;
}


/* This function stores n along with the reals for reading,
 * but on reading it assumes that n matches the value in the checkpoint file,
 * a fatal error is generated when this is not the case.
 */
static int do_cpte_reals(XDR *xd,int cptp,int ecpt,int sflags,
                         int n,real **v,FILE *list)
{
    return do_cpte_reals_low(xd,cptp,ecpt,sflags,n,NULL,v,list,ecprREAL);
}

/* This function does the same as do_cpte_reals,
 * except that on reading it ignores the passed value of *n
 * and stored the value read from the checkpoint file in *n.
 */
static int do_cpte_n_reals(XDR *xd,int cptp,int ecpt,int sflags,
                           int *n,real **v,FILE *list)
{
    return do_cpte_reals_low(xd,cptp,ecpt,sflags,-1,n,v,list,ecprREAL);
}

static int do_cpte_real(XDR *xd,int cptp,int ecpt,int sflags,
                        real *r,FILE *list)
{
    int n;

    return do_cpte_reals_low(xd,cptp,ecpt,sflags,1,NULL,&r,list,ecprREAL);
}

static int do_cpte_ints(XDR *xd,int cptp,int ecpt,int sflags,
                        int n,int **v,FILE *list)
{
    bool_t res=0;
    int  dtc=xdr_datatype_int;
    int *vp,*va=NULL;
    int  nf,dt,i;
    
    nf = n;
    res = xdr_int(xd,&nf);
    if (res == 0)
    {
        return -1;
    }
    if (list == NULL && v != NULL && nf != n)
    {
        gmx_fatal(FARGS,"Count mismatch for state entry %s, code count is %d, file count is %d\n",st_names(cptp,ecpt),n,nf);
    }
    dt = dtc;
    res = xdr_int(xd,&dt);
    if (res == 0)
    {
        return -1;
    }
    if (dt != dtc)
    {
        gmx_fatal(FARGS,"Type mismatch for state entry %s, code type is %s, file type is %s\n",
                  st_names(cptp,ecpt),xdr_datatype_names[dtc],
                  xdr_datatype_names[dt]);
    }
    if (list || !(sflags & (1<<ecpt)) || v == NULL)
    {
        snew(va,nf);
        vp = va;
    }
    else
    {
        if (*v == NULL)
        {
            snew(*v,nf);
        }
        vp = *v;
    }
    res = xdr_vector(xd,(char *)vp,nf,
                     (unsigned int)sizeof(int),(xdrproc_t)xdr_int);
    if (res == 0)
    {
        return -1;
    }
    if (list)
    {
        pr_ivec(list,0,st_names(cptp,ecpt),vp,nf,TRUE);
    }
    if (va)
    {
        sfree(va);
    }

    return 0;
}

static int do_cpte_int(XDR *xd,int cptp,int ecpt,int sflags,
                       int *i,FILE *list)
{
    return do_cpte_ints(xd,cptp,ecpt,sflags,1,&i,list);
}

static int do_cpte_doubles(XDR *xd,int cptp,int ecpt,int sflags,
                           int n,double **v,FILE *list)
{
    bool_t res=0;
    int  dtc=xdr_datatype_double;
    double *vp,*va=NULL;
    int  nf,dt,i;
    
    nf = n;
    res = xdr_int(xd,&nf);
    if (res == 0)
    {
        return -1;
    }
    if (list == NULL && nf != n)
    {
        gmx_fatal(FARGS,"Count mismatch for state entry %s, code count is %d, file count is %d\n",st_names(cptp,ecpt),n,nf);
    }
    dt = dtc;
    res = xdr_int(xd,&dt);
    if (res == 0)
    {
        return -1;
    }
    if (dt != dtc)
    {
        gmx_fatal(FARGS,"Precision mismatch for state entry %s, code precision is %s, file precision is %s\n",
                  st_names(cptp,ecpt),xdr_datatype_names[dtc],
                  xdr_datatype_names[dt]);
    }
    if (list || !(sflags & (1<<ecpt)))
    {
        snew(va,nf);
        vp = va;
    }
    else
    {
        if (*v == NULL)
        {
            snew(*v,nf);
        }
        vp = *v;
    }
    res = xdr_vector(xd,(char *)vp,nf,
                     (unsigned int)sizeof(double),(xdrproc_t)xdr_double);
    if (res == 0)
    {
        return -1;
    }
    if (list)
    {
        pr_doubles(list,0,st_names(cptp,ecpt),vp,nf);
    }
    if (va)
    {
        sfree(va);
    }

    return 0;
}

static int do_cpte_double(XDR *xd,int cptp,int ecpt,int sflags,
                          double *r,FILE *list)
{
    return do_cpte_doubles(xd,cptp,ecpt,sflags,1,&r,list);
}


static int do_cpte_rvecs(XDR *xd,int cptp,int ecpt,int sflags,
                         int n,rvec **v,FILE *list)
{
    int n3;

    return do_cpte_reals_low(xd,cptp,ecpt,sflags,
                             n*DIM,NULL,(real **)v,list,ecprRVEC);
}

static int do_cpte_matrix(XDR *xd,int cptp,int ecpt,int sflags,
                          matrix v,FILE *list)
{
    real *vr;
    real ret;

    vr = (real *)&(v[0][0]);
    ret = do_cpte_reals_low(xd,cptp,ecpt,sflags,
                            DIM*DIM,NULL,&vr,NULL,ecprMATRIX);
    
    if (list && ret == 0)
    {
        pr_rvecs(list,0,st_names(cptp,ecpt),v,DIM);
    }
    
    return ret;
}

static int do_cpte_matrices(XDR *xd,int cptp,int ecpt,int sflags,
                            int n,matrix **v,FILE *list)
{
    bool_t res=0;
    matrix *vp,*va=NULL;
    real *vr;
    int  nf,i,j,k;
    int  ret;

    nf = n;
    res = xdr_int(xd,&nf);
    if (res == 0)
    {
        return -1;
    }
    if (list == NULL && nf != n)
    {
        gmx_fatal(FARGS,"Count mismatch for state entry %s, code count is %d, file count is %d\n",st_names(cptp,ecpt),n,nf);
    }
    if (list || !(sflags & (1<<ecpt)))
    {
        snew(va,nf);
        vp = va;
    }
    else
    {
        if (*v == NULL)
        {
            snew(*v,nf);
        }
        vp = *v;
    }
    snew(vr,nf*DIM*DIM);
    for(i=0; i<nf; i++)
    {
        for(j=0; j<DIM; j++)
        {
            for(k=0; k<DIM; k++)
            {
                vr[(i*DIM+j)*DIM+k] = vp[i][j][k];
            }
        }
    }
    ret = do_cpte_reals_low(xd,cptp,ecpt,sflags,
                            nf*DIM*DIM,NULL,&vr,NULL,ecprMATRIX);
    for(i=0; i<nf; i++)
    {
        for(j=0; j<DIM; j++)
        {
            for(k=0; k<DIM; k++)
            {
                vp[i][j][k] = vr[(i*DIM+j)*DIM+k];
            }
        }
    }
    sfree(vr);
    
    if (list && ret == 0)
    {
        for(i=0; i<nf; i++)
        {
            pr_rvecs(list,0,st_names(cptp,ecpt),vp[i],DIM);
        }
    }
    if (va)
    {
        sfree(va);
    }
    
    return ret;
}

static void do_cpt_header(XDR *xd,gmx_bool bRead,int *file_version,
                          char **version,char **btime,char **buser,char **bmach,
                          char **fprog,char **ftime,
                          int *eIntegrator,int *simulation_part,
                          gmx_large_int_t *step,double *t,
                          int *nnodes,int *dd_nc,int *npme,
                          int *natoms,int *ngtc, int *nnhpres, int *nhchainlength,
                          int *flags_state,int *flags_eks,int *flags_enh,
                          FILE *list)
{
    bool_t res=0;
    int  magic;
    int  idum=0;
    int  i;
    char *fhost;

    if (bRead)
    {
        magic = -1;
    }
    else
    {
        magic = CPT_MAGIC1;
    }
    res = xdr_int(xd,&magic);
    if (res == 0)
    {
        gmx_fatal(FARGS,"The checkpoint file is empty/corrupted, or maybe you are out of quota?");
    }
    if (magic != CPT_MAGIC1)
    {
        gmx_fatal(FARGS,"Start of file magic number mismatch, checkpoint file has %d, should be %d\n"
                  "The checkpoint file is corrupted or not a checkpoint file",
                  magic,CPT_MAGIC1);
    }
    if (!bRead)
    {
        snew(fhost,255);
#ifdef HAVE_UNISTD_H
        if (gethostname(fhost,255) != 0)
        {
            sprintf(fhost,"unknown");
        }
#else
        sprintf(fhost,"unknown");
#endif  
    }
    do_cpt_string_err(xd,bRead,"GROMACS version"           ,version,list);
    do_cpt_string_err(xd,bRead,"GROMACS build time"        ,btime,list);
    do_cpt_string_err(xd,bRead,"GROMACS build user"        ,buser,list);
    do_cpt_string_err(xd,bRead,"GROMACS build machine"     ,bmach,list);
    do_cpt_string_err(xd,bRead,"generating program"        ,fprog,list);
    do_cpt_string_err(xd,bRead,"generation time"           ,ftime,list);
    *file_version = cpt_version;
    do_cpt_int_err(xd,"checkpoint file version",file_version,list);
    if (*file_version > cpt_version)
    {
        gmx_fatal(FARGS,"Attempting to read a checkpoint file of version %d with code of version %d\n",*file_version,cpt_version);
    }
    if (*file_version >= 12)
    {
        do_cpt_string_err(xd,bRead,"generating host"           ,&fhost,list);
        if (list == NULL)
        {
            sfree(fhost);
        }
    }
    do_cpt_int_err(xd,"#atoms"            ,natoms     ,list);
    do_cpt_int_err(xd,"#T-coupling groups",ngtc       ,list);
    if (*file_version >= 10) 
    {
        do_cpt_int_err(xd,"#Nose-Hoover T-chains",nhchainlength,list);
    }
    else
    {
        *nhchainlength = 1;
    }
    if (*file_version >= 11)
    {
        do_cpt_int_err(xd,"#Nose-Hoover T-chains for barostat ",nnhpres,list);
    }
    else
    {
        *nnhpres = 0;
    }
    do_cpt_int_err(xd,"integrator"        ,eIntegrator,list);
	if (*file_version >= 3)
	{
		do_cpt_int_err(xd,"simulation part #", simulation_part,list);
	}
	else
	{
		*simulation_part = 1;
	}
    if (*file_version >= 5)
    {
        do_cpt_step_err(xd,"step"         ,step       ,list);
    }
    else
    {
        do_cpt_int_err(xd,"step"          ,&idum      ,list);
        *step = idum;
    }
    do_cpt_double_err(xd,"t"              ,t          ,list);
    do_cpt_int_err(xd,"#PP-nodes"         ,nnodes     ,list);
    idum = 1;
    do_cpt_int_err(xd,"dd_nc[x]",dd_nc ? &(dd_nc[0]) : &idum,list);
    do_cpt_int_err(xd,"dd_nc[y]",dd_nc ? &(dd_nc[1]) : &idum,list);
    do_cpt_int_err(xd,"dd_nc[z]",dd_nc ? &(dd_nc[2]) : &idum,list);
    do_cpt_int_err(xd,"#PME-only nodes",npme,list);
    do_cpt_int_err(xd,"state flags",flags_state,list);
	if (*file_version >= 4)
    {
        do_cpt_int_err(xd,"ekin data flags",flags_eks,list);
        do_cpt_int_err(xd,"energy history flags",flags_enh,list);
    }
    else
    {
        *flags_eks  = 0;
        *flags_enh   = (*flags_state >> (estORIRE_DTAV+1));
        *flags_state = (*flags_state & ~((1<<(estORIRE_DTAV+1)) |
                                         (1<<(estORIRE_DTAV+2)) |
                                         (1<<(estORIRE_DTAV+3))));
    }
}

static int do_cpt_footer(XDR *xd,gmx_bool bRead,int file_version)
{
    bool_t res=0;
    int  magic;
    
    if (file_version >= 2)
    {
        magic = CPT_MAGIC2;
        res = xdr_int(xd,&magic);
        if (res == 0)
        {
            cp_error();
        }
        if (magic != CPT_MAGIC2)
        {
            return -1;
        }
    }

    return 0;
}

static int do_cpt_state(XDR *xd,gmx_bool bRead,
                        int fflags,t_state *state,
                        gmx_bool bReadRNG,FILE *list)
{
    int  sflags;
    int  **rng_p,**rngi_p;
    int  i;
    int  ret;
    int  nnht,nnhtp;

    ret = 0;
    
    nnht = state->nhchainlength*state->ngtc;
    nnhtp = state->nhchainlength*state->nnhpres;

    if (bReadRNG)
    {
        rng_p  = (int **)&state->ld_rng;
        rngi_p = &state->ld_rngi;
    }
    else
    {
        /* Do not read the RNG data */
        rng_p  = NULL;
        rngi_p = NULL;
    }

    sflags = state->flags;
    for(i=0; (i<estNR && ret == 0); i++)
    {
        if (fflags & (1<<i))
        {
            switch (i)
            {
            case estLAMBDA:  ret = do_cpte_real(xd,0,i,sflags,&state->lambda,list); break;
            case estBOX:     ret = do_cpte_matrix(xd,0,i,sflags,state->box,list); break;
            case estBOX_REL: ret = do_cpte_matrix(xd,0,i,sflags,state->box_rel,list); break;
            case estBOXV:    ret = do_cpte_matrix(xd,0,i,sflags,state->boxv,list); break;
            case estPRES_PREV: ret = do_cpte_matrix(xd,0,i,sflags,state->pres_prev,list); break;
            case estSVIR_PREV:  ret = do_cpte_matrix(xd,0,i,sflags,state->svir_prev,list); break;
            case estFVIR_PREV:  ret = do_cpte_matrix(xd,0,i,sflags,state->fvir_prev,list); break;
            case estNH_XI:   ret = do_cpte_doubles(xd,0,i,sflags,nnht,&state->nosehoover_xi,list); break;
            case estNH_VXI:  ret = do_cpte_doubles(xd,0,i,sflags,nnht,&state->nosehoover_vxi,list); break;
            case estNHPRES_XI:   ret = do_cpte_doubles(xd,0,i,sflags,nnhtp,&state->nhpres_xi,list); break;
            case estNHPRES_VXI:  ret = do_cpte_doubles(xd,0,i,sflags,nnhtp,&state->nhpres_vxi,list); break;
            case estTC_INT:  ret = do_cpte_doubles(xd,0,i,sflags,state->ngtc,&state->therm_integral,list); break;
            case estVETA:    ret = do_cpte_real(xd,0,i,sflags,&state->veta,list); break;
            case estVOL0:    ret = do_cpte_real(xd,0,i,sflags,&state->vol0,list); break;
            case estX:       ret = do_cpte_rvecs(xd,0,i,sflags,state->natoms,&state->x,list); break;
            case estV:       ret = do_cpte_rvecs(xd,0,i,sflags,state->natoms,&state->v,list); break;
            case estSDX:     ret = do_cpte_rvecs(xd,0,i,sflags,state->natoms,&state->sd_X,list); break;
            case estLD_RNG:  ret = do_cpte_ints(xd,0,i,sflags,state->nrng,rng_p,list); break;
            case estLD_RNGI: ret = do_cpte_ints(xd,0,i,sflags,state->nrngi,rngi_p,list); break;
            case estDISRE_INITF:  ret = do_cpte_real (xd,0,i,sflags,&state->hist.disre_initf,list); break;
            case estDISRE_RM3TAV: ret = do_cpte_n_reals(xd,0,i,sflags,&state->hist.ndisrepairs,&state->hist.disre_rm3tav,list); break;
            case estORIRE_INITF:  ret = do_cpte_real (xd,0,i,sflags,&state->hist.orire_initf,list); break;
            case estORIRE_DTAV:   ret = do_cpte_n_reals(xd,0,i,sflags,&state->hist.norire_Dtav,&state->hist.orire_Dtav,list); break;
            default:
                gmx_fatal(FARGS,"Unknown state entry %d\n"
                          "You are probably reading a new checkpoint file with old code",i);
            }
        }
    }
    
    return ret;
}

static int do_cpt_ekinstate(XDR *xd,gmx_bool bRead,
                            int fflags,ekinstate_t *ekins,
                            FILE *list)
{
    int  i;
    int  ret;

    ret = 0;

    for(i=0; (i<eeksNR && ret == 0); i++)
    {
        if (fflags & (1<<i))
        {
            switch (i)
            {
                
			case eeksEKIN_N:     ret = do_cpte_int(xd,1,i,fflags,&ekins->ekin_n,list); break;
			case eeksEKINH :     ret = do_cpte_matrices(xd,1,i,fflags,ekins->ekin_n,&ekins->ekinh,list); break;
			case eeksEKINF:      ret = do_cpte_matrices(xd,1,i,fflags,ekins->ekin_n,&ekins->ekinf,list); break;
			case eeksEKINO:      ret = do_cpte_matrices(xd,1,i,fflags,ekins->ekin_n,&ekins->ekinh_old,list); break;
            case eeksEKINTOTAL:  ret = do_cpte_matrix(xd,1,i,fflags,ekins->ekin_total,list); break;
            case eeksEKINSCALEF: ret = do_cpte_doubles(xd,1,i,fflags,ekins->ekin_n,&ekins->ekinscalef_nhc,list); break;
            case eeksVSCALE:     ret = do_cpte_doubles(xd,1,i,fflags,ekins->ekin_n,&ekins->vscale_nhc,list); break;
            case eeksEKINSCALEH: ret = do_cpte_doubles(xd,1,i,fflags,ekins->ekin_n,&ekins->ekinscaleh_nhc,list); break;
 			case eeksDEKINDL :   ret = do_cpte_real(xd,1,i,fflags,&ekins->dekindl,list); break;
            case eeksMVCOS:      ret = do_cpte_real(xd,1,i,fflags,&ekins->mvcos,list); break;			
            default:
                gmx_fatal(FARGS,"Unknown ekin data state entry %d\n"
                          "You are probably reading a new checkpoint file with old code",i);
            }
        }
    }
    
    return ret;
}


static int do_cpt_enerhist(XDR *xd,gmx_bool bRead,
                           int fflags,energyhistory_t *enerhist,
                           FILE *list)
{
    int  i;
    int  j;
    int  ret;

    ret = 0;

    if (bRead)
    {
        enerhist->nsteps     = 0;
        enerhist->nsum       = 0;
        enerhist->nsteps_sim = 0;
        enerhist->nsum_sim   = 0;
        enerhist->dht        = NULL;

        if (fflags & (1<< eenhENERGY_DELTA_H_NN) )
        {
            snew(enerhist->dht,1);
            enerhist->dht->ndh = NULL;
            enerhist->dht->dh = NULL;
            enerhist->dht->start_lambda_set=FALSE;
        }
    }

    for(i=0; (i<eenhNR && ret == 0); i++)
    {
        if (fflags & (1<<i))
        {
            switch (i)
            {
                case eenhENERGY_N:     ret = do_cpte_int(xd,2,i,fflags,&enerhist->nener,list); break;
                case eenhENERGY_AVER:  ret = do_cpte_doubles(xd,2,i,fflags,enerhist->nener,&enerhist->ener_ave,list); break;
                case eenhENERGY_SUM:   ret = do_cpte_doubles(xd,2,i,fflags,enerhist->nener,&enerhist->ener_sum,list); break;
                case eenhENERGY_NSUM:  do_cpt_step_err(xd,eenh_names[i],&enerhist->nsum,list); break;
                case eenhENERGY_SUM_SIM: ret = do_cpte_doubles(xd,2,i,fflags,enerhist->nener,&enerhist->ener_sum_sim,list); break;
                case eenhENERGY_NSUM_SIM:   do_cpt_step_err(xd,eenh_names[i],&enerhist->nsum_sim,list); break;
                case eenhENERGY_NSTEPS:     do_cpt_step_err(xd,eenh_names[i],&enerhist->nsteps,list); break;
                case eenhENERGY_NSTEPS_SIM: do_cpt_step_err(xd,eenh_names[i],&enerhist->nsteps_sim,list); break;
                case eenhENERGY_DELTA_H_NN: do_cpt_int_err(xd,eenh_names[i], &(enerhist->dht->nndh), list); 
                    if (bRead) /* now allocate memory for it */
                    {
                        snew(enerhist->dht->dh, enerhist->dht->nndh);
                        snew(enerhist->dht->ndh, enerhist->dht->nndh);
                        for(j=0;j<enerhist->dht->nndh;j++)
                        {
                            enerhist->dht->ndh[j] = 0;
                            enerhist->dht->dh[j] = NULL;
                        }
                    }
                break;
                case eenhENERGY_DELTA_H_LIST: 
                    for(j=0;j<enerhist->dht->nndh;j++)
                    {
                        ret=do_cpte_n_reals(xd, 2, i, fflags, &enerhist->dht->ndh[j], &(enerhist->dht->dh[j]), list); 
                    }
                    break;
                case eenhENERGY_DELTA_H_STARTTIME: 
                    ret=do_cpte_double(xd, 2, i, fflags, &(enerhist->dht->start_time), list); break;
                case eenhENERGY_DELTA_H_STARTLAMBDA: 
                    enerhist->dht->start_lambda_set=TRUE;
                    ret=do_cpte_double(xd, 2, i, fflags, &(enerhist->dht->start_lambda), list); break;
                default:
                    gmx_fatal(FARGS,"Unknown energy history entry %d\n"
                              "You are probably reading a new checkpoint file with old code",i);
            }
        }
    }

    if ((fflags & (1<<eenhENERGY_SUM)) && !(fflags & (1<<eenhENERGY_SUM_SIM)))
    {
        /* Assume we have an old file format and copy sum to sum_sim */
        srenew(enerhist->ener_sum_sim,enerhist->nener);
        for(i=0; i<enerhist->nener; i++)
        {
            enerhist->ener_sum_sim[i] = enerhist->ener_sum[i];
        }
        fflags |= (1<<eenhENERGY_SUM_SIM);
    }
    
    if ( (fflags & (1<<eenhENERGY_NSUM)) &&
        !(fflags & (1<<eenhENERGY_NSTEPS)))
    {
        /* Assume we have an old file format and copy nsum to nsteps */
        enerhist->nsteps = enerhist->nsum;
        fflags |= (1<<eenhENERGY_NSTEPS);
    }
    if ( (fflags & (1<<eenhENERGY_NSUM_SIM)) &&
        !(fflags & (1<<eenhENERGY_NSTEPS_SIM)))
    {
        /* Assume we have an old file format and copy nsum to nsteps */
        enerhist->nsteps_sim = enerhist->nsum_sim;
        fflags |= (1<<eenhENERGY_NSTEPS_SIM);
    }

    return ret;
}

static int do_cpt_files(XDR *xd, gmx_bool bRead, 
                        gmx_file_position_t **p_outputfiles, int *nfiles, 
                        FILE *list, int file_version)
{
    int    i,j;
    gmx_off_t  offset;
    gmx_off_t  mask = 0xFFFFFFFFL;
    int    offset_high,offset_low;
    char   *buf;
    gmx_file_position_t *outputfiles;

    if (do_cpt_int(xd,"number of output files",nfiles,list) != 0)
    {
        return -1;
    }

    if(bRead)
    {
        snew(*p_outputfiles,*nfiles);
    }

    outputfiles = *p_outputfiles;

    for(i=0;i<*nfiles;i++)
    {
        /* 64-bit XDR numbers are not portable, so it is stored as separate high/low fractions */
        if(bRead)
        {
            do_cpt_string_err(xd,bRead,"output filename",&buf,list);
            strncpy(outputfiles[i].filename,buf,CPTSTRLEN-1);
            if(list==NULL)
            {
                sfree(buf);			
            }

            if (do_cpt_int(xd,"file_offset_high",&offset_high,list) != 0)
            {
                return -1;
            }
            if (do_cpt_int(xd,"file_offset_low",&offset_low,list) != 0)
            {
                return -1;
            }
#if (SIZEOF_GMX_OFF_T > 4)
            outputfiles[i].offset = ( ((gmx_off_t) offset_high) << 32 ) | ( (gmx_off_t) offset_low & mask );
#else
            outputfiles[i].offset = offset_low;
#endif
        }
        else
        {
            buf = outputfiles[i].filename;
            do_cpt_string_err(xd,bRead,"output filename",&buf,list);
            /* writing */
            offset      = outputfiles[i].offset;
            if (offset == -1)
            {
                offset_low  = -1;
                offset_high = -1;
            }
            else
            {
#if (SIZEOF_GMX_OFF_T > 4)
                offset_low  = (int) (offset & mask);
                offset_high = (int) ((offset >> 32) & mask);
#else
                offset_low  = offset;
                offset_high = 0;
#endif
            }
            if (do_cpt_int(xd,"file_offset_high",&offset_high,list) != 0)
            {
                return -1;
            }
            if (do_cpt_int(xd,"file_offset_low",&offset_low,list) != 0)
            {
                return -1;
            }
        }
        if (file_version >= 8)
        {
            if (do_cpt_int(xd,"file_checksum_size",&(outputfiles[i].chksum_size),
                           list) != 0)
            {
                return -1;
            }
            if (do_cpt_u_chars(xd,"file_checksum",16,outputfiles[i].chksum,list) != 0)
            {
                return -1;
            }
        } 
        else 
        {
            outputfiles[i].chksum_size = -1;
        }
    }
    return 0;
}


void write_checkpoint(const char *fn,gmx_bool bNumberAndKeep,
                      FILE *fplog,t_commrec *cr,
                      int eIntegrator,int simulation_part,
                      gmx_large_int_t step,double t,t_state *state)
{
    t_fileio *fp;
    int  file_version;
    char *version;
    char *btime;
    char *buser;
    char *bmach;
    char *fprog;
    char *fntemp; /* the temporary checkpoint file name */
    time_t now;
    char timebuf[STRLEN];
    int  nppnodes,npmenodes,flag_64bit;
    char buf[1024],suffix[5+STEPSTRSIZE],sbuf[STEPSTRSIZE];
    gmx_file_position_t *outputfiles;
    int  noutputfiles;
    char *ftime;
    int  flags_eks,flags_enh,i;
    t_fileio *ret;
		
    if (PAR(cr))
    {
        if (DOMAINDECOMP(cr))
        {
            nppnodes  = cr->dd->nnodes;
            npmenodes = cr->npmenodes;
        }
        else
        {
            nppnodes  = cr->nnodes;
            npmenodes = 0;
        }
    }
    else
    {
        nppnodes  = 1;
        npmenodes = 0;
    }

    /* make the new temporary filename */
    snew(fntemp, strlen(fn)+5+STEPSTRSIZE);
    strcpy(fntemp,fn);
    fntemp[strlen(fn) - strlen(ftp2ext(fn2ftp(fn))) - 1] = '\0';
    sprintf(suffix,"_%s%s","step",gmx_step_str(step,sbuf));
    strcat(fntemp,suffix);
    strcat(fntemp,fn+strlen(fn) - strlen(ftp2ext(fn2ftp(fn))) - 1);
   
    time(&now);
    gmx_ctime_r(&now,timebuf,STRLEN);

    if (fplog)
    { 
        fprintf(fplog,"Writing checkpoint, step %s at %s\n\n",
                gmx_step_str(step,buf),timebuf);
    }
    
    /* Get offsets for open files */
    gmx_fio_get_output_file_positions(&outputfiles, &noutputfiles);

    fp = gmx_fio_open(fntemp,"w");
	
    if (state->ekinstate.bUpToDate)
    {
        flags_eks =
            ((1<<eeksEKIN_N) | (1<<eeksEKINH) | (1<<eeksEKINF) | 
             (1<<eeksEKINO) | (1<<eeksEKINSCALEF) | (1<<eeksEKINSCALEH) | 
             (1<<eeksVSCALE) | (1<<eeksDEKINDL) | (1<<eeksMVCOS));
    }
    else
    {
        flags_eks = 0;
    }

    flags_enh = 0;
    if (state->enerhist.nsum > 0 || state->enerhist.nsum_sim > 0)
    {
        flags_enh |= (1<<eenhENERGY_N);
        if (state->enerhist.nsum > 0)
        {
            flags_enh |= ((1<<eenhENERGY_AVER) | (1<<eenhENERGY_SUM) |
                          (1<<eenhENERGY_NSTEPS) | (1<<eenhENERGY_NSUM));
        }
        if (state->enerhist.nsum_sim > 0)
        {
            flags_enh |= ((1<<eenhENERGY_SUM_SIM) | (1<<eenhENERGY_NSTEPS_SIM) |
                          (1<<eenhENERGY_NSUM_SIM));
        }
        if (state->enerhist.dht)
        {
            flags_enh |= ( (1<< eenhENERGY_DELTA_H_NN) |
                           (1<< eenhENERGY_DELTA_H_LIST) | 
                           (1<< eenhENERGY_DELTA_H_STARTTIME) |
                           (1<< eenhENERGY_DELTA_H_STARTLAMBDA) );
        }
    }

    
    version = gmx_strdup(VERSION);
    btime   = gmx_strdup(BUILD_TIME);
    buser   = gmx_strdup(BUILD_USER);
    bmach   = gmx_strdup(BUILD_MACHINE);
    fprog   = gmx_strdup(Program());

    ftime   = &(timebuf[0]);
    
    do_cpt_header(gmx_fio_getxdr(fp),FALSE,&file_version,
                  &version,&btime,&buser,&bmach,&fprog,&ftime,
                  &eIntegrator,&simulation_part,&step,&t,&nppnodes,
                  DOMAINDECOMP(cr) ? cr->dd->nc : NULL,&npmenodes,
                  &state->natoms,&state->ngtc,&state->nnhpres,
                  &state->nhchainlength, &state->flags,&flags_eks,&flags_enh,
                  NULL);
    
    sfree(version);
    sfree(btime);
    sfree(buser);
    sfree(bmach);
    sfree(fprog);

    if((do_cpt_state(gmx_fio_getxdr(fp),FALSE,state->flags,state,TRUE,NULL) < 0)        ||
       (do_cpt_ekinstate(gmx_fio_getxdr(fp),FALSE,flags_eks,&state->ekinstate,NULL) < 0)||
       (do_cpt_enerhist(gmx_fio_getxdr(fp),FALSE,flags_enh,&state->enerhist,NULL) < 0)  ||
       (do_cpt_files(gmx_fio_getxdr(fp),FALSE,&outputfiles,&noutputfiles,NULL,
                     file_version) < 0))
    {
        gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of quota?");
    }

    do_cpt_footer(gmx_fio_getxdr(fp),FALSE,file_version);

    /* we really, REALLY, want to make sure to physically write the checkpoint, 
       and all the files it depends on, out to disk. Because we've
       opened the checkpoint with gmx_fio_open(), it's in our list
       of open files.  */
    ret=gmx_fio_all_output_fsync();

    if (ret)
    {
        char buf[STRLEN];
        sprintf(buf,
                "Cannot fsync '%s'; maybe you are out of disk space or quota?",
                gmx_fio_getname(ret));

        if (getenv(GMX_IGNORE_FSYNC_FAILURE_ENV)==NULL)
        {
            gmx_file(buf);
        }
        else
        {
            gmx_warning(buf);
        }
    }

    if( gmx_fio_close(fp) != 0)
    {
        gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of quota?");
    }

    /* we don't move the checkpoint if the user specified they didn't want it,
       or if the fsyncs failed */
    if (!bNumberAndKeep && !ret)
    {
        if (gmx_fexist(fn))
        {
            /* Rename the previous checkpoint file */
            strcpy(buf,fn);
            buf[strlen(fn) - strlen(ftp2ext(fn2ftp(fn))) - 1] = '\0';
            strcat(buf,"_prev");
            strcat(buf,fn+strlen(fn) - strlen(ftp2ext(fn2ftp(fn))) - 1);
#ifndef GMX_FAHCORE
            /* we copy here so that if something goes wrong between now and
             * the rename below, there's always a state.cpt.
             * If renames are atomic (such as in POSIX systems),
             * this copying should be unneccesary.
             */
            gmx_file_copy(fn, buf, FALSE);
            /* We don't really care if this fails: 
             * there's already a new checkpoint.
             */
#else
            gmx_file_rename(fn, buf);
#endif
        }
        if (gmx_file_rename(fntemp, fn) != 0)
        {
            gmx_file("Cannot rename checkpoint file; maybe you are out of quota?");
        }
    }

    sfree(outputfiles);
    sfree(fntemp);

#ifdef GMX_FAHCORE
    /*code for alternate checkpointing scheme.  moved from top of loop over 
      steps */
    fcRequestCheckPoint();
    if ( fcCheckPointParallel( cr->nodeid, NULL,0) == 0 ) {
        gmx_fatal( 3,__FILE__,__LINE__, "Checkpoint error on step %d\n", step );
    }
#endif /* end GMX_FAHCORE block */
}

static void print_flag_mismatch(FILE *fplog,int sflags,int fflags)
{
    int i;
    
    fprintf(fplog,"\nState entry mismatch between the simulation and the checkpoint file\n");
    fprintf(fplog,"Entries which are not present in the checkpoint file will not be updated\n");
    fprintf(fplog,"  %24s    %11s    %11s\n","","simulation","checkpoint");
    for(i=0; i<estNR; i++)
    {
        if ((sflags & (1<<i)) || (fflags & (1<<i)))
        {
            fprintf(fplog,"  %24s    %11s    %11s\n",
                    est_names[i],
                    (sflags & (1<<i)) ? "  present  " : "not present",
                    (fflags & (1<<i)) ? "  present  " : "not present");
        }
    }
}

static void check_int(FILE *fplog,const char *type,int p,int f,gmx_bool *mm)
{
	FILE *fp = fplog ? fplog : stderr;

    if (p != f)
    {
		fprintf(fp,"  %s mismatch,\n",type);
		fprintf(fp,"    current program: %d\n",p);
		fprintf(fp,"    checkpoint file: %d\n",f);
		fprintf(fp,"\n");
        *mm = TRUE;
    }
}

static void check_string(FILE *fplog,const char *type,const char *p,
                         const char *f,gmx_bool *mm)
{
	FILE *fp = fplog ? fplog : stderr;
	
    if (strcmp(p,f) != 0)
    {
		fprintf(fp,"  %s mismatch,\n",type);
		fprintf(fp,"    current program: %s\n",p);
		fprintf(fp,"    checkpoint file: %s\n",f);
		fprintf(fp,"\n");
        *mm = TRUE;
    }
}

static void check_match(FILE *fplog,
                        char *version,
                        char *btime,char *buser,char *bmach,char *fprog,
                        t_commrec *cr,gmx_bool bPartDecomp,int npp_f,int npme_f,
                        ivec dd_nc,ivec dd_nc_f)
{
    int  npp;
    gmx_bool mm;
    
    mm = FALSE;
    
    check_string(fplog,"Version"      ,VERSION      ,version,&mm);
    check_string(fplog,"Build time"   ,BUILD_TIME   ,btime  ,&mm);
    check_string(fplog,"Build user"   ,BUILD_USER   ,buser  ,&mm);
    check_string(fplog,"Build machine",BUILD_MACHINE,bmach  ,&mm);
    check_string(fplog,"Program name" ,Program()    ,fprog  ,&mm);
    
    check_int   (fplog,"#nodes"       ,cr->nnodes   ,npp_f+npme_f ,&mm);
    if (bPartDecomp)
    {
        dd_nc[XX] = 1;
        dd_nc[YY] = 1;
        dd_nc[ZZ] = 1;
    }
    if (cr->nnodes > 1)
    {
        check_int (fplog,"#PME-nodes"  ,cr->npmenodes,npme_f     ,&mm);

        npp = cr->nnodes;
        if (cr->npmenodes >= 0)
        {
            npp -= cr->npmenodes;
        }
        if (npp == npp_f)
        {
            check_int (fplog,"#DD-cells[x]",dd_nc[XX]    ,dd_nc_f[XX],&mm);
            check_int (fplog,"#DD-cells[y]",dd_nc[YY]    ,dd_nc_f[YY],&mm);
            check_int (fplog,"#DD-cells[z]",dd_nc[ZZ]    ,dd_nc_f[ZZ],&mm);
        }
    }
    
    if (mm)
    {
		fprintf(stderr,
				"Gromacs binary or parallel settings not identical to previous run.\n"
				"Continuation is exact, but is not guaranteed to be binary identical%s.\n\n",
				fplog ? ",\n see the log file for details" : "");
		
        if (fplog)
        {
			fprintf(fplog,
					"Gromacs binary or parallel settings not identical to previous run.\n"
					"Continuation is exact, but is not guaranteed to be binary identical.\n\n");
		}
    }
}

static void read_checkpoint(const char *fn,FILE **pfplog,
                            t_commrec *cr,gmx_bool bPartDecomp,ivec dd_nc,
                            int eIntegrator,gmx_large_int_t *step,double *t,
                            t_state *state,gmx_bool *bReadRNG,gmx_bool *bReadEkin,
                            int *simulation_part,gmx_bool bAppendOutputFiles)
{
    t_fileio *fp;
    int  i,j,rc;
    int  file_version;
    char *version,*btime,*buser,*bmach,*fprog,*ftime;
	char filename[STRLEN],buf[STEPSTRSIZE];
    int  nppnodes,eIntegrator_f,nppnodes_f,npmenodes_f;
    ivec dd_nc_f;
    int  natoms,ngtc,nnhpres,nhchainlength,fflags,flags_eks,flags_enh;
    int  d;
    int  ret;
    gmx_file_position_t *outputfiles;
    int  nfiles;
    t_fileio *chksum_file;
    FILE* fplog = *pfplog;
    unsigned char digest[16];
#if !((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
    struct flock fl;  /* don't initialize here: the struct order is OS 
                         dependent! */
#endif

    const char *int_warn=
              "WARNING: The checkpoint file was generated with integrator %s,\n"
              "         while the simulation uses integrator %s\n\n";
    const char *sd_note=
        "NOTE: The checkpoint file was for %d nodes doing SD or BD,\n"
        "      while the simulation uses %d SD or BD nodes,\n"
        "      continuation will be exact, except for the random state\n\n";
    
#if !((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__) 
    fl.l_type=F_WRLCK;
    fl.l_whence=SEEK_SET;
    fl.l_start=0;
    fl.l_len=0;
    fl.l_pid=0;
#endif

    if (PARTDECOMP(cr))
    {
        gmx_fatal(FARGS,
                  "read_checkpoint not (yet) supported with particle decomposition");
    }
    
    fp = gmx_fio_open(fn,"r");
    do_cpt_header(gmx_fio_getxdr(fp),TRUE,&file_version,
                  &version,&btime,&buser,&bmach,&fprog,&ftime,
                  &eIntegrator_f,simulation_part,step,t,
                  &nppnodes_f,dd_nc_f,&npmenodes_f,
                  &natoms,&ngtc,&nnhpres,&nhchainlength,
                  &fflags,&flags_eks,&flags_enh,NULL);
    
    if (cr == NULL || MASTER(cr))
    {
        fprintf(stderr,"\nReading checkpoint file %s generated: %s\n\n",
                fn,ftime);
    }
	
	/* This will not be written if we do appending, since fplog is still NULL then */
    if (fplog)
    {
        fprintf(fplog,"\n");
        fprintf(fplog,"Reading checkpoint file %s\n",fn);
        fprintf(fplog,"  file generated by:     %s\n",fprog);  
        fprintf(fplog,"  file generated at:     %s\n",ftime);  
        fprintf(fplog,"  GROMACS build time:    %s\n",btime);  
        fprintf(fplog,"  GROMACS build user:    %s\n",buser);  
        fprintf(fplog,"  GROMACS build machine: %s\n",bmach);  
        fprintf(fplog,"  simulation part #:     %d\n",*simulation_part);
        fprintf(fplog,"  step:                  %s\n",gmx_step_str(*step,buf));
        fprintf(fplog,"  time:                  %f\n",*t);  
        fprintf(fplog,"\n");
    }
    
    if (natoms != state->natoms)
    {
        gmx_fatal(FARGS,"Checkpoint file is for a system of %d atoms, while the current system consists of %d atoms",natoms,state->natoms);
    }
    if (ngtc != state->ngtc)
    {
        gmx_fatal(FARGS,"Checkpoint file is for a system of %d T-coupling groups, while the current system consists of %d T-coupling groups",ngtc,state->ngtc);
    }
    if (nnhpres != state->nnhpres)
    {
        gmx_fatal(FARGS,"Checkpoint file is for a system of %d NH-pressure-coupling variables, while the current system consists of %d NH-pressure-coupling variables",nnhpres,state->nnhpres);
    }

    init_gtc_state(state,state->ngtc,state->nnhpres,nhchainlength); /* need to keep this here to keep the tpr format working */
    /* write over whatever was read; we use the number of Nose-Hoover chains from the checkpoint */
    
    if (eIntegrator_f != eIntegrator)
    {
        if (MASTER(cr))
        {
            fprintf(stderr,int_warn,EI(eIntegrator_f),EI(eIntegrator));
        }
		if(bAppendOutputFiles)
		{
			gmx_fatal(FARGS,
					  "Output file appending requested, but input/checkpoint integrators do not match.\n"
					  "Stopping the run to prevent you from ruining all your data...\n"
					  "If you _really_ know what you are doing, try with the -noappend option.\n");
		}
        if (fplog)
        {
            fprintf(fplog,int_warn,EI(eIntegrator_f),EI(eIntegrator));
        }
    }

    if (!PAR(cr))
    {
        nppnodes = 1;
        cr->npmenodes = 0;
    }
    else if (bPartDecomp)
    {
        nppnodes = cr->nnodes;
        cr->npmenodes = 0;
    }
    else if (cr->nnodes == nppnodes_f + npmenodes_f)
    {
        if (cr->npmenodes < 0)
        {
            cr->npmenodes = npmenodes_f;
        }
        nppnodes = cr->nnodes - cr->npmenodes;
        if (nppnodes == nppnodes_f)
        {
            for(d=0; d<DIM; d++)
            {
                if (dd_nc[d] == 0)
                {
                    dd_nc[d] = dd_nc_f[d];
                }
            }
        }
    }
    else
    {
        /* The number of PP nodes has not been set yet */
        nppnodes = -1;
    }

    if ((EI_SD(eIntegrator) || eIntegrator == eiBD) && nppnodes > 0)
    {
        /* Correct the RNG state size for the number of PP nodes.
         * Such assignments should all be moved to one central function.
         */
        state->nrng  = nppnodes*gmx_rng_n();
        state->nrngi = nppnodes;
    }
    
    *bReadRNG = TRUE;
    if (fflags != state->flags)
    {
		
        if (MASTER(cr))
        {
			if(bAppendOutputFiles)
			{
				gmx_fatal(FARGS,
						  "Output file appending requested, but input and checkpoint states are not identical.\n"
						  "Stopping the run to prevent you from ruining all your data...\n"
						  "You can try with the -noappend option, and get more info in the log file.\n");
			}
			
            if (getenv("GMX_ALLOW_CPT_MISMATCH") == NULL)
            {
                gmx_fatal(FARGS,"You seem to have switched ensemble, integrator, T and/or P-coupling algorithm between the cpt and tpr file. The recommended way of doing this is passing the cpt file to grompp (with option -t) instead of to mdrun. If you know what you are doing, you can override this error by setting the env.var. GMX_ALLOW_CPT_MISMATCH");
            }
            else
            {
                fprintf(stderr,
                        "WARNING: The checkpoint state entries do not match the simulation,\n"
                        "         see the log file for details\n\n");
            }
        }
		
		if(fplog)
		{
			print_flag_mismatch(fplog,state->flags,fflags);
		}
    }
    else
    {
        if ((EI_SD(eIntegrator) || eIntegrator == eiBD) &&
            nppnodes != nppnodes_f)
        {
            *bReadRNG = FALSE;
            if (MASTER(cr))
            {
                fprintf(stderr,sd_note,nppnodes_f,nppnodes);
            }
            if (fplog)
            {
                fprintf(fplog ,sd_note,nppnodes_f,nppnodes);
            }
        }
        if (MASTER(cr))
        {
            check_match(fplog,version,btime,buser,bmach,fprog,
                        cr,bPartDecomp,nppnodes_f,npmenodes_f,dd_nc,dd_nc_f);
        }
    }
    ret = do_cpt_state(gmx_fio_getxdr(fp),TRUE,fflags,state,*bReadRNG,NULL);
    if (ret)
    {
        cp_error();
    }
    ret = do_cpt_ekinstate(gmx_fio_getxdr(fp),TRUE,
                           flags_eks,&state->ekinstate,NULL);
    if (ret)
    {
        cp_error();
    }
    *bReadEkin = ((flags_eks & (1<<eeksEKINH)) || (flags_eks & (1<<eeksEKINF)) || (flags_eks & (1<<eeksEKINO)) ||
                  (flags_eks & (1<<eeksEKINSCALEF)) | (flags_eks & (1<<eeksEKINSCALEH)) | (flags_eks & (1<<eeksVSCALE)));
    
    ret = do_cpt_enerhist(gmx_fio_getxdr(fp),TRUE,
                          flags_enh,&state->enerhist,NULL);
    if (ret)
    {
        cp_error();
    }

    if (file_version < 6)
    {
        const char *warn="Reading checkpoint file in old format, assuming that the run that generated this file started at step 0, if this is not the case the averages stored in the energy file will be incorrect.";

        fprintf(stderr,"\nWARNING: %s\n\n",warn);
        if (fplog)
        {
            fprintf(fplog,"\nWARNING: %s\n\n",warn);
        }
        state->enerhist.nsum     = *step;
        state->enerhist.nsum_sim = *step;
    }

	ret = do_cpt_files(gmx_fio_getxdr(fp),TRUE,&outputfiles,&nfiles,NULL,file_version);
	if (ret)
	{
		cp_error();
	}
					   
    ret = do_cpt_footer(gmx_fio_getxdr(fp),TRUE,file_version);
    if (ret)
    {
        cp_error();
    }
    if( gmx_fio_close(fp) != 0)
	{
		gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of quota?");
	}
    
    sfree(fprog);
    sfree(ftime);
    sfree(btime);
    sfree(buser);
    sfree(bmach);
	
	/* If the user wants to append to output files,
     * we use the file pointer positions of the output files stored
     * in the checkpoint file and truncate the files such that any frames
     * written after the checkpoint time are removed.
     * All files are md5sum checked such that we can be sure that
     * we do not truncate other (maybe imprortant) files.
	 */
    if (bAppendOutputFiles)
    {
        if (fn2ftp(outputfiles[0].filename)!=efLOG)
        {
            /* make sure first file is log file so that it is OK to use it for 
             * locking
             */
            gmx_fatal(FARGS,"The first output file should always be the log "
                      "file but instead is: %s", outputfiles[0].filename);
        }
        for(i=0;i<nfiles;i++)
        {
            if (outputfiles[i].offset < 0)
            {
                gmx_fatal(FARGS,"The original run wrote a file called '%s' which "
                    "is larger than 2 GB, but mdrun did not support large file"
                    " offsets. Can not append. Run mdrun with -noappend",
                    outputfiles[i].filename);
            }
#ifdef GMX_FAHCORE
            chksum_file=gmx_fio_open(outputfiles[i].filename,"a");

#else
            chksum_file=gmx_fio_open(outputfiles[i].filename,"r+");

            /* lock log file */                
            if (i==0)
            {
#if !((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__) 
                if (fcntl(fileno(gmx_fio_getfp(chksum_file)), F_SETLK, &fl)
                    ==-1)
#else
                if (_locking(fileno(gmx_fio_getfp(chksum_file)), _LK_NBLCK, LONG_MAX)==-1)
#endif
                {
                    if (errno!=EACCES && errno!=EAGAIN)
                    {
                        gmx_fatal(FARGS,"Failed to lock: %s. %s.",
                                  outputfiles[i].filename, strerror(errno));
                    }
                    else 
                    {
                        gmx_fatal(FARGS,"Failed to lock: %s. Already running "
                                  "simulation?", outputfiles[i].filename);
                    }
                }
            }
            
            /* compute md5 chksum */ 
            if (outputfiles[i].chksum_size != -1)
            {
                if (gmx_fio_get_file_md5(chksum_file,outputfiles[i].offset,
                                     digest) != outputfiles[i].chksum_size)  /*at the end of the call the file position is at the end of the file*/
                {
                    gmx_fatal(FARGS,"Can't read %d bytes of '%s' to compute checksum. The file has been replaced or its contents has been modified.",
                              outputfiles[i].chksum_size, 
                              outputfiles[i].filename);
                }
            } 
            if (i==0)  /*log file needs to be seeked in case we need to truncate (other files are truncated below)*/
            {
                if (gmx_fio_seek(chksum_file,outputfiles[i].offset))
                {
                	gmx_fatal(FARGS,"Seek error! Failed to truncate log-file: %s.", strerror(errno));
                }
            }
#endif
            
            if (i==0) /*open log file here - so that lock is never lifted 
                        after chksum is calculated */
            {
                *pfplog = gmx_fio_getfp(chksum_file);
            }
            else
            {
                gmx_fio_close(chksum_file);
            }
#ifndef GMX_FAHCORE            
            /* compare md5 chksum */
            if (outputfiles[i].chksum_size != -1 &&
                memcmp(digest,outputfiles[i].chksum,16)!=0) 
            {
                if (debug)
                {
                    fprintf(debug,"chksum for %s: ",outputfiles[i].filename);
                    for (j=0; j<16; j++)
                    {
                        fprintf(debug,"%02x",digest[j]);
                    }
                    fprintf(debug,"\n");
                }
                gmx_fatal(FARGS,"Checksum wrong for '%s'. The file has been replaced or its contents has been modified.",
                          outputfiles[i].filename);
            }
#endif        

              
            if (i!=0) /*log file is already seeked to correct position */
            {
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
                rc = gmx_wintruncate(outputfiles[i].filename,outputfiles[i].offset);
#else            
                rc = truncate(outputfiles[i].filename,outputfiles[i].offset);
#endif
                if(rc!=0)
                {
                    gmx_fatal(FARGS,"Truncation of file %s failed.",outputfiles[i].filename);
                }
            }
        }
    }

    sfree(outputfiles);
}


void load_checkpoint(const char *fn,FILE **fplog,
                     t_commrec *cr,gmx_bool bPartDecomp,ivec dd_nc,
                     t_inputrec *ir,t_state *state,
                     gmx_bool *bReadRNG,gmx_bool *bReadEkin,gmx_bool bAppend)
{
    gmx_large_int_t step;
    double t;

    if (SIMMASTER(cr)) {
      /* Read the state from the checkpoint file */
      read_checkpoint(fn,fplog,
                      cr,bPartDecomp,dd_nc,
                      ir->eI,&step,&t,state,bReadRNG,bReadEkin,
                      &ir->simulation_part,bAppend);
    }
    if (PAR(cr)) {
      gmx_bcast(sizeof(cr->npmenodes),&cr->npmenodes,cr);
      gmx_bcast(DIM*sizeof(dd_nc[0]),dd_nc,cr);
      gmx_bcast(sizeof(step),&step,cr);
      gmx_bcast(sizeof(*bReadRNG),bReadRNG,cr);
      gmx_bcast(sizeof(*bReadEkin),bReadEkin,cr);
    }
    ir->bContinuation    = TRUE;
    if (ir->nsteps >= 0)
    {
        ir->nsteps          += ir->init_step - step;
    }
    ir->init_step        = step;
	ir->simulation_part += 1;
}

static void read_checkpoint_data(t_fileio *fp,int *simulation_part,
                                 gmx_large_int_t *step,double *t,t_state *state,
                                 gmx_bool bReadRNG,
                                 int *nfiles,gmx_file_position_t **outputfiles)
{
    int  file_version;
    char *version,*btime,*buser,*bmach,*fprog,*ftime;
    int  eIntegrator;
    int  nppnodes,npme;
    ivec dd_nc;
    int  flags_eks,flags_enh;
    int  nfiles_loc;
    gmx_file_position_t *files_loc=NULL;
    int  ret;
	
    do_cpt_header(gmx_fio_getxdr(fp),TRUE,&file_version,
                  &version,&btime,&buser,&bmach,&fprog,&ftime,
                  &eIntegrator,simulation_part,step,t,&nppnodes,dd_nc,&npme,
                  &state->natoms,&state->ngtc,&state->nnhpres,&state->nhchainlength,
                  &state->flags,&flags_eks,&flags_enh,NULL);
    ret =
        do_cpt_state(gmx_fio_getxdr(fp),TRUE,state->flags,state,bReadRNG,NULL);
    if (ret)
    {
        cp_error();
    }
    ret = do_cpt_ekinstate(gmx_fio_getxdr(fp),TRUE,
                           flags_eks,&state->ekinstate,NULL);
    if (ret)
    {
        cp_error();
    }
    ret = do_cpt_enerhist(gmx_fio_getxdr(fp),TRUE,
                          flags_enh,&state->enerhist,NULL);
    if (ret)
    {
        cp_error();
    }

    ret = do_cpt_files(gmx_fio_getxdr(fp),TRUE,
                       outputfiles != NULL ? outputfiles : &files_loc,
                       outputfiles != NULL ? nfiles : &nfiles_loc,
                       NULL,file_version);
    if (files_loc != NULL)
    {
        sfree(files_loc);
    }
	
    if (ret)
    {
        cp_error();
    }
	
    ret = do_cpt_footer(gmx_fio_getxdr(fp),TRUE,file_version);
    if (ret)
    {
        cp_error();
    }

    sfree(fprog);
    sfree(ftime);
    sfree(btime);
    sfree(buser);
    sfree(bmach);
}

void 
read_checkpoint_state(const char *fn,int *simulation_part,
                      gmx_large_int_t *step,double *t,t_state *state)
{
    t_fileio *fp;
    
    fp = gmx_fio_open(fn,"r");
    read_checkpoint_data(fp,simulation_part,step,t,state,FALSE,NULL,NULL);
    if( gmx_fio_close(fp) != 0)
	{
		gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of quota?");
	}
}

void read_checkpoint_trxframe(t_fileio *fp,t_trxframe *fr)
{
    t_state state;
    int simulation_part;
    gmx_large_int_t step;
    double t;
    
    init_state(&state,0,0,0,0);
    
    read_checkpoint_data(fp,&simulation_part,&step,&t,&state,FALSE,NULL,NULL);
    
    fr->natoms  = state.natoms;
    fr->bTitle  = FALSE;
    fr->bStep   = TRUE;
    fr->step    = gmx_large_int_to_int(step,
                                    "conversion of checkpoint to trajectory");
    fr->bTime   = TRUE;
    fr->time    = t;
    fr->bLambda = TRUE;
    fr->lambda  = state.lambda;
    fr->bAtoms  = FALSE;
    fr->bX      = (state.flags & (1<<estX));
    if (fr->bX)
    {
        fr->x     = state.x;
        state.x   = NULL;
    }
    fr->bV      = (state.flags & (1<<estV));
    if (fr->bV)
    {
        fr->v     = state.v;
        state.v   = NULL;
    }
    fr->bF      = FALSE;
    fr->bBox    = (state.flags & (1<<estBOX));
    if (fr->bBox)
    {
        copy_mat(state.box,fr->box);
    }
    done_state(&state);
}

void list_checkpoint(const char *fn,FILE *out)
{
    t_fileio *fp;
    int  file_version;
    char *version,*btime,*buser,*bmach,*fprog,*ftime;
    int  eIntegrator,simulation_part,nppnodes,npme;
    gmx_large_int_t step;
    double t;
    ivec dd_nc;
    t_state state;
    int  flags_eks,flags_enh;
    int  indent;
    int  i,j;
    int  ret;
    gmx_file_position_t *outputfiles;
	int  nfiles;
	
    init_state(&state,-1,-1,-1,-1);

    fp = gmx_fio_open(fn,"r");
    do_cpt_header(gmx_fio_getxdr(fp),TRUE,&file_version,
                  &version,&btime,&buser,&bmach,&fprog,&ftime,
                  &eIntegrator,&simulation_part,&step,&t,&nppnodes,dd_nc,&npme,
                  &state.natoms,&state.ngtc,&state.nnhpres,&state.nhchainlength,
                  &state.flags,&flags_eks,&flags_enh,out);
    ret = do_cpt_state(gmx_fio_getxdr(fp),TRUE,state.flags,&state,TRUE,out);
    if (ret)
    {
        cp_error();
    }
    ret = do_cpt_ekinstate(gmx_fio_getxdr(fp),TRUE,
                           flags_eks,&state.ekinstate,out);
    if (ret)
    {
        cp_error();
    }
    ret = do_cpt_enerhist(gmx_fio_getxdr(fp),TRUE,
                          flags_enh,&state.enerhist,out);

    if (ret == 0)
    {
		do_cpt_files(gmx_fio_getxdr(fp),TRUE,&outputfiles,&nfiles,out,file_version);
	}
	
    if (ret == 0)
    {
        ret = do_cpt_footer(gmx_fio_getxdr(fp),TRUE,file_version);
    }
	
    if (ret)
    {
        cp_warning(out);
    }
    if( gmx_fio_close(fp) != 0)
	{
		gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of quota?");
	}
    
    done_state(&state);
}


static gmx_bool exist_output_file(const char *fnm_cp,int nfile,const t_filenm fnm[])
{
    int i;

    /* Check if the output file name stored in the checkpoint file
     * is one of the output file names of mdrun.
     */
    i = 0;
    while (i < nfile &&
           !(is_output(&fnm[i]) && strcmp(fnm_cp,fnm[i].fns[0]) == 0))
    {
        i++;
    }
    
    return (i < nfile && gmx_fexist(fnm_cp));
}

/* This routine cannot print tons of data, since it is called before the log file is opened. */
gmx_bool read_checkpoint_simulation_part(const char *filename, int *simulation_part,
                                     gmx_large_int_t *cpt_step,t_commrec *cr,
                                     gmx_bool bAppendReq,
                                     int nfile,const t_filenm fnm[],
                                     const char *part_suffix,gmx_bool *bAddPart)
{
    t_fileio *fp;
    gmx_large_int_t step=0;
	double t;
    t_state state;
    int  nfiles;
    gmx_file_position_t *outputfiles;
    int  nexist,f;
    gmx_bool bAppend;
    char *fn,suf_up[STRLEN];

    bAppend = FALSE;

    if (SIMMASTER(cr)) {
        if(!gmx_fexist(filename) || (!(fp = gmx_fio_open(filename,"r")) ))
        {
            *simulation_part = 0;
        }
        else 
        {
            init_state(&state,0,0,0,0);

            read_checkpoint_data(fp,simulation_part,&step,&t,&state,FALSE,
                                 &nfiles,&outputfiles);
            if( gmx_fio_close(fp) != 0)
            {
                gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of quota?");
            }
            done_state(&state);

            if (bAppendReq)
            {
                nexist = 0;
                for(f=0; f<nfiles; f++)
                {
                    if (exist_output_file(outputfiles[f].filename,nfile,fnm))
                    {
                        nexist++;
                    }
                }
                if (nexist == nfiles)
                {
                    bAppend = bAppendReq;
                }
                else if (nexist > 0)
                {
                    fprintf(stderr,
                            "Output file appending has been requested,\n"
                            "but some output files listed in the checkpoint file %s\n"
                            "are not present or are named differently by the current program:\n",
                            filename);
                    fprintf(stderr,"output files present:");
                    for(f=0; f<nfiles; f++)
                    {
                        if (exist_output_file(outputfiles[f].filename,
                                              nfile,fnm))
                        {
                            fprintf(stderr," %s",outputfiles[f].filename);
                        }
                    }
                    fprintf(stderr,"\n");
                    fprintf(stderr,"output files not present or named differently:");
                    for(f=0; f<nfiles; f++)
                    {
                        if (!exist_output_file(outputfiles[f].filename,
                                               nfile,fnm))
                        {
                            fprintf(stderr," %s",outputfiles[f].filename);
                        }
                    }
                    fprintf(stderr,"\n");
                    
                    gmx_fatal(FARGS,"File appending requested, but only %d of the %d output files are present",nexist,nfiles);
                }
            }
            
            if (bAppend)
            {
                if (nfiles == 0)
                {
                    gmx_fatal(FARGS,"File appending requested, but no output file information is stored in the checkpoint file");
                }
                fn = outputfiles[0].filename;
                if (strlen(fn) < 4 ||
                    gmx_strcasecmp(fn+strlen(fn)-4,ftp2ext(efLOG)) == 0)
                {
                    gmx_fatal(FARGS,"File appending requested, but the log file is not the first file listed in the checkpoint file");
                }
                /* Set bAddPart to whether the suffix string '.part' is present
                 * in the log file name.
                 */
                strcpy(suf_up,part_suffix);
                upstring(suf_up);
                *bAddPart = (strstr(fn,part_suffix) != NULL ||
                             strstr(fn,suf_up) != NULL);
            }

            sfree(outputfiles);
        }
    }
    if (PAR(cr))
    {
        gmx_bcast(sizeof(*simulation_part),simulation_part,cr);

        if (*simulation_part > 0 && bAppendReq)
        {
            gmx_bcast(sizeof(bAppend),&bAppend,cr);
            gmx_bcast(sizeof(*bAddPart),bAddPart,cr);
        }
    }
    if (NULL != cpt_step)
    {
        *cpt_step = step;
    }

    return bAppend;
}
