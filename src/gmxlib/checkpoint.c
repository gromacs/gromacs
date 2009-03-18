/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * $Id$
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <time.h>
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

#ifdef GMX_FAHCORE
#include "corewrap.h"
#endif

#define CPT_MAGIC1 171817
#define CPT_MAGIC2 171819

/* cpt_version should normally only be changed
 * when the header of footer format changes.
 * The state data format itself is backward and forward compatible.
 * But old code can not read a new entry that is present in the file
 * (but can read a new format when new entries are not present).
 */
static const int cpt_version = 6;

enum { ecpdtINT, ecpdtFLOAT, ecpdtDOUBLE, ecpdtNR };

const char *ecpdt_names[ecpdtNR] = { "int", "float", "double" };

const char *est_names[estNR]=
{
    "FE-lambda",
    "box", "box-rel", "box-v", "pres_prev",
    "nosehoover-xi", "thermostat-integral",
    "x", "v", "SDx", "CGp", "LD-rng", "LD-rng-i",
    "disre_initf", "disre_rm3tav",
    "orire_initf", "orire_Dtav",
};

enum { eeksEKINH_N, eeksEKINH, eeksDEKINDL, eeksMVCOS, eeksNR };

const char *eeks_names[eeksNR]=
{
    "Ekinh_n", "Ekinh", "dEkindlambda", "mv_cos"
};

enum { eenhENERGY_N, eenhENERGY_AVER, eenhENERGY_SUM, eenhENERGY_NSUM,
       eenhENERGY_SUM_SIM, eenhENERGY_NSUM_SIM, eenhNR };

const char *eenh_names[eenhNR]=
{
    "energy_n", "energy_aver", "energy_sum", "energy_nsum",
    "energy_sum_sim", "energy_nsum_sim"
};

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

static void do_cpt_string_err(XDR *xd,bool bRead,char *desc,char **s,FILE *list)
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

static int do_cpt_int(XDR *xd,char *desc,int *i,FILE *list)
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

static void do_cpt_int_err(XDR *xd,char *desc,int *i,FILE *list)
{
    if (do_cpt_int(xd,desc,i,list) < 0)
    {
        cp_error();
    }
}

static void do_cpt_step_err(XDR *xd,const char *desc,gmx_step_t *i,FILE *list)
{
    bool_t res=0;
    char   buf[22];

    res = xdr_gmx_step_t(xd,i,"reading checkpoint file");
    if (res == 0)
    {
        cp_error();
    }
    if (list)
    {
        fprintf(list,"%s = %s\n",desc,gmx_step_str(*i,buf));
    }
}

static void do_cpt_double_err(XDR *xd,char *desc,double *f,FILE *list)
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

enum { ecprREAL, ecprRVEC, ecprMATRIX };

static int do_cpte_reals_low(XDR *xd,int cptp,int ecpt,int sflags,
                             int n,real **v,
                             FILE *list,int erealtype)
{
    bool_t res=0;
#ifndef GMX_DOUBLE
    int  dtc=ecpdtFLOAT;
#else
    int  dtc=ecpdtDOUBLE;
#endif
    real *vp,*va=NULL;
    float  *vf;
    double *vd;
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
        fprintf(stderr,"Precision mismatch for state entry %s, code precision is %s, file precision is %s\n",
                st_names(cptp,ecpt),ecpdt_names[dtc],ecpdt_names[dt]);
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
    if (dt == ecpdtFLOAT)
    {
        if (dtc == ecpdtFLOAT)
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
        if (dtc != ecpdtFLOAT)
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
        if (dtc == ecpdtDOUBLE)
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
        if (dtc != ecpdtDOUBLE)
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



static int do_cpte_reals(XDR *xd,int cptp,int ecpt,int sflags,
                         int n,real **v,FILE *list)
{
    return do_cpte_reals_low(xd,cptp,ecpt,sflags,n,v,list,ecprREAL);
}

static int do_cpte_real(XDR *xd,int cptp,int ecpt,int sflags,
                        real *r,FILE *list)
{
    return do_cpte_reals_low(xd,cptp,ecpt,sflags,1,&r,list,ecprREAL);
}

static int do_cpte_ints(XDR *xd,int cptp,int ecpt,int sflags,
                        int n,int **v,FILE *list)
{
    bool_t res=0;
    int  dtc=ecpdtINT;
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
                  st_names(cptp,ecpt),ecpdt_names[dtc],ecpdt_names[dt]);
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
    int  dtc=ecpdtDOUBLE;
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
                  st_names(cptp,ecpt),ecpdt_names[dtc],ecpdt_names[dt]);
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


static int do_cpte_rvecs(XDR *xd,int cptp,int ecpt,int sflags,
                         int n,rvec **v,FILE *list)
{
   return do_cpte_reals_low(xd,cptp,ecpt,sflags,
                            n*DIM,(real **)v,list,ecprRVEC);
}

static int do_cpte_matrix(XDR *xd,int cptp,int ecpt,int sflags,
                          matrix v,FILE *list)
{
    real *vr;
    real ret;

    vr = (real *)&(v[0][0]);
    ret = do_cpte_reals_low(xd,cptp,ecpt,sflags,DIM*DIM,&vr,NULL,ecprMATRIX);
    
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
                            nf*DIM*DIM,&vr,NULL,ecprMATRIX);
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

static void do_cpt_header(XDR *xd,bool bRead,int *file_version,
                          char **version,char **btime,char **buser,char **bmach,
                          char **fprog,char **ftime,
                          int *eIntegrator,int *simulation_part,
                          gmx_step_t *step,double *t,
                          int *nnodes,int *dd_nc,int *npme,
                          int *natoms,int *ngtc,
                          int *flags_state,int *flags_eks,int *flags_enh,
                          FILE *list)
{
    bool_t res=0;
    int  magic;
    int  idum=0;
    int  i;
    
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
    do_cpt_int_err(xd,"#atoms"            ,natoms     ,list);
    do_cpt_int_err(xd,"#T-coupling groups",ngtc       ,list);
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

static int do_cpt_footer(XDR *xd,bool bRead,int file_version)
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

static int do_cpt_state(XDR *xd,bool bRead,
                        int fflags,t_state *state,
                        bool bReadRNG,FILE *list)
{
    int  sflags;
    int  **rng_p,**rngi_p;
    int  i;
    int  ret;
    
    ret = 0;

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
            case estLAMBDA:  ret = do_cpte_real  (xd,0,i,sflags,&state->lambda,list); break;
            case estBOX:     ret = do_cpte_matrix(xd,0,i,sflags,state->box,list); break;
            case estBOX_REL: ret = do_cpte_matrix(xd,0,i,sflags,state->box_rel,list); break;
            case estBOXV:    ret = do_cpte_matrix(xd,0,i,sflags,state->boxv,list); break;
            case estPRES_PREV: ret = do_cpte_matrix(xd,0,i,sflags,state->pres_prev,list); break;
            case estNH_XI:   ret = do_cpte_reals (xd,0,i,sflags,state->ngtc,&state->nosehoover_xi,list); break;
            case estTC_INT:  ret = do_cpte_doubles(xd,0,i,sflags,state->ngtc,&state->therm_integral,list); break;
            case estX:       ret = do_cpte_rvecs (xd,0,i,sflags,state->natoms,&state->x,list); break;
            case estV:       ret = do_cpte_rvecs (xd,0,i,sflags,state->natoms,&state->v,list); break;
            case estSDX:     ret = do_cpte_rvecs (xd,0,i,sflags,state->natoms,&state->sd_X,list); break;
            case estLD_RNG:  ret = do_cpte_ints  (xd,0,i,sflags,state->nrng,rng_p,list); break;
            case estLD_RNGI: ret = do_cpte_ints (xd,0,i,sflags,state->nrngi,rngi_p,list); break;
            case estDISRE_INITF:  ret = do_cpte_real (xd,0,i,sflags,&state->hist.disre_initf,list); break;
            case estDISRE_RM3TAV: ret = do_cpte_reals(xd,0,i,sflags,state->hist.ndisrepairs,&state->hist.disre_rm3tav,list); break;
            case estORIRE_INITF:  ret = do_cpte_real (xd,0,i,sflags,&state->hist.orire_initf,list); break;
            case estORIRE_DTAV:   ret = do_cpte_reals(xd,0,i,sflags,state->hist.norire_Dtav,&state->hist.orire_Dtav,list); break;
            default:
                gmx_fatal(FARGS,"Unknown state entry %d\n"
                          "You are probably reading a new checkpoint file with old code",i);
            }
        }
    }
    
    return ret;
}

static int do_cpt_ekinstate(XDR *xd,bool bRead,
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

			case eeksEKINH_N: ret = do_cpte_int(xd,1,i,fflags,&ekins->ekinh_n,list); break;
			case eeksEKINH:   ret = do_cpte_matrices(xd,1,i,fflags,ekins->ekinh_n,&ekins->ekinh,list); break;
 			case eeksDEKINDL: ret = do_cpte_real(xd,1,i,fflags,&ekins->dekindl,list); break;
            case eeksMVCOS:   ret = do_cpte_real(xd,1,i,fflags,&ekins->mvcos,list); break;			
            default:
                gmx_fatal(FARGS,"Unknown ekin data state entry %d\n"
                          "You are probably reading a new checkpoint file with old code",i);
            }
        }
    }
    
    return ret;
}

static int do_cpt_enerhist(XDR *xd,bool bRead,
                           int fflags,energyhistory_t *enerhist,
                           FILE *list)
{
    int  i;
    int  ret;

    ret = 0;

    if (bRead)
    {
        enerhist->nsum     = 0;
        enerhist->nsum_sim = 0;
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
            case eenhENERGY_NSUM_SIM: do_cpt_step_err(xd,eenh_names[i],&enerhist->nsum_sim,list); break;
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

    return ret;
}

static int
do_cpt_files(XDR *xd, bool bRead, gmx_file_position_t **p_outputfiles, int *nfiles, FILE *list)
{
	int    i;
	off_t  offset;
	off_t  mask = 0xFFFFFFFFL;
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
#if (SIZEOF_OFF_T > SIZEOF_INT)
			outputfiles[i].offset = ( ((off_t) offset_high) << 32 ) | ( (off_t) offset_low );
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
#if (SIZEOF_OFF_T > SIZEOF_INT)
			offset_low  = (int) (offset & mask);
			offset_high = (int) ((offset >> 32) & mask);
#else
			offset_low  = offset;
			offset_high = 0;
#endif
			if (do_cpt_int(xd,"file_offset_high",&offset_high,list) != 0)
			{
				return -1;
			}
			if (do_cpt_int(xd,"file_offset_low",&offset_low,list) != 0)
			{
				return -1;
			}
		}
	}
	return 0;
}


void write_checkpoint(char *fn,FILE *fplog,t_commrec *cr,
                      int eIntegrator,int simulation_part,
                      gmx_step_t step,double t,t_state *state)
{
    int  fp;
    int  file_version;
    char *version=VERSION;
    char *btime=BUILD_TIME;
    char *buser=BUILD_USER;
    char *bmach=BUILD_MACHINE;
    char *fprog;
    char *ftime;
    time_t now;
    int  nppnodes,npmenodes,flag_64bit;
    char buf[1024];
	gmx_file_position_t *outputfiles;
	int  noutputfiles;
    int  flags_eks,flags_enh,i;
		
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
    
    if (gmx_fexist(fn))
    {
        /* Rename the previous checkpoint file */
        strcpy(buf,fn);
        buf[strlen(fn) - strlen(ftp2ext(fn2ftp(fn))) - 1] = '\0';
        strcat(buf,"_prev");
        strcat(buf,fn+strlen(fn) - strlen(ftp2ext(fn2ftp(fn))) - 1);
		(void)remove(buf); /* Unix will overwrite buf if it exists, but for windows we need to remove first */
        if(rename(fn,buf) != 0)
		{
			gmx_file("Cannot rename checkpoint file; maybe you are out of quota?");
		}
    }
    
    fprog = Program();
    now = time(NULL);
    ftime = strdup(ctime(&now));
    ftime[strlen(ftime)-1] = '\0';

	/* No need to pollute stderr every time we write a checkpoint file */
    /* fprintf(stderr,"\nWriting checkpoint, step %d at %s\n",step,ftime); */
    if (fplog)
    { 
        fprintf(fplog,"Writing checkpoint, step %s at %s\n\n",
                gmx_step_str(step,buf),ftime);
    }
    
    /* Get offsets for open files */
	gmx_fio_get_output_file_positions(&outputfiles, &noutputfiles);

    fp = gmx_fio_open(fn,"w");
	
    flags_eks =
        ((1<<eeksEKINH_N) | (1<<eeksEKINH) | (1<<eeksDEKINDL) |
         (1<<eeksMVCOS));

    flags_enh = 0;
    if (state->enerhist.nsum > 0 || state->enerhist.nsum_sim > 0)
    {
        flags_enh |= (1<<eenhENERGY_N);
        if (state->enerhist.nsum > 0)
        {
            flags_enh |= ((1<<eenhENERGY_AVER) | (1<<eenhENERGY_SUM) |
                          (1<<eenhENERGY_NSUM));
        }
        if (state->enerhist.nsum_sim > 0)
        {
            flags_enh |= ((1<<eenhENERGY_SUM_SIM) | (1<<eenhENERGY_NSUM_SIM));
        }
    }

    do_cpt_header(gmx_fio_getxdr(fp),FALSE,&file_version,
                  &version,&btime,&buser,&bmach,&fprog,&ftime,
                  &eIntegrator,&simulation_part,&step,&t,&nppnodes,
                  DOMAINDECOMP(cr) ? cr->dd->nc : NULL,&npmenodes,
                  &state->natoms,&state->ngtc,
                  &state->flags,&flags_eks,&flags_enh,NULL);

    if( (do_cpt_state(gmx_fio_getxdr(fp),FALSE,state->flags,state,TRUE,NULL) < 0)          ||
		(do_cpt_ekinstate(gmx_fio_getxdr(fp),FALSE,flags_eks,&state->ekinstate,NULL) < 0)  ||
		(do_cpt_enerhist(gmx_fio_getxdr(fp),FALSE,flags_enh,&state->enerhist,NULL) < 0)    ||
	    (do_cpt_files(gmx_fio_getxdr(fp),FALSE,&outputfiles,&noutputfiles,NULL) < 0))
	{
		gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of quota?");
	}

    do_cpt_footer(gmx_fio_getxdr(fp),FALSE,file_version);

    if( gmx_fio_close(fp) != 0)
	{
		gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of quota?");
	}
    
    sfree(ftime);
	sfree(outputfiles);
	#ifdef GMX_FAHCORE
    /*code for alternate checkpointing scheme.  moved from top of loop over steps */
      if ( fcCheckPointParallel( (MASTER(cr)==0), NULL) == 0 ) {
        gmx_fatal( 3,__FILE__,__LINE__, "Checkpoint error on step %d\n", step );
      }

	#endif /* end FAHCORE block */

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

static void check_int(FILE *fplog,char *type,int p,int f,bool *mm)
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

static void check_string(FILE *fplog,char *type,char *p,char *f,bool *mm)
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
                        t_commrec *cr,bool bPartDecomp,int npp_f,int npme_f,
                        ivec dd_nc,ivec dd_nc_f)
{
    int  npp;
    bool mm;
    
    mm = FALSE;
    
    check_string(fplog,"Version"      ,VERSION      ,version,&mm);
    check_string(fplog,"Build time"   ,BUILD_TIME   ,btime  ,&mm);
    check_string(fplog,"Build user"   ,BUILD_USER   ,buser  ,&mm);
    check_string(fplog,"Build machine",BUILD_MACHINE,bmach  ,&mm);
    check_string(fplog,"Program name" ,Program()    ,fprog  ,&mm);
    
    npp = cr->nnodes - cr->npmenodes;
    check_int   (fplog,"#nodes"       ,cr->nnodes   ,npp_f+npme_f ,&mm);
    if (bPartDecomp)
    {
        dd_nc[XX] = 1;
        dd_nc[YY] = 1;
        dd_nc[ZZ] = 1;
    }
    if (npp > 1)
    {
        check_int (fplog,"#PME-nodes"  ,cr->npmenodes,npme_f     ,&mm);
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

static void 
read_checkpoint(char *fn,FILE *fplog,
                t_commrec *cr,bool bPartDecomp,ivec dd_nc,
                int eIntegrator,gmx_step_t *step,double *t,
                t_state *state,bool *bReadRNG,bool *bReadEkin,
                int *simulation_part,bool bAppendOutputFiles)
{
    int  fp,i,j;
    int  file_version;
    char *version,*btime,*buser,*bmach,*fprog,*ftime;
	char filename[STRLEN],buf[22];
    int  nppnodes,eIntegrator_f,nppnodes_f,npmenodes_f;
    ivec dd_nc_f;
    int  natoms,ngtc,fflags,flags_eks,flags_enh;
    int  d;
    int  ret;
	gmx_file_position_t *outputfiles;
	int  nfiles;
	
    char *int_warn=
        "WARNING: The checkpoint file was generator with integrator %s,\n"
        "         while the simulation uses integrator %s\n\n";
    char *sd_note=
        "NOTE: The checkpoint file was for %d nodes doing SD or BD,\n"
        "      while the simulation uses %d SD or BD nodes,\n"
        "      continuation will be exact, except for the random state\n\n";
    
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
                  &natoms,&ngtc,
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
					  "If you _really_ know what you are doing, try without the -append option.\n");
		}
        if (fplog)
        {
            fprintf(fplog,int_warn,EI(eIntegrator_f),EI(eIntegrator));
        }
    }

    if (!PAR(cr))
    {
        nppnodes = 1;
    }
    else if (bPartDecomp)
    {
        nppnodes = cr->nnodes;
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
						  "You can try without the -append option, and get more info in the log file.\n");
			}
			
            fprintf(stderr,
                    "WARNING: The checkpoint state entries do not match the simulation,\n"
                    "         see the log file for details\n\n");
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
    *bReadEkin = (flags_eks & (1<<eeksEKINH));

    ret = do_cpt_enerhist(gmx_fio_getxdr(fp),TRUE,
                          flags_enh,&state->enerhist,NULL);
    if (ret)
    {
        cp_error();
    }

    if (file_version < 6)
    {
        fprintf(stderr,"\nWARNING: Reading checkpoint file in old format, assuming that the run that generated this file started at step 0, if this is not the case the energy averages will be incorrect.\n\n");
        if (fplog)
        {
            fprintf(fplog,"\nWARNING: Reading checkpoint file in old format, assuming that the run that generated this file started at step 0, if this is not the case the energy averages will be incorrect.\n\n");
        }
        state->enerhist.nsum     = *step;
        state->enerhist.nsum_sim = *step;
    }

	ret = do_cpt_files(gmx_fio_getxdr(fp),TRUE,&outputfiles,&nfiles,NULL);
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
	
	/* If the user wants to append to output files (danger), we use the file pointer
	 * positions of the output files stored in the checkpoint file and truncate the
	 * files such that any frames written after the checkpoint time are removed.
	 *
	 * You will get REALLY fun problems if you use the -append option by provide
	 * mdrun with other input files (in-frame truncation in the wrong places). Suit yourself!
	 */
	if(bAppendOutputFiles)
	{
		for(i=0;i<nfiles;i++)
		{
			if(0 != gmx_truncatefile(outputfiles[i].filename,outputfiles[i].offset) )
            {
                gmx_fatal(FARGS,"Truncation of file %s failed.",outputfiles[i].filename);
            }
		}
	}
	
	sfree(outputfiles);
}


void load_checkpoint(char *fn,FILE *fplog,
                     t_commrec *cr,bool bPartDecomp,ivec dd_nc,
                     t_inputrec *ir,t_state *state,
                     bool *bReadRNG,bool *bReadEkin,bool bAppend)
{
    gmx_step_t step;
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
      gmx_bcast(sizeof(&bReadRNG),bReadRNG,cr);
      gmx_bcast(sizeof(&bReadEkin),bReadEkin,cr);
    }
    ir->bContinuation    = TRUE;
    ir->nsteps          += ir->init_step - step;
    ir->init_step        = step;
	ir->simulation_part += 1;
}

static void low_read_checkpoint_state(int fp,int *simulation_part,
                                      gmx_step_t *step,double *t,t_state *state,
                                      bool bReadRNG)
{
    int  file_version;
    char *version,*btime,*buser,*bmach,*fprog,*ftime;
    int  eIntegrator;
    int  nppnodes,npme;
    ivec dd_nc;
    int  flags_eks,flags_enh;
    int  ret;
    gmx_file_position_t *outputfiles;
	int  nfiles;
	
    do_cpt_header(gmx_fio_getxdr(fp),TRUE,&file_version,
                  &version,&btime,&buser,&bmach,&fprog,&ftime,
                  &eIntegrator,simulation_part,step,t,&nppnodes,dd_nc,&npme,
                  &state->natoms,&state->ngtc,
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

	ret = do_cpt_files(gmx_fio_getxdr(fp),TRUE,&outputfiles,&nfiles,NULL);
	sfree(outputfiles);
	
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
read_checkpoint_state(char *fn,int *simulation_part,
                      gmx_step_t *step,double *t,t_state *state)
{
    int  fp;
    
    fp = gmx_fio_open(fn,"r");
    low_read_checkpoint_state(fp,simulation_part,step,t,state,TRUE);
    if( gmx_fio_close(fp) != 0)
	{
		gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of quota?");
	}
}

void read_checkpoint_trxframe(int fp,t_trxframe *fr)
{
    t_state state;
    int simulation_part;
    gmx_step_t step;
    double t;
    
    init_state(&state,0,0);
    
    low_read_checkpoint_state(fp,&simulation_part,&step,&t,&state,FALSE);
    
    fr->natoms  = state.natoms;
    fr->bTitle  = FALSE;
    fr->bStep   = TRUE;
    fr->step    = gmx_step_t_to_int(step,
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

void list_checkpoint(char *fn,FILE *out)
{
    int  fp;
    int  file_version;
    char *version,*btime,*buser,*bmach,*fprog,*ftime;
    int  eIntegrator,simulation_part,nppnodes,npme;
    gmx_step_t step;
    double t;
    ivec dd_nc;
    t_state state;
    int  flags_eks,flags_enh;
    int  indent;
    int  i,j;
    int  ret;
    gmx_file_position_t *outputfiles;
	int  nfiles;
	
    init_state(&state,-1,-1);

    fp = gmx_fio_open(fn,"r");
    do_cpt_header(gmx_fio_getxdr(fp),TRUE,&file_version,
                  &version,&btime,&buser,&bmach,&fprog,&ftime,
                  &eIntegrator,&simulation_part,&step,&t,&nppnodes,dd_nc,&npme,
                  &state.natoms,&state.ngtc,
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
		do_cpt_files(gmx_fio_getxdr(fp),TRUE,&outputfiles,&nfiles,out);
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


/* This routine cannot print tons of data, since it is called before the log file is opened. */
int
read_checkpoint_simulation_part(char *filename)
{
    int  fp;
	int  file_version;
    char *version,*btime,*buser,*bmach,*fprog,*ftime;
    int  eIntegrator_f,nppnodes_f,npmenodes_f;
    ivec dd_nc_f;
    gmx_step_t step;
	double t;
    int  natoms,ngtc,fflags,flags_eks,flags_enh,simulation_part;
		
	if(!gmx_fexist(filename) || ( (fp = gmx_fio_open(filename,"r")) < 0 ))
	{
		return 0;
	}
	
    do_cpt_header(gmx_fio_getxdr(fp),
                  TRUE,&file_version,
                  &version,&btime,&buser,&bmach,&fprog,&ftime,
                  &eIntegrator_f,&simulation_part,
                  &step,&t,&nppnodes_f,dd_nc_f,&npmenodes_f,
                  &natoms,&ngtc,
                  &fflags,&flags_eks,&flags_enh,NULL);
    if( gmx_fio_close(fp) != 0)
	{
		gmx_file("Cannot read/write checkpoint; corrupt file, or maybe you are out of quota?");
	}
	
	sfree(version);
	sfree(btime);
	sfree(buser);
	sfree(bmach);
	sfree(fprog);
	sfree(ftime);
	
	return simulation_part;
}

