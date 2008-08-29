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
#include "checkpoint.h"

#define CPT_MAGIC1 171817
#define CPT_MAGIC2 171819

/* cpt_version should normally only be changed
 * when the header of footer format changes.
 * The state data format itself is backward and forward compatible.
 * But old code can not read a new entry that is present in the file
 * (but can read a new format when new entries are not present).
 */
static const int cpt_version = 2;

enum { ecpdtINT, ecpdtFLOAT, ecpdtDOUBLE, ecpdtNR };

const char *ecpdt_names[ecpdtNR] = { "int", "float", "double" };

const char *est_names[estNR]=
{
    "FE-lambda",
    "box", "box-rel", "box-v", "pres_prev",
    "nosehoover-xi", "thermostat-integral",
    "x", "v", "SDx", "CGp", "LD-rng", "LD-rng-i",
    "disre_initf", "disre_rm3tav",
    "orire_initf", "orire_Dtav"
};

static void cp_warning(FILE *fp)
{
    fprintf(fp,"\nWARNING: Checkpoint file is corrupted or truncated\n\n");
}

static void cp_error()
{
    gmx_fatal(FARGS,"Checkpoint file is corrupted or truncated");
}

static void do_cpt_string_err(XDR *xd,bool bRead,char *desc,char **s,FILE *list)
{
#define CPTSTRLEN 1024
    bool_t res=0;
    int dt=ecpdtINT;
    
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
    if (do_cpt_int(xd,desc,i,list))
    {
        cp_error();
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

static int do_cpte_reals_low(XDR *xd,int ecpt,int sflags,int n,real **v,
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
        gmx_fatal(FARGS,"Count mismatch for state entry %s, code count is %d, file count is %d\n",est_names[ecpt],n,nf);
    }
    dt = dtc;
    res = xdr_int(xd,&dt);
    if (res == 0)
    {
        return -1;
    }
    if (dt != dtc)
    {
        fprintf(stderr,"Precision mismatch for state entry %s, code precision is %s, file precision is %s\n",est_names[ecpt],
                dtc==ecpdtFLOAT ? "float" : "double",
                dt ==ecpdtFLOAT ? "float" : "double");
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
            pr_reals(list,0,est_names[ecpt],vp,nf);
            break;
        case ecprRVEC:
            pr_rvecs(list,0,est_names[ecpt],(rvec *)vp,nf/3);
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

static int do_cpte_reals(XDR *xd,int ecpt,int sflags,int n,real **v,
                         FILE *list)
{
    return do_cpte_reals_low(xd,ecpt,sflags,n,v,list,ecprREAL);
}

static int do_cpte_real(XDR *xd,int ecpt,int sflags,real *r,FILE *list)
{
    return do_cpte_reals_low(xd,ecpt,sflags,1,&r,list,ecprREAL);
}

static int do_cpte_ints(XDR *xd,int ecpt,int sflags,int n,int **v,FILE *list)
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
        gmx_fatal(FARGS,"Count mismatch for state entry %s, code count is %d, file count is %d\n",est_names[ecpt],n,nf);
    }
    dt = dtc;
    res = xdr_int(xd,&dt);
    if (res == 0)
    {
        return -1;
    }
    if (dt != dtc)
    {
        gmx_fatal(FARGS,"Type mismatch for state entry %s, code type is %s, file type is %s\n",est_names[ecpt],ecpdt_names[dtc],ecpdt_names[dt]);
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
        pr_ivec(list,0,est_names[ecpt],vp,nf,TRUE);
    }
    if (va)
    {
        sfree(va);
    }

    return 0;
}

static int do_cpte_int(XDR *xd,int ecpt,int sflags,int *i,FILE *list)
{
    return do_cpte_ints(xd,ecpt,sflags,1,&i,list);
}

static int do_cpte_doubles(XDR *xd,int ecpt,int sflags,int n,double **v,
                           FILE *list)
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
        gmx_fatal(FARGS,"Count mismatch for state entry %s, code count is %d, file count is %d\n",est_names[ecpt],n,nf);
    }
    dt = dtc;
    res = xdr_int(xd,&dt);
    if (res == 0)
    {
        return -1;
    }
    if (dt != dtc)
    {
        gmx_fatal(FARGS,"Precision mismatch for state entry %s, code precision is %s, file precision is %s\n",est_names[ecpt],
                  dtc==ecpdtFLOAT ? "float" : "double",
                  dt ==ecpdtFLOAT ? "float" : "double");
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
        pr_doubles(list,0,est_names[ecpt],vp,nf);
    }
    if (va)
    {
        sfree(va);
    }

    return 0;
}

static int do_cpte_rvecs(XDR *xd,int ecpt,int sflags,int n,rvec **v,
                         FILE *list)
{
   return do_cpte_reals_low(xd,ecpt,sflags,n*DIM,(real **)v,list,ecprRVEC);
}

static int do_cpte_matrix(XDR *xd,int ecpt,int sflags,matrix v,FILE *list)
{
    real *vr;
    real ret;

    vr = (real *)&(v[0][0]);
    ret = do_cpte_reals_low(xd,ecpt,sflags,DIM*DIM,&vr,NULL,ecprMATRIX);
    
    if (list && ret == 0)
    {
        pr_rvecs(list,0,est_names[ecpt],v,DIM);
    }
    
    return ret;
}

static void do_cpt_header(XDR *xd,bool bRead,int *file_version,
                          char **version,char **btime,char **buser,char **bmach,
                          char **fprog,char **ftime,
                          int *eIntegrator,int *step,double *t,
                          int *nnodes,int *dd_nc,int *npme,
                          int *natoms,int *ngtc,int *flags,
                          FILE *list)
{
    bool_t res=0;
    int  magic;
    int  idum;
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
        gmx_fatal(FARGS,"The checkpoint file is empty, corrupted or not a checkpoint file");
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
    do_cpt_int_err(xd,"step"              ,step       ,list);
    do_cpt_double_err(xd,"t"              ,t          ,list);
    do_cpt_int_err(xd,"#PP-nodes"         ,nnodes     ,list);
    idum = 1;
    do_cpt_int_err(xd,"dd_nc[x]",dd_nc ? &(dd_nc[0]) : &idum,list);
    do_cpt_int_err(xd,"dd_nc[y]",dd_nc ? &(dd_nc[1]) : &idum,list);
    do_cpt_int_err(xd,"dd_nc[z]",dd_nc ? &(dd_nc[2]) : &idum,list);
    do_cpt_int_err(xd,"#PME-only nodes",npme,list);
    do_cpt_int_err(xd,"state flags",flags,list);
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
            case estLAMBDA:  ret = do_cpte_real  (xd,i,sflags,&state->lambda,list); break;
            case estBOX:     ret = do_cpte_matrix(xd,i,sflags,state->box,list); break;
            case estBOX_REL: ret = do_cpte_matrix(xd,i,sflags,state->box_rel,list); break;
            case estBOXV:    ret = do_cpte_matrix(xd,i,sflags,state->boxv,list); break;
            case estPRES_PREV: ret = do_cpte_matrix(xd,i,sflags,state->pres_prev,list); break;
            case estNH_XI:   ret = do_cpte_reals (xd,i,sflags,state->ngtc,&state->nosehoover_xi,list); break;
            case estTC_INT:  ret = do_cpte_doubles(xd,i,sflags,state->ngtc,&state->therm_integral,list); break;
            case estX:       ret = do_cpte_rvecs (xd,i,sflags,state->natoms,&state->x,list); break;
            case estV:       ret = do_cpte_rvecs (xd,i,sflags,state->natoms,&state->v,list); break;
            case estSDX:     ret = do_cpte_rvecs (xd,i,sflags,state->natoms,&state->sd_X,list); break;
            case estLD_RNG:  ret = do_cpte_ints  (xd,i,sflags,state->nrng,rng_p,list); break;
            case estLD_RNGI: ret = do_cpte_ints (xd,i,sflags,state->nrngi,rngi_p,list); break;
            case estDISRE_INITF:  ret = do_cpte_real (xd,i,sflags,&state->hist.disre_initf,list); break;
            case estDISRE_RM3TAV: ret = do_cpte_reals(xd,i,sflags,state->hist.ndisrepairs,&state->hist.disre_rm3tav,list); break;
            case estORIRE_INITF:  ret = do_cpte_real (xd,i,sflags,&state->hist.orire_initf,list); break;
            case estORIRE_DTAV:   ret = do_cpte_reals(xd,i,sflags,state->hist.norire_Dtav,&state->hist.orire_Dtav,list); break;
            default:
                gmx_fatal(FARGS,"Unknown state entry %d\n"
                          "You are probably reading a new checkpoint file with old code",i);
            }
        }
    }
    
    return ret;
}

void write_checkpoint(char *fn,FILE *fplog,t_commrec *cr,
                      int eIntegrator,int step,double t,t_state *state)
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
    int  nppnodes,npmenodes;
    char buf[1024];
    
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
    
    if (fexist(fn))
    {
        /* Rename the previous checkpoint file */
        strcpy(buf,fn);
        buf[strlen(fn) - strlen(ftp2ext(fn2ftp(fn))) - 1] = '\0';
        strcat(buf,"_prev");
        strcat(buf,fn+strlen(fn) - strlen(ftp2ext(fn2ftp(fn))) - 1);
        rename(fn,buf);
    }
    
    fprog = Program();
    now = time(NULL);
    ftime = strdup(ctime(&now));
    ftime[strlen(ftime)-1] = '\0';
    
    fprintf(stderr,"\nWriting checkpoint, step %d at %s\n",step,ftime);
    if (fplog)
    {
        fprintf(fplog,"Writing checkpoint, step %d at %s\n\n",step,ftime);
    }
    
    fp = gmx_fio_open(fn,"w");
    do_cpt_header(gmx_fio_getxdr(fp),FALSE,&file_version,
                  &version,&btime,&buser,&bmach,&fprog,&ftime,
                  &eIntegrator,&step,&t,&nppnodes,
                  DOMAINDECOMP(cr) ? cr->dd->nc : NULL,&npmenodes,
                  &state->natoms,&state->ngtc,&state->flags,NULL);
    do_cpt_state(gmx_fio_getxdr(fp),FALSE,state->flags,state,TRUE,NULL);
    do_cpt_footer(gmx_fio_getxdr(fp),FALSE,file_version);
    

    gmx_fio_close(fp);
    
    sfree(ftime);
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
    if (p != f)
    {
        if (fplog)
        {
            fprintf(fplog,"  %s mismatch,\n",type);
            fprintf(fplog,"    current program: %d\n",p);
            fprintf(fplog,"    checkpoint file: %d\n",f);
            fprintf(fplog,"\n");
        }
        *mm = TRUE;
    }
}

static void check_string(FILE *fplog,char *type,char *p,char *f,bool *mm)
{
    if (strcmp(p,f) != 0)
    {
        if (fplog)
        {
            fprintf(fplog,"  %s mismatch,\n",type);
            fprintf(fplog,"    current program: %s\n",p);
            fprintf(fplog,"    checkpoint file: %s\n",f);
            fprintf(fplog,"\n");
        }
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
    
    check_string(fplog,"Version"     ,VERSION      ,version,&mm);
    check_string(fplog,"Build time"  ,BUILD_TIME   ,btime  ,&mm);
    check_string(fplog,"Build user"  ,BUILD_USER   ,buser  ,&mm);
    check_string(fplog,"Build time"  ,BUILD_MACHINE,bmach  ,&mm);
    check_string(fplog,"Program name",Program()    ,fprog  ,&mm);
    
    npp = cr->nnodes - cr->npmenodes;
    check_int   (fplog,"#PP-nodes"   ,npp          ,npp_f      ,&mm);
    if (bPartDecomp)
    {
        dd_nc[XX] = 1;
        dd_nc[YY] = 1;
        dd_nc[ZZ] = 1;
    }
    if (npp == npp_f && npp > 1)
    {
        check_int (fplog,"#PME-nodes"  ,cr->npmenodes,npme_f     ,&mm);
        check_int (fplog,"#DD-cells[x]",dd_nc[XX]    ,dd_nc_f[XX],&mm);
        check_int (fplog,"#DD-cells[y]",dd_nc[YY]    ,dd_nc_f[YY],&mm);
        check_int (fplog,"#DD-cells[z]",dd_nc[ZZ]    ,dd_nc_f[ZZ],&mm);
    }
    
    if (mm)
    {
        if (fplog)
        {
            fprintf(fplog,"Continuation is exact, but is not guaranteed to be binary identical\n\n");
            fprintf(stderr,"Continuation is exact, but is not guaranteed to be binary identical,\n"
                    "see the log file for details\n\n");
        }
    }
}

void read_checkpoint(char *fn,FILE *fplog,
                     t_commrec *cr,bool bPartDecomp,ivec dd_nc,
                     int eIntegrator,int *step,double *t,
                     t_state *state,bool *bReadRNG)
{
    int  fp;
    int  file_version;
    char *version,*btime,*buser,*bmach,*fprog,*ftime;
    int  nppnodes,eIntegrator_f,nppnodes_f,npmenodes_f;
    ivec dd_nc_f;
    int  natoms,ngtc,fflags;
    int  d;
    int  ret;
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
                  "read_checkpoint not supported with particle decomposition");
    }
    
    fp = gmx_fio_open(fn,"r");
    do_cpt_header(gmx_fio_getxdr(fp),TRUE,&file_version,
                  &version,&btime,&buser,&bmach,&fprog,&ftime,
                  &eIntegrator_f,step,t,&nppnodes_f,dd_nc_f,&npmenodes_f,
                  &natoms,&ngtc,&fflags,NULL);
    
    if (cr == NULL || MASTER(cr))
    {
        fprintf(stderr,"\nReading checkpoint file %s generated: %s\n\n",
                fn,ftime);
    }
    if (fplog)
    {
        fprintf(fplog,"\n");
        fprintf(fplog,"Reading checkpoint file %s\n",fn);
        fprintf(fplog,"  file generated by:     %s\n",fprog);  
        fprintf(fplog,"  file generated at:     %s\n",ftime);  
        fprintf(fplog,"  GROMACS build time:    %s\n",btime);  
        fprintf(fplog,"  GROMACS build user:    %s\n",buser);  
        fprintf(fplog,"  GROMACS build machine: %s\n",bmach);  
        fprintf(fplog,"  step %d\n",*step);  
        fprintf(fplog,"  time %f\n",*t);  
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
        if (fplog)
        {
            fprintf(fplog,int_warn,EI(eIntegrator_f),EI(eIntegrator));
        }
    }

    if (!PAR(cr))
    {
        nppnodes = 1;
    } else if (bPartDecomp)
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
    
    *bReadRNG = TRUE;
    if (fflags != state->flags)
    {
        if (MASTER(cr))
        {
            fprintf(stderr,
                    "WARNING: The checkpoint state entries do not match the simulation,\n"
                    "         see the log file for details\n\n");
        }
        print_flag_mismatch(fplog,state->flags,fflags);
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
    ret = do_cpt_footer(gmx_fio_getxdr(fp),TRUE,file_version);
    if (ret)
    {
        cp_error();
    }
    gmx_fio_close(fp);
    
    sfree(fprog);
    sfree(ftime);
    sfree(btime);
    sfree(buser);
    sfree(bmach);
}

static void low_read_checkpoint_state(int fp,
                                      int *step,double *t,t_state *state,
                                      bool bReadRNG)
{
    int  file_version;
    char *version,*btime,*buser,*bmach,*fprog,*ftime;
    int  eIntegrator;
    int  nppnodes,npme;
    ivec dd_nc;
    int  ret;
    
    do_cpt_header(gmx_fio_getxdr(fp),TRUE,&file_version,
                  &version,&btime,&buser,&bmach,&fprog,&ftime,
                  &eIntegrator,step,t,&nppnodes,dd_nc,&npme,
                  &state->natoms,&state->ngtc,&state->flags,NULL);
    ret =
        do_cpt_state(gmx_fio_getxdr(fp),TRUE,state->flags,state,bReadRNG,NULL);
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

void read_checkpoint_state(char *fn,int *step,double *t,t_state *state)
{
    int  fp;
    
    fp = gmx_fio_open(fn,"r");
    low_read_checkpoint_state(fp,step,t,state,TRUE);
    gmx_fio_close(fp);
}

void read_checkpoint_trxframe(int fp,t_trxframe *fr)
{
    t_state state;
    int step;
    double t;
    
    init_state(&state,0,0);
    
    low_read_checkpoint_state(fp,&step,&t,&state,FALSE);
    
    fr->natoms  = state.natoms;
    fr->bTitle  = FALSE;
    fr->bStep   = TRUE;
    fr->step    = step;
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
    int  eIntegrator,step,nppnodes,npme;
    double t;
    ivec dd_nc;
    t_state state;
    int  indent;
    int  i;
    int  ret;
    
    init_state(&state,-1,-1);

    fp = gmx_fio_open(fn,"r");
    do_cpt_header(gmx_fio_getxdr(fp),TRUE,&file_version,
                  &version,&btime,&buser,&bmach,&fprog,&ftime,
                  &eIntegrator,&step,&t,&nppnodes,dd_nc,&npme,
                  &state.natoms,&state.ngtc,&state.flags,out);
    ret = do_cpt_state(gmx_fio_getxdr(fp),TRUE,state.flags,&state,TRUE,out);
    if (ret == 0)
    {
        ret = do_cpt_footer(gmx_fio_getxdr(fp),TRUE,file_version);
    }
    if (ret)
    {
        cp_warning(out);
    }
    gmx_fio_close(fp);
    
    done_state(&state);
}
