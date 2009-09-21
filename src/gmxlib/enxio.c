/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "futil.h"
#include "string2.h"
#include "gmx_fatal.h"
#include "smalloc.h"
#include "gmxfio.h"
#include "enxio.h"
#include "vec.h"

/* The source code in this file should be thread-safe. 
         Please keep it that way. */

/* This number should be increased whenever the file format changes! */
static const int enx_version = 3;


/* Stuff for reading pre 4.1 energy files */
typedef struct {
    bool     bOldFileOpen;   /* Is this an open old file? */
    bool     bReadFirstStep; /* Did we read the first step? */
    int      first_step;     /* First step in the energy file */
    int      step_prev;      /* Previous step */
    int      nsum_prev;      /* Previous step sum length */
    t_energy *ener_prev;     /* Previous energy sums */
} ener_old_t;

struct ener_file
{
    ener_old_t eo;
    int fp;
    int framenr;
    real frametime;
};


void free_enxframe(t_enxframe *fr)
{
  int b;

  if (fr->e_alloc)
    sfree(fr->ener);
  if (fr->d_alloc) {
    sfree(fr->disre_rm3tav);
    sfree(fr->disre_rt);
  }
  for(b=0; b<fr->nblock; b++)
    sfree(fr->block[b]);
  sfree(fr->block);
  sfree(fr->b_alloc);
  sfree(fr->nr);
}

static void edr_strings(XDR *xdr,bool bRead,int file_version,
                        int n,gmx_enxnm_t **nms)
{
    int  i;
    gmx_enxnm_t *nm;

    if (*nms == NULL)
    {
        snew(*nms,n);
    }
    for(i=0; i<n; i++)
    {
        nm = &(*nms)[i];
        if (bRead)
        {
            if (nm->name)
            {
                sfree(nm->name);
                nm->name = NULL;
            }
            if (nm->unit)
            {
                sfree(nm->unit);
                nm->unit = NULL;
            }
        }
        if(!xdr_string(xdr,&(nm->name),STRLEN))
        {
            gmx_file("Cannot write energy names to file; maybe you are out of quota?");
        }
        if (file_version >= 2)
        {
            if(!xdr_string(xdr,&(nm->unit),STRLEN))
            {
                gmx_file("Cannot write energy names to file; maybe you are out of quota?");
            }
        }
        else
        {
            nm->unit = strdup("kJ/mol");
        }
    }
}

static void gen_units(int n,char ***units)
{
    int i;

    snew(*units,n);
    for(i=0; i<n; i++)
    {
        (*units)[i] = strdup("kJ/mol");
    }
}

void do_enxnms(ener_file_t ef,int *nre,gmx_enxnm_t **nms)
{
    int  magic=-55555;
    XDR  *xdr;
    bool bRead = gmx_fio_getread(ef->fp);
    int  file_version;
    int  i;
    
    gmx_fio_select(ef->fp);

    xdr = gmx_fio_getxdr(ef->fp);
    
    if (!xdr_int(xdr,&magic))
    {
        if(!bRead)
        {
            gmx_file("Cannot write energy names to file; maybe you are out of quota?");
        }
        *nre=0;
        return;
    }
    if (magic > 0)
    {
        /* Assume this is an old edr format */
        file_version = 1;
        *nre = magic;
        ef->eo.bOldFileOpen = TRUE;
        ef->eo.bReadFirstStep = FALSE;
        srenew(ef->eo.ener_prev,*nre);
    }
    else
    {
        ef->eo.bOldFileOpen=FALSE;

        if (magic != -55555)
        {
            gmx_fatal(FARGS,"Energy names magic number mismatch, this is not a GROMACS edr file");
        }
        file_version = enx_version;
        xdr_int(xdr,&file_version);
        if (file_version > enx_version)
        {
            gmx_fatal(FARGS,"reading tpx file (%s) version %d with version %d program",gmx_fio_getname(ef->fp),file_version,enx_version);
        }
        xdr_int(xdr,nre);
    }
    if (file_version != enx_version)
    {
        fprintf(stderr,"Note: enx file_version %d, software version %d\n",
                file_version,enx_version);
    }

    edr_strings(xdr,bRead,file_version,*nre,nms);
}

static bool do_eheader(ener_file_t ef,int *file_version,t_enxframe *fr,
                       int nre_test,bool *bWrongPrecision,bool *bOK)
{
    int  magic=-7777777;
    real r;
    int  block,i,zero=0,dum=0;
    bool bRead = gmx_fio_getread(ef->fp);
    int  tempfix_nr=0;
    
    if (nre_test >= 0)
    {
        *bWrongPrecision = FALSE;
    }

    *bOK=TRUE;
    /* The original energy frame started with a real,
     * so we have to use a real for compatibility.
     * This is VERY DIRTY code, since do_eheader can be called
     * with the wrong precision set and then we could read r > -1e10,
     * while actually the intention was r < -1e10.
     * When nre_test >= 0, do_eheader should therefore terminate
     * before the number of i/o calls starts depending on what has been read
     * (which is the case for for instance the block sizes for variable
     * number of blocks, where this number is read before).
     */
    r = -2e10;
    if (!do_real(r))
    {
        return FALSE;
    }
    if (r > -1e10)
    {
        /* Assume we are reading an old format */
        *file_version = 1;
        fr->t = r;
        if (!do_int(dum))   *bOK = FALSE;
        fr->step = dum;
    }
    else
    {
        if (!do_int (magic))       *bOK = FALSE;
        if (magic != -7777777)
        {
            gmx_fatal(FARGS,"Energy header magic number mismatch, this is not a GROMACS edr file");
        }
        *file_version = enx_version;
        if (!do_int (*file_version)) *bOK = FALSE;
        if (*bOK && *file_version > enx_version)
        {
            gmx_fatal(FARGS,"reading tpx file (%s) version %d with version %d program",gmx_fio_getname(ef->fp),file_version,enx_version);
        }
        if (!do_double(fr->t))       *bOK = FALSE;
        if (!do_gmx_large_int(fr->step)) *bOK = FALSE;
        if (!bRead && fr->nsum == 1) {
            /* Do not store sums of length 1,
             * since this does not add information.
             */
            if (!do_int (zero))      *bOK = FALSE;
        } else {
            if (!do_int (fr->nsum))  *bOK = FALSE;
        }
        if (*file_version >= 3)
        {
            do_gmx_large_int(fr->nsteps);
        }
        else
        {
            fr->nsteps = max(1,fr->nsum);
        }
    }
    if (!do_int (fr->nre))     *bOK = FALSE;
    if (!do_int (fr->ndisre))  *bOK = FALSE;
    if (!do_int (fr->nblock))  *bOK = FALSE;

    /* Frames could have nre=0, so we can not rely only on the fr->nre check */
    if (bRead && nre_test >= 0 &&
        ((fr->nre > 0 && fr->nre != nre_test) ||
         fr->nre < 0 || fr->ndisre < 0 || fr->nblock < 0))
    {
        *bWrongPrecision = TRUE;
        return *bOK;
    }
	
    if (*bOK && bRead && fr->nblock>fr->nr_alloc)
    {
        srenew(fr->nr,fr->nblock);
        srenew(fr->b_alloc,fr->nblock);
        srenew(fr->block,fr->nblock);
        for(i=fr->nr_alloc; i<fr->nblock; i++) {
            fr->block[i]   = NULL;
            fr->b_alloc[i] = 0;
        }
        fr->nr_alloc = fr->nblock;
    }
    for(block=0; block<fr->nblock; block++)
    {
        if (!do_int (fr->nr[block])) 
        {
            *bOK = FALSE;
        }
    }
    if (!do_int (fr->e_size))  *bOK = FALSE;
    if (!do_int (fr->d_size))  *bOK = FALSE;
    /* Do a dummy int to keep the format compatible with the old code */
    if (!do_int (dum))         *bOK = FALSE;
    
    if (*bOK && *file_version == 1 && nre_test < 0)
    {
#if 0
        if (fp >= ener_old_nalloc)
        {
            gmx_incons("Problem with reading old format energy files");
        }
#endif
        
        if (!ef->eo.bReadFirstStep)
        {
            ef->eo.bReadFirstStep = TRUE;
            ef->eo.first_step     = fr->step;
            ef->eo.nsum_prev      = 0;
        }
        
        fr->nsum = fr->step - ef->eo.first_step + 1;
        fr->nsteps = fr->nsum;
    }
	
    return *bOK;
}

void free_enxnms(int n,gmx_enxnm_t *nms)
{
    int i;

    for(i=0; i<n; i++)
    {
        sfree(nms[i].name);
        sfree(nms[i].unit);
    }

    sfree(nms);
}

void close_enx(ener_file_t ef)
{
    if(gmx_fio_close(ef->fp) != 0)
    {
        gmx_file("Cannot close energy file; it might be corrupt, or maybe you are out of quota?");  
    }
}

static bool empty_file(const char *fn)
{
    FILE *fp;
    char dum;
    int  ret;
    bool bEmpty;
    
    fp = gmx_fio_fopen(fn,"r");
    ret = fread(&dum,sizeof(dum),1,fp);
    bEmpty = feof(fp);
    gmx_fio_fclose(fp);
    
    return bEmpty;
}


ener_file_t open_enx(const char *fn,const char *mode)
{
    int        nre,i;
    gmx_enxnm_t *nms=NULL;
    int        file_version=-1;
    t_enxframe *fr;
    bool       bWrongPrecision,bDum=TRUE;
    struct ener_file *ef;

    snew(ef,1);

    if (mode[0]=='r') {
        ef->fp=gmx_fio_open(fn,mode);
        gmx_fio_select(ef->fp);
        gmx_fio_setprecision(ef->fp,FALSE);
        do_enxnms(ef,&nre,&nms);
        snew(fr,1);
        do_eheader(ef,&file_version,fr,nre,&bWrongPrecision,&bDum);
        if(!bDum)
        {
            gmx_file("Cannot read energy file header. Corrupt file?");
        }

        /* Now check whether this file is in single precision */
        if (!bWrongPrecision &&
            ((fr->e_size && (fr->nre == nre) && 
              (nre*4*(long int)sizeof(float) == fr->e_size)) ||
             (fr->d_size && 
              (fr->ndisre*(long int)sizeof(float)*2+(long int)sizeof(int) 
               == fr->d_size))))
        {
            fprintf(stderr,"Opened %s as single precision energy file\n",fn);
            free_enxnms(nre,nms);
        }
        else
        {
            gmx_fio_rewind(ef->fp);
            gmx_fio_select(ef->fp);
            gmx_fio_setprecision(ef->fp,TRUE);
            do_enxnms(ef,&nre,&nms);
            do_eheader(ef,&file_version,fr,nre,&bWrongPrecision,&bDum);
            if(!bDum)
            {
                gmx_file("Cannot write energy file header; maybe you are out of quota?");
            }

            if (((fr->e_size && (fr->nre == nre) && 
                            (nre*4*(long int)sizeof(double) == fr->e_size)) ||
                        (fr->d_size && 
                         (fr->ndisre*(long int)sizeof(double)*2+
                                     (long int)sizeof(int) == 
                          fr->d_size))))
                fprintf(stderr,"Opened %s as double precision energy file\n",
                        fn);
            else {
                if (empty_file(fn))
                    gmx_fatal(FARGS,"File %s is empty",fn);
                else
                    gmx_fatal(FARGS,"Energy file %s not recognized, maybe different CPU?",
                              fn);
            }
            free_enxnms(nre,nms);
        }
        free_enxframe(fr);
        sfree(fr);
        gmx_fio_rewind(ef->fp);
    }
    else 
        ef->fp = gmx_fio_open(fn,mode);

    ef->framenr=0;
    ef->frametime=0;

    return ef;
}

int enx_file_pointer(const ener_file_t ef)
{
    return ef->fp;
}

static void convert_full_sums(ener_old_t *ener_old,t_enxframe *fr)
{
    int nstep_all;
    int ne,ns,i;
    double esum_all,eav_all;
    
    if (fr->nsum > 0)
    {
        ne = 0;
        ns = 0;
        for(i=0; i<fr->nre; i++)
        {
            if (fr->ener[i].e    != 0) ne++;
            if (fr->ener[i].esum != 0) ns++;
        }
        if (ne > 0 && ns == 0)
        {
            /* We do not have all energy sums */
            fr->nsum = 0;
        }
    }
    
    /* Convert old full simulation sums to sums between energy frames */
    nstep_all = fr->step - ener_old->first_step + 1;
    if (fr->nsum > 1 && fr->nsum == nstep_all && ener_old->nsum_prev > 0)
    {
        /* Set the new sum length: the frame step difference */
        fr->nsum = fr->step - ener_old->step_prev;
        for(i=0; i<fr->nre; i++)
        {
            esum_all = fr->ener[i].esum;
            eav_all  = fr->ener[i].eav;
            fr->ener[i].esum = esum_all - ener_old->ener_prev[i].esum;
            fr->ener[i].eav  = eav_all  - ener_old->ener_prev[i].eav
                - dsqr(ener_old->ener_prev[i].esum/(nstep_all - fr->nsum)
                       - esum_all/nstep_all)*
                (nstep_all - fr->nsum)*nstep_all/(double)fr->nsum;
            ener_old->ener_prev[i].esum = esum_all;
            ener_old->ener_prev[i].eav  = eav_all;
        }
        ener_old->nsum_prev = nstep_all;
    }
    else if (fr->nsum > 0)
    {
        if (fr->nsum != nstep_all)
        {
            fprintf(stderr,"\nWARNING: something is wrong with the energy sums, will not use exact averages\n");
            ener_old->nsum_prev = 0;
        }
        else
        {
            ener_old->nsum_prev = nstep_all;
        }
        /* Copy all sums to ener_prev */
        for(i=0; i<fr->nre; i++)
        {
            ener_old->ener_prev[i].esum = fr->ener[i].esum;
            ener_old->ener_prev[i].eav  = fr->ener[i].eav;
        }
    }
    
    ener_old->step_prev = fr->step;
}

bool do_enx(ener_file_t ef,t_enxframe *fr)
{
    int       file_version=-1;
    int       i,block;
    bool      bRead,bOK,bOK1,bSane;
    real      tmp1,tmp2,rdum;
    char      buf[22];
    
    bOK = TRUE;
    bRead = gmx_fio_getread(ef->fp);
    if (!bRead)
    {  
        fr->e_size = fr->nre*sizeof(fr->ener[0].e)*4;
        fr->d_size = fr->ndisre*(sizeof(fr->disre_rm3tav[0]) + 
                                 sizeof(fr->disre_rt[0]));
    }
    gmx_fio_select(ef->fp);
    
    if (!do_eheader(ef,&file_version,fr,-1,NULL,&bOK))
    {
        if (bRead)
        {
            fprintf(stderr,"\rLast energy frame read %d time %8.3f           ",
                    ef->framenr-1,ef->frametime);
            if (!bOK)
            {
                fprintf(stderr,
                        "\nWARNING: Incomplete energy frame: nr %d time %8.3f\n",
                        ef->framenr,fr->t);
            }
        }
        else
        {
            gmx_file("Cannot write energy file header; maybe you are out of quota?");
        }
        return FALSE;
    }
    if (bRead)
    {
        if ((ef->framenr <   20 || ef->framenr %   10 == 0) &&
            (ef->framenr <  200 || ef->framenr %  100 == 0) &&
            (ef->framenr < 2000 || ef->framenr % 1000 == 0))
        {
            fprintf(stderr,"\rReading energy frame %6d time %8.3f           ",
                    ef->framenr,fr->t);
        }
        ef->framenr++;
        ef->frametime = fr->t;
    }
    /* Check sanity of this header */
    bSane = (fr->nre > 0 || fr->ndisre > 0);
    for(block=0; block<fr->nblock; block++)
    {
        bSane = bSane || (fr->nr[block] > 0);
    }
    if (!((fr->step >= 0) && bSane))
    {
        fprintf(stderr,"\nWARNING: there may be something wrong with energy file %s\n",
                gmx_fio_getname(ef->fp));
        fprintf(stderr,"Found: step=%s, nre=%d, ndisre=%d, nblock=%d, time=%g.\n"
                "Trying to skip frame expect a crash though\n",
                gmx_step_str(fr->step,buf),fr->nre,fr->ndisre,fr->nblock,fr->t);
    }
    if (bRead && fr->nre > fr->e_alloc)
    {
        srenew(fr->ener,fr->nre);
        for(i=fr->e_alloc; (i<fr->nre); i++)
        {
            fr->ener[i].e    = 0;
            fr->ener[i].eav  = 0;
            fr->ener[i].esum = 0;
        }
        fr->e_alloc = fr->nre;
    }
    
    for(i=0; i<fr->nre; i++)
    {
        bOK = bOK && do_real(fr->ener[i].e);
        
        /* Do not store sums of length 1,
         * since this does not add information.
         */
        if (file_version == 1 ||
            (bRead && fr->nsum > 0) || fr->nsum > 1)
        {
            tmp1 = fr->ener[i].eav;
            bOK = bOK && do_real(tmp1);
            if (bRead)
                fr->ener[i].eav = tmp1;
            
            /* This is to save only in single precision (unless compiled in DP) */
            tmp2 = fr->ener[i].esum;
            bOK = bOK && do_real(tmp2);
            if (bRead)
                fr->ener[i].esum = tmp2;
            
            if (file_version == 1)
            {
                /* Old, unused real */
                rdum = 0;
                bOK = bOK && do_real(rdum);
            }
        }
    }
    
    /* Here we can not check for file_version==1, since one could have
     * continued an old format simulation with a new one with mdrun -append.
     */
    if (bRead && ef->eo.bOldFileOpen)
    {
        /* Convert old full simulation sums to sums between energy frames */
        convert_full_sums(&(ef->eo),fr);
    }
    if (fr->ndisre)
    {
        if (bRead && fr->ndisre>fr->d_alloc)
        {
            srenew(fr->disre_rm3tav,fr->ndisre);
            srenew(fr->disre_rt,fr->ndisre);
            fr->d_alloc = fr->ndisre;
        }
        ndo_real(fr->disre_rm3tav,fr->ndisre,bOK1);
        bOK = bOK && bOK1;
        ndo_real(fr->disre_rt,fr->ndisre,bOK1);
        bOK = bOK && bOK1;
    }
    for(block=0; block<fr->nblock; block++)
    {
        if (bRead && fr->nr[block]>fr->b_alloc[block])
        {
            srenew(fr->block[block],fr->nr[block]);
            fr->b_alloc[block] = fr->nr[block];
        }
        ndo_real(fr->block[block],fr->nr[block],bOK1);
        bOK = bOK && bOK1;
    }
    
    if(!bRead)
    {
        if( gmx_fio_flush(ef->fp) != 0)
        {
            gmx_file("Cannot write energy file; maybe you are out of quota?");
        }
    }
    
    if (!bOK)
    {
        if (bRead)
        {
            fprintf(stderr,"\nLast energy frame read %d",
                    ef->framenr-1);
            fprintf(stderr,"\nWARNING: Incomplete energy frame: nr %d time %8.3f\n",
                    ef->framenr,fr->t);
        }
        else
        {
            gmx_fatal(FARGS,"could not write energies");
        }
        return FALSE; 
    }
    
    return TRUE;
}

static real find_energy(const char *name, int nre, gmx_enxnm_t *enm,
                        t_enxframe *fr)
{
    int i;
    
    for(i=0; i<nre; i++)
    {
        if (strcmp(enm[i].name,name) == 0)
        {
            return  fr->ener[i].e;
        }
    }
    
    gmx_fatal(FARGS,"Could not find energy term named '%s'",name);
    
    return 0;
}


void get_enx_state(const char *fn, real t, gmx_groups_t *groups, t_inputrec *ir,
                   t_state *state)
{
  /* Should match the names in mdebin.c */
  static const char *boxvel_nm[] = {
  "Box-Vel-XX", "Box-Vel-YY", "Box-Vel-ZZ",
  "Box-Vel-YX", "Box-Vel-ZX", "Box-Vel-ZY"
  };
  
  static const char *pcouplmu_nm[] = {
    "Pcoupl-Mu-XX", "Pcoupl-Mu-YY", "Pcoupl-Mu-ZZ",
    "Pcoupl-Mu-YX", "Pcoupl-Mu-ZX", "Pcoupl-Mu-ZY"
  };
  int ind0[] = { XX,YY,ZZ,YY,ZZ,ZZ };
  int ind1[] = { XX,YY,ZZ,XX,XX,YY };

  int nre,nfr,i,ni,npcoupl;
  char       buf[STRLEN];
  gmx_enxnm_t *enm;
  t_enxframe *fr;
  ener_file_t in;

  in = open_enx(fn,"r");
  do_enxnms(in,&nre,&enm);
  snew(fr,1);
  nfr = 0;
  while ((nfr==0 || fr->t != t) && do_enx(in,fr)) {
    nfr++;
  }
  close_enx(in);
  fprintf(stderr,"\n");

  if (nfr == 0 || fr->t != t)
    gmx_fatal(FARGS,"Could not find frame with time %f in '%s'",t,fn);
  
  npcoupl = TRICLINIC(ir->compress) ? 6 : 3;
  if (ir->epc == epcPARRINELLORAHMAN) {
    clear_mat(state->boxv);
    for(i=0; i<npcoupl; i++) {
      state->boxv[ind0[i]][ind1[i]] =
	find_energy(boxvel_nm[i],nre,enm,fr);
    }
    fprintf(stderr,"\nREAD %d BOX VELOCITIES FROM %s\n\n",npcoupl,fn);
  }

  if (ir->etc == etcNOSEHOOVER) {
    for(i=0; i<state->ngtc; i++) {
      ni = groups->grps[egcTC].nm_ind[i];
      sprintf(buf,"Xi-%s",*(groups->grpname[ni]));
      state->nosehoover_xi[i] = find_energy(buf,nre,enm,fr);
    }
    fprintf(stderr,"\nREAD %d NOSE-HOOVER Xi's FROM %s\n\n",state->ngtc,fn);
  }
  
  free_enxnms(nre,enm);
  free_enxframe(fr);
  sfree(fr);
}

