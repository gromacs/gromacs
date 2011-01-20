/*
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 *
 * Gromacs is a library for molecular simulation and trajectory analysis,
 * written by Erik Lindahl, David van der Spoel, Berk Hess, and others - for
 * a full list of developers and information, check out http://www.gromacs.org
 *
 * This program is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the Free 
 * Software Foundation; either version 2 of the License, or (at your option) any 
 * later version.
 * As a special exception, you may use this file as part of a free software
 * library without restriction.  Specifically, if other files instantiate
 * templates or use macros or inline functions from this file, or you compile
 * this file and link it with other files to produce an executable, this
 * file does not by itself cause the resulting executable to be covered by
 * the GNU Lesser General Public License.  
 *
 * In plain-speak: do not worry about classes/macros/templates either - only
 * changes to the library have to be LGPL, not an application linking with it.
 *
 * To help fund GROMACS development, we humbly ask that you cite
 * the papers people have written on it - you can find them on the website!
 */

/*! \file  nb_kernel211_power6.c
 *  \brief Wrapper for fortran nonbonded kernel 211
 *
 *  \internal
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <types/simple.h>

#ifndef F77_FUNC
#define F77_FUNC(name,NAME) name ## _
#endif

/* Declarations of Fortran routines. 
 * We avoid using underscores in F77 identifiers for portability! 
 */
void
F77_FUNC(pwr6kernel211,PWR6KERNEL211)
     (int *         nri,        int *         iinr,     
      int *         jindex,     int *         jjnr,   
      int *         shift,      real *        shiftvec,
      real *        fshift,     int *         gid, 
      real *        pos,        real *        faction,
      real *        charge,     real *        facel,
      real *        krf,        real *        crf,  
      real *        Vc,         int *         type,   
      int *         ntype,      real *        vdwparam,
      real *        Vvdw,       real *        tabscale,
      real *        VFtab,      real *        invsqrta, 
      real *        dvda,       real *        gbtabscale,
      real *        GBtab,      int *         nthreads, 
      int *         count,      void *        mtx,
      int *         outeriter,  int *         inneriter,
      real *        work);
     
void
F77_FUNC(pwr6kernel211nf,PWR6KERNEL211NF)
     (int *         nri,        int *         iinr,     
      int *         jindex,     int *         jjnr,   
      int *         shift,      real *        shiftvec,
      real *        fshift,     int *         gid, 
      real *        pos,        real *        faction,
      real *        charge,     real *        facel,
      real *        krf,        real *        crf,  
      real *        Vc,         int *         type,   
      int *         ntype,      real *        vdwparam,
      real *        Vvdw,       real *        tabscale,
      real *        VFtab,      real *        invsqrta, 
      real *        dvda,       real *        gbtabscale,
      real *        GBtab,      int *         nthreads, 
      int *         count,      void *        mtx,
      int *         outeriter,  int *         inneriter,
      real *        work);
     

void
nb_kernel211_power6
     (int *         nri,        int *         iinr,     
      int *         jindex,     int *         jjnr,   
      int *         shift,      real *        shiftvec,
      real *        fshift,     int *         gid, 
      real *        pos,        real *        faction,
      real *        charge,     real *        facel,
      real *        krf,        real *        crf,  
      real *        Vc,         int *         type,   
      int *         ntype,      real *        vdwparam,
      real *        Vvdw,       real *        tabscale,
      real *        VFtab,      real *        invsqrta, 
      real *        dvda,       real *        gbtabscale,
      real *        GBtab,      int *         nthreads, 
      int *         count,      void *        mtx,
      int *         outeriter,  int *         inneriter,
      real *        work)
{
  F77_FUNC(pwr6kernel211,PWR6KERNEL211)
    (nri,iinr,jindex,jjnr,shift,shiftvec,fshift,gid,pos,faction,
     charge,facel,krf,crf,Vc,type,ntype,vdwparam,Vvdw,tabscale,
     VFtab,invsqrta,dvda,gbtabscale,GBtab,nthreads,count,mtx,
     outeriter,inneriter,work);
}



void
nb_kernel211nf_power6
     (int *         nri,        int           iinr[],     
      int           jindex[],   int           jjnr[],   
      int           shift[],    real          shiftvec[],
      real          fshift[],   int           gid[], 
      real          pos[],      real          faction[],
      real          charge[],   real *        facel,
      real *        krf,        real *        crf,  
      real          Vc[],       int           type[],   
      int *         ntype,      real          vdwparam[],
      real          Vvdw[],     real *        tabscale,
      real          VFtab[],    real          invsqrta[], 
      real          dvda[],     real *        gbtabscale,
      real          GBtab[],    int *         nthreads, 
      int *         count,      void *        mtx,
      int *         outeriter,  int *         inneriter,
      real *        work)
{
  F77_FUNC(pwr6kernel211nf,PWR6KERNEL211NF)
    (nri,iinr,jindex,jjnr,shift,shiftvec,fshift,gid,pos,faction,
     charge,facel,krf,crf,Vc,type,ntype,vdwparam,Vvdw,tabscale,
     VFtab,invsqrta,dvda,gbtabscale,GBtab,nthreads,count,mtx,
     outeriter,inneriter,work);
}


