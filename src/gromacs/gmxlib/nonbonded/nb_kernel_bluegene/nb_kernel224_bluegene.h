/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2004
 * David van der Spoel, Erik Lindahl, University of Groningen.
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
#ifndef _NB_KERNEL224_H_
#define _NB_KERNEL224_H_

/* This header is never installed, so we can use config.h */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <types/simple.h>

/*! \file  nb_kernel224.h
 *  \brief Nonbonded kernel 224 (RF Coul + Bham, TIP4p-TIP4p)
 *
 *  \internal
 */

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif


/*! \brief Nonbonded kernel 224 with forces.
 *
 *  \internal  Generated at compile time in either C or Fortran
 *             depending on configuration settings. The name of
 *             the function in C is nb_kernel224. For Fortran the
 *             name mangling depends on the compiler, but in Gromacs
 *             you can handle it automatically with the macro
 *             F77_OR_C_FUNC_(nb_kernel224,NB_KERNEL224), which
 *             expands to the correct identifier.
 *
 *  <b>Coulomb interaction:</b> Reaction-Field <br>
 *  <b>VdW interaction:</b> Buckingham <br>
 *  <b>Water optimization:</b> TIP4p - TIP4p <br>
 *  <b>Forces calculated:</b> Yes <br>
 *
 *  \note All level1 and level2 nonbonded kernels use the same
 *        call sequence. Parameters are documented in nb_kernel.h
 */
void
nb_kernel224_bluegene
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
                 real          work[]);


/*! \brief Nonbonded kernel 224 without forces.
 *
 *  \internal  Generated at compile time in either C or Fortran
 *             depending on configuration settings. The name of
 *             the function in C is nb_kernel224. For Fortran the
 *             name mangling depends on the compiler, but in Gromacs
 *             you can handle it automatically with the macro
 *             F77_OR_C_FUNC_(nb_kernel224,NB_KERNEL224), which
 *             expands to the correct identifier.
 *
 *  <b>Coulomb interaction:</b> Reaction-Field <br>
 *  <b>VdW interaction:</b> Buckingham <br>
 *  <b>Water optimization:</b> TIP4p - TIP4p <br>
 *  <b>Forces calculated:</b> No <br>
 *
 *  \note All level1 and level2 nonbonded kernels use the same
 *        call sequence. Parameters are documented in nb_kernel.h
 */
void
nb_kernel224nf_bluegene
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
                 real          work[]);


#ifdef __cplusplus
}
#endif

#endif /* _NB_KERNEL224_H_ */
