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
#ifndef _NBKERNEL232_POWER6_H_
#define _NBKERNEL232_POWER6_H_

/*! \file  nb_kernel232_power6.h
 *  \brief Nonbonded kernel 232 (RF Coul + Tab VdW, SPC-SPC)
 *
 *  \internal
 */

#include "types/simple.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif


/*! \brief Nonbonded kernel 232 with forces.
 *
 *  <b>Coulomb interaction:</b> Reaction-Field <br>
 *  <b>VdW interaction:</b> Tabulated  <br>
 *  <b>Water optimization:</b> SPC - SPC <br>
 *  <b>Forces calculated:</b> Yes <br>
 *
 *  \note All level1 and level2 nonbonded kernels use the same
 *        call sequence. Parameters are documented in nb_kernel_power6.h
 */
void
nb_kernel232_power6
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


/*! \brief Nonbonded kernel 232 without forces.
 *
 *  <b>Coulomb interaction:</b> Reaction-Field <br>
 *  <b>VdW interaction:</b> Tabulated  <br>
 *  <b>Water optimization:</b> SPC - SPC <br>
 *  <b>Forces calculated:</b> No <br>
 *
 *  \note All level1 and level2 nonbonded kernels use the same
 *        call sequence. Parameters are documented in nb_kernel_power6.h
 */
void
nb_kernel232nf_power6
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

#endif /* _NBKERNEL232_POWER6_H_ */
