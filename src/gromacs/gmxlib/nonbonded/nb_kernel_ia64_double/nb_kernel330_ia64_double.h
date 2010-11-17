/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-3304
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
#ifndef _NB_KERNEL330_IA64D_H_
#define _NB_KERNEL330_IA64D_H_

/*! \file  nb_kernel330_ia64_double.h
 *  \brief ia64(double)-optimized versions of nonbonded kernel 330
 *
 *  \internal
 */

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/*! \brief Nonbonded kernel 330 with forces, optimized for ia64 double precision assembly.
 *
 *  \internal
 *
 *  <b>Coulomb interaction:</b> Tabulated <br>
 *  <b>VdW interaction:</b> Tabulated <br>
 *  <b>Water optimization:</b> No <br>
 *  <b>Forces calculated:</b> Yes <br>
 *
 *  \note All level1 and level2 nonbonded kernels use the same
 *        call sequence. Parameters are documented in nb_kernel.h
 */
void
nb_kernel330_ia64_double  (int *    nri,        int      iinr[],    int      jindex[],
                     int      jjnr[],     int      shift[],   double   shiftvec[],
                     double   fshift[],   int      gid[],     double   pos[],
                     double   faction[],  double   charge[],  double * facel,
                     double * krf,        double * crf,       double   Vc[],
                     int      type[],     int *    ntype,     double   vdwparam[],
                     double   Vvdw[],     double * tabscale,  double   VFtab[],
                     double   invsqrta[], double   dvda[],    double * gbtabscale,
                     double   GBtab[],    int *    nthreads,  int *    count,
                     void *   mtx,        int *    outeriter, int *    inneriter,
                     double * work);

#ifdef __cplusplus
}
#endif

#endif /* _NB_KERNEL330_IA64D_H_ */
