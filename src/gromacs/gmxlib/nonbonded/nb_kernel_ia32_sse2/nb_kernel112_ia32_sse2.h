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
#ifndef _NB_KERNEL112_IA32_SSE2_H_
#define _NB_KERNEL112_IA32_SSE2_H_


/*! \file  nb_kernel112_ia32_sse2.h
 *  \brief ia32 SSE2-optimized versions of nonbonded kernel 112
 *
 *  \internal
 */

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/*! \brief Nonbonded kernel 112 with forces, optimized for ia32 sse2.
 *
 *  \internal
 *
 *  <b>Coulomb interaction:</b> Standard 1/r <br>
 *  <b>VdW interaction:</b> Lennard-Jones <br>
 *  <b>Water optimization:</b> Pairs of SPC/TIP3P waters interaction <br>
 *  <b>Forces calculated:</b> Yes <br>
 *
 *  \note All level1 and level2 nonbonded kernels use the same
 *        call sequence. Parameters are documented in nb_kernel.h
 */
void
nb_kernel112_ia32_sse2  (int *    nri,        int      iinr[],    int      jindex[],
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




/*! \brief Nonbonded kernel 112 without forces, optimized for ia32 sse2.
 *
 *  \internal
 *
 *  <b>Coulomb interaction:</b> Standard 1/r <br>
 *  <b>VdW interaction:</b> Lennard-Jones <br>
 *  <b>Water optimization:</b> Pairs of SPC/TIP3P waters interaction <br>
 *  <b>Forces calculated:</b> No <br>
 *
 *  \note All level1 and level2 nonbonded kernels use the same
 *        call sequence. Parameters are documented in nb_kernel.h
 */
void
nb_kernel112nf_ia32_sse2(int *    nri,        int      iinr[],    int      jindex[],
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


#endif /* _NB_KERNEL112_IA32_SSE2_H_ */
