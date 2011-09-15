/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2011
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

#ifndef _INTERACTION_CONST_
#define _INTERACTION_CONST_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    /* VdW */
    real rvdw;

    /* Cut-off */
    real rlist;

    /* PME/Ewald */
    real ewaldcoeff;
  
    /* Dielectric constant resp. multiplication factor for charges */
    real epsilon_r,epsilon_rf,epsfac;  
  
    /* Constants for reaction fields */
    real k_rf,c_rf;

    /* type of electrostatics (defined in enums.h) */
    int  eeltype;

    /* Force/energy interpolation tables, linear in force, quadratic in V */
    real tabq_scale;
    int  tabq_size;
    /* Coulomb table, size of array is tabsize*4 (when used) */
    real *tabq_coul_FDV0;
} interaction_const_t;

#ifdef __cplusplus
}
#endif

#endif /* _INTERACTION_CONST_ */
