/*
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
 * GRoups of Organic Molecules in ACtion for Science
 */


#include "simple.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

typedef struct {
    real    Th;             /* Temperature at half step        */
    real    T;              /* Temperature at full step        */
    tensor  ekinh;          /* Kinetic energy at half step     */
    tensor  ekinh_old;      /* Kinetic energy at old half step */
    tensor  ekinf;          /* Kinetic energy at full step     */
    real    lambda;         /* Berendsen coupling lambda       */
    double  ekinscalef_nhc; /* Scaling factor for NHC- full step */
    double  ekinscaleh_nhc; /* Scaling factor for NHC- half step */
    double  vscale_nhc;     /* Scaling factor for NHC- velocity */
} gmx_temperature_coupling_group_outputs_t;

typedef struct {
    int                                       ngroups;    /* The number of T-coupling groups */
    gmx_temperature_coupling_group_outputs_t *group_data; /* T-coupling data */
} gmx_temperature_coupling_outputs_t;

typedef struct {
    int     nat;    /* Number of atoms in this group		*/
    rvec    u;      /* Mean velocities of home particles        */
    rvec    uold;   /* Previous mean velocities of home particles   */
    double  mA;     /* Mass for topology A		                */
    double  mB;     /* Mass for topology B		                */
} t_grp_acc;

typedef struct {
    gmx_bool   bDoAcceleration; /* Whether acceleration will occur */
    int        ngroups;         /* The number of acceleration groups */
    t_grp_acc *group_data;      /* Acceleration data			*/
} gmx_constant_acceleration_t;

typedef struct {
    real    accel;  /* The acceleration for the cosine profile      */
    real    mvcos;  /* The cos momenta of home particles            */
    real    vcos;   /* The velocity of the cosine profile           */
} t_cosine_acceleration;

typedef struct {
    tensor         **ekin_work_alloc; /* Allocated locations for *_work members */
    tensor         **ekin_work;       /* Per-thread work buffers for accumulating tcstat */
    real           **dekindl_work;    /* Per-thread work buffers for accumulating dekindl */
    tensor           ekin;            /* overall kinetic energy               */
    tensor           ekinh;           /* overall 1/2 step kinetic energy      */
    real             dekindl;         /* dEkin/dlambda at half step           */
    real             dekindl_old;     /* dEkin/dlambda at old half step       */
} gmx_ekindata_t;

#define GID(igid, jgid, gnr) ((igid < jgid) ? (igid*gnr+jgid) : (jgid*gnr+igid))

#ifdef __cplusplus
}
#endif
