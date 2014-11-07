/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef GMX_LEGACYHEADERS_TYPES_GROUP_H
#define GMX_LEGACYHEADERS_TYPES_GROUP_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
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
} t_grp_tcstat;

typedef struct {
    int     nat;    /* Number of atoms in this group		*/
    rvec    u;      /* Mean velocities of home particles        */
    rvec    uold;   /* Previous mean velocities of home particles   */
    double  mA;     /* Mass for topology A		                */
    double  mB;     /* Mass for topology B		                */
} t_grp_acc;

typedef struct {
    real    cos_accel;  /* The acceleration for the cosine profile      */
    real    mvcos;      /* The cos momenta of home particles            */
    real    vcos;       /* The velocity of the cosine profile           */
} t_cos_acc;

typedef struct {
    gmx_bool         bNEMD;
    int              ngtc;            /* The number of T-coupling groups      */
    t_grp_tcstat    *tcstat;          /* T-coupling data            */
    tensor         **ekin_work_alloc; /* Allocated locations for *_work members */
    tensor         **ekin_work;       /* Work arrays for tcstat per thread    */
    real           **dekindl_work;    /* Work location for dekindl per thread */
    int              ngacc;           /* The number of acceleration groups    */
    t_grp_acc       *grpstat;         /* Acceleration data			*/
    tensor           ekin;            /* overall kinetic energy               */
    tensor           ekinh;           /* overall 1/2 step kinetic energy      */
    real             dekindl;         /* dEkin/dlambda at half step           */
    real             dekindl_old;     /* dEkin/dlambda at old half step       */
    t_cos_acc        cosacc;          /* Cosine acceleration data             */
} gmx_ekindata_t;

#define GID(igid, jgid, gnr) ((igid < jgid) ? (igid*gnr+jgid) : (jgid*gnr+igid))

#ifdef __cplusplus
}
#endif

#endif
