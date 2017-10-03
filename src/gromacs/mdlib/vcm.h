/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_VCM_H
#define GMX_MDLIB_VCM_H

#include <stdio.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_groups_t;
struct t_inputrec;

typedef struct {
    rvec      p;        /* Linear momentum                     */
    rvec      x;        /* Center of mass                      */
    rvec      j;        /* Angular momentum                    */
    tensor    i;        /* Moment of inertia                   */
    real      mass;     /* Mass                                */
} t_vcm_thread;

typedef struct {
    int           nr;          /* Number of groups                    */
    int           size;        /* Size of group arrays                */
    int           stride;      /* Stride for thread data              */
    int           mode;        /* One of the enums above              */
    gmx_bool      ndim;        /* The number of dimensions for corr.  */
    real          timeStep;    /* The time step for COMM removal      */
    real         *group_ndf;   /* Number of degrees of freedom        */
    rvec         *group_p;     /* Linear momentum per group           */
    rvec         *group_v;     /* Linear velocity per group           */
    rvec         *group_x;     /* Center of mass per group            */
    rvec         *group_j;     /* Angular momentum per group          */
    rvec         *group_w;     /* Angular velocity (omega)            */
    tensor       *group_i;     /* Moment of inertia per group         */
    real         *group_mass;  /* Mass per group                      */
    char        **group_name;  /* These two are copies to pointers in */
    t_vcm_thread *thread_vcm;  /* Temporary data per thread and group */
} t_vcm;

t_vcm *init_vcm(FILE *fp, gmx_groups_t *groups, const t_inputrec *ir);

/* Do a per group center of mass things */
void calc_vcm_grp(int start, int homenr, t_mdatoms *md,
                  rvec x[], rvec v[], t_vcm *vcm);

/* Set the COM velocity to zero and potentially correct the COM position.
 *
 * With linear modes nullptr can be passed for x.
 */
void do_stopcm_grp(int homenr,
                   const unsigned short *group_id,
                   rvec x[], rvec v[], const t_vcm &vcm);

void check_cm_grp(FILE *fp, t_vcm *vcm, t_inputrec *ir, real Temp_Max);

#endif
