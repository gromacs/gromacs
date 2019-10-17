/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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
#ifndef _nb_kernel_h_
#define _nb_kernel_h_

#include <stdio.h>

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/utility/real.h"

struct t_blocka;

/* Structure to collect kernel data not available in forcerec or mdatoms structures.
 * This is only used inside the nonbonded module.
 */
typedef struct
{
    int                    flags;
    const struct t_blocka* exclusions;
    real*                  lambda;
    real*                  dvdl;

    /* pointers to tables */
    t_forcetable* table_elec;
    t_forcetable* table_vdw;
    t_forcetable* table_elec_vdw;

    /* potentials */
    real* energygrp_elec;
    real* energygrp_vdw;
} nb_kernel_data_t;


typedef void nb_kernel_t(t_nblist* gmx_restrict nlist,
                         rvec* gmx_restrict x,
                         rvec* gmx_restrict f,
                         struct t_forcerec* gmx_restrict fr,
                         t_mdatoms* gmx_restrict mdatoms,
                         nb_kernel_data_t* gmx_restrict kernel_data,
                         t_nrnb* gmx_restrict nrnb);

#endif /* _nb_kernel_h_ */
