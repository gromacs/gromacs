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

/* The gmx_bools indicate whether a field was read from the trajectory.
 * Do not try to use a pointer when its gmx_bool is FALSE, as memory might
 * not be allocated.
 */

#ifndef GMX_FILEIO_TRX_H
#define GMX_FILEIO_TRX_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

struct t_atoms;

typedef struct gmxvmdplugin t_gmxvmdplugin;

typedef struct t_trxframe
{
    int      flags;            /* flags for read_first/next_frame  */
    int      not_ok;           /* integrity flags                  */
    gmx_bool bDouble;          /* Double precision?                */
    int      natoms;           /* number of atoms (atoms, x, v, f) */
    real     t0;               /* time of the first frame, needed  *
                                * for skipping frames with -dt     */
    real     tf;               /* internal frame time - DO NOT CHANGE */
    real     tpf;              /* time of the previous frame, not  */
                               /* the read, but real file frames   */
    real     tppf;             /* time of two frames ago           */
                               /* tpf and tppf are needed to       */
                               /* correct rounding errors for -e   */
    gmx_bool        bTitle;
    const char     *title;     /* title of the frame            */
    gmx_bool        bStep;
    int             step;      /* MD step number                   */
    gmx_bool        bTime;
    real            time;      /* time of the frame                */
    gmx_bool        bLambda;
    gmx_bool        bFepState; /* does it contain fep_state?       */
    real            lambda;    /* free energy perturbation lambda  */
    int             fep_state; /* which fep state are we in? */
    gmx_bool        bAtoms;
    struct t_atoms *atoms;     /* atoms struct (natoms)            */
    gmx_bool        bPrec;
    real            prec;      /* precision of x, fraction of 1 nm */
    gmx_bool        bX;
    rvec           *x;         /* coordinates (natoms)             */
    gmx_bool        bV;
    rvec           *v;         /* velocities (natoms)              */
    gmx_bool        bF;
    rvec           *f;         /* forces (natoms)                  */
    gmx_bool        bBox;
    matrix          box;       /* the 3 box vectors                */
    gmx_bool        bPBC;
    int             ePBC;      /* the type of pbc                  */
    t_gmxvmdplugin* vmdplugin;
} t_trxframe;

#ifdef __cplusplus
}
#endif

#endif
