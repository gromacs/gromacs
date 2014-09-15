/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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
#ifndef _fcdata_h
#define _fcdata_h

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef real rvec5[5];

/* Distance restraining stuff */
typedef struct {
    int      dr_weighting; /* Weighting of pairs in one restraint              */
    gmx_bool dr_bMixed;    /* Use sqrt of the instantaneous times              *
                            * the time averaged violation                      */
    real     dr_fc;        /* Force constant for disres,                       *
                            * which is multiplied by a (possibly)              *
                            * different factor for each restraint              */
    real  dr_tau;          /* Time constant for disres		          */
    real  ETerm;           /* multiplication factor for time averaging         */
    real  ETerm1;          /* 1 - ETerm1                                       */
    real  exp_min_t_tau;   /* Factor for slowly switching on the force         */
    int   nres;            /* The number of distance restraints                */
    int   npair;           /* The number of distance restraint pairs           */
    real  sumviol;         /* The sum of violations                            */
    real *rt;              /* The calculated instantaneous distance (npr)      */
    real *rm3tav;          /* The calculated time averaged distance (npr)      */
    real *Rtl_6;           /* The calculated instantaneous r^-6 (nr)           */
    real *Rt_6;            /* The calculated inst. ens. averaged r^-6 (nr)     */
    real *Rtav_6;          /* The calculated time and ens. averaged r^-6 (nr)  */
    int   nsystems;        /* The number of systems for ensemble averaging     */
} t_disresdata;


/* Orientation restraining stuff */
typedef struct {
    real      fc;            /* Force constant for the restraints                  */
    real      edt;           /* Multiplication factor for time averaging           */
    real      edt_1;         /* 1 - edt                                            */
    real      exp_min_t_tau; /* Factor for slowly switching on the force         */
    int       nr;            /* The number of orientation restraints               */
    int       nex;           /* The number of experiments                          */
    int       nref;          /* The number of atoms for the fit                    */
    real     *mref;          /* The masses of the reference atoms                  */
    rvec     *xref;          /* The reference coordinates for the fit (nref)       */
    rvec     *xtmp;          /* Temporary array for fitting (nref)                 */
    matrix    R;             /* Rotation matrix to rotate to the reference coor.   */
    tensor   *S;             /* Array of order tensors for each experiment (nexp)  */
    rvec5    *Dinsl;         /* The order matrix D for all restraints (nr x 5)     */
    rvec5    *Dins;          /* The ensemble averaged D (nr x 5)                   */
    rvec5    *Dtav;          /* The time and ensemble averaged D (nr x 5)          */
    real     *oinsl;         /* The calculated instantaneous orientations          */
    real     *oins;          /* The calculated emsemble averaged orientations      */
    real     *otav;          /* The calculated time and ensemble averaged orient.  */
    real      rmsdev;        /* The weighted (using kfac) RMS deviation            */
    rvec5    *tmp;           /* An array of temporary 5-vectors (nex);             */
    real   ***TMP;           /* An array of temporary 5x5 matrices (nex);          */
    real     *eig;           /* Eigenvalues/vectors, for output only (nex x 12)    */

    /* variables for diagonalization with diagonalize_orires_tensors()*/
    double **M;
    double  *eig_diag;
    double **v;
} t_oriresdata;

typedef struct {
    int   n;      /* n+1 is the number of points */
    real  scale;  /* distance between two points */
    real *data;   /* the actual table data, per point there are 4 numbers */
} bondedtable_t;

/*
 * Data struct used in the force calculation routines
 * for storing the tables for bonded interactions and
 * for storing information which is needed in following steps
 * (for instance for time averaging in distance retraints)
 * or for storing output, since force routines only return the potential.
 */
typedef struct {
    bondedtable_t *bondtab;
    bondedtable_t *angletab;
    bondedtable_t *dihtab;

    t_disresdata   disres;
    t_oriresdata   orires;
} t_fcdata;

#ifdef __cplusplus
}
#endif

#endif /* _fcdata_h */
