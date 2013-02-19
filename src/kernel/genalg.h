/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#ifndef _genalg_h
#define _genalg_h

typedef struct {
    int  np;         /* Number of points           */
    int  atype;      /* Atom type                  */
    int  ptype;      /* Parameter type             */
    real rmin, rmax; /* Minimum and maximum value  */
    real dr;
    real rval;       /* Current value              */
} t_range;

typedef struct {
    int     NP, D;
    int     strategy;
    int     seed;
    int     ipop, gen;    /* Current population member and generation */
    int     imin;         /* Member with lowest energy */
    real    CR, FF;
    real  **pold, **pnew; /* Old and new populations */
    real   *best, *bestit, *cost, *tmp, *msf, *energy;
    tensor *pres;
    rvec   *scale;
} t_genalg;

enum {
    eseSIGMA, eseEPSILON, eseBHAMA, eseBHAMB, eseBHAMC,
    eseCELLX, eseCELLY, eseCELLZ, eseNR
};

extern real value_rand(t_range *r, int *seed);

extern t_genalg *init_ga(FILE *fplog, const char *infile, int D, t_range range[]);

extern void update_ga(FILE *fpout_ptr, t_range range[], t_genalg *ga);

extern gmx_bool print_ga(FILE *fp, t_genalg *ga, real msf, tensor pres, rvec scale,
                         real energy, t_range range[], real tol);

extern real cost(tensor P, real MSF, real energy);

#endif
