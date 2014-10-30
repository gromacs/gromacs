/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2005,2006,2007,2008,2009,2010,2011,2012,2013,2014, by the GROMACS development team, led by
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

#include "gmxpre.h"

#include "domdec_utility.h"

#include <assert.h>
#include <string.h>

#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "domdec_internal.h"


void dd_set_cginfo_threaded(int *index_gl, int cg0, int cg1,
                            t_forcerec *fr, char *bLocalCG)
{
    cginfo_mb_t *cginfo_mb;
    int         *cginfo;
    int          cg;

    cginfo_mb = fr->cginfo_mb;
    cginfo    = fr->cginfo;

    /* Looking up the cginfo for all local/non-local atoms is expensive,
     * even when there is only one molecule block, so we multi-thread it.
     */
    int nthread gmx_unused;

    // cppcheck-suppress unreadVariable
    nthread = gmx_omp_nthreads_get(emntDomdec);
#pragma omp parallel for num_threads(nthread) schedule(static)
    for (cg = cg0; cg < cg1; cg++)
    {
        cginfo[cg] = ddcginfo(cginfo_mb, index_gl[cg]);
    }

    if (bLocalCG != NULL)
    {
        for (cg = cg0; cg < cg1; cg++)
        {
            bLocalCG[index_gl[cg]] = TRUE;
        }
    }
}
