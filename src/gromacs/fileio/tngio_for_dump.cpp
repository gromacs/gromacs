/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
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
#include "tngio_for_tools.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "tngio.h"

#ifdef GMX_USE_TNG
#include "../../external/tng_io/include/tng_io.h"
#endif

#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/tools/dump.h"

void list_tng_for_gmx_dump(const char *fn, gmx_bool bXVG)
{
#ifdef GMX_USE_TNG
    tng_trajectory_t     tng;
    gmx_int64_t          nframe = 0;
    gmx_int64_t          i, *block_ids = NULL, step, ndatablocks;
    gmx_int64_t          pos_block_id = TNG_TRAJ_POSITIONS;
    gmx_bool             bOK;

    gmx_tng_open(fn, 'r', &tng);
    gmx_print_tng_molecule_system(tng, stdout);

    bOK    = gmx_get_tng_data_block_types_of_next_frame(tng, -1,
                                                        bXVG ? 1 : 0,
                                                        bXVG ? &pos_block_id : NULL,
                                                        &step, &ndatablocks,
                                                        &block_ids);
    do
    {
        for (i = 0; i < ndatablocks; i++)
        {
            double               frame_time;
            real                 prec, *values = NULL;
            gmx_int64_t          n_values_per_frame, n_atoms;
            char                 block_name[STRLEN];

            gmx_get_tng_data_next_frame_of_block_type(tng, block_ids[i], &values,
                                                      &step, &frame_time,
                                                      &n_values_per_frame, &n_atoms,
                                                      &prec,
                                                      block_name, STRLEN, &bOK);
            if (!bOK)
            {
                /* Can't write any output because we don't know what
                   arrays are valid. */
                fprintf(stderr, "\nWARNING: Incomplete frame at time %g, will not write output\n", frame_time);
                list_tng_inner(fn, (0 == i), bXVG, values, step, frame_time,
                               n_values_per_frame, n_atoms, prec, nframe, block_name);
            }
        }
        nframe++;
    }
    while (gmx_get_tng_data_block_types_of_next_frame(tng, step,
                                                      bXVG ? 1 : 0,
                                                      bXVG ? &pos_block_id : NULL,
                                                      &step,
                                                      &ndatablocks,
                                                      &block_ids));

    if (block_ids)
    {
        sfree(block_ids);
    }

    gmx_tng_close(&tng);
#else
    GMX_UNUSED_VALUE(fn);
    GMX_UNUSED_VALUE(bXVG);
#endif
}
