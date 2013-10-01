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

#include <vector>

#include "tngio.h"
#include "tng_io.h"
#include "trx.h"

#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/legacyheaders/physics.h"
#include "gromacs/legacyheaders/gmx_fatal.h"

struct tng_trajectory_output
{
    tng_trajectory_t in, out;
};

void prepare_tng_writing(const char              *filename,
                         char                     mode,
                         tng_trajectory_output_t *tngptr,
                         int                      natoms)
{
#ifdef GMX_USE_TNG

    snew(*tngptr, 1);
    tng_trajectory_output_t tng = *tngptr;

    tng_open(filename, mode, &tng->out);
    tng_trajectory_t output = tng->out;

    tng_implicit_num_particles_set(output, natoms);

    tng_trajectory_t input = tng->in;
    if (input)
    {
        /* Set parameters (time per frame, n frames per frame set and
         * writing intervals of positions, box shape and lambdas) of
         * the output tng container based on their respective values
         * int the input tng container */

        // TODO copy the topology from tngin if it exists
        int64_t n_frames_per_frame_set, interval = -1;
        double  time;
        // TODO make this configurable
        //    char compression = bUseLossyCompression ? TNG_TNG_COMPRESSION : TNG_GZIP_COMPRESSION;
        char compression = TNG_TNG_COMPRESSION;

        tng_time_per_frame_get(input, &time);
        tng_time_per_frame_set(output, time);

        if (tng_data_get_stride_length(input, TNG_TRAJ_POSITIONS, -1, &interval)
            == TNG_SUCCESS)
        {
            tng_util_generic_write_interval_set(output, interval, 3, TNG_TRAJ_POSITIONS,
                                                "POSITIONS", TNG_PARTICLE_BLOCK_DATA,
                                                compression);
            tng_util_generic_write_interval_set(output, interval, 9, TNG_TRAJ_BOX_SHAPE,
                                                "BOX SHAPE", TNG_NON_PARTICLE_BLOCK_DATA,
                                                TNG_GZIP_COMPRESSION);
            tng_util_generic_write_interval_set(output, interval, 1, TNG_GMX_LAMBDA,
                                                "LAMBDAS", TNG_NON_PARTICLE_BLOCK_DATA,
                                                TNG_GZIP_COMPRESSION);
        }
        tng_num_frames_per_frame_set_get(input, &n_frames_per_frame_set);
        tng_num_frames_per_frame_set_set(output, n_frames_per_frame_set);
    }
    else
    {
        tng_num_frames_per_frame_set_set(output, 1);
    }
#endif
    // TODO support more functionality
}

void open_tng_for_reading(const char              *fn,
                          tng_trajectory_output_t *tng)
{
    snew(*tng, 1);
    tng_open(fn, 'r', &(*tng)->in);
}

void write_tng_from_trxframe(tng_trajectory_output_t tng,
                             t_trxframe             *frame)
{
#ifdef GMX_USE_TNG
    if(frame->step > 0)
    {
        double time_per_frame = frame->time * PICO / frame->step;
        tng_time_per_frame_set(tng->out, time_per_frame);
    }
    // TODO what is gc?
    fwrite_tng(tng->out,
               // TODO make compression configurable
               TRUE,
               frame->step,
               frame->time,
               0,
               (const rvec *) frame->box,
               frame->natoms,
               (const rvec *) frame->x,
               (const rvec *) frame->v,
               (const rvec *) frame->f);
#endif
}

void tng_tools_close(tng_trajectory_output_t tng)
{
#ifdef GMX_USE_TNG
    tng_close(&tng->in);
    tng_close(&tng->out);
#endif
}

static void
tng_convert_array_to_real_array(void       *from,
                                real       *to,
                                const float fact,
                                const int   natoms,
                                const int   nvalues,
                                const char  datatype)
{
    int i, j;

    switch (datatype)
    {
        case TNG_FLOAT_DATA:
            if (sizeof(real) == sizeof(float))
            {
                if (fact == 1)
                {
                    memcpy(to, from, nvalues * sizeof(real) * natoms);
                }
                else
                {
                    for (i = 0; i < natoms; i++)
                    {
                        for (j = 0; j < nvalues; j++)
                        {
                            to[i*nvalues+j] = (real)((float *)from)[i*nvalues+j] * fact;
                        }
                    }
                }
            }
            else
            {
                for (i = 0; i < natoms; i++)
                {
                    for (j = 0; j < nvalues; j++)
                    {
                        to[i*nvalues+j] = (real)((float *)from)[i*nvalues+j] * fact;
                    }
                }
            }
            break;
        case TNG_INT_DATA:
            for (i = 0; i < natoms; i++)
            {
                for (j = 0; j < nvalues; j++)
                {
                    to[i*nvalues+j] = (real)((int64_t *)from)[i*nvalues+j] * fact;
                }
            }
            break;
        case TNG_DOUBLE_DATA:
            if (sizeof(real) == sizeof(double))
            {
                if (fact == 1)
                {
                    memcpy(to, from, nvalues * sizeof(real) * natoms);
                }
                else
                {
                    for (i = 0; i < natoms; i++)
                    {
                        for (j = 0; j < nvalues; j++)
                        {
                            to[i*nvalues+j] = (real)((double *)from)[i*nvalues+j] * fact;
                        }
                    }
                }
            }
            else
            {
                for (i = 0; i < natoms; i++)
                {
                    for (j = 0; j < nvalues; j++)
                    {
                        to[i*nvalues+j] = (real)((double *)from)[i*nvalues+j] * fact;
                    }
                }
            }
            break;
        default:
            gmx_incons("Illegal datatype when converting values to a real array!");
            return;
    }
    return;
}

real getDistanceScaleFactor(tng_trajectory_t in)
{
    int64_t exp = -1;
    real    distanceScaleFactor;

    // TODO Ideally, TNG can do this for us
    tng_distance_unit_exponential_get(in, &exp);

    // GROMACS expects distances in nm
    switch (exp)
    {
        case 9:
            distanceScaleFactor = NANO/NANO;
            break;
        case 10:
            distanceScaleFactor = NANO/ANGSTROM;
            break;
        default:
            distanceScaleFactor = pow(10.0, exp + 9.0);
    }

    return distanceScaleFactor;
}

gmx_bool read_next_tng_frame(tng_trajectory_output_t tng,
                             t_trxframe             *fr)
{
#ifdef GMX_USE_TNG
    gmx_bool            bOK = TRUE;
    tng_function_status stat;
    int64_t             numberOfAtoms = -1, frameNumber = -1;
    char                datatype      = -1;
    void               *values        = NULL;
    double              frameTime     = -1.0;

    stat = tng_num_particles_get(tng->in, &numberOfAtoms);
    if (stat != TNG_SUCCESS)
    {
        gmx_file("Cannot determine number of atoms from TNG file.");
    }
    fr->natoms = numberOfAtoms;

    stat = tng_util_particle_data_next_frame_read(tng->in,
                                                  TNG_TRAJ_POSITIONS,
                                                  &values,
                                                  &datatype,
                                                  &frameNumber,
                                                  &frameTime);
    if (stat == TNG_CRITICAL)
    {
        if (values)
        {
            free(values);
        }
        gmx_file("Cannot read positions from TNG file.");
        return FALSE;
    }
    else if (stat == TNG_FAILURE)
    {
        // TODO what should we be doing here?
        return FALSE;
    }

    fr->step = (int) frameNumber;
    // Convert the time to ps
    fr->time = frameTime / PICO;

    srenew(fr->x, fr->natoms);
    tng_convert_array_to_real_array(values,
                                    (real *) fr->x,
                                    getDistanceScaleFactor(tng->in),
                                    fr->natoms,
                                    DIM,
                                    datatype);

    /* values does not have to be freed before reading next frame. It will
     * be reallocated if it is not NULL. */

    // TODO how do we fill lambda, etc. (if we need to? check GROMACS
    // client code)

    stat = tng_util_non_particle_data_next_frame_read(tng->in,
                                                      TNG_TRAJ_BOX_SHAPE,
                                                      &values,
                                                      &datatype,
                                                      &frameNumber,
                                                      &frameTime);
    if (stat != TNG_SUCCESS)
    {
        gmx_file("Cannot read box shape from TNG file.");
    }
    tng_convert_array_to_real_array(values,
                                    (real *) fr->box,
                                    getDistanceScaleFactor(tng->in),
                                    DIM,
                                    DIM,
                                    datatype);

    /* values must be freed before leaving this function */
    free(values);

    return bOK;
#else
    return FALSE;
#endif
}
