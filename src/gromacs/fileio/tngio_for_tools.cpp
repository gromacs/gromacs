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
#include "trx.h"

#ifdef GMX_USE_TNG
#include "tng_io.h"
#endif

#include "gromacs/legacyheaders/types/atoms.h"
#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/legacyheaders/physics.h"
#include "gromacs/legacyheaders/gmx_fatal.h"
#include "gromacs/utility/common.h"

struct tng_trajectory_output
{
    tng_trajectory_t in, out;
};

void prepare_tng_writing(const char              *filename,
                         char                     mode,
                         tng_trajectory_t        *input,
                         tng_trajectory_t        *output,
                         int                      nAtoms)
{
#ifdef GMX_USE_TNG

    tng_open(filename, mode, output);

    /* Do we have an input file in TNG format? If so, then there's
       more data we can copy over, rather than having to improvise. */
    if (*input)
    {
        /* Set parameters (time per frame, n frames per frame set and
         * writing intervals of positions, box shape and lambdas) of
         * the output tng container based on their respective values
         * int the input tng container */

        tng_molecule_system_copy(*input, *output);

        int64_t n_frames_per_frame_set, interval = -1;
        double  time;
        // TODO make this configurable
        //    char compression = bUseLossyCompression ? TNG_TNG_COMPRESSION : TNG_GZIP_COMPRESSION;
        char compression = TNG_TNG_COMPRESSION;

        tng_time_per_frame_get(*input, &time);
        tng_time_per_frame_set(*output, time);

        tng_num_frames_per_frame_set_get(*input, &n_frames_per_frame_set);
        tng_num_frames_per_frame_set_set(*output, n_frames_per_frame_set);

        if (tng_data_get_stride_length(*input, TNG_TRAJ_POSITIONS, -1, &interval)
            == TNG_SUCCESS)
        {
            tng_util_generic_write_interval_set(*output, interval, 3, TNG_TRAJ_POSITIONS,
                                                "POSITIONS", TNG_PARTICLE_BLOCK_DATA,
                                                compression);
            tng_util_generic_write_interval_set(*output, interval, 9, TNG_TRAJ_BOX_SHAPE,
                                                "BOX SHAPE", TNG_NON_PARTICLE_BLOCK_DATA,
                                                TNG_GZIP_COMPRESSION);
            tng_util_generic_write_interval_set(*output, interval, 1, TNG_GMX_LAMBDA,
                                                "LAMBDAS", TNG_NON_PARTICLE_BLOCK_DATA,
                                                TNG_GZIP_COMPRESSION);
        }
    }
    else
    {
        tng_implicit_num_particles_set(*output, nAtoms);

        tng_num_frames_per_frame_set_set(*output, 1);
    }
#else
    GMX_UNUSED_VALUE(filename);
    GMX_UNUSED_VALUE(mode);
    GMX_UNUSED_VALUE(input);
    GMX_UNUSED_VALUE(output);
    GMX_UNUSED_VALUE(nAtoms);
#endif
    // TODO support more functionality
}

void open_tng_for_reading(const char              *fn,
                          tng_trajectory_t        *input)
{
#ifdef GMX_USE_TNG
    tng_open(fn, 'r', input);
#else
    gmx_file("GROMACS was compiled without TNG support, cannot handle this file type");
    GMX_UNUSED_VALUE(fn);
    GMX_UNUSED_VALUE(input);
#endif
}

void write_tng_from_trxframe(tng_trajectory_t        output,
                             t_trxframe             *frame)
{
#ifdef GMX_USE_TNG
    if (frame->step > 0)
    {
        double timePerFrame = frame->time * PICO / frame->step;
        tng_time_per_frame_set(output, timePerFrame);
    }
    // TODO what is gc?
    fwrite_tng(output,
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
#else
    GMX_UNUSED_VALUE(output);
    GMX_UNUSED_VALUE(frame);
#endif
}

void tng_tools_close(tng_trajectory_t *tng)
{
#ifdef GMX_USE_TNG
    tng_close(tng);
#else
    GMX_UNUSED_VALUE(tng);
#endif
}

#ifdef GMX_USE_TNG
static void
tng_convert_array_to_real_array(void       *from,
                                real       *to,
                                const float fact,
                                const int   nAtoms,
                                const int   nValues,
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
                    memcpy(to, from, nValues * sizeof(real) * nAtoms);
                }
                else
                {
                    for (i = 0; i < nAtoms; i++)
                    {
                        for (j = 0; j < nValues; j++)
                        {
                            to[i*nValues+j] = (real)((float *)from)[i*nValues+j] * fact;
                        }
                    }
                }
            }
            else
            {
                for (i = 0; i < nAtoms; i++)
                {
                    for (j = 0; j < nValues; j++)
                    {
                        to[i*nValues+j] = (real)((float *)from)[i*nValues+j] * fact;
                    }
                }
            }
            break;
        case TNG_INT_DATA:
            for (i = 0; i < nAtoms; i++)
            {
                for (j = 0; j < nValues; j++)
                {
                    to[i*nValues+j] = (real)((int64_t *)from)[i*nValues+j] * fact;
                }
            }
            break;
        case TNG_DOUBLE_DATA:
            if (sizeof(real) == sizeof(double))
            {
                if (fact == 1)
                {
                    memcpy(to, from, nValues * sizeof(real) * nAtoms);
                }
                else
                {
                    for (i = 0; i < nAtoms; i++)
                    {
                        for (j = 0; j < nValues; j++)
                        {
                            to[i*nValues+j] = (real)((double *)from)[i*nValues+j] * fact;
                        }
                    }
                }
            }
            else
            {
                for (i = 0; i < nAtoms; i++)
                {
                    for (j = 0; j < nValues; j++)
                    {
                        to[i*nValues+j] = (real)((double *)from)[i*nValues+j] * fact;
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
#endif

gmx_bool read_next_tng_frame(tng_trajectory_t        input,
                             t_trxframe             *fr)
{
#ifdef GMX_USE_TNG
    gmx_bool            bOK = TRUE;
    tng_function_status stat;
    int64_t             numberOfAtoms = -1, frameNumber = -1;
    int64_t             requestedIds[5] = {TNG_TRAJ_BOX_SHAPE, TNG_TRAJ_POSITIONS,
                                           TNG_TRAJ_VELOCITIES, TNG_TRAJ_FORCES,
                                           TNG_GMX_LAMBDA};
    int64_t             nBlocks, blockId, *blockIds = NULL;
    char                datatype      = -1, codecId;
    void               *values        = NULL;
    double              frameTime     = -1.0;
    int                 size, blockDependency;
    float               prec;


    fr->bStep     = FALSE;
    fr->bTime     = FALSE;
    fr->bLambda   = FALSE;
    fr->bAtoms    = FALSE;
    fr->bPrec     = FALSE;
    fr->bX        = FALSE;
    fr->bV        = FALSE;
    fr->bF        = FALSE;
    fr->bBox      = FALSE;

    stat = tng_num_particles_get(input, &numberOfAtoms);
    if (stat != TNG_SUCCESS)
    {
        gmx_file("Cannot determine number of atoms from TNG file.");
    }
    fr->natoms = numberOfAtoms;

    if (!get_tng_data_block_types_of_next_frame(input,
                                                fr->step,
                                                5, requestedIds,
                                                &frameNumber,
                                                &nBlocks,
                                                &blockIds))
    {
        return FALSE;
    }

    if(nBlocks == 0)
    {
        return FALSE;
    }

    for(int64_t i = 0; i < nBlocks; i++)
    {
        blockId = blockIds[i];
        tng_data_block_dependency_get(input, blockId, &blockDependency);
        if(blockDependency & TNG_PARTICLE_DEPENDENT)
        {
            stat = tng_util_particle_data_next_frame_read(input,
                                                          blockId,
                                                          &values,
                                                          &datatype,
                                                          &frameNumber,
                                                          &frameTime);
        }
        else
        {
            stat = tng_util_non_particle_data_next_frame_read(input,
                                                              blockId,
                                                              &values,
                                                              &datatype,
                                                              &frameNumber,
                                                              &frameTime);
        }
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
            continue;
        }
        switch(blockId)
        {
            case TNG_TRAJ_BOX_SHAPE:
                switch (datatype)
                {
                    case TNG_INT_DATA:
                        size = sizeof(int64_t);
                        break;
                    case TNG_FLOAT_DATA:
                        size = sizeof(float);
                        break;
                    case TNG_DOUBLE_DATA:
                        size = sizeof(double);
                        break;
                    default:
                        size = 0; /* Just to make the compiler happy. */
                        gmx_incons("Illegal datatype of box shape values!");
                }
                for (int i = 0; i < DIM; i++)
                {
                    tng_convert_array_to_real_array((char *)(values) + size * i * DIM,
                                                    (real *) fr->box[i],
                                                    getDistanceScaleFactor(input),
                                                    1,
                                                    DIM,
                                                    datatype);
                }
                fr->bBox = TRUE;
                break;
            case TNG_TRAJ_POSITIONS:
                srenew(fr->x, fr->natoms);
                tng_convert_array_to_real_array(values,
                                                (real *) fr->x,
                                                getDistanceScaleFactor(input),
                                                fr->natoms,
                                                DIM,
                                                datatype);
                fr->bX = TRUE;
                tng_util_frame_current_compression_get(input, blockId, &codecId, &prec);
                /* This must be updated if/when more lossy compression methods are added */
                if(codecId == TNG_TNG_COMPRESSION)
                {
                    fr->prec = prec;
                    fr->bPrec = TRUE;
                }
                break;
            case TNG_TRAJ_VELOCITIES:
                srenew(fr->v, fr->natoms);
                tng_convert_array_to_real_array(values,
                                                (real *) fr->v,
                                                getDistanceScaleFactor(input),
                                                fr->natoms,
                                                DIM,
                                                datatype);
                fr->bV = TRUE;
                break;
            case TNG_TRAJ_FORCES:
                srenew(fr->f, fr->natoms);
                tng_convert_array_to_real_array(values,
                                                (real *) fr->f,
                                                getDistanceScaleFactor(input),
                                                fr->natoms,
                                                DIM,
                                                datatype);
                fr->bF = TRUE;
                break;
            case TNG_GMX_LAMBDA:
                switch (datatype)
                {
                    case TNG_FLOAT_DATA:
                        fr->lambda = (*(float *)values);
                        break;
                    case TNG_DOUBLE_DATA:
                        fr->lambda = (*(double *)values);
                        break;
                    default:
                        gmx_incons("Illegal datatype lambda value!");
                }
                fr->bLambda = TRUE;
                break;
            default:
                gmx_incons("Illegal block type!");
        }
    /* values does not have to be freed before reading next frame. It will
     * be reallocated if it is not NULL. */
    }

    fr->step = (int) frameNumber;
    fr->bStep = TRUE;
    // Convert the time to ps
    fr->time = frameTime / PICO;
    fr->bTime = TRUE;

    /* values must be freed before leaving this function */
    sfree(values);

    return bOK;
#else
    GMX_UNUSED_VALUE(input);
    GMX_UNUSED_VALUE(fr);
    return FALSE;
#endif
}

void tng_set_compression_precision(tng_trajectory_t tng,
                                   real prec)
{
#ifdef GMX_USE_TNG
    tng_compression_precision_set(tng, prec);
#else
    GMX_UNUSED_VALUE(tng);
    GMX_UNUSED_VALUE(prec);
#endif
}

void print_tng_molecule_system(tng_trajectory_t input,
                               FILE *stream)
{
#ifdef GMX_USE_TNG
    int64_t nMolecules, nChains, nResidues, nAtoms, *molCntList;
    tng_molecule_t molecule;
    tng_chain_t chain;
    tng_residue_t residue;
    tng_atom_t atom;
    char str[256], varNAtoms;

    tng_num_molecule_types_get(input, &nMolecules);
    tng_molecule_cnt_list_get(input, &molCntList);
    tng_num_particles_variable_get(input, &varNAtoms);

    for (int64_t i = 0; i < nMolecules; i++)
    {
        tng_molecule_of_index_get(input, i, &molecule);
        tng_molecule_name_get(input, molecule, str, 256);
        if(varNAtoms == TNG_CONSTANT_N_ATOMS)
        {
            if((int)molCntList[i] == 0)
            {
                continue;
            }
            fprintf(stream, "Molecule: %s, count: %d\n", str, (int)molCntList[i]);
        }
        else
        {
            fprintf(stream, "Molecule: %s\n", str);
        }
        tng_molecule_num_chains_get(input, molecule, &nChains);
        if(nChains > 0)
        {
            for(int64_t j = 0; j < nChains; j++)
            {
                tng_molecule_chain_of_index_get(input, molecule, j, &chain);
                tng_chain_name_get(input, chain, str, 256);
                fprintf(stream, "\tChain: %s\n", str);
                tng_chain_num_residues_get(input, chain, &nResidues);
                for(int64_t k = 0; k < nResidues; k++)
                {
                    tng_chain_residue_of_index_get(input, chain, k, &residue);
                    tng_residue_name_get(input, residue, str, 256);
                    fprintf(stream, "\t\tResidue: %s\n", str);
                    tng_residue_num_atoms_get(input, residue, &nAtoms);
                    for(int64_t l = 0; l < nAtoms; l++)
                    {
                        tng_residue_atom_of_index_get(input, residue, l, &atom);
                        tng_atom_name_get(input, atom, str, 256);
                        fprintf(stream, "\t\t\tAtom: %s", str);
                        tng_atom_type_get(input, atom, str, 256);
                        fprintf(stream, " (%s)\n", str);
                    }
                }
            }
        }
        /* It is possible to have a molecule without chains, in which case
         * residues in the molecule can be iterated through without going
         * through chains. */
        else
        {
            tng_molecule_num_residues_get(input, molecule, &nResidues);
            if(nResidues > 0)
            {
                for(int64_t k = 0; k < nResidues; k++)
                {
                    tng_molecule_residue_of_index_get(input, molecule, k, &residue);
                    tng_residue_name_get(input, residue, str, 256);
                    fprintf(stream, "\t\tResidue: %s\n", str);
                    tng_residue_num_atoms_get(input, residue, &nAtoms);
                    for(int64_t l = 0; l < nAtoms; l++)
                    {
                        tng_residue_atom_of_index_get(input, residue, l, &atom);
                        tng_atom_name_get(input, atom, str, 256);
                        fprintf(stream, "\t\t\tAtom: %s", str);
                        tng_atom_type_get(input, atom, str, 256);
                        fprintf(stream, " (%s)\n", str);
                    }
                }
            }
            else
            {
                tng_molecule_num_atoms_get(input, molecule, &nAtoms);
                for(int64_t l = 0; l < nAtoms; l++)
                {
                    tng_molecule_atom_of_index_get(input, molecule, l, &atom);
                    tng_atom_name_get(input, atom, str, 256);
                    fprintf(stream, "\t\t\tAtom: %s", str);
                    tng_atom_type_get(input, atom, str, 256);
                    fprintf(stream, " (%s)\n", str);
                }
            }
        }
    }
#else
    GMX_UNUSED_VALUE(input);
    GMX_UNUSED_VALUE(stream);
#endif
}

gmx_bool get_tng_data_block_types_of_next_frame(tng_trajectory_t input,
                                                int frame,
                                                int nRequestedIds,
                                                int64_t *requestedIds,
                                                int64_t *nextFrame,
                                                int64_t *nBlocks,
                                                int64_t **blockIds)
{
#ifdef GMX_USE_TNG
    tng_function_status stat;

    stat = tng_util_trajectory_next_frame_present_data_blocks_find(input, frame,
                                                                   nRequestedIds, requestedIds,
                                                                   nextFrame,
                                                                   nBlocks, blockIds);

    if(stat == TNG_CRITICAL)
    {
        gmx_file("Cannot read TNG file. Cannot find data blocks of next frame.");
    }
    else if(stat == TNG_FAILURE)
    {
        return FALSE;
    }
    return TRUE;
#else
    GMX_UNUSED_VALUE(input);
    GMX_UNUSED_VALUE(frame);
    GMX_UNUSED_VALUE(nRequestedIds);
    GMX_UNUSED_VALUE(requestedIds);
    GMX_UNUSED_VALUE(nextFrame);
    GMX_UNUSED_VALUE(nBlocks);
    GMX_UNUSED_VALUE(blockIds);
    return FALSE;
#endif
}

gmx_bool get_tng_data_next_frame_of_block_type(tng_trajectory_t input,
                                               int64_t blockId,
                                               real **values,
                                               int64_t *frameNumber,
                                               double *frameTime,
                                               int64_t *nValuesPerFrame,
                                               int64_t *nAtoms,
                                               real *prec,
                                               char *name,
                                               int maxLen,
                                               gmx_bool *bOK)
{
#ifdef GMX_USE_TNG
    tng_function_status stat;
    char datatype = -1, codecId;
    int blockDependency;
    void *data = 0;
    float localPrec;

    stat = tng_data_block_name_get(input, blockId, name, maxLen);
    if(stat != TNG_SUCCESS)
    {
        gmx_file("Cannot read next frame of TNG file");
    }
    stat = tng_data_block_dependency_get(input, blockId, &blockDependency);
    if(stat != TNG_SUCCESS)
    {
        gmx_file("Cannot read next frame of TNG file");
    }
    if(blockDependency & TNG_PARTICLE_DEPENDENT)
    {
        tng_num_particles_get(input, nAtoms);
        stat = tng_util_particle_data_next_frame_read(input,
                                                      blockId,
                                                      &data,
                                                      &datatype,
                                                      frameNumber,
                                                      frameTime);
    }
    else
    {
        *nAtoms = 1; /* There are not actually any atoms, but it is used for
                        allocating memory */
        stat = tng_util_non_particle_data_next_frame_read(input,
                                                          blockId,
                                                          &data,
                                                          &datatype,
                                                          frameNumber,
                                                          frameTime);
    }
    if(stat == TNG_CRITICAL)
    {
        gmx_file("Cannot read next frame of TNG file");
    }
    if(stat == TNG_FAILURE)
    {
        *bOK = TRUE;
        return FALSE;
    }

    stat = tng_data_block_num_values_per_frame_get(input, blockId, nValuesPerFrame);
    if(stat != TNG_SUCCESS)
    {
        gmx_file("Cannot read next frame of TNG file");
    }
    snew(*values, sizeof(real) * *nValuesPerFrame * *nAtoms);
    tng_convert_array_to_real_array(data,
                                    *values,
                                    getDistanceScaleFactor(input),
                                    *nAtoms,
                                    *nValuesPerFrame,
                                    datatype);

    tng_util_frame_current_compression_get(input, blockId, &codecId, &localPrec);

    /* This must be updated if/when more lossy compression methods are added */
    if(codecId != TNG_TNG_COMPRESSION)
    {
        *prec = -1.0;
    }
    else
    {
        *prec = localPrec;
    }

    *bOK = TRUE;
    return TRUE;
#else
    GMX_UNUSED_VALUE(input);
    GMX_UNUSED_VALUE(blockId);
    GMX_UNUSED_VALUE(values);
    GMX_UNUSED_VALUE(frameNumber);
    GMX_UNUSED_VALUE(frameTime);
    GMX_UNUSED_VALUE(nValuesPerFrame);
    GMX_UNUSED_VALUE(nAtoms);
    GMX_UNUSED_VALUE(prec);
    GMX_UNUSED_VALUE(name);
    GMX_UNUSED_VALUE(maxLen);
    GMX_UNUSED_VALUE(bOK);
    return FALSE;
#endif
}