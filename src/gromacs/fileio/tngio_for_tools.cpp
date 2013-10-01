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
#else
typedef struct tng_trajectory *tng_trajectory_t;
#endif

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
                         int                      natoms)
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
        tng_implicit_num_particles_set(*output, natoms);

        tng_num_frames_per_frame_set_set(*output, 1);
    }
#else
    GMX_UNUSED_VALUE(filename);
    GMX_UNUSED_VALUE(mode);
    GMX_UNUSED_VALUE(input);
    GMX_UNUSED_VALUE(output);
    GMX_UNUSED_VALUE(natoms);
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
        double time_per_frame = frame->time * PICO / frame->step;
        tng_time_per_frame_set(output, time_per_frame);
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
#endif

gmx_bool read_next_tng_frame(tng_trajectory_t        input,
                             t_trxframe             *fr)
{
#ifdef GMX_USE_TNG
    gmx_bool            bOK = TRUE;
    tng_function_status stat;
    int64_t             numberOfAtoms = -1, frameNumber = -1;
    char                datatype      = -1;
    void               *values        = NULL;
    double              frameTime     = -1.0;
    int                 size;

    stat = tng_num_particles_get(input, &numberOfAtoms);
    if (stat != TNG_SUCCESS)
    {
        gmx_file("Cannot determine number of atoms from TNG file.");
    }
    fr->natoms = numberOfAtoms;

    stat = tng_util_particle_data_next_frame_read(input,
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
                                    getDistanceScaleFactor(input),
                                    fr->natoms,
                                    DIM,
                                    datatype);

    /* values does not have to be freed before reading next frame. It will
     * be reallocated if it is not NULL. */

    // TODO how do we fill lambda, etc. (if we need to? check GROMACS
    // client code)

    stat = tng_util_non_particle_data_next_frame_read(input,
                                                      TNG_TRAJ_BOX_SHAPE,
                                                      &values,
                                                      &datatype,
                                                      &frameNumber,
                                                      &frameTime);
    if (stat != TNG_SUCCESS)
    {
        gmx_file("Cannot read box shape from TNG file.");
    }
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

    /* values must be freed before leaving this function */
    free(values);

    return bOK;
#else
    GMX_UNUSED_VALUE(input);
    GMX_UNUSED_VALUE(fr);
    return FALSE;
#endif
}

void print_tng_molecule_system(tng_trajectory_t input,
                               FILE *stream)
{
#ifdef GMX_USE_TNG
    int64_t n_molecules, n_chains, n_residues, n_atoms, *mol_cnt_list;
    tng_molecule_t molecule;
    tng_chain_t chain;
    tng_residue_t residue;
    tng_atom_t atom;
    char str[256], var_n_atoms;

    tng_num_molecule_types_get(input, &n_molecules);
    tng_molecule_cnt_list_get(input, &mol_cnt_list);
    tng_num_particles_variable_get(input, &var_n_atoms);

    for (int64_t i = 0; i < n_molecules; i++)
    {
        tng_molecule_of_index_get(input, i, &molecule);
        tng_molecule_name_get(input, molecule, str, 256);
        if(var_n_atoms == TNG_CONSTANT_N_ATOMS)
        {
            if((int)mol_cnt_list[i] == 0)
            {
                continue;
            }
            fprintf(stream, "Molecule: %s, count: %d\n", str, (int)mol_cnt_list[i]);
        }
        else
        {
            fprintf(stream, "Molecule: %s\n", str);
        }
        tng_molecule_num_chains_get(input, molecule, &n_chains);
        if(n_chains > 0)
        {
            for(int64_t j = 0; j < n_chains; j++)
            {
                tng_molecule_chain_of_index_get(input, molecule, j, &chain);
                tng_chain_name_get(input, chain, str, 256);
                fprintf(stream, "\tChain: %s\n", str);
                tng_chain_num_residues_get(input, chain, &n_residues);
                for(int64_t k = 0; k < n_residues; k++)
                {
                    tng_chain_residue_of_index_get(input, chain, k, &residue);
                    tng_residue_name_get(input, residue, str, 256);
                    fprintf(stream, "\t\tResidue: %s\n", str);
                    tng_residue_num_atoms_get(input, residue, &n_atoms);
                    for(int64_t l = 0; l < n_atoms; l++)
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
        else
        {
            tng_molecule_num_residues_get(input, molecule, &n_residues);
            if(n_residues > 0)
            {
                for(int64_t k = 0; k < n_residues; k++)
                {
                    tng_molecule_residue_of_index_get(input, molecule, k, &residue);
                    tng_residue_name_get(input, residue, str, 256);
                    fprintf(stream, "\t\tResidue: %s\n", str);
                    tng_residue_num_atoms_get(input, residue, &n_atoms);
                    for(int64_t l = 0; l < n_atoms; l++)
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
                tng_molecule_num_atoms_get(input, molecule, &n_atoms);
                for(int64_t l = 0; l < n_atoms; l++)
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
                                                int n_requested_ids,
                                                int64_t *requested_ids,
                                                int64_t *next_frame,
                                                int64_t *n_blocks,
                                                int64_t **block_ids)
{
#ifdef GMX_USE_TNG
    tng_function_status stat;

    stat = tng_util_trajectory_next_frame_present_data_blocks_find(input, frame,
                                                                   n_requested_ids, requested_ids,
                                                                   next_frame,
                                                                   n_blocks, block_ids);

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
    GMX_UNUSED_VALUE(fr);
    GMX_UNUSED_VALUE(n_blocks);
    GMX_UNUSED_VALUE(block_ids);
    return FALSE;
#endif
}

gmx_bool get_tng_data_next_frame_of_block_type(tng_trajectory_t input,
                                               int64_t blockId,
                                               real **values,
                                               int64_t *frameNumber,
                                               double *frameTime,
                                               int64_t *n_values_per_frame,
                                               int64_t *n_atoms,
                                               real *prec,
                                               char *name,
                                               int max_len,
                                               gmx_bool *bOK)
{
    tng_function_status stat;
    char datatype = -1, codec_id;
    int block_dependency;
    void *data = 0;

    stat = tng_data_block_name_get(input, blockId, name, max_len);
    if(stat != TNG_SUCCESS)
    {
        gmx_file("Cannot read next frame of TNG file");
    }
    stat = tng_data_block_dependency_get(input, blockId, &block_dependency);
    if(stat != TNG_SUCCESS)
    {
        gmx_file("Cannot read next frame of TNG file");
    }
    if(block_dependency & TNG_PARTICLE_DEPENDENT)
    {
        tng_num_particles_get(input, n_atoms);
        stat = tng_util_particle_data_next_frame_read(input,
                                                      blockId,
                                                      &data,
                                                      &datatype,
                                                      frameNumber,
                                                      frameTime);
    }
    else
    {
        *n_atoms = 1; /* There are not actually any atoms, but it is used for
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

    stat = tng_data_block_num_values_per_frame_get(input, blockId, n_values_per_frame);
    if(stat != TNG_SUCCESS)
    {
        gmx_file("Cannot read next frame of TNG file");
    }
    snew(*values, sizeof(real) * *n_values_per_frame * *n_atoms);
    tng_convert_array_to_real_array(data,
                                    *values,
                                    getDistanceScaleFactor(input),
                                    *n_atoms,
                                    *n_values_per_frame,
                                    datatype);
    
    tng_util_frame_current_compression_get(input, blockId, &codec_id, prec);
    
    /* This must be updated if/when more lossy compression methods are added */
    if(codec_id != TNG_TNG_COMPRESSION)
    {
        *prec = -1.0;
    }

    *bOK = TRUE;
    return TRUE;
}