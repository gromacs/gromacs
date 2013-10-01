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
#include "tngio.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include "gromacs/legacyheaders/domdec.h"
#include "gromacs/legacyheaders/gmx_fatal.h"
#include "gromacs/utility/gmxassert.h"
#ifdef GMX_USE_TNG
#include "tng_io.h"
#else
typedef tng_trajectory_t int;
#endif

const int defaultFramesPerFrameSet = 100;

void tng_open(const char       *filename,
              char              mode,
              tng_trajectory_t *tng)
{
#ifdef GMX_USE_TNG
    /* tng must not be pointing at already allocated memory.
     * Memory will be allocated by tng_util_trajectory_open() and must
     * later on be freed by tng_util_trajectory_close(). */
    tng_util_trajectory_open(filename, mode, tng);
#endif
}

void tng_close(tng_trajectory_t *tng)
{
    /* TODO: What behaviour is defined by TNG if clients try to close
       before open? Or close after an error condition arises? The
       existence of a segfault when we called without the GROMACS
       check of MASTER(cr) suggests TNG should be more defensive
       here. */
#ifdef GMX_USE_TNG
    tng_util_trajectory_close(tng);
#endif
}

void tng_add_top(tng_trajectory_t  tng,
                 const gmx_mtop_t *mtop)
{
#ifdef GMX_USE_TNG
    int                  mol_it, at_it;
    char                 chain_name[2];
    const gmx_moltype_t *mol_type;
    const t_atoms       *atoms;
    const t_atom        *at;
    const t_resinfo     *resin;
    tng_molecule_t       tng_mol;
    tng_chain_t          tng_chain;
    tng_residue_t        tng_res;
    tng_atom_t           tng_atom;

    for (mol_it = 0; mol_it < mtop->nmoltype; mol_it++)
    {
        mol_type = &mtop->moltype[mol_it];

        /* Add a molecule to the TNG trajectory with the same name as the
         * current molecule. */
        tng_molecule_add(tng, *mol_type->name, &tng_mol);

        atoms = &mol_type->atoms;

        /* FIXME: The TNG atoms should contain mass and atomB info (for free
         * energy calculations */
        for (at_it = 0; at_it < atoms->nr; at_it++)
        {
            at = &atoms->atom[at_it];
            /* FIXME: Currently the TNG API can only add atoms belonging to a
             * residue and chain. */
            if (atoms->nres > 0)
            {
                resin         = &atoms->resinfo[at->resind];
                chain_name[0] = resin->chainid;
                chain_name[1] = 0;
                if (tng_molecule_chain_find (tng, tng_mol, chain_name,
                                             (int64_t)-1, &tng_chain) !=
                    TNG_SUCCESS)
                {
                    tng_molecule_chain_add (tng, tng_mol, chain_name,
                                            &tng_chain);
                }

                /* FIXME: When TNG supports both residue index and residue
                 * number the latter should be used. */
                if (tng_chain_residue_find(tng, tng_chain, *resin->name,
                                           at->resind + 1, &tng_res)
                    != TNG_SUCCESS)
                {
                    tng_chain_residue_add(tng, tng_chain, *resin->name, &tng_res);
                }
                tng_residue_atom_add(tng, tng_res, *(atoms->atomname[at_it]), *(atoms->atomtype[at_it]), &tng_atom);
            }
        }
        tng_molecule_cnt_set(tng, tng_mol, mtop->molblock[mol_it].nmol);
    }
#endif
}

/* Compute least common denominator of n1 and n2 if they are positive.
 *
 * If only one of n1 and n2 is positive, then return it.
 * If neither n1 or n2 is positive, then return -1. */
static int
lcd_if_positive(int n1, int n2)
{
    if (0 >= n1)
    {
        return (0 >= n2) ? -1 : n2;
    }
    if (0 >= n2)
    {
        return n1;
    }

    /* We have a non-trivial least common denominator to compute. */
    return lcd(n1, n2);
}

void tng_set_writing_frequency(tng_trajectory_t  tng,
                               const t_inputrec *ir)
{
#ifdef GMX_USE_TNG
    int     lcd;
    int64_t numFramesPerFrameSet;

    /* Set the number of frames per frame set to contain at least
     * defaultFramesPerFrameSet of the lowest common denominator of
     * the writing frequency of positions and velocities. */
    /* FIXME later: consider nstenergy also? */
    lcd = lcd_if_positive(ir->nstxout, ir->nstvout);
    lcd = lcd_if_positive(lcd, ir->nstfout);
    if (0 >= lcd)
    {
        return;
    }

    numFramesPerFrameSet = std::min((int) ir->nsteps,
                                    lcd * defaultFramesPerFrameSet);

    tng_num_frames_per_frame_set_set(tng, numFramesPerFrameSet);

    if (ir->nstxout)
    {
        tng_util_pos_write_frequency_set(tng, ir->nstxout);
        /* The design of TNG makes it awkward to try to write a box
         * with multiple periodicities, which might be co-prime. Since
         * the use cases for the box with a frame consisting only of
         * velocities seem low, for now we associate box writing with
         * position writing. */
        tng_util_box_shape_write_frequency_set(tng, ir->nstxout);
        /* TODO: if/when we write energies to TNG also, reconsider how
         * and when box information is written, because GROMACS
         * behaviour pre-5.0 was to write the box with every
         * trajectory frame and every energy frame, and probably
         * people depend on this. */
    }
    /* TODO: we are not passing a frequency to the
       tng_util_*_write_frequency_set() functions, we are passing an
       interval between events. Should this be fixed in the TNG
       API? */
    if (ir->nstvout)
    {
        tng_util_vel_write_frequency_set(tng, ir->nstvout);
    }
    if (ir->nstfout)
    {
        tng_util_force_write_frequency_set(tng, ir->nstfout);
    }
#endif
}

void tng_set_time_step(tng_trajectory_t  tng,
                               const t_inputrec *ir)
{
#ifdef GMX_USE_TNG
    /* In TNG the time is stored in seconds */
    tng_time_per_frame_set(tng, ir->delta_t*0.000000000001);
#endif
}
/* FIXME: Write a reading function. The need to do allocation will be
 * messy (and the other GROMACS trajectory reading code is not a good
 * example here), but hopefully reading the TNG headers will help set
 * this up.
 * - Comment: The TNG functions handle the memory allocation. But it
 * needs to be explicitly freed afterwards.
 */

void fwrite_tng(tng_trajectory_t tng,
                int step, real t, real lambda,
                const rvec *box, int natoms,
                const rvec *x, const rvec *v, const rvec *f)
{
#ifdef GMX_USE_TNG
    bool                bOK = true;
    tng_function_status stat;

    GMX_ASSERT(box, "Need a non-NULL box");
    stat = tng_util_box_shape_with_time_write(tng, step, (double)t*0.000000000001, reinterpret_cast<const float *>(box));
    bOK  = bOK && stat == TNG_SUCCESS;

    /* TODO: what does TNG do if you pass it a NULL pointer? Do we
     * need to even check? */
    if (x)
    {
        stat = tng_util_pos_with_time_write(tng, step, (double)t*0.000000000001, reinterpret_cast<const float *>(x));
        bOK  = bOK && stat == TNG_SUCCESS;
    }
    /* TODO: should we be giving up now if !bOK? */
    if (v)
    {
        stat = tng_util_vel_with_time_write(tng, step, (double)t*0.000000000001, reinterpret_cast<const float *>(v));
        bOK  = bOK && stat == TNG_SUCCESS;
    }
    if (f)
    {
        stat = tng_util_force_with_time_write(tng, step, (double)t*0.000000000001, reinterpret_cast<const float *>(f));
        bOK  = bOK && stat == TNG_SUCCESS;
    }

    /* TODO: do something with time and lambda - they don't get
     * written with every frame, but when should they get written? */

    if (!bOK)
    {
        gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
    }
#endif
}
