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

#include "gromacs/legacyheaders/domdec.h"
#include "gromacs/legacyheaders/gmx_fatal.h"
#include "gromacs/utility/gmxassert.h"
#ifdef GMX_USE_TNG
#include "tng_io.h"
#else
typedef tng_trajectory_t int;
#endif

/* TODO: I'm assuming TNG copes with being passed a NULL
   tng_trajectory_t */

void tng_open(const char *filename,
              char mode,
              tng_trajectory_t *tng)
{
#ifdef GMX_USE_TNG
  tng_util_trajectory_open(filename, mode, tng);
#endif
}

void tng_close(tng_trajectory_t *tng)
{
#ifdef GMX_USE_TNG
    tng_util_trajectory_close(tng);
#endif
}

void tng_add_top(tng_trajectory_t tng,
                 const gmx_mtop_t *mtop)
{
#ifdef GMX_USE_TNG
    int mol_it, at_it;
    char chain_name[2];
    const gmx_moltype_t *mol_type;
    const t_atoms *atoms;
    const t_atom *at;
    const t_resinfo *resin;
    tng_molecule_t tng_mol;
    tng_chain_t tng_chain;
    tng_residue_t tng_res;
    tng_atom_t tng_atom;

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
                resin = &atoms->resinfo[at->resind];
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

void tng_set_writing_frequency(tng_trajectory_t tng,
                               const t_inputrec *ir)
{
#ifdef GMX_USE_TNG
    int lowcd;
    int64_t n_frames_per_frame_set;

    /* Set the number of frames per frame set to contain at least 100
     * of the lowest common denominator of the writing frequency of
     * positions and velocities.
     * FIXME: Forces should also be taken into account.
     * FIXME later: energies also?
     */
    if(ir->nstxout && ir->nstvout)
    {
        lowcd = lcd(ir->nstxout, ir->nstvout);
    }
    else if(ir->nstxout)
    {
        lowcd = ir->nstxout;
    }
    else
    {
        lowcd = ir->nstvout;
    }
    if (lowcd <= 0)
    {
        lowcd = 1;
    }
    n_frames_per_frame_set = lowcd * 100;
    if (n_frames_per_frame_set > ir->nsteps)
    {
        n_frames_per_frame_set = ir->nsteps;
    }

    tng_num_frames_per_frame_set_set(tng, n_frames_per_frame_set);

    if (ir->nstxout)
    {
        tng_util_pos_write_frequency_set(tng, ir->nstxout);
    }
    if (ir->nstvout)
    {
        tng_util_vel_write_frequency_set(tng, ir->nstvout);
    }
    if (ir->nstfout)
    {
        tng_util_force_write_frequency_set(tng, ir->nstfout);
    }
    if (ir->nstxout || ir->nstvout || ir->nstfout)
    {
        tng_util_box_shape_write_frequency_set(tng, lowcd);
        /* FIXME: Currently lcd(3,5) means we write a box every 1
           step, but the trajectory writing routine will not (and
           should not) get called every lcd(3,5) steps. We should be
           able to write a box (and perhaps other junk, like time, or
           lambda) every time we write output, because that functions
           as a label for the output. And this kind of data might be
           missing for a given frame, too. We should also be able to
           just write the box at a known frequency, too. I hope the
           spec allowed for all that ;-) In real-world use, lcd() will
           tend to be more sensible, but specs and code need either to
           cater for all use cases, or declare their constraints. */
        /* TODO if this use of lcd eventually goes away, revert
           changes to the lcd() function */
    }
#endif
}

/* FIXME: Write a reading function. The need to do allocation will be
   messy (and the other GROMACS trajectory reading code is not a good
   example here), but hopefully reading the TNG headers will help set
   this up. */

void fwrite_tng(tng_trajectory_t tng, int step, real t, real lambda,
                const rvec *box, int natoms, const rvec *x, const rvec *v, const rvec *f)
{
#ifdef GMX_USE_TNG
    bool bOK = true;
    tng_function_status stat;

    GMX_ASSERT(box, "Need a non-NULL box");
    stat = tng_util_box_shape_write(tng, step, reinterpret_cast<const float *>(box));
    bOK = bOK && stat == TNG_SUCCESS;

    if (x)
    {
        stat = tng_util_pos_write(tng, step, reinterpret_cast<const float *>(x));
        bOK = bOK && stat == TNG_SUCCESS;
    }
    if (v)
    {
        stat = tng_util_vel_write(tng, step, reinterpret_cast<const float *>(v));
        bOK = bOK && stat == TNG_SUCCESS;
    }
    if (f)
    {
        stat = tng_util_force_write(tng, step, reinterpret_cast<const float *>(f));
        bOK = bOK && stat == TNG_SUCCESS;
    }

    /* TODO: do something with time and lambda */

    if (!bOK)
    {
        gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
    }
#endif
}

