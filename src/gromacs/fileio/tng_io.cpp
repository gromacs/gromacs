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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "tng_io.h"

#ifdef GMX_USE_TNG
void init_tng_top(tng_trajectory_t tng, const gmx_mtop_t *mtop)
/* This function converts the current topology to TNG molecular data */
{
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
}

void init_tng_writing_frequency(tng_trajectory_t tng, t_inputrec *ir)
{
    int lowcd;
    int64_t n_frames_per_frame_set;

    /* Set the number of frames per frame set to contain at least 100
     * of the lowest common denominator of the writing frequency of
     * positions and velocities.
     * FIXME: Forces should also be taken into account.
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

    if(ir->nstxout)
    {
        tng_util_pos_write_frequency_set(tng, ir->nstxout);
    }
    if(ir->nstvout)
    {
        tng_util_vel_write_frequency_set(tng, ir->nstvout);
    }
    if(ir->nstfout)
    {
        tng_util_force_write_frequency_set(tng, ir->nstfout);
    }
    if(ir->nstxout || ir->nstvout || ir->nstfout)
    {
        tng_util_box_shape_write_frequency_set(tng, lowcd);
    }
}

#endif

#ifdef GMX_USE_TNG
/* TNG: STUB */
/*static gmx_bool do_htng(tng_trajectory_t tng, gmx_bool bRead, t_trnheader *sh,
                        rvec *box, rvec *x, rvec *v, rvec *f)
{
    gmx_bool bOK;
    tng_function_status stat;
    int64_t n_particles;

    bOK = TRUE;

    *//* FIXME: The box shape has to be written with a stride length of the
     * lowest common denominator of the other written data. *//*

    *//* FIXME: Only writing floats *//*

    if (sh->box_size != 0)
    {
        stat = tng_util_box_shape_write(tng, sh->step, (float *)box);
        bOK = bOK && stat == TNG_SUCCESS;
    }
    if (sh->vir_size != 0)
    {
        *//* FIXME: Not implemented yet *//*
    }
    if (sh->pres_size != 0)
    {
        *//* FIXME: Not implemented yet *//*
    }
    if (sh->x_size != 0)
    {
        stat = tng_util_pos_write(tng, sh->step, (float *)x);
        bOK = bOK && stat == TNG_SUCCESS;
    }
    if (sh->v_size != 0)
    {
        stat = tng_util_vel_write(tng, sh->step, (float *)v);
        bOK = bOK && stat == TNG_SUCCESS;
    }
    if (sh->f_size != 0)
    {
        stat = tng_util_force_write(tng, sh->step, (float *)f);
        bOK = bOK && stat == TNG_SUCCESS;
    }

    return bOK;
}*/
#endif

#ifdef GMX_USE_TNG
/* TNG: STUB */
/*static gmx_bool do_tng(tng_trajectory_t tng, gmx_bool bRead, int *step, real *t, real *lambda,
                       rvec *box, int *natoms, rvec *x, rvec *v, rvec *f)
{
    t_trnheader *sh;
    gmx_bool     bOK;

    snew(sh, 1);
    if (!bRead)
    {
        sh->box_size = (box) ? sizeof(matrix) : 0;
        sh->x_size   = ((x) ? (*natoms*sizeof(x[0])) : 0);
        sh->v_size   = ((v) ? (*natoms*sizeof(v[0])) : 0);
        sh->f_size   = ((f) ? (*natoms*sizeof(f[0])) : 0);
        sh->natoms   = *natoms;
        sh->step     = *step;
        sh->nre      = 0;
        sh->t        = *t;
        sh->lambda   = *lambda;
    }
    *//* FIXME: Reading does not work yet. *//*
    if (bRead)
    {
        *natoms = sh->natoms;
        *step   = sh->step;
        *t      = sh->t;
        *lambda = sh->lambda;
        if (sh->ir_size)
        {
            gmx_file("inputrec in trn file");
        }
        if (sh->e_size)
        {
            gmx_file("energies in trn file");
        }
        if (sh->top_size)
        {
            gmx_file("topology in trn file");
        }
        if (sh->sym_size)
        {
            gmx_file("symbol table in trn file");
        }
    }
    bOK = do_htng(tng, bRead, sh, box, x, v, f);

    sfree(sh);

    return bOK;
}*/
#endif

#ifdef GMX_USE_TNG
/* TNG: STUB */
/*void fwrite_tng(tng_trajectory_t tng, int step, real t, real lambda,
                rvec *box, int natoms, rvec *x, rvec *v, rvec *f)
{
    if (do_tng(tng, FALSE, &step, &t, &lambda, box, &natoms, x, v, f) == FALSE)
    {
        gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
    }
}*/
#endif

#ifdef GMX_USE_TNG
/* TNG: STUB */
/*static gmx_bool do_htng(tng_trajectory_t tng, gmx_bool bRead, t_trnheader *sh,
                        rvec *box, rvec *x, rvec *v, rvec *f)
{
    gmx_bool bOK;
    tng_function_status stat;
    int64_t n_particles;

    bOK = TRUE;

    *//* FIXME: The box shape has to be written with a stride length of the
     * lowest common denominator of the other written data. *//*

    *//* FIXME: Only writing floats *//*

    if (sh->box_size != 0)
    {
        stat = tng_util_box_shape_write(tng, sh->step, (float *)box);
        bOK = bOK && stat == TNG_SUCCESS;
    }
    if (sh->vir_size != 0)
    {
        *//* FIXME: Not implemented yet *//*
    }
    if (sh->pres_size != 0)
    {
        *//* FIXME: Not implemented yet *//*
    }
    if (sh->x_size != 0)
    {
        stat = tng_util_pos_write(tng, sh->step, (float *)x);
        bOK = bOK && stat == TNG_SUCCESS;
    }
    if (sh->v_size != 0)
    {
        stat = tng_util_vel_write(tng, sh->step, (float *)v);
        bOK = bOK && stat == TNG_SUCCESS;
    }
    if (sh->f_size != 0)
    {
        stat = tng_util_force_write(tng, sh->step, (float *)f);
        bOK = bOK && stat == TNG_SUCCESS;
    }

    return bOK;
}*/
#endif

#ifdef GMX_USE_TNG
/* TNG: STUB */
/*static gmx_bool do_tng(tng_trajectory_t tng, gmx_bool bRead, int *step, real *t, real *lambda,
                       rvec *box, int *natoms, rvec *x, rvec *v, rvec *f)
{
    t_trnheader *sh;
    gmx_bool     bOK;

    snew(sh, 1);
    if (!bRead)
    {
        sh->box_size = (box) ? sizeof(matrix) : 0;
        sh->x_size   = ((x) ? (*natoms*sizeof(x[0])) : 0);
        sh->v_size   = ((v) ? (*natoms*sizeof(v[0])) : 0);
        sh->f_size   = ((f) ? (*natoms*sizeof(f[0])) : 0);
        sh->natoms   = *natoms;
        sh->step     = *step;
        sh->nre      = 0;
        sh->t        = *t;
        sh->lambda   = *lambda;
    }
    *//* FIXME: Reading does not work yet. *//*
    if (bRead)
    {
        *natoms = sh->natoms;
        *step   = sh->step;
        *t      = sh->t;
        *lambda = sh->lambda;
        if (sh->ir_size)
        {
            gmx_file("inputrec in trn file");
        }
        if (sh->e_size)
        {
            gmx_file("energies in trn file");
        }
        if (sh->top_size)
        {
            gmx_file("topology in trn file");
        }
        if (sh->sym_size)
        {
            gmx_file("symbol table in trn file");
        }
    }
    bOK = do_htng(tng, bRead, sh, box, x, v, f);

    sfree(sh);

    return bOK;
}*/
#endif

#ifdef GMX_USE_TNG
/* TNG: STUB */
/*void fwrite_tng(tng_trajectory_t tng, int step, real t, real lambda,
                rvec *box, int natoms, rvec *x, rvec *v, rvec *f)
{
    if (do_tng(tng, FALSE, &step, &t, &lambda, box, &natoms, x, v, f) == FALSE)
    {
        gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
    }
}*/
#endif


