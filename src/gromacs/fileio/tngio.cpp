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

#include "tng_io.h"
#include <algorithm>
#include "gromacs/legacyheaders/gmx_fatal.h"
#include "gromacs/legacyheaders/physics.h"
#include "gromacs/utility/gmxassert.h"
#include "gmxfio.h"

const int defaultFramesPerFrameSet = 100;

/* Define pointers to specific writing functions depending on if we
 * write float or double data */
typedef tng_function_status (*tng_set_writing_interval_func_pointer)(tng_trajectory_t,
                                                                     const int64_t,
                                                                     const int64_t,
                                                                     const int64_t,
                                                                     const char*,
                                                                     const char,
                                                                     const char);
#ifdef GMX_DOUBLE
tng_set_writing_interval_func_pointer tng_set_writing_interval = tng_util_generic_write_interval_double_set;
typedef tng_function_status (*tng_write_data_func_pointer)(tng_trajectory_t,
                                                           const int64_t,
                                                           const double,
                                                           const double*,
                                                           const int64_t,
                                                           const int64_t,
                                                           const char  *,
                                                           const char,
                                                           const char);
tng_write_data_func_pointer           tng_write_data = tng_util_generic_with_time_double_write;
#else
tng_set_writing_interval_func_pointer tng_set_writing_interval = tng_util_generic_write_interval_set;
typedef tng_function_status (*tng_write_data_func_pointer)(tng_trajectory_t,
                                                           const int64_t,
                                                           const double,
                                                           const float*,
                                                           const int64_t,
                                                           const int64_t,
                                                           const char *,
                                                           const char,
                                                           const char);
tng_write_data_func_pointer tng_write_data = tng_util_generic_with_time_write;
#endif

static const char *modeToVerb(char mode)
{
    switch(mode)
    {
    case 'r':
        return "reading";
        break;
    case 'w':
        return "writing";
        break;
    case 'a':
        return "appending";
        break;
    default:
        gmx_fatal(FARGS, "Invalid file opening mode %c", mode);
        return "";
    }
}

/*! \libinternal \brief Compute least common denominator of n1 and n2
 *
 * This duplicates a routine in domdec_setup.c, but that's probably
 * best until we make a general maths utility module. */
static int lcd(int n1, int n2)
{
    int d, i;

    d = 1;
    for (i = 2; (i <= n1 && i <= n2); i++)
    {
        if (n1 % i == 0 && n2 % i == 0)
        {
            d = i;
        }
    }

    return d;
}

/*! \libinternal \brief Compute least common denominator of n1 and n2
 * if they are positive.
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

// TODO Should file operations be mutexed?

void tng_open(const char       gmx_unused *filename,
              char             gmx_unused  mode,
              tng_trajectory_t gmx_unused *tng)
{
#ifdef GMX_USE_TNG
    /* First check whether we have to make a backup,
     * only for writing, not for read or append.
     */
    if (mode == 'w')
    {
#ifndef GMX_FAHCORE
        /* only make backups for normal gromacs */
        make_backup(filename);
#endif
    }

    /* tng must not be pointing at already allocated memory.
     * Memory will be allocated by tng_util_trajectory_open() and must
     * later on be freed by tng_util_trajectory_close(). */
    if (TNG_SUCCESS != tng_util_trajectory_open(filename, mode, tng))
    {
        /* TNG does return more than one degree of error, but there is
           no use case for GROMACS handling the non-fatal errors
           gracefully. */
        gmx_fatal(FARGS,
                  "%s while opening %s for %s",
                  gmx_strerror("file"),
                  filename,
                  modeToVerb(mode));
    }

    /* TODO: Here we should add a time stamp (of time of creation)
     * if mode is 'w'. We should also add info about user and program. */
#endif
}

void tng_close(tng_trajectory_t gmx_unused *tng)
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

static void addTngMoleculeFromTopology(tng_trajectory_t tng,
                                       const char *moleculeName,
                                       const t_atoms *atoms,
                                       int64_t numMolecules,
                                       tng_molecule_t *tng_mol)
{
    // TODO can these TNG calls fail and can/should we do anything
    // about it?
    tng_molecule_add(tng, moleculeName, tng_mol);

    /* FIXME: The TNG atoms should contain mass and atomB info (for free
     * energy calculations) */
    for (int at_it = 0; at_it < atoms->nr; at_it++)
    {
        const t_atom *at = &atoms->atom[at_it];
        /* FIXME: Currently the TNG API can only add atoms belonging to a
         * residue and chain. */
        if (atoms->nres > 0)
        {
            const t_resinfo *resin = &atoms->resinfo[at->resind];
            char chainName[2] = {resin->chainid, 0};
            tng_chain_t tng_chain = NULL;
            tng_residue_t tng_res  = NULL;
            tng_atom_t tng_atom    = NULL;

            if (tng_molecule_chain_find (tng, *tng_mol, chainName,
                                         (int64_t)-1, &tng_chain) !=
                TNG_SUCCESS)
            {
                tng_molecule_chain_add (tng, *tng_mol, chainName,
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
    tng_molecule_cnt_set(tng, *tng_mol, numMolecules);
}

static void tng_add_mtop(tng_trajectory_t gmx_unused  tng,
                         const gmx_mtop_t gmx_unused *mtop)
{
#ifdef GMX_USE_TNG
    int                  i, j;
    const t_ilist       *ilist;
    tng_bond_t           tng_bond;

    for (int mol_it = 0; mol_it < mtop->nmolblock; mol_it++)
    {
        tng_molecule_t tng_mol = NULL;
        const gmx_moltype_t *mol_type =
            &mtop->moltype[mtop->molblock[mol_it].type];

        /* Add a molecule to the TNG trajectory with the same name as the
         * current molecule. */
        addTngMoleculeFromTopology(tng,
                                   *(mol_type->name),
                                   &mol_type->atoms,
                                   mtop->molblock[mol_it].nmol,
                                   &tng_mol);

        /* Bonds have to be deduced from interactions (constraints etc). Different
         * interactions have different sets of parameters. */
        /* Constraints are specified using two atoms */
        for (i = F_CONSTR; i <= F_CONSTRNC; i++)
        {
            ilist = &mol_type->ilist[i];
            if (ilist)
            {
                j = 1;
                while (j < ilist->nr)
                {
                    tng_molecule_bond_add(tng, tng_mol, ilist->iatoms[j], ilist->iatoms[j+1], &tng_bond);
                    j += 3;
                }
            }
        }
        /* Settle is described using three atoms */
        ilist = &mol_type->ilist[F_SETTLE];
        if (ilist)
        {
            j = 1;
            while (j < ilist->nr)
            {
                tng_molecule_bond_add(tng, tng_mol, ilist->iatoms[j], ilist->iatoms[j+1], &tng_bond);
                tng_molecule_bond_add(tng, tng_mol, ilist->iatoms[j], ilist->iatoms[j+2], &tng_bond);
                j += 4;
            }
        }
        /* FIXME: Are there more interactions that should be used to generate bonds? Virtual sites? */
    }
#endif
}

// TODO this is never called - use g_tool -s information here?
static void tng_add_topology(tng_trajectory_t gmx_unused  tng,
                             const t_topology gmx_unused *top)
{
#ifdef GMX_USE_TNG
    for (int mol_it = 0; mol_it < top->mols.nr; mol_it++)
    {
        tng_molecule_t tng_mol = NULL;

        /* Add a molecule to the TNG trajectory with the same name as the
         * current molecule. */
        addTngMoleculeFromTopology(tng,
                                   "unknown molecule",
                                   &top->atoms,
                                   top->mols.nr,
                                   &tng_mol);
    }
#endif
}

void tng_set_writing_intervals(tng_trajectory_t gmx_unused  tng,
                               const gmx_bool gmx_unused    bUseLossyCompression,
                               const t_inputrec gmx_unused *ir)
{
#ifdef GMX_USE_TNG
    int  xout, vout, fout;
    int low_denom = -1, lowest = -1;
    char compression;

    tng_set_frames_per_frame_set(tng, bUseLossyCompression, ir);

    if (bUseLossyCompression)
    {
        xout        = ir->nstxtcout;
        vout        = 0;
        fout        = 0;
        compression = TNG_TNG_COMPRESSION;
    }
    else
    {
        xout        = ir->nstxout;
        vout        = ir->nstvout;
        fout        = ir->nstfout;
        compression = TNG_GZIP_COMPRESSION;
    }
    if (xout)
    {
        tng_set_writing_interval(tng, xout, 3, TNG_TRAJ_POSITIONS,
                                 "POSITIONS", TNG_PARTICLE_BLOCK_DATA,
                                 compression);
        /* The design of TNG makes it awkward to try to write a box
         * with multiple periodicities, which might be co-prime. Since
         * the use cases for the box with a frame consisting only of
         * velocities seem low, for now we associate box writing with
         * position writing. */
        tng_set_writing_interval(tng, xout, 9, TNG_TRAJ_BOX_SHAPE,
                                 "BOX SHAPE", TNG_NON_PARTICLE_BLOCK_DATA,
                                 TNG_GZIP_COMPRESSION);
        /* TODO: if/when we write energies to TNG also, reconsider how
         * and when box information is written, because GROMACS
         * behaviour pre-5.0 was to write the box with every
         * trajectory frame and every energy frame, and probably
         * people depend on this. */

        low_denom = lcd_if_positive(low_denom, xout);
        if (lowest < 0 || xout < lowest)
        {
            lowest = xout;
        }
    }
    if (vout)
    {
        tng_set_writing_interval(tng, ir->nstvout, 3, TNG_TRAJ_VELOCITIES,
                                 "VELOCITIES", TNG_PARTICLE_BLOCK_DATA,
                                 compression);
        low_denom = lcd_if_positive(low_denom, vout);
        if (lowest < 0 || vout < lowest)
        {
            lowest = vout;
        }
    }
    if (fout)
    {
        tng_set_writing_interval(tng, ir->nstfout, 3, TNG_TRAJ_FORCES,
                                 "FORCES", TNG_PARTICLE_BLOCK_DATA,
                                 TNG_GZIP_COMPRESSION);
        low_denom = lcd_if_positive(low_denom, fout);
        if (lowest < 0 || fout < lowest)
        {
            lowest = fout;
        }
    }
    if (low_denom > 0)
    {
        /* Lambdas written at an interval of the lowest common denominator
         * of other output */
        tng_set_writing_interval(tng, low_denom, 1, TNG_GMX_LAMBDA,
                                 "LAMBDAS", TNG_NON_PARTICLE_BLOCK_DATA,
                                 TNG_GZIP_COMPRESSION);

        if(low_denom < lowest / 10)
        {
            gmx_warning("The lowest common denominator of trajectory output is "
                        "every %d step(s), whereas the shortest output interval "
                        "is every %d steps.", low_denom, lowest);
        }
    }
#endif
}

void tng_prepare_md_writing(tng_trajectory_t gmx_unused  tng,
                            const gmx_mtop_t gmx_unused *mtop,
                            const t_inputrec gmx_unused *ir)
{
#ifdef GMX_USE_TNG
    tng_add_mtop(tng, mtop);
    tng_set_writing_intervals(tng, FALSE, ir);
    tng_time_per_frame_set(tng, ir->delta_t * PICO);
#endif
}

void tng_prepare_low_prec_writing(tng_trajectory_t gmx_unused  tng,
                                  const gmx_mtop_t gmx_unused *mtop,
                                  const t_inputrec gmx_unused *ir)
{
#ifdef GMX_USE_TNG
    tng_add_mtop(tng, mtop);
    tng_add_selection_groups(tng, mtop);
    tng_set_writing_intervals(tng, TRUE, ir);
    tng_time_per_frame_set(tng, ir->delta_t * PICO);
    tng_compression_precision_set(tng, 1.0/ir->xtcprec);
#endif
}

void tng_add_selection_groups(tng_trajectory_t gmx_unused  tng,
                              const gmx_mtop_t gmx_unused *mtop)
{
#ifdef GMX_USE_TNG
    const gmx_moltype_t *mol_type;
    const t_atoms       *atoms;
    const t_atom        *at;
    const t_resinfo     *resin;
    const t_ilist       *ilist;
    int                  n_atoms = 0, i = 0, j, mol_it, at_it, name_index;
    int                  atom_offset = 0;
    tng_molecule_t       mol, iter_mol;
    tng_chain_t          chain;
    tng_residue_t        res;
    tng_atom_t           atom;
    tng_bond_t           tng_bond;
    int64_t              n_mols;
    char                *group_name;

    /* TODO: Set the molecule name according to the selection name */
    name_index = *mtop->groups.grps[egcXTC].nm_ind;
    group_name = *mtop->groups.grpname[name_index];

    tng_molecule_alloc(tng, &mol);
    tng_molecule_name_set(tng, mol, group_name);
    tng_molecule_chain_add(tng, mol, "", &chain);
    for (mol_it = 0; mol_it < mtop->nmoltype; mol_it++)
    {
        mol_type = &mtop->moltype[mtop->molblock[mol_it].type];

        atoms = &mol_type->atoms;

        for (j = 0; j < mtop->molblock[mol_it].nmol; j++)
        {
            bool bAtomsAdded = FALSE;
            for (at_it = 0; at_it < atoms->nr; at_it++, i++)
            {
                char *res_name;
                int   res_id;

                if (ggrpnr(&mtop->groups, egcXTC, i) != 0)
                {
                    continue;
                }
                at = &atoms->atom[at_it];
                if (atoms->nres > 0)
                {
                    resin = &atoms->resinfo[at->resind];
                    /* FIXME: When TNG supports both residue index and residue
                     * number the latter should be used. */
                    res_name = *resin->name;
                    res_id   = at->resind + 1;
                }
                else
                {
                    res_name = (char *)"";
                    res_id   = 0;
                }
                if (tng_chain_residue_find(tng, chain, res_name, res_id, &res)
                    != TNG_SUCCESS)
                {
                    tng_chain_residue_w_id_add(tng, chain, res_name, res_id, &res);
                }
                tng_residue_atom_w_id_add(tng, res, *(atoms->atomname[at_it]),
                                          *(atoms->atomtype[at_it]),
                                          atom_offset + at_it, &atom);
                n_atoms++;
                bAtomsAdded = TRUE;
            }
            /* Add bonds. */
            if(bAtomsAdded)
            {
                for (int k = F_CONSTR; k <= F_CONSTRNC; k++)
                {
                    ilist = &mol_type->ilist[k];
                    if (ilist)
                    {
                        int l = 1;
                        while (l < ilist->nr)
                        {
                            int atom1, atom2;
                            atom1 = ilist->iatoms[l] + atom_offset;
                            atom2 = ilist->iatoms[l+1] + atom_offset;
                            if(ggrpnr(&mtop->groups, egcXTC, atom1) == 0 &&
                               ggrpnr(&mtop->groups, egcXTC, atom2) == 0)
                            {
                                tng_molecule_bond_add(tng, mol, ilist->iatoms[l],
                                                      ilist->iatoms[l+1], &tng_bond);
                            }
                            l += 3;
                        }
                    }
                }
                /* Settle is described using three atoms */
                ilist = &mol_type->ilist[F_SETTLE];
                if (ilist)
                {
                    int l = 1;
                    while (l < ilist->nr)
                    {
                        int atom1, atom2, atom3;
                        atom1 = ilist->iatoms[l] + atom_offset;
                        atom2 = ilist->iatoms[l+1] + atom_offset;
                        atom3 = ilist->iatoms[l+2] + atom_offset;
                        if(ggrpnr(&mtop->groups, egcXTC, atom1) == 0)
                        {
                            if(ggrpnr(&mtop->groups, egcXTC, atom2) == 0)
                            {
                                tng_molecule_bond_add(tng, mol, atom1,
                                                      atom2, &tng_bond);
                            }
                            if(ggrpnr(&mtop->groups, egcXTC, atom3) == 0)
                            {
                                tng_molecule_bond_add(tng, mol, atom1,
                                                      atom3, &tng_bond);
                            }
                        }
                        l += 4;
                    }
                }
            }
            atom_offset += atoms->nr;
        }
    }
    if (n_atoms != i)
    {
        tng_molecule_existing_add(tng, &mol);
        tng_molecule_cnt_set(tng, mol, 1);
        tng_num_molecule_types_get(tng, &n_mols);
        for (int64_t k = 0; k < n_mols; k++)
        {
            tng_molecule_of_index_get(tng, k, &iter_mol);
            if (iter_mol == mol)
            {
                continue;
            }
            tng_molecule_cnt_set(tng, iter_mol, 0);
        }
    }
    else
    {
        tng_molecule_free(tng, &mol);
        return;
    }
#endif
}

void tng_set_frames_per_frame_set(tng_trajectory_t gmx_unused  tng,
                                  const gmx_bool gmx_unused    bUseLossyCompression,
                                  const t_inputrec gmx_unused *ir)
{
#ifdef GMX_USE_TNG
    int     lcd;
    int64_t numFramesPerFrameSet;

    /* Set the number of frames per frame set to contain at least
     * defaultFramesPerFrameSet of the lowest common denominator of
     * the writing frequency of positions and velocities. */
    /* FIXME later: consider nstenergy also? */
    if (bUseLossyCompression)
    {
        lcd = ir->nstxtcout;
    }
    else
    {
        lcd = lcd_if_positive(ir->nstxout, ir->nstvout);
        lcd = lcd_if_positive(lcd, ir->nstfout);
    }
    if (0 >= lcd)
    {
        return;
    }

    numFramesPerFrameSet = std::min((int) ir->nsteps,
                                    lcd * defaultFramesPerFrameSet);

    tng_num_frames_per_frame_set_set(tng, numFramesPerFrameSet);
#endif
}

void fwrite_tng(tng_trajectory_t tng,
                const gmx_bool   bUseLossyCompression,
                int              step,
                real             elapsedPicoSeconds,
                real             lambda,
                const rvec      *box,
                int              natoms,
                const rvec      *x,
                const rvec      *v,
                const rvec      *f)
{
#ifdef GMX_USE_TNG
    double              elapsedSeconds = elapsedPicoSeconds * PICO;
    int64_t             tng_n_particles;
    char                compression;

    GMX_ASSERT(box, "Need a non-NULL box");

    if (!tng)
    {
        return;
    }

    tng_num_particles_get(tng, &tng_n_particles);
    if (natoms != (int)tng_n_particles)
    {
        tng_implicit_num_particles_set(tng, natoms);
    }

    if (bUseLossyCompression)
    {
        compression = TNG_TNG_COMPRESSION;
    }
    else
    {
        compression = TNG_GZIP_COMPRESSION;
    }

    /* TODO: what does TNG do if you pass it a NULL pointer? Do we
     * need to even check?
     * - Comment: The functions return TNG_FAILURE without performing any writing if the data pointer is NULL.
     * It would be possible to ignore TNG_FAILURE and just stop if receiving TNG_CRITICAL and thereby not require
     * any checks here. But it can be good to catch failures as well.
     */
    /* The writing is done using tng_write_data, which writes float or double
     * depending on the GROMACS compilation. */
    if (x)
    {
        if(tng_write_data(tng, step, elapsedSeconds,
                          reinterpret_cast<const real *>(x),
                          3, TNG_TRAJ_POSITIONS, "POSITIONS",
                          TNG_PARTICLE_BLOCK_DATA,
                          compression) != TNG_SUCCESS)
        {
            gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
            return;
        }
        /* TNG-MF1 compression only compresses positions and velocities. Use lossless
         * compression for box shape regardless of output mode */
        if(tng_write_data(tng, step, elapsedSeconds,
                          reinterpret_cast<const real *>(box),
                          9, TNG_TRAJ_BOX_SHAPE, "BOX SHAPE",
                          TNG_NON_PARTICLE_BLOCK_DATA,
                          TNG_GZIP_COMPRESSION) != TNG_SUCCESS)
        {
            gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
            return;
        }
    }

    if (v)
    {
        if(tng_write_data(tng, step, elapsedSeconds,
                          reinterpret_cast<const real *>(v),
                          3, TNG_TRAJ_VELOCITIES, "VELOCITIES",
                          TNG_PARTICLE_BLOCK_DATA,
                          compression) != TNG_SUCCESS)
        {
            gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
            return;
        }
    }
    
    if (f)
    {
        /* TNG-MF1 compression only compresses positions and velocities. Use lossless
         * compression for forces regardless of output mode */
        if(tng_write_data(tng, step, elapsedSeconds,
                          reinterpret_cast<const real *>(f),
                          3, TNG_TRAJ_FORCES, "FORCES",
                          TNG_PARTICLE_BLOCK_DATA,
                          TNG_GZIP_COMPRESSION) != TNG_SUCCESS)
        {
            gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
            return;
        }
    }

    /* TNG-MF1 compression only compresses positions and velocities. Use lossless
     * compression for lambdas regardless of output mode */
    if(tng_write_data(tng, step, elapsedSeconds,
                      reinterpret_cast<const real *>(&lambda),
                      1, TNG_GMX_LAMBDA, "LAMBDAS",
                      TNG_NON_PARTICLE_BLOCK_DATA,
                      TNG_GZIP_COMPRESSION) != TNG_SUCCESS)
    {
        gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
        return;
    }
#endif
}
