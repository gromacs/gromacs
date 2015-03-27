/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, 2015 by the GROMACS development team, led by
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

#include "tngio.h"

#include "config.h"

#ifdef GMX_USE_TNG
#include "tng/tng_io.h"
#endif

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/types/ifunc.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/sysinfo.h"

static const char *modeToVerb(char mode)
{
    const char *p;
    switch (mode)
    {
        case 'r':
            p = "reading";
            break;
        case 'w':
            p = "writing";
            break;
        case 'a':
            p = "appending";
            break;
        default:
            gmx_fatal(FARGS, "Invalid file opening mode %c", mode);
            p = "";
            break;
    }
    return p;
}

void gmx_tng_open(const char       *filename,
                  char              mode,
                  tng_trajectory_t *tng)
{
#ifdef GMX_USE_TNG
    /* First check whether we have to make a backup,
     * only for writing, not for read or append.
     */
    if (mode == 'w')
    {
        make_backup(filename);
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

    if (mode == 'w' || mode == 'a')
    {
        char hostname[256];
        gmx_gethostname(hostname, 256);
        if (mode == 'w')
        {
            tng_first_computer_name_set(*tng, hostname);
        }
        else
        {
            tng_last_computer_name_set(*tng, hostname);
        }

        char        programInfo[256];
        const char *precisionString = "";
#ifdef GMX_DOUBLE
        precisionString = " (double precision)";
#endif
        sprintf(programInfo, "%.100s, %.128s%.24s",
                gmx::getProgramContext().displayName(),
                GromacsVersion(), precisionString);
        if (mode == 'w')
        {
            tng_first_program_name_set(*tng, programInfo);
        }
        else
        {
            tng_last_program_name_set(*tng, programInfo);
        }

        char username[256];
        if (!gmx_getusername(username, 256))
        {
            if (mode == 'w')
            {
                tng_first_user_name_set(*tng, username);
            }
            else
            {
                tng_last_user_name_set(*tng, username);
                tng_file_headers_write(*tng, TNG_USE_HASH);
            }
        }
    }
#else
    gmx_file("GROMACS was compiled without TNG support, cannot handle this file type");
    GMX_UNUSED_VALUE(filename);
    GMX_UNUSED_VALUE(mode);
    GMX_UNUSED_VALUE(tng);
#endif
}

void gmx_tng_close(tng_trajectory_t *tng)
{
    /* We have to check that tng is set because
     * tng_util_trajectory_close wants to return a NULL in it, and
     * gives a fatal error if it is NULL. */
#ifdef GMX_USE_TNG
    if (tng)
    {
        tng_util_trajectory_close(tng);
    }
#else
    GMX_UNUSED_VALUE(tng);
#endif
}

#ifdef GMX_USE_TNG
static void addTngMoleculeFromTopology(tng_trajectory_t     tng,
                                       const char          *moleculeName,
                                       const t_atoms       *atoms,
                                       gmx_int64_t          numMolecules,
                                       tng_molecule_t      *tngMol)
{
    tng_chain_t      tngChain = NULL;
    tng_residue_t    tngRes   = NULL;

    if (tng_molecule_add(tng, moleculeName, tngMol) != TNG_SUCCESS)
    {
        gmx_file("Cannot add molecule to TNG molecular system.");
    }

    /* FIXME: The TNG atoms should contain mass and atomB info (for free
     * energy calculations), i.e. in when it's available in TNG (2.0). */
    for (int atomIndex = 0; atomIndex < atoms->nr; atomIndex++)
    {
        const t_atom *at = &atoms->atom[atomIndex];
        /* FIXME: Currently the TNG API can only add atoms belonging to a
         * residue and chain. Wait for TNG 2.0*/
        if (atoms->nres > 0)
        {
            const t_resinfo *resInfo        = &atoms->resinfo[at->resind];
            char             chainName[2]   = {resInfo->chainid, 0};
            tng_atom_t       tngAtom        = NULL;
            t_atom          *prevAtom;

            if (atomIndex > 0)
            {
                prevAtom = &atoms->atom[atomIndex - 1];
            }
            else
            {
                prevAtom = 0;
            }

            /* If this is the first atom or if the residue changed add the
             * residue to the TNG molecular system. */
            if (!prevAtom || resInfo != &atoms->resinfo[prevAtom->resind])
            {
                /* If this is the first atom or if the chain changed add
                 * the chain to the TNG molecular system. */
                if (!prevAtom || resInfo->chainid !=
                    atoms->resinfo[prevAtom->resind].chainid)
                {
                    tng_molecule_chain_add(tng, *tngMol, chainName,
                                           &tngChain);
                }
                /* FIXME: When TNG supports both residue index and residue
                 * number the latter should be used. Wait for TNG 2.0*/
                tng_chain_residue_add(tng, tngChain, *resInfo->name, &tngRes);
            }
            tng_residue_atom_add(tng, tngRes, *(atoms->atomname[atomIndex]), *(atoms->atomtype[atomIndex]), &tngAtom);
        }
    }
    tng_molecule_cnt_set(tng, *tngMol, numMolecules);
}

void gmx_tng_add_mtop(tng_trajectory_t  tng,
                      const gmx_mtop_t *mtop)
{
    int                  i, j;
    const t_ilist       *ilist;
    tng_bond_t           tngBond;

    if (!mtop)
    {
        /* No topology information available to add. */
        return;
    }

    for (int molIndex = 0; molIndex < mtop->nmolblock; molIndex++)
    {
        tng_molecule_t       tngMol  = NULL;
        const gmx_moltype_t *molType =
            &mtop->moltype[mtop->molblock[molIndex].type];

        /* Add a molecule to the TNG trajectory with the same name as the
         * current molecule. */
        addTngMoleculeFromTopology(tng,
                                   *(molType->name),
                                   &molType->atoms,
                                   mtop->molblock[molIndex].nmol,
                                   &tngMol);

        /* Bonds have to be deduced from interactions (constraints etc). Different
         * interactions have different sets of parameters. */
        /* Constraints are specified using two atoms */
        for (i = 0; i < F_NRE; i++)
        {
            if (IS_CHEMBOND(i))
            {
                ilist = &molType->ilist[i];
                if (ilist)
                {
                    j = 1;
                    while (j < ilist->nr)
                    {
                        tng_molecule_bond_add(tng, tngMol, ilist->iatoms[j], ilist->iatoms[j+1], &tngBond);
                        j += 3;
                    }
                }
            }
        }
        /* Settle is described using three atoms */
        ilist = &molType->ilist[F_SETTLE];
        if (ilist)
        {
            j = 1;
            while (j < ilist->nr)
            {
                tng_molecule_bond_add(tng, tngMol, ilist->iatoms[j], ilist->iatoms[j+1], &tngBond);
                tng_molecule_bond_add(tng, tngMol, ilist->iatoms[j], ilist->iatoms[j+2], &tngBond);
                j += 4;
            }
        }
    }
}

/*! \libinternal \brief Compute greatest common divisor of n1 and n2
 * if they are positive.
 *
 * If only one of n1 and n2 is positive, then return it.
 * If neither n1 or n2 is positive, then return -1. */
static int
greatest_common_divisor_if_positive(int n1, int n2)
{
    if (0 >= n1)
    {
        return (0 >= n2) ? -1 : n2;
    }
    if (0 >= n2)
    {
        return n1;
    }

    /* We have a non-trivial greatest common divisor to compute. */
    return gmx_greatest_common_divisor(n1, n2);
}

/* By default try to write 100 frames (of actual output) in each frame set.
 * This number is the number of outputs of the most frequently written data
 * type per frame set.
 * TODO for 5.1: Verify that 100 frames per frame set is efficient for most
 * setups regarding compression efficiency and compression time. Make this
 * a hidden command-line option? */
const int defaultFramesPerFrameSet = 100;

/*! \libinternal \brief  Set the number of frames per frame
 * set according to output intervals.
 * The default is that 100 frames are written of the data
 * that is written most often. */
static void tng_set_frames_per_frame_set(tng_trajectory_t  tng,
                                         const gmx_bool    bUseLossyCompression,
                                         const t_inputrec *ir)
{
    int     gcd = -1;

    /* Set the number of frames per frame set to contain at least
     * defaultFramesPerFrameSet of the lowest common denominator of
     * the writing interval of positions and velocities. */
    /* FIXME after 5.0: consider nstenergy also? */
    if (bUseLossyCompression)
    {
        gcd = ir->nstxout_compressed;
    }
    else
    {
        gcd = greatest_common_divisor_if_positive(ir->nstxout, ir->nstvout);
        gcd = greatest_common_divisor_if_positive(gcd, ir->nstfout);
    }
    if (0 >= gcd)
    {
        return;
    }

    tng_num_frames_per_frame_set_set(tng, gcd * defaultFramesPerFrameSet);
}

/*! \libinternal \brief Set the data-writing intervals, and number of
 * frames per frame set */
static void set_writing_intervals(tng_trajectory_t  tng,
                                  const gmx_bool    bUseLossyCompression,
                                  const t_inputrec *ir)
{
    /* Define pointers to specific writing functions depending on if we
     * write float or double data */
    typedef tng_function_status (*set_writing_interval_func_pointer)(tng_trajectory_t,
                                                                     const gmx_int64_t,
                                                                     const gmx_int64_t,
                                                                     const gmx_int64_t,
                                                                     const char*,
                                                                     const char,
                                                                     const char);
#ifdef GMX_DOUBLE
    set_writing_interval_func_pointer set_writing_interval = tng_util_generic_write_interval_double_set;
#else
    set_writing_interval_func_pointer set_writing_interval = tng_util_generic_write_interval_set;
#endif
    int  xout, vout, fout;
//     int  gcd = -1, lowest = -1;
    char compression;

    tng_set_frames_per_frame_set(tng, bUseLossyCompression, ir);

    if (bUseLossyCompression)
    {
        xout        = ir->nstxout_compressed;
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
        set_writing_interval(tng, xout, 3, TNG_TRAJ_POSITIONS,
                             "POSITIONS", TNG_PARTICLE_BLOCK_DATA,
                             compression);
        /* The design of TNG makes it awkward to try to write a box
         * with multiple periodicities, which might be co-prime. Since
         * the use cases for the box with a frame consisting only of
         * velocities seem low, for now we associate box writing with
         * position writing. */
        set_writing_interval(tng, xout, 9, TNG_TRAJ_BOX_SHAPE,
                             "BOX SHAPE", TNG_NON_PARTICLE_BLOCK_DATA,
                             TNG_GZIP_COMPRESSION);
        /* TODO: if/when we write energies to TNG also, reconsider how
         * and when box information is written, because GROMACS
         * behaviour pre-5.0 was to write the box with every
         * trajectory frame and every energy frame, and probably
         * people depend on this. */

        /* TODO: If we need to write lambda values at steps when
         * positions (or other data) are not also being written, then
         * code in mdoutf.c will need to match however that is
         * organized here. */
        set_writing_interval(tng, xout, 1, TNG_GMX_LAMBDA,
                             "LAMBDAS", TNG_NON_PARTICLE_BLOCK_DATA,
                             TNG_GZIP_COMPRESSION);

        /* FIXME: gcd and lowest currently not used. */
//         gcd = greatest_common_divisor_if_positive(gcd, xout);
//         if (lowest < 0 || xout < lowest)
//         {
//             lowest = xout;
//         }
    }
    if (vout)
    {
        set_writing_interval(tng, ir->nstvout, 3, TNG_TRAJ_VELOCITIES,
                             "VELOCITIES", TNG_PARTICLE_BLOCK_DATA,
                             compression);

        /* FIXME: gcd and lowest currently not used. */
//         gcd = greatest_common_divisor_if_positive(gcd, vout);
//         if (lowest < 0 || vout < lowest)
//         {
//             lowest = vout;
//         }
    }
    if (fout)
    {
        set_writing_interval(tng, ir->nstfout, 3, TNG_TRAJ_FORCES,
                             "FORCES", TNG_PARTICLE_BLOCK_DATA,
                             TNG_GZIP_COMPRESSION);

        /* FIXME: gcd and lowest currently not used. */
//         gcd = greatest_common_divisor_if_positive(gcd, fout);
//         if (lowest < 0 || fout < lowest)
//         {
//             lowest = fout;
//         }
    }
    /* FIXME: See above. gcd interval for lambdas is disabled. */
//     if (gcd > 0)
//     {
//         /* Lambdas written at an interval of the lowest common denominator
//          * of other output */
//         set_writing_interval(tng, gcd, 1, TNG_GMX_LAMBDA,
//                                  "LAMBDAS", TNG_NON_PARTICLE_BLOCK_DATA,
//                                  TNG_GZIP_COMPRESSION);
//
//         if (gcd < lowest / 10)
//         {
//             gmx_warning("The lowest common denominator of trajectory output is "
//                         "every %d step(s), whereas the shortest output interval "
//                         "is every %d steps.", gcd, lowest);
//         }
//     }
}
#endif

void gmx_tng_prepare_md_writing(tng_trajectory_t  tng,
                                const gmx_mtop_t *mtop,
                                const t_inputrec *ir)
{
#ifdef GMX_USE_TNG
    gmx_tng_add_mtop(tng, mtop);
    set_writing_intervals(tng, FALSE, ir);
    tng_time_per_frame_set(tng, ir->delta_t * PICO);
#else
    GMX_UNUSED_VALUE(tng);
    GMX_UNUSED_VALUE(mtop);
    GMX_UNUSED_VALUE(ir);
#endif
}

#ifdef GMX_USE_TNG
/* Check if all atoms in the molecule system are selected
 * by a selection group of type specified by gtype. */
static gmx_bool all_atoms_selected(const gmx_mtop_t *mtop,
                                   const int         gtype)
{
    const gmx_moltype_t     *molType;
    const t_atoms           *atoms;

    /* Iterate through all atoms in the molecule system and
     * check if they belong to a selection group of the
     * requested type. */
    for (int molIndex = 0, i = 0; molIndex < mtop->nmoltype; molIndex++)
    {
        molType = &mtop->moltype[mtop->molblock[molIndex].type];

        atoms = &molType->atoms;

        for (int j = 0; j < mtop->molblock[molIndex].nmol; j++)
        {
            for (int atomIndex = 0; atomIndex < atoms->nr; atomIndex++, i++)
            {
                if (ggrpnr(&mtop->groups, gtype, i) != 0)
                {
                    return FALSE;
                }
            }
        }
    }
    return TRUE;
}

/* Create TNG molecules which will represent each of the selection
 * groups for writing. But do that only if there is actually a
 * specified selection group and it is not the whole system.
 * TODO: Currently the only selection that is taken into account
 * is egcCompressedX, but other selections should be added when
 * e.g. writing energies is implemented.
 */
static void add_selection_groups(tng_trajectory_t  tng,
                                 const gmx_mtop_t *mtop)
{
    const gmx_moltype_t     *molType;
    const t_atoms           *atoms;
    const t_atom            *at;
    const t_resinfo         *resInfo;
    const t_ilist           *ilist;
    int                      nameIndex;
    int                      atom_offset = 0;
    tng_molecule_t           mol, iterMol;
    tng_chain_t              chain;
    tng_residue_t            res;
    tng_atom_t               atom;
    tng_bond_t               tngBond;
    gmx_int64_t              nMols;
    char                    *groupName;

    /* TODO: When the TNG molecules block is more flexible TNG selection
     * groups should not need all atoms specified. It should be possible
     * just to specify what molecules are selected (and/or which atoms in
     * the molecule) and how many (if applicable). */

    /* If no atoms are selected we do not need to create a
     * TNG selection group molecule. */
    if (mtop->groups.ngrpnr[egcCompressedX] == 0)
    {
        return;
    }

    /* If all atoms are selected we do not have to create a selection
     * group molecule in the TNG molecule system. */
    if (all_atoms_selected(mtop, egcCompressedX))
    {
        return;
    }

    /* The name of the TNG molecule containing the selection group is the
     * same as the name of the selection group. */
    nameIndex = *mtop->groups.grps[egcCompressedX].nm_ind;
    groupName = *mtop->groups.grpname[nameIndex];

    tng_molecule_alloc(tng, &mol);
    tng_molecule_name_set(tng, mol, groupName);
    tng_molecule_chain_add(tng, mol, "", &chain);
    for (int molIndex = 0, i = 0; molIndex < mtop->nmoltype; molIndex++)
    {
        molType = &mtop->moltype[mtop->molblock[molIndex].type];

        atoms = &molType->atoms;

        for (int j = 0; j < mtop->molblock[molIndex].nmol; j++)
        {
            bool bAtomsAdded = FALSE;
            for (int atomIndex = 0; atomIndex < atoms->nr; atomIndex++, i++)
            {
                char *res_name;
                int   res_id;

                if (ggrpnr(&mtop->groups, egcCompressedX, i) != 0)
                {
                    continue;
                }
                at = &atoms->atom[atomIndex];
                if (atoms->nres > 0)
                {
                    resInfo = &atoms->resinfo[at->resind];
                    /* FIXME: When TNG supports both residue index and residue
                     * number the latter should be used. */
                    res_name = *resInfo->name;
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
                    /* Since there is ONE chain for selection groups do not keep the
                     * original residue IDs - otherwise there might be conflicts. */
                    tng_chain_residue_add(tng, chain, res_name, &res);
                }
                tng_residue_atom_w_id_add(tng, res, *(atoms->atomname[atomIndex]),
                                          *(atoms->atomtype[atomIndex]),
                                          atom_offset + atomIndex, &atom);
                bAtomsAdded = TRUE;
            }
            /* Add bonds. */
            if (bAtomsAdded)
            {
                for (int k = 0; k < F_NRE; k++)
                {
                    if (IS_CHEMBOND(k))
                    {
                        ilist = &molType->ilist[k];
                        if (ilist)
                        {
                            int l = 1;
                            while (l < ilist->nr)
                            {
                                int atom1, atom2;
                                atom1 = ilist->iatoms[l] + atom_offset;
                                atom2 = ilist->iatoms[l+1] + atom_offset;
                                if (ggrpnr(&mtop->groups, egcCompressedX, atom1) == 0 &&
                                    ggrpnr(&mtop->groups, egcCompressedX, atom2) == 0)
                                {
                                    tng_molecule_bond_add(tng, mol, ilist->iatoms[l],
                                                          ilist->iatoms[l+1], &tngBond);
                                }
                                l += 3;
                            }
                        }
                    }
                }
                /* Settle is described using three atoms */
                ilist = &molType->ilist[F_SETTLE];
                if (ilist)
                {
                    int l = 1;
                    while (l < ilist->nr)
                    {
                        int atom1, atom2, atom3;
                        atom1 = ilist->iatoms[l] + atom_offset;
                        atom2 = ilist->iatoms[l+1] + atom_offset;
                        atom3 = ilist->iatoms[l+2] + atom_offset;
                        if (ggrpnr(&mtop->groups, egcCompressedX, atom1) == 0)
                        {
                            if (ggrpnr(&mtop->groups, egcCompressedX, atom2) == 0)
                            {
                                tng_molecule_bond_add(tng, mol, atom1,
                                                      atom2, &tngBond);
                            }
                            if (ggrpnr(&mtop->groups, egcCompressedX, atom3) == 0)
                            {
                                tng_molecule_bond_add(tng, mol, atom1,
                                                      atom3, &tngBond);
                            }
                        }
                        l += 4;
                    }
                }
            }
            atom_offset += atoms->nr;
        }
    }
    tng_molecule_existing_add(tng, &mol);
    tng_molecule_cnt_set(tng, mol, 1);
    tng_num_molecule_types_get(tng, &nMols);
    for (gmx_int64_t k = 0; k < nMols; k++)
    {
        tng_molecule_of_index_get(tng, k, &iterMol);
        if (iterMol == mol)
        {
            continue;
        }
        tng_molecule_cnt_set(tng, iterMol, 0);
    }
}
#endif

void gmx_tng_set_compression_precision(tng_trajectory_t tng,
                                       real             prec)
{
#ifdef GMX_USE_TNG
    tng_compression_precision_set(tng, prec);
#else
    GMX_UNUSED_VALUE(tng);
    GMX_UNUSED_VALUE(prec);
#endif
}

void gmx_tng_prepare_low_prec_writing(tng_trajectory_t  tng,
                                      const gmx_mtop_t *mtop,
                                      const t_inputrec *ir)
{
#ifdef GMX_USE_TNG
    gmx_tng_add_mtop(tng, mtop);
    add_selection_groups(tng, mtop);
    set_writing_intervals(tng, TRUE, ir);
    tng_time_per_frame_set(tng, ir->delta_t * PICO);
    gmx_tng_set_compression_precision(tng, ir->x_compression_precision);
#else
    GMX_UNUSED_VALUE(tng);
    GMX_UNUSED_VALUE(mtop);
    GMX_UNUSED_VALUE(ir);
#endif
}

void gmx_fwrite_tng(tng_trajectory_t tng,
                    const gmx_bool   bUseLossyCompression,
                    int              step,
                    real             elapsedPicoSeconds,
                    real             lambda,
                    const rvec      *box,
                    int              nAtoms,
                    const rvec      *x,
                    const rvec      *v,
                    const rvec      *f)
{
#ifdef GMX_USE_TNG
    typedef tng_function_status (*write_data_func_pointer)(tng_trajectory_t,
                                                           const gmx_int64_t,
                                                           const double,
                                                           const real*,
                                                           const gmx_int64_t,
                                                           const gmx_int64_t,
                                                           const char*,
                                                           const char,
                                                           const char);
#ifdef GMX_DOUBLE
    static write_data_func_pointer           write_data           = tng_util_generic_with_time_double_write;
#else
    static write_data_func_pointer           write_data           = tng_util_generic_with_time_write;
#endif
    double                                   elapsedSeconds = elapsedPicoSeconds * PICO;
    gmx_int64_t                              nParticles;
    char                                     compression;


    if (!tng)
    {
        /* This function might get called when the type of the
           compressed trajectory is actually XTC. So we exit and move
           on. */
        return;
    }

    tng_num_particles_get(tng, &nParticles);
    if (nAtoms != (int)nParticles)
    {
        tng_implicit_num_particles_set(tng, nAtoms);
    }

    if (bUseLossyCompression)
    {
        compression = TNG_TNG_COMPRESSION;
    }
    else
    {
        compression = TNG_GZIP_COMPRESSION;
    }

    /* The writing is done using write_data, which writes float or double
     * depending on the GROMACS compilation. */
    if (x)
    {
        GMX_ASSERT(box, "Need a non-NULL box if positions are written");

        if (write_data(tng, step, elapsedSeconds,
                       reinterpret_cast<const real *>(x),
                       3, TNG_TRAJ_POSITIONS, "POSITIONS",
                       TNG_PARTICLE_BLOCK_DATA,
                       compression) != TNG_SUCCESS)
        {
            gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
        }
        /* TNG-MF1 compression only compresses positions and velocities. Use lossless
         * compression for box shape regardless of output mode */
        if (write_data(tng, step, elapsedSeconds,
                       reinterpret_cast<const real *>(box),
                       9, TNG_TRAJ_BOX_SHAPE, "BOX SHAPE",
                       TNG_NON_PARTICLE_BLOCK_DATA,
                       TNG_GZIP_COMPRESSION) != TNG_SUCCESS)
        {
            gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
        }
    }

    if (v)
    {
        if (write_data(tng, step, elapsedSeconds,
                       reinterpret_cast<const real *>(v),
                       3, TNG_TRAJ_VELOCITIES, "VELOCITIES",
                       TNG_PARTICLE_BLOCK_DATA,
                       compression) != TNG_SUCCESS)
        {
            gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
        }
    }

    if (f)
    {
        /* TNG-MF1 compression only compresses positions and velocities. Use lossless
         * compression for forces regardless of output mode */
        if (write_data(tng, step, elapsedSeconds,
                       reinterpret_cast<const real *>(f),
                       3, TNG_TRAJ_FORCES, "FORCES",
                       TNG_PARTICLE_BLOCK_DATA,
                       TNG_GZIP_COMPRESSION) != TNG_SUCCESS)
        {
            gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
        }
    }

    /* TNG-MF1 compression only compresses positions and velocities. Use lossless
     * compression for lambdas regardless of output mode */
    if (write_data(tng, step, elapsedSeconds,
                   reinterpret_cast<const real *>(&lambda),
                   1, TNG_GMX_LAMBDA, "LAMBDAS",
                   TNG_NON_PARTICLE_BLOCK_DATA,
                   TNG_GZIP_COMPRESSION) != TNG_SUCCESS)
    {
        gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
    }
#else
    GMX_UNUSED_VALUE(tng);
    GMX_UNUSED_VALUE(bUseLossyCompression);
    GMX_UNUSED_VALUE(step);
    GMX_UNUSED_VALUE(elapsedPicoSeconds);
    GMX_UNUSED_VALUE(lambda);
    GMX_UNUSED_VALUE(box);
    GMX_UNUSED_VALUE(nAtoms);
    GMX_UNUSED_VALUE(x);
    GMX_UNUSED_VALUE(v);
    GMX_UNUSED_VALUE(f);
#endif
}

void fflush_tng(tng_trajectory_t tng)
{
#ifdef GMX_USE_TNG
    if (!tng)
    {
        return;
    }
    tng_frame_set_premature_write(tng, TNG_USE_HASH);
#else
    GMX_UNUSED_VALUE(tng);
#endif
}

float gmx_tng_get_time_of_final_frame(tng_trajectory_t tng)
{
#ifdef GMX_USE_TNG
    gmx_int64_t nFrames;
    double      time;
    float       fTime;

    tng_num_frames_get(tng, &nFrames);
    tng_util_time_of_frame_get(tng, nFrames - 1, &time);

    fTime = time / PICO;
    return fTime;
#else
    GMX_UNUSED_VALUE(tng);
    return -1.0;
#endif
}
