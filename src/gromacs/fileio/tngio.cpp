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

#include <cmath>

#include <algorithm>
#include <iterator>
#include <memory>
#include <vector>

#if GMX_USE_TNG
#include "tng/tng_io.h"
#endif

#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/sysinfo.h"
#include "gromacs/utility/unique_cptr.h"

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
#if GMX_USE_TNG
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
                  "File I/O error while opening %s for %s",
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
#if GMX_DOUBLE
        precisionString = " (double precision)";
#endif
        sprintf(programInfo, "%.100s %.128s%.24s",
                gmx::getProgramContext().displayName(),
                gmx_version(), precisionString);
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
#if GMX_USE_TNG
    if (tng)
    {
        tng_util_trajectory_close(tng);
    }
#else
    GMX_UNUSED_VALUE(tng);
#endif
}

#if GMX_USE_TNG
static void addTngMoleculeFromTopology(tng_trajectory_t     tng,
                                       const char          *moleculeName,
                                       const t_atoms       *atoms,
                                       gmx_int64_t          numMolecules,
                                       tng_molecule_t      *tngMol)
{
    tng_chain_t      tngChain = nullptr;
    tng_residue_t    tngRes   = nullptr;

    if (tng_molecule_add(tng, moleculeName, tngMol) != TNG_SUCCESS)
    {
        gmx_file("Cannot add molecule to TNG molecular system.");
    }

    for (int atomIndex = 0; atomIndex < atoms->nr; atomIndex++)
    {
        const t_atom *at = &atoms->atom[atomIndex];
        /* FIXME: Currently the TNG API can only add atoms belonging to a
         * residue and chain. Wait for TNG 2.0*/
        if (atoms->nres > 0)
        {
            const t_resinfo *resInfo        = &atoms->resinfo[at->resind];
            char             chainName[2]   = {resInfo->chainid, 0};
            tng_atom_t       tngAtom        = nullptr;
            t_atom          *prevAtom;

            if (atomIndex > 0)
            {
                prevAtom = &atoms->atom[atomIndex - 1];
            }
            else
            {
                prevAtom = nullptr;
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
    int                i;
    int                j;
    std::vector<real>  atomCharges;
    std::vector<real>  atomMasses;
    const t_ilist     *ilist;
    tng_bond_t         tngBond;
    char               datatype;

    if (!mtop)
    {
        /* No topology information available to add. */
        return;
    }

#if GMX_DOUBLE
    datatype = TNG_DOUBLE_DATA;
#else
    datatype = TNG_FLOAT_DATA;
#endif

    atomCharges.reserve(mtop->natoms);
    atomMasses.reserve(mtop->natoms);

    for (int molIndex = 0; molIndex < mtop->nmolblock; molIndex++)
    {
        tng_molecule_t       tngMol  = nullptr;
        const gmx_moltype_t *molType = &mtop->moltype[mtop->molblock[molIndex].type];

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
        /* First copy atom charges and masses, first atom by atom and then copy the memory for the molecule instances.
         * FIXME: Atom B state data should also be written to TNG (v 2.0?) */
        for (int atomCounter = 0; atomCounter < molType->atoms.nr; atomCounter++)
        {
            atomCharges.push_back(molType->atoms.atom[atomCounter].q);
            atomMasses.push_back(molType->atoms.atom[atomCounter].m);
        }
        for (int molCounter = 1; molCounter < mtop->molblock[molIndex].nmol; molCounter++)
        {
            std::copy_n(atomCharges.end() - molType->atoms.nr, molType->atoms.nr, std::back_inserter(atomCharges));
            std::copy_n(atomMasses.end()  - molType->atoms.nr, molType->atoms.nr, std::back_inserter(atomMasses));
        }
    }
    /* Write the TNG data blocks. */
    tng_particle_data_block_add(tng, TNG_TRAJ_PARTIAL_CHARGES, "PARTIAL CHARGES",
                                datatype, TNG_NON_TRAJECTORY_BLOCK,
                                1, 1, 1, 0, mtop->natoms,
                                TNG_GZIP_COMPRESSION, atomCharges.data());
    tng_particle_data_block_add(tng, TNG_TRAJ_MASSES, "ATOM MASSES",
                                datatype, TNG_NON_TRAJECTORY_BLOCK,
                                1, 1, 1, 0, mtop->natoms,
                                TNG_GZIP_COMPRESSION, atomMasses.data());
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
#if GMX_DOUBLE
    set_writing_interval_func_pointer set_writing_interval = tng_util_generic_write_interval_double_set;
#else
    set_writing_interval_func_pointer set_writing_interval = tng_util_generic_write_interval_set;
#endif
    int  xout, vout, fout;
    int  gcd = -1, lowest = -1;
    char compression;

    tng_set_frames_per_frame_set(tng, bUseLossyCompression, ir);

    if (bUseLossyCompression)
    {
        xout        = ir->nstxout_compressed;

        /* If there is no uncompressed coordinate output write forces
           and velocities to the compressed tng file. */
        if (ir->nstxout)
        {
            vout        = 0;
            fout        = 0;
        }
        else
        {
            vout        = ir->nstvout;
            fout        = ir->nstfout;
        }
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
        /* TODO: if/when we write energies to TNG also, reconsider how
         * and when box information is written, because GROMACS
         * behaviour pre-5.0 was to write the box with every
         * trajectory frame and every energy frame, and probably
         * people depend on this. */

        gcd = greatest_common_divisor_if_positive(gcd, xout);
        if (lowest < 0 || xout < lowest)
        {
            lowest = xout;
        }
    }
    if (vout)
    {
        set_writing_interval(tng, vout, 3, TNG_TRAJ_VELOCITIES,
                             "VELOCITIES", TNG_PARTICLE_BLOCK_DATA,
                             compression);

        gcd = greatest_common_divisor_if_positive(gcd, vout);
        if (lowest < 0 || vout < lowest)
        {
            lowest = vout;
        }
    }
    if (fout)
    {
        set_writing_interval(tng, fout, 3, TNG_TRAJ_FORCES,
                             "FORCES", TNG_PARTICLE_BLOCK_DATA,
                             TNG_GZIP_COMPRESSION);

        gcd = greatest_common_divisor_if_positive(gcd, fout);
        if (lowest < 0 || fout < lowest)
        {
            lowest = fout;
        }
    }
    if (gcd > 0)
    {
        /* Lambdas and box shape written at an interval of the lowest common
           denominator of other output */
        set_writing_interval(tng, gcd, 1, TNG_GMX_LAMBDA,
                             "LAMBDAS", TNG_NON_PARTICLE_BLOCK_DATA,
                             TNG_GZIP_COMPRESSION);

        set_writing_interval(tng, gcd, 9, TNG_TRAJ_BOX_SHAPE,
                             "BOX SHAPE", TNG_NON_PARTICLE_BLOCK_DATA,
                             TNG_GZIP_COMPRESSION);
        if (gcd < lowest / 10)
        {
            gmx_warning("The lowest common denominator of trajectory output is "
                        "every %d step(s), whereas the shortest output interval "
                        "is every %d steps.", gcd, lowest);
        }
    }
}
#endif

void gmx_tng_prepare_md_writing(tng_trajectory_t  tng,
                                const gmx_mtop_t *mtop,
                                const t_inputrec *ir)
{
#if GMX_USE_TNG
    gmx_tng_add_mtop(tng, mtop);
    set_writing_intervals(tng, FALSE, ir);
    tng_time_per_frame_set(tng, ir->delta_t * PICO);
#else
    GMX_UNUSED_VALUE(tng);
    GMX_UNUSED_VALUE(mtop);
    GMX_UNUSED_VALUE(ir);
#endif
}

#if GMX_USE_TNG
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
#if GMX_USE_TNG
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
#if GMX_USE_TNG
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
                    gmx_int64_t      step,
                    real             elapsedPicoSeconds,
                    real             lambda,
                    const rvec      *box,
                    int              nAtoms,
                    const rvec      *x,
                    const rvec      *v,
                    const rvec      *f)
{
#if GMX_USE_TNG
    typedef tng_function_status (*write_data_func_pointer)(tng_trajectory_t,
                                                           const gmx_int64_t,
                                                           const double,
                                                           const real*,
                                                           const gmx_int64_t,
                                                           const gmx_int64_t,
                                                           const char*,
                                                           const char,
                                                           const char);
#if GMX_DOUBLE
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
     * compression for lambdas and box shape regardless of output mode */
    if (write_data(tng, step, elapsedSeconds,
                   reinterpret_cast<const real *>(box),
                   9, TNG_TRAJ_BOX_SHAPE, "BOX SHAPE",
                   TNG_NON_PARTICLE_BLOCK_DATA,
                   TNG_GZIP_COMPRESSION) != TNG_SUCCESS)
    {
        gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
    }

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
#if GMX_USE_TNG
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
#if GMX_USE_TNG
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

void gmx_prepare_tng_writing(const char              *filename,
                             char                     mode,
                             tng_trajectory_t        *input,
                             tng_trajectory_t        *output,
                             int                      nAtoms,
                             const gmx_mtop_t        *mtop,
                             const int               *index,
                             const char              *indexGroupName)
{
#if GMX_USE_TNG
    /* FIXME after 5.0: Currently only standard block types are read */
    const int           defaultNumIds              = 5;
    static gmx_int64_t  fallbackIds[defaultNumIds] =
    {
        TNG_TRAJ_BOX_SHAPE, TNG_TRAJ_POSITIONS,
        TNG_TRAJ_VELOCITIES, TNG_TRAJ_FORCES,
        TNG_GMX_LAMBDA
    };
    static char         fallbackNames[defaultNumIds][32] =
    {
        "BOX SHAPE", "POSITIONS", "VELOCITIES",
        "FORCES", "LAMBDAS"
    };

    typedef tng_function_status (*set_writing_interval_func_pointer)(tng_trajectory_t,
                                                                     const gmx_int64_t,
                                                                     const gmx_int64_t,
                                                                     const gmx_int64_t,
                                                                     const char*,
                                                                     const char,
                                                                     const char);
#if GMX_DOUBLE
    set_writing_interval_func_pointer set_writing_interval = tng_util_generic_write_interval_double_set;
#else
    set_writing_interval_func_pointer set_writing_interval = tng_util_generic_write_interval_set;
#endif

    gmx_tng_open(filename, mode, output);

    /* Do we have an input file in TNG format? If so, then there's
       more data we can copy over, rather than having to improvise. */
    if (*input)
    {
        /* Set parameters (compression, time per frame, molecule
         * information, number of frames per frame set and writing
         * intervals of positions, box shape and lambdas) of the
         * output tng container based on their respective values int
         * the input tng container */
        double      time, compression_precision;
        gmx_int64_t n_frames_per_frame_set, interval = -1;

        tng_compression_precision_get(*input, &compression_precision);
        tng_compression_precision_set(*output, compression_precision);
        // TODO make this configurable in a future version
        char compression_type = TNG_TNG_COMPRESSION;

        tng_molecule_system_copy(*input, *output);

        tng_time_per_frame_get(*input, &time);
        tng_time_per_frame_set(*output, time);

        tng_num_frames_per_frame_set_get(*input, &n_frames_per_frame_set);
        tng_num_frames_per_frame_set_set(*output, n_frames_per_frame_set);

        for (int i = 0; i < defaultNumIds; i++)
        {
            if (tng_data_get_stride_length(*input, fallbackIds[i], -1, &interval)
                == TNG_SUCCESS)
            {
                switch (fallbackIds[i])
                {
                    case TNG_TRAJ_POSITIONS:
                    case TNG_TRAJ_VELOCITIES:
                        set_writing_interval(*output, interval, 3, fallbackIds[i],
                                             fallbackNames[i], TNG_PARTICLE_BLOCK_DATA,
                                             compression_type);
                        break;
                    case TNG_TRAJ_FORCES:
                        set_writing_interval(*output, interval, 3, fallbackIds[i],
                                             fallbackNames[i], TNG_PARTICLE_BLOCK_DATA,
                                             TNG_GZIP_COMPRESSION);
                        break;
                    case TNG_TRAJ_BOX_SHAPE:
                        set_writing_interval(*output, interval, 9, fallbackIds[i],
                                             fallbackNames[i], TNG_NON_PARTICLE_BLOCK_DATA,
                                             TNG_GZIP_COMPRESSION);
                        break;
                    case TNG_GMX_LAMBDA:
                        set_writing_interval(*output, interval, 1, fallbackIds[i],
                                             fallbackNames[i], TNG_NON_PARTICLE_BLOCK_DATA,
                                             TNG_GZIP_COMPRESSION);
                        break;
                    default:
                        continue;
                }
            }
        }

    }
    else
    {
        /* TODO after trjconv is modularized: fix this so the user can
           change precision when they are doing an operation where
           this makes sense, and not otherwise.

           char compression = bUseLossyCompression ? TNG_TNG_COMPRESSION : TNG_GZIP_COMPRESSION;
           gmx_tng_set_compression_precision(*output, ndec2prec(nDecimalsOfPrecision));
         */
        gmx_tng_add_mtop(*output, mtop);
        tng_num_frames_per_frame_set_set(*output, 1);
    }

    if (index && nAtoms > 0)
    {
        gmx_tng_setup_atom_subgroup(*output, nAtoms, index, indexGroupName);
    }

    /* If for some reason there are more requested atoms than there are atoms in the
     * molecular system create a number of implicit atoms (without atom data) to
     * compensate for that. */
    if (nAtoms >= 0)
    {
        tng_implicit_num_particles_set(*output, nAtoms);
    }
#else
    GMX_UNUSED_VALUE(filename);
    GMX_UNUSED_VALUE(mode);
    GMX_UNUSED_VALUE(input);
    GMX_UNUSED_VALUE(output);
    GMX_UNUSED_VALUE(nAtoms);
    GMX_UNUSED_VALUE(mtop);
    GMX_UNUSED_VALUE(index);
    GMX_UNUSED_VALUE(indexGroupName);
#endif
}

void gmx_write_tng_from_trxframe(tng_trajectory_t        output,
                                 const t_trxframe       *frame,
                                 int                     natoms)
{
#if GMX_USE_TNG
    if (frame->step > 0)
    {
        double timePerFrame = frame->time * PICO / frame->step;
        tng_time_per_frame_set(output, timePerFrame);
    }
    if (natoms < 0)
    {
        natoms = frame->natoms;
    }
    gmx_fwrite_tng(output,
                   TRUE,
                   frame->step,
                   frame->time,
                   0,
                   frame->box,
                   natoms,
                   frame->x,
                   frame->v,
                   frame->f);
#else
    GMX_UNUSED_VALUE(output);
    GMX_UNUSED_VALUE(frame);
    GMX_UNUSED_VALUE(natoms);
#endif
}

namespace
{

#if GMX_USE_TNG
void
convert_array_to_real_array(void       *from,
                            real       *to,
                            const float fact,
                            const int   nAtoms,
                            const int   nValues,
                            const char  datatype)
{
    int        i, j;

    const bool useDouble = GMX_DOUBLE;
    switch (datatype)
    {
        case TNG_FLOAT_DATA:
            if (!useDouble)
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
                            to[i*nValues+j] = reinterpret_cast<float *>(from)[i*nValues+j] * fact;
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
                        to[i*nValues+j] = reinterpret_cast<float *>(from)[i*nValues+j] * fact;
                    }
                }
            }
            break;
        case TNG_INT_DATA:
            for (i = 0; i < nAtoms; i++)
            {
                for (j = 0; j < nValues; j++)
                {
                    to[i*nValues+j] = reinterpret_cast<gmx_int64_t *>(from)[i*nValues+j] * fact;
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
                            to[i*nValues+j] = reinterpret_cast<double *>(from)[i*nValues+j] * fact;
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
                        to[i*nValues+j] = reinterpret_cast<double *>(from)[i*nValues+j] * fact;
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
    gmx_int64_t exp = -1;
    real        distanceScaleFactor;

    // TODO Hopefully, TNG 2.0 will do this kind of thing for us
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

} // namespace

void gmx_tng_setup_atom_subgroup(tng_trajectory_t tng,
                                 const int        nind,
                                 const int       *ind,
                                 const char      *name)
{
#if GMX_USE_TNG
    gmx_int64_t              nAtoms, cnt, nMols;
    tng_molecule_t           mol, iterMol;
    tng_chain_t              chain;
    tng_residue_t            res;
    tng_atom_t               atom;
    tng_function_status      stat;

    tng_num_particles_get(tng, &nAtoms);

    if (nAtoms == nind)
    {
        return;
    }

    stat = tng_molecule_find(tng, name, -1, &mol);
    if (stat == TNG_SUCCESS)
    {
        tng_molecule_num_atoms_get(tng, mol, &nAtoms);
        tng_molecule_cnt_get(tng, mol, &cnt);
        if (nAtoms == nind)
        {
            stat = TNG_SUCCESS;
        }
        else
        {
            stat = TNG_FAILURE;
        }
    }
    if (stat == TNG_FAILURE)
    {
        /* The indexed atoms are added to one separate molecule. */
        tng_molecule_alloc(tng, &mol);
        tng_molecule_name_set(tng, mol, name);
        tng_molecule_chain_add(tng, mol, "", &chain);

        for (int i = 0; i < nind; i++)
        {
            char        temp_name[256], temp_type[256];

            /* Try to retrieve the residue name of the atom */
            stat = tng_residue_name_of_particle_nr_get(tng, ind[i], temp_name, 256);
            if (stat != TNG_SUCCESS)
            {
                temp_name[0] = '\0';
            }
            /* Check if the molecule of the selection already contains this residue */
            if (tng_chain_residue_find(tng, chain, temp_name, -1, &res)
                != TNG_SUCCESS)
            {
                tng_chain_residue_add(tng, chain, temp_name, &res);
            }
            /* Try to find the original name and type of the atom */
            stat = tng_atom_name_of_particle_nr_get(tng, ind[i], temp_name, 256);
            if (stat != TNG_SUCCESS)
            {
                temp_name[0] = '\0';
            }
            stat = tng_atom_type_of_particle_nr_get(tng, ind[i], temp_type, 256);
            if (stat != TNG_SUCCESS)
            {
                temp_type[0] = '\0';
            }
            tng_residue_atom_w_id_add(tng, res, temp_name, temp_type, ind[i], &atom);
        }
        tng_molecule_existing_add(tng, &mol);
    }
    /* Set the count of the molecule containing the selected atoms to 1 and all
     * other molecules to 0 */
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
#else
    GMX_UNUSED_VALUE(tng);
    GMX_UNUSED_VALUE(nind);
    GMX_UNUSED_VALUE(ind);
    GMX_UNUSED_VALUE(name);
#endif
}

/* TODO: If/when TNG acquires the ability to copy data blocks without
 * uncompressing them, then this implemenation should be reconsidered.
 * Ideally, gmx trjconv -f a.tng -o b.tng -b 10 -e 20 would be fast
 * and lose no information. */
gmx_bool gmx_read_next_tng_frame(tng_trajectory_t            input,
                                 t_trxframe                 *fr,
                                 gmx_int64_t                *requestedIds,
                                 int                         numRequestedIds)
{
#if GMX_USE_TNG
    gmx_bool                bOK = TRUE;
    tng_function_status     stat;
    gmx_int64_t             numberOfAtoms = -1, frameNumber = -1;
    gmx_int64_t             nBlocks, blockId, *blockIds = nullptr, codecId;
    char                    datatype      = -1;
    void                   *values        = nullptr;
    double                  frameTime     = -1.0;
    int                     size, blockDependency;
    double                  prec;
    const int               defaultNumIds = 5;
    static gmx_int64_t      fallbackRequestedIds[defaultNumIds] =
    {
        TNG_TRAJ_BOX_SHAPE, TNG_TRAJ_POSITIONS,
        TNG_TRAJ_VELOCITIES, TNG_TRAJ_FORCES,
        TNG_GMX_LAMBDA
    };


    fr->bStep     = FALSE;
    fr->bTime     = FALSE;
    fr->bLambda   = FALSE;
    fr->bAtoms    = FALSE;
    fr->bPrec     = FALSE;
    fr->bX        = FALSE;
    fr->bV        = FALSE;
    fr->bF        = FALSE;
    fr->bBox      = FALSE;

    /* If no specific IDs were requested read all block types that can
     * currently be interpreted */
    if (!requestedIds || numRequestedIds == 0)
    {
        numRequestedIds = defaultNumIds;
        requestedIds    = fallbackRequestedIds;
    }

    stat = tng_num_particles_get(input, &numberOfAtoms);
    if (stat != TNG_SUCCESS)
    {
        gmx_file("Cannot determine number of atoms from TNG file.");
    }
    fr->natoms = numberOfAtoms;

    bool nextFrameExists = gmx_get_tng_data_block_types_of_next_frame(input,
                                                                      fr->step,
                                                                      numRequestedIds,
                                                                      requestedIds,
                                                                      &frameNumber,
                                                                      &nBlocks,
                                                                      &blockIds);
    gmx::unique_cptr<gmx_int64_t, gmx::free_wrapper> blockIdsGuard(blockIds);
    if (!nextFrameExists)
    {
        return FALSE;
    }

    if (nBlocks == 0)
    {
        return FALSE;
    }

    for (gmx_int64_t i = 0; i < nBlocks; i++)
    {
        blockId = blockIds[i];
        tng_data_block_dependency_get(input, blockId, &blockDependency);
        if (blockDependency & TNG_PARTICLE_DEPENDENT)
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
            gmx_file("Cannot read positions from TNG file.");
            return FALSE;
        }
        else if (stat == TNG_FAILURE)
        {
            continue;
        }
        switch (blockId)
        {
            case TNG_TRAJ_BOX_SHAPE:
                switch (datatype)
                {
                    case TNG_INT_DATA:
                        size = sizeof(gmx_int64_t);
                        break;
                    case TNG_FLOAT_DATA:
                        size = sizeof(float);
                        break;
                    case TNG_DOUBLE_DATA:
                        size = sizeof(double);
                        break;
                    default:
                        gmx_incons("Illegal datatype of box shape values!");
                }
                for (int i = 0; i < DIM; i++)
                {
                    convert_array_to_real_array(reinterpret_cast<char *>(values) + size * i * DIM,
                                                reinterpret_cast<real *>(fr->box[i]),
                                                getDistanceScaleFactor(input),
                                                1,
                                                DIM,
                                                datatype);
                }
                fr->bBox = TRUE;
                break;
            case TNG_TRAJ_POSITIONS:
                srenew(fr->x, fr->natoms);
                convert_array_to_real_array(values,
                                            reinterpret_cast<real *>(fr->x),
                                            getDistanceScaleFactor(input),
                                            fr->natoms,
                                            DIM,
                                            datatype);
                fr->bX = TRUE;
                tng_util_frame_current_compression_get(input, blockId, &codecId, &prec);
                /* This must be updated if/when more lossy compression methods are added */
                if (codecId == TNG_TNG_COMPRESSION)
                {
                    fr->prec  = prec;
                    fr->bPrec = TRUE;
                }
                break;
            case TNG_TRAJ_VELOCITIES:
                srenew(fr->v, fr->natoms);
                convert_array_to_real_array(values,
                                            (real *) fr->v,
                                            getDistanceScaleFactor(input),
                                            fr->natoms,
                                            DIM,
                                            datatype);
                fr->bV = TRUE;
                tng_util_frame_current_compression_get(input, blockId, &codecId, &prec);
                /* This must be updated if/when more lossy compression methods are added */
                if (codecId == TNG_TNG_COMPRESSION)
                {
                    fr->prec  = prec;
                    fr->bPrec = TRUE;
                }
                break;
            case TNG_TRAJ_FORCES:
                srenew(fr->f, fr->natoms);
                convert_array_to_real_array(values,
                                            reinterpret_cast<real *>(fr->f),
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
                        fr->lambda = *(reinterpret_cast<float *>(values));
                        break;
                    case TNG_DOUBLE_DATA:
                        fr->lambda = *(reinterpret_cast<double *>(values));
                        break;
                    default:
                        gmx_incons("Illegal datatype lambda value!");
                }
                fr->bLambda = TRUE;
                break;
            default:
                gmx_warning("Illegal block type! Currently GROMACS tools can only handle certain data types. Skipping block.");
        }
        /* values does not have to be freed before reading next frame. It will
         * be reallocated if it is not NULL. */
    }

    fr->step  = frameNumber;
    fr->bStep = TRUE;
    // Convert the time to ps
    fr->time  = frameTime / PICO;
    fr->bTime = TRUE;

    // TODO This does not leak, but is not exception safe.
    /* values must be freed before leaving this function */
    sfree(values);

    return bOK;
#else
    GMX_UNUSED_VALUE(input);
    GMX_UNUSED_VALUE(fr);
    GMX_UNUSED_VALUE(requestedIds);
    GMX_UNUSED_VALUE(numRequestedIds);
    return FALSE;
#endif
}

void gmx_print_tng_molecule_system(tng_trajectory_t input,
                                   FILE            *stream)
{
#if GMX_USE_TNG
    gmx_int64_t         nMolecules, nChains, nResidues, nAtoms, nFramesRead;
    gmx_int64_t         strideLength, nParticlesRead, nValuesPerFrameRead, *molCntList;
    tng_molecule_t      molecule;
    tng_chain_t         chain;
    tng_residue_t       residue;
    tng_atom_t          atom;
    tng_function_status stat;
    char                str[256];
    char                varNAtoms;
    char                datatype;
    void               *data = nullptr;
    std::vector<real>   atomCharges;
    std::vector<real>   atomMasses;

    tng_num_molecule_types_get(input, &nMolecules);
    tng_molecule_cnt_list_get(input, &molCntList);
    /* Can the number of particles change in the trajectory or is it constant? */
    tng_num_particles_variable_get(input, &varNAtoms);

    for (gmx_int64_t i = 0; i < nMolecules; i++)
    {
        tng_molecule_of_index_get(input, i, &molecule);
        tng_molecule_name_get(input, molecule, str, 256);
        if (varNAtoms == TNG_CONSTANT_N_ATOMS)
        {
            if ((int)molCntList[i] == 0)
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
        if (nChains > 0)
        {
            for (gmx_int64_t j = 0; j < nChains; j++)
            {
                tng_molecule_chain_of_index_get(input, molecule, j, &chain);
                tng_chain_name_get(input, chain, str, 256);
                fprintf(stream, "\tChain: %s\n", str);
                tng_chain_num_residues_get(input, chain, &nResidues);
                for (gmx_int64_t k = 0; k < nResidues; k++)
                {
                    tng_chain_residue_of_index_get(input, chain, k, &residue);
                    tng_residue_name_get(input, residue, str, 256);
                    fprintf(stream, "\t\tResidue: %s\n", str);
                    tng_residue_num_atoms_get(input, residue, &nAtoms);
                    for (gmx_int64_t l = 0; l < nAtoms; l++)
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
            if (nResidues > 0)
            {
                for (gmx_int64_t k = 0; k < nResidues; k++)
                {
                    tng_molecule_residue_of_index_get(input, molecule, k, &residue);
                    tng_residue_name_get(input, residue, str, 256);
                    fprintf(stream, "\t\tResidue: %s\n", str);
                    tng_residue_num_atoms_get(input, residue, &nAtoms);
                    for (gmx_int64_t l = 0; l < nAtoms; l++)
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
                for (gmx_int64_t l = 0; l < nAtoms; l++)
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

    tng_num_particles_get(input, &nAtoms);
    stat = tng_particle_data_vector_get(input, TNG_TRAJ_PARTIAL_CHARGES, &data, &nFramesRead,
                                        &strideLength, &nParticlesRead,
                                        &nValuesPerFrameRead, &datatype);
    if (stat == TNG_SUCCESS)
    {
        atomCharges.resize(nAtoms);
        convert_array_to_real_array(data,
                                    atomCharges.data(),
                                    1,
                                    nAtoms,
                                    1,
                                    datatype);

        fprintf(stream, "Atom Charges (%d):\n", int(nAtoms));
        for (gmx_int64_t i = 0; i < nAtoms; i += 10)
        {
            fprintf(stream, "Atom Charges [%8d-]=[", int(i));
            for (gmx_int64_t j = 0; (j < 10 && i + j < nAtoms); j++)
            {
                fprintf(stream, " %12.5e", atomCharges[i + j]);
            }
            fprintf(stream, "]\n");
        }
    }

    stat = tng_particle_data_vector_get(input, TNG_TRAJ_MASSES, &data, &nFramesRead,
                                        &strideLength, &nParticlesRead,
                                        &nValuesPerFrameRead, &datatype);
    if (stat == TNG_SUCCESS)
    {
        atomMasses.resize(nAtoms);
        convert_array_to_real_array(data,
                                    atomMasses.data(),
                                    1,
                                    nAtoms,
                                    1,
                                    datatype);

        fprintf(stream, "Atom Masses (%d):\n", int(nAtoms));
        for (gmx_int64_t i = 0; i < nAtoms; i += 10)
        {
            fprintf(stream, "Atom Masses [%8d-]=[", int(i));
            for (gmx_int64_t j = 0; (j < 10 && i + j < nAtoms); j++)
            {
                fprintf(stream, " %12.5e", atomMasses[i + j]);
            }
            fprintf(stream, "]\n");
        }
    }

    sfree(data);
#else
    GMX_UNUSED_VALUE(input);
    GMX_UNUSED_VALUE(stream);
#endif
}

gmx_bool gmx_get_tng_data_block_types_of_next_frame(tng_trajectory_t     input,
                                                    int                  frame,
                                                    int                  nRequestedIds,
                                                    gmx_int64_t         *requestedIds,
                                                    gmx_int64_t         *nextFrame,
                                                    gmx_int64_t         *nBlocks,
                                                    gmx_int64_t        **blockIds)
{
#if GMX_USE_TNG
    tng_function_status stat;

    stat = tng_util_trajectory_next_frame_present_data_blocks_find(input, frame,
                                                                   nRequestedIds, requestedIds,
                                                                   nextFrame,
                                                                   nBlocks, blockIds);

    if (stat == TNG_CRITICAL)
    {
        gmx_file("Cannot read TNG file. Cannot find data blocks of next frame.");
    }
    else if (stat == TNG_FAILURE)
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

gmx_bool gmx_get_tng_data_next_frame_of_block_type(tng_trajectory_t     input,
                                                   gmx_int64_t          blockId,
                                                   real               **values,
                                                   gmx_int64_t         *frameNumber,
                                                   double              *frameTime,
                                                   gmx_int64_t         *nValuesPerFrame,
                                                   gmx_int64_t         *nAtoms,
                                                   real                *prec,
                                                   char                *name,
                                                   int                  maxLen,
                                                   gmx_bool            *bOK)
{
#if GMX_USE_TNG
    tng_function_status stat;
    char                datatype = -1;
    gmx_int64_t         codecId;
    int                 blockDependency;
    void               *data = nullptr;
    double              localPrec;

    stat = tng_data_block_name_get(input, blockId, name, maxLen);
    if (stat != TNG_SUCCESS)
    {
        gmx_file("Cannot read next frame of TNG file");
    }
    stat = tng_data_block_dependency_get(input, blockId, &blockDependency);
    if (stat != TNG_SUCCESS)
    {
        gmx_file("Cannot read next frame of TNG file");
    }
    if (blockDependency & TNG_PARTICLE_DEPENDENT)
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
        stat    = tng_util_non_particle_data_next_frame_read(input,
                                                             blockId,
                                                             &data,
                                                             &datatype,
                                                             frameNumber,
                                                             frameTime);
    }
    if (stat == TNG_CRITICAL)
    {
        gmx_file("Cannot read next frame of TNG file");
    }
    if (stat == TNG_FAILURE)
    {
        *bOK = TRUE;
        return FALSE;
    }

    stat = tng_data_block_num_values_per_frame_get(input, blockId, nValuesPerFrame);
    if (stat != TNG_SUCCESS)
    {
        gmx_file("Cannot read next frame of TNG file");
    }
    srenew(*values, sizeof(real) * *nValuesPerFrame * *nAtoms);
    convert_array_to_real_array(data,
                                *values,
                                getDistanceScaleFactor(input),
                                *nAtoms,
                                *nValuesPerFrame,
                                datatype);

    tng_util_frame_current_compression_get(input, blockId, &codecId, &localPrec);

    /* This must be updated if/when more lossy compression methods are added */
    if (codecId != TNG_TNG_COMPRESSION)
    {
        *prec = -1.0;
    }
    else
    {
        *prec = localPrec;
    }

    sfree(data);

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
