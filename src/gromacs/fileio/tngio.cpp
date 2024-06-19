/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "tngio.h"

#include "config.h"

#include <cmath>
#include <cstdint>
#include <cstring>

#include <algorithm>
#include <iterator>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"

#if GMX_USE_TNG
#    include "tng/tng_io.h"
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

#if !GMX_USE_TNG
using tng_trajectory_t = void*;
#endif

/*! \brief Gromacs Wrapper around tng datatype
 *
 * This could in principle hold any GROMACS-specific requirements not yet
 * implemented in or not relevant to the TNG library itself. However, for now
 * we only use it to handle some shortcomings we have discovered, where the TNG
 * API itself is a bit fragile and can end up overwriting data if called several
 * times with the same frame number.
 * The logic to determine the time per step was also a bit fragile. This is not
 * critical, but since we anyway need a wrapper for ensuring unique frame
 * numbers, we can also use it to store the time of the first step and use that
 * to derive a slightly better/safer estimate of the time per step.
 *
 * At some future point where we have a second-generation TNG API we should
 * consider removing this again.
 */
struct gmx_tng_trajectory
{
    tng_trajectory_t tng;                  //!< Actual TNG handle (pointer)
    bool             lastStepDataIsValid;  //!< True if lastStep has been set
    std::int64_t     lastStep;             //!< Index/step used for last frame
    bool             lastTimeDataIsValid;  //!< True if lastTime has been set
    double           lastTime;             //!< Time of last frame (TNG unit is seconds)
    bool             timePerFrameIsSet;    //!< True if we have set the time per frame
    int              boxOutputInterval;    //!< Number of steps between the output of box size
    int              lambdaOutputInterval; //!< Number of steps between the output of lambdas
};

#if GMX_USE_TNG
static const char* modeToVerb(char mode)
{
    const char* p;
    switch (mode)
    {
        case 'r': p = "reading"; break;
        case 'w': p = "writing"; break;
        case 'a': p = "appending"; break;
        default: gmx_fatal(FARGS, "Invalid file opening mode %c", mode);
    }
    return p;
}
#endif

#if !GMX_USE_TNG && defined(__clang__)
/* gmx_tng_open will not return when the build configuration did not
 * support TNG at all. That leads to warnings that suggest using
 * [[noreturn]]. But that attribute also has to be used in the header
 * file, and we do not want the header file to depend on the build
 * configuration. So we suppress the warning that gives the
 * suggestion. */
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wmissing-noreturn"
#endif

void gmx_tng_open(const std::filesystem::path& filename, char mode, gmx_tng_trajectory_t* gmx_tng)
{
#if GMX_USE_TNG
    /* First check whether we have to make a backup,
     * only for writing, not for read or append.
     */
    if (mode == 'w')
    {
        make_backup(filename);
    }

    *gmx_tng                        = new gmx_tng_trajectory;
    (*gmx_tng)->lastStepDataIsValid = false;
    (*gmx_tng)->lastTimeDataIsValid = false;
    (*gmx_tng)->timePerFrameIsSet   = false;
    tng_trajectory_t* tng           = &(*gmx_tng)->tng;

    /* tng must not be pointing at already allocated memory.
     * Memory will be allocated by tng_util_trajectory_open() and must
     * later on be freed by tng_util_trajectory_close(). */
    if (TNG_SUCCESS != tng_util_trajectory_open(filename.string().c_str(), mode, tng))
    {
        /* TNG does return more than one degree of error, but there is
           no use case for GROMACS handling the non-fatal errors
           gracefully. */
        gmx_fatal(FARGS, "File I/O error while opening %s for %s", filename.string().c_str(), modeToVerb(mode));
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
        const char* precisionString = "";
#    if GMX_DOUBLE
        precisionString = " (double precision)";
#    endif
        sprintf(programInfo, "%.100s %.128s%.24s", gmx::getProgramContext().displayName(), gmx_version(), precisionString);
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
    GMX_UNUSED_VALUE(gmx_tng);
#endif
}

#if !GMX_USE_TNG && defined(__clang__)
#    pragma clang diagnostic pop
#endif

void gmx_tng_close(gmx_tng_trajectory_t* gmx_tng)
{
    /* We have to check that tng is set because
     * tng_util_trajectory_close wants to return a NULL in it, and
     * gives a fatal error if it is NULL. */
#if GMX_USE_TNG
    if (gmx_tng == nullptr || *gmx_tng == nullptr)
    {
        return;
    }
    tng_trajectory_t* tng = &(*gmx_tng)->tng;

    if (tng)
    {
        tng_util_trajectory_close(tng);
    }
    delete *gmx_tng;
    *gmx_tng = nullptr;

#else
    GMX_UNUSED_VALUE(gmx_tng);
#endif
}

#if GMX_USE_TNG
static void addTngMoleculeFromTopology(gmx_tng_trajectory_t gmx_tng,
                                       const char*          moleculeName,
                                       const t_atoms*       atoms,
                                       int64_t              numMolecules,
                                       tng_molecule_t*      tngMol)
{
    tng_trajectory_t tng      = gmx_tng->tng;
    tng_chain_t      tngChain = nullptr;
    tng_residue_t    tngRes   = nullptr;

    if (tng_molecule_add(tng, moleculeName, tngMol) != TNG_SUCCESS)
    {
        gmx_file("Cannot add molecule to TNG molecular system.");
    }

    for (int atomIndex = 0; atomIndex < atoms->nr; atomIndex++)
    {
        const t_atom* at = &atoms->atom[atomIndex];
        /* FIXME: Currently the TNG API can only add atoms belonging to a
         * residue and chain. Wait for TNG 2.0*/
        if (atoms->nres > 0)
        {
            const t_resinfo* resInfo      = &atoms->resinfo[at->resind];
            char             chainName[2] = { resInfo->chainid, 0 };
            tng_atom_t       tngAtom      = nullptr;
            t_atom*          prevAtom;

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
                if (!prevAtom || resInfo->chainid != atoms->resinfo[prevAtom->resind].chainid)
                {
                    tng_molecule_chain_add(tng, *tngMol, chainName, &tngChain);
                }
                /* FIXME: When TNG supports both residue index and residue
                 * number the latter should be used. Wait for TNG 2.0*/
                tng_chain_residue_add(tng, tngChain, *resInfo->name, &tngRes);
            }
            tng_residue_atom_add(
                    tng, tngRes, *(atoms->atomname[atomIndex]), *(atoms->atomtype[atomIndex]), &tngAtom);
        }
    }
    tng_molecule_cnt_set(tng, *tngMol, numMolecules);
}

void gmx_tng_add_mtop(gmx_tng_trajectory_t gmx_tng, const gmx_mtop_t* mtop)
{
    int               i;
    int               j;
    std::vector<real> atomCharges;
    std::vector<real> atomMasses;
    tng_bond_t        tngBond;
    char              datatype;

    tng_trajectory_t tng = gmx_tng->tng;

    if (!mtop)
    {
        /* No topology information available to add. */
        return;
    }

#    if GMX_DOUBLE
    datatype = TNG_DOUBLE_DATA;
#    else
    datatype                                               = TNG_FLOAT_DATA;
#    endif

    atomCharges.reserve(mtop->natoms);
    atomMasses.reserve(mtop->natoms);

    for (const gmx_molblock_t& molBlock : mtop->molblock)
    {
        tng_molecule_t       tngMol  = nullptr;
        const gmx_moltype_t* molType = &mtop->moltype[molBlock.type];

        /* Add a molecule to the TNG trajectory with the same name as the
         * current molecule. */
        addTngMoleculeFromTopology(gmx_tng, *(molType->name), &molType->atoms, molBlock.nmol, &tngMol);

        /* Bonds have to be deduced from interactions (constraints etc). Different
         * interactions have different sets of parameters. */
        /* Constraints are specified using two atoms */
        for (i = 0; i < F_NRE; i++)
        {
            if (IS_CHEMBOND(i))
            {
                const InteractionList& ilist = molType->ilist[i];
                j                            = 1;
                while (j < ilist.size())
                {
                    tng_molecule_bond_add(tng, tngMol, ilist.iatoms[j], ilist.iatoms[j + 1], &tngBond);
                    j += 3;
                }
            }
        }
        /* Settle is described using three atoms */
        const InteractionList& ilist = molType->ilist[F_SETTLE];
        j                            = 1;
        while (j < ilist.size())
        {
            tng_molecule_bond_add(tng, tngMol, ilist.iatoms[j], ilist.iatoms[j + 1], &tngBond);
            tng_molecule_bond_add(tng, tngMol, ilist.iatoms[j], ilist.iatoms[j + 2], &tngBond);
            j += 4;
        }
        /* First copy atom charges and masses, first atom by atom and then copy the memory for the molecule instances.
         * FIXME: Atom B state data should also be written to TNG (v 2.0?) */
        for (int atomCounter = 0; atomCounter < molType->atoms.nr; atomCounter++)
        {
            atomCharges.push_back(molType->atoms.atom[atomCounter].q);
            atomMasses.push_back(molType->atoms.atom[atomCounter].m);
        }
        for (int molCounter = 1; molCounter < molBlock.nmol; molCounter++)
        {
            std::copy_n(atomCharges.end() - molType->atoms.nr,
                        molType->atoms.nr,
                        std::back_inserter(atomCharges));
            std::copy_n(atomMasses.end() - molType->atoms.nr, molType->atoms.nr, std::back_inserter(atomMasses));
        }
    }
    /* Write the TNG data blocks. */
    tng_particle_data_block_add(tng,
                                TNG_TRAJ_PARTIAL_CHARGES,
                                "PARTIAL CHARGES",
                                datatype,
                                TNG_NON_TRAJECTORY_BLOCK,
                                1,
                                1,
                                1,
                                0,
                                mtop->natoms,
                                TNG_GZIP_COMPRESSION,
                                atomCharges.data());
    tng_particle_data_block_add(tng,
                                TNG_TRAJ_MASSES,
                                "ATOM MASSES",
                                datatype,
                                TNG_NON_TRAJECTORY_BLOCK,
                                1,
                                1,
                                1,
                                0,
                                mtop->natoms,
                                TNG_GZIP_COMPRESSION,
                                atomMasses.data());
}

/*! \libinternal \brief Compute greatest common divisor of n1 and n2
 * if they are positive.
 *
 * If only one of n1 and n2 is positive, then return it.
 * If neither n1 or n2 is positive, then return -1. */
static int greatest_common_divisor_if_positive(int n1, int n2)
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
    return std::gcd(n1, n2);
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
static void tng_set_frames_per_frame_set(gmx_tng_trajectory_t gmx_tng,
                                         const gmx_bool       bUseLossyCompression,
                                         const t_inputrec*    ir)
{
    int              gcd = -1;
    tng_trajectory_t tng = gmx_tng->tng;

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
static void set_writing_intervals(gmx_tng_trajectory_t gmx_tng,
                                  const gmx_bool       bUseLossyCompression,
                                  const t_inputrec*    ir)
{
    tng_trajectory_t tng = gmx_tng->tng;

    /* Define pointers to specific writing functions depending on if we
     * write float or double data */
    typedef tng_function_status (*set_writing_interval_func_pointer)(
            tng_trajectory_t, const int64_t, const int64_t, const int64_t, const char*, const char, const char);
#    if GMX_DOUBLE
    set_writing_interval_func_pointer set_writing_interval = tng_util_generic_write_interval_double_set;
#    else
    set_writing_interval_func_pointer set_writing_interval = tng_util_generic_write_interval_set;
#    endif
    int  xout, vout, fout;
    int  gcd = -1, lowest = -1;
    char compression;

    tng_set_frames_per_frame_set(gmx_tng, bUseLossyCompression, ir);

    if (bUseLossyCompression)
    {
        xout = ir->nstxout_compressed;

        /* If there is no uncompressed coordinate output write forces
           and velocities to the compressed tng file. */
        if (ir->nstxout)
        {
            vout = 0;
            fout = 0;
        }
        else
        {
            vout = ir->nstvout;
            fout = ir->nstfout;
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
        set_writing_interval(
                tng, xout, 3, TNG_TRAJ_POSITIONS, "POSITIONS", TNG_PARTICLE_BLOCK_DATA, compression);
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
        set_writing_interval(
                tng, vout, 3, TNG_TRAJ_VELOCITIES, "VELOCITIES", TNG_PARTICLE_BLOCK_DATA, compression);

        gcd = greatest_common_divisor_if_positive(gcd, vout);
        if (lowest < 0 || vout < lowest)
        {
            lowest = vout;
        }
    }
    if (fout)
    {
        set_writing_interval(
                tng, fout, 3, TNG_TRAJ_FORCES, "FORCES", TNG_PARTICLE_BLOCK_DATA, TNG_GZIP_COMPRESSION);

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
        set_writing_interval(
                tng, gcd, 1, TNG_GMX_LAMBDA, "LAMBDAS", TNG_NON_PARTICLE_BLOCK_DATA, TNG_GZIP_COMPRESSION);

        set_writing_interval(
                tng, gcd, 9, TNG_TRAJ_BOX_SHAPE, "BOX SHAPE", TNG_NON_PARTICLE_BLOCK_DATA, TNG_GZIP_COMPRESSION);
        gmx_tng->lambdaOutputInterval = gcd;
        gmx_tng->boxOutputInterval    = gcd;
        if (gcd < lowest / 10)
        {
            gmx_warning(
                    "The lowest common denominator of trajectory output is "
                    "every %d step(s), whereas the shortest output interval "
                    "is every %d steps.",
                    gcd,
                    lowest);
        }
    }
}
#endif

void gmx_tng_prepare_md_writing(gmx_tng_trajectory_t gmx_tng, const gmx_mtop_t* mtop, const t_inputrec* ir)
{
#if GMX_USE_TNG
    gmx_tng_add_mtop(gmx_tng, mtop);
    set_writing_intervals(gmx_tng, FALSE, ir);
    tng_time_per_frame_set(gmx_tng->tng, ir->delta_t * gmx::c_pico);
    gmx_tng->timePerFrameIsSet = true;
#else
    GMX_UNUSED_VALUE(gmx_tng);
    GMX_UNUSED_VALUE(mtop);
    GMX_UNUSED_VALUE(ir);
#endif
}

#if GMX_USE_TNG
/* Check if all atoms in the molecule system are selected
 * by a selection group of type specified by gtype. */
static gmx_bool all_atoms_selected(const gmx_mtop_t* mtop, const SimulationAtomGroupType gtype)
{
    /* Iterate through all atoms in the molecule system and
     * check if they belong to a selection group of the
     * requested type. */
    int i = 0;
    for (const gmx_molblock_t& molBlock : mtop->molblock)
    {
        const gmx_moltype_t& molType = mtop->moltype[molBlock.type];
        const t_atoms&       atoms   = molType.atoms;

        for (int j = 0; j < molBlock.nmol; j++)
        {
            for (int atomIndex = 0; atomIndex < atoms.nr; atomIndex++, i++)
            {
                if (getGroupType(mtop->groups, gtype, i) != 0)
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
static void add_selection_groups(gmx_tng_trajectory_t gmx_tng, const gmx_mtop_t* mtop)
{
    const t_atoms*   atoms;
    const t_atom*    at;
    const t_resinfo* resInfo;
    int              nameIndex;
    int              atom_offset = 0;
    tng_molecule_t   mol, iterMol;
    tng_chain_t      chain;
    tng_residue_t    res;
    tng_atom_t       atom;
    tng_bond_t       tngBond;
    int64_t          nMols;
    char*            groupName;
    tng_trajectory_t tng = gmx_tng->tng;

    /* TODO: When the TNG molecules block is more flexible TNG selection
     * groups should not need all atoms specified. It should be possible
     * just to specify what molecules are selected (and/or which atoms in
     * the molecule) and how many (if applicable). */

    /* If no atoms are selected we do not need to create a
     * TNG selection group molecule. */
    if (mtop->groups.numberOfGroupNumbers(SimulationAtomGroupType::CompressedPositionOutput) == 0)
    {
        return;
    }

    /* If all atoms are selected we do not have to create a selection
     * group molecule in the TNG molecule system. */
    if (all_atoms_selected(mtop, SimulationAtomGroupType::CompressedPositionOutput))
    {
        return;
    }

    /* The name of the TNG molecule containing the selection group is the
     * same as the name of the selection group. */
    nameIndex = mtop->groups.groups[SimulationAtomGroupType::CompressedPositionOutput][0];
    groupName = *mtop->groups.groupNames[nameIndex];

    tng_molecule_alloc(tng, &mol);
    tng_molecule_name_set(tng, mol, groupName);
    tng_molecule_chain_add(tng, mol, "", &chain);
    int i = 0;
    for (const gmx_molblock_t& molBlock : mtop->molblock)
    {
        const gmx_moltype_t& molType = mtop->moltype[molBlock.type];

        atoms = &molType.atoms;

        for (int j = 0; j < molBlock.nmol; j++)
        {
            bool bAtomsAdded = FALSE;
            for (int atomIndex = 0; atomIndex < atoms->nr; atomIndex++, i++)
            {
                char* res_name;
                int   res_id;

                if (getGroupType(mtop->groups, SimulationAtomGroupType::CompressedPositionOutput, i) != 0)
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
                    res_name = const_cast<char*>("");
                    res_id   = 0;
                }
                if (tng_chain_residue_find(tng, chain, res_name, res_id, &res) != TNG_SUCCESS)
                {
                    /* Since there is ONE chain for selection groups do not keep the
                     * original residue IDs - otherwise there might be conflicts. */
                    tng_chain_residue_add(tng, chain, res_name, &res);
                }
                tng_residue_atom_w_id_add(tng,
                                          res,
                                          *(atoms->atomname[atomIndex]),
                                          *(atoms->atomtype[atomIndex]),
                                          atom_offset + atomIndex,
                                          &atom);
                bAtomsAdded = TRUE;
            }
            /* Add bonds. */
            if (bAtomsAdded)
            {
                for (int k = 0; k < F_NRE; k++)
                {
                    if (IS_CHEMBOND(k))
                    {
                        const InteractionList& ilist = molType.ilist[k];
                        for (int l = 1; l < ilist.size(); l += 3)
                        {
                            int atom1, atom2;
                            atom1 = ilist.iatoms[l] + atom_offset;
                            atom2 = ilist.iatoms[l + 1] + atom_offset;
                            if (getGroupType(mtop->groups, SimulationAtomGroupType::CompressedPositionOutput, atom1)
                                        == 0
                                && getGroupType(mtop->groups, SimulationAtomGroupType::CompressedPositionOutput, atom2)
                                           == 0)
                            {
                                tng_molecule_bond_add(
                                        tng, mol, ilist.iatoms[l], ilist.iatoms[l + 1], &tngBond);
                            }
                        }
                    }
                }
                /* Settle is described using three atoms */
                const InteractionList& ilist = molType.ilist[F_SETTLE];
                for (int l = 1; l < ilist.size(); l += 4)
                {
                    int atom1, atom2, atom3;
                    atom1 = ilist.iatoms[l] + atom_offset;
                    atom2 = ilist.iatoms[l + 1] + atom_offset;
                    atom3 = ilist.iatoms[l + 2] + atom_offset;
                    if (getGroupType(mtop->groups, SimulationAtomGroupType::CompressedPositionOutput, atom1)
                        == 0)
                    {
                        if (getGroupType(mtop->groups, SimulationAtomGroupType::CompressedPositionOutput, atom2)
                            == 0)
                        {
                            tng_molecule_bond_add(tng, mol, atom1, atom2, &tngBond);
                        }
                        if (getGroupType(mtop->groups, SimulationAtomGroupType::CompressedPositionOutput, atom3)
                            == 0)
                        {
                            tng_molecule_bond_add(tng, mol, atom1, atom3, &tngBond);
                        }
                    }
                }
            }
            atom_offset += atoms->nr;
        }
    }
    tng_molecule_existing_add(tng, &mol);
    tng_molecule_cnt_set(tng, mol, 1);
    tng_num_molecule_types_get(tng, &nMols);
    for (int64_t k = 0; k < nMols; k++)
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

void gmx_tng_set_compression_precision(gmx_tng_trajectory_t gmx_tng, real prec)
{
#if GMX_USE_TNG
    tng_compression_precision_set(gmx_tng->tng, prec);
#else
    GMX_UNUSED_VALUE(gmx_tng);
    GMX_UNUSED_VALUE(prec);
#endif
}

void gmx_tng_prepare_low_prec_writing(gmx_tng_trajectory_t gmx_tng, const gmx_mtop_t* mtop, const t_inputrec* ir)
{
#if GMX_USE_TNG
    gmx_tng_add_mtop(gmx_tng, mtop);
    add_selection_groups(gmx_tng, mtop);
    set_writing_intervals(gmx_tng, TRUE, ir);
    tng_time_per_frame_set(gmx_tng->tng, ir->delta_t * gmx::c_pico);
    gmx_tng->timePerFrameIsSet = true;
    gmx_tng_set_compression_precision(gmx_tng, ir->x_compression_precision);
#else
    GMX_UNUSED_VALUE(gmx_tng);
    GMX_UNUSED_VALUE(mtop);
    GMX_UNUSED_VALUE(ir);
#endif
}

void gmx_fwrite_tng(gmx_tng_trajectory_t gmx_tng,
                    const gmx_bool       bUseLossyCompression,
                    int64_t              step,
                    real                 elapsedPicoSeconds,
                    real                 lambda,
                    const rvec*          box,
                    int                  nAtoms,
                    const rvec*          x,
                    const rvec*          v,
                    const rvec*          f)
{
#if GMX_USE_TNG
    typedef tng_function_status (*write_data_func_pointer)(tng_trajectory_t,
                                                           const int64_t,
                                                           const double,
                                                           const real*,
                                                           const int64_t,
                                                           const int64_t,
                                                           const char*,
                                                           const char,
                                                           const char);
#    if GMX_DOUBLE
    static write_data_func_pointer write_data = tng_util_generic_with_time_double_write;
#    else
    static write_data_func_pointer    write_data           = tng_util_generic_with_time_write;
#    endif
    double  elapsedSeconds = elapsedPicoSeconds * gmx::c_pico;
    int64_t nParticles;
    char    compression;


    if (!gmx_tng)
    {
        /* This function might get called when the type of the
           compressed trajectory is actually XTC. So we exit and move
           on. */
        return;
    }
    tng_trajectory_t tng = gmx_tng->tng;

    // While the GROMACS interface to this routine specifies 'step', TNG itself
    // only uses 'frame index' internally, although it suggests that it's a good
    // idea to let this represent the step rather than just frame numbers.
    //
    // However, the way the GROMACS interface works, there are cases where
    // the step is simply not set, so all frames rather have step=0.
    // If we call the lower-level TNG routines with the same frame index
    // (which is set from the step) multiple times, things will go very
    // wrong (overwritten data).
    //
    // To avoid this, we require the step value to always be larger than the
    // previous value, and if this is not true we simply set it to a value
    // one unit larger than the previous step.
    //
    // Note: We do allow the user to specify any crazy value the want for the
    // first step. If it is "not set", GROMACS will have used the value 0.
    if (gmx_tng->lastStepDataIsValid && step <= gmx_tng->lastStep)
    {
        step = gmx_tng->lastStep + 1;
    }

    // Now that we have fixed the potentially incorrect step, we can also set
    // the time per frame if it was not already done, and we have
    // step/time information from the last frame so it is possible to calculate it.
    // Note that we do allow the time to be the same for all steps, or even
    // decreasing. In that case we will simply say the time per step is 0
    // or negative. A bit strange, but a correct representation of the data :-)
    if (!gmx_tng->timePerFrameIsSet && gmx_tng->lastTimeDataIsValid && gmx_tng->lastStepDataIsValid)
    {
        double       deltaTime = elapsedSeconds - gmx_tng->lastTime;
        std::int64_t deltaStep = step - gmx_tng->lastStep;
        tng_time_per_frame_set(tng, deltaTime / deltaStep);
        gmx_tng->timePerFrameIsSet = true;
    }

    tng_num_particles_get(tng, &nParticles);
    if (nAtoms != static_cast<int>(nParticles))
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

        if (write_data(tng,
                       step,
                       elapsedSeconds,
                       reinterpret_cast<const real*>(x),
                       3,
                       TNG_TRAJ_POSITIONS,
                       "POSITIONS",
                       TNG_PARTICLE_BLOCK_DATA,
                       compression)
            != TNG_SUCCESS)
        {
            gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
        }
    }

    if (v)
    {
        if (write_data(tng,
                       step,
                       elapsedSeconds,
                       reinterpret_cast<const real*>(v),
                       3,
                       TNG_TRAJ_VELOCITIES,
                       "VELOCITIES",
                       TNG_PARTICLE_BLOCK_DATA,
                       compression)
            != TNG_SUCCESS)
        {
            gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
        }
    }

    if (f)
    {
        /* TNG-MF1 compression only compresses positions and velocities. Use lossless
         * compression for forces regardless of output mode */
        if (write_data(tng,
                       step,
                       elapsedSeconds,
                       reinterpret_cast<const real*>(f),
                       3,
                       TNG_TRAJ_FORCES,
                       "FORCES",
                       TNG_PARTICLE_BLOCK_DATA,
                       TNG_GZIP_COMPRESSION)
            != TNG_SUCCESS)
        {
            gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
        }
    }

    if (box)
    {
        /* TNG-MF1 compression only compresses positions and velocities. Use lossless
         * compression for box shape regardless of output mode */
        if (write_data(tng,
                       step,
                       elapsedSeconds,
                       reinterpret_cast<const real*>(box),
                       9,
                       TNG_TRAJ_BOX_SHAPE,
                       "BOX SHAPE",
                       TNG_NON_PARTICLE_BLOCK_DATA,
                       TNG_GZIP_COMPRESSION)
            != TNG_SUCCESS)
        {
            gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
        }
    }

    if (lambda >= 0)
    {
        /* TNG-MF1 compression only compresses positions and velocities. Use lossless
         * compression for lambda regardless of output mode */
        if (write_data(tng,
                       step,
                       elapsedSeconds,
                       reinterpret_cast<const real*>(&lambda),
                       1,
                       TNG_GMX_LAMBDA,
                       "LAMBDAS",
                       TNG_NON_PARTICLE_BLOCK_DATA,
                       TNG_GZIP_COMPRESSION)
            != TNG_SUCCESS)
        {
            gmx_file("Cannot write TNG trajectory frame; maybe you are out of disk space?");
        }
    }

    // Update the data in the wrapper

    gmx_tng->lastStepDataIsValid = true;
    gmx_tng->lastStep            = step;
    gmx_tng->lastTimeDataIsValid = true;
    gmx_tng->lastTime            = elapsedSeconds;
#else
    GMX_UNUSED_VALUE(gmx_tng);
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

void fflush_tng(gmx_tng_trajectory_t gmx_tng)
{
#if GMX_USE_TNG
    if (!gmx_tng)
    {
        return;
    }
    tng_frame_set_premature_write(gmx_tng->tng, TNG_USE_HASH);
#else
    GMX_UNUSED_VALUE(gmx_tng);
#endif
}

float gmx_tng_get_time_of_final_frame(gmx_tng_trajectory_t gmx_tng)
{
#if GMX_USE_TNG
    int64_t          nFrames;
    double           time;
    float            fTime;
    tng_trajectory_t tng = gmx_tng->tng;

    tng_num_frames_get(tng, &nFrames);
    tng_util_time_of_frame_get(tng, nFrames - 1, &time);

    fTime = time / gmx::c_pico;
    return fTime;
#else
    GMX_UNUSED_VALUE(gmx_tng);
    return -1.0;
#endif
}

void gmx_prepare_tng_writing(const std::filesystem::path& filename,
                             char                         mode,
                             gmx_tng_trajectory_t*        gmx_tng_input,
                             gmx_tng_trajectory_t*        gmx_tng_output,
                             int                          nAtoms,
                             const gmx_mtop_t*            mtop,
                             gmx::ArrayRef<const int>     index,
                             const char*                  indexGroupName)
{
#if GMX_USE_TNG
    tng_trajectory_t* input = (gmx_tng_input && *gmx_tng_input) ? &(*gmx_tng_input)->tng : nullptr;
    /* FIXME after 5.0: Currently only standard block types are read */
    const int      defaultNumIds              = 5;
    static int64_t fallbackIds[defaultNumIds] = {
        TNG_TRAJ_BOX_SHAPE, TNG_TRAJ_POSITIONS, TNG_TRAJ_VELOCITIES, TNG_TRAJ_FORCES, TNG_GMX_LAMBDA
    };
    static char fallbackNames[defaultNumIds][32] = {
        "BOX SHAPE", "POSITIONS", "VELOCITIES", "FORCES", "LAMBDAS"
    };

    typedef tng_function_status (*set_writing_interval_func_pointer)(
            tng_trajectory_t, const int64_t, const int64_t, const int64_t, const char*, const char, const char);
#    if GMX_DOUBLE
    set_writing_interval_func_pointer set_writing_interval = tng_util_generic_write_interval_double_set;
#    else
    set_writing_interval_func_pointer set_writing_interval = tng_util_generic_write_interval_set;
#    endif

    gmx_tng_open(filename, mode, gmx_tng_output);
    tng_trajectory_t* output = &(*gmx_tng_output)->tng;

    /* Do we have an input file in TNG format? If so, then there's
       more data we can copy over, rather than having to improvise. */
    if (gmx_tng_input && *gmx_tng_input)
    {
        /* Set parameters (compression, time per frame, molecule
         * information, number of frames per frame set and writing
         * intervals of positions, box shape and lambdas) of the
         * output tng container based on their respective values in
         * the input tng container */
        double  time, compression_precision;
        int64_t n_frames_per_frame_set, interval = -1;

        tng_compression_precision_get(*input, &compression_precision);
        tng_compression_precision_set(*output, compression_precision);
        // TODO make this configurable in a future version
        char compression_type = TNG_TNG_COMPRESSION;

        tng_molecule_system_copy(*input, *output);
        /* This should be done before adding the trajectory data blocks below since
         * adding them forces the TNG header, including the molecule data, to be written. */
        if ((!index.empty()) && nAtoms > 0)
        {
            gmx_tng_setup_atom_subgroup(*gmx_tng_output, index, indexGroupName);
        }

        tng_time_per_frame_get(*input, &time);
        /* Only write the time per frame if it was written (and valid). E.g., single
         * frame files do not usually contain any time per frame information. */
        if (time >= 0)
        {
            (*gmx_tng_input)->timePerFrameIsSet = true;
            tng_time_per_frame_set(*output, time);
            // Since we have copied the value from the input TNG we should not change it again
            (*gmx_tng_output)->timePerFrameIsSet = true;
        }

        tng_num_frames_per_frame_set_get(*input, &n_frames_per_frame_set);
        tng_num_frames_per_frame_set_set(*output, n_frames_per_frame_set);

        for (int i = 0; i < defaultNumIds; i++)
        {
            if (tng_data_get_stride_length(*input, fallbackIds[i], -1, &interval) == TNG_SUCCESS)
            {
                switch (fallbackIds[i])
                {
                    case TNG_TRAJ_POSITIONS:
                    case TNG_TRAJ_VELOCITIES:
                        set_writing_interval(*output,
                                             interval,
                                             3,
                                             fallbackIds[i],
                                             fallbackNames[i],
                                             TNG_PARTICLE_BLOCK_DATA,
                                             compression_type);
                        break;
                    case TNG_TRAJ_FORCES:
                        set_writing_interval(*output,
                                             interval,
                                             3,
                                             fallbackIds[i],
                                             fallbackNames[i],
                                             TNG_PARTICLE_BLOCK_DATA,
                                             TNG_GZIP_COMPRESSION);
                        break;
                    case TNG_TRAJ_BOX_SHAPE:
                        set_writing_interval(*output,
                                             interval,
                                             9,
                                             fallbackIds[i],
                                             fallbackNames[i],
                                             TNG_NON_PARTICLE_BLOCK_DATA,
                                             TNG_GZIP_COMPRESSION);
                        (*gmx_tng_output)->boxOutputInterval = interval;
                        break;
                    case TNG_GMX_LAMBDA:
                        set_writing_interval(*output,
                                             interval,
                                             1,
                                             fallbackIds[i],
                                             fallbackNames[i],
                                             TNG_NON_PARTICLE_BLOCK_DATA,
                                             TNG_GZIP_COMPRESSION);
                        (*gmx_tng_output)->lambdaOutputInterval = interval;
                        break;
                    default: continue;
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
        gmx_tng_add_mtop(*gmx_tng_output, mtop);
        if ((!index.empty()) && nAtoms > 0)
        {
            gmx_tng_setup_atom_subgroup(*gmx_tng_output, index, indexGroupName);
        }
        tng_num_frames_per_frame_set_set(*output, 1);
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
    GMX_UNUSED_VALUE(gmx_tng_input);
    GMX_UNUSED_VALUE(gmx_tng_output);
    GMX_UNUSED_VALUE(nAtoms);
    GMX_UNUSED_VALUE(mtop);
    GMX_UNUSED_VALUE(index);
    GMX_UNUSED_VALUE(indexGroupName);
#endif
}

void gmx_write_tng_from_trxframe(gmx_tng_trajectory_t gmx_tng_output, const t_trxframe* frame, int natoms)
{
#if GMX_USE_TNG
    if (natoms < 0)
    {
        natoms = frame->natoms;
    }
    gmx_fwrite_tng(
            gmx_tng_output, TRUE, frame->step, frame->time, 0, frame->box, natoms, frame->x, frame->v, frame->f);
#else
    GMX_UNUSED_VALUE(gmx_tng_output);
    GMX_UNUSED_VALUE(frame);
    GMX_UNUSED_VALUE(natoms);
#endif
}

namespace
{

#if GMX_USE_TNG
void convert_array_to_real_array(void*       from,
                                 real*       to,
                                 const float fact,
                                 const int   nAtoms,
                                 const int   nValues,
                                 const char  datatype)
{
    int i, j;

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
                            to[i * nValues + j] = reinterpret_cast<float*>(from)[i * nValues + j] * fact;
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
                        to[i * nValues + j] = reinterpret_cast<float*>(from)[i * nValues + j] * fact;
                    }
                }
            }
            break;
        case TNG_INT_DATA:
            for (i = 0; i < nAtoms; i++)
            {
                for (j = 0; j < nValues; j++)
                {
                    to[i * nValues + j] = reinterpret_cast<int64_t*>(from)[i * nValues + j] * fact;
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
                            to[i * nValues + j] = reinterpret_cast<double*>(from)[i * nValues + j] * fact;
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
                        to[i * nValues + j] = reinterpret_cast<double*>(from)[i * nValues + j] * fact;
                    }
                }
            }
            break;
        default: gmx_incons("Illegal datatype when converting values to a real array!");
    }
}

real getDistanceScaleFactor(gmx_tng_trajectory_t in)
{
    int64_t exp = -1;
    real    distanceScaleFactor;

    // TODO Hopefully, TNG 2.0 will do this kind of thing for us
    tng_distance_unit_exponential_get(in->tng, &exp);

    // GROMACS expects distances in nm
    switch (exp)
    {
        case 9: distanceScaleFactor = gmx::c_nano / gmx::c_nano; break;
        case 10: distanceScaleFactor = gmx::c_nano / gmx::c_angstrom; break;
        default: distanceScaleFactor = std::pow(10.0, exp + 9.0);
    }

    return distanceScaleFactor;
}
#endif

} // namespace

void gmx_tng_setup_atom_subgroup(gmx_tng_trajectory_t gmx_tng, gmx::ArrayRef<const int> ind, const char* name)
{
#if GMX_USE_TNG
    int64_t             nAtoms, cnt, nMols;
    tng_molecule_t      mol, iterMol;
    tng_chain_t         chain;
    tng_residue_t       res;
    tng_atom_t          atom;
    tng_function_status stat;
    tng_trajectory_t    tng = gmx_tng->tng;

    tng_num_particles_get(tng, &nAtoms);

    if (nAtoms == ind.ssize())
    {
        return;
    }

    /* Check if there is one molecule type matching the selection name. */
    stat = tng_molecule_find(tng, name, -1, &mol);
    if (stat == TNG_SUCCESS)
    {
        tng_molecule_num_atoms_get(tng, mol, &nAtoms);
        tng_molecule_cnt_get(tng, mol, &cnt);
        if (nAtoms * cnt == ind.ssize())
        {
            stat = TNG_SUCCESS;
        }
        else
        {
            /* Since the molecule that matched the name of the selection did
             * not match the number of atoms in the selection set the number
             * of molecules of that type to 0. Below a new molecule type will
             * be added matching that of the selection. */
            tng_molecule_cnt_set(tng, mol, 0);
            stat = TNG_FAILURE;
        }
    }
    if (stat == TNG_FAILURE)
    {
        /* The indexed atoms are added to one separate molecule. */
        tng_molecule_alloc(tng, &mol);
        tng_molecule_name_set(tng, mol, name);
        tng_molecule_chain_add(tng, mol, "", &chain);

        for (gmx::Index i = 0; i < ind.ssize(); i++)
        {
            char temp_name[256], temp_type[256];

            /* Try to retrieve the residue name of the atom */
            stat = tng_residue_name_of_particle_nr_get(tng, ind[i], temp_name, 256);
            if (stat != TNG_SUCCESS)
            {
                temp_name[0] = '\0';
            }
            /* Check if the molecule of the selection already contains this residue */
            if (tng_chain_residue_find(tng, chain, temp_name, -1, &res) != TNG_SUCCESS)
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
        /* Set the count of the add molecule containing the selected atoms to 1 */
        tng_molecule_cnt_set(tng, mol, 1);
    }
    /* Set the count of other molecule types to 0 */
    tng_num_molecule_types_get(tng, &nMols);
    for (int64_t k = 0; k < nMols; k++)
    {
        tng_molecule_of_index_get(tng, k, &iterMol);
        if (iterMol == mol)
        {
            continue;
        }
        tng_molecule_cnt_set(tng, iterMol, 0);
    }
#else
    GMX_UNUSED_VALUE(gmx_tng);
    GMX_UNUSED_VALUE(ind);
    GMX_UNUSED_VALUE(name);
#endif
}

/* TODO: If/when TNG acquires the ability to copy data blocks without
 * uncompressing them, then this implementation should be reconsidered.
 * Ideally, gmx trjconv -f a.tng -o b.tng -b 10 -e 20 would be fast
 * and lose no information. */
gmx_bool gmx_read_next_tng_frame(gmx_tng_trajectory_t gmx_tng_input,
                                 t_trxframe*          fr,
                                 int64_t*             requestedIds,
                                 int                  numRequestedIds)
{
#if GMX_USE_TNG
    tng_trajectory_t    input = gmx_tng_input->tng;
    gmx_bool            bOK   = TRUE;
    tng_function_status stat;
    int64_t             numberOfAtoms = -1, frameNumber = -1;
    int64_t             nBlocks, blockId, *blockIds = nullptr, codecId;
    char                datatype  = -1;
    void*               values    = nullptr;
    double              frameTime = -1.0;
    int                 size, blockDependency;
    double              prec;
    const int           defaultNumIds                       = 5;
    static int64_t      fallbackRequestedIds[defaultNumIds] = {
        TNG_TRAJ_BOX_SHAPE, TNG_TRAJ_POSITIONS, TNG_TRAJ_VELOCITIES, TNG_TRAJ_FORCES, TNG_GMX_LAMBDA
    };


    fr->bStep   = FALSE;
    fr->bTime   = FALSE;
    fr->bLambda = FALSE;
    fr->bAtoms  = FALSE;
    fr->bPrec   = FALSE;
    fr->bX      = FALSE;
    fr->bV      = FALSE;
    fr->bF      = FALSE;
    fr->bBox    = FALSE;

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

    /* If the current (last read/written) TNG step is recorded use that. Otherwise retrieve it from the frame data. */
    std::int64_t step;
    if (gmx_tng_input->lastStepDataIsValid)
    {
        step = gmx_tng_input->lastStep;
    }
    else
    {
        step = fr->step;
    }

    bool nextFrameExists = gmx_get_tng_data_block_types_of_next_frame(
            gmx_tng_input, step, numRequestedIds, requestedIds, &frameNumber, &nBlocks, &blockIds);
    gmx::unique_cptr<int64_t, gmx::free_wrapper> blockIdsGuard(blockIds);
    if (!nextFrameExists)
    {
        return FALSE;
    }

    if (nBlocks == 0)
    {
        return FALSE;
    }

    for (int64_t i = 0; i < nBlocks; i++)
    {
        blockId = blockIds[i];
        tng_data_block_dependency_get(input, blockId, &blockDependency);
        if (blockDependency & TNG_PARTICLE_DEPENDENT)
        {
            stat = tng_util_particle_data_next_frame_read(
                    input, blockId, &values, &datatype, &frameNumber, &frameTime);
        }
        else
        {
            stat = tng_util_non_particle_data_next_frame_read(
                    input, blockId, &values, &datatype, &frameNumber, &frameTime);
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
                    case TNG_INT_DATA: size = sizeof(int64_t); break;
                    case TNG_FLOAT_DATA: size = sizeof(float); break;
                    case TNG_DOUBLE_DATA: size = sizeof(double); break;
                    default: gmx_incons("Illegal datatype of box shape values!");
                }
                for (int i = 0; i < DIM; i++)
                {
                    convert_array_to_real_array(reinterpret_cast<char*>(values) + size * i * DIM,
                                                reinterpret_cast<real*>(fr->box[i]),
                                                getDistanceScaleFactor(gmx_tng_input),
                                                1,
                                                DIM,
                                                datatype);
                }
                fr->bBox = TRUE;
                break;
            case TNG_TRAJ_POSITIONS:
                srenew(fr->x, fr->natoms);
                convert_array_to_real_array(values,
                                            reinterpret_cast<real*>(fr->x),
                                            getDistanceScaleFactor(gmx_tng_input),
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
                                            reinterpret_cast<real*>(fr->v),
                                            getDistanceScaleFactor(gmx_tng_input),
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
                                            reinterpret_cast<real*>(fr->f),
                                            getDistanceScaleFactor(gmx_tng_input),
                                            fr->natoms,
                                            DIM,
                                            datatype);
                fr->bF = TRUE;
                break;
            case TNG_GMX_LAMBDA:
                switch (datatype)
                {
                    case TNG_FLOAT_DATA: fr->lambda = *(reinterpret_cast<float*>(values)); break;
                    case TNG_DOUBLE_DATA: fr->lambda = *(reinterpret_cast<double*>(values)); break;
                    default: gmx_incons("Illegal datatype lambda value!");
                }
                fr->bLambda = TRUE;
                break;
            default:
                gmx_warning(
                        "Illegal block type! Currently GROMACS tools can only handle certain data "
                        "types. Skipping block.");
        }
        /* values does not have to be freed before reading next frame. It will
         * be reallocated if it is not NULL. */
    }

    fr->step  = frameNumber;
    fr->bStep = TRUE;

    // Convert the time to ps
    fr->time  = frameTime / gmx::c_pico;
    fr->bTime = (frameTime > 0);

    /* Update the data in the wrapper */
    gmx_tng_input->lastStepDataIsValid = true;
    gmx_tng_input->lastStep            = frameNumber;
    gmx_tng_input->lastTimeDataIsValid = true;
    gmx_tng_input->lastTime            = frameTime;

    // TODO This does not leak, but is not exception safe.
    /* values must be freed before leaving this function */
    sfree(values);

    return bOK;
#else
    GMX_UNUSED_VALUE(gmx_tng_input);
    GMX_UNUSED_VALUE(fr);
    GMX_UNUSED_VALUE(requestedIds);
    GMX_UNUSED_VALUE(numRequestedIds);
    return FALSE;
#endif
}

void gmx_print_tng_molecule_system(gmx_tng_trajectory_t gmx_tng_input, FILE* stream)
{
#if GMX_USE_TNG
    int64_t             nMolecules, nChains, nResidues, nAtoms, nFramesRead;
    int64_t             strideLength, nParticlesRead, nValuesPerFrameRead, *molCntList;
    tng_molecule_t      molecule;
    tng_chain_t         chain;
    tng_residue_t       residue;
    tng_atom_t          atom;
    tng_function_status stat;
    char                str[256];
    char                varNAtoms;
    char                datatype;
    void*               data = nullptr;
    std::vector<real>   atomCharges;
    std::vector<real>   atomMasses;
    tng_trajectory_t    input = gmx_tng_input->tng;

    tng_num_molecule_types_get(input, &nMolecules);
    tng_molecule_cnt_list_get(input, &molCntList);
    /* Can the number of particles change in the trajectory or is it constant? */
    tng_num_particles_variable_get(input, &varNAtoms);

    for (int64_t i = 0; i < nMolecules; i++)
    {
        tng_molecule_of_index_get(input, i, &molecule);
        tng_molecule_name_get(input, molecule, str, 256);
        if (varNAtoms == TNG_CONSTANT_N_ATOMS)
        {
            if (static_cast<int>(molCntList[i]) == 0)
            {
                continue;
            }
            fprintf(stream, "Molecule: %s, count: %d\n", str, static_cast<int>(molCntList[i]));
        }
        else
        {
            fprintf(stream, "Molecule: %s\n", str);
        }
        tng_molecule_num_chains_get(input, molecule, &nChains);
        if (nChains > 0)
        {
            for (int64_t j = 0; j < nChains; j++)
            {
                tng_molecule_chain_of_index_get(input, molecule, j, &chain);
                tng_chain_name_get(input, chain, str, 256);
                fprintf(stream, "\tChain: %s\n", str);
                tng_chain_num_residues_get(input, chain, &nResidues);
                for (int64_t k = 0; k < nResidues; k++)
                {
                    tng_chain_residue_of_index_get(input, chain, k, &residue);
                    tng_residue_name_get(input, residue, str, 256);
                    fprintf(stream, "\t\tResidue: %s\n", str);
                    tng_residue_num_atoms_get(input, residue, &nAtoms);
                    for (int64_t l = 0; l < nAtoms; l++)
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
                for (int64_t k = 0; k < nResidues; k++)
                {
                    tng_molecule_residue_of_index_get(input, molecule, k, &residue);
                    tng_residue_name_get(input, residue, str, 256);
                    fprintf(stream, "\t\tResidue: %s\n", str);
                    tng_residue_num_atoms_get(input, residue, &nAtoms);
                    for (int64_t l = 0; l < nAtoms; l++)
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
                for (int64_t l = 0; l < nAtoms; l++)
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
    stat = tng_particle_data_vector_get(input,
                                        TNG_TRAJ_PARTIAL_CHARGES,
                                        &data,
                                        &nFramesRead,
                                        &strideLength,
                                        &nParticlesRead,
                                        &nValuesPerFrameRead,
                                        &datatype);
    if (stat == TNG_SUCCESS)
    {
        atomCharges.resize(nAtoms);
        convert_array_to_real_array(data, atomCharges.data(), 1, nAtoms, 1, datatype);

        fprintf(stream, "Atom Charges (%d):\n", int(nAtoms));
        for (int64_t i = 0; i < nAtoms; i += 10)
        {
            fprintf(stream, "Atom Charges [%8d-]=[", int(i));
            for (int64_t j = 0; (j < 10 && i + j < nAtoms); j++)
            {
                fprintf(stream, " %12.5e", atomCharges[i + j]);
            }
            fprintf(stream, "]\n");
        }
    }

    stat = tng_particle_data_vector_get(
            input, TNG_TRAJ_MASSES, &data, &nFramesRead, &strideLength, &nParticlesRead, &nValuesPerFrameRead, &datatype);
    if (stat == TNG_SUCCESS)
    {
        atomMasses.resize(nAtoms);
        convert_array_to_real_array(data, atomMasses.data(), 1, nAtoms, 1, datatype);

        fprintf(stream, "Atom Masses (%d):\n", int(nAtoms));
        for (int64_t i = 0; i < nAtoms; i += 10)
        {
            fprintf(stream, "Atom Masses [%8d-]=[", int(i));
            for (int64_t j = 0; (j < 10 && i + j < nAtoms); j++)
            {
                fprintf(stream, " %12.5e", atomMasses[i + j]);
            }
            fprintf(stream, "]\n");
        }
    }

    sfree(data);
#else
    GMX_UNUSED_VALUE(gmx_tng_input);
    GMX_UNUSED_VALUE(stream);
#endif
}

gmx_bool gmx_get_tng_data_block_types_of_next_frame(gmx_tng_trajectory_t gmx_tng_input,
                                                    int                  frame,
                                                    int                  nRequestedIds,
                                                    int64_t*             requestedIds,
                                                    int64_t*             nextFrame,
                                                    int64_t*             nBlocks,
                                                    int64_t**            blockIds)
{
#if GMX_USE_TNG
    tng_function_status stat;
    tng_trajectory_t    input = gmx_tng_input->tng;

    stat = tng_util_trajectory_next_frame_present_data_blocks_find(
            input, frame, nRequestedIds, requestedIds, nextFrame, nBlocks, blockIds);

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
    GMX_UNUSED_VALUE(gmx_tng_input);
    GMX_UNUSED_VALUE(frame);
    GMX_UNUSED_VALUE(nRequestedIds);
    GMX_UNUSED_VALUE(requestedIds);
    GMX_UNUSED_VALUE(nextFrame);
    GMX_UNUSED_VALUE(nBlocks);
    GMX_UNUSED_VALUE(blockIds);
    return FALSE;
#endif
}

gmx_bool gmx_get_tng_data_next_frame_of_block_type(gmx_tng_trajectory_t gmx_tng_input,
                                                   int64_t              blockId,
                                                   real**               values,
                                                   int64_t*             frameNumber,
                                                   double*              frameTime,
                                                   int64_t*             nValuesPerFrame,
                                                   int64_t*             nAtoms,
                                                   real*                prec,
                                                   char*                name,
                                                   int                  maxLen,
                                                   gmx_bool*            bOK)
{
#if GMX_USE_TNG
    tng_function_status stat;
    char                datatype = -1;
    int64_t             codecId;
    int                 blockDependency;
    void*               data = nullptr;
    double              localPrec;
    tng_trajectory_t    input = gmx_tng_input->tng;

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
        stat = tng_util_particle_data_next_frame_read(
                input, blockId, &data, &datatype, frameNumber, frameTime);
    }
    else
    {
        *nAtoms = 1; /* There are not actually any atoms, but it is used for
                        allocating memory */
        stat = tng_util_non_particle_data_next_frame_read(
                input, blockId, &data, &datatype, frameNumber, frameTime);
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
    convert_array_to_real_array(
            data, *values, getDistanceScaleFactor(gmx_tng_input), *nAtoms, *nValuesPerFrame, datatype);

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
    GMX_UNUSED_VALUE(gmx_tng_input);
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

int gmx_tng_get_box_output_interval(gmx_tng_trajectory_t gmx_tng)
{
#if GMX_USE_TNG
    return gmx_tng->boxOutputInterval;
#else
    GMX_UNUSED_VALUE(gmx_tng);
    return -1;
#endif
}

int gmx_tng_get_lambda_output_interval(gmx_tng_trajectory_t gmx_tng)
{
#if GMX_USE_TNG
    return gmx_tng->lambdaOutputInterval;
#else
    GMX_UNUSED_VALUE(gmx_tng);
    return -1;
#endif
}
