/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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

/*! \internal \file
 *
 * \brief This file defines functions used by the domdec module
 * while managing the construction, use and error checking for
 * topologies local to a DD rank.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "gromacs/domdec/computemultibodycutoffs.h"

#include <cmath>

#include <memory>
#include <vector>

#include "gromacs/domdec/options.h"
#include "gromacs/domdec/reversetopology.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/logger.h"

using gmx::ArrayRef;
using gmx::DDBondedChecking;
using gmx::RVec;

typedef struct
{
    real r2;
    int  ftype;
    int  a1;
    int  a2;
} bonded_distance_t;

/*! \brief Compare distance^2 \p r2 against the distance in \p bd and if larger store it along with \p ftype and atom indices \p a1 and \p a2 */
static void update_max_bonded_distance(real r2, int ftype, int a1, int a2, bonded_distance_t* bd)
{
    if (r2 > bd->r2)
    {
        bd->r2    = r2;
        bd->ftype = ftype;
        bd->a1    = a1;
        bd->a2    = a2;
    }
}

//! Returns the squared distance between two atoms, with or without PBC correction
template<bool usePbc>
static real inline distanceSquared(const t_pbc* pbc, const RVec& x1, const RVec& x2)
{
    if constexpr (usePbc)
    {
        rvec dx;
        pbc_dx(pbc, x1, x2, dx);
        return norm2(dx);
    }
    else
    {
        return distance2(x1, x2);

        GMX_UNUSED_VALUE(pbc);
    }
}

/*! \brief Set the distance, function type and atom indices for the longest distance between atoms of molecule type \p molt for two-body and multi-body bonded interactions */
template<bool usePbc>
static void bonded_cg_distance_mol(const gmx_moltype_t*   molt,
                                   const DDBondedChecking ddBondedChecking,
                                   gmx_bool               bExcl,
                                   const t_pbc*           pbc,
                                   ArrayRef<const RVec>   x,
                                   bonded_distance_t*     bd_2b,
                                   bonded_distance_t*     bd_mb)
{
    const ReverseTopOptions rtOptions(ddBondedChecking);

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (dd_check_ftype(ftype, rtOptions))
        {
            const auto& il   = molt->ilist[ftype];
            int         nral = NRAL(ftype);
            if (nral > 1)
            {
                for (int i = 0; i < il.size(); i += 1 + nral)
                {
                    for (int ai = 0; ai < nral; ai++)
                    {
                        int atomI = il.iatoms[i + 1 + ai];
                        for (int aj = ai + 1; aj < nral; aj++)
                        {
                            int atomJ = il.iatoms[i + 1 + aj];
                            if (atomI != atomJ)
                            {
                                const real rij2 = distanceSquared<usePbc>(pbc, x[atomI], x[atomJ]);

                                update_max_bonded_distance(
                                        rij2, ftype, atomI, atomJ, (nral == 2) ? bd_2b : bd_mb);
                            }
                        }
                    }
                }
            }
        }
    }
    if (bExcl)
    {
        const auto& excls = molt->excls;
        for (gmx::Index ai = 0; ai < excls.ssize(); ai++)
        {
            for (const int aj : excls[ai])
            {
                if (ai != aj)
                {
                    const real rij2 = distanceSquared<usePbc>(pbc, x[ai], x[aj]);

                    /* There is no function type for exclusions, use -1 */
                    update_max_bonded_distance(rij2, -1, ai, aj, bd_2b);
                }
            }
        }
    }
}

/*! \brief Set the distance, function type and atom indices for the longest atom distance involved in intermolecular interactions for two-body and multi-body bonded interactions */
static void bonded_distance_intermol(const InteractionLists& ilists_intermol,
                                     const DDBondedChecking  ddBondedChecking,
                                     ArrayRef<const RVec>    x,
                                     PbcType                 pbcType,
                                     const matrix            box,
                                     bonded_distance_t*      bd_2b,
                                     bonded_distance_t*      bd_mb)
{
    t_pbc pbc;

    set_pbc(&pbc, pbcType, box);

    const ReverseTopOptions rtOptions(ddBondedChecking);

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (dd_check_ftype(ftype, rtOptions))
        {
            const auto& il   = ilists_intermol[ftype];
            int         nral = NRAL(ftype);

            /* No nral>1 check here, since intermol interactions always
             * have nral>=2 (and the code is also correct for nral=1).
             */
            for (int i = 0; i < il.size(); i += 1 + nral)
            {
                for (int ai = 0; ai < nral; ai++)
                {
                    int atom_i = il.iatoms[i + 1 + ai];

                    for (int aj = ai + 1; aj < nral; aj++)
                    {
                        rvec dx;

                        int atom_j = il.iatoms[i + 1 + aj];

                        pbc_dx(&pbc, x[atom_i], x[atom_j], dx);

                        const real rij2 = norm2(dx);

                        update_max_bonded_distance(rij2, ftype, atom_i, atom_j, (nral == 2) ? bd_2b : bd_mb);
                    }
                }
            }
        }
    }
}

//! Returns whether \p molt has at least one virtual site
static bool moltypeHasVsite(const gmx_moltype_t& molt)
{
    bool hasVsite = false;
    for (int i = 0; i < F_NRE; i++)
    {
        if ((interaction_function[i].flags & IF_VSITE) && !molt.ilist[i].empty())
        {
            hasVsite = true;
        }
    }

    return hasVsite;
}

//! Returns coordinates not broken over PBC for a molecule
static void getWholeMoleculeCoordinates(const gmx_moltype_t*  molt,
                                        const gmx_ffparams_t* ffparams,
                                        PbcType               pbcType,
                                        t_graph*              graph,
                                        const matrix          box,
                                        ArrayRef<const RVec>  x,
                                        ArrayRef<RVec>        xs)
{
    if (pbcType != PbcType::No)
    {
        mk_mshift(nullptr, graph, pbcType, box, as_rvec_array(x.data()));

        shift_x(graph, box, as_rvec_array(x.data()), as_rvec_array(xs.data()));
        /* By doing an extra mk_mshift the molecules that are broken
         * because they were e.g. imported from another software
         * will be made whole again. Such are the healing powers
         * of GROMACS.
         */
        mk_mshift(nullptr, graph, pbcType, box, as_rvec_array(xs.data()));
    }
    else
    {
        /* We copy the coordinates so the original coordinates remain
         * unchanged, just to be 100% sure that we do not affect
         * binary reproducibility of simulations.
         */
        for (int i = 0; i < molt->atoms.nr; i++)
        {
            copy_rvec(x[i], xs[i]);
        }
    }

    if (moltypeHasVsite(*molt))
    {
        gmx::constructVirtualSites(xs, ffparams->iparams, molt->ilist);
    }
}

void dd_bonded_cg_distance(const gmx::MDLogger&   mdlog,
                           const gmx_mtop_t&      mtop,
                           const t_inputrec&      inputrec,
                           ArrayRef<const RVec>   x,
                           const matrix           box,
                           const DDBondedChecking ddBondedChecking,
                           real*                  r_2b,
                           real*                  r_mb)
{
    bonded_distance_t bd_2b = { 0, -1, -1, -1 };
    bonded_distance_t bd_mb = { 0, -1, -1, -1 };

    bool bExclRequired = inputrecExclForces(&inputrec);

    t_pbc pbc;
    if (inputrec.bPeriodicMols)
    {
        set_pbc(&pbc, inputrec.pbcType, box);
    }

    *r_2b         = 0;
    *r_mb         = 0;
    int at_offset = 0;
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        const gmx_moltype_t& molt = mtop.moltype[molb.type];
        if (molt.atoms.nr == 1 || molb.nmol == 0)
        {
            at_offset += molb.nmol * molt.atoms.nr;
        }
        else
        {
            t_graph graph;
            if (inputrec.pbcType != PbcType::No)
            {
                graph = mk_graph_moltype(molt);
            }

            std::vector<RVec> xs(molt.atoms.nr);
            for (int mol = 0; mol < molb.nmol; mol++)
            {
                ArrayRef<const RVec> xMolWithPbc = x.subArray(at_offset, molt.atoms.nr);

                bonded_distance_t bd_mol_2b = { 0, -1, -1, -1 };
                bonded_distance_t bd_mol_mb = { 0, -1, -1, -1 };

                if (!inputrec.bPeriodicMols)
                {
                    getWholeMoleculeCoordinates(
                            &molt, &mtop.ffparams, inputrec.pbcType, &graph, box, xMolWithPbc, xs);

                    bonded_cg_distance_mol<false>(
                            &molt, ddBondedChecking, bExclRequired, nullptr, xs, &bd_mol_2b, &bd_mol_mb);
                }
                else
                {
                    // Compute distances using PBC corrections
                    bonded_cg_distance_mol<true>(
                            &molt, ddBondedChecking, bExclRequired, &pbc, xMolWithPbc, &bd_mol_2b, &bd_mol_mb);
                }

                /* Process the mol data adding the atom index offset */
                update_max_bonded_distance(bd_mol_2b.r2,
                                           bd_mol_2b.ftype,
                                           at_offset + bd_mol_2b.a1,
                                           at_offset + bd_mol_2b.a2,
                                           &bd_2b);
                update_max_bonded_distance(bd_mol_mb.r2,
                                           bd_mol_mb.ftype,
                                           at_offset + bd_mol_mb.a1,
                                           at_offset + bd_mol_mb.a2,
                                           &bd_mb);

                at_offset += molt.atoms.nr;
            }
        }
    }

    if (mtop.bIntermolecularInteractions)
    {
        GMX_RELEASE_ASSERT(mtop.intermolecular_ilist,
                           "We should have an ilist when intermolecular interactions are on");

        bonded_distance_intermol(
                *mtop.intermolecular_ilist, ddBondedChecking, x, inputrec.pbcType, box, &bd_2b, &bd_mb);
    }

    *r_2b = std::sqrt(bd_2b.r2);
    *r_mb = std::sqrt(bd_mb.r2);

    if (*r_2b > 0 || *r_mb > 0)
    {
        GMX_LOG(mdlog.info).appendText("Initial maximum distances in bonded interactions:");
        if (*r_2b > 0)
        {
            GMX_LOG(mdlog.info)
                    .appendTextFormatted(
                            "    two-body bonded interactions: %5.3f nm, %s, atoms %d %d",
                            *r_2b,
                            (bd_2b.ftype >= 0) ? interaction_function[bd_2b.ftype].longname : "Exclusion",
                            bd_2b.a1 + 1,
                            bd_2b.a2 + 1);
        }
        if (*r_mb > 0)
        {
            GMX_LOG(mdlog.info)
                    .appendTextFormatted(
                            "  multi-body bonded interactions: %5.3f nm, %s, atoms %d %d",
                            *r_mb,
                            interaction_function[bd_mb.ftype].longname,
                            bd_mb.a1 + 1,
                            bd_mb.a2 + 1);
        }
    }
}
