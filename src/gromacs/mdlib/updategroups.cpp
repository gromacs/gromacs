/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/* \internal \file
 *
 * \brief Implements functions for generating update groups
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include "updategroups.h"

#include <cmath>

#include "gromacs/mdlib/constr.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"

namespace gmx
{

/*! \brief Returns whether \p moltype contains flexible constraints */
static bool hasFlexibleConstraints(const gmx_moltype_t            &moltype,
                                   gmx::ArrayRef<const t_iparams>  iparams)
{
    for (auto &ilist : extractILists(moltype.ilist, IF_CONSTRAINT))
    {
        if (ilist.functionType != F_SETTLE)
        {
            for (size_t i = 0; i < ilist.iatoms.size(); i += ilistStride(ilist))
            {
                if (isConstraintFlexible(iparams.data(), ilist.iatoms[i]))
                {
                    return true;
                }
            }
        }
    }

    return false;
}

/*! \brief Returns whether moltype has incompatible vsites.
 *
 * For simplicity the only compatible vsites are linear 2 or 3 atom sites
 * that are constructed in between the 2 or 3 contructing atoms,
 */
static bool hasIncompatibleVsites(const gmx_moltype_t            &moltype,
                                  gmx::ArrayRef<const t_iparams>  iparams)
{
    bool hasIncompatibleVsites = false;

    for (auto &ilist : extractILists(moltype.ilist, IF_VSITE))
    {
        if (ilist.functionType == F_VSITE2 || ilist.functionType == F_VSITE3)
        {
            for (size_t i = 0; i < ilist.iatoms.size(); i += ilistStride(ilist))
            {
                const t_iparams &iparam = iparams[ilist.iatoms[i]];
                real             coeffMin;
                real             coeffSum;
                if (ilist.functionType == F_VSITE2)
                {
                    coeffMin = iparam.vsite.a;
                    coeffSum = iparam.vsite.a;
                }
                else
                {
                    coeffMin = std::min(iparam.vsite.a, iparam.vsite.b);
                    coeffSum = iparam.vsite.a + iparam.vsite.b;
                }
                if (coeffMin < 0 || coeffSum > 1)
                {
                    hasIncompatibleVsites = true;
                    break;
                }
            }
        }
        else
        {
            hasIncompatibleVsites = true;
            break;
        }
    }

    return hasIncompatibleVsites;
}

/*! \brief Returns a merged list with constraints of all types */
static InteractionList jointConstraintList(const gmx_moltype_t &moltype)
{
    InteractionList   ilistCombined;
    std::vector<int> &iatoms = ilistCombined.iatoms;

    for (auto &ilist : extractILists(moltype.ilist, IF_CONSTRAINT))
    {
        if (ilist.functionType == F_SETTLE)
        {
            for (size_t i = 0; i < ilist.iatoms.size(); i += ilistStride(ilist))
            {
                iatoms.push_back(-1);
                iatoms.push_back(ilist.iatoms[i + 1]);
                iatoms.push_back(ilist.iatoms[i + 2]);
                iatoms.push_back(-1);
                iatoms.push_back(ilist.iatoms[i + 1]);
                iatoms.push_back(ilist.iatoms[i + 3]);
                iatoms.push_back(-1);
                iatoms.push_back(ilist.iatoms[i + 2]);
                iatoms.push_back(ilist.iatoms[i + 3]);
            }
        }
        else
        {
            GMX_RELEASE_ASSERT(NRAL(ilist.functionType) == 2, "Can only handle two-atom non-SETTLE constraints");

            iatoms.insert(iatoms.end(),
                          ilist.iatoms.begin(), ilist.iatoms.end());
        }
    }

    return ilistCombined;
}

/*! \brief Struct for returning an atom range */
struct AtomIndexExtremes
{
    int minAtom; //!< The minimum atom index
    int maxAtom; //!< The maximum atom index
};

/*! \brief Returns the range of constructing atom for vsite with atom index \p a */
static AtomIndexExtremes
vsiteConstructRange(int                  a,
                    const gmx_moltype_t &moltype)
{
    AtomIndexExtremes extremes = { -1, -1 };

    for (auto &ilist : extractILists(moltype.ilist, IF_VSITE))
    {
        for (size_t i = 0; i < ilist.iatoms.size(); i += ilistStride(ilist))
        {
            if (ilist.iatoms[i + 1] == a)
            {
                extremes.minAtom = ilist.iatoms[i + 2];
                extremes.maxAtom = ilist.iatoms[i + 2];
                for (size_t j = i + 3; j < i + ilistStride(ilist); j++)
                {
                    extremes.minAtom = std::min(extremes.minAtom, ilist.iatoms[j]);
                    extremes.maxAtom = std::max(extremes.maxAtom, ilist.iatoms[j]);
                }
                return extremes;
            }
        }
    }

    GMX_RELEASE_ASSERT(false, "If a is a vsite, we should have found constructing atoms");

    return extremes;
}

/*! \brief Returns the range of atoms constrained to atom \p a (including \p a itself) */
static AtomIndexExtremes
constraintAtomRange(int                    a,
                    const t_blocka        &at2con,
                    const InteractionList &ilistConstraints)
{
    AtomIndexExtremes extremes = { a, a };

    for (int i = at2con.index[a]; i < at2con.index[a + 1]; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            int atomJ        = ilistConstraints.iatoms[at2con.a[i]*3 + 1 + j];
            extremes.minAtom = std::min(extremes.minAtom, atomJ);
            extremes.maxAtom = std::max(extremes.maxAtom, atomJ);
        }
    }

    return extremes;
}

/*! \brief Returns a list that tells whether atoms in \p moltype are vsites */
static std::vector<bool> buildIsParticleVsite(const gmx_moltype_t &moltype)
{
    std::vector<bool> isVsite(moltype.atoms.nr);

    for (auto &ilist : extractILists(moltype.ilist, IF_VSITE))
    {
        for (size_t i = 0; i < ilist.iatoms.size(); i += ilistStride(ilist))
        {
            int vsiteAtom      = ilist.iatoms[i + 1];
            isVsite[vsiteAtom] = true;
        }
    }

    return isVsite;
}

/*! \brief Returns the size of the update group starting at \p firstAtom or 0 when criteria (see updategroups.h) are not met */
static int detectGroup(int                    firstAtom,
                       const gmx_moltype_t   &moltype,
                       const t_blocka        &at2con,
                       const InteractionList &ilistConstraints)
{
    /* We should be using moltype.atoms.atom[].ptype for checking whether
     * a particle is a vsite. But the test code can't fill t_atoms,
     * because it uses C pointers which get double freed.
     */
    std::vector<bool> isParticleVsite = buildIsParticleVsite(moltype);

    /* A non-vsite atom without constraints is an update group by itself */
    if (!isParticleVsite[firstAtom] &&
        at2con.index[firstAtom + 1] - at2con.index[firstAtom] == 0)
    {
        return 1;
    }

    /* A (potential) update group starts at firstAtom and should consist
     * of two or more atoms and possibly vsites. At least one atom should
     * have constraints with all other atoms and vsites should have all
     * constructing atoms inside the group. Here we increase lastAtom until
     * the criteria are fulfilled or exit when criteria cannot be met.
     */
    int numAtomsWithConstraints = 0;
    int maxConstraintsPerAtom   = 0;
    int lastAtom                = firstAtom;
    int a                       = firstAtom;
    while (a <= lastAtom)
    {
        if (isParticleVsite[a])
        {
            AtomIndexExtremes extremes = vsiteConstructRange(a, moltype);
            if (extremes.minAtom < firstAtom)
            {
                /* A constructing atom is outside the group,
                 * we can not use update groups.
                 */
                return 0;
            }
            lastAtom = std::max(lastAtom, extremes.maxAtom);
        }
        else
        {
            int numConstraints = at2con.index[a + 1] - at2con.index[a];
            if (numConstraints == 0)
            {
                /* We can not have unconstrained atoms in an update group */
                return 0;
            }
            /* This atom has at least one constraint.
             * Check whether all constraints are within the group
             * and whether we need to extend the group.
             */
            numAtomsWithConstraints += 1;
            maxConstraintsPerAtom    = std::max(maxConstraintsPerAtom, numConstraints);

            AtomIndexExtremes extremes =
                constraintAtomRange(a, at2con, ilistConstraints);
            if (extremes.minAtom < firstAtom)
            {
                /* Constraint to atom outside the "group" */
                return 0;
            }
            lastAtom = std::max(lastAtom, extremes.maxAtom);
        }

        a++;
    }

    /* lastAtom might be followed by a vsite that is constructed from atoms
     * with index <= lastAtom. Handle this case.
     */
    if (lastAtom + 1 < moltype.atoms.nr &&
        isParticleVsite[lastAtom + 1])
    {
        AtomIndexExtremes extremes = vsiteConstructRange(lastAtom + 1, moltype);
        if (extremes.minAtom < firstAtom)
        {
            /* Constructing atom precedes the group */
            return 0;
        }
        else if (extremes.maxAtom <= lastAtom)
        {
            /* All constructing atoms are in the group, add the vsite to the group */
            lastAtom++;
        }
        else if (extremes.minAtom <= lastAtom)
        {
            /* Some, but not all constructing atoms are in the group */
            return 0;
        }
    }

    GMX_RELEASE_ASSERT(maxConstraintsPerAtom < numAtomsWithConstraints,
                       "We have checked that atoms are only constrained to atoms within the group,"
                       "so each atom should have fewer constraints than the number of atoms");
    /* Check that at least one atom is constrained to all others */
    if (maxConstraintsPerAtom != numAtomsWithConstraints - 1)
    {
        return 0;
    }

    return lastAtom - firstAtom + 1;
}

/*! \brief Returns a list of update groups for \p moltype */
static RangePartitioning
makeUpdateGroups(const gmx_moltype_t            &moltype,
                 gmx::ArrayRef<const t_iparams>  iparams)
{
    RangePartitioning groups;

    /* Update groups are not compatible with flexible constraints.
     * Without dynamics the flexible constraints are ignored,
     * but since performance for EM/NM is less critical, we do not
     * use update groups to keep the code here simpler.
     */
    if (hasFlexibleConstraints(moltype, iparams) ||
        hasIncompatibleVsites(moltype, iparams))
    {
        return groups;
    }

    /* Combine all constraint ilists into a single one */
    InteractionList constraintsCombined = jointConstraintList(moltype);
    t_ilist         ilistsCombined[F_NRE];
    ilistsCombined[F_CONSTR].nr         = constraintsCombined.iatoms.size();
    ilistsCombined[F_CONSTR].iatoms     = constraintsCombined.iatoms.data();
    ilistsCombined[F_CONSTRNC].nr       = 0;
    /* We "include" flexible constraints, but none are present (checked above) */
    t_blocka             at2con         = make_at2con(moltype.atoms.nr,
                                                      ilistsCombined, iparams.data(),
                                                      FlexibleConstraintTreatment::Include);

    bool satisfiesCriteria = true;

    int  firstAtom         = 0;
    while (satisfiesCriteria && firstAtom < moltype.atoms.nr)
    {
        int numAtomsInGroup =
            detectGroup(firstAtom, moltype, at2con, constraintsCombined);

        if (numAtomsInGroup == 0)
        {
            satisfiesCriteria = false;
        }
        else
        {
            groups.appendBlock(numAtomsInGroup);
        }
        firstAtom += numAtomsInGroup;
    }

    if (!satisfiesCriteria)
    {
        /* Make groups empty, to signal not satisfying the criteria */
        groups.clear();
    }

    done_blocka(&at2con);

    return groups;
}

std::vector<RangePartitioning> makeUpdateGroups(const gmx_mtop_t &mtop)
{
    std::vector<RangePartitioning> updateGroups;

    bool                           systemSatisfiesCriteria = true;
    for (const gmx_moltype_t &moltype : mtop.moltype)
    {
        updateGroups.push_back(makeUpdateGroups(moltype, mtop.ffparams.iparams));

        if (updateGroups.back().numBlocks() == 0)
        {
            systemSatisfiesCriteria = false;
        }
    }

    if (!systemSatisfiesCriteria)
    {
        updateGroups.clear();
    }

    return updateGroups;
}

/*! \brief Returns the maximum update group radius for \p moltype */
static real
computeMaxUpdateGroupRadius(const gmx_moltype_t            &moltype,
                            gmx::ArrayRef<const t_iparams>  iparams,
                            const RangePartitioning        &updateGroups)
{
    GMX_RELEASE_ASSERT(!hasFlexibleConstraints(moltype, iparams),
                       "Flexible constraints are not supported here");

    const InteractionList &settles = moltype.ilist[F_SETTLE];

    t_blocka               at2con =
        make_at2con(moltype, iparams, FlexibleConstraintTreatment::Include);

    const int  stride = 1 + NRAL(F_CONSTR);
    const real half   = 0.5;

    real       maxRadius = 0;
    for (int group = 0; group < updateGroups.numBlocks(); group++)
    {
        if (updateGroups.block(group).size() == 1)
        {
            /* Single atom group, radius is zero */
            continue;
        }

        /* Find the atom maxAtom with the maximum number of constraints */
        int maxNumConstraints = 0;
        int maxAtom           = -1;
        for (int a : updateGroups.block(group))
        {
            int numConstraints = at2con.index[a + 1] - at2con.index[a];
            if (numConstraints > maxNumConstraints)
            {
                maxNumConstraints = numConstraints;
                maxAtom           = a;
            }
        }
        GMX_ASSERT(maxAtom >= 0 || settles.size() > 0,
                   "We should have at least two atoms in the group with constraints");
        if (maxAtom < 0)
        {
            continue;
        }

        for (int i = at2con.index[maxAtom]; i < at2con.index[maxAtom + 1]; i++)
        {
            int conIndex = at2con.a[i]*stride;
            int iparamsIndex;
            if (conIndex < moltype.ilist[F_CONSTR].size())
            {
                iparamsIndex = moltype.ilist[F_CONSTR].iatoms[conIndex];
            }
            else
            {
                iparamsIndex = moltype.ilist[F_CONSTRNC].iatoms[conIndex - moltype.ilist[F_CONSTR].size()];
            }
            /* Here we take the maximum constraint length of the A and B
             * topology, which assumes lambda is between 0 and 1 for
             * free-energy runs.
             */
            real constraintLength = std::max(iparams[iparamsIndex].constr.dA,
                                             iparams[iparamsIndex].constr.dB);
            if (at2con.index[maxAtom + 1] == at2con.index[maxAtom] + 1)
            {
                /* Single constraint: use the distance from the midpoint */
                maxRadius = std::max(maxRadius, half*constraintLength);
            }
            else
            {
                /* Multiple constraints: use the distance from the central
                 * atom. We can do better than this if we would assume
                 * that the angles between bonds do not stretch a lot.
                 */
                maxRadius = std::max(maxRadius, constraintLength);
            }
        }
    }

    for (int i = 0; i < settles.size(); i += 1 + NRAL(F_SETTLE))
    {
        const real dOH   = iparams[settles.iatoms[i]].settle.doh;
        const real dHH   = iparams[settles.iatoms[i]].settle.dhh;
        /* Compute distance^2 from center of geometry to O and H */
        const real dCO2  = (4*dOH*dOH - dHH*dHH)/9;
        const real dCH2  = (dOH*dOH + 2*dHH*dHH)/9;
        const real dCAny = std::sqrt(std::max(dCO2, dCH2));
        maxRadius        = std::max(maxRadius, dCAny);
    }

    done_blocka(&at2con);

    return maxRadius;
}

real computeMaxUpdateGroupRadius(const gmx_mtop_t                       &mtop,
                                 gmx::ArrayRef<const RangePartitioning>  updateGroups)
{
    if (updateGroups.empty())
    {
        return 0;
    }

    GMX_RELEASE_ASSERT(static_cast<size_t>(updateGroups.size()) == mtop.moltype.size(),
                       "We need one update group entry per moleculetype");

    real maxRadius = 0;

    for (size_t moltype = 0; moltype < mtop.moltype.size(); moltype++)
    {
        maxRadius
            = std::max(maxRadius,
                       computeMaxUpdateGroupRadius(mtop.moltype[moltype],
                                                   mtop.ffparams.iparams,
                                                   updateGroups[moltype]));
    }

    return maxRadius;
}

} // namespace gmx
