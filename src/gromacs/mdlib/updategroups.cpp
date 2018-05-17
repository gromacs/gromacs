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

/*! \brief The stride of the constraint ilist iatoms list */
static int c_constraintStride = 1 + NRAL(F_CONSTR);

/*! \brief Returns whether \p moltype contains flexible constraints */
static bool hasFlexibleConstraints(const gmx_moltype_t &moltype,
                                   const t_iparams     *iparams)
{
    bool hasFlexibleConstraints = false;

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if ((interaction_function[ftype].flags & IF_CONSTRAINT) &&
            ftype != F_SETTLE)
        {
            const InteractionList &ilist = moltype.ilist[ftype];

            for (int i = 0; i < ilist.size(); i += c_constraintStride)
            {
                if (isConstraintFlexible(iparams, ilist.iatoms[i]))
                {
                    hasFlexibleConstraints = true;
                }
            }
        }
    }

    return hasFlexibleConstraints;
}

/*! \brief Returns whether moltype has incompatible vsites.
 *
 * For simpliticy the only compatible vsites are linear 2 or 3 atom sites
 * that are constructed in between the 2 or 3 contructing atoms,
 */
static bool hasIncompatibleVsites(const gmx_moltype_t &moltype,
                                  const t_iparams     *iparams)
{
    bool hasIncompatibleVsites = false;

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        const InteractionList &ilist = moltype.ilist[ftype];

        if ((interaction_function[ftype].flags & IF_VSITE) &&
            ilist.size() > 0)
        {
            if (ftype == F_VSITE2 || ftype == F_VSITE3)
            {
                for (int i = 0; i < ilist.size(); i += 1 + NRAL(ftype))
                {
                    const t_iparams &iparam   = iparams[ilist.iatoms[i]];
                    real             coeffMin = iparam.vsite.a;
                    real             coeffSum = iparam.vsite.a;
                    if (ftype == F_VSITE3)
                    {
                        coeffMin  = std::min(coeffMin, iparam.vsite.b);
                        coeffSum += iparam.vsite.b;
                    }
                    if (coeffMin < 0 || coeffSum > 1)
                    {
                        hasIncompatibleVsites = true;
                    }
                }
            }
            else
            {
                hasIncompatibleVsites = true;
            }
        }
    }

    return hasIncompatibleVsites;
}

/*! \brief Returns the other atom in a constraint involving \p ourAtomIndex */
static int getOtherConstraintAtom(int            constraintIndex,
                                  const t_ilist &ilistConstraints,
                                  int            ourAtomIndex)
{
    const t_iatom *iatoms = ilistConstraints.iatoms + constraintIndex*c_constraintStride;

    if (iatoms[1] == ourAtomIndex)
    {
        return iatoms[2];
    }
    else
    {
        return iatoms[1];
    }
}

/*! \brief Returns a merged list with constraints of all types */
static InteractionList jointConstraintList(const gmx_moltype_t &moltype)
{
    InteractionList   ilistCombined;
    std::vector<int> &iatoms = ilistCombined.iatoms;

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_CONSTRAINT)
        {
            const InteractionList &ilist = moltype.ilist[ftype];

            if (ftype == F_SETTLE)
            {
                const int settleStride = 1 + NRAL(ftype);

                for (int i = 0; i < ilist.size(); i += settleStride)
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
                GMX_RELEASE_ASSERT(NRAL(ftype) == 2, "Can only handle two atom non-SETTLE constraints");

                iatoms.insert(iatoms.end(),
                              ilist.iatoms.begin(), ilist.iatoms.end());
            }
        }
    }

    return ilistCombined;
}

/*! \brief Struct for returning an atom range */
struct AtomExtremes
{
    int minAtom; //!< The minimum atom index
    int maxAtom; //!< The maximum atom index
};

/*! \brief Returns the range of constructing atom for vsite with atom index \p a */
static AtomExtremes
vsiteConstructRange(int                  a,
                    const gmx_moltype_t &moltype)
{
    AtomExtremes extremes = { -1, -1 };

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            const InteractionList &ilist = moltype.ilist[ftype];
            for (int i = 0; i < ilist.size(); i += 1 + NRAL(ftype))
            {
                if (ilist.iatoms[i + 1] == a)
                {
                    extremes.minAtom = ilist.iatoms[i + 2];
                    extremes.maxAtom = ilist.iatoms[i + 2];
                    for (int j = i + 3; j < i + 1 + NRAL(ftype); j++)
                    {
                        extremes.minAtom = std::min(extremes.minAtom, ilist.iatoms[j]);
                        extremes.maxAtom = std::max(extremes.maxAtom, ilist.iatoms[j]);
                    }
                    return extremes;
                }
            }
        }
    }

    GMX_RELEASE_ASSERT(false, "If a is a vsite, we should have found constructing atoms");

    return extremes;
}

/*! \brief Returns a list that tells whether atoms in \p moltype are vsites */
static std::vector<bool> buildIsParticleVsite(const gmx_moltype_t &moltype)
{
    std::vector<bool> isVsite(moltype.atoms.nr);

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            const InteractionList &ilist = moltype.ilist[ftype];
            for (int i = 0; i < ilist.size(); i += 1 + NRAL(ftype))
            {
                isVsite[ilist.iatoms[i + 1]] = true;
            }
        }
    }

    return isVsite;
}

/*! \brief Returns the size of the update group starting at \p firstAtom or 0 when criteria are not met */
static int detectGroup(int                  firstAtom,
                       const gmx_moltype_t &moltype,
                       const t_blocka      &at2con,
                       const t_ilist       &ilistConstraints)
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

    int numAtomsWithConstraints = 0;
    int maxConstraintsPerAtom   = 0;
    int lastAtom                = firstAtom;
    int a                       = firstAtom;
    while (a <= lastAtom)
    {
        if (isParticleVsite[a])
        {
            AtomExtremes extremes = vsiteConstructRange(a, moltype);
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
            numAtomsWithConstraints += 1;
            maxConstraintsPerAtom    = std::max(maxConstraintsPerAtom, numConstraints);

            for (int i = at2con.index[a]; i < at2con.index[a + 1]; i++)
            {
                int otherAtom =
                    getOtherConstraintAtom(at2con.a[i], ilistConstraints, a);
                if (otherAtom < firstAtom)
                {
                    /* Constraint to atom outside the "group" */
                    return 0;
                }
                lastAtom = std::max(lastAtom, otherAtom);
            }
        }

        a++;
    }

    /* lastAtom might be followed by a vsite that is constructed from atoms
     * with index <= lastAtom. Handle this case.
     */
    if (lastAtom + 1 < moltype.atoms.nr &&
        isParticleVsite[lastAtom + 1])
    {
        AtomExtremes extremes = vsiteConstructRange(lastAtom + 1, moltype);
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
                       "so each atom should have fewer contraints than the number of atoms");
    /* Check that at least one atom is constrained to all others */
    if (maxConstraintsPerAtom < numAtomsWithConstraints - 1)
    {
        return 0;
    }

    return lastAtom - firstAtom + 1;
}

/*! \brief Returns a list of update groups for \p moltype */
static RangePartitioning makeUpdateGroups(const gmx_moltype_t &moltype,
                                          const t_iparams     *iparams)
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
    /* We "include" flexible constraints, but none a present (checked above) */
    t_blocka             at2con         = make_at2con(moltype.atoms.nr,
                                                      ilistsCombined, iparams,
                                                      FlexibleConstraintTreatment::Include);

    bool satisfiesCriteria = true;

    int  firstAtom         = 0;
    while (satisfiesCriteria && firstAtom < moltype.atoms.nr)
    {
        int numAtomsInGroup =
            detectGroup(firstAtom, moltype, at2con, ilistsCombined[F_CONSTR]);

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
static real computeMaxUpdateGroupRadius(const gmx_moltype_t     &moltype,
                                        const t_iparams         *iparams,
                                        const RangePartitioning &updateGroups)
{
    GMX_RELEASE_ASSERT(!hasFlexibleConstraints(moltype, iparams), "Flexible constraints are not supported here");

    t_blocka   at2con = make_at2con(moltype, iparams,
                                    FlexibleConstraintTreatment::Include);

    const int  stride = 1 + NRAL(F_CONSTR);
    const real half   = 0.5;

    real       maxRadius = 0;
    for (int group = 0; group < updateGroups.numBlocks(); group++)
    {
        int a = updateGroups.block(group).begin();
        for (int i = at2con.index[a]; i < at2con.index[a + 1]; i++)
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
            if (at2con.index[a + 1] == at2con.index[a] + 1)
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

    const InteractionList &settles = moltype.ilist[F_SETTLE];
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

    GMX_RELEASE_ASSERT(static_cast<size_t>(updateGroups.size()) == mtop.moltype.size(), "We need one update group entry per moleculetype");

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
