/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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

#include <unordered_map>
#include <variant>

#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/message_string_collector.h"

namespace gmx
{

//! Reasons why the system can be incompatible with update groups
enum class IncompatibilityReasons
{
    FlexibleConstraint,
    IncompatibleVsite,
    VsiteConstructingAtomsSplit,
    ConstrainedAtomOrder,
    NoCentralConstraintAtom,
    Count
};

//! Strings explaining why the system is incompatible with update groups
static const EnumerationArray<IncompatibilityReasons, std::string> reasonStrings = {
    "flexible constraints are present",
    "an incompatible virtual site type is used",
    "the construction atoms of a virtual site are only partly with a group of constrained atoms",
    "atoms that are (in)directly constrained together are interdispersed with other atoms",
    "there are three or more consecutively coupled constraints"
};

/*! \brief Returns whether \p moltype contains flexible constraints */
static bool hasFlexibleConstraints(const gmx_moltype_t& moltype, gmx::ArrayRef<const t_iparams> iparams)
{
    for (auto& ilist : extractILists(moltype.ilist, IF_CONSTRAINT))
    {
        if (ilist.functionType != F_SETTLE)
        {
            for (size_t i = 0; i < ilist.iatoms.size(); i += ilistStride(ilist))
            {
                if (isConstraintFlexible(iparams, ilist.iatoms[i]))
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
static bool hasIncompatibleVsites(const gmx_moltype_t& moltype, gmx::ArrayRef<const t_iparams> iparams)
{
    bool hasIncompatibleVsites = false;

    for (auto& ilist : extractILists(moltype.ilist, IF_VSITE))
    {
        if (ilist.functionType == F_VSITE2 || ilist.functionType == F_VSITE3)
        {
            for (size_t i = 0; i < ilist.iatoms.size(); i += ilistStride(ilist))
            {
                const t_iparams& iparam = iparams[ilist.iatoms[i]];
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
static InteractionList jointConstraintList(const gmx_moltype_t& moltype)
{
    InteractionList   ilistCombined;
    std::vector<int>& iatoms = ilistCombined.iatoms;

    for (auto& ilist : extractILists(moltype.ilist, IF_CONSTRAINT))
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
            GMX_RELEASE_ASSERT(NRAL(ilist.functionType) == 2,
                               "Can only handle two-atom non-SETTLE constraints");

            iatoms.insert(iatoms.end(), ilist.iatoms.begin(), ilist.iatoms.end());
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
static AtomIndexExtremes vsiteConstructRange(int a, const gmx_moltype_t& moltype)
{
    AtomIndexExtremes extremes = { -1, -1 };

    for (auto& ilist : extractILists(moltype.ilist, IF_VSITE))
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
static AtomIndexExtremes constraintAtomRange(int                     a,
                                             const ListOfLists<int>& at2con,
                                             const InteractionList&  ilistConstraints)
{
    AtomIndexExtremes extremes = { a, a };

    for (const int constraint : at2con[a])
    {
        for (int j = 0; j < 2; j++)
        {
            int atomJ        = ilistConstraints.iatoms[constraint * 3 + 1 + j];
            extremes.minAtom = std::min(extremes.minAtom, atomJ);
            extremes.maxAtom = std::max(extremes.maxAtom, atomJ);
        }
    }

    return extremes;
}

/*! \brief Returns a list that tells whether atoms in \p moltype are vsites */
static std::vector<bool> buildIsParticleVsite(const gmx_moltype_t& moltype)
{
    std::vector<bool> isVsite(moltype.atoms.nr);

    for (auto& ilist : extractILists(moltype.ilist, IF_VSITE))
    {
        for (size_t i = 0; i < ilist.iatoms.size(); i += ilistStride(ilist))
        {
            int vsiteAtom      = ilist.iatoms[i + 1];
            isVsite[vsiteAtom] = true;
        }
    }

    return isVsite;
}

/*! \brief Returns the size of the update group starting at \p firstAtom or an incompatibility reason */
static std::variant<int, IncompatibilityReasons> detectGroup(int                     firstAtom,
                                                             const gmx_moltype_t&    moltype,
                                                             const ListOfLists<int>& at2con,
                                                             const InteractionList& ilistConstraints)
{
    /* We should be using moltype.atoms.atom[].ptype for checking whether
     * a particle is a vsite. But the test code can't fill t_atoms,
     * because it uses C pointers which get double freed.
     */
    std::vector<bool> isParticleVsite = buildIsParticleVsite(moltype);

    /* A non-vsite atom without constraints is an update group by itself */
    if (!isParticleVsite[firstAtom] && at2con[firstAtom].empty())
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
                return IncompatibilityReasons::VsiteConstructingAtomsSplit;
            }
            lastAtom = std::max(lastAtom, extremes.maxAtom);
        }
        else
        {
            const int numConstraints = at2con[a].ssize();
            if (numConstraints == 0)
            {
                /* We can not have unconstrained atoms in an update group */
                return IncompatibilityReasons::ConstrainedAtomOrder;
            }
            /* This atom has at least one constraint.
             * Check whether all constraints are within the group
             * and whether we need to extend the group.
             */
            numAtomsWithConstraints += 1;
            maxConstraintsPerAtom = std::max(maxConstraintsPerAtom, numConstraints);

            AtomIndexExtremes extremes = constraintAtomRange(a, at2con, ilistConstraints);
            if (extremes.minAtom < firstAtom)
            {
                /* Constraint to atom outside the "group" */
                return IncompatibilityReasons::ConstrainedAtomOrder;
            }
            lastAtom = std::max(lastAtom, extremes.maxAtom);
        }

        a++;
    }

    /* lastAtom might be followed by a vsite that is constructed from atoms
     * with index <= lastAtom. Handle this case.
     */
    if (lastAtom + 1 < moltype.atoms.nr && isParticleVsite[lastAtom + 1])
    {
        AtomIndexExtremes extremes = vsiteConstructRange(lastAtom + 1, moltype);
        if (extremes.minAtom < firstAtom)
        { // NOLINT bugprone-branch-clone
            /* Constructing atom precedes the group */
            return IncompatibilityReasons::VsiteConstructingAtomsSplit;
        }
        else if (extremes.maxAtom <= lastAtom)
        {
            /* All constructing atoms are in the group, add the vsite to the group */
            lastAtom++;
        }
        else if (extremes.minAtom <= lastAtom)
        {
            /* Some, but not all constructing atoms are in the group */
            return IncompatibilityReasons::VsiteConstructingAtomsSplit;
        }
    }

    GMX_RELEASE_ASSERT(maxConstraintsPerAtom < numAtomsWithConstraints,
                       "We have checked that atoms are only constrained to atoms within the group,"
                       "so each atom should have fewer constraints than the number of atoms");
    /* Check that at least one atom is constrained to all others */
    if (maxConstraintsPerAtom != numAtomsWithConstraints - 1)
    {
        return IncompatibilityReasons::NoCentralConstraintAtom;
    }

    return lastAtom - firstAtom + 1;
}

/*! \brief Returns a list of update groups for \p moltype or an incompatibility reason */
static std::variant<RangePartitioning, IncompatibilityReasons>
makeUpdateGroupingsPerMoleculeType(const gmx_moltype_t& moltype, gmx::ArrayRef<const t_iparams> iparams)
{
    RangePartitioning groups;

    /* Update groups are not compatible with flexible constraints.
     * Without dynamics the flexible constraints are ignored,
     * but since performance for EM/NM is less critical, we do not
     * use update groups to keep the code here simpler.
     */
    if (hasFlexibleConstraints(moltype, iparams))
    {
        return IncompatibilityReasons::FlexibleConstraint;
    }
    if (hasIncompatibleVsites(moltype, iparams))
    {
        return IncompatibilityReasons::IncompatibleVsite;
    }

    /* Combine all constraint ilists into a single one */
    std::array<InteractionList, F_NRE> ilistsCombined;
    ilistsCombined[F_CONSTR] = jointConstraintList(moltype);
    /* We "include" flexible constraints, but none are present (checked above) */
    const ListOfLists<int> at2con = make_at2con(
            moltype.atoms.nr, ilistsCombined, iparams, FlexibleConstraintTreatment::Include);

    int firstAtom = 0;
    while (firstAtom < moltype.atoms.nr)
    {
        const auto detectionResult = detectGroup(firstAtom, moltype, at2con, ilistsCombined[F_CONSTR]);

        if (std::holds_alternative<IncompatibilityReasons>(detectionResult))
        {
            // Can not use update groups, return the reason for incompatiblity
            return std::get<IncompatibilityReasons>(detectionResult);
        }

        const int numAtomsInGroup = std::get<int>(detectionResult);

        groups.appendBlock(numAtomsInGroup);

        firstAtom += numAtomsInGroup;
    }

    return groups;
}

std::variant<std::vector<RangePartitioning>, std::string> makeUpdateGroupingsPerMoleculeType(const gmx_mtop_t& mtop)
{
    std::vector<RangePartitioning> updateGroupingsPerMoleculeType;

    for (const gmx_moltype_t& moltype : mtop.moltype)
    {
        const auto detectionResult = makeUpdateGroupingsPerMoleculeType(moltype, mtop.ffparams.iparams);

        if (std::holds_alternative<IncompatibilityReasons>(detectionResult))
        {
            // Can not use update groups, return a string with the reason for incompatiblity
            return reasonStrings[std::get<IncompatibilityReasons>(detectionResult)];
        }

        updateGroupingsPerMoleculeType.push_back(std::get<RangePartitioning>(detectionResult));
    }

    return updateGroupingsPerMoleculeType;
}

/*! \brief Returns a map of angles ilist.iatoms indices with the middle atom as key */
static std::unordered_multimap<int, int> getAngleIndices(const gmx_moltype_t& moltype)
{
    const InteractionList& angles = moltype.ilist[F_ANGLES];

    std::unordered_multimap<int, int> indices(angles.size());

    for (int i = 0; i < angles.size(); i += 1 + NRAL(F_ANGLES))
    {
        indices.insert({ angles.iatoms[i + 2], i });
    }

    return indices;
}

/*! \brief When possible, computes the maximum radius of constrained atom in an update group
 *
 * Supports groups with 2 or 3 atoms where all partner atoms are connected to
 * each other by angle potentials. The temperature is used to compute a radius
 * that is not exceeded with a chance of 10^-9. Note that this computation
 * assumes there are no other strong forces working on these angular
 * degrees of freedom.
 * The return value is -1 when all partners are not connected to each other
 * by one angle potential, when a potential is perturbed or when an angle
 * could reach more than 180 degrees.
 */
template<int numPartnerAtoms>
static real constraintGroupRadius(const gmx_moltype_t&                     moltype,
                                  gmx::ArrayRef<const t_iparams>           iparams,
                                  const int                                centralAtom,
                                  const ListOfLists<int>&                  at2con,
                                  const std::unordered_multimap<int, int>& angleIndices,
                                  const real                               constraintLength,
                                  const real                               temperature)
{
    const int numConstraints = at2con[centralAtom].ssize();
    GMX_RELEASE_ASSERT(numConstraints == numPartnerAtoms,
                       "We expect as many constraints as partner atoms here");

    std::array<int, numPartnerAtoms> partnerAtoms;
    for (int i = 0; i < numPartnerAtoms; i++)
    {
        const int ind = at2con[centralAtom][i] * 3;
        if (ind >= moltype.ilist[F_CONSTR].size())
        {
            /* This is a flexible constraint, we don't optimize for that */
            return -1;
        }
        const int a1    = moltype.ilist[F_CONSTR].iatoms[ind + 1];
        const int a2    = moltype.ilist[F_CONSTR].iatoms[ind + 2];
        partnerAtoms[i] = (a1 == centralAtom ? a2 : a1);
    }

    const InteractionList&           angles      = moltype.ilist[F_ANGLES];
    auto                             range       = angleIndices.equal_range(centralAtom);
    int                              angleType   = -1;
    std::array<int, numPartnerAtoms> numAngles   = { 0 };
    bool                             areSameType = true;
    for (auto it = range.first; it != range.second; ++it)
    {
        /* Check if the outer atoms in the angle are both partner atoms */
        int numAtomsFound = 0;
        for (int ind = it->second + 1; ind < it->second + 4; ind += 2)
        {
            for (const int& partnerAtom : partnerAtoms)
            {
                if (angles.iatoms[ind] == partnerAtom)
                {
                    numAtomsFound++;
                }
            }
        }
        if (numAtomsFound == 2)
        {
            /* Check that the angle potentials have the same type */
            if (angleType == -1)
            {
                angleType = angles.iatoms[it->second];
            }
            else if (angles.iatoms[it->second] != angleType)
            {
                areSameType = false;
            }
            /* Count the number of angle interactions per atoms */
            for (int ind = it->second + 1; ind < it->second + 4; ind += 2)
            {
                for (size_t i = 0; i < partnerAtoms.size(); i++)
                {
                    if (angles.iatoms[ind] == partnerAtoms[i])
                    {
                        numAngles[i]++;
                    }
                }
            }
        }
    }

    bool criteriaSatisfied = areSameType;
    for (int i = 0; i < numPartnerAtoms; i++)
    {
        if (numAngles[i] != numPartnerAtoms - 1)
        {
            criteriaSatisfied = false;
        }
    }

    /* We don't bother optimizing the perturbed angle case */
    const t_iparams& angleParams = iparams[angleType];
    if (criteriaSatisfied && angleParams.harmonic.rB == angleParams.harmonic.rA
        && angleParams.harmonic.krB == angleParams.harmonic.krA)
    {
        /* Set number of stddevs such that change of exceeding < 10^-9 */
        constexpr real c_numSigma = 6.0;
        /* Compute the maximally stretched angle */
        const real eqAngle = angleParams.harmonic.rA * gmx::c_deg2Rad;
        const real fc      = angleParams.harmonic.krA;
        const real maxAngle =
                eqAngle + c_numSigma * gmx::c_boltz * temperature / ((numPartnerAtoms - 1) * fc);
        if (maxAngle >= M_PI)
        {
            return -1;
        }

        if (numPartnerAtoms == 2)
        {
            /* With two atoms constrainted to a cental atom we have a triangle
             * with two equal sides because the constraint type is equal.
             * Return the distance from the COG to the farthest two corners,
             * i.e. the partner atoms.
             */
            real distMidPartner = std::sin(0.5 * maxAngle) * constraintLength;
            real distCentralMid = std::cos(0.5 * maxAngle) * constraintLength;
            real distCogMid     = distCentralMid * numPartnerAtoms / (numPartnerAtoms + 1);
            real distCogPartner = std::sqrt(gmx::square(distMidPartner) + gmx::square(distCogMid));

            return distCogPartner;
        }
        else if (numPartnerAtoms == 3)
        {
            /* With three atoms constrainted to a cental atom we have the
             * largest COG-atom distance when one partner atom (the outlier)
             * moves out by stretching its two angles by an equal amount.
             * The math here gets a bit more involved, but it is still
             * rather straightforward geometry.
             * We first compute distances in the plane of the three partners.
             * Here we have two "equilibrium" partners and one outlier.
             * We make use of the "Mid" point between the two "Eq" partners.
             * We project the central atom on this plane.
             * Then we compute the distance of the central atom to the plane.
             * The distance of the COG to the ourlier is returned.
             */
            real halfDistEqPartner  = std::sin(0.5 * eqAngle) * constraintLength;
            real distPartnerOutlier = 2 * std::sin(0.5 * maxAngle) * constraintLength;
            real distMidOutlier =
                    std::sqrt(gmx::square(distPartnerOutlier) - gmx::square(halfDistEqPartner));
            real distMidCenterInPlane =
                    0.5 * (distMidOutlier - gmx::square(halfDistEqPartner) / distMidOutlier);
            real distCentralToPlane =
                    std::sqrt(gmx::square(constraintLength) - gmx::square(halfDistEqPartner)
                              - gmx::square(distMidCenterInPlane));
            real distCogOutlierH = distCentralToPlane / (numPartnerAtoms + 1);
            real distCogOutlierP =
                    distMidOutlier - (distMidOutlier + distMidCenterInPlane) / (numPartnerAtoms + 1);
            real distCogOutlier = std::sqrt(gmx::square(distCogOutlierH) + gmx::square(distCogOutlierP));

            return distCogOutlier;
        }
        else
        {
            GMX_RELEASE_ASSERT(false, "Only 2 or 3 constraints are supported here");
        }
    }

    return -1;
}

/*! \brief Returns the maximum update group radius for \p moltype */
static real computeMaxUpdateGroupRadius(const gmx_moltype_t&           moltype,
                                        gmx::ArrayRef<const t_iparams> iparams,
                                        const RangePartitioning&       updateGrouping,
                                        real                           temperature)
{
    GMX_RELEASE_ASSERT(!hasFlexibleConstraints(moltype, iparams),
                       "Flexible constraints are not supported here");

    const InteractionList& settles = moltype.ilist[F_SETTLE];

    const ListOfLists<int> at2con = make_at2con(moltype, iparams, FlexibleConstraintTreatment::Include);

    const auto angleIndices = getAngleIndices(moltype);

    real maxRadius = 0;
    for (int group = 0; group < updateGrouping.numBlocks(); group++)
    {
        if (updateGrouping.block(group).size() == 1)
        {
            /* Single atom group, radius is zero */
            continue;
        }

        /* Find the atom maxAtom with the maximum number of constraints */
        int maxNumConstraints = 0;
        int maxAtom           = -1;
        for (int a : updateGrouping.block(group))
        {
            const int numConstraints = at2con[a].ssize();
            if (numConstraints > maxNumConstraints)
            {
                maxNumConstraints = numConstraints;
                maxAtom           = a;
            }
        }
        GMX_ASSERT(maxAtom >= 0 || !settles.empty(),
                   "We should have at least two atoms in the group with constraints");
        if (maxAtom < 0)
        {
            continue;
        }

        bool allTypesAreEqual     = true;
        int  constraintType       = -1;
        real maxConstraintLength  = 0;
        real sumConstraintLengths = 0;
        bool isFirstConstraint    = true;
        for (const int constraint : at2con[maxAtom])
        {
            int conIndex = constraint * (1 + NRAL(F_CONSTR));
            int iparamsIndex;
            if (conIndex < moltype.ilist[F_CONSTR].size())
            {
                iparamsIndex = moltype.ilist[F_CONSTR].iatoms[conIndex];
            }
            else
            {
                iparamsIndex =
                        moltype.ilist[F_CONSTRNC].iatoms[conIndex - moltype.ilist[F_CONSTR].size()];
            }
            if (isFirstConstraint)
            {
                constraintType    = iparamsIndex;
                isFirstConstraint = false;
            }
            else if (iparamsIndex != constraintType)
            {
                allTypesAreEqual = false;
            }
            /* Here we take the maximum constraint length of the A and B
             * topology, which assumes lambda is between 0 and 1 for
             * free-energy runs.
             */
            real constraintLength =
                    std::max(iparams[iparamsIndex].constr.dA, iparams[iparamsIndex].constr.dB);
            maxConstraintLength = std::max(maxConstraintLength, constraintLength);
            sumConstraintLengths += constraintLength;
        }

        int  numConstraints = at2con[maxAtom].ssize();
        real radius;
        if (numConstraints == 1)
        {
            /* Single constraint: the radius is the distance from the midpoint */
            radius = 0.5_real * maxConstraintLength;
        }
        else
        {
            radius = -1;

            /* With 2 constraints the maximum possible radius is the
             * constraint length, so we can use that for the 0 K case.
             */
            if (numConstraints == 2 && allTypesAreEqual && temperature > 0)
            {
                radius = constraintGroupRadius<2>(
                        moltype, iparams, maxAtom, at2con, angleIndices, maxConstraintLength, temperature);
            }
            /* With 3 constraints the maximum possible radius is 1.4 times
             * the constraint length, so it is worth computing a smaller
             * radius to enable update groups for systems in a small box.
             */
            if (numConstraints == 3 && allTypesAreEqual && temperature >= 0)
            {
                radius = constraintGroupRadius<3>(
                        moltype, iparams, maxAtom, at2con, angleIndices, maxConstraintLength, temperature);
                if (temperature == 0 && radius >= 0)
                {
                    /* Add a 10% margin for deviation at 0 K */
                    radius *= 1.1_real;
                }
            }

            if (radius < 0)
            {
                /* Worst case: atom with the longest constraint on one side
                 * of the center, all others on the opposite side
                 */
                radius = maxConstraintLength
                         + (sumConstraintLengths - 2 * maxConstraintLength) / (numConstraints + 1);
            }
        }
        maxRadius = std::max(maxRadius, radius);
    }

    for (int i = 0; i < settles.size(); i += 1 + NRAL(F_SETTLE))
    {
        const real dOH = iparams[settles.iatoms[i]].settle.doh;
        const real dHH = iparams[settles.iatoms[i]].settle.dhh;
        /* Compute distance^2 from center of geometry to O and H */
        const real dCO2  = (4 * dOH * dOH - dHH * dHH) / 9;
        const real dCH2  = (dOH * dOH + 2 * dHH * dHH) / 9;
        const real dCAny = std::sqrt(std::max(dCO2, dCH2));
        maxRadius        = std::max(maxRadius, dCAny);
    }

    return maxRadius;
}

real computeMaxUpdateGroupRadius(const gmx_mtop_t&                      mtop,
                                 gmx::ArrayRef<const RangePartitioning> updateGroupingsPerMoleculeType,
                                 real                                   temperature)
{
    if (updateGroupingsPerMoleculeType.empty())
    {
        return 0;
    }

    GMX_RELEASE_ASSERT(updateGroupingsPerMoleculeType.size() == mtop.moltype.size(),
                       "We need one update group entry per moleculetype");

    real maxRadius = 0;

    for (size_t moltype = 0; moltype < mtop.moltype.size(); moltype++)
    {
        const real radiusOfThisMoleculeType = computeMaxUpdateGroupRadius(
                mtop.moltype[moltype], mtop.ffparams.iparams, updateGroupingsPerMoleculeType[moltype], temperature);
        maxRadius = std::max(maxRadius, radiusOfThisMoleculeType);
    }

    return maxRadius;
}

UpdateGroups::UpdateGroups(std::vector<RangePartitioning>&& updateGroupingPerMoleculeType,
                           const real                       maxUpdateGroupRadius) :
    useUpdateGroups_(true),
    updateGroupingPerMoleculeType_(std::move(updateGroupingPerMoleculeType)),
    maxUpdateGroupRadius_(maxUpdateGroupRadius)
{
}

bool systemHasConstraintsOrVsites(const gmx_mtop_t& mtop)
{
    IListRange ilistRange(mtop);
    return std::any_of(ilistRange.begin(), ilistRange.end(), [](const auto& ilists) {
        return !extractILists(ilists.list(), IF_CONSTRAINT | IF_VSITE).empty();
    });
}

UpdateGroups makeUpdateGroups(const gmx::MDLogger&             mdlog,
                              std::vector<RangePartitioning>&& updateGroupingPerMoleculeType,
                              const real                       maxUpdateGroupRadius,
                              const bool                       doRerun,
                              const bool                       useDomainDecomposition,
                              const bool                       systemHasConstraintsOrVsites,
                              const real                       cutoffMargin)
{
    GMX_RELEASE_ASSERT(!updateGroupingPerMoleculeType.empty(), "We need the update grouping");

    MessageStringCollector messages;

    messages.startContext("When checking whether update groups are usable:");

    messages.appendIf(doRerun, "Rerun does not support update groups");

    messages.appendIf(!useDomainDecomposition,
                      "Domain decomposition is not active, so there is no need for update groups");

    messages.appendIf(!systemHasConstraintsOrVsites,
                      "No constraints or virtual sites are in use, so it is best not to use update "
                      "groups");

    messages.appendIf(
            getenv("GMX_NO_UPDATEGROUPS") != nullptr,
            "Environment variable GMX_NO_UPDATEGROUPS prohibited the use of update groups");

    // To use update groups, the large domain-to-domain cutoff
    // distance should be compatible with the box size.
    messages.appendIf(2 * maxUpdateGroupRadius >= cutoffMargin,
                      "The combination of rlist and box size prohibits the use of update groups");

    if (!messages.isEmpty())
    {
        // Log why we can't use update groups
        GMX_LOG(mdlog.info).appendText(messages.toString());
        return UpdateGroups();
    }

    // Success!
    return UpdateGroups(std::move(updateGroupingPerMoleculeType), maxUpdateGroupRadius);
}


} // namespace gmx
