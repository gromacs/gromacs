/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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

#include "gpp_atomtype.h"

#include <climits>
#include <cmath>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <optional>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

struct AtomTypeData
{
    //! Explicit constructor.
    AtomTypeData(const t_atom&            a,
                 const std::string&       name,
                 const InteractionOfType& nb,
                 const int                bondAtomType,
                 const int                atomNumber) :
        atom_(a), name_(name), nb_(nb), bondAtomType_(bondAtomType), atomNumber_(atomNumber)
    {
    }
    //! Actual atom data.
    t_atom atom_;
    //! Atom name.
    std::string name_;
    //! Nonbonded data.
    InteractionOfType nb_;
    //! Bonded atomtype for the type.
    int bondAtomType_;
    //! Atom number for the atom type.
    int atomNumber_;
};

class PreprocessingAtomTypes::Impl
{
public:
    //! The number for currently loaded entries.
    size_t size() const { return types.size(); }
    //! The actual atom type data.
    std::vector<AtomTypeData> types;
    //! Map from \c types[i].name to \c i for quick look-up in \ref atomTypeFromName. Ref #3974.
    std::unordered_map<std::string, int> nameToAtomType;
};

bool PreprocessingAtomTypes::isSet(int nt) const
{
    return ((nt >= 0) && (nt < gmx::ssize(*this)));
}

std::optional<int> PreprocessingAtomTypes::atomTypeFromName(const std::string& str) const
{
    /* Atom types are always case sensitive */
    const auto found = impl_->nameToAtomType.find(str);
    if (found == impl_->nameToAtomType.end())
    {
        return std::nullopt;
    }
    else
    {
        GMX_ASSERT(str == impl_->types[found->second].name_,
                   "Invalid data in atomTypeFromName lookup table");
        return std::make_optional(found->second);
    }
}

size_t PreprocessingAtomTypes::size() const
{
    return impl_->size();
}

std::optional<const std::string> PreprocessingAtomTypes::atomNameFromAtomType(int nt) const
{
    return isSet(nt) ? std::make_optional(impl_->types[nt].name_) : std::nullopt;
}

std::optional<real> PreprocessingAtomTypes::atomMassFromAtomType(int nt) const
{
    return isSet(nt) ? std::make_optional(impl_->types[nt].atom_.m) : std::nullopt;
}

std::optional<real> PreprocessingAtomTypes::atomChargeFromAtomType(int nt) const
{
    return isSet(nt) ? std::make_optional(impl_->types[nt].atom_.q) : std::nullopt;
}

std::optional<ParticleType> PreprocessingAtomTypes::atomParticleTypeFromAtomType(int nt) const
{
    return isSet(nt) ? std::make_optional(impl_->types[nt].atom_.ptype) : std::nullopt;
}

std::optional<int> PreprocessingAtomTypes::bondAtomTypeFromAtomType(int nt) const
{
    return isSet(nt) ? std::make_optional(impl_->types[nt].bondAtomType_) : std::nullopt;
}

std::optional<int> PreprocessingAtomTypes::atomNumberFromAtomType(int nt) const
{
    return isSet(nt) ? std::make_optional(impl_->types[nt].atomNumber_) : std::nullopt;
}

std::optional<real> PreprocessingAtomTypes::atomNonBondedParamFromAtomType(int nt, int param) const
{
    if (!isSet(nt))
    {
        return std::nullopt;
    }
    gmx::ArrayRef<const real> forceParam = impl_->types[nt].nb_.forceParam();
    if ((param < 0) || (param >= MAXFORCEPARAM))
    {
        return std::nullopt;
    }
    return std::make_optional(forceParam[param]);
}

PreprocessingAtomTypes::PreprocessingAtomTypes() : impl_(new Impl) {}

PreprocessingAtomTypes::PreprocessingAtomTypes(PreprocessingAtomTypes&& old) noexcept :
    impl_(std::move(old.impl_))
{
}

PreprocessingAtomTypes& PreprocessingAtomTypes::operator=(PreprocessingAtomTypes&& old) noexcept
{
    impl_ = std::move(old.impl_);
    return *this;
}

PreprocessingAtomTypes::~PreprocessingAtomTypes() {}

int PreprocessingAtomTypes::addType(const t_atom&            a,
                                    const std::string&       name,
                                    const InteractionOfType& nb,
                                    int                      bondAtomType,
                                    int                      atomNumber)
{
    auto position = atomTypeFromName(name);
    if (!position.has_value())
    {
        impl_->types.emplace_back(a, name, nb, bondAtomType, atomNumber);
        const int newType           = impl_->types.size() - 1;
        impl_->nameToAtomType[name] = newType;
        return newType;
    }
    else
    {
        return *position;
    }
}

std::optional<int> PreprocessingAtomTypes::setType(int                      nt,
                                                   const t_atom&            a,
                                                   const std::string&       name,
                                                   const InteractionOfType& nb,
                                                   int                      bondAtomType,
                                                   int                      atomNumber)
{
    if (!isSet(nt))
    {
        return std::nullopt;
    }

    impl_->types[nt].atom_         = a;
    impl_->types[nt].name_         = name;
    impl_->types[nt].nb_           = nb;
    impl_->types[nt].bondAtomType_ = bondAtomType;
    impl_->types[nt].atomNumber_   = atomNumber;

    return std::make_optional(nt);
}

static int search_atomtypes(const PreprocessingAtomTypes*          ga,
                            int*                                   n,
                            gmx::ArrayRef<int>                     typelist,
                            int                                    thistype,
                            gmx::ArrayRef<const InteractionOfType> interactionTypes,
                            int                                    ftype)
{
    int nn    = *n;
    int nrfp  = NRFP(ftype);
    int ntype = ga->size();

    int i;
    for (i = 0; (i < nn); i++)
    {
        if (typelist[i] == thistype)
        {
            /* This type number has already been added */
            break;
        }

        /* Otherwise, check if the parameters are identical to any previously added type */

        bool bFound = true;
        for (int j = 0; j < ntype && bFound; j++)
        {
            /* Check nonbonded parameters */
            gmx::ArrayRef<const real> forceParam1 =
                    interactionTypes[ntype * typelist[i] + j].forceParam();
            gmx::ArrayRef<const real> forceParam2 = interactionTypes[ntype * thistype + j].forceParam();
            for (int k = 0; (k < nrfp) && bFound; k++)
            {
                bFound = forceParam1[k] == forceParam2[k];
            }

            /* Check atomnumber */
            int tli = typelist[i];
            bFound  = bFound
                     && (*ga->atomNumberFromAtomType(tli) == *ga->atomNumberFromAtomType(thistype));
        }
        if (bFound)
        {
            break;
        }
    }

    if (i == nn)
    {
        if (nn == ntype)
        {
            gmx_fatal(FARGS, "Atomtype horror n = %d, %s, %d", nn, __FILE__, __LINE__);
        }
        typelist[nn] = thistype;
        nn++;
    }
    *n = nn;

    return i;
}

void PreprocessingAtomTypes::renumberTypes(gmx::ArrayRef<InteractionsOfType> plist,
                                           gmx_mtop_t*                       mtop,
                                           int*                              wall_atomtype,
                                           bool                              bVerbose)
{
    int nat, ftype, ntype;

    ntype = size();
    std::vector<int> typelist(ntype);

    if (bVerbose)
    {
        fprintf(stderr, "renumbering atomtypes...\n");
    }

    /* Since the bonded interactions have been assigned now,
     * we want to reduce the number of atom types by merging
     * ones with identical nonbonded interactions, in addition
     * to removing unused ones.
     *
     * With QM/MM we also check that the atom numbers match
     */

    /* Get nonbonded interaction type */
    if (plist[F_LJ].size() > 0)
    {
        ftype = F_LJ;
    }
    else
    {
        ftype = F_BHAM;
    }

    /* Renumber atomtypes by first making a list of which ones are actually used.
     * We provide the list of nonbonded parameters so search_atomtypes
     * can determine if two types should be merged.
     */
    nat = 0;
    for (const gmx_moltype_t& moltype : mtop->moltype)
    {
        const t_atoms* atoms = &moltype.atoms;
        for (int i = 0; (i < atoms->nr); i++)
        {
            atoms->atom[i].type = search_atomtypes(
                    this, &nat, typelist, atoms->atom[i].type, plist[ftype].interactionTypes, ftype);
            atoms->atom[i].typeB = search_atomtypes(
                    this, &nat, typelist, atoms->atom[i].typeB, plist[ftype].interactionTypes, ftype);
        }
    }

    for (int i = 0; i < 2; i++)
    {
        if (wall_atomtype[i] >= 0)
        {
            wall_atomtype[i] = search_atomtypes(
                    this, &nat, typelist, wall_atomtype[i], plist[ftype].interactionTypes, ftype);
        }
    }

    std::vector<AtomTypeData> new_types;
    /* We now have a list of unique atomtypes in typelist */

    /* Renumber nlist */
    std::vector<InteractionOfType> nbsnew;

    // Reset the map used for fast lookups, and refill it below
    impl_->nameToAtomType.clear();

    for (int i = 0; (i < nat); i++)
    {
        int mi = typelist[i];
        for (int j = 0; (j < nat); j++)
        {
            int                      mj              = typelist[j];
            const InteractionOfType& interactionType = plist[ftype].interactionTypes[ntype * mi + mj];
            nbsnew.emplace_back(interactionType.atoms(),
                                interactionType.forceParam(),
                                interactionType.interactionTypeName());
        }
        new_types.push_back(impl_->types[mi]);
        impl_->nameToAtomType[std::string(impl_->types[mi].name_)] = new_types.size() - 1;
    }

    mtop->ffparams.atnr = nat;

    impl_->types                  = new_types;
    plist[ftype].interactionTypes = nbsnew;
}
