/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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

#include "gpp_atomtype.h"

#include <climits>
#include <cmath>
#include <cstring>

#include <algorithm>

#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

struct AtomTypeData
{
    //! Explicit constructor.
    AtomTypeData(const t_atom& a, char** name, const InteractionOfType& nb, const int bondAtomType, const int atomNumber) :
        atom_(a),
        name_(name),
        nb_(nb),
        bondAtomType_(bondAtomType),
        atomNumber_(atomNumber)
    {
    }
    //! Actual atom data.
    t_atom atom_;
    //! Atom name.
    char** name_;
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
};

bool PreprocessingAtomTypes::isSet(int nt) const
{
    return ((nt >= 0) && (nt < gmx::ssize(*this)));
}

int PreprocessingAtomTypes::atomTypeFromName(const std::string& str) const
{
    /* Atom types are always case sensitive */
    auto found = std::find_if(impl_->types.begin(), impl_->types.end(),
                              [&str](const auto& type) { return str == *type.name_; });
    if (found == impl_->types.end())
    {
        return NOTSET;
    }
    else
    {
        return std::distance(impl_->types.begin(), found);
    }
}

size_t PreprocessingAtomTypes::size() const
{
    return impl_->size();
}

const char* PreprocessingAtomTypes::atomNameFromAtomType(int nt) const
{
    return isSet(nt) ? *(impl_->types[nt].name_) : nullptr;
}

real PreprocessingAtomTypes::atomMassFromAtomType(int nt) const
{
    return isSet(nt) ? impl_->types[nt].atom_.m : NOTSET;
}

real PreprocessingAtomTypes::atomChargeFromAtomType(int nt) const
{
    return isSet(nt) ? impl_->types[nt].atom_.q : NOTSET;
}

int PreprocessingAtomTypes::atomParticleTypeFromAtomType(int nt) const
{
    return isSet(nt) ? impl_->types[nt].atom_.ptype : NOTSET;
}

int PreprocessingAtomTypes::bondAtomTypeFromAtomType(int nt) const
{
    return isSet(nt) ? impl_->types[nt].bondAtomType_ : NOTSET;
}

int PreprocessingAtomTypes::atomNumberFromAtomType(int nt) const
{
    return isSet(nt) ? impl_->types[nt].atomNumber_ : NOTSET;
}

real PreprocessingAtomTypes::atomNonBondedParamFromAtomType(int nt, int param) const
{
    if (!isSet(nt))
    {
        return NOTSET;
    }
    gmx::ArrayRef<const real> forceParam = impl_->types[nt].nb_.forceParam();
    if ((param < 0) || (param >= MAXFORCEPARAM))
    {
        return NOTSET;
    }
    return forceParam[param];
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

int PreprocessingAtomTypes::addType(t_symtab*                tab,
                                    const t_atom&            a,
                                    const std::string&       name,
                                    const InteractionOfType& nb,
                                    int                      bondAtomType,
                                    int                      atomNumber)
{
    int position = atomTypeFromName(name);
    if (position == NOTSET)
    {
        impl_->types.emplace_back(a, put_symtab(tab, name.c_str()), nb, bondAtomType, atomNumber);
        return atomTypeFromName(name);
    }
    else
    {
        return position;
    }
}

int PreprocessingAtomTypes::setType(int                      nt,
                                    t_symtab*                tab,
                                    const t_atom&            a,
                                    const std::string&       name,
                                    const InteractionOfType& nb,
                                    int                      bondAtomType,
                                    int                      atomNumber)
{
    if (!isSet(nt))
    {
        return NOTSET;
    }

    impl_->types[nt].atom_         = a;
    impl_->types[nt].name_         = put_symtab(tab, name.c_str());
    impl_->types[nt].nb_           = nb;
    impl_->types[nt].bondAtomType_ = bondAtomType;
    impl_->types[nt].atomNumber_   = atomNumber;

    return nt;
}

void PreprocessingAtomTypes::printTypes(FILE* out)
{
    fprintf(out, "[ %s ]\n", dir2str(Directive::d_atomtypes));
    fprintf(out, "; %6s  %8s  %8s  %8s  %12s  %12s\n", "type", "mass", "charge", "particle", "c6",
            "c12");
    for (auto& entry : impl_->types)
    {
        fprintf(out, "%8s  %8.3f  %8.3f  %8s  %12e  %12e\n", *(entry.name_), entry.atom_.m,
                entry.atom_.q, "A", entry.nb_.c0(), entry.nb_.c1());
    }

    fprintf(out, "\n");
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
            bFound = bFound && (ga->atomNumberFromAtomType(tli) == ga->atomNumberFromAtomType(thistype));
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
            atoms->atom[i].type  = search_atomtypes(this, &nat, typelist, atoms->atom[i].type,
                                                   plist[ftype].interactionTypes, ftype);
            atoms->atom[i].typeB = search_atomtypes(this, &nat, typelist, atoms->atom[i].typeB,
                                                    plist[ftype].interactionTypes, ftype);
        }
    }

    for (int i = 0; i < 2; i++)
    {
        if (wall_atomtype[i] >= 0)
        {
            wall_atomtype[i] = search_atomtypes(this, &nat, typelist, wall_atomtype[i],
                                                plist[ftype].interactionTypes, ftype);
        }
    }

    std::vector<AtomTypeData> new_types;
    /* We now have a list of unique atomtypes in typelist */

    /* Renumber nlist */
    std::vector<InteractionOfType> nbsnew;

    for (int i = 0; (i < nat); i++)
    {
        int mi = typelist[i];
        for (int j = 0; (j < nat); j++)
        {
            int                      mj              = typelist[j];
            const InteractionOfType& interactionType = plist[ftype].interactionTypes[ntype * mi + mj];
            nbsnew.emplace_back(interactionType.atoms(), interactionType.forceParam(),
                                interactionType.interactionTypeName());
        }
        new_types.push_back(impl_->types[mi]);
    }

    mtop->ffparams.atnr = nat;

    impl_->types                  = new_types;
    plist[ftype].interactionTypes = nbsnew;
}

void PreprocessingAtomTypes::copyTot_atomtypes(t_atomtypes* atomtypes) const
{
    /* Copy the atomtype data to the topology atomtype list */
    int ntype     = size();
    atomtypes->nr = ntype;
    snew(atomtypes->atomnumber, ntype);

    for (int i = 0; i < ntype; i++)
    {
        atomtypes->atomnumber[i] = impl_->types[i].atomNumber_;
    }
}
