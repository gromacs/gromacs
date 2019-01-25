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
    AtomTypeData(const t_atom  &a,
                 char         **name,
                 const t_param &nb,
                 const int      bondAtomType,
                 const int      atomNumber) :
        atom_(a), name_(name), nb_(nb),
        bondAtomType_(bondAtomType),
        atomNumber_(atomNumber)
    {}
    //! Actual atom data.
    t_atom   atom_;
    //! Atom name.
    char   **name_;
    //! Nonbonded data.
    t_param  nb_;
    //! Bonded atomtype for the type.
    int      bondAtomType_;
    //! Atom number for the atom type.
    int      atomNumber_;
};

class PreprocessingAtomType::Impl
{
    public:
        //! The number for currently loaded entries.
        int nr() const { return types.size(); }
        //! The actual atom type data.
        std::vector<AtomTypeData> types;
};


int PreprocessingAtomType::atomTypeFromString(const std::string &str) const
{
    /* Atom types are always case sensitive */
    for (auto it = impl_->types.begin(); it != impl_->types.end(); it++)
    {
        if (strcmp(str.c_str(), *it->name_) == 0)
        {
            return std::distance(impl_->types.begin(), it);
        }
    }

    return NOTSET;
}

int PreprocessingAtomType::nr() const
{
    return impl_->nr();
}

const char *PreprocessingAtomType::atomNameFromType(int nt) const
{
    if ((nt < 0) || (nt >= impl_->nr()))
    {
        return nullptr;
    }

    return *(impl_->types[nt].name_);
}

real PreprocessingAtomType::atomMassAFromType(int nt) const
{
    if ((nt < 0) || (nt >= impl_->nr()))
    {
        return NOTSET;
    }

    return impl_->types[nt].atom_.m;
}

real PreprocessingAtomType::atomMassBFromType(int nt) const
{
    if ((nt < 0) || (nt >= impl_->nr()))
    {
        return NOTSET;
    }

    return impl_->types[nt].atom_.mB;
}

real PreprocessingAtomType::atomChargeAFromType(int nt) const
{
    if ((nt < 0) || (nt >= impl_->nr()))
    {
        return NOTSET;
    }

    return impl_->types[nt].atom_.q;
}

real PreprocessingAtomType::atomChargeBFromType(int nt) const
{
    if ((nt < 0) || (nt >= impl_->nr()))
    {
        return NOTSET;
    }

    return impl_->types[nt].atom_.qB;
}

int PreprocessingAtomType::atomParameterFromType(int nt) const
{
    if ((nt < 0) || (nt >= impl_->nr()))
    {
        return NOTSET;
    }

    return impl_->types[nt].atom_.ptype;
}

int PreprocessingAtomType::bondAtomParameterFromType(int nt) const
{
    if ((nt < 0) || (nt >= impl_->nr()))
    {
        return NOTSET;
    }

    return impl_->types[nt].bondAtomType_;
}

int PreprocessingAtomType::atomNumberFromType(int nt) const
{
    if ((nt < 0) || (nt >= impl_->nr()))
    {
        return NOTSET;
    }

    return impl_->types[nt].atomNumber_;
}

real PreprocessingAtomType::atomParameter(int nt, int param) const
{
    if ((nt < 0) || (nt >= impl_->nr()))
    {
        return NOTSET;
    }
    if ((param < 0) || (param >= MAXFORCEPARAM))
    {
        return NOTSET;
    }
    return impl_->types[nt].nb_.c[param];
}

PreprocessingAtomType::PreprocessingAtomType()
    : impl_(new Impl)
{}

PreprocessingAtomType::~PreprocessingAtomType()
{}

int PreprocessingAtomType::addType(t_symtab          *tab,
                                   const t_atom      &a,
                                   const char        *name,
                                   const t_param     &nb,
                                   int                bondAtomType,
                                   int                atomNumber)
{
    auto found = std::find_if(impl_->types.begin(), impl_->types.end(),
                              [&name](const AtomTypeData &data)
                              { return strcmp(name, *data.name_) == 0; });

    if (found == impl_->types.end())
    {
        impl_->types.push_back(AtomTypeData(a,
                                            put_symtab(tab, name),
                                            nb,
                                            bondAtomType,
                                            atomNumber));
    }
    return impl_->nr();
}

int PreprocessingAtomType::setType(int                nt,
                                   t_symtab          *tab,
                                   const t_atom      &a,
                                   const char        *name,
                                   const t_param     &nb,
                                   int                bondAtomType,
                                   int                atomNumber)
{
    if ((nt < 0) || (nt >= impl_->nr()))
    {
        return NOTSET;
    }

    impl_->types[nt].atom_         = a;
    impl_->types[nt].name_         = put_symtab(tab, name);
    impl_->types[nt].nb_           = nb;
    impl_->types[nt].bondAtomType_ = bondAtomType;
    impl_->types[nt].atomNumber_   = atomNumber;

    return nt;
}

void PreprocessingAtomType::printTypes(FILE * out)
{
    fprintf (out, "[ %s ]\n", dir2str(Directive::d_atomtypes));
    fprintf (out, "; %6s  %8s  %8s  %8s  %12s  %12s\n",
             "type", "mass", "charge", "particle", "c6", "c12");
    for (auto &entry : impl_->types)
    {
        fprintf(out, "%8s  %8.3f  %8.3f  %8s  %12e  %12e\n",
                *(entry.name_), entry.atom_.m, entry.atom_.q, "A",
                entry.nb_.c0(), entry.nb_.c1());
    }

    fprintf (out, "\n");
}

static int search_atomtypes(const PreprocessingAtomType *ga,
                            int                         *n,
                            gmx::ArrayRef<int>           typelist,
                            int                          thistype,
                            t_param                      param[],
                            int                          ftype)
{
    int      i, nn, nrfp, j, k, ntype, tli;
    bool     bFound = FALSE;

    nn    = *n;
    nrfp  = NRFP(ftype);
    ntype = ga->nr();

    for (i = 0; (i < nn); i++)
    {
        if (typelist[i] == thistype)
        {
            /* This type number has already been added */
            break;
        }

        /* Otherwise, check if the parameters are identical to any previously added type */

        bFound = TRUE;
        for (j = 0; j < ntype && bFound; j++)
        {
            /* Check nonbonded parameters */
            for (k = 0; k < nrfp && bFound; k++)
            {
                bFound = (param[ntype*typelist[i]+j].c[k] == param[ntype*thistype+j].c[k]);
            }

            /* Check atomnumber */
            tli    = typelist[i];
            bFound = bFound &&
                (ga->atomNumberFromType(tli) == ga->atomNumberFromType(thistype));
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

void PreprocessingAtomType::renumberTypes(gmx::ArrayRef<SystemParameters> plist,
                                          gmx_mtop_t                     *mtop,
                                          int                            *wall_atomtype,
                                          bool                            bVerbose)
{
    int         nat, ftype, ntype;

    ntype = nr();
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
    if (plist[F_LJ].nr > 0)
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
    for (const gmx_moltype_t &moltype : mtop->moltype)
    {
        const t_atoms *atoms = &moltype.atoms;
        for (int i = 0; (i < atoms->nr); i++)
        {
            atoms->atom[i].type =
                search_atomtypes(this, &nat, typelist, atoms->atom[i].type,
                                 plist[ftype].param, ftype);
            atoms->atom[i].typeB =
                search_atomtypes(this, &nat, typelist, atoms->atom[i].typeB,
                                 plist[ftype].param, ftype);
        }
    }

    for (int i = 0; i < 2; i++)
    {
        if (wall_atomtype[i] >= 0)
        {
            wall_atomtype[i] = search_atomtypes(this, &nat, typelist, wall_atomtype[i],
                                                plist[ftype].param, ftype);
        }
    }

    std::vector<AtomTypeData> new_types;
    /* We now have a list of unique atomtypes in typelist */

    /* Renumber nlist */
    /* Renumber nlist */
    t_param *nbsnew = nullptr;
    snew(nbsnew, plist[ftype].nr);

    int nrfp  = NRFP(ftype);

    int k = 0;
    for (int i = 0; (i < nat); i++)
    {
        int mi = typelist[i];
        for (int j = 0; (j < nat); j++, k++)
        {
            int mj = typelist[j];
            for (int l = 0; (l < nrfp); l++)
            {
                nbsnew[k].c[l] = plist[ftype].param[ntype*mi+mj].c[l];
            }
            new_types.push_back(impl_->types[mi]);
        }
    }

    int i;
    for (i = 0; (i < nat*nat); i++)
    {
        for (int l = 0; (l < nrfp); l++)
        {
            plist[ftype].param[i].c[l] = nbsnew[i].c[l];
        }
    }
    plist[ftype].nr     = i;
    mtop->ffparams.atnr = nat;

    impl_->types = new_types;

    sfree(nbsnew);
}

void PreprocessingAtomType::copyTot_atomtypes(t_atomtypes *atomtypes) const
{
    /* Copy the atomtype data to the topology atomtype list */
    int ntype         = nr();
    atomtypes->nr = ntype;
    snew(atomtypes->atomnumber, ntype);

    for (int i = 0; i < ntype; i++)
    {
        atomtypes->atomnumber[i] = impl_->types[i].atomNumber_;
    }
}
