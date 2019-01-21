/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements classes from atomsbuilder.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 */
#include "gmxpre.h"

#include "atomsbuilder.h"

#include <algorithm>

#include "gromacs/math/vec.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

/********************************************************************
 * AtomsBuilder
 */

AtomsBuilder::AtomsBuilder(std::vector<AtomInfo> *atoms,
       std::vector<Residue>  *resinfo,
       std::vector<PdbEntry> *pdb, t_symtab *symtab)
    : atoms_(*atoms), resinfo_(*resinfo), pdb_(*pdb), symtab_(symtab),
      currentResidueIndex_(resinfo->size()), nextResidueNumber_(-1)
{
    if (!resinfo->empty())
    {
        nextResidueNumber_ = resinfo->back().nr_ + 1;
    }
}

AtomsBuilder::~AtomsBuilder()
{
}

char **AtomsBuilder::symtabString(char **source)
{
    if (symtab_ != nullptr)
    {
        return put_symtab(symtab_, *source);
    }
    return source;
}

void AtomsBuilder::clearAtoms()
{
    atoms_.clear();
    resinfo_.clear();
    currentResidueIndex_ = 0;
    nextResidueNumber_   = -1;
}

int AtomsBuilder::currentAtomCount() const
{
    return atoms_.size();
}

void AtomsBuilder::setNextResidueNumber(int number)
{
    nextResidueNumber_ = number;
}

void AtomsBuilder::addAtom(const AtomInfo &atom, const PdbEntry &pdb)
{
    atoms_.push_back(atom);
    atoms_.back().resind_ = currentResidueIndex_;
    if (allAtomsHavePdbInfo(pdb_))
    {
        if (pdb.isSet_)
        {
            pdb_.push_back(pdb);
        }
        else
        {
            pdb_.push_back(PdbEntry());
        }
    }
}

void AtomsBuilder::startResidue(const Residue &resinfo)
{
    if (nextResidueNumber_ == -1)
    {
        nextResidueNumber_ = resinfo.nr_;
    }
    const int index = resinfo_.size();
    resinfo_.push_back(resinfo);
    resinfo_.back().nr_ = nextResidueNumber_;
    ++nextResidueNumber_;
    currentResidueIndex_      = index;
}

void AtomsBuilder::finishResidue(const Residue &resinfo)
{
    if (nextResidueNumber_ == -1)
    {
        nextResidueNumber_ = resinfo.nr_;
    }
    const int index = currentResidueIndex_;
    resinfo_[index]      = resinfo;
    resinfo_[index].nr_   = nextResidueNumber_;
    resinfo_[index].name_ = symtabString(resinfo.name_);
    ++nextResidueNumber_;
    currentResidueIndex_      = index + 1;
}

void AtomsBuilder::discardCurrentResidue()
{
    int index = atoms_.size() - 1;
    while (index > 0 && atoms_.back().resind_ == currentResidueIndex_)
    {
        atoms_.erase(atoms_.end() -1);
    }
}

void AtomsBuilder::mergeAtoms(gmx::ArrayRef<const AtomInfo> atoms,
                              gmx::ArrayRef<const Residue> resinfo,
                              gmx::ArrayRef<const PdbEntry> pdb)
{
    int prevResInd = -1;
    if (!pdb.empty() && (atoms.size() != pdb.size()))
    {
        GMX_THROW(InternalError("Size of atom and pdb entries are not the same"));
    }

    for (int i = 0; i < atoms.size(); ++i)
    {
        const int resind = atoms[i].resind_;
        if (resind != prevResInd)
        {
            startResidue(resinfo[resind]);
            prevResInd = resind;
        }
        addAtom(atoms[i], pdb[i]);
    }
}

/********************************************************************
 * AtomsRemover
 */

AtomsRemover::AtomsRemover(int size)
    : removed_(size, 0)
{
}

AtomsRemover::~AtomsRemover()
{
}

void AtomsRemover::refreshAtomCount(gmx::ArrayRef<const AtomInfo> atoms)
{
    removed_.resize(atoms.size(), 0);
}

void AtomsRemover::markAll()
{
    std::fill(removed_.begin(), removed_.end(), 1);
}

void AtomsRemover::markResidue(gmx::ArrayRef<const AtomInfo> atoms, int atomIndex, bool bStatus)
{
    const int resind = atoms[atomIndex].resind_;
    while (atomIndex > 0 && resind == atoms[atomIndex - 1].resind_)
    {
        --atomIndex;
    }
    while (atomIndex < atoms.size() && resind == atoms[atomIndex].resind_)
    {
        removed_[atomIndex] = (bStatus ? 1 : 0);
        ++atomIndex;
    }
}

void AtomsRemover::removeMarkedElements(std::vector<RVec> *container) const
{
    GMX_RELEASE_ASSERT(container->size() == removed_.size(),
                       "Mismatching contained passed for removing values");
    int j = 0;
    for (size_t i = 0; i < removed_.size(); ++i)
    {
        if (!removed_[i])
        {
            (*container)[j] = (*container)[i];
            ++j;
        }
    }
    container->resize(j);
}

void AtomsRemover::removeMarkedElements(std::vector<real> *container) const
{
    GMX_RELEASE_ASSERT(container->size() == removed_.size(),
                       "Mismatching contained passed for removing values");
    int j = 0;
    for (size_t i = 0; i < removed_.size(); ++i)
    {
        if (!removed_[i])
        {
            (*container)[j] = (*container)[i];
            ++j;
        }
    }
    container->resize(j);
}

void AtomsRemover::removeMarkedAtoms(std::vector<AtomInfo> *atoms,
                                     std::vector<Residue>  *resinfo,
                                     std::vector<PdbEntry> *pdb) const
{
    const int    originalAtomCount = atoms->size();
    AtomsBuilder builder(atoms, resinfo, pdb, nullptr);
    if (atoms->size() > 0)
    {
        builder.setNextResidueNumber(resinfo->at(0).nr_);
    }
    builder.clearAtoms();
    if (!pdb->empty() && (atoms->size() != pdb->size()))
    {
        GMX_THROW(InternalError("Size of atom and pdb entries are not the same"));
    }
    int          prevResInd = -1;
    for (int i = 0; i < originalAtomCount; ++i)
    {
        if (!removed_[i])
        {
            const int resind = atoms->at(i).resind_;
            if (resind != prevResInd)
            {
                builder.startResidue(resinfo->at(resind));
                prevResInd = resind;
            }
            PdbEntry newPdb;
            if (!pdb->empty())
            {
                newPdb = pdb->at(i);
            }
            builder.addAtom(atoms->at(i), newPdb);
        }
    }
}

} // namespace gmx
