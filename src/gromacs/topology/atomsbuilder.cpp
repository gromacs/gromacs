/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

/********************************************************************
 * AtomsBuilder
 */

AtomsBuilder::AtomsBuilder(t_atoms *atoms, SymbolTable *symtab)
    : atoms_(atoms), symtab_(symtab),
      nrAlloc_(atoms->getNatoms()), nresAlloc_(atoms->getNresidues()),
      currentResidueIndex_(atoms->getNresidues()), nextResidueNumber_(-1)
{
    if (atoms->getNresidues() > 0)
    {
        nextResidueNumber_ = atoms->resinfo[atoms->getNresidues() - 1].nr + 1;
    }
}

AtomsBuilder::~AtomsBuilder()
{
}

SymbolPtr AtomsBuilder::symtabString(const char *source)
{
    GMX_RELEASE_ASSERT(symtab_, "Need valid symtab to add entry");

    return put_symtab(symtab_, source);
}

void AtomsBuilder::reserve(int atomCount, int residueCount)
{
    /*
       atoms_->atom.resize(atomCount);
       atoms_->atomname.resize(atomCount);
       atoms_->resinfo.resize(residueCount);
       if (!atoms_->pdbinfo.empty())
       {
        atoms_->pdbinfo.resize(atomCount);
       }
     */
    nrAlloc_   = atomCount;
    nresAlloc_ = residueCount;
}

void AtomsBuilder::clearAtoms()
{
    atoms_->atom.empty();
    atoms_->resinfo.empty();
    currentResidueIndex_ = 0;
    nextResidueNumber_   = -1;
}

int AtomsBuilder::currentAtomCount() const
{
    return atoms_->getNatoms();
}

void AtomsBuilder::setNextResidueNumber(int number)
{
    nextResidueNumber_ = number;
}

void AtomsBuilder::addAtom(const t_atoms &atoms, int i)
{
    atoms_->atom.emplace_back(atoms.atom[i]);
    atoms_->atomname.emplace_back(symtabString(atoms.atomname[i]->c_str()));
    atoms_->atom.end()->resind = currentResidueIndex_;
    if (!atoms_->pdbinfo.empty())
    {
        if (!atoms.pdbinfo.empty())
        {
            atoms_->pdbinfo.emplace_back(atoms.pdbinfo[i]);
        }
        else
        {
            t_pdbinfo newPdbinfo;
            gmx_pdbinfo_init_default(&newPdbinfo);
            atoms_->pdbinfo.emplace_back(newPdbinfo);
        }
    }
}

void AtomsBuilder::startResidue(const t_resinfo &resinfo)
{
    if (nextResidueNumber_ == -1)
    {
        nextResidueNumber_ = resinfo.nr;
    }
    const int index = atoms_->getNresidues();
    atoms_->resinfo.emplace_back(resinfo);
    atoms_->resinfo.end()->nr   = nextResidueNumber_;
    atoms_->resinfo.end()->name = symtabString(resinfo.name->c_str());
    ++nextResidueNumber_;
    currentResidueIndex_      = index;
}

void AtomsBuilder::finishResidue(const t_resinfo &resinfo)
{
    if (nextResidueNumber_ == -1)
    {
        nextResidueNumber_ = resinfo.nr;
    }
    const int index = currentResidueIndex_;
    atoms_->resinfo[index]      = resinfo;
    atoms_->resinfo[index].nr   = nextResidueNumber_;
    atoms_->resinfo[index].name = symtabString(resinfo.name->c_str());
    ++nextResidueNumber_;
    currentResidueIndex_      = index + 1;
//    if (index >= atoms_->nres)
//    {
//        ++atoms_->nres;
//    }
}

void AtomsBuilder::discardCurrentResidue()
{
    int index = atoms_->getNatoms() - 1;
    while (index > 0 && atoms_->atom[index - 1].resind == currentResidueIndex_)
    {
        atoms_->atom.erase(atoms_->atom.end());
        --index;
    }
    atoms_->resinfo.erase(atoms_->resinfo.end());
}

void AtomsBuilder::mergeAtoms(const t_atoms &atoms)
{
    int prevResInd = -1;
    for (int i = 0; i < atoms.getNatoms(); ++i)
    {
        const int resind = atoms.atom[i].resind;
        if (resind != prevResInd)
        {
            startResidue(atoms.resinfo[resind]);
            prevResInd = resind;
        }
        addAtom(atoms, i);
    }
}

/********************************************************************
 * AtomsRemover
 */

AtomsRemover::AtomsRemover(const t_atoms &atoms)
    : removed_(atoms.getNatoms(), 0)
{
}

AtomsRemover::~AtomsRemover()
{
}

void AtomsRemover::refreshAtomCount(const t_atoms &atoms)
{
    removed_.resize(atoms.getNatoms(), 0);
}

void AtomsRemover::markAll()
{
    std::fill(removed_.begin(), removed_.end(), 1);
}

void AtomsRemover::markResidue(const t_atoms &atoms, int atomIndex, bool bStatus)
{
    const int resind = atoms.atom[atomIndex].resind;
    while (atomIndex > 0 && resind == atoms.atom[atomIndex - 1].resind)
    {
        --atomIndex;
    }
    while (atomIndex < atoms.getNatoms() && resind == atoms.atom[atomIndex].resind)
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

void AtomsRemover::removeMarkedAtoms(t_atoms *atoms) const
{
    const int    originalAtomCount = atoms->getNatoms();
    AtomsBuilder builder(atoms, nullptr);
    if (atoms->getNatoms() > 0)
    {
        builder.setNextResidueNumber(atoms->resinfo[0].nr);
    }
    builder.clearAtoms();
    int          prevResInd = -1;
    for (int i = 0; i < originalAtomCount; ++i)
    {
        if (!removed_[i])
        {
            const int resind = atoms->atom[i].resind;
            if (resind != prevResInd)
            {
                builder.startResidue(atoms->resinfo[resind]);
                prevResInd = resind;
            }
            builder.addAtom(*atoms, i);
        }
    }
}

} // namespace gmx
