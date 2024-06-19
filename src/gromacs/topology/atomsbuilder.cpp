/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
 * \brief
 * Implements classes from atomsbuilder.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/topology/atomsbuilder.h"

#include <cstddef>

#include <algorithm>
#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

/********************************************************************
 * AtomsBuilder
 */

AtomsBuilder::AtomsBuilder(t_atoms* atoms, t_symtab* symtab) :
    atoms_(atoms),
    symtab_(symtab),
    nrAlloc_(atoms->nr),
    nresAlloc_(atoms->nres),
    currentResidueIndex_(atoms->nres),
    nextResidueNumber_(-1)
{
    if (atoms->nres > 0)
    {
        nextResidueNumber_ = atoms->resinfo[atoms->nres - 1].nr + 1;
    }
}

AtomsBuilder::~AtomsBuilder() {}

char** AtomsBuilder::symtabString(char** source)
{
    if (symtab_ != nullptr)
    {
        return put_symtab(symtab_, *source);
    }
    return source;
}

void AtomsBuilder::reserve(int atomCount, int residueCount)
{
    srenew(atoms_->atom, atomCount);
    srenew(atoms_->atomname, atomCount);
    srenew(atoms_->resinfo, residueCount);
    if (atoms_->pdbinfo != nullptr)
    {
        srenew(atoms_->pdbinfo, atomCount);
    }
    nrAlloc_   = atomCount;
    nresAlloc_ = residueCount;
}

void AtomsBuilder::clearAtoms()
{
    atoms_->nr           = 0;
    atoms_->nres         = 0;
    currentResidueIndex_ = 0;
    nextResidueNumber_   = -1;
}

int AtomsBuilder::currentAtomCount() const
{
    return atoms_->nr;
}

void AtomsBuilder::setNextResidueNumber(int number)
{
    nextResidueNumber_ = number;
}

void AtomsBuilder::addAtom(const t_atoms& atoms, int i)
{
    const int index            = atoms_->nr;
    atoms_->atom[index]        = atoms.atom[i];
    atoms_->atomname[index]    = symtabString(atoms.atomname[i]);
    atoms_->atom[index].resind = currentResidueIndex_;
    if (atoms_->pdbinfo != nullptr)
    {
        if (atoms.pdbinfo != nullptr)
        {
            atoms_->pdbinfo[index] = atoms.pdbinfo[i];
        }
        else
        {
            gmx_pdbinfo_init_default(&atoms_->pdbinfo[index]);
        }
    }
    ++atoms_->nr;
}

void AtomsBuilder::startResidue(const t_resinfo& resinfo)
{
    if (nextResidueNumber_ == -1)
    {
        nextResidueNumber_ = resinfo.nr;
    }
    const int index             = atoms_->nres;
    atoms_->resinfo[index]      = resinfo;
    atoms_->resinfo[index].nr   = nextResidueNumber_;
    atoms_->resinfo[index].name = symtabString(resinfo.name);
    ++nextResidueNumber_;
    currentResidueIndex_ = index;
    ++atoms_->nres;
}

void AtomsBuilder::finishResidue(const t_resinfo& resinfo)
{
    if (nextResidueNumber_ == -1)
    {
        nextResidueNumber_ = resinfo.nr;
    }
    const int index             = currentResidueIndex_;
    atoms_->resinfo[index]      = resinfo;
    atoms_->resinfo[index].nr   = nextResidueNumber_;
    atoms_->resinfo[index].name = symtabString(resinfo.name);
    ++nextResidueNumber_;
    currentResidueIndex_ = index + 1;
    if (index >= atoms_->nres)
    {
        ++atoms_->nres;
    }
}

void AtomsBuilder::discardCurrentResidue()
{
    int index = atoms_->nr - 1;
    while (index > 0 && atoms_->atom[index - 1].resind == currentResidueIndex_)
    {
        --index;
    }
    atoms_->nr   = index;
    atoms_->nres = currentResidueIndex_;
}

void AtomsBuilder::mergeAtoms(const t_atoms& atoms)
{
    if (atoms_->nr + atoms.nr > nrAlloc_ || atoms_->nres + atoms.nres > nresAlloc_)
    {
        reserve(atoms_->nr + atoms.nr, atoms_->nres + atoms.nres);
    }
    int prevResInd = -1;
    for (int i = 0; i < atoms.nr; ++i)
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

AtomsRemover::AtomsRemover(const t_atoms& atoms) : removed_(atoms.nr, 0) {}

AtomsRemover::~AtomsRemover() {}

void AtomsRemover::refreshAtomCount(const t_atoms& atoms)
{
    removed_.resize(atoms.nr, 0);
}

void AtomsRemover::markAll()
{
    std::fill(removed_.begin(), removed_.end(), 1);
}

void AtomsRemover::markResidue(const t_atoms& atoms, int atomIndex, bool bStatus)
{
    const int resind = atoms.atom[atomIndex].resind;
    while (atomIndex > 0 && resind == atoms.atom[atomIndex - 1].resind)
    {
        --atomIndex;
    }
    while (atomIndex < atoms.nr && resind == atoms.atom[atomIndex].resind)
    {
        removed_[atomIndex] = (bStatus ? 1 : 0);
        ++atomIndex;
    }
}

void AtomsRemover::removeMarkedElements(std::vector<RVec>* container) const
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

void AtomsRemover::removeMarkedElements(std::vector<real>* container) const
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

void AtomsRemover::removeMarkedAtoms(t_atoms* atoms) const
{
    const int    originalAtomCount = atoms->nr;
    AtomsBuilder builder(atoms, nullptr);
    if (atoms->nr > 0)
    {
        builder.setNextResidueNumber(atoms->resinfo[0].nr);
    }
    builder.clearAtoms();
    int prevResInd = -1;
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
