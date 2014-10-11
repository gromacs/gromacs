/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Utility classes for manipulating \c t_atoms structures.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 */
#ifndef GMX_TOPOLOGY_ATOMSBUILDER_H
#define GMX_TOPOLOGY_ATOMSBUILDER_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/real.h"

struct t_atoms;
struct t_resinfo;

namespace gmx
{

class AtomsBuilder
{
    public:
        explicit AtomsBuilder(t_atoms *atoms);
        ~AtomsBuilder();

        void reserve(int atomCount, int residueCount);
        void clearAtoms();

        int currentAtomCount() const;

        void setNextResidueNumber(int number);
        void addAtom(const t_atoms &atoms, int i);
        void startResidue(const t_resinfo &resinfo);
        void finishResidue(const t_resinfo &resinfo);
        void discardCurrentResidue();

        void mergeAtoms(const t_atoms &atoms);

    private:
        t_atoms *atoms_;
        int      nrAlloc_;
        int      nresAlloc_;
        int      currentResidueIndex_;
        int      nextResidueNumber_;

        GMX_DISALLOW_COPY_AND_ASSIGN(AtomsBuilder);
};

class AtomsRemover
{
    public:
        explicit AtomsRemover(const t_atoms &atoms);
        ~AtomsRemover();

        void markAll();
        void markResidue(const t_atoms &atoms, int atomIndex, bool bStatus);
        bool isMarked(int atomIndex) const { return removed_[atomIndex] != 0; }

        void removeMarkedVectors(rvec array[]) const;
        void removeMarkedValues(real array[]) const;
        void removeMarkedAtoms(t_atoms *atoms) const;

    private:
        std::vector<char> removed_;

        GMX_DISALLOW_COPY_AND_ASSIGN(AtomsRemover);
};

} // namespace gmx

#endif
