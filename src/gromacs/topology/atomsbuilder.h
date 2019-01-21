/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/real.h"

struct AtomInfo;
struct Residue;
struct PdbEntry;
struct t_symtab;

namespace gmx
{

class AtomsBuilder
{
    public:
        AtomsBuilder(std::vector<AtomInfo> *atoms,
                     std::vector<Residue> *resinfo,
                     std::vector<PdbEntry> *pdb,
                     t_symtab *symtab);
        ~AtomsBuilder();

        void clearAtoms();

        int currentAtomCount() const;

        void setNextResidueNumber(int number);
        void addAtom(const AtomInfo &atom, const PdbEntry &pdb);
        void startResidue(const Residue &resinfo);
        void finishResidue(const Residue &resinfo);
        void discardCurrentResidue();

        void mergeAtoms(gmx::ArrayRef<const AtomInfo> atoms,
                        gmx::ArrayRef<const Residue> resinfo,
                        gmx::ArrayRef<const PdbEntry> pdb);

        gmx::ArrayRef<const AtomInfo> atoms() const { return atoms_; }
        gmx::ArrayRef<const Residue> resinfo() const { return resinfo_; }
        gmx::ArrayRef<const PdbEntry> pdb() const { return pdb_; }

    private:
        char **symtabString(char **source);

        std::vector<AtomInfo> atoms_;
        std::vector<Residue> resinfo_;
        std::vector<PdbEntry> pdb_;
        t_symtab *symtab_;
        int       nrAlloc_;
        int       nresAlloc_;
        int       currentResidueIndex_;
        int       nextResidueNumber_;

        GMX_DISALLOW_COPY_AND_ASSIGN(AtomsBuilder);
};

class AtomsRemover
{
    public:
        explicit AtomsRemover(int size);
        ~AtomsRemover();

        void refreshAtomCount(gmx::ArrayRef<const AtomInfo> atoms);

        void markAll();
        void markResidue(gmx::ArrayRef<const AtomInfo> atoms, int atomIndex, bool bStatus);
        bool isMarked(int atomIndex) const { return removed_[atomIndex] != 0; }

        void removeMarkedElements(std::vector<RVec> *container) const;
        void removeMarkedElements(std::vector<real> *container) const;
        void removeMarkedAtoms(std::vector<AtomInfo> *atoms,
                               std::vector<Residue>  *resinfo,
                               std::vector<PdbEntry> *pdb) const;

    private:
        std::vector<char> removed_;

        GMX_DISALLOW_COPY_AND_ASSIGN(AtomsRemover);
};

} // namespace gmx

#endif
