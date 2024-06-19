/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
/*!\file
 * \internal
 * \brief
 * Implements helper methods to place particle COM in boxes.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_pbcutil
 */
#include "gmxpre.h"

#include "gromacs/pbcutil/com.h"

#include <algorithm>
#include <iterator>
#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbcenums.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/range.h"
#include "gromacs/utility/real.h"

enum class PbcType : int;

namespace gmx
{

namespace
{

/*! \brief
 * Calculates shift to place COM into a box.
 *
 * \param[in] position     Unshifted COM position.
 * \param[in] box          The current box to place the COM in.
 * \param[in] pbcType      What kind of box is being used.
 * \param[in] unitCellType Type of unitcell being used.
 * \param[in] centerType   How things should be centered.
 * \returns The shift needed to place the COM into the box.
 */
RVec evaluateShiftToBox(const RVec&          position,
                        const matrix         box,
                        const PbcType&       pbcType,
                        const UnitCellType&  unitCellType,
                        const CenteringType& centerType)
{
    RVec newCOMPosition(position);
    auto comArrayRef = arrayRefFromArray(&newCOMPosition, 1);
    switch (unitCellType)
    {
        case (UnitCellType::Rectangular): put_atoms_in_box(pbcType, box, comArrayRef); break;
        case (UnitCellType::Triclinic):
            put_atoms_in_triclinic_unitcell(static_cast<int>(centerType), box, comArrayRef);
            break;
        case (UnitCellType::Compact):
            put_atoms_in_compact_unitcell(pbcType, static_cast<int>(centerType), box, comArrayRef);
            break;
        default: GMX_RELEASE_ASSERT(false, "Unhandled type of unit cell");
    }
    RVec shift(0, 0, 0);
    rvec_sub(newCOMPosition, position, shift);
    return shift;
}

/*! \brief
 * Calculates the COM for each collection of atoms.
 *
 * \param[in] x       View on coordinates of the molecule.
 * \param[in] moltype Which molecule type to calculate for.
 * \param[in] atomOffset If needed, point from where to count the first atom to process.
 * \returns The center of mass for the molecule.
 */
RVec calculateCOM(ArrayRef<const RVec> x, const gmx_moltype_t& moltype, const int atomOffset = 0)
{
    RVec   com(0, 0, 0);
    double totalMass   = 0;
    int    currentAtom = atomOffset;
    for (const auto& coord : x)
    {
        real mass = moltype.atoms.atom[currentAtom].m;
        for (int d = 0; d < DIM; d++)
        {
            com[d] += mass * coord[d];
        }
        totalMass += mass;
        currentAtom++;
    }
    svmul(1.0 / totalMass, com, com);
    return com;
}

} // namespace

void shiftAtoms(const RVec& shift, ArrayRef<RVec> x)
{
    std::transform(
            std::begin(x), std::end(x), std::begin(x), [shift](RVec elemX) { return elemX + shift; });
}

void placeCoordinatesWithCOMInBox(const PbcType&      pbcType,
                                  const UnitCellType  unitCellType,
                                  const CenteringType centerType,
                                  const matrix        box,
                                  ArrayRef<RVec>      x,
                                  const gmx_mtop_t&   mtop,
                                  const COMShiftType  comShiftType)
{
    GMX_RELEASE_ASSERT(comShiftType != COMShiftType::Count, "Using COUNT of enumeration");
    // loop over all molecule blocks, then over all molecules in these blocks
    int molb = 0;
    for (const auto& molblock : mtop.molblock)
    {
        const MoleculeBlockIndices& ind              = mtop.moleculeBlockIndices[molb];
        const gmx_moltype_t&        moltype          = mtop.moltype[molblock.type];
        const int                   atomsPerMolecule = ind.numAtomsPerMolecule;
        int                         atomStart        = ind.globalAtomStart;
        for (int mol = 0; mol < molblock.nmol; mol++)
        {
            std::vector<Range<int>> atomRanges;
            if (comShiftType == COMShiftType::Molecule)
            {
                atomRanges.emplace_back(0, atomsPerMolecule);
            }
            else if (comShiftType == COMShiftType::Residue)
            {
                atomRanges = atomRangeOfEachResidue(moltype);
            }
            for (const auto& atomRange : atomRanges)
            {
                const int      firstAtomGlobalPosition = atomStart + *atomRange.begin();
                ArrayRef<RVec> coords = x.subArray(firstAtomGlobalPosition, atomRange.size());
                const RVec     com    = calculateCOM(coords, moltype, *atomRange.begin());
                const RVec shift = evaluateShiftToBox(com, box, pbcType, unitCellType, centerType);
                shiftAtoms(shift, coords);
            }
            atomStart += atomsPerMolecule;
        }
        molb++;
    }
}

} // namespace gmx
