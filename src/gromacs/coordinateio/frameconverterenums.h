/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
/*! \file
 * \brief
 * Defines enum class defining the guarantees provided by different frameconverters
 * for the coordiante file manipulations done by them.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \libinternal
 * \ingroup module_coordinateio
 */
#ifndef GMX_COORDINATEIO_FRAMECONVERTERENUMS_H
#define GMX_COORDINATEIO_FRAMECONVERTERENUMS_H

namespace gmx
{

/*!\brief
 * The enums here define the guarantees provided by frameconverters concerning the
 * modifications they provide.
 *
 * A method can specify different kind of guarantees for the operation, and aggregate
 * methods can combine those flags to provide an overview about what has happened to a
 * coordinate frame during processing.
 *
 * \libinternal
 * \ingroup module_coordinateio
 *
 */
enum class FrameConverterFlags : unsigned long
{
    /*! \brief
     * Base setting means no guarantee is done by a module.
     */
    NoGuarantee = 1 << 0,
    /*! \brief
     * Tells us that molecules have been made whole.
     */
    MoleculesAreWhole = 1 << 1,
    /*! \brief
     * Tells us that no jumps over periodic boundary are present. Implies molecules are made whole.
     */
    NoPBCJumps = 1 << 2,
    /*! \brief
     * Tells us that COM of all molecules is in the box.
     */
    MoleculeCOMInBox = 1 << 3,
    /*! \brief
     * Tells us that COM of residues is in the box.
     */
    ResidueCOMInBox = 1 << 4,
    /*! \brief
     * Tells us that all atoms are in the box.
     */
    AtomsInBox = 1 << 5,
    /*! \brief
     * Tells us a converter changed the unit cell to rectangular.
     */
    UnitCellIsRectangular = 1 << 6,
    /*! \brief
     * Tells us a converter changed the unit cell to triclinic.
     */
    UnitCellIsTriclinic = 1 << 7,
    /*! \brief
     * Tells us a converter changed the unit cell to compact.
     */
    UnitCellIsCompact = 1 << 8,
    /*! \brief
     * Tells us that converter centered system in a box.
     *
     * Invalidated by calling a method that changes the type of box.
     */
    SystemIsCenteredInBox = 1 << 9,
    /*! \brief
     * Tells us that converter has fit coordinate data to reference.
     *
     * Specific for rotational and translational fit.
     */
    FitToReferenceRotTrans = 1 << 10,
    /*! \brief
     * Tells us that converter fit coordinate data to reference.
     *
     * Specific for rotation and translation in XY plane.
     */
    FitToReferenceRotTransXY = 1 << 11,
    /*! \brief
     * Tells us that converter fit coordinate data to reference.
     *
     * Specific for translational fit.
     */
    FitToReferenceTranslation = 1 << 12,
    /*! \brief
     * Tells us that converter fit coordinate data to reference.
     *
     * Specific for translational fit in XY plane.
     */
    FitToReferenceTranslationXY = 1 << 13,
    /*! \brief
     * Tells us that converter fit coordinate data to reference.
     *
     * Specific for progressive fit.
     */
    FitToReferenceProgressive = 1 << 14,
    /*! \brief
     * Tells us that converter has set a new center for the system.
     *
     * This affects the routines that place atoms in the box and the unit cell changing routines,
     * so they have to be invalidated by any method setting this flag.
     */
    NewSystemCenter = 1 << 15,
    //! Final entry to get number of operations.
    Count
};

//! Conversion of flag to its corresponding unsigned long value.
inline unsigned long convertFlag(FrameConverterFlags flag)
{
    return static_cast<unsigned long>(flag);
}

} // namespace gmx

#endif
