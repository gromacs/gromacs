/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * Declares input generator class for CP2K QMMM
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_QMMMINPUTGENERATOR_H
#define GMX_APPLIED_FORCES_QMMMINPUTGENERATOR_H

#include <set>
#include <string>

#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "qmmmtypes.h"

enum class PbcType : int;

namespace gmx
{

/*! \internal \brief
 * Class that takes QMMMParameters, Coordinates, Point charges, Box dimensions, pbcType.
 * Generates QM/MM sample input parameters and pdb-style coordinates for CP2K.
 * Input are generated as std::string objects which can be stored in tpr KVT
 * and/or flushed into the files.
 */
class QMMMInputGenerator
{
public:
    /*! \brief Construct QMMMInputGenerator from its parameters
     *
     * \param[in] parameters Structure with QMMM parameters
     * \param[in] pbcType Periodic boundary conditions
     * \param[in] box Matrix with full box of the system
     * \param[in] q Point charges of each atom in the system
     * \param[in] x Coordinates of each atom in the system
     */
    QMMMInputGenerator(const QMMMParameters& parameters,
                       PbcType               pbcType,
                       const matrix          box,
                       ArrayRef<const real>  q,
                       ArrayRef<const RVec>  x);

    /*! \brief Generates sample CP2K input file
     *
     */
    std::string generateCP2KInput() const;

    /*! \brief Generates PDB file suitable for usage with CP2K.
     *  In that PDB file Point Charges of MM atoms are provided with Extended Beta field
     */
    std::string generateCP2KPdb() const;

    //! \brief Returns computed QM box dimensions
    const matrix& qmBox() const;

    //! \brief Returns computed translation vector in order to center QM atoms inside QM box
    const RVec& qmTrans() const;

private:
    //! \brief Check if atom belongs to the global index of qmAtoms_
    bool isQMAtom(Index globalAtomIndex) const;

    /*!\brief Calculates dimensions and center of the QM box.
     *  Also evaluates translation for the system in order to center QM atoms inside QM box
     *
     *  \param[in] scale Factor of how much QM box would be bigger than the radius of QM system
     *  \param[in] minNorm Minimum norm of the QM box vector
     */
    void computeQMBox(real scale, real minNorm);

    //! \brief Generates &GLOBAL section of CP2K Input
    static std::string generateGlobalSection();

    //! \brief Generates &DFT section of CP2K Input
    std::string generateDFTSection() const;

    //! \brief Generates &QMMM section of CP2K Input
    std::string generateQMMMSection() const;

    //! \brief Generates &MM section of CP2K Input
    static std::string generateMMSection();

    //! \brief Generates &SUBSYS section of CP2K Input
    std::string generateSubsysSection() const;

    //! QMMM Parameters structure
    const QMMMParameters& parameters_;
    //! Simulation PbcType
    PbcType pbc_;
    //! Simulation Box
    matrix box_;
    //! QM box
    matrix qmBox_;
    //! PBC-aware center of QM subsystem
    RVec qmCenter_;
    //! Translation that shifts qmCenter_ to the center of qmBox_
    RVec qmTrans_;
    //! Set containing indexes of all QM atoms
    std::set<Index> qmAtoms_;
    //! Atoms point charges
    ArrayRef<const real> q_;
    //! Atoms coordinates
    ArrayRef<const RVec> x_;

    //! Scale of the generated QM box with respect to the QM-subsystem size
    static constexpr real sc_qmBoxScale = 1.5;
    //! Minimum length of the generated QM box vectors in nm
    static constexpr real sc_qmBoxMinLength = 1.0;
};

/*! \brief Transforms vector a such as distance from it to the plane defined by vectors b and c
 * will be h minimum length will be milL and maximum length maxL
 *
 * \param[in] a Vector which should be scaled
 * \param[in] b First vector that forms the plane
 * \param[in] c Second vector that forms the plane
 * \param[in] h Distance from the end of a to the plane of (b,c)
 * \param[in] minNorm Minimum norm of vector
 * \param[in] maxNorm Maximum norm of vector
 */
RVec computeQMBoxVec(const RVec& a, const RVec& b, const RVec& c, real h, real minNorm, real maxNorm);

} // namespace gmx

#endif // GMX_APPLIED_FORCES_QMMMINPUTGENERATOR_H
