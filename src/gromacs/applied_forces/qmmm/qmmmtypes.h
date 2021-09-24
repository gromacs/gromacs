/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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
 * Declares structers and types needed to evaluate forces and energies for QM/MM
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_QMMMTYPES_H
#define GMX_APPLIED_FORCES_QMMMTYPES_H

#include <memory>
#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

namespace gmx
{

/*! \internal
 * \brief Helper structure with indexes of broken bonds between QM and MM
 * Used to determine and store pair of QM and MM atoms between which chemical bond is broken
 */
struct LinkFrontier
{
    //! Global index of QM atom at Frontier
    index qm;
    //! Global index of MM atom at Frontier
    index mm;
};

/*! \brief Enumerator for supported QM methods
 * Also could be INPUT which means external input file provided
 * with the name determined by QMMMParameters::qminputfilename_
 */
enum class QMMMQMMethod
{
    PBE,   //!< DFT with PBE functional
    BLYP,  //!< DFT with BLYP functional
    INPUT, //!< User provides suitable input file for QM package
    Count
};

//! The names of the supported QM methods
static const EnumerationArray<QMMMQMMethod, const char*> c_qmmmQMMethodNames = {
    { "PBE", "BLYP", "INPUT" }
};

//! symbols of the elements in periodic table
const std::vector<std::string> periodic_system = {
    "X  ", "H  ", "He ", "Li ", "Be ", "B  ", "C  ", "N  ", "O  ", "F  ", "Ne ", "Na ",
    "Mg ", "Al ", "Si ", "P  ", "S  ", "Cl ", "Ar ", "K  ", "Ca ", "Sc ", "Ti ", "V  ",
    "Cr ", "Mn ", "Fe ", "Co ", "Ni ", "Cu ", "Zn ", "Ga ", "Ge ", "As ", "Se ", "Br ",
    "Kr ", "Rb ", "Sr ", "Y  ", "Zr ", "Nb ", "Mo ", "Tc ", "Ru ", "Rh ", "Pd ", "Ag ",
    "Cd ", "In ", "Sn ", "Sb ", "Te ", "I  ", "Xe ", "Cs ", "Ba ", "La ", "Ce ", "Pr ",
    "Nd ", "Pm ", "Sm ", "Eu ", "Gd ", "Tb ", "Dy ", "Ho ", "Er ", "Tm ", "Yb ", "Lu ",
    "Hf ", "Ta ", "W  ", "Re ", "Os ", "Ir ", "Pt ", "Au ", "Hg ", "Tl ", "Pb ", "Bi ",
    "Po ", "At ", "Rn ", "Fr ", "Ra ", "Ac ", "Th ", "Pa ", "U  ", "Np ", "Pu ", "Am ",
    "Cm ", "Bk ", "Cf ", "Es ", "Fm ", "Md ", "No ", "Lr ", "Rf ", "Db ", "Sg ", "Bh ",
    "Hs ", "Mt ", "Ds ", "Rg ", "Cn ", "Nh ", "Fl ", "Mc ", "Lv ", "Ts ", "Og "
};

/*! \internal
 * \brief Holding all parameters needed for QM/MM simulation.
 * Also used for setting all default parameter values.
 */
struct QMMMParameters
{
    //! Indicate if QM/MM is active (default false)
    bool active_ = false;
    //! Indices of the atoms that are part of the QM region (default whole System)
    std::vector<index> qmIndices_;
    //! Indices of the atoms that are part of the MM region (default no MM atoms)
    std::vector<index> mmIndices_;
    //! Vector with pairs of indicies defining broken bonds in QMMM (default determined from topology)
    std::vector<LinkFrontier> link_;
    //! Vector with atomic numbers of all atoms in the system (default determined from topology)
    std::vector<int> atomNumbers_;
    //! Total charge of QM system (default 0)
    int qmCharge_ = 0;
    //! Total multiplicity of QM system (default 1)
    int qmMultiplicity_ = 1;
    //! Method used for QM calculation (default DFT with PBE functional)
    QMMMQMMethod qmMethod_ = QMMMQMMethod::PBE;
    /*! \brief String containing name of the CP2K files (*.inp, *.out, *.pdb)
     * default value empty, means will be deduced from *.tpr name during mdrun
     */
    std::string qmFileNameBase_;
    //! String containing whole CP2K input which can be stored inside *.tpr
    std::string qmInput_;
    //! String containing PDB file for CP2K input which can be stored inside *.tpr
    std::string qmPdb_;
    //! Matrix that contains vectors defining QM box
    matrix qmBox_;
    //! Translation vector to center QM subsystem inside the QM Box
    RVec qmTrans_;

    //! Constructor with default initializers for arrays
    QMMMParameters() :
        qmBox_{ { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } }, qmTrans_{ 0.0, 0.0, 0.0 }
    {
    }

    GMX_DISALLOW_COPY_AND_ASSIGN(QMMMParameters);
};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_QMMMTYPES_H
