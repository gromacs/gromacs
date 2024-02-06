/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * Declares helper functions for reading structure factors from datafile
 *
 * \author Alexey Shvetsov <alexxyum@gmail.com>
 * \ingroup module_trajectoryanalysis
 */

#ifndef GMX_TRAJECTORYANALYSIS_MODULES_SCATTERINGFACTORS_H
#define GMX_TRAJECTORYANALYSIS_MODULES_SCATTERINGFACTORS_H

#include <array>
#include <string>
#include <vector>

namespace gmx
{

/*! \internal \brief
 * Cromer-Mann scattering factor parameters to compute structure factor dependent on Q
 *
 * \f[s(q) = \sum_{i}^{4}  a_i * \exp(- b_i * (\frac{q}{4\pi})^2) + c\f]
 * s(q) units is number of electrons so s(0) ~ Z (atomic number)
 */
struct CromerMannParameters
{
    //! parameter a
    std::array<double, 4> a;
    //! parameter b
    std::array<double, 4> b;
    //! parameter c
    double c;
};

/*! \internal \brief
 * Neutron scattering factor parameters for an atom type.
 */
struct AtomicStructureFactor
{
    //! atomic isotope
    std::string isotope;
    //! Number of Protons
    int atomicNumber;
    //! neutron scattering length (1e-15m)
    double neutronIsotropicScatteringLength;
    //! Parameters for the Cromer Mann fit
    CromerMannParameters xrayCromerMannParameters;
};

//! Helper function to read in atomic scattering data from file.
std::vector<AtomicStructureFactor> readAtomicStructureFactors();


} // namespace gmx


#endif // GMX_TRAJECTORYANALYSIS_MODULES_SCATTERINGFACTORS_H
