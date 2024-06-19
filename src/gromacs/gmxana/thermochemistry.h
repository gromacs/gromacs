/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
 * Code for computing entropy and heat capacity from eigenvalues
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef GMXANA_THERMOCHEMISTRY_H
#define GMXANA_THERMOCHEMISTRY_H

#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

namespace gmx
{
template<typename>
class ArrayRef;
}

/*! \brief Compute zero point energy from an array of eigenvalues.
 *
 * This routine first converts the eigenvalues from a normal mode
 * analysis to frequencies and then computes the zero point energy.
 *
 * \param[in] eigval       The eigenvalues
 * \param[in] scale_factor Factor to scale frequencies by before computing cv
 * \return The zero point energy (kJ/mol)
 */
double calcZeroPointEnergy(gmx::ArrayRef<const real> eigval, real scale_factor);

/*! \brief Compute heat capacity due to vibrational motion
 *
 * \param[in] eigval       The eigenvalues
 * \param[in] temperature  Temperature (K)
 * \param[in] linear       TRUE if this is a linear molecule
 * \param[in] scale_factor Factor to scale frequencies by before computing cv
 * \return The heat capacity at constant volume (J/mol K)
 */
double calcVibrationalHeatCapacity(gmx::ArrayRef<const real> eigval, real temperature, bool linear, real scale_factor);

/*! \brief Compute entropy due to translational motion
 *
 * Following the equations in J. W. Ochterski,
 * Thermochemistry in Gaussian, Gaussian, Inc., 2000
 * Pitssburg PA
 *
 * \param[in] mass         Molecular mass (Dalton)
 * \param[in] temperature  Temperature (K)
 * \param[in] pressure     Pressure (bar) at which to compute
 * \returns The translational entropy (J/mol K)
 */
double calcTranslationalEntropy(real mass, real temperature, real pressure);

/*! \brief Compute entropy due to rotational motion
 *
 * Following the equations in J. W. Ochterski,
 * Thermochemistry in Gaussian, Gaussian, Inc., 2000
 * Pitssburg PA
 *
 * \param[in] temperature  Temperature (K)
 * \param[in] natom        Number of atoms
 * \param[in] linear       TRUE if this is a linear molecule
 * \param[in] theta        The principal moments of inertia (unit of Energy)
 * \param[in] sigma_r      Symmetry factor, should be >= 1
 * \returns The rotational entropy (J/mol K)
 */
double calcRotationalEntropy(real temperature, int natom, bool linear, const rvec theta, real sigma_r);

/*! \brief Compute internal energy due to vibrational motion
 *
 * \param[in] eigval       The eigenvalues
 * \param[in] temperature  Temperature (K)
 * \param[in] linear       TRUE if this is a linear molecule
 * \param[in] scale_factor Factor to scale frequencies by before computing E
 * \return The internal energy (J/mol K)
 */
double calcVibrationalInternalEnergy(gmx::ArrayRef<const real> eigval, real temperature, bool linear, real scale_factor);

/*! \brief Compute entropy using Schlitter formula
 *
 * Computes entropy for a molecule / molecular system using the
 * algorithm due to Schlitter (Chem. Phys. Lett. 215 (1993)
 * 617-621).
 * The input should be eigenvalues from a covariance analysis,
 * the units of the eigenvalues are those of energy.
 *
 * \param[in] eigval       The eigenvalues
 * \param[in] temperature  Temperature (K)
 * \param[in] linear       True if this is a linear molecule (typically a diatomic molecule).
 * \return the entropy (J/mol K)
 */
double calcSchlitterEntropy(gmx::ArrayRef<const real> eigval, real temperature, bool linear);

/*! \brief Compute entropy using Quasi-Harmonic formula
 *
 * Computes entropy for a molecule / molecular system using the
 * Quasi-harmonic algorithm (Macromolecules 1984, 17, 1370).
 * The input should be eigenvalues from a normal mode analysis.
 * In both cases the units of the eigenvalues are those of energy.
 *
 * \param[in] eigval       The eigenvalues
 * \param[in] temperature  Temperature (K)
 * \param[in] linear       True if this is a linear molecule (typically a diatomic molecule).
 * \param[in] scale_factor Factor to scale frequencies by before computing S0
 * \return the entropy (J/mol K)
 */
double calcQuasiHarmonicEntropy(gmx::ArrayRef<const real> eigval, real temperature, bool linear, real scale_factor);

#endif
