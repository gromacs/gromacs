/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
/*! \libinternal \file
 * \brief
 * Declares module creation function for applied electric fields
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_applied_forces
 * \inlibraryapi
 */
#ifndef GMX_APPLIED_FORCES_ELECTRICFIELD_H
#define GMX_APPLIED_FORCES_ELECTRICFIELD_H

#include <memory>
#include <string_view>

namespace gmx
{

class IMDModule;

/*! \internal
    \brief Information about the electric-field module.
 *
 * Provides name and method to create a electric-field module.
 */
struct ElectricFieldModuleInfo
{
    /*! \brief
     * Creates a module for an external electric field.
     *
     * The returned class describes the time dependent electric field that can
     * be applied to all charges in a simulation. The field is described
     * by the following:
     *     E(t) = A cos(omega*(t-t0))*exp(-sqr(t-t0)/(2.0*sqr(sigma)));
     * If sigma = 0 there is no pulse and we have instead
     *     E(t) = A cos(omega*t)
     *
     * force is kJ mol^-1 nm^-1 = e * kJ mol^-1 nm^-1 / e
     *
     * WARNING:
     * There can be problems with the virial.
     * Since the field is not self-consistent this is unavoidable.
     * For neutral molecules the virial is correct within this approximation.
     * For neutral systems with many charged molecules the error is small.
     * But for systems with a net charge or a few charged molecules
     * the error can be significant when the field is high.
     * Solution: implement a self-consistent electric field into PME.
     */
    static std::unique_ptr<IMDModule> create();
    //! The name of the module
    static constexpr std::string_view sc_name = "electric-field";
};

} // namespace gmx

#endif
