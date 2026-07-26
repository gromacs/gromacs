/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * \brief Bonded force test runners.
 *
 * Declares test runner interface for bonded force evaluation. The test runners
 * abstract class is used to unify the interfaces for CPU and GPU implementations,
 * so the same test can run on the same data using different backends.
 *
 * \ingroup module_listed_forces
 */
#ifndef GMX_LISTED_FORCES_TESTS_BONDEDTESTRUNNERS_H
#define GMX_LISTED_FORCES_TESTS_BONDEDTESTRUNNERS_H

#include <string>
#include <vector>

#include "gromacs/listed_forces/bonded.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/test_device.h"

#include "bondedtestdata.h"

struct t_pbc;

namespace gmx
{
namespace test
{

/*! \brief Bonded force test runner interface.
 *
 * Wraps the actual implementation (CPU or GPU) of bonded force evaluation
 * into a common interface.
 */
class IBondedTestRunner
{
public:
    virtual ~IBondedTestRunner() = default;

    /*! \brief Run the bonded kernel for the given flavor and fill output.
     *
     * \param[in]  input   Input parameters and interaction type
     * \param[in]  x       Coordinates
     * \param[in]  pbc     Periodic boundary conditions
     * \param[in]  lambda  Free-energy lambda
     * \param[in]  iatoms  Interaction atom indices
     * \param[in]  flavor  Bonded kernel flavor to run
     * \param[out] output  Output energies, forces, and shift forces
     */
    virtual void run(const iListInput&           input,
                     const PaddedVector<RVec>&   x,
                     const t_pbc&                pbc,
                     real                        lambda,
                     const std::vector<t_iatom>& iatoms,
                     BondedKernelFlavor          flavor,
                     OutputQuantities*           output) = 0;

    /*! \brief Whether this runner supports the given flavor for this input.
     *
     * \param[in] flavor  Bonded kernel flavor
     * \param[in] input   Input (used to check FEP, interaction type, etc.)
     * \param[in] lambda  Current lambda (for FEP)
     */
    virtual bool supportsFlavor(BondedKernelFlavor flavor, const iListInput& input, real lambda) const = 0;

    /*! \brief Human-readable description of the hardware used by this runner.
     */
    virtual std::string hardwareDescription() const = 0;
};

/*! \brief CPU implementation of the bonded force test runner.
 */
class BondedHostTestRunner : public IBondedTestRunner
{
public:
    void run(const iListInput&           input,
             const PaddedVector<RVec>&   x,
             const t_pbc&                pbc,
             real                        lambda,
             const std::vector<t_iatom>& iatoms,
             BondedKernelFlavor          flavor,
             OutputQuantities*           output) override;
    bool supportsFlavor(BondedKernelFlavor flavor, const iListInput& input, real lambda) const override;
    std::string hardwareDescription() const override { return "CPU"; }
};

/*! \brief GPU implementation of the bonded force test runner.
 */
class BondedDeviceTestRunner : public IBondedTestRunner
{
public:
    explicit BondedDeviceTestRunner(const TestDevice& testDevice) : testDevice_(testDevice) {}
    void run(const iListInput&           input,
             const PaddedVector<RVec>&   x,
             const t_pbc&                pbc,
             real                        lambda,
             const std::vector<t_iatom>& iatoms,
             BondedKernelFlavor          flavor,
             OutputQuantities*           output) override;
    bool supportsFlavor(BondedKernelFlavor flavor, const iListInput& input, real lambda) const override;
    std::string hardwareDescription() const override { return testDevice_.description(); }

private:
    const TestDevice& testDevice_;
};

} // namespace test
} // namespace gmx

#endif // GMX_LISTED_FORCES_TESTS_BONDEDTESTRUNNERS_H
