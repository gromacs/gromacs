/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief Shared test data and utilities for bonded force tests.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_listed_forces
 */
#ifndef GMX_LISTED_FORCES_TESTS_BONDEDTESTDATA_H
#define GMX_LISTED_FORCES_TESTS_BONDEDTESTDATA_H

#include <array>
#include <optional>
#include <ostream>
#include <vector>

#include "gromacs/listed_forces/bonded.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"

struct t_pbc;

namespace gmx
{
namespace test
{

//! Number of atoms used in these tests.
constexpr int c_numAtomsBondedTest = 4;

//! Charges used for bonded test
static constexpr std::array<real, c_numAtomsBondedTest> c_bondedTestCharges = { { 1.5, -2.0, 1.5, -1.0 } };

/*! \brief Output from bonded kernels
 *
 * \todo Later this might turn into the actual output struct. */
struct OutputQuantities
{
    //! Energy of this interaction
    real energy = 0;
    //! Derivative with respect to lambda
    real dvdlambda = 0;
    //! Shift vectors
    rvec fshift[c_numShiftVectors] = { { 0 } };
    //! Forces
    alignas(GMX_REAL_MAX_SIMD_WIDTH * sizeof(real)) rvec4 f[c_numAtomsBondedTest] = { { 0 } };
};

/*! \brief Input structure for listed forces tests
 */
struct iListInput
{
public:
    //! Function type
    std::optional<InteractionFunction> ftype;
    //! Tolerance for float evaluation
    float ftoler = 1e-6;
    //! Tolerance for double evaluation
    double dtoler = 1e-8;
    //! Do free energy perturbation?
    bool fep = false;
    //! Interaction parameters
    t_iparams iparams = { { 0 } };

    friend std::ostream& operator<<(std::ostream& out, const iListInput& input);

    //! Constructor
    iListInput() {}

    /*! \brief Constructor with tolerance
     *
     * \param[in] ftol Single precision tolerance
     * \param[in] dtol Double precision tolerance
     */
    iListInput(float ftol, double dtol)
    {
        ftoler = ftol;
        dtoler = dtol;
    }
    /*! \brief Set parameters for harmonic potential
     */
    iListInput setHarmonic(InteractionFunction ft, real rA, real krA, real rB, real krB);
    iListInput setHarmonic(InteractionFunction ft, real rA, real krA);
    iListInput setCubic(real b0, real kb, real kcub);
    iListInput setMorse(real b0A, real cbA, real betaA, real b0B, real cbB, real betaB);
    iListInput setMorse(real b0A, real cbA, real betaA);
    iListInput setFene(real bm, real kb);
    iListInput setLinearAngle(real klinA, real aA, real klinB, real aB);
    iListInput setLinearAngle(real klinA, real aA);
    iListInput
    setUreyBradley(real thetaA, real kthetaA, real r13A, real kUBA, real thetaB, real kthetaB, real r13B, real kUBB);
    iListInput setUreyBradley(real thetaA, real kthetaA, real r13A, real kUBA);
    iListInput setCrossBondBonds(real r1e, real r2e, real krr);
    iListInput setCrossBondAngles(real r1e, real r2e, real r3e, real krt);
    iListInput setQuarticAngles(real theta, const real c[5]);
    iListInput setPDihedrals(InteractionFunction ft, real phiA, real cpA, int mult, real phiB, real cpB);
    iListInput setPDihedrals(InteractionFunction ft, real phiA, real cpA, int mult);
    iListInput setRbDihedrals(const real rbcA[NR_RBDIHS], const real rbcB[NR_RBDIHS]);
    iListInput setRbDihedrals(const real rbc[NR_RBDIHS]);
    iListInput setPolarization(real alpha);
    iListInput setAnharmPolarization(real alpha, real drcut, real khyp);
    iListInput setTholePolarization(real a, real alpha1, real alpha2);
    iListInput setWaterPolarization(real alpha_x, real alpha_y, real alpha_z, real rOH, real rHH, real rOD);
};

/*! \brief Utility to fill iatoms struct
 *
 * \param[in]  ftype  Function type
 * \param[out] iatoms Pointer to iatoms struct
 */
void fillIatoms(std::optional<InteractionFunction> ftype, std::vector<t_iatom>* iatoms);

class TestReferenceChecker;

/*! \brief Utility to check the output from bonded tests
 *
 * \param[in] checker Reference checker
 * \param[in] output  The output from the test to check
 * \param[in] bondedKernelFlavor  Flavor for determining what output to check
 */
void checkOutput(TestReferenceChecker*   checker,
                 const OutputQuantities& output,
                 BondedKernelFlavor      bondedKernelFlavor);

} // namespace test
} // namespace gmx

#endif // GMX_LISTED_FORCES_TESTS_BONDEDTESTDATA_H
