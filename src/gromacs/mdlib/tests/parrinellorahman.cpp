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
#include "gmxpre.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>

#include <array>
#include <initializer_list>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/matrix.h"
#include "gromacs/math/multidimarray.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/coupling.h"
#include "gromacs/mdspan/extents.h"
#include "gromacs/mdspan/layouts.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/boxutilities.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/loggerbuilder.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/naming.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

//! GoogleTest expectations about whether matrices are diagonal or not
void isAlmostDiagonalMatrix(const char* name, const Matrix3x3& m)
{
    SCOPED_TRACE(formatString("Testing that %s is almost diagonal", name));
    FloatingPointTolerance tolerance = absoluteTolerance(2e-5);
    EXPECT_REAL_EQ_TOL(0, m(XX, YY), tolerance);
    EXPECT_REAL_EQ_TOL(0, m(XX, ZZ), tolerance);
    EXPECT_REAL_EQ_TOL(0, m(YY, XX), tolerance);
    EXPECT_REAL_EQ_TOL(0, m(YY, ZZ), tolerance);
    EXPECT_REAL_EQ_TOL(0, m(ZZ, XX), tolerance);
    EXPECT_REAL_EQ_TOL(0, m(ZZ, YY), tolerance);
}

/*! \brief Some kinds of matrices of simulation box pressure
 *
 * The symmetric ones are useful for testing e.g. that anisotropic
 * coupling preserves symmetry if the inputs are symmetric. */
enum class PressureMatrixType
{
    UniformDiagonal,       //! all three diagonal entries are equal
    PartlyUniformDiagonal, //! two diagonal entries are equal
    Diagonal,              //! diagonal with unequal entries
    General,               //! All elements non-zero
    Extreme,               //! All elements non-zero and of large magnitude
    Count
};

/*! \brief Helper array for describing tests
 *
 * These are used in naming test reference data files, and we have
 * hard limits on file lengths, so we must abbreviate. */
constexpr gmx::EnumerationArray<PressureMatrixType, const char*> sc_pressureMatrixTypeNames = {
    "pmt_unif diag",
    "pmt_part unif diag",
    "pmt_diagonal",
    "pmt_general",
    "pmt_extreme"
};

//! Some matrices of simulation box pressure
const EnumerationArray<PressureMatrixType, Matrix3x3> c_pressureMatrices = {
    // clang-format off
     Matrix3x3{{ -2.1,    0,    0,
                    0, -2.1,    0,
                    0,    0, -2.1 }},
     Matrix3x3{{ -2.1,    0,    0,
                    0, -2.1,    0,
                    0,    0,  3.4 }},
     Matrix3x3{{ -2.1,    0,    0,
                    0, -1.7,    0,
                    0,    0,  3.4 }},
     Matrix3x3{{ -2.1,  1.3, -1.1,
                  1.4, -1.7,  2.3,
                 -1.0,  2.2,  3.4}},
     // One with very large pressure to ensure warnings are triggered
     Matrix3x3{{  2e5,  1e5,  3e5,
                  3e5,  1e5,  2e5,
                  1e5,  2e5,  3e5}}
    // clang-format on
};

//! Simulation box types
enum class BoxShape : int
{
    Cubic,
    Rectilinear,
    RhombicDodecahedronXYSquare,
    RhombicDodecahedronXYHexagon,
    TruncatedOctahedron,
    Other,
    Count
};

constexpr gmx::EnumerationArray<BoxShape, const char*> sc_boxShapeNames = {
    "shape_cubic",
    "shape_rect", // abbreviates "rectilinear"
    "shape_rds",  // abbreviates "rhombic dodecahedron XY square"
    "shape_rdh",  // abbreviates "rhombic dodecahedron XY hexagon"
    "shape_to",   // abbreviates "truncated octahedron",
    "shape_other"
};

/*! \brief Convenience typedef of the test input parameters
 *
 * Parameters:
 * - pressure-coupling MDP options
 * - the matrix properties of the input pressure
 * - the shape of the box
 * - the box matrix
 * - the box velocity matrix
 */
using ParrinelloRahmanTestParameters =
        std::tuple<PressureCouplingOptions, PressureMatrixType, BoxShape, Matrix3x3, Matrix3x3>;

//! Helper function for naming tests
std::string pressureCouplingOptionsToString(const PressureCouplingOptions options)
{
    return enumValueToString(options.epct);
}

//! Helper functor for naming tests
struct Matrix3x3ToString
{
    std::string prefix;
    std::string operator()(const Matrix3x3& m) const
    {
        return prefix + ' ' + doubleToString(m(XX, XX));
    }
};

//! Tuple of formatters to name the parameterized test cases
const NameOfTestFromTuple<ParrinelloRahmanTestParameters> sc_testNamer{ std::make_tuple(
        pressureCouplingOptionsToString,
        sc_pressureMatrixTypeNames,
        sc_boxShapeNames,
        Matrix3x3ToString{ "box" },
        Matrix3x3ToString{ "boxv" }) };

//! Test fixture - abbreviated ParrinelloRahman to ParrRahm for shorter refdata filenames
using ParrRahmTest = ::testing::TestWithParam<ParrinelloRahmanTestParameters>;

TEST_P(ParrRahmTest, Works)
{
    // In this test, we
    // * set up simulation conditions (a box of a given shape with size and velocity),
    // * perform the coupling one time to update velocity, M, and mu,
    // * do a consistency check
    // * do a regression check

    // Unpack the parameters
    auto [options, pressureMatrixType, boxShape, box, boxVelocity] = GetParam();

    // Fill legacy marix data structures
    matrix legacyPressure;
    fillLegacyMatrix(c_pressureMatrices[pressureMatrixType], legacyPressure);
    matrix legacyBox;
    fillLegacyMatrix(box, legacyBox);
    matrix legacyBoxVelocity;
    fillLegacyMatrix(boxVelocity, legacyBoxVelocity);

    // Do not test non-equilibrium box-deformation functionality yet.
    // This will be easier once parrinellorahman_pcoupl has been
    // modernized.
    matrix deform = { { 0._real } };

    // Make the relative box vectors that are used in preserving the
    // shape of the box during pressure coupling.
    matrix     legacyBoxRel;
    const int  ndim       = options.epct == PressureCouplingType::SemiIsotropic ? 2 : 3;
    const bool initBoxRel = true;
    clear_mat(legacyBoxRel);
    do_box_rel(ndim, deform, legacyBoxRel, legacyBox, initBoxRel);

    // Prepare to test the logging
    StringOutputStream logStream;
    LoggerBuilder      builder;
    builder.addTargetStream(MDLogger::LogLevel::Warning, &logStream);
    LoggerOwner     logOwner = builder.build();
    const MDLogger& mdlog    = logOwner.logger();

    // Call the Parrinello-Rahman pressure-coupling function to produce
    // new values in legacyBoxVelocity, M, and mu
    Matrix3x3 M;
    Matrix3x3 mu;
    parrinellorahman_pcoupl(
            mdlog, 0, options, deform, 0.02, legacyPressure, legacyBox, legacyBoxRel, legacyBoxVelocity, &M, &mu);

    // Test the outputs
    if (options.epct == PressureCouplingType::Isotropic)
    {
        isAlmostDiagonalMatrix("M", M);
        isAlmostDiagonalMatrix("mu", mu);
    }

    if (options.epct == PressureCouplingType::Isotropic)
    // The following should also work with non-isotropic couplings,
    // but does not because the numerical accuracy is not preserved
    // properly.
    {
        SCOPED_TRACE("Check the diagonal values of mu are equal");
        const uint64_t         singleUlpDiff = 10;
        const uint64_t         doubleUlpDiff = boxShape == BoxShape::Other ? 10 : 5;
        FloatingPointTolerance tolerance(0, 0, 0, 0, singleUlpDiff, doubleUlpDiff, false);
        EXPECT_REAL_EQ_TOL(mu(XX, XX), mu(YY, YY), tolerance);
        EXPECT_REAL_EQ_TOL(mu(XX, XX), mu(ZZ, ZZ), tolerance);
        SCOPED_TRACE("Check the diagonal values of M are equal");
        EXPECT_REAL_EQ_TOL(M(XX, XX), M(YY, YY), tolerance);
        EXPECT_REAL_EQ_TOL(M(XX, XX), M(ZZ, ZZ), tolerance);
    }

    // The "other" box is just a set of numbers we use for ensuring
    // that we aren't tying the implementation too closely to
    // conventional box shapes. It sometimes issues the warning, but
    // we don't care about testing whether and when it does.
    if (boxShape != BoxShape::Other)
    {
        // If the pressure is extreme, then the warning will be issued, so
        // we test that. Also, if the initial box velocity is large
        // enough relative to the box, then the scaling will large enough
        // to trigger the warning.
        if (pressureMatrixType == PressureMatrixType::Extreme
            || (std::abs(boxVelocity(0, 0)) >= 100._real && box(0, 0) < 10000._real))
        {
            EXPECT_NE(logStream.toString(), "");
        }
        else
        {
            EXPECT_EQ(logStream.toString(), "");
        }
    }

    // Finally run some regression tests on the actual numbers produced
    checkTestNameLength({});
    // We don't need to check the pressure values when using the very
    // extreme pressure that we used to test that the warning about
    // large box scaling can be issued, and they would need very
    // different tolerances also.
    if (pressureMatrixType != PressureMatrixType::Extreme)
    {
        TestReferenceData    refData;
        TestReferenceChecker checker(refData.rootChecker());
        const uint64_t       singleUlpDiff = 20, doubleUlpDiff = 10;
        checker.setDefaultTolerance(FloatingPointTolerance(0, 0, 0, 0, singleUlpDiff, doubleUlpDiff, false));
        std::vector<real> newBoxVelocities;
        for (int d = 0; d < DIM; d++)
        {
            for (int n = 0; n < DIM; n++)
            {
                newBoxVelocities.push_back(legacyBoxVelocity[d][n]);
            }
        }
        checker.checkSequence(newBoxVelocities.begin(), newBoxVelocities.end(), "Box Velocities");
    }
}

//! Template for box-making function
template<BoxShape boxShape>
Matrix3x3 makeBox(real spacing);

//! Template specialization to make a cubic box
template<>
Matrix3x3 makeBox<BoxShape::Cubic>(real spacing)
{
    // clang-format off
    return Matrix3x3
    {{ spacing, 0._real, 0._real,
       0._real, spacing, 0._real,
       0._real, 0._real, spacing }};
    // clang-format on
}

//! Template specialization to make a rectilinear box
template<>
Matrix3x3 makeBox<BoxShape::Rectilinear>(real spacing)
{
    // clang-format off
    return Matrix3x3
    {{ spacing, 0._real,           0._real,
       0._real, spacing / 2._real, 0._real,
       0._real, 0._real,           spacing * 3._real }};
    // clang-format on
}

//! Template specialization to make a rhombic dodecahedral box where the XY face is a square
template<>
Matrix3x3 makeBox<BoxShape::RhombicDodecahedronXYSquare>(real spacing)
{
    // clang-format off
    return Matrix3x3
    {{ spacing,           0._real,           0._real,
       0._real,           spacing,           0._real,
       spacing / 2._real, spacing / 2._real, spacing * std::sqrt(2._real) / 2._real }};
    // clang-format on
}

//! Template specialization to make a rhombic dodecahedral box where the XY face is a hexagon
template<>
Matrix3x3 makeBox<BoxShape::RhombicDodecahedronXYHexagon>(real spacing)
{
    // clang-format off
    return Matrix3x3
    {{ spacing,           0._real,                                0._real,
       spacing / 2._real, spacing * std::sqrt(3._real) / 2._real, 0._real,
       spacing / 2._real, spacing * std::sqrt(3._real) / 6._real, spacing * std::sqrt(6._real) / 3._real }};
    // clang-format on
}

//! Template specialization to make a truncated octahedral box
template<>
Matrix3x3 makeBox<BoxShape::TruncatedOctahedron>(real spacing)
{
    // clang-format off
    return Matrix3x3
    {{  spacing,            0._real,                                          0._real,
        spacing / 3._real,  spacing * std::sqrt(2._real) * 2._real / 3._real, 0._real,
       -spacing / 3._real,  spacing * std::sqrt(2._real)           / 3._real, spacing * std::sqrt(6._real) / 3._real }};
    // clang-format on
}

//! Template specialization to make an arbitrary conformant box
template<>
Matrix3x3 makeBox<BoxShape::Other>(real spacing)
{
    // clang-format off
    return Matrix3x3
    {{  spacing,           0._real,            0._real,
        spacing / 2._real, spacing * 1.1_real, 0._real,
       -spacing / 2._real, spacing * 1.2_real, spacing * 0.9_real }};
    // clang-format on
}

//! Make a vector of boxes whose type is \p boxShape and leading dimension one of the \p boxSizes
template<BoxShape boxShape>
std::vector<Matrix3x3> makeBoxes(std::initializer_list<real> boxSizes)
{
    std::vector<Matrix3x3> boxes;
    for (const real boxSize : boxSizes)
    {
        boxes.push_back(makeBox<boxShape>(boxSize));
    }
    return boxes;
}

//! Sets of box vectors to use in tests
const EnumerationArray<BoxShape, std::vector<Matrix3x3>> c_boxVectors = {
    makeBoxes<BoxShape::Cubic>({ 2.5, 1e4 }),
    makeBoxes<BoxShape::Rectilinear>({ 2.5, 1e4 }),
    makeBoxes<BoxShape::RhombicDodecahedronXYSquare>({ 2.5, 1e4 }),
    makeBoxes<BoxShape::RhombicDodecahedronXYHexagon>({ 2.5, 1e4 }),
    makeBoxes<BoxShape::TruncatedOctahedron>({ 2.5, 1e4 }),
    makeBoxes<BoxShape::Other>({ 2.5, 1e4 }),
};

//! Sets of box velocity vectors to use in tests
const EnumerationArray<BoxShape, std::vector<Matrix3x3>> c_boxVelocities = {
    makeBoxes<BoxShape::Cubic>({ -1e2, -1., -1e-7, 0., 1e-7, 1., 1e2 }),
    makeBoxes<BoxShape::Rectilinear>({ -1e2, -1., -1e-7, 0., 1e-7, 1., 1e2 }),
    makeBoxes<BoxShape::RhombicDodecahedronXYSquare>({ -1e2, -1., -1e-7, 0., 1e-7, 1., 1e2 }),
    makeBoxes<BoxShape::RhombicDodecahedronXYHexagon>({ -1e2, -1., -1e-7, 0., 1e-7, 1., 1e2 }),
    makeBoxes<BoxShape::TruncatedOctahedron>({ -1e2, -1., -1e-7, 0., 1e-7, 1., 1e2 }),
    makeBoxes<BoxShape::Other>({ -1e2, -1., -1e-7, 0., 1e-7, 1., 1e2 }),
};

//! Sets of pressure-coupling MDP options to use in tests
const std::vector<PressureCouplingOptions> c_options = []() {
    PressureCouplingOptions options;
    options.epc               = PressureCoupling::ParrinelloRahman;
    options.tau_p             = 1.;
    Matrix3x3 compressibility = diagonalMatrix<real, 3, 3>(4.5e-5);
    fillLegacyMatrix(compressibility, options.compress);

    Matrix3x3 referencePressure = identityMatrix<real, 3, 3>();
    fillLegacyMatrix(referencePressure, options.ref_p);

    std::vector<PressureCouplingOptions> optionsVector;
    for (const PressureCouplingType pressureCouplingType : { PressureCouplingType::Isotropic,
                                                             // PressureCouplingType::SemiIsotropic,
                                                             PressureCouplingType::Anisotropic })
    {
        options.epct = pressureCouplingType;
        optionsVector.push_back(options);
    }

    return optionsVector;
}();

using testing::Combine;
using testing::Values;
using testing::ValuesIn;
INSTANTIATE_TEST_SUITE_P(Cubic,
                         ParrRahmTest,
                         Combine(ValuesIn(c_options),
                                 ValuesIn(EnumerationWrapper<PressureMatrixType>{}),
                                 Values(BoxShape::Cubic),
                                 ValuesIn(c_boxVectors[BoxShape::Cubic]),
                                 ValuesIn(c_boxVelocities[BoxShape::Cubic])),
                         sc_testNamer);

INSTANTIATE_TEST_SUITE_P(Rectilinear,
                         ParrRahmTest,
                         Combine(ValuesIn(c_options),
                                 ValuesIn(EnumerationWrapper<PressureMatrixType>{}),
                                 Values(BoxShape::Rectilinear),
                                 ValuesIn(c_boxVectors[BoxShape::Rectilinear]),
                                 ValuesIn(c_boxVelocities[BoxShape::Rectilinear])),
                         sc_testNamer);

INSTANTIATE_TEST_SUITE_P(RhombDodecXYSquare,
                         ParrRahmTest,
                         Combine(ValuesIn(c_options),
                                 ValuesIn(EnumerationWrapper<PressureMatrixType>{}),
                                 Values(BoxShape::RhombicDodecahedronXYSquare),
                                 ValuesIn(c_boxVectors[BoxShape::RhombicDodecahedronXYSquare]),
                                 ValuesIn(c_boxVelocities[BoxShape::RhombicDodecahedronXYSquare])),
                         sc_testNamer);

INSTANTIATE_TEST_SUITE_P(RhombDodecXYHex,
                         ParrRahmTest,
                         Combine(ValuesIn(c_options),
                                 ValuesIn(EnumerationWrapper<PressureMatrixType>{}),
                                 Values(BoxShape::RhombicDodecahedronXYHexagon),
                                 ValuesIn(c_boxVectors[BoxShape::RhombicDodecahedronXYHexagon]),
                                 ValuesIn(c_boxVelocities[BoxShape::RhombicDodecahedronXYHexagon])),
                         sc_testNamer);

INSTANTIATE_TEST_SUITE_P(TruncOct,
                         ParrRahmTest,
                         Combine(ValuesIn(c_options),
                                 ValuesIn(EnumerationWrapper<PressureMatrixType>{}),
                                 Values(BoxShape::TruncatedOctahedron),
                                 ValuesIn(c_boxVectors[BoxShape::TruncatedOctahedron]),
                                 ValuesIn(c_boxVelocities[BoxShape::TruncatedOctahedron])),
                         sc_testNamer);

INSTANTIATE_TEST_SUITE_P(Other,
                         ParrRahmTest,
                         Combine(ValuesIn(c_options),
                                 ValuesIn(EnumerationWrapper<PressureMatrixType>{}),
                                 Values(BoxShape::Other),
                                 ValuesIn(c_boxVectors[BoxShape::Other]),
                                 ValuesIn(c_boxVelocities[BoxShape::Other])),
                         sc_testNamer);

} // namespace
} // namespace test
} // namespace gmx
