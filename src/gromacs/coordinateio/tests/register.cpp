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
/*!\file
 * \internal
 * \brief
 * Tests for frameconverter registration method.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */

#include "gmxpre.h"

#include "gromacs/coordinateio/frameconverters/register.h"

#include <memory>
#include <numeric>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/coordinateio/frameconverterenums.h"
#include "gromacs/coordinateio/iframeconverter.h"
#include "gromacs/coordinateio/tests/frameconverter.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/trajectory/trajectoryframe.h"

namespace gmx
{

namespace test
{

class RegisterFrameConverterTest : public FrameConverterTestBase
{
public:
    //! Run the test.
    void runTest();
    /*! \brief Checks that the memory stored is consistent.
     *
     * \param[in] expectedOutcome If results should match or not.
     */
    void checkMemory(bool expectedOutcome);
    //! Access underlying storage object.
    ProcessFrameConversion* method() { return &method_; }
    //! Access to none owning frame object.
    t_trxframe* newFrame() { return newFrame_; }

private:
    //! Storage object.
    ProcessFrameConversion method_;
    //! Non owning pointer to new coordinate frame.
    t_trxframe* newFrame_;
};

void RegisterFrameConverterTest::runTest()
{
    newFrame_ = method()->prepareAndTransformCoordinates(frame());
}
void RegisterFrameConverterTest::checkMemory(bool expectedOutcome)
{
    EXPECT_EQ(expectedOutcome, x() == newFrame()->x);
    EXPECT_EQ(expectedOutcome, v() == newFrame()->v);
    EXPECT_EQ(expectedOutcome, f() == newFrame()->f);
}

TEST_F(RegisterFrameConverterTest, NoConverterWorks)
{
    EXPECT_EQ(0, method()->getNumberOfConverters());
    runTest();
    EXPECT_TRUE((method()->guarantee() & convertFlag(FrameConverterFlags::NoGuarantee)) != 0U);
    EXPECT_FALSE((method()->guarantee() & convertFlag(FrameConverterFlags::AtomsInBox)) != 0U);
    checkMemory(false);
}

TEST_F(RegisterFrameConverterTest, RegistrationWorks)
{
    EXPECT_EQ(0, method()->getNumberOfConverters());
    method()->addFrameConverter(std::make_unique<DummyConverter>(FrameConverterFlags::AtomsInBox));
    EXPECT_EQ(1, method()->getNumberOfConverters());
    method()->addFrameConverter(
            std::make_unique<DummyConverter>(FrameConverterFlags::FitToReferenceProgressive));
    EXPECT_EQ(2, method()->getNumberOfConverters());

    runTest();
    EXPECT_TRUE((method()->guarantee() & convertFlag(FrameConverterFlags::AtomsInBox)) != 0U);
    checkMemory(false);
}

TEST_F(RegisterFrameConverterTest, NewConverterCanInvalidateGuarantees)
{
    method()->addFrameConverter(std::make_unique<DummyConverter>(FrameConverterFlags::AtomsInBox));
    method()->addFrameConverter(std::make_unique<DummyConverter>(FrameConverterFlags::NewSystemCenter));
    method()->addFrameConverter(std::make_unique<DummyConverter>(FrameConverterFlags::MoleculeCOMInBox));
    runTest();
    EXPECT_FALSE((method()->guarantee() & convertFlag(FrameConverterFlags::AtomsInBox)) != 0U);
    EXPECT_TRUE((method()->guarantee() & convertFlag(FrameConverterFlags::MoleculeCOMInBox)) != 0U);
}

} // namespace test

} // namespace gmx
