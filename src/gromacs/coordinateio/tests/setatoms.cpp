/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
/*!\internal
 * \file
 * \brief
 * Tests for gmx::SetAtoms
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */


#include "gmxpre.h"

#include "gromacs/coordinateio/outputadapters/setatoms.h"

#include <filesystem>
#include <memory>
#include <string>
#include <utility>

#include <gtest/gtest.h>

#include "gromacs/coordinateio/coordinatefileenums.h"
#include "gromacs/coordinateio/tests/coordinate_test.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace gmx
{

namespace test
{

/*!\brief
 * Test fixture to prepare a setatoms object to pass data through.
 */
class SetAtomsTest : public gmx::test::CommandLineTestBase
{
public:
    //! Prepare test fixture topology to have atoms available.
    SetAtomsTest()
    {
        topology_.fillFromInputFile(TestFileManager::getInputFilePath("lysozyme.pdb").string());
    }

    /*! \brief
     * Get access to the method for changing atom information.
     *
     * \param[in] type Type to use for setting up method.
     * \param[in] haveAtomsAvailable Decide if method shouöd have atoms stored or not.
     */
    SetAtoms* setAtoms(ChangeAtomsType type, bool haveAtomsAvailable)
    {
        if (!setAtoms_)
        {
            if (haveAtomsAvailable)
            {
                setAtoms_ = std::make_unique<SetAtoms>(type, topology_.copyAtoms());
            }
            else
            {
                setAtoms_ = std::make_unique<SetAtoms>(type, nullptr);
            }
        }
        return setAtoms_.get();
    }
    //! Get access to a t_atoms struct to use in the dummy t_trxframe.
    t_atoms* getAtomsCopy() { return const_cast<t_atoms*>(topology_.atoms()); }

private:
    //! Object to use for tests
    SetAtomsPointer setAtoms_;
    //! Local topology to get atoms from
    TopologyInformation topology_;
};

TEST_F(SetAtomsTest, RemovesExistingAtoms)
{
    t_trxframe frame;
    clear_trxframe(&frame, true);

    frame.bAtoms = true;
    frame.atoms  = getAtomsCopy();

    EXPECT_NE(frame.atoms, nullptr);

    SetAtoms* method = setAtoms(ChangeAtomsType::Never, true);
    EXPECT_NO_THROW(method->processFrame(0, &frame));

    EXPECT_FALSE(frame.bAtoms);
    EXPECT_EQ(frame.atoms, nullptr);
}

TEST_F(SetAtomsTest, AddsNewAtoms)
{
    t_trxframe frame;
    clear_trxframe(&frame, true);

    frame.bAtoms = false;
    frame.atoms  = nullptr;

    SetAtoms* method = setAtoms(ChangeAtomsType::AlwaysFromStructure, true);
    EXPECT_NO_THROW(method->processFrame(0, &frame));

    EXPECT_TRUE(frame.bAtoms);
    EXPECT_NE(frame.atoms, nullptr);
}

TEST_F(SetAtomsTest, ThrowsOnRequiredAtomsNotAvailable)
{
    t_trxframe frame;
    clear_trxframe(&frame, true);

    frame.bAtoms = false;
    frame.atoms  = nullptr;

    SetAtoms* method = setAtoms(ChangeAtomsType::Always, false);
    EXPECT_THROW(method->processFrame(0, &frame), InconsistentInputError);
}

TEST_F(SetAtomsTest, WillUseOldAtomsWhenNoNewAvailable)
{
    t_trxframe frame;
    clear_trxframe(&frame, true);

    frame.bAtoms = true;
    frame.atoms  = getAtomsCopy();

    EXPECT_NE(frame.atoms, nullptr);

    SetAtoms* method = setAtoms(ChangeAtomsType::Always, false);
    EXPECT_NO_THROW(method->processFrame(0, &frame));
}

TEST_F(SetAtomsTest, ThrowsWhenUserAtomReplacementNotPossible)
{
    t_trxframe frame;
    clear_trxframe(&frame, true);

    frame.bAtoms = true;
    frame.atoms  = getAtomsCopy();

    EXPECT_NE(frame.atoms, nullptr);

    SetAtoms* method = setAtoms(ChangeAtomsType::AlwaysFromStructure, false);
    EXPECT_THROW(method->processFrame(0, &frame), InconsistentInputError);
}

} // namespace test

} // namespace gmx
