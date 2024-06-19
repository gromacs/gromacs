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
/*!\file
 * \libinternal
 * \brief
 * Helpers and data for outputadapter module tests.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inlibraryapi
 * \ingroup module_coordinateio
 */

#ifndef GMX_COORDINATEIO_TESTS_MODULE_H
#define GMX_COORDINATEIO_TESTS_MODULE_H

#include "gmxpre.h"

#include "config.h"

#include <memory>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/coordinateio/coordinatefileenums.h"
#include "gromacs/coordinateio/ioutputadapter.h"
#include "gromacs/coordinateio/outputadaptercontainer.h"
#include "gromacs/coordinateio/outputadapters/outputselector.h"
#include "gromacs/coordinateio/outputadapters/setatoms.h"
#include "gromacs/coordinateio/outputadapters/setbox.h"
#include "gromacs/coordinateio/outputadapters/setforces.h"
#include "gromacs/coordinateio/outputadapters/setprecision.h"
#include "gromacs/coordinateio/outputadapters/setstarttime.h"
#include "gromacs/coordinateio/outputadapters/settimestep.h"
#include "gromacs/coordinateio/outputadapters/setvelocities.h"
#include "gromacs/coordinateio/requirements.h"
#include "gromacs/coordinateio/tests/coordinate_test.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/selection/selection.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

/*!\libinternal \brief  Helper to test supported file names. */
class SetAtomsSupportedFiles : public ModuleTest
{
public:
    void prepareTest(const char* filename)
    {
        addTopology();
        //! Storage for requirements.
        OutputRequirements requirements;

        requirements.atoms = ChangeAtomsType::AlwaysFromStructure;

        EXPECT_NO_THROW(runTest(filename, requirements));
    }
};

/*!\libinternal \brief  Helper to test supported file names. */
class SetAtomsUnSupportedFiles : public ModuleTest
{
public:
    void prepareTest(const char* filename)
    {
        //! Storage for requirements.
        OutputRequirements requirements;

        requirements.atoms = ChangeAtomsType::AlwaysFromStructure;

        EXPECT_THROW(runTest(filename, requirements), InconsistentInputError);
    }
};

/*!\libinternal \brief  Helper to test supported file names. */
class AnyOutputSupportedFiles : public ModuleTest, public ModuleSelection
{
public:
    void prepareTest(const char* filename)
    {
        addTopology();
        //! Storage for requirements.
        OutputRequirements requirements;
        //! Local atoms
        Selection sel;
        //! Local box
        matrix box;

        clear_mat(box);

        addOptionForSelection(&dummySelection_, true);
        setSelectionOptionValues(getOption(), &dummySelection_, true);

        copy_mat(requirements.newBox, box);
        requirements.box       = ChangeFrameInfoType::Always;
        requirements.frameTime = ChangeFrameTimeType::Both;

        EXPECT_NO_THROW(runTest(filename, requirements));
    }
};

/*!\libinternal \brief  Helper to test support for disabling all output options. */
class NoOptionalOutput : public ModuleTest
{
public:
    void prepareTest(const char* filename)
    {
        addTopology();
        //! Storage for requirements.
        OutputRequirements requirements;

        requirements.atoms     = ChangeAtomsType::Never;
        requirements.velocity  = ChangeSettingType::Never;
        requirements.force     = ChangeSettingType::Never;
        requirements.precision = ChangeFrameInfoType::Always;
        requirements.prec      = 3;

        EXPECT_NO_THROW(runTest(filename, requirements));
    }
};

/*!\libinternal \brief  Helper to test that invalid selection is rejected */
class OutputSelectorDeathTest : public ModuleTest, public ModuleSelection
{
public:
    void prepareTest()
    {
        //! Storage for frameadapters.
        OutputAdapterContainer adapters(CoordinateFileFlags::Base);
        //! Local atoms
        Selection sel;

        addOptionForSelection(&sel, false);
        setSelectionOptionValues(getOption(), &sel, false);

        GMX_ASSERT_DEATH_IF_SUPPORTED(adapters.addAdapter(std::make_unique<OutputSelector>(sel),
                                                          CoordinateFileFlags::RequireCoordinateSelection),
                                      "Need a valid selection out of simple atom indices");
    }
};

/*!\libinternal \brief  Helper to test supported file names. */
class SetVelocitySupportedFiles : public ModuleTest
{
public:
    void prepareTest(const char* filename)
    {
        addTopology();
        //! Storage for requirements.
        OutputRequirements requirements;

        requirements.velocity = ChangeSettingType::Always;

        EXPECT_NO_THROW(runTest(filename, requirements));
    }
};

/*!\libinternal \brief  Helper to test supported file names. */
class SetVelocityUnSupportedFiles : public ModuleTest
{
public:
    void prepareTest(const char* filename)
    {
        //! Storage for requirements.
        OutputRequirements requirements;

        requirements.velocity = ChangeSettingType::Always;

        EXPECT_THROW(runTest(filename, requirements), InconsistentInputError);
    }
};

/*!\libinternal \brief  Helper to test supported file names. */
class SetForceSupportedFiles : public ModuleTest
{
public:
    void prepareTest(const char* filename)
    {
        addTopology();
        //! Storage for requirements.
        OutputRequirements requirements;

        requirements.force = ChangeSettingType::Always;

        EXPECT_NO_THROW(runTest(filename, requirements));
    }
};

/*!\libinternal \brief  Helper to test supported file names. */
class SetForceUnSupportedFiles : public ModuleTest
{
public:
    void prepareTest(const char* filename)
    {
        //! Storage for requirements.
        OutputRequirements requirements;

        requirements.force = ChangeSettingType::Always;

        EXPECT_THROW(runTest(filename, requirements), InconsistentInputError);
    }
};

/*!\libinternal \brief  Helper to test supported file names. */
class SetPrecisionSupportedFiles : public ModuleTest
{
public:
    /*! \brief Set up the test case before running.
     *
     * \param[in] filename Name of the outputfile used to specify file type.
     */
    void prepareTest(const char* filename)
    {
        addTopology();
        OutputRequirements requirements;

        requirements.precision = ChangeFrameInfoType::Always;
        requirements.prec      = 5;

        EXPECT_NO_THROW(runTest(filename, requirements));
    }
};

/*!\libinternal \brief  Helper to test supported file names. */
class SetPrecisionUnSupportedFiles : public ModuleTest
{
public:
    /*! \brief Set up the test case before running.
     *
     * \param[in] filename Name of the outputfile used to specify file type.
     */
    void prepareTest(const char* filename)
    {
        OutputRequirements requirements;

        requirements.precision = ChangeFrameInfoType::Always;
        requirements.prec      = 5;

        EXPECT_THROW(runTest(filename, requirements), InconsistentInputError);
    }
};

//! Names here work for setAtoms module
const char* const setAtomsSupported[] = {
#if GMX_USE_TNG
    "spc2-traj.tng",
#endif
    "spc2-traj.gro",
    "spc2-traj.pdb",
};

//! Names here don't work for setAtoms module
const char* const setAtomsUnSupported[] = { "spc2-traj.trr", "spc2-traj.xtc", "spc2-traj.g96" };

/*! \brief
 *  Names here work for stuff that has no specific requirements.
 *
 *  PDB and GRO format are not tested here because they also require atoms information
 *  that is incompatible with the other output formats.
 */
const char* const anySupported[] = { "spc2-traj.trr",
#if GMX_USE_TNG
                                     "spc2-traj.tng",
#endif
                                     "spc2-traj.xtc",
                                     "spc2-traj.g96" };

//! Names here work for setVelocity module
const char* const setVelocitySupported[] = {
#if GMX_USE_TNG
    "spc2-traj.tng",
#endif
    "spc2-traj.gro",
    "spc2-traj.trr",
};

//! Names here don't work for setVelocity module
const char* const setVelocityUnSupported[] = { "spc2-traj.xtc", "spc2-traj.pdb", "spc2-traj.g96" };

//! Names here work for setForce module
const char* const setForceSupported[] = {
#if GMX_USE_TNG
    "spc2-traj.tng",
#endif
    "spc2-traj.trr",
};

//! Names here don't work for setForce module
const char* const setForceUnSupported[] = { "spc2-traj.xtc",
                                            "spc2-traj.pdb",
                                            "spc2-traj.gro",
                                            "spc2-traj.g96" };

//! Names here work for setPrecision module
const char* const setPrecisionSupported[] = {
#if GMX_USE_TNG
    "spc2-traj.tng",
#endif
    "spc2-traj.xtc",
};

//! Names here don't work for setPrecision module
const char* const setPrecisionUnSupported[] = { "spc2-traj.trr",
                                                "spc2-traj.pdb",
                                                "spc2-traj.gro",
                                                "spc2-traj.g96" };


} // namespace test

} // namespace gmx

#endif
