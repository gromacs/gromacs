/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Tests for gmx clustsize
 *
 * \todo These will be superseded by tests of the new style analysis
 * modules.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */

#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/modules/clustsize.h"

#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"
#include "testutils/textblockmatchers.h"
#include "testutils/xvgtest.h"

#include "moduletest.h"

namespace
{

using gmx::test::CommandLine;
using gmx::test::ExactTextMatch;
using gmx::test::NoTextMatch;
using gmx::test::XvgMatch;

/********************************************************************
 * Tests for gmx::analysismodules::ClustSize.
 */

//! Test fixture for the `clustsize` analysis module.
typedef gmx::test::TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::ClustSizeInfo>
    ClustsizeTest;

TEST_F(ClustsizeTest, NoMolDefaultCutoff)
{
    const char *const command[] = { "clustsize", "-selection", "0" };

    setInputFile("-f", "clustsize.pdb");
    setInputFile("-s", "clustsize.tpr");
    setInputFile("-n", "clustsize.ndx");
    setOutputFile("-mc", ".xvg", XvgMatch());
    setOutputFile("-nc", ".xvg", XvgMatch());
    setOutputFile("-ac", ".xvg", XvgMatch());
    setOutputFile("-hc", ".xvg", XvgMatch());
    setDatasetTolerance("avclust", gmx::test::relativeToleranceAsFloatingPoint(4, 1e-8));
    runTest(CommandLine(command));
}

TEST_F(ClustsizeTest, NoMolShortCutoff)
{
    const char *const command[] = { "clustsize", "-cut", "0.3", "-selection", "0" };

    setInputFile("-f", "clustsize.pdb");
    setInputFile("-s", "clustsize.tpr");
    setInputFile("-n", "clustsize.ndx");
    setOutputFile("-mc", ".xvg", XvgMatch());
    setOutputFile("-nc", ".xvg", XvgMatch());
    setOutputFile("-ac", ".xvg", XvgMatch());
    setOutputFile("-hc", ".xvg", XvgMatch());
    setDatasetTolerance("avclust", gmx::test::relativeToleranceAsFloatingPoint(4, 1e-8));
    runTest(CommandLine(command));
}

TEST_F(ClustsizeTest, MolDefaultCutoff)
{
    const char *const command[] = { "clustsize", "-mol", "-selection", "0" };

    setInputFile("-f", "clustsize.pdb");
    setInputFile("-s", "clustsize.tpr");
    setOutputFile("-mc", ".xvg", XvgMatch());
    setOutputFile("-nc", ".xvg", XvgMatch());
    setOutputFile("-ac", ".xvg", XvgMatch());
    setOutputFile("-hc", ".xvg", XvgMatch());
    setDatasetTolerance("avclust", gmx::test::relativeToleranceAsFloatingPoint(4, 1e-8));
    runTest(CommandLine(command));
}

TEST_F(ClustsizeTest, MolShortCutoff)
{
    const char *const command[] = { "clustsize", "-mol", "-cut", "0.3", "-selection", "0" };

    setInputFile("-f", "clustsize.pdb");
    setInputFile("-s", "clustsize.tpr");
    setOutputFile("-mc", ".xvg", XvgMatch());
    setOutputFile("-nc", ".xvg", XvgMatch());
    setOutputFile("-ac", ".xvg", XvgMatch());
    setOutputFile("-hc", ".xvg", XvgMatch());
    setDatasetTolerance("avclust", gmx::test::relativeToleranceAsFloatingPoint(4, 1e-8));
    runTest(CommandLine(command));
}

} // namespace
