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
 * Declares function to add the content of an xvg file to a checker.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_XVGTESTS_H
#define GMX_TESTUTILS_XVGTESTS_H

#include <string>

#include "testutils/testasserts.h"
#include "testutils/textblockmatchers.h"

namespace gmx
{

class TextInputStream;

namespace test
{

class TestReferenceChecker;

struct XvgMatchSettings
{
    XvgMatchSettings() : tolerance(defaultRealTolerance()), testData(true) {}

    FloatingPointTolerance tolerance;
    bool                   testData;
};

/*! \brief
 * Adds content of xvg file to TestReferenceChecker object.
 *
 * \param[in] input       Stream that provides the xvg content.
 * \param[in,out] checker Checker to use.
 * \param[in] settings    Settings to use for matching.
 *
 * Parses an xvg file from the input stream, and checks the contents against
 * reference data.  \p settings can be used to customize the matching.
 * Only a single data set is supported (but multiple columns work).
 * A subset of xmgrace formatting is also checked; static content that is
 * nearly always the same is skipped.
 *
 * \see XvgMatch
 */
void checkXvgFile(TextInputStream* input, TestReferenceChecker* checker, const XvgMatchSettings& settings);

/*! \libinternal \brief
 * Match the contents as an xvg file.
 *
 * \see checkXvgFile()
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class XvgMatch : public ITextBlockMatcherSettings
{
public:
    //! Sets the tolerance for matching data point values.
    XvgMatch& tolerance(const FloatingPointTolerance& tolerance)
    {
        settings_.tolerance = tolerance;
        return *this;
    }
    /*! \brief
     * Sets whether the actual data is checked.
     *
     * If set to `false`, only the legends are checked.  Use this if the
     * data is already tested using different means.
     */
    XvgMatch& testData(bool test)
    {
        settings_.testData = test;
        return *this;
    }

    TextBlockMatcherPointer createMatcher() const override;

private:
    XvgMatchSettings settings_;
};

} // namespace test

} // namespace gmx

#endif
