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
/*! \libinternal \file
 * \brief
 * Declares function to add the content of a conf file to a checker.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_CONFTEST_H
#define GMX_TESTUTILS_CONFTEST_H

#include <string>

#include "testutils/testasserts.h"
#include "testutils/textblockmatchers.h"

namespace gmx
{

class TextInputStream;

namespace test
{

class TestReferenceChecker;

struct ConfMatchSettings
{
    FloatingPointTolerance tolerance              = defaultRealTolerance();
    bool                   matchFullConfiguration = false;
};

/*! \brief
 * Adds content of a gro file to TestReferenceChecker object.
 *
 * \param[in] input       Stream that provides the gro content.
 * \param[in,out] checker Checker to use.
 * \param[in] settings    Settings to use for matching.
 *
 * Parses a gro file from the input stream, and checks the contents against
 * reference data (only first two lines for now).
 *
 * \see ConfMatch
 */
void checkConfFile(TextInputStream* input, TestReferenceChecker* checker, const ConfMatchSettings& settings);

/*! \libinternal \brief
 * Match the contents as an gro file.
 *
 * \see checkGroFile()
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class ConfMatch : public ITextBlockMatcherSettings
{
public:
    //! Sets the tolerance for matching floating point values.
    ConfMatch& tolerance(const FloatingPointTolerance& tolerance)
    {
        settings_.tolerance = tolerance;
        return *this;
    }
    //! Sets whether the full file contents should be matched
    ConfMatch& matchFullConfiguration(const bool flag)
    {
        settings_.matchFullConfiguration = flag;
        return *this;
    }

    TextBlockMatcherPointer createMatcher() const override;

private:
    ConfMatchSettings settings_;
};

} // namespace test

} // namespace gmx

#endif
