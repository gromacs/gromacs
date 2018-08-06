/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018, by the GROMACS development team, led by
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

#include "gromacs/fileio/filetypes.h"

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
    ConfMatchSettings() : tolerance_(defaultRealTolerance()),
                          filetype_(efGRO), fullCheck_(false)
    {
    }

    FloatingPointTolerance  tolerance_;
    int                     filetype_;
    bool                    fullCheck_;
};

/*! \brief
 * Adds content of a configuration file to TestReferenceChecker object.
 *
 * \param[in] input       Stream that provides the gro content.
 * \param[in,out] checker Checker to use.
 * \param[in] settings    Settings to use for matching.
 *
 * The exact form of the checker is determined by the file type used when
 * declaring the checker.
 */
void checkConfFile(TextInputStream         *input,
                   TestReferenceChecker    *checker,
                   const ConfMatchSettings &settings);

/*! \brief
 * Performs comparison checks on a GRO file.
 *
 * \param[in] input       Stream that provides the gro content.
 * \param[in,out] checker Checker to use.
 * \param[in] settings    Settings to use for matching.
 *
 * The checker parses the file from input and compares against the stored
 * reference data. If the fullCheck_ setting is used, the complete file
 * will be compared field by field. Otherwise only the first two lines are compared.
 */
void checkGroFile(TextInputStream         *input,
                  TestReferenceChecker    *checker,
                  const ConfMatchSettings &settings);

/*! \brief
 * Performs comparison checks on a PDB file.
 *
 * \param[in] input       Stream that provides the gro content.
 * \param[in,out] checker Checker to use.
 * \param[in] settings    Settings to use for matching.
 *
 * The checker parses the file from input and compares against the stored
 * reference data. If the fullCheck_ setting is used, the complete file
 * will be compared field by field. Otherwise only the title line is checked.
 */
void checkPdbFile(TextInputStream         *input,
                  TestReferenceChecker    *checker,
                  const ConfMatchSettings &settings);

/*! \libinternal \brief
 * Match the contents as an gro file.
 *
 * \see checkConfFile()
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class ConfMatch : public ITextBlockMatcherSettings
{
    public:
        //! Sets the tolerance for matching floating point values.
        ConfMatch &tolerance(const FloatingPointTolerance &tolerance)
        {
            settings_.tolerance_ = tolerance;
            return *this;
        }
        //! Sets the filetype of the output file, used to determine the checker that should be used.
        ConfMatch &filetype(const int &filetype)
        {
            settings_.filetype_ = filetype;
            return *this;
        }
        //! Sets the flag to determine if a complete file comparison should be attempted or not.
        ConfMatch &fullCheck(const bool &fullCheck)
        {
            settings_.fullCheck_ = fullCheck;
            return *this;
        }

        virtual TextBlockMatcherPointer createMatcher() const;

    private:
        ConfMatchSettings  settings_;
};

} // namespace test

} // namespace gmx

#endif
