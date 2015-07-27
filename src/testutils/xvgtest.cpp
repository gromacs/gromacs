/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Implements routine to check the content of xvg files.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "xvgtest.h"

#include <cerrno>
#include <cstdlib>

#include <vector>

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textstream.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

void checkXvgFile(TextInputStream      *input,
                  TestReferenceChecker *checker)
{
    std::string line;
    int         nrow = 0;

    while (input->readLine(&line))
    {
        if (!((line.find("#") != line.npos) ||
              (line.find("@") != line.npos)))
        {
            std::vector<std::string> split = splitString(line);
            std::vector<real>        row;

            for (std::vector<std::string>::iterator si = split.begin();
                 (si < split.end()); ++si)
            {
                const char *ptr = si->c_str();
                char       *endptr;
                errno = 0;
                double      dval = std::strtod(ptr, &endptr);
                if (errno == ERANGE)
                {
                    GMX_THROW(InvalidInputError("Invalid value: '" + *si
                                                + "'; it causes an overflow/underflow"));
                }
                if (*ptr == '\0' || *endptr != '\0')
                {
                    GMX_THROW(InvalidInputError("Invalid value: '" + *si
                                                + "'; expected a number"));
                }
                row.push_back(dval);
            }
            std::string buf = formatString("Row%d", nrow++);
            checker->checkSequence(row.begin(), row.end(), buf.c_str());
        }
    }
    std::string buf = formatString("Row%d", nrow++);
    checker->checkPresent(false, buf.c_str());
}

} // namespace test

} // namespace gmx
