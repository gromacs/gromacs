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

#include "testxvg.h"

#include <cstdlib>

#include <vector>

#include "gromacs/utility/filestream.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

void checkXvgFile(std::string                &fileName,
                  int                         nColumn,
                  double                      tolerance,
                  test::TestReferenceChecker &checker)
{
    if (fileName.size() == 0)
    {
        return;
    }
    TextInputFile                     fp(fileName);
    bool                              bOK;
    std::vector<std::vector<double> > result;

    if (nColumn < 1)
    {
        return;
    }
    result.resize(nColumn);

    checker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1, tolerance));
    do
    {

        std::string        line;

        bOK = fp.readLine(&line);

        if (bOK)
        {
            std::vector<std::string> split = splitString(line);
            int i = 0;

            for (std::vector<std::string>::iterator si = split.begin(); (si < split.end()) && (i < nColumn); ++si, ++i)
            {
                if ((si->rfind("#") != si->npos) ||
                    (si->rfind("@") != si->npos))
                {
                    break;
                }
                result[i].push_back(atof(si->c_str()));
            }
        }
    }
    while (bOK);

    for (int i = 0; (i < nColumn); i++)
    {
        if (result[i].size() > 0)
        {
            std::string buf = formatString("result%d", i);
            checker.checkSequenceArray(result[i].size(), &result[i][0], buf.c_str());
        }
    }
}

} // namespace test

} // namespace gmx
