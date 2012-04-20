/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Implements functions in datapath.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_testutils
 */
#include "datapath.h"

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/path.h"

static const char *g_testDataPath = NULL;

namespace gmx
{
namespace test
{

std::string getTestFilePath(const char *filename)
{
    return Path::join(getTestDataPath(), filename);
}

const char *getTestDataPath()
{
    GMX_RELEASE_ASSERT(g_testDataPath != NULL, "Test data path not set");
    return g_testDataPath;
}

void setTestDataPath(const char *path)
{
    GMX_RELEASE_ASSERT(Directory::exists(path),
                       "Test data directory does not exist");
    g_testDataPath = path;
}

} // namespace test
} // namespace gmx
