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
 * Implements functions in testoptions.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_testutils
 */
#include "testoptions.h"

#include <cstdio>
#include <cstdlib>

#include <new>
#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <gmock/gmock.h>

#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/errorcodes.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programinfo.h"

#include "datapath.h"
#include "refdata.h"
#include "testexceptions.h"

namespace
{

//! Stored command line for gmx::test::parseTestOptions().
boost::scoped_ptr<std::vector<std::string> > s_commandLine;

} // namespace

namespace gmx
{
namespace test
{

void initTestUtils(const char *dataPath, int *argc, char *argv[])
{
    try
    {
        ProgramInfo::init(*argc, argv);
        ::testing::InitGoogleMock(argc, argv);
        if (dataPath != NULL)
        {
            TestFileManager::setTestDataPath(dataPath);
        }
        initReferenceData(argc, argv);
        boost::scoped_ptr<std::vector<std::string> > commandLine(
                new std::vector<std::string>());
        for (int i = 0; i < *argc; ++i)
        {
            commandLine->push_back(argv[i]);
        }
        swap(commandLine, s_commandLine);
    }
    catch (const std::exception &ex)
    {
        printFatalErrorMessage(stderr, ex);
        std::exit(1);
    }
    ::gmx::setFatalErrorHandler(NULL);
}

void parseTestOptions(Options *options)
{
    GMX_RELEASE_ASSERT(s_commandLine.get() != NULL,
                       "Test options not initialized");
    try
    {
        CommandLineParser(options).parse(s_commandLine.get());
        options->finish();
    }
    catch (const GromacsException &ex)
    {
        GMX_THROW_WRAPPER_TESTEXCEPTION(ex);
    }
}

} // namespace test
} // namespace gmx
