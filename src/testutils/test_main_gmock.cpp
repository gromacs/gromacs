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
/*! \libinternal \file
 * \brief
 * main() for unit tests that use Google C++ Mocking Framework.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 */
#include <gmock/gmock.h>

#include "gromacs/fatalerror/fatalerror.h"
#include "testutils/datapath.h"
#include "testutils/refdata.h"

/*! \brief
 * Initializes unit testing with Google C++ Mocking Framework.
 */
int main(int argc, char *argv[])
{
    ::testing::InitGoogleMock(&argc, argv);
#ifdef TEST_DATA_PATH
    ::gmx::test::setTestDataPath(TEST_DATA_PATH);
    if (::gmx::test::initReferenceData(&argc, argv) != 0)
    {
        return 1;
    }
#endif
    ::gmx::setFatalErrorHandler(NULL);
    return RUN_ALL_TESTS();
}
