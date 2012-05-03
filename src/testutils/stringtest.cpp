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
 * Implements gmx::test::StringTestBase.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_testutils
 */
#include "testutils/stringtest.h"

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/file.h"
#include "testutils/testoptions.h"

namespace gmx
{
namespace test
{

bool StringTestBase::s_bWriteToStdOut = false;

void StringTestBase::SetUpTestCase()
{
    Options options(NULL, NULL);
    options.addOption(BooleanOption("stdout").store(&s_bWriteToStdOut));
    parseTestOptions(&options);
}

StringTestBase::StringTestBase()
{
}

StringTestBase::~StringTestBase()
{
}

TestReferenceChecker &
StringTestBase::checker()
{
    if (checker_.get() == NULL)
    {
        checker_.reset(new TestReferenceChecker(data_.rootChecker()));
    }
    return *checker_;
}

void
StringTestBase::checkText(const std::string &text, const char *id)
{
    if (s_bWriteToStdOut)
    {
        printf("%s:\n", id);
        printf("%s[END]\n", text.c_str());
    }
    else
    {
        checker().checkStringBlock(text, id);
    }
}

void
StringTestBase::checkFileContents(const std::string &filename, const char *id)
{
    std::string text = File::readToString(filename);
    checkText(text, id);
}

} // namespace test
} // namespace gmx
