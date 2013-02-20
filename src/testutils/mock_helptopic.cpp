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
 * Implements classes in mock_helptopic.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "mock_helptopic.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace gmx
{
namespace test
{

/********************************************************************
 * MockHelpTopic
 */

// static
MockHelpTopic &
MockHelpTopic::addSubTopic(gmx::AbstractCompositeHelpTopic *parent,
                           const char *name, const char *title,
                           const char *text)
{
    MockHelpTopic *topic = new MockHelpTopic(name, title, text);
    parent->addSubTopic(gmx::HelpTopicPointer(topic));
    return *topic;
}

MockHelpTopic::MockHelpTopic(const char *name, const char *title, const char *text)
    : name_(name), title_(title), text_(text)
{
}

MockHelpTopic::~MockHelpTopic()
{
}

const char *MockHelpTopic::name() const
{
    return name_;
}

const char *MockHelpTopic::title() const
{
    return title_;
}

std::string MockHelpTopic::helpText() const
{
    return text_;
}

MockHelpTopic &
MockHelpTopic::addSubTopic(const char *name, const char *title,
                           const char *text)
{
    return addSubTopic(this, name, title, text);
}

} // namespace test
} // namespace gmx
