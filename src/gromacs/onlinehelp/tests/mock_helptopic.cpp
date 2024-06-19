/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Implements classes in mock_helptopic.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_onlinehelp
 */
#include "gmxpre.h"

#include "mock_helptopic.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/onlinehelp/ihelptopic.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
namespace test
{

/********************************************************************
 * MockHelpTopic
 */

// static
MockHelpTopic& MockHelpTopic::addSubTopic(gmx::AbstractCompositeHelpTopic* parent,
                                          const char*                      name,
                                          const char*                      title,
                                          const char*                      text)
{
    MockHelpTopic* topic = new MockHelpTopic(name, title, text);
    parent->addSubTopic(gmx::HelpTopicPointer(topic));
    return *topic;
}

MockHelpTopic::MockHelpTopic(const char* name, const char* title, const char* text) :
    name_(name), title_(title), text_(text != nullptr ? text : "")
{
    if (!isNullOrEmpty(text))
    {
        using ::testing::_;
        using ::testing::Invoke;
        ON_CALL(*this, writeHelp(_)).WillByDefault(Invoke(this, &MockHelpTopic::writeHelpBase));
    }
}

MockHelpTopic::~MockHelpTopic() {}

const char* MockHelpTopic::name() const
{
    return name_;
}

const char* MockHelpTopic::title() const
{
    return title_;
}

std::string MockHelpTopic::helpText() const
{
    return text_;
}

MockHelpTopic& MockHelpTopic::addSubTopic(const char* name, const char* title, const char* text)
{
    return addSubTopic(this, name, title, text);
}

} // namespace test
} // namespace gmx
