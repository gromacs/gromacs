/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014,2015, by the GROMACS development team, led by
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
 * Declares mock implementation of gmx::HelpTopicInterface.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_onlinehelp
 */
#ifndef GMX_ONLINEHELP_TESTS_MOCK_HELPTOPIC_H
#define GMX_ONLINEHELP_TESTS_MOCK_HELPTOPIC_H

#include <gmock/gmock.h>

#include "gromacs/onlinehelp/helptopic.h"
#include "gromacs/onlinehelp/helpwritercontext.h"

namespace gmx
{
namespace test
{

class MockHelpTopic : public AbstractCompositeHelpTopic
{
    public:
        static MockHelpTopic &addSubTopic(
            gmx::AbstractCompositeHelpTopic *parent,
            const char *name, const char *title, const char *text);

        MockHelpTopic(const char *name, const char *title, const char *text);
        virtual ~MockHelpTopic();

        virtual const char *name() const;
        virtual const char *title() const;

        MOCK_CONST_METHOD1(writeHelp, void(const HelpWriterContext &context));

        MockHelpTopic &addSubTopic(const char *name, const char *title,
                                   const char *text);
        using AbstractCompositeHelpTopic::addSubTopic;

        /*! \brief
         * Calls base class writeHelp() method.
         *
         * This provides the possibility for the mock to do the actual help
         * writing.
         */
        void writeHelpBase(const HelpWriterContext &context)
        {
            AbstractCompositeHelpTopic::writeHelp(context);
        }

    private:
        virtual std::string helpText() const;

        const char             *name_;
        const char             *title_;
        const char             *text_;
};

} // namespace test
} // namespace gmx

#endif
