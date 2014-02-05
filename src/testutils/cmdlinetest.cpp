/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
 * Implements gmx::test::CommandLine.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "cmdlinetest.h"

#include <cstdlib>
#include <cstring>
#include <sstream>

#include <new>
#include <vector>

#include "gromacs/commandline/cmdlineprogramcontext.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{
namespace test
{

/********************************************************************
 * CommandLine::Impl
 */

class CommandLine::Impl
{
    public:
        Impl(const char *const cmdline[], size_t count);
        ~Impl();

        std::vector<char *>     args_;
        std::vector<char *>     argv_;
        int                     argc_;
};

CommandLine::Impl::Impl(const char *const cmdline[], size_t count)
{
    args_.reserve(count);
    argv_.reserve(count + 1);
    argc_ = static_cast<int>(count);
    for (size_t i = 0; i < count; ++i)
    {
        char *arg = strdup(cmdline[i]);
        if (arg == NULL)
        {
            throw std::bad_alloc();
        }
        args_.push_back(arg);
        argv_.push_back(arg);
    }
    argv_.push_back(NULL);
}

CommandLine::Impl::~Impl()
{
    for (size_t i = 0; i < args_.size(); ++i)
    {
        std::free(args_[i]);
    }
}

/********************************************************************
 * CommandLine
 */

CommandLine::CommandLine()
    : impl_(new Impl(NULL, 0))
{
}

CommandLine::CommandLine(const ConstArrayRef<const char *> &cmdline)
    : impl_(new Impl(cmdline.data(), cmdline.size()))
{
}

CommandLine::CommandLine(const CommandLine &other)
    : impl_(new Impl(other.argv(), other.argc()))
{
}

CommandLine::~CommandLine()
{
}

void CommandLine::initFromArray(const ConstArrayRef<const char *> &cmdline)
{
    impl_.reset(new Impl(cmdline.data(), cmdline.size()));
}

void CommandLine::append(const char *arg)
{
    GMX_RELEASE_ASSERT(impl_->argc_ == static_cast<int>(impl_->args_.size()),
                       "Command-line has been modified externally");
    size_t newSize = impl_->args_.size() + 1;
    impl_->args_.reserve(newSize);
    impl_->argv_.reserve(newSize + 1);
    char *newArg = strdup(arg);
    if (newArg == NULL)
    {
        throw std::bad_alloc();
    }
    impl_->args_.push_back(newArg);
    impl_->argv_.pop_back(); // Remove the trailing NULL.
    impl_->argv_.push_back(newArg);
    impl_->argv_.push_back(NULL);
    impl_->argc_ = static_cast<int>(newSize);
}

namespace
{

//! Helper function for converting values to strings
template <typename T>
std::string value2string(T value)
{
    std::stringstream ss;
    ss << value;
    return ss.str();
}

}

void CommandLine::addOption(const char *name,
                            const char *value)
{
    append(name);
    append(value);
}

void CommandLine::addOption(const char        *name,
                            const std::string &value)
{
    addOption(name, value.c_str());
}

void CommandLine::addOption(const char *name,
                            int         value)
{
    append(name);
    append(value2string(value));
}

void CommandLine::addOption(const char *name,
                            double      value)
{
    append(name);
    append(value2string(value));
}

int &CommandLine::argc()
{
    return impl_->argc_;
}
char **CommandLine::argv()
{
    return &impl_->argv_[0];
}
int CommandLine::argc() const
{
    return impl_->argc_;
}
const char *const *CommandLine::argv() const
{
    return &impl_->argv_[0];
}
const char *CommandLine::arg(int i) const
{
    return impl_->argv_[i];
}

std::string CommandLine::toString() const
{
    return CommandLineProgramContext(argc(), argv()).commandLine();
}

} // namespace test
} // namespace gmx
