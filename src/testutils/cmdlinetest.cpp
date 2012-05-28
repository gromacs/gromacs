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
 * Implements gmx::test::CommandLine.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_testutils
 */
#include "cmdlinetest.h"

#include <cstdlib>
#include <cstring>

#include <new>
#include <vector>

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
    argv_.reserve(count);
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

CommandLine::CommandLine(const char *const cmdline[], size_t count)
    : impl_(new Impl(cmdline, count))
{
}

CommandLine::CommandLine(const CommandLine &other)
    : impl_(new Impl(other.argv(), other.argc()))
{
}

CommandLine::~CommandLine()
{
}

void CommandLine::append(const char *arg)
{
    size_t newSize = impl_->args_.size() + 1;
    impl_->args_.reserve(newSize);
    impl_->argv_.reserve(newSize);
    char *newArg = strdup(arg);
    if (newArg == NULL)
    {
        throw std::bad_alloc();
    }
    impl_->args_.push_back(newArg);
    impl_->argv_.push_back(newArg);
    impl_->argc_ = static_cast<int>(newSize);
}

int &CommandLine::argc() { return impl_->argc_; }
char **CommandLine::argv() { return &impl_->argv_[0]; }
int CommandLine::argc() const { return impl_->argc_; }
const char *const *CommandLine::argv() const { return &impl_->argv_[0]; }
const char *CommandLine::arg(int i) const { return impl_->argv_[i]; }

std::string CommandLine::toString() const
{
    std::string result;
    for (size_t i = 0; i < impl_->args_.size(); ++i)
    {
        if (i > 0)
        {
            result.append(" ");
        }
        const char *arg = impl_->args_[i];
        bool bSpaces = (std::strchr(arg, ' ') != NULL);
        if (bSpaces)
        {
            result.append("'");
        }
        result.append(arg);
        if (bSpaces)
        {
            result.append("'");
        }
    }
    return result;
}

} // namespace test
} // namespace gmx
