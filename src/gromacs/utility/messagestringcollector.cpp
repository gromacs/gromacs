/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2014, by the GROMACS development team, led by
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
 * Implements gmx::MessageStringCollector.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "messagestringcollector.h"

#include <vector>

#include "gromacs/utility/gmxassert.h"

namespace gmx
{

class MessageStringCollector::Impl
{
    public:
        Impl() : prevContext_(0) {}

        std::vector<std::string> contexts_;
        std::string              text_;
        size_t                   prevContext_;
};

MessageStringCollector::MessageStringCollector()
    : impl_(new Impl)
{
}

MessageStringCollector::~MessageStringCollector()
{
}

void MessageStringCollector::startContext(const char *name)
{
    impl_->contexts_.push_back(name);
}

void MessageStringCollector::append(const std::string &message)
{
    int indent = static_cast<int>(impl_->prevContext_ * 2);
    if (!impl_->contexts_.empty())
    {
        std::vector<std::string>::const_iterator ci;
        for (ci = impl_->contexts_.begin() + impl_->prevContext_;
             ci != impl_->contexts_.end(); ++ci)
        {
            impl_->text_.append(indent, ' ');
            impl_->text_.append(*ci);
            impl_->text_.append("\n");
            indent += 2;
        }
    }
    impl_->prevContext_ = impl_->contexts_.size();

    // TODO: Put this into a more generic helper, could be useful elsewhere
    size_t pos = 0;
    while (pos < message.size())
    {
        size_t nextpos = message.find_first_of('\n', pos);
        impl_->text_.append(indent, ' ');
        impl_->text_.append(message.substr(pos, nextpos - pos));
        impl_->text_.append("\n");
        if (nextpos == std::string::npos)
        {
            break;
        }
        pos = nextpos + 1;
    }
}

void MessageStringCollector::finishContext()
{
    GMX_RELEASE_ASSERT(!impl_->contexts_.empty(),
                       "finishContext() called without context");
    impl_->contexts_.pop_back();
    if (impl_->prevContext_ > impl_->contexts_.size())
    {
        impl_->prevContext_ = impl_->contexts_.size();
    }
}

void MessageStringCollector::clear()
{
    impl_->contexts_.clear();
    impl_->text_.clear();
    impl_->prevContext_ = 0;
}

bool MessageStringCollector::isEmpty() const
{
    return impl_->text_.empty();
}

std::string MessageStringCollector::toString() const
{
    return impl_->text_;
}

} // namespace gmx
