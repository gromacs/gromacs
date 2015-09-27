/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Declares internal data structures for the reference data framework.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_REFDATA_IMPL_H
#define GMX_TESTUTILS_REFDATA_IMPL_H

#include <list>
#include <memory>
#include <string>

#include "gromacs/utility/gmxassert.h"

namespace gmx
{
namespace test
{

class ReferenceDataEntry
{
    public:
        typedef std::unique_ptr<ReferenceDataEntry> EntryPointer;
        typedef std::list<EntryPointer> ChildList;
        typedef ChildList::const_iterator ChildIterator;

        static EntryPointer createRoot()
        {
            return EntryPointer(new ReferenceDataEntry("", ""));
        }

        ReferenceDataEntry(const char *type, const char *id)
            : type_(type), id_(id != NULL ? id : ""), isTextBlock_(false)
        {
        }

        const std::string &type() const { return type_; }
        const std::string &id() const { return id_; }
        bool isCompound() const { return !children_.empty(); }
        bool isTextBlock() const { return isTextBlock_; }
        const std::string &value() const { return value_; }
        const ChildList &children() const { return children_; }

        void setValue(const std::string &value)
        {
            GMX_RELEASE_ASSERT(!isCompound(),
                               "Cannot have a value for a compound entry");
            value_ = value;
        }
        void setTextBlockValue(const std::string &value)
        {
            GMX_RELEASE_ASSERT(!isCompound(),
                               "Cannot have a value for a compound entry");
            value_       = value;
            isTextBlock_ = true;
        }
        ChildIterator addChild(EntryPointer child)
        {
            GMX_RELEASE_ASSERT(isCompound() || value_.empty(),
                               "Cannot make an entry with a value to a compound");
            return children_.insert(children_.end(), move(child));
        }

    private:
        std::string  type_;
        std::string  id_;
        std::string  value_;
        bool         isTextBlock_;
        ChildList    children_;
};

} // namespace test
} // namespace gmx

#endif
