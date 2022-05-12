/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
    typedef std::list<EntryPointer>             ChildList;
    typedef ChildList::const_iterator           ChildIterator;

    static EntryPointer createRoot() { return std::make_unique<ReferenceDataEntry>("", ""); }

    ReferenceDataEntry(const char* type, const char* id) :
        type_(type),
        id_(id != nullptr ? id : ""),
        isTextBlock_(false),
        hasBeenChecked_(false),
        correspondingOutputEntry_(nullptr)
    {
    }

    const std::string&  type() const { return type_; }
    const std::string&  id() const { return id_; }
    bool                isCompound() const { return !children_.empty(); }
    bool                isTextBlock() const { return isTextBlock_; }
    const std::string&  value() const { return value_; }
    const ChildList&    children() const { return children_; }
    ReferenceDataEntry* correspondingOutputEntry() const { return correspondingOutputEntry_; }

    bool idMatches(const char* id) const
    {
        return (id == nullptr && id_.empty()) || (id != nullptr && id_ == id);
    }

    ChildIterator findChild(const char* id, const ChildIterator& prev) const
    {
        if (children_.empty())
        {
            return children_.end();
        }
        ChildIterator child          = prev;
        bool          wrappingSearch = true;
        if (child != children_.end())
        {
            if (id == nullptr && (*child)->id().empty())
            {
                wrappingSearch = false;
                ++child;
                if (child == children_.end())
                {
                    return children_.end();
                }
            }
        }
        else
        {
            child          = children_.begin();
            wrappingSearch = false;
        }
        do
        {
            if ((*child)->idMatches(id))
            {
                return child;
            }
            ++child;
            if (wrappingSearch && child == children_.end())
            {
                child = children_.begin();
            }
        } while (child != children_.end() && child != prev);
        return children_.end();
    }
    bool isValidChild(const ChildIterator& prev) const { return prev != children_.end(); }

    bool hasBeenChecked() const { return hasBeenChecked_; }
    void setChecked() { hasBeenChecked_ = true; }

    void setCheckedIncludingChildren()
    {
        setChecked();
        for (const auto& child : children_)
        {
            child->setCheckedIncludingChildren();
        }
    }

    EntryPointer cloneToOutputEntry()
    {
        EntryPointer entry(new ReferenceDataEntry(type_.c_str(), id_.c_str()));
        setCorrespondingOutputEntry(entry.get());
        entry->setValue(value());
        entry->isTextBlock_ = isTextBlock_;
        return entry;
    }

    void setValue(const std::string& value)
    {
        GMX_RELEASE_ASSERT(!isCompound(), "Cannot have a value for a compound entry");
        value_ = value;
    }
    void setTextBlockValue(const std::string& value)
    {
        GMX_RELEASE_ASSERT(!isCompound(), "Cannot have a value for a compound entry");
        value_       = value;
        isTextBlock_ = true;
    }
    void makeCompound(const char* type)
    {
        type_.assign(type);
        value_.clear();
        isTextBlock_ = false;
    }
    void setCorrespondingOutputEntry(ReferenceDataEntry* entry)
    {
        GMX_RELEASE_ASSERT(correspondingOutputEntry_ == nullptr, "Output entry already exists");
        correspondingOutputEntry_ = entry;
    }
    ChildIterator addChild(EntryPointer child)
    {
        GMX_RELEASE_ASSERT(isCompound() || value_.empty(),
                           "Cannot make an entry with a value to a compound");
        return children_.insert(children_.end(), std::move(child));
    }

private:
    std::string         type_;
    std::string         id_;
    std::string         value_;
    bool                isTextBlock_;
    ChildList           children_;
    bool                hasBeenChecked_;
    ReferenceDataEntry* correspondingOutputEntry_;
};

} // namespace test
} // namespace gmx

#endif
