/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares gmx::FlagsTemplate.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_FLAGS_H
#define GMX_UTILITY_FLAGS_H

namespace gmx
{

/*! \brief
 * Template class for typesafe handling of combination of flags.
 *
 * \tparam FlagType An enumerated type that holds the possible single flags.
 *
 * This class is not used publicly, but is present in an installed header
 * because it is used internally in public template classes.
 *
 * Currently, it is not completely transparent, since or'ing together two
 * \p FlagType flags does not automatically create a FlagsTemplate object.
 * Also, some operators and more complex operations (like testing for multiple
 * flags at the same time) are missing, but can be added if the need arises.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
template <typename FlagType>
class FlagsTemplate
{
    public:
        //! Creates a flags object with no flags set.
        FlagsTemplate() : flags_(0) {}
        //! Creates a flags object from a single flag.
        FlagsTemplate(FlagType flag) : flags_(flag) {}

        /*! \brief
         * Tests if the given flag is set.
         *
         * Note that if \p flag has more than a single bit set, then returns
         * true if any of them is set.
         */
        bool test(FlagType flag) const { return (flags_ & flag) != 0; }
        //! Clears all flags.
        void clearAll() { flags_ = 0; }
        //! Sets the given flag.
        void set(FlagType flag) { flags_ |= flag; }
        //! Clears the given flag.
        void clear(FlagType flag) { flags_ &= ~flag; }
        //! Sets or clears the given flag.
        void set(FlagType flag, bool bSet)
        {
            if (bSet)
            {
                set(flag);
            }
            else
            {
                clear(flag);
            }
        }

        //! Combines flags from two flags objects.
        FlagsTemplate<FlagType>
        operator|(const FlagsTemplate<FlagType> &other) const
        {
            return FlagsTemplate<FlagType>(flags_ | other.flags_);
        }
        //! Combines flags from another flag object.
        FlagsTemplate<FlagType> &
        operator|=(const FlagsTemplate<FlagType> &other)
        {
            flags_ |= other.flags_;
            return *this;
        }
        //! Combined flags from two flags objects.
        FlagsTemplate<FlagType>
        operator&(const FlagsTemplate<FlagType> &other) const
        {
            return FlagsTemplate<FlagType>(flags_ & other.flags_);
        }
        //! Returns an object with all flags flipped.
        FlagsTemplate<FlagType> operator~() const
        {
            return FlagsTemplate<FlagType>(~flags_);
        }

    private:
        //! Creates a flags object with the given flags.
        explicit FlagsTemplate(unsigned long flags) : flags_(flags) {}

        unsigned long           flags_;
};

} // namespace gmx

#endif
