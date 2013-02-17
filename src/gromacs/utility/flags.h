/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
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
 * \tparam T An enumerated type that holds the possible single flags.
 *
 * This class is not used publicly, but is present in an installed header
 * because it is used internally in public template classes.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
template <typename T>
class FlagsTemplate
{
    public:
        //! Creates a flags object with no flags set.
        FlagsTemplate() : flags_(0) {}
        //! Creates a flags object from a single flag.
        FlagsTemplate(T flag) : flags_(flag) {}

        //! Returns true if the given flag is set.
        bool test(T flag) const { return flags_ & flag; }
        //! Clears all flags.
        void clearAll() { flags_ = 0; }
        //! Sets the given flag.
        void set(T flag) { flags_ |= flag; }
        //! Clears the given flag.
        void clear(T flag) { flags_ &= ~flag; }
        //! Sets or clears the given flag.
        void set(T flag, bool bSet)
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
        FlagsTemplate<T> operator|(const FlagsTemplate<T> &other) const
        {
            return FlagsTemplate<T>(flags_ | other.flags_);
        }
        //! Combines flags from another flag object.
        FlagsTemplate<T> &operator|=(const FlagsTemplate<T> &other)
        {
            flags_ |= other.flags_;
            return *this;
        }
        //! Combined flags from two flags objects.
        FlagsTemplate<T> operator&(const FlagsTemplate<T> &other) const
        {
            return FlagsTemplate<T>(flags_ & other.flags_);
        }
        //! Returns an object with all flags flipped.
        FlagsTemplate<T> operator~() const
        {
            return FlagsTemplate<T>(~flags_);
        }

    private:
        //! Creates a flags object with the given flags.
        explicit FlagsTemplate(unsigned long flags) : flags_(flags) {}

        unsigned long           flags_;
};

} // namespace gmx

#endif
