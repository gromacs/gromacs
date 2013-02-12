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
 *
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
