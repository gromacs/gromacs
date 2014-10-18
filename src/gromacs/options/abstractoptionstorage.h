/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2014, by the GROMACS development team, led by
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
 * Declares gmx::AbstractOptionStorage.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_ABSTRACTOPTIONSTORAGE_H
#define GMX_OPTIONS_ABSTRACTOPTIONSTORAGE_H

#include <string>

#include "gromacs/options/optionflags.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class AbstractOption;
class OptionInfo;
class Options;

/*! \libinternal \brief
 * Abstract base class for converting, validating, and storing option values.
 *
 * This class should normally not be subclassed directly, but the
 * OptionStorageTemplate should be used instead.  The templated class provides
 * basic functionality for most of the pure virtual methods, and also
 * integrates well with option setting objects derived from OptionTemplate.
 *
 * \inlibraryapi
 * \ingroup module_options
 *
 * \internal
 * This class really consists of two parts: the public interface that is
 * used by the internal implementation of the options module, and the
 * interface that derived classes use to provide type-dependent functionality.
 * The latter consists of a few pure virtual methods, of which a few simple
 * query methods are also part of the module-internal interface, others are
 * protected and called by the non-virtual methods when needed.
 * The reason why these two roles are in one class is twofold:
 *  -# Both the derived classes and the internal module implementation may need
 *     access to the same information like the allowed number of values and the
 *     name of the option.
 *  -# Having only one class is consistent with the structure used for options
 *     settings objects: there is very direct correspondence between
 *     AbstractOption and AbstractOptionStorage and between OptionTemplate and
 *     OptionStorageTemplate.
 */
class AbstractOptionStorage
{
    public:
        virtual ~AbstractOptionStorage();

        //! Returns true if the option has been set.
        bool isSet() const { return hasFlag(efOption_Set); }
        /*! \brief
         * Returns true if the option is a boolean option.
         *
         * This is used to optionally support an alternative syntax where an
         * option provided with no value sets the value to true and an
         * option prefixed with "no" clears the value.
         */
        bool isBoolean() const;
        //! Returns true if the option is a hidden option.
        bool isHidden() const { return hasFlag(efOption_Hidden); }
        //! Returns true if the option is required.
        bool isRequired() const { return hasFlag(efOption_Required); }
        //! Returns true if the option is vector-valued.
        bool isVector() const { return hasFlag(efOption_Vector); }
        //! Returns the name of the option.
        const std::string &name() const { return name_; }
        //! Returns the description of the option set by the calling code.
        const std::string &description() const { return descr_; }

        //! Returns true if defaultValueIfSet() value is specified.
        bool defaultValueIfSetExists() const
        { return hasFlag(efOption_DefaultValueIfSetExists); }
        //! Returns the minimum number of values required in one set.
        int minValueCount() const { return minValueCount_; }
        //! Returns the maximum allowed number of values in one set (-1 = no limit).
        int maxValueCount() const { return maxValueCount_; }

        /*! \brief
         * Returns an option info object corresponding to this option.
         */
        virtual OptionInfo &optionInfo() = 0;
        /*! \brief
         * Returns a short string describing the type of the option.
         */
        virtual std::string typeString() const = 0;
        /*! \brief
         * Formats additional description for the option.
         *
         * If this method returns a non-empty string, it is appended to the
         * plain description when printing help texts.
         * The default implementation returns an empty string.
         */
        virtual std::string formatExtraDescription() const
        { return std::string(); }
        /*! \brief
         * Returns the number of option values added so far.
         */
        virtual int valueCount() const = 0;
        /*! \brief
         * Returns the i'th value formatted as a string.
         *
         * If \p i is DefaultValueIfSetIndex, should format the default value
         * if set (see OptionTemplate::defaultValueIfSet()).
         */
        virtual std::string formatValue(int i) const = 0;
        //! \copydoc OptionInfo::formatDefaultValueIfSet()
        std::string formatDefaultValueIfSet() const
        { return formatValue(DefaultValueIfSetIndex); }

        /*! \brief
         * Starts adding values from a new source for the option.
         *
         * This marks the vurrent value of the option as a default value,
         * causing next call to startSet() to clear it.  This allows values
         * from the new source to overwrite old values.
         *
         * This method does not throw.
         */
        void startSource();
        /*! \brief
         * Starts adding a new set of values for the option.
         *
         * \throws  InvalidInputError if option is specified multiple times,
         *      but is not specified to accept it.
         *
         * If the parameter is specified multiple times, startSet() should be
         * called before the values for each instance.
         *
         * Strong exception safety guarantee.
         */
        void startSet();
        /*! \brief
         * Adds a new value for the option, converting it from a string.
         *
         * \param[in] value  String value to convert.
         * \throws  InvalidInputError if value cannot be converted, or
         *      if there are too many values.
         *
         * This method should only be called between startSet() and
         * finishSet().
         */
        void appendValue(const std::string &value);
        /*! \brief
         * Performs validation and/or actions once a set of values has been
         * added.
         *
         * \throws  InvalidInputError if too few values have been provided, or
         *      if the valid values since previous startSet() are invalid as a
         *      set.
         *
         * If the parameter is specified multiple times, finishSet() should be
         * called after the values for each instance.
         */
        void finishSet();
        /*! \brief
         * Performs validation and/or actions once all values have been added.
         *
         * \throws InvalidInputError if the option is required but not set, or
         *      if all valid values together are invalid as a set.
         *
         * This method should be called after all values have been provided
         * with appendValue().
         */
        void finish();

    protected:
        //! Index used with formatValue() for formatting default value if set.
        static const int DefaultValueIfSetIndex = -1;

        /*! \brief
         * Initializes the storage object from the settings object.
         *
         * \param[in] settings  Option settings.
         * \param[in] staticFlags Option flags that are always set and specify
         *      generic behavior of the option.
         * \throws  APIError if invalid settings have been provided.
         */
        AbstractOptionStorage(const AbstractOption &settings,
                              OptionFlags           staticFlags);

        //! Marks the option as set.
        void markAsSet() { flags_.set(efOption_Set); }
        //! Returns true if the given flag is set.
        bool hasFlag(OptionFlag flag) const { return flags_.test(flag); }
        //! Sets the given flag.
        void setFlag(OptionFlag flag) { return flags_.set(flag); }
        //! Clears the given flag.
        void clearFlag(OptionFlag flag) { return flags_.clear(flag); }

        /*! \brief
         * Sets a new minimum number of values required in one set.
         *
         * \param[in] count  New minimum number of values (must be > 0).
         * \throws InvalidInputError if already provided values violate the limit.
         *
         * If values have already been provided, it is checked that there are
         * enough.
         *
         * Cannot be called for options with ::efOption_MultipleTimes set,
         * because it is impossible to check the requirement after the values
         * have been set.
         * If attempted, will assert.
         */
        void setMinValueCount(int count);
        /*! \brief
         * Sets a new maximum number of values required in one set.
         *
         * \param[in] count  New maximum number of values
         *                   (must be > 0, or -1 for no limit).
         * \throws InvalidInputError if already provided values violate the limit.
         *
         * If values have already been provided, it is checked that there are
         * not too many.
         *
         * Cannot be called for options with ::efOption_MultipleTimes set,
         * because it is impossible to check the requirement after the values
         * have been set.
         * If attempted, will assert.
         */
        void setMaxValueCount(int count);

        /*! \brief
         * Removes all values from temporary storage for a set.
         *
         * This function is always called before starting to add values to
         * a set, allowing the storage to clear its internal buffers.
         *
         * Should not throw.
         */
        virtual void clearSet() = 0;
        /*! \brief
         * Adds a new value, converting it from a string.
         *
         * \param[in] value  String value to convert.
         * \throws  InvalidInputError if \p value is not valid for this option
         *      or if there have been too many values in the set.
         *
         * This method may be called multiple times if the underlying
         * option is defined to accept multiple values.
         *
         * \see OptionStorageTemplate::convertValue()
         */
        virtual void convertValue(const std::string &value) = 0;
        /*! \brief
         * Performs validation and/or actions once a set of values has been
         * added.
         *
         * \throws  InvalidInputError if the values in the set are not valid
         *      as a whole.
         *
         * This method may be called multiple times if the underlying option
         * can be specified multiple times.
         * This method is not currently called if one of the convertValue()
         * calls throwed.
         *
         * \todo
         * Improve the call semantics.
         *
         * \see OptionStorageTemplate::processSetValues()
         */
        virtual void processSet() = 0;
        /*! \brief
         * Performs validation and/or actions once all values have been added.
         *
         * \throws  InvalidInputError if all provided values are not valid as
         *      a set.
         *
         * This method is always called once.
         *
         * If the method throws, implementation should take care to leave the
         * option in a consistent, meaningful state.  However, currently none
         * of the implementations actually throw in any situation where the
         * option may be left in an inconsistent state.
         */
        virtual void processAll() = 0;

    private:
        std::string             name_;
        std::string             descr_;
        //! Flags for the option.
        OptionFlags             flags_;
        //! Minimum number of values required (in one set).
        int                     minValueCount_;
        //! Maximum allowed number of values (in one set), or -1 if no limit.
        int                     maxValueCount_;
        //! Whether we are currently assigning values to a set.
        bool                    bInSet_;
        //! Whether there were errors in set values.
        bool                    bSetValuesHadErrors_;

        GMX_DISALLOW_COPY_AND_ASSIGN(AbstractOptionStorage);
};

} // namespace gmx

#endif
