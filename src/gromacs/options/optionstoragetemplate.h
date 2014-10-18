/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014, by the GROMACS development team, led by
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
 * Defines gmx::OptionStorageTemplate template.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_OPTIONSTORAGETEMPLATE_H
#define GMX_OPTIONS_OPTIONSTORAGETEMPLATE_H

#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "gromacs/options/abstractoption.h"
#include "gromacs/options/abstractoptionstorage.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

class Options;

/*! \libinternal \brief
 * Templated base class for constructing option value storage classes.
 *
 * \tparam T Assignable type that stores a single option value.
 *
 * Provides an implementation of the clearSet(), valueCount(), and processSet()
 * methods of AbstractOptionStorage, as well as a basic no-action
 * implementation of processAll().  Two new virtual methods are added:
 * processSetValues() and refreshValues().  The default implementation of
 * processSetValues() does nothing, and refreshValues() is used to update
 * secondary storage after values have been added/changed.
 * This leaves typeString(), formatValue(), and convertValue() to be
 * implemented in derived classes.  processSetValues() and processAll() can
 * also be implemented if necessary.
 *
 * Implements transaction support for adding values within a set: all calls to
 * addValue() add the value to a temporary storage, processSetValues() operates
 * on this temporary storage, and commitValues() then copies these values to
 * the real storage.  commitValues() provides a strong exception safety
 * guarantee for the process (and it only throws if it runs out of memory).
 *
 * \inlibraryapi
 * \ingroup module_options
 */
template <typename T>
class OptionStorageTemplate : public AbstractOptionStorage
{
    public:
        //! Alias for the template class for use in base classes.
        typedef OptionStorageTemplate<T> MyBase;
        //! Type of the container that contains the current values.
        typedef std::vector<T> ValueList;

        virtual ~OptionStorageTemplate();

        // No implementation in this class for the pure virtual methods, but
        // the declarations are still included for clarity.
        virtual std::string typeString() const = 0;
        virtual int valueCount() const { return static_cast<int>(values_->size()); }
        /*! \copydoc gmx::AbstractOptionStorage::formatValue()
         *
         * OptionStorageTemplate implements handling of DefaultValueIfSetIndex
         * in this method, as well as checking that \p i is a valid index.
         * Derived classes must implement formatSingleValue() to provide the
         * actual formatting for a value of type \p T.
         */
        virtual std::string formatValue(int i) const;

    protected:
        /*! \brief
         * Initializes the storage from option settings.
         *
         * \param[in] settings  Option settings.
         * \param[in] staticFlags Option flags that are always set and specify
         *      generic behavior of the option.
         * \throws  APIError if invalid settings have been provided.
         */
        template <class U>
        explicit OptionStorageTemplate(const OptionTemplate<T, U> &settings,
                                       OptionFlags staticFlags = OptionFlags());


        virtual void clearSet();
        /*! \copydoc gmx::AbstractOptionStorage::convertValue()
         *
         * Derived classes should call addValue() after they have converted
         * \p value to the storage type.  It is allowed to call addValue()
         * more than once, or not at all.  OptionsAssigner::appendValue()
         * provides the same exception safety guarantee as this method, so it
         * should be considered whether the implementation can be made strongly
         * exception safe.
         */
        virtual void convertValue(const std::string &value) = 0;
        /*! \brief
         * Processes values for a set after all have been converted.
         *
         * \param[in,out] values Valid values in the set.
         * \throws InvalidInputError if the values do not form a valid set.
         *
         * This method is called after all convertValue() calls for a set.
         * \p values contains all values that were validly converted by
         * convertValue().  The derived class may alter the values, but should
         * in such a case ensure that a correct number of values is produced.
         * If the derived class throws, all values in \p values are discarded.
         */
        virtual void processSetValues(ValueList *values)
        {
            GMX_UNUSED_VALUE(values);
        }
        /*! \copydoc gmx::AbstractOptionStorage::processSet()
         *
         * OptionStorageTemplate implements transaction support for a set of
         * values in this method (see the class description), and provides a
         * more detailed processSetValues() method that can be overridden in
         * subclasses to process the actual values.  Derived classes should
         * override that method instead of this one if set value processing is
         * necessary.
         */
        virtual void processSet();
        /*! \copydoc gmx::AbstractOptionStorage::processAll()
         *
         * The implementation in OptionStorageTemplate does nothing.
         */
        virtual void processAll()
        {
        }
        /*! \brief
         * Formats a single value as a string.
         *
         * \param[in] value  Value to format.
         * \returns   \p value formatted as a string.
         *
         * The derived class must provide this method to format values a
         * strings.  Called by formatValue() to do the actual formatting.
         */
        virtual std::string formatSingleValue(const T &value) const = 0;

        /*! \brief
         * Removes all values from the storage.
         *
         * Does not throw.
         */
        void clear() { values_->clear(); }
        /*! \brief
         * Adds a value to a temporary storage.
         *
         * \param[in] value  Value to add. A copy is made.
         * \throws std::bad_alloc if out of memory.
         * \throws InvalidInputError if the maximum value count has been reached.
         *
         * Derived classes should call this function from the convertValue()
         * implementation to add converted values to the storage.
         * If the maximum value count has been reached, the value is discarded
         * and an exception is thrown.
         *
         * If adding values outside convertValue() (e.g., to set a custom
         * default value), derived classes should call clearSet() before adding
         * values (unless in the constructor) and commitValues() once all
         * values are added.
         */
        void addValue(const T &value);
        /*! \brief
         * Commits values added with addValue().
         *
         * \throws std::bad_alloc if out of memory.
         *
         * If this function succeeds, values added with addValue() since the
         * previous clearSet() are added to the storage for the option.
         * Only throws in out-of-memory conditions, and provides the strong
         * exception safety guarantee.
         *
         * See addValue() for cases where this method should be used in derived
         * classes.
         *
         * Calls refreshValues() and clearSet() if it is successful.
         */
        void commitValues();
        /*! \brief
         * Updates alternative store locations.
         *
         * Derived classes should override this method if they implement
         * alternative store locations, and copy/translate values from the
         * values() vector to these alternative storages.  They should also
         * call the base class implementation as part of their implementation.
         *
         * Should be called in derived classes if values are modified directly
         * through the values() method, e.g., in processAll().  Does not need
         * to be called if commitValues() is used.
         *
         * Does not throw, and derived classes should not change that.
         */
        virtual void refreshValues();

        /*! \brief
         * Sets the default value for the option.
         *
         * \param[in] value  Default value to set.
         * \throws    std::bad_alloc if out of memory.
         *
         * This method can be used from the derived class constructor to
         * programmatically set a default value.
         */
        void setDefaultValue(const T &value);
        /*! \brief
         * Sets the default value if set for the option.
         *
         * \param[in] value  Default value to set.
         * \throws    std::bad_alloc if out of memory.
         *
         * This method can be used from the derived class constructor to
         * programmatically set a default value.
         */
        void setDefaultValueIfSet(const T &value);

        /*! \brief
         * Provides derived classes access to the current list of values.
         *
         * The non-const variant should only be used from processAll() in
         * derived classes if necessary, and refreshValues() should be called
         * if any changes are made.
         */
        ValueList       &values() { return *values_; }
        //! Provides derived classes access to the current list of values.
        const ValueList &values() const { return *values_; }

    private:
        /*! \brief
         * Vector for temporary storage of values before commitSet() is called.
         */
        ValueList               setValues_;
        /*! \brief
         * Vector for primary storage of option values.
         *
         * Is never NULL; points either to externally provided vector, or an
         * internally allocated one.  The allocation is performed by the
         * constructor.
         *
         * Primarily, modifications to values are done only to this storage,
         * and other storage locations are updated only when refreshValues() is
         * called.
         */
        ValueList                   *values_;
        T                           *store_;
        int                         *countptr_;
        boost::scoped_ptr<ValueList> ownedValues_;
        boost::scoped_ptr<T>         defaultValueIfSet_;

        // Copy and assign disallowed by base.
};


template <typename T>
template <class U>
OptionStorageTemplate<T>::OptionStorageTemplate(const OptionTemplate<T, U> &settings,
                                                OptionFlags staticFlags)
    : AbstractOptionStorage(settings, staticFlags),
      values_(settings.storeVector_),
      store_(settings.store_),
      countptr_(settings.countptr_)
{
    // If the maximum number of values is not known, storage to
    // caller-allocated memory is unsafe.
    if (store_ != NULL && (maxValueCount() < 0 || hasFlag(efOption_MultipleTimes)))
    {
        GMX_THROW(APIError("Cannot set user-allocated storage for arbitrary number of values"));
    }
    if (values_ == NULL)
    {
        ownedValues_.reset(new std::vector<T>);
        values_ = ownedValues_.get();
    }
    if (hasFlag(efOption_NoDefaultValue)
        && (settings.defaultValue_ != NULL
            || settings.defaultValueIfSet_ != NULL))
    {
        GMX_THROW(APIError("Option does not support default value, but one is set"));
    }
    if (store_ != NULL && countptr_ == NULL && !isVector()
        && minValueCount() != maxValueCount())
    {
        GMX_THROW(APIError("Count storage is not set, although the number of produced values is not known"));
    }
    if (!hasFlag(efOption_NoDefaultValue))
    {
        setFlag(efOption_HasDefaultValue);
        if (settings.defaultValue_ != NULL)
        {
            setDefaultValue(*settings.defaultValue_);
        }
        else if (ownedValues_.get() != NULL && store_ != NULL)
        {
            values_->clear();
            int count = (settings.isVector() ?
                         settings.maxValueCount_ : settings.minValueCount_);
            for (int i = 0; i < count; ++i)
            {
                values_->push_back(store_[i]);
            }
        }
        if (settings.defaultValueIfSet_ != NULL)
        {
            setDefaultValueIfSet(*settings.defaultValueIfSet_);
        }
    }
}


template <typename T>
OptionStorageTemplate<T>::~OptionStorageTemplate()
{
}


template <typename T>
std::string OptionStorageTemplate<T>::formatValue(int i) const
{
    GMX_RELEASE_ASSERT(i == DefaultValueIfSetIndex || (i >= 0 && i < valueCount()),
                       "Invalid value index");
    if (i == DefaultValueIfSetIndex)
    {
        if (defaultValueIfSet_.get() != NULL)
        {
            return formatSingleValue(*defaultValueIfSet_);
        }
        return std::string();
    }
    return formatSingleValue(values()[i]);
}


template <typename T>
void OptionStorageTemplate<T>::clearSet()
{
    setValues_.clear();
}


template <typename T>
void OptionStorageTemplate<T>::processSet()
{
    processSetValues(&setValues_);
    if (setValues_.empty() && defaultValueIfSet_.get() != NULL)
    {
        addValue(*defaultValueIfSet_);
        setFlag(efOption_HasDefaultValue);
    }
    else
    {
        clearFlag(efOption_HasDefaultValue);
    }
    if (!hasFlag(efOption_DontCheckMinimumCount)
        && setValues_.size() < static_cast<size_t>(minValueCount()))
    {
        GMX_THROW(InvalidInputError("Too few (valid) values"));
    }
    commitValues();
}


template <typename T>
void OptionStorageTemplate<T>::addValue(const T &value)
{
    if (maxValueCount() >= 0
        && setValues_.size() >= static_cast<size_t>(maxValueCount()))
    {
        GMX_THROW(InvalidInputError("Too many values"));
    }
    setValues_.push_back(value);
}


template <typename T>
void OptionStorageTemplate<T>::commitValues()
{
    if (hasFlag(efOption_ClearOnNextSet))
    {
        values_->swap(setValues_);
    }
    else
    {
        values_->insert(values_->end(), setValues_.begin(), setValues_.end());
    }
    clearSet();
    refreshValues();
}


template <typename T>
void OptionStorageTemplate<T>::refreshValues()
{
    if (countptr_ != NULL)
    {
        *countptr_ = static_cast<int>(values_->size());
    }
    if (store_ != NULL)
    {
        for (size_t i = 0; i < values_->size(); ++i)
        {
            store_[i] = (*values_)[i];
        }
    }
}


template <typename T>
void OptionStorageTemplate<T>::setDefaultValue(const T &value)
{
    if (hasFlag(efOption_NoDefaultValue))
    {
        GMX_THROW(APIError("Option does not support default value, but one is set"));
    }
    if (hasFlag(efOption_HasDefaultValue))
    {
        setFlag(efOption_ExplicitDefaultValue);
        clear();
        clearSet();
        addValue(value);
        // TODO: As this is called from the constructor, it should not call
        // virtual functions.
        commitValues();
    }
}


template <typename T>
void OptionStorageTemplate<T>::setDefaultValueIfSet(const T &value)
{
    if (hasFlag(efOption_NoDefaultValue))
    {
        GMX_THROW(APIError("Option does not support default value, but one is set"));
    }
    if (hasFlag(efOption_MultipleTimes))
    {
        GMX_THROW(APIError("defaultValueIfSet() is not supported with allowMultiple()"));
    }
    setFlag(efOption_DefaultValueIfSetExists);
    defaultValueIfSet_.reset(new T(value));
}

} // namespace gmx

#endif
