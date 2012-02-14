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
/*! \libinternal \file
 * \brief
 * Defines gmx::OptionStorageTemplate template.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_OPTIONSTORAGETEMPLATE_H
#define GMX_OPTIONS_OPTIONSTORAGETEMPLATE_H

#include <memory>
#include <string>
#include <vector>

#include "../fatalerror/exceptions.h"
#include "../fatalerror/gmxassert.h"

#include "abstractoption.h"
#include "abstractoptionstorage.h"

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
        virtual const char *typeString() const = 0;
        virtual int valueCount() const { return static_cast<int>(_values->size()); }
        virtual std::string formatValue(int i) const = 0;

    protected:
        /*! \brief
         * Initializes the storage from option settings.
         *
         * \param[in] settings  Option settings.
         * \param[in] options   Option collection that will contain the
         *     option.
         * \param[in] staticFlags Option flags that are always set and specify
         *      generic behavior of the option.
         * \throws  APIError if invalid settings have been provided.
         */
        template <class U>
        OptionStorageTemplate(const OptionTemplate<T, U> &settings, Options *options,
                              OptionFlags staticFlags = OptionFlags());


        virtual void clearSet();
        /*! \copydoc AbstractOptionStorage::convertValue()
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
        }
        /*! \copydoc AbstractOptionStorage::processSet()
         *
         * OptionStorage template implements transaction support for a set of
         * values in this method (see the class description), and provides a
         * more detailed processSetValues() method that can be overridden in
         * subclasses to process the actual values.  Derived classes should
         * override that method instead of this one if set value processing is
         * necessary.
         */
        virtual void processSet();
        /*! \copydoc AbstractOptionStorage::processAll()
         *
         * The implementation in OptionStorageTemplate does nothing.
         */
        virtual void processAll()
        {
        }

        /*! \brief
         * Removes all values from the storage.
         *
         * Does not throw.
         */
        void clear() { _values->clear(); }
        /*! \brief
         * Adds a value to a temporary storage.
         *
         * \param[in] value  Value to add. A copy is made.
         * \throws InvalidInputError if the maximum value count has been reached.
         *
         * Derived classes should call this function from the convertValue()
         * implementation to add converted values to the storage.
         * If the maximum value cont has been reached, the value is discarded
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
         * Provides derived classes access to the current list of values.
         *
         * The non-const variant should only be used from processAll() in
         * derived classes if necessary, and refreshValues() should be called
         * if any changes are made.
         */
        ValueList &values() { return *_values; }
        //! Provides derived classes access to the current list of values.
        const ValueList &values() const { return *_values; }

    private:
        /*! \brief
         * Vector for temporary storage of values before commitSet() is called.
         */
        ValueList               _setValues;
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
        ValueList              *_values;
        T                      *_store;
        int                    *_countptr;
        // Could be scoped_ptr
        std::auto_ptr<T>        _defaultValueIfSet;

        // Copy and assign disallowed by base.
};


template <typename T>
template <class U>
OptionStorageTemplate<T>::OptionStorageTemplate(const OptionTemplate<T, U> &settings,
                                                Options *options,
                                                OptionFlags staticFlags)
    : AbstractOptionStorage(settings, options, staticFlags),
      _values(settings._storeVector),
      _store(settings._store),
      _countptr(settings._countptr)
{
    std::auto_ptr<std::vector<T> > valueGuard;
    if (!_values)
    {
        // The flag should be set for proper error checking.
        GMX_RELEASE_ASSERT(!hasFlag(efExternalValueVector),
                           "Internal inconsistency");
        valueGuard.reset(new std::vector<T>);
        _values = valueGuard.get();
    }
    if (hasFlag(efNoDefaultValue)
        && (settings._defaultValue != NULL
            || settings._defaultValueIfSet != NULL))
    {
        GMX_THROW(APIError("Option does not support default value, but one is set"));
    }
    if (_store != NULL && _countptr == NULL && !hasFlag(efVector)
        && minValueCount() != maxValueCount())
    {
        GMX_THROW(APIError("Count storage is not set, although the number of produced values is not known"));
    }
    if (!hasFlag(efNoDefaultValue))
    {
        setFlag(efHasDefaultValue);
        if (settings._defaultValue != NULL)
        {
            _values->clear();
            addValue(*settings._defaultValue);
            // TODO: This is a bit hairy, as it indirectly calls a virtual function.
            commitValues();
        }
        else if (!hasFlag(efExternalValueVector) && _store != NULL)
        {
            _values->clear();
            int count = (settings.isVector() ?
                            settings._maxValueCount : settings._minValueCount);
            for (int i = 0; i < count; ++i)
            {
                _values->push_back(_store[i]);
            }
        }
        if (settings._defaultValueIfSet != NULL)
        {
            if (hasFlag(efMulti))
            {
                GMX_THROW(APIError("defaultValueIfSet() is not supported with allowMultiple()"));
            }
            _defaultValueIfSet.reset(new T(*settings._defaultValueIfSet));
        }
    }
    setFlag(efClearOnNextSet);
    valueGuard.release();
}


template <typename T>
OptionStorageTemplate<T>::~OptionStorageTemplate()
{
    if (!hasFlag(efExternalValueVector))
    {
        delete _values;
    }
}


template <typename T>
void OptionStorageTemplate<T>::clearSet()
{
    _setValues.clear();
}


template <typename T>
void OptionStorageTemplate<T>::processSet()
{
    processSetValues(&_setValues);
    if (_setValues.empty() && _defaultValueIfSet.get() != NULL)
    {
        addValue(*_defaultValueIfSet);
        setFlag(efHasDefaultValue);
    }
    else
    {
        clearFlag(efHasDefaultValue);
    }
    if (!hasFlag(efDontCheckMinimumCount)
        && _setValues.size() < static_cast<size_t>(minValueCount()))
    {
        clearSet();
        GMX_THROW(InvalidInputError("Too few (valid) values"));
    }
    commitValues();
}


template <typename T>
void OptionStorageTemplate<T>::addValue(const T &value)
{
    if (maxValueCount() >= 0
        && _setValues.size() >= static_cast<size_t>(maxValueCount()))
    {
        GMX_THROW(InvalidInputError("Too many values"));
    }
    _setValues.push_back(value);
}


template <typename T>
void OptionStorageTemplate<T>::commitValues()
{
    if (hasFlag(efClearOnNextSet))
    {
        _values->swap(_setValues);
        clearFlag(efClearOnNextSet);
    }
    else
    {
        _values->insert(_values->end(), _setValues.begin(), _setValues.end());
    }
    clearSet();
    refreshValues();
}


template <typename T>
void OptionStorageTemplate<T>::refreshValues()
{
    if (_countptr != NULL)
    {
        *_countptr = static_cast<int>(_values->size());
    }
    if (_store != NULL)
    {
        for (size_t i = 0; i < _values->size(); ++i)
        {
            _store[i] = (*_values)[i];
        }
    }
}

} // namespace gmx

#endif
