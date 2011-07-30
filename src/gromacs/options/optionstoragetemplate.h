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

#include <string>
#include <vector>

#include "../fatalerror/exceptions.h"
#include "../fatalerror/gmxassert.h"

#include "abstractoption.h"
#include "abstractoptionstorage.h"

namespace gmx
{

class Options;

/*! \brief
 * Templated base class for constructing option value storage classes.
 *
 * \tparam T Assignable type that stores a single option value.
 *
 * Provides an implementation of the clearSet() and valueCount() methods of
 * AbstractOptionStorage, as well as a basic implementation of processSet() and
 * processAll().  This leaves typeString(), formatValue(), and convertValue()
 * to be implemented in derived classes.
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
        virtual int valueCount() const { return _values->size(); }
        virtual std::string formatValue(int i) const = 0;

    protected:
        /*! \brief
         * Initializes the storage from option settings.
         *
         * \retval 0 on success.
         *
         * \see OptionTemplate::createDefaultStorage()
         */
        template <class U>
        OptionStorageTemplate(const OptionTemplate<T, U> &settings, Options *options,
                              OptionFlags staticFlags = OptionFlags());


        virtual void clearSet();
        /*! \copydoc AbstractOptionStorage::convertValue()
         *
         * Derived classes should call addValue() after they have converted
         * \p value to the storage type.
         */
        virtual void convertValue(const std::string &value) = 0;
        /*! \brief
         * Processes values for a set after all have been converted.
         *
         * \param[in,out] values Valid values in the set.
         *
         * This method is called after all convertValue() calls for a set.
         * \p values contains all values that were validly converted by
         * convertValue().  The derived class may alter the values.
         */
        virtual void processSetValues(ValueList *values)
        {
        }
        /*! \copydoc AbstractOptionStorage::processSet()
         *
         * OptionStorage template implements this method, and provides a more
         * detailed processSetValues() method that can be overridden in
         * subclasses.
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
         */
        void clear() { _values->clear(); }
        /*! \brief
         * Adds a value to the storage.
         *
         * \param[in] value  Value to add. A copy is made.
         * \retval 0 on success.
         * \retval ::eeInvalidInput if the maximum value count has been reached.
         *
         * Derived classes should call this function from the convertValue()
         * implementation to add converted values to the storage.
         */
        void addValue(const T &value);
        /*! \brief
         * Commits values added with addValue().
         *
         * Derived classes should call this method if they use addValue()
         * outside convertValue(), e.g., to set a default value.
         */
        void commitValues();
        /*! \brief
         * Store values in alternate locations.
         */
        virtual void refreshValues();

        //! Provides derived classes access to the current list of values.
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
         * internally allocated one.  The allocation is performed by init().
         *
         * addValue() adds values only to this storage.  Other memory locations
         * are updated only when processAll() is called.
         */
        ValueList              *_values;
        T                      *_store;
        int                    *_countptr;
        T                      *_defaultValueIfSet;

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
      _countptr(settings._countptr),
      _defaultValueIfSet(NULL)
{
    if (!_values)
    {
        // The flag should be set for proper error checking.
        GMX_RELEASE_ASSERT(!hasFlag(efExternalValueVector),
                           "Internal inconsistency");
        _values = new std::vector<T>;
    }
    GMX_RELEASE_ASSERT(!hasFlag(efNoDefaultValue)
                       || (settings._defaultValue == NULL
                           && settings._defaultValueIfSet == NULL),
                       "Option does not support default value, but one is set");
    if (_store != NULL && _countptr == NULL && !hasFlag(efVector)
        && minValueCount() != maxValueCount())
    {
        GMX_THROW(APIError("Count storage is not set, although the number of produced values is not known"));
    }
    if (!hasFlag(efNoDefaultValue))
    {
        if (settings._defaultValue != NULL)
        {
            _values->clear();
            addValue(*settings._defaultValue);
            setFlag(efHasDefaultValue);
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
            setFlag(efHasDefaultValue);
        }
        if (settings._defaultValueIfSet != NULL)
        {
            _defaultValueIfSet = new T(*settings._defaultValueIfSet);
        }
    }
}


template <typename T>
OptionStorageTemplate<T>::~OptionStorageTemplate()
{
    if (!hasFlag(efExternalValueVector))
    {
        delete _values;
    }
    delete _defaultValueIfSet;
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
    if (_setValues.empty() && _defaultValueIfSet != NULL)
    {
        addValue(*_defaultValueIfSet);
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
    if (hasFlag(efHasDefaultValue))
    {
        _values->swap(_setValues);
        clearFlag(efHasDefaultValue);
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
        *_countptr = _values->size();
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
