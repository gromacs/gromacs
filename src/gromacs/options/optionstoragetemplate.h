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
 * Provides an implementation of the clear() and valueCount() methods of
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


        virtual void clear();
        /*! \copydoc AbstractOptionStorage::convertValue()
         *
         * Derived classes should call addValue() after they have converted
         * \p value to the storage type.
         */
        virtual void convertValue(const std::string &value) = 0;
        /*! \copydoc AbstractOptionStorage::processSet()
         *
         * The implementation in OptionStorageTemplate copies the values
         * from the main storage vector to alternate locations, and always
         * succeeds.  Derived classes should always call the base class
         * implementation if they override this method.
         */
        virtual void processSet(int nvalues)
        {
            processValues(nvalues);
        }
        /*! \copydoc AbstractOptionStorage::processAll()
         *
         * The implementation in OptionStorageTemplate does nothing.
         * Derived classes should still always call the base class
         * implementation if they override this method.
         */
        virtual void processAll()
        {
        }

        /*! \brief
         * Adds a value to the storage.
         *
         * \param[in] value  Value to add. A copy is made.
         * \retval 0 on success.
         * \retval ::eeInvalidInput if the maximum value count has been reached.
         *
         * Derived classes should call this function from the convertValue()
         * implementation to add converted values to the storage.
         * It is only necessary to check the return value if addValue() is
         * called more than once from one convertValue() invocation, or if
         * ::efConversionMayNotAddValues is specified.
         */
        void addValue(const T &value);
        /*! \brief
         * Store values in alternate locations.
         *
         * \param[in] nvalues  Number of values to process.
         *
         * Stores the last \p nvalues values added with addValue() to the
         * alternate storage locations.
         *
         * Derived classes should call this method if they use addValue()
         * outside convertValue(), e.g., to set a default value.
         */
        void processValues(int nvalues);

        //! Provides derived classes access to the current list of values.
        ValueList &values() { return *_values; }
        //! Provides derived classes access to the current list of values.
        const ValueList &values() const { return *_values; }

    private:
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
void OptionStorageTemplate<T>::clear()
{
    _values->clear();
}


template <typename T>
void OptionStorageTemplate<T>::addValue(const T &value)
{
    incrementValueCount();
    _values->push_back(value);
    if (_store != NULL)
    {
        _store[_values->size() - 1] = value;
    }
    if (_countptr != NULL)
    {
        *_countptr = _values->size();
    }
}


template <typename T>
void OptionStorageTemplate<T>::processValues(int nvalues)
{
    if (nvalues == 0 && _defaultValueIfSet != NULL)
    {
        addValue(*_defaultValueIfSet);
        nvalues = 1;
    }
    if (_countptr != NULL)
    {
        *_countptr = _values->size();
    }
    if (_store != NULL)
    {
        for (size_t i = _values->size() - nvalues; i < _values->size(); ++i)
        {
            _store[i] = (*_values)[i];
        }
    }
}

} // namespace gmx

#endif
