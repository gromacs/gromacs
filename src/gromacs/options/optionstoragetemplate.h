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

#include <cassert>

#include <string>
#include <vector>

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

        /*! \brief
         * Initializes the storage from option settings.
         *
         * \retval 0 on success.
         *
         * \see OptionTemplate::createDefaultStorage()
         */
        template <class U>
        int init(const OptionTemplate<T, U> &settings, Options *options);

        // No implementation in this class for the pure virtual methods, but
        // the declarations are still included for clarity.
        virtual const char *typeString() const = 0;
        virtual int valueCount() const { return _values->size(); }
        virtual std::string formatValue(int i) const = 0;

    protected:
        //! Initializes default values (no storage).
        OptionStorageTemplate();

        virtual void clear();
        /*! \copydoc AbstractOptionStorage::convertValue()
         *
         * Derived classes should call addValue() after they have converted
         * \p value to the storage type.
         */
        virtual int convertValue(const std::string &value,
                                 AbstractErrorReporter *errors) = 0;
        /*! \copydoc AbstractOptionStorage::processSet()
         *
         * The implementation in OptionStorageTemplate copies the values
         * from the main storage vector to alternate locations, and always
         * succeeds.  Derived classes should always call the base class
         * implementation if they override this method.
         */
        virtual int processSet(int nvalues,
                               AbstractErrorReporter * /*errors*/)
        {
            processValues(nvalues, true);
            return 0;
        }
        /*! \copydoc AbstractOptionStorage::processAll()
         *
         * The implementation in OptionStorageTemplate does nothing, and always
         * returns zero.  Derived classes should still always call the base
         * class implementation if they override this method.
         */
        virtual int processAll(AbstractErrorReporter * /*errors*/)
        { return 0; }

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
        int addValue(const T &value);
        /*! \brief
         * Store values in alternate locations.
         *
         * \param[in] nvalues  Number of values to process.
         * \param[in] bDoArray Whether to put values in the array storage as
         *      well.
         *
         * Stores the last \p nvalues values added with addValue() to the
         * alternate storage locations.
         * The current implementation asserts if it is called more than once
         * with \p bDoArray set to true.
         *
         * Derived classes should call this method if they use addValue()
         * outside convertValue(), e.g., to set a default value.  In such
         * cases, \p bDoArray should be set to false.
         */
        void processValues(int nvalues, bool bDoArray);

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
        T                     **_storeArray;
        int                    *_nvalptr;

        // Copy and assign disallowed by base.
};

/*! \brief
 * Helper function for creating storage objects.
 *
 * \tparam     T  Type of the settings object (derived from OptionTemplate).
 * \tparam     S  Type of the storage object (derived from OptionStorageTemplate).
 * \param[in]  settings  Settings object to pass to S::init().
 * \param[in]  options   Options object to pass to S::init().
 * \param[out] output    Pointer to the created storage object.
 * \returns    The return value of S::init().
 *
 * Creates a new instance of S and calls the init() method with the provided
 * parameters.  If the initialization fails, destroys the partially constructed
 * object.
 *
 * \inlibraryapi
 */
template <class T, class S> int
createOptionStorage(const T *settings, Options *options,
                    AbstractOptionStorage **output)
{
    S *storage = new S;
    int rc = storage->init(*settings, options);
    if (rc != 0)
    {
        delete storage;
        return rc;
    }
    *output = storage;
    return 0;
}


template <typename T>
OptionStorageTemplate<T>::OptionStorageTemplate()
    : _values(NULL), _store(NULL), _storeArray(NULL), _nvalptr(NULL)
{
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
template <class U>
int OptionStorageTemplate<T>::init(const OptionTemplate<T, U> &settings, Options *options)
{
    // It's impossible for the caller to do proper memory management if
    // the provided memory is not initialized as NULL.
    assert(settings._storeArray == NULL || *settings._storeArray == NULL);
    int rc = AbstractOptionStorage::init(settings, options);
    if (rc != 0)
    {
        return rc;
    }
    _store      = settings._store;
    _storeArray = settings._storeArray;
    _nvalptr    = settings._nvalptr;
    _values     = settings._storeVector;
    if (!_values)
    {
        // The flag should be set for proper error checking.
        assert(!hasFlag(efExternalValueVector));
        _values = new std::vector<T>;
    }
    // If the option does not support default values, one should not be set.
    assert(!hasFlag(efNoDefaultValue) || settings._defaultValue == NULL);
    if (!hasFlag(efNoDefaultValue))
    {
        if (settings._defaultValue != NULL)
        {
            _values->clear();
            addValue(*settings._defaultValue);
            processValues(1, false);
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
    }
    return 0;
}


template <typename T>
void OptionStorageTemplate<T>::clear()
{
    _values->clear();
}


template <typename T>
int OptionStorageTemplate<T>::addValue(const T &value)
{
    int rc = incrementValueCount();
    if (rc == 0)
    {
        _values->push_back(value);
    }
    return rc;
}


template <typename T>
void OptionStorageTemplate<T>::processValues(int nvalues, bool bDoArray)
{
    if (_nvalptr)
    {
        *_nvalptr = _values->size();
    }
    if (bDoArray && _storeArray != NULL)
    {
        assert(*_storeArray == NULL);
        *_storeArray = new T[_values->size()];
    }
    for (size_t i = _values->size() - nvalues; i < _values->size(); ++i)
    {
        if (_store)
        {
            _store[i] = (*_values)[i];
        }
        if (bDoArray && _storeArray != NULL)
        {
            (*_storeArray)[i] = (*_values)[i];
        }
    }
}

} // namespace gmx

#endif
