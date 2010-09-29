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
 * Provides an implementation of the clear(), finish(), valueCount(), and
 * formatValues() methods of the AbstractOptionStorage interface, leaving
 * typeString(), appendValue(), finishSet(), and formatValue() to be
 * implemented in derived classes.
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
        virtual void clear();
        /*! \copydoc AbstractOptionStorage::appendValue()
         *
         * Derived classes should call addValue() after they have converted
         * \p value to the storage type.
         */
        virtual int appendValue(const std::string &value,
                                AbstractErrorReporter *errors) = 0;
        virtual int finishSet(int nvalues,
                              AbstractErrorReporter *errors) = 0;
        /*! \copydoc AbstractOptionStorage::finish()
         *
         * The implementation in AbstractStorage copies the values
         * from the main storage vector to alternate locations, and always
         * succeeds.  Derived classes should always call the base class
         * implementation if they override this method.
         */
        virtual int finish(AbstractErrorReporter *errors);
        virtual int valueCount() const { return _values->size(); }
        virtual std::string formatValue(int i) const = 0;
        virtual std::string formatValues() const;

    protected:
        //! Initializes default values (no storage).
        OptionStorageTemplate();

        /*! \brief
         * Adds a value to the storage.
         *
         * \param[in] value  Value to add. A copy is made.
         *
         * Derived classes should call this function from the appendValue()
         * implementation to add converted values to the storage.
         */
        void addValue(const T &value);

        //! Returns the Options object that houses the parameter.
        Options &hostOptions() { return *_options; }
        //! \copydoc hostOptions()
        const Options &hostOptions() const { return *_options; }

        //! Provides derived classes access to the current list of values.
        ValueList &values() { return *_values; }
        //! Provides derived classes access to the current list of values.
        const ValueList &values() const { return *_values; }

    private:
        //! Internal flags for managing the storage.
        enum Flags {
            //! Memory for \p _values is managed by this object.
            efOwnValueVector   = 1<<0,
        };

        /*! \brief
         * Vector for primary storage of option values.
         *
         * Is never NULL; points either to externally provided vector, or an
         * internally allocated one.  The allocation is performed by init().
         *
         * addValue() adds values only to this storage.  Other memory locations
         * are updated only when finish() is called.
         */
        ValueList              *_values;
        T                      *_store;
        T                     **_storeArray;
        int                    *_nvalptr;
        int                     _flags;
        Options                *_options;

        // Disallow copy and assign.
        OptionStorageTemplate(const OptionStorageTemplate<T> &);
        void operator =(const OptionStorageTemplate<T> &);
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
    : _values(NULL), _store(NULL), _storeArray(NULL), _nvalptr(NULL),
      _flags(0), _options(NULL)
{
}

template <typename T>
OptionStorageTemplate<T>::~OptionStorageTemplate()
{
    if (_flags & efOwnValueVector)
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
    _store      = settings._store;
    _storeArray = settings._storeArray;
    _nvalptr    = settings._nvalptr;
    _values     = settings._storeVector;
    _options    = options;
    if (!_values)
    {
        _values = new std::vector<T>;
        _flags |= efOwnValueVector;
    }
    if (settings._defaultValue)
    {
        _values->clear();
        addValue(*settings._defaultValue);
    }
    else if ((_flags & efOwnValueVector) && _store)
    {
        _values->clear();
        int count = (settings.isVector() ?
                        settings._maxValueCount : settings._minValueCount);
        for (int i = 0; i < count; ++i)
        {
            _values->push_back(_store[i]);
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
int OptionStorageTemplate<T>::finish(AbstractErrorReporter * /*errors*/)
{
    if (_nvalptr)
    {
        *_nvalptr = _values->size();
    }
    if (_storeArray)
    {
        *_storeArray = new T[_values->size()];
    }
    for (size_t i = 0; i < _values->size(); ++i)
    {
        if (_store)
        {
            _store[i] = (*_values)[i];
        }
        if (_storeArray)
        {
            (*_storeArray)[i] = (*_values)[i];
        }
    }
    return 0;
}

template <typename T>
void OptionStorageTemplate<T>::addValue(const T &value)
{
    // We store the value here to allow limited dependence between option
    // values.
    if (_store != NULL)
    {
        _store[_values->size()] = value;
    }
    _values->push_back(value);
}

template <typename T>
std::string OptionStorageTemplate<T>::formatValues() const
{
    std::string result;
    int count = valueCount();
    for (int i = 0; i < count; ++i)
    {
        if (i != 0)
        {
            result.append(" ");
        }
        result.append(formatValue(i));
    }
    return result;
}

} // namespace gmx

#endif
