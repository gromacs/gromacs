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
/*! \file
 * \brief
 * Defines gmx::AbstractOption and a related template.
 *
 * This header defines base classes for option settings that are used with
 * Options::addOption().  These classes implement the "named parameter"
 * idiom for specifying option properties.
 *
 * These classes also take care of creating and setting up the actual option
 * objects.
 *
 * This header is needed directly only when implementing new option types,
 * but methods of OptionTemplate are visible even to the normal user through
 * its subclasses.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_ABSTRACTOPTION_H
#define GMX_OPTIONS_ABSTRACTOPTION_H

#include <string>
#include <vector>

#include "../utility/uniqueptr.h"

#include "optionflags.h"

namespace gmx
{

class AbstractOptionStorage;
template <typename T> class OptionStorageTemplate;
class Options;

//! Smart pointer for managing an AbstractOptionStorage object.
typedef gmx_unique_ptr<AbstractOptionStorage>::type
        AbstractOptionStoragePointer;

/*! \brief
 * Abstract base class for specifying option properties.
 *
 * Concrete classes should normally not derive directly from this class,
 * but from OptionTemplate instead.  Classes derived from this class
 * are mainly designed to implement the "named parameter" idiom.  For
 * efficiency and clarity, these classes should only store values provided to
 * them.  All error checking and memory management should be postponed to the
 * point when the actual option is created.
 *
 * Subclasses should override createStorage() to create the correct type
 * of storage object.
 *
 * \ingroup module_options
 */
class AbstractOption
{
    public:
        // Virtual only for completeness, in normal use should not be needed.
        virtual ~AbstractOption() { }

    protected:
        /*! \cond libapi */
        //! Initializes the name and default values for an option.
        explicit AbstractOption(const char *name)
            : _minValueCount(1), _maxValueCount(1),
              _name(name), _descr(NULL)
        { }

        /*! \brief
         * Creates a default storage object for the option.
         *
         * \returns The created storage object.
         * \throws  APIError if invalid option settings have been provided.
         *
         * This method is called by Options::addOption() when initializing an
         * option from the settings.
         *
         * Derived classes should implement the method to create an actual
         * storage object and populate it with correct values.
         * They should also throw APIError if they detect problems.
         *
         * Should only be called by Options::addOption().
         */
        virtual AbstractOptionStoragePointer createStorage() const = 0;

        /*! \brief
         * Creates the description string for the option.
         *
         * \returns Description string for the option.
         *
         * This function is virtual to allow derived classes to customize the
         * description programmatically, e.g., by adding the list of allowed
         * values.
         * The default implementation simply returns the user-provided
         * description.
         */
        virtual std::string createDescription() const
        { return _descr ? _descr : ""; }

        //! Sets the description for the option.
        void setDescription(const char *descr) { _descr = descr; }
        //! Sets a flag for the option.
        void setFlag(OptionFlag flag) { _flags.set(flag); }
        //! Clears a flag for the option.
        void clearFlag(OptionFlag flag) { _flags.clear(flag); }
        //! Sets or clears a flag for the option.
        void setFlag(OptionFlag flag, bool bSet) { _flags.set(flag, bSet); }
        //! Returns true if the option is vector-valued.
        bool isVector() const { return hasFlag(efVector); }
        //! Sets the option to be vector-valued.
        void setVector()
        {
            setFlag(efVector);
            _minValueCount = 1;
            if (_maxValueCount == 1)
            {
                _maxValueCount = 3;
            }
        }
        //! Sets the required number of values for the option.
        void setValueCount(int count)
        {
            if (!hasFlag(efVector))
            {
                _minValueCount = count;
            }
            _maxValueCount = count;
        }

        //! Minimum number of values required for the option.
        int                     _minValueCount;
        //! Maximum number of values allowed for the option.
        int                     _maxValueCount;
        //! \endcond

    private:
        //! Returns true if a flag has been set.
        bool hasFlag(OptionFlag flag) const { return _flags.test(flag); }

        const char             *_name;
        //! Pointer to description of the option.
        const char             *_descr;
        OptionFlags             _flags;

        /*! \brief
         * Needed to initialize an AbstractOptionStorage object from this class
         * without otherwise unnecessary accessors.
         */
        friend class AbstractOptionStorage;
        /*! \brief
         * Needed to be able to call createStorage().
         */
        friend class Options;
};

/*! \brief
 * Templated base class for constructing concrete option settings classes.
 *
 * \tparam T Assignable type that stores a single option value.
 * \tparam U Type of the derived class.
 *
 * This template is used as a base class like this:
 * \code
class ConcreteOption : public OptionTemplate<int, ConcreteOption>
{
 * \endcode
 *
 * All public functions in this class return \c *this casted to a reference to
 * \p U.  They do not throw.
 *
 * For examples of how to use classes derived from this class, see the class
 * documentation for Options.
 *
 * \inlibraryapi
 * \ingroup module_options
 */
template <typename T, class U>
class OptionTemplate : public AbstractOption
{
    public:
        //! Type that stores a single option value.
        typedef T ValueType;
        //! Alias for the derived class type.
        typedef U MyClass;

        /*! \brief
         * Sets a description for the option.
         *
         * \param[in] descr Description to set.
         *
         * String in \p descr is copied when the option is created.
         */
        MyClass &description(const char *descr)
        { setDescription(descr); return me(); }
        //! Hides the option from normal help output.
        MyClass &hidden(bool bHidden = true)
        { setFlag(efHidden, bHidden); return me(); }
        //! Requires the option to be specified explicitly.
        MyClass &required(bool bRequired = true)
        { setFlag(efRequired, bRequired); return me(); }
        //! Allows the option to be specified multiple times.
        MyClass &allowMultiple(bool bMulti = true)
        { setFlag(efMulti, bMulti); return me(); }
        //! Requires exactly \p count values for the option.
        MyClass &valueCount(int count) { setValueCount(count); return me(); }
        //! Allows any number of values for the option.
        MyClass &multiValue() { _maxValueCount = -1; return me(); }

        /*! \brief
         * Sets a default value for the option.
         *
         * \param[in] defaultValue Default value.
         *
         * If the option is never set, the default value is copied to the
         * assigned storage.  Note that if the option is not set and there
         * is no default value, the storage is not altered, which can also be
         * used to provide a default value.  The latter method has to be used
         * if the option can take multiple values.
         * If required() is specified, only affects the default value shown in
         * help output.
         *
         * \p defaultValue is copied when the option is created.
         */
        MyClass &defaultValue(const T &defaultValue)
        { _defaultValue = &defaultValue; return me(); }
        /*! \brief
         * Sets a default value for the option when it is set.
         *
         * \param[in] defaultValue Default value.
         *
         * This value is used if the option is set, but no value is provided.
         * If the option is never set, the value set with defaultValue() is
         * used.  Can only be used for options that accept a single value.
         *
         * \p defaultValue is copied when the option is created.
         */
        MyClass &defaultValueIfSet(const T &defaultValue)
        { _defaultValueIfSet = &defaultValue; return me(); }
        /*! \brief
         * Stores value(s) in memory pointed by \p store.
         *
         * \param[in] store  Storage for option value(s).
         *
         * The caller is responsible for allocating enough memory such that
         * the any allowed number of values fits into the array pointed by
         * \p store.  If there is no maximum allowed number or if the maximum
         * is inconveniently large, storeVector() should be used.
         *
         * For information on when values are available in the storage, see
         * storeVector().
         *
         * The pointer provided should remain valid as long as the associated
         * Options object exists.
         */
        MyClass &store(T *store)
        { setFlag(efExternalStore); _store = store; return me(); }
        /*! \brief
         * Stores number of values in the value pointed by \p countptr.
         *
         * \param[in] countptr Storage for the number of values.
         *
         * For information on when values are available in the storage, see
         * storeVector().
         *
         * The pointers provided should remain valid as long as the associated
         * Options object exists.
         */
        MyClass &storeCount(int *countptr)
        { _countptr = countptr; return me(); }
        /*! \brief
         * Stores option values in the provided vector.
         *
         * \param[in] store  Vector to store option values in.
         *
         * Values are added to the vector after each successful set of values
         * is parsed.  Note that for some options, the value may be changed
         * later, and is only guaranteed to be correct after Options::finish()
         * has been called.
         *
         * The pointer provided should remain valid as long as the associated
         * Options object exists.
         */
        MyClass &storeVector(std::vector<T> *store)
        { setFlag(efExternalValueVector); _storeVector = store; return me(); }

    protected:
        /*! \cond libapi */
        //! Alias for the template class for use in base classes.
        typedef OptionTemplate<T, U> MyBase;

        //! Initializes the name and default values for an option.
        explicit OptionTemplate(const char *name)
            : AbstractOption(name),
              _defaultValue(NULL), _defaultValueIfSet(NULL), _store(NULL),
              _countptr(NULL), _storeVector(NULL)
        { }

        /*! \brief
         * Returns a pointer to user-specified default value, or NULL if there
         * is none.
         */
        const T *defaultValue() const { return _defaultValue; }
        /*! \brief
         * Returns a pointer to user-specified default value, or NULL if there
         * is none.
         */
        const T *defaultValueIfSet() const { return _defaultValueIfSet; }
        //! Returns \p *this casted into MyClass to reduce typing.
        MyClass &me() { return static_cast<MyClass &>(*this); }
        //! \endcond

    private:
        const T                *_defaultValue;
        const T                *_defaultValueIfSet;
        T                      *_store;
        int                    *_countptr;
        std::vector<T>         *_storeVector;

        /*! \brief
         * Needed to initialize storage from this class without otherwise
         * unnecessary accessors.
         */
        friend class OptionStorageTemplate<T>;
};

} // namespace gmx

#endif
