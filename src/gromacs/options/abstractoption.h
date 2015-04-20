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
/*! \file
 * \brief
 * Defines gmx::AbstractOption, gmx::OptionTemplate and gmx::OptionInfo.
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
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_ABSTRACTOPTION_H
#define GMX_OPTIONS_ABSTRACTOPTION_H

#include <string>
#include <vector>

#include "gromacs/options/optionflags.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class AbstractOptionStorage;
template <typename T> class OptionStorageTemplate;
class OptionManagerContainer;
class Options;

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
 * of storage object.  If they use their own info type derived from OptionInfo,
 * they should also have a public typedef \c InfoType that specifies that
 * info type.  This is required for Options::addOption() to return the correct
 * info type.
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
            : minValueCount_(1), maxValueCount_(1),
              name_(name), descr_(NULL)
        { }

        /*! \brief
         * Creates a default storage object for the option.
         *
         * \param[in] managers  Manager container (unused if the option does
         *     not use a manager).
         * \returns   The created storage object.
         * \throws    APIError if invalid option settings have been provided.
         *
         * This method is called by Options::addOption() when initializing an
         * option from the settings.
         *
         * Derived classes should implement the method to create an actual
         * storage object and populate it with correct values.
         * They should also throw APIError if they detect problems.
         *
         * Should only be called by Options::addOption().
         *
         * The ownership of the return value is passed, but is not using a
         * smart pointer to avoid introducing such a dependency in an installed
         * header.  The implementation will always consist of a single `new`
         * call and returning that value, and the caller always immediately
         * wraps the pointer in a smart pointer, so there is not exception
         * safety issue.
         */
        virtual AbstractOptionStorage *createStorage(
            const OptionManagerContainer &managers) const = 0;

        //! Sets the description for the option.
        void setDescription(const char *descr) { descr_ = descr; }
        //! Sets a flag for the option.
        void setFlag(OptionFlag flag) { flags_.set(flag); }
        //! Clears a flag for the option.
        void clearFlag(OptionFlag flag) { flags_.clear(flag); }
        //! Sets or clears a flag for the option.
        void setFlag(OptionFlag flag, bool bSet) { flags_.set(flag, bSet); }
        //! Returns true if the option is vector-valued.
        bool isVector() const { return hasFlag(efOption_Vector); }
        /*! \brief
         * Sets the option to be vector-valued.
         *
         * This method is provided for convenience to make management of value
         * counts easier.  In order to implement a vector-valued option, the
         * class derived from AbstractOption should expose a method that calls
         * this method, and the storage object derived from
         * AbstractOptionStorage should check isVector().
         * If only a single value is provided, the storage object should fill
         * the whole vector with that value.
         *
         * The length of the vector (the value of maxValueCount_) must be
         * fixed.  The default length is 3 elements.
         */
        void setVector()
        {
            setFlag(efOption_Vector);
            minValueCount_ = 1;
            if (maxValueCount_ == 1)
            {
                maxValueCount_ = 3;
            }
        }
        //! Sets the required number of values for the option.
        void setValueCount(int count)
        {
            if (!hasFlag(efOption_Vector))
            {
                minValueCount_ = count;
            }
            maxValueCount_ = count;
        }

        //! Minimum number of values required for the option.
        int                     minValueCount_;
        //! Maximum number of values allowed for the option.
        int                     maxValueCount_;
        //! \endcond

    private:
        //! Returns true if a flag has been set.
        bool hasFlag(OptionFlag flag) const { return flags_.test(flag); }

        const char             *name_;
        //! Pointer to description of the option.
        const char             *descr_;
        OptionFlags             flags_;

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
        { setFlag(efOption_Hidden, bHidden); return me(); }
        /*! \brief
         * Requires the option to be specified explicitly.
         *
         * Note that if you specify defaultValue() together with required(),
         * the user is not required to explicitly provide the option.
         * In this case, required() only affects possible help output.
         */
        MyClass &required(bool bRequired = true)
        { setFlag(efOption_Required, bRequired); return me(); }
        //! Allows the option to be specified multiple times.
        MyClass &allowMultiple(bool bMulti = true)
        { setFlag(efOption_MultipleTimes, bMulti); return me(); }
        //! Requires exactly \p count values for the option.
        MyClass &valueCount(int count) { setValueCount(count); return me(); }
        //! Allows any number of values for the option.
        MyClass &multiValue(bool bMulti = true)
        { if (bMulti) { maxValueCount_ = -1; } return me(); }

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
         *
         * \p defaultValue is copied when the option is created.
         */
        MyClass &defaultValue(const T &defaultValue)
        { defaultValue_ = &defaultValue; return me(); }
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
        { defaultValueIfSet_ = &defaultValue; return me(); }
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
        { store_ = store; return me(); }
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
        { countptr_ = countptr; return me(); }
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
        { storeVector_ = store; return me(); }

    protected:
        /*! \cond libapi */
        //! Alias for the template class for use in base classes.
        typedef OptionTemplate<T, U> MyBase;

        //! Initializes the name and default values for an option.
        explicit OptionTemplate(const char *name)
            : AbstractOption(name),
              defaultValue_(NULL), defaultValueIfSet_(NULL), store_(NULL),
              countptr_(NULL), storeVector_(NULL)
        { }

        /*! \brief
         * Returns a pointer to user-specified default value, or NULL if there
         * is none.
         */
        const T *defaultValue() const { return defaultValue_; }
        /*! \brief
         * Returns a pointer to user-specified default value, or NULL if there
         * is none.
         */
        const T *defaultValueIfSet() const { return defaultValueIfSet_; }
        //! Returns \p *this casted into MyClass to reduce typing.
        MyClass &me() { return static_cast<MyClass &>(*this); }
        //! \endcond

    private:
        const T                *defaultValue_;
        const T                *defaultValueIfSet_;
        T                      *store_;
        int                    *countptr_;
        std::vector<T>         *storeVector_;

        /*! \brief
         * Needed to initialize storage from this class without otherwise
         * unnecessary accessors.
         */
        friend class OptionStorageTemplate<T>;
};

/*! \brief
 * Gives information and allows modifications to an option after creation.
 *
 * When an option is added with Options::addOption(), an object of a subclass
 * of OptionInfo is returned.  This object can be later used to access
 * information about the option.  Non-const methods also allow later changing
 * (some of) the option settings provided at initialization time.
 * The properties accessible/modifiable through this interface are implemented
 * based on need, and may not be implemented for all cases.
 *
 * \if libapi
 * This class is also used by OptionsVisitor and OptionsModifyingVisitor as
 * the interface that allows querying/modifying each visited option.
 * \endif
 *
 * This class isolates the details of the internal option implementation from
 * callers.  Although this class is a simple reference to the underlying
 * implementation, it is implemented as non-copyable to allow const/non-const
 * status of a reference to this class to indicate whether modifications are
 * allowed.  Otherwise, separate classes would be needed for access and
 * modification, complicating the implementation.  In the implementation,
 * there is always a single OptionInfo instance referring to one option.
 * The underlying implementation object always owns this instance, and only
 * references are passed to callers.
 *
 * \see Options::addOption()
 * \if libapi
 * \see OptionsVisitor
 * \see OptionsModifyingVisitor
 * \endif
 *
 * \inpublicapi
 * \ingroup module_options
 */
class OptionInfo
{
    public:
        virtual ~OptionInfo();

        /*! \brief
         * Test whether the option is of a particular type.
         *
         * \tparam InfoType  Option type to test for. Should be a class derived
         *      from OptionInfo.
         */
        template <class InfoType>
        bool isType() const
        {
            return toType<InfoType>() != NULL;
        }
        /*! \brief
         * Convert the info object to a particular type if the type is correct.
         *
         * \tparam InfoType  Option type to convert to. Should be a class
         *      derived from OptionInfo.
         * \retval this converted to a pointer to \p InfoType, or NULL if the
         *      conversion is not possible.
         */
        template <class InfoType>
        InfoType *toType()
        {
            return dynamic_cast<InfoType *>(this);
        }
        //! \copydoc toType()
        template <class InfoType>
        const InfoType *toType() const
        {
            return dynamic_cast<const InfoType *>(this);
        }

        //! Returns true if the option has been set.
        bool isSet() const;
        //! Returns true if the option is a hidden option.
        bool isHidden() const;
        //! Returns true if the option is required.
        bool isRequired() const;
        //! Returns the minimum number of values that this option accepts.
        int minValueCount() const;
        //! Returns the maximum number of values that this option accepts.
        int maxValueCount() const;
        //! Returns the name of the option.
        const std::string &name() const;
        //! Returns the type of the option as a string.
        std::string type() const;
        //! Returns the description of the option.
        std::string formatDescription() const;
        /*! \brief
         * Returns the default value if set for the option as a string.
         *
         * \see OptionTemplate::defaultValueIfSet()
         */
        std::string formatDefaultValueIfSet() const;

        //! Returns the number of values given for the option.
        int valueCount() const;
        //! Returns the i'th value of the option as a string.
        std::string formatValue(int i) const;

    protected:
        /*! \cond libapi */
        /*! \brief
         * Wraps a given option object.
         *
         * Does not throw.
         */
        explicit OptionInfo(AbstractOptionStorage *option);

        //! Returns the wrapped option storage object.
        AbstractOptionStorage       &option() { return option_; }
        //! Returns the wrapped option storage object.
        const AbstractOptionStorage &option() const { return option_; }
        //! \endcond

    private:
        //! The wrapped option.
        AbstractOptionStorage  &option_;

        GMX_DISALLOW_COPY_AND_ASSIGN(OptionInfo);
};

} // namespace gmx

#endif
