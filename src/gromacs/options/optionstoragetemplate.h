/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
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

#include <memory>
#include <string>
#include <vector>

#include "gromacs/options/abstractoption.h"
#include "gromacs/options/abstractoptionstorage.h"
#include "gromacs/utility/any.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "valueconverter.h"
#include "valuestore.h"

namespace gmx
{

class Options;

/*! \libinternal \brief
 * Templated base class for constructing option value storage classes.
 *
 * \tparam T Assignable type that stores a single option value.
 *
 * Provides an implementation of the clearSet(), valueCount(), processSet(),
 * and defaultValuesAsStrings() methods of AbstractOptionStorage, as well as a
 * basic no-action implementation of processAll().  Two new virtual methods are
 * added: processSetValues() and formatSingleValue().
 * This leaves typeString(), convertValue() and formatStringValue() to be
 * implemented in derived classes.
 * processSetValues() and processAll() can also be implemented if necessary.
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
template<typename T>
class OptionStorageTemplate : public AbstractOptionStorage
{
public:
    //! Alias for the template class for use in base classes.
    typedef OptionStorageTemplate<T> MyBase;
    //! Type of the container that contains the current values.
    typedef std::vector<T> ValueList;

    // No implementation in this class for the pure virtual methods, but
    // the declarations are still included for clarity.
    // The various copydoc calls are needed with Doxygen 1.8.10, although
    // things work without with 1.8.5...
    std::string typeString() const override = 0;
    //! \copydoc gmx::AbstractOptionStorage::valueCount()
    int valueCount() const override { return store_->valueCount(); }
    //! \copydoc gmx::AbstractOptionStorage::defaultValues()
    std::vector<Any> defaultValues() const override;
    /*! \copydoc gmx::AbstractOptionStorage::defaultValuesAsStrings()
     *
     * OptionStorageTemplate implements handling of defaultValueIfSet()
     * cases and composing the vector.
     * Derived classes must implement formatSingleValue() to provide the
     * actual formatting for a value of type \p T.
     */
    std::vector<std::string> defaultValuesAsStrings() const override;

protected:
    //! Smart pointer for managing the final storage interface.
    typedef std::unique_ptr<IOptionValueStore<T>> StorePointer;

    /*! \brief
     * Initializes the storage from option settings.
     *
     * \param[in] settings  Option settings.
     * \param[in] staticFlags Option flags that are always set and specify
     *      generic behavior of the option.
     * \throws  APIError if invalid settings have been provided.
     */
    template<class U>
    explicit OptionStorageTemplate(const OptionTemplate<T, U>& settings,
                                   OptionFlags                 staticFlags = OptionFlags());
    /*! \brief
     * Initializes the storage from base option settings.
     *
     * \param[in] settings  Option settings.
     * \param[in] store     Final storage location.
     * \throws  APIError if invalid settings have been provided.
     *
     * This constructor works for cases where there is no matching
     * OptionTemplate (e.g., EnumOption).
     */
    OptionStorageTemplate(const AbstractOption& settings, StorePointer store);

    //! \copydoc gmx::AbstractOptionStorage::clearSet()
    void clearSet() override;
    /*! \copydoc gmx::AbstractOptionStorage::convertValue()
     *
     * Derived classes should call addValue() after they have converted
     * \p value to the storage type.  It is allowed to call addValue()
     * more than once, or not at all.  OptionsAssigner::appendValue()
     * provides the same exception safety guarantee as this method, so it
     * should be considered whether the implementation can be made strongly
     * exception safe.
     */
    void convertValue(const Any& value) override = 0;
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
    virtual void processSetValues(ValueList* values) { GMX_UNUSED_VALUE(values); }
    /*! \copydoc gmx::AbstractOptionStorage::processSet()
     *
     * OptionStorageTemplate implements transaction support for a set of
     * values in this method (see the class description), and provides a
     * more detailed processSetValues() method that can be overridden in
     * subclasses to process the actual values.  Derived classes should
     * override that method instead of this one if set value processing is
     * necessary.
     */
    void processSet() override;
    /*! \copydoc gmx::AbstractOptionStorage::processAll()
     *
     * The implementation in OptionStorageTemplate does nothing.
     */
    void processAll() override {}
    /*! \brief
     * Formats a single value as a string.
     *
     * \param[in] value  Value to format.
     * \returns   \p value formatted as a string.
     *
     * The derived class must provide this method to format values a
     * strings.  Called by defaultValuesAsStrings() to do the actual
     * formatting.
     */
    virtual std::string formatSingleValue(const T& value) const = 0;

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
    void addValue(const T& value);
    /*! \brief
     * Commits values added with addValue().
     *
     * \throws std::bad_alloc if out of memory.
     *
     * If this function succeeds, values added with addValue() since the
     * previous clearSet() are added to the storage for the option.
     * Only throws in out-of-memory conditions, and provides the strong
     * exception safety guarantee as long as the copy constructor of `T`
     * does not throw.
     *
     * See addValue() for cases where this method should be used in derived
     * classes.
     *
     * Calls clearSet() if it is successful.
     */
    void commitValues();

    /*! \brief
     * Sets the default value for the option.
     *
     * \param[in] value  Default value to set.
     * \throws    std::bad_alloc if out of memory.
     *
     * This method can be used from the derived class constructor to
     * programmatically set a default value.
     */
    void setDefaultValue(const T& value);
    /*! \brief
     * Sets the default value if set for the option.
     *
     * \param[in] value  Default value to set.
     * \throws    std::bad_alloc if out of memory.
     *
     * This method can be used from the derived class constructor to
     * programmatically set a default value.
     */
    void setDefaultValueIfSet(const T& value);

    /*! \brief
     * Provides derived classes access to the current list of values.
     *
     * The non-const any should only be used from processAll() in
     * derived classes if necessary.
     */
    ArrayRef<T> values() { return store_->values(); }
    //! Provides derived classes access to the current list of values.
    ArrayRef<const T> values() const { return store_->values(); }

private:
    //! Creates the internal storage object for final values..
    StorePointer createStore(ValueList* storeVector, T* store, int* storeCount, int initialCount);

    /*! \brief
     * Vector for temporary storage of values before commitSet() is called.
     */
    ValueList setValues_;
    //! Final storage for option values.
    const StorePointer store_;
    // This never releases ownership.
    std::unique_ptr<T> defaultValueIfSet_;

    // Copy and assign disallowed by base.
};


/*! \libinternal \brief
 * Simplified option storage template for options that have one-to-one value
 * conversion.
 *
 * \tparam T Assignable type that stores a single option value.
 *
 * To implement an option that always map a single input value to a single
 * output value, derive from this class instead of OptionStorageTemplate.
 * This class implements convertValue() in terms of two new virtual methods:
 * initConverter() and processValue().
 *
 * To specify how different types of values need to be converted, implement
 * initConverter().
 * To do common post-processing of the values after conversion, but before they
 * are added to the underlying storage, override processValue().
 *
 * \inlibraryapi
 * \ingroup module_options
 */
template<typename T>
class OptionStorageTemplateSimple : public OptionStorageTemplate<T>
{
public:
    //! Alias for the template class for use in base classes.
    typedef OptionStorageTemplateSimple<T> MyBase;

protected:
    //! Alias for the converter parameter type for initConverter().
    typedef OptionValueConverterSimple<T> ConverterType;

    //! Initializes the storage.
    template<class U>
    explicit OptionStorageTemplateSimple(const OptionTemplate<T, U>& settings,
                                         OptionFlags                 staticFlags = OptionFlags()) :
        OptionStorageTemplate<T>(settings, staticFlags), initialized_(false)
    {
    }
    //! Initializes the storage.
    OptionStorageTemplateSimple(const AbstractOption&                           settings,
                                typename OptionStorageTemplate<T>::StorePointer store) :
        OptionStorageTemplate<T>(settings, std::move(store)), initialized_(false)
    {
    }

    std::vector<Any> normalizeValues(const std::vector<Any>& values) const override
    {
        const_cast<MyBase*>(this)->ensureConverterInitialized();
        std::vector<Any> result;
        result.reserve(values.size());
        for (const auto& value : values)
        {
            result.push_back(normalizeValue(converter_.convert(value)));
        }
        return result;
    }

    /*! \brief
     * Specifies how different types are converted.
     *
     * See OptionValueConverterSimple for more details.
     */
    virtual void initConverter(ConverterType* converter) = 0;
    /*! \brief
     * Post-processes a value after conversion to the output type.
     *
     * \param[in] value  Value after conversion.
     * \returns   Value to store for the option.
     *
     * The default implementation only provides an identity mapping.
     */
    virtual T processValue(const T& value) const { return value; }
    /*! \brief
     * Converts a single value to normalized type.
     *
     * \param[in] value  Value after conversion.
     * \returns   Value to store for the option.
     *
     * This can be overridden to serialize a different type than `T`
     * when using the option with KeyValueTreeObject.
     */
    virtual Any normalizeValue(const T& value) const { return Any::create<T>(processValue(value)); }

private:
    void convertValue(const Any& any) override
    {
        ensureConverterInitialized();
        this->addValue(processValue(converter_.convert(any)));
    }
    void ensureConverterInitialized()
    {
        if (!initialized_)
        {
            initConverter(&converter_);
            initialized_ = true;
        }
    }

    ConverterType converter_;
    bool          initialized_;
};


/********************************************************************
 * OptionStorageTemplate implementation
 */

template<typename T>
template<class U>
OptionStorageTemplate<T>::OptionStorageTemplate(const OptionTemplate<T, U>& settings,
                                                OptionFlags                 staticFlags) :
    AbstractOptionStorage(settings, staticFlags),
    store_(createStore(settings.storeVector_,
                       settings.store_,
                       settings.countptr_,
                       (settings.isVector() ? settings.maxValueCount_ : settings.minValueCount_)))
{
    if (hasFlag(efOption_NoDefaultValue)
        && (settings.defaultValue_ != nullptr || settings.defaultValueIfSet_ != nullptr))
    {
        GMX_THROW(APIError("Option does not support default value, but one is set"));
    }
    if (!hasFlag(efOption_NoDefaultValue))
    {
        setFlag(efOption_HasDefaultValue);
        if (settings.defaultValue_ != nullptr)
        {
            setDefaultValue(*settings.defaultValue_);
        }
        if (settings.defaultValueIfSet_ != nullptr)
        {
            setDefaultValueIfSet(*settings.defaultValueIfSet_);
        }
    }
}


template<typename T>
OptionStorageTemplate<T>::OptionStorageTemplate(const AbstractOption& settings, StorePointer store) :
    AbstractOptionStorage(settings, OptionFlags()), store_(std::move(store))
{
}


template<typename T>
std::unique_ptr<IOptionValueStore<T>> OptionStorageTemplate<T>::createStore(ValueList* storeVector,
                                                                            T*         store,
                                                                            int* storeCount, // NOLINT(readability-non-const-parameter) passed non-const to OptionValueStorePlain
                                                                            int initialCount)
{
    if (storeVector != nullptr)
    {
        GMX_RELEASE_ASSERT(store == nullptr && storeCount == nullptr,
                           "Cannot specify more than one storage location");
        return StorePointer(new OptionValueStoreVector<T>(storeVector));
    }
    else if (store != nullptr)
    {
        // If the maximum number of values is not known, storage to
        // caller-allocated memory is unsafe.
        if (maxValueCount() < 0 || hasFlag(efOption_MultipleTimes))
        {
            GMX_THROW(APIError("Cannot set user-allocated storage for arbitrary number of values"));
        }
        if (storeCount == nullptr && !isVector() && minValueCount() != maxValueCount())
        {
            GMX_THROW(
                    APIError("Count storage is not set, although the number of produced values is "
                             "not known"));
        }
        if (hasFlag(efOption_NoDefaultValue))
        {
            initialCount = 0;
        }
        return StorePointer(new OptionValueStorePlain<T>(store, storeCount, initialCount));
    }
    GMX_RELEASE_ASSERT(storeCount == nullptr, "Cannot specify count storage without value storage");
    return StorePointer(new OptionValueStoreNull<T>());
}


template<typename T>
std::vector<Any> OptionStorageTemplate<T>::defaultValues() const
{
    std::vector<Any> result;
    if (hasFlag(efOption_NoDefaultValue))
    {
        return result;
    }
    GMX_RELEASE_ASSERT(
            hasFlag(efOption_HasDefaultValue),
            "Current option implementation can only provide default values before assignment");
    for (const auto& value : values())
    {
        result.push_back(Any::create<T>(value));
    }
    return normalizeValues(result);
}


template<typename T>
std::vector<std::string> OptionStorageTemplate<T>::defaultValuesAsStrings() const
{
    std::vector<std::string> result;
    if (hasFlag(efOption_NoDefaultValue))
    {
        return result;
    }
    GMX_RELEASE_ASSERT(
            hasFlag(efOption_HasDefaultValue),
            "Current option implementation can only provide default values before assignment");
    for (const auto& value : values())
    {
        result.push_back(formatSingleValue(value));
    }
    if (result.empty() || (result.size() == 1 && result[0].empty()))
    {
        result.clear();
        if (defaultValueIfSet_ != nullptr)
        {
            result.push_back(formatSingleValue(*defaultValueIfSet_));
        }
    }
    return result;
}


template<typename T>
void OptionStorageTemplate<T>::clearSet()
{
    setValues_.clear();
}


template<typename T>
void OptionStorageTemplate<T>::processSet()
{
    processSetValues(&setValues_);
    if (setValues_.empty() && defaultValueIfSet_ != nullptr)
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


template<typename T>
void OptionStorageTemplate<T>::addValue(const T& value)
{
    if (maxValueCount() >= 0 && setValues_.size() >= static_cast<size_t>(maxValueCount()))
    {
        GMX_THROW(InvalidInputError("Too many values"));
    }
    setValues_.push_back(value);
}

template<typename T>
void OptionStorageTemplate<T>::commitValues()
{
    if (hasFlag(efOption_ClearOnNextSet))
    {
        store_->clear();
    }
    store_->reserve(setValues_.size());
    // For bool the loop variable isn't a reference (it's its special reference type)
    CLANG_DIAGNOSTIC_IGNORE("-Wrange-loop-analysis");
    for (const auto& value : setValues_)
    {
        store_->append(value);
    }
    CLANG_DIAGNOSTIC_RESET;
    clearSet();
}

template<typename T>
void OptionStorageTemplate<T>::setDefaultValue(const T& value)
{
    if (hasFlag(efOption_NoDefaultValue))
    {
        GMX_THROW(APIError("Option does not support default value, but one is set"));
    }
    if (hasFlag(efOption_HasDefaultValue))
    {
        setFlag(efOption_ExplicitDefaultValue);
        store_->clear();
        store_->append(value);
    }
}


template<typename T>
void OptionStorageTemplate<T>::setDefaultValueIfSet(const T& value)
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
    defaultValueIfSet_ = std::make_unique<T>(value);
}

} // namespace gmx

#endif
