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
/*! \file
 * \brief
 * Declares option objects for basic option types.
 *
 * Together with options.h, this header forms the part of the public API
 * that most classes will use to provide options.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_BASICOPTIONS_H
#define GMX_OPTIONS_BASICOPTIONS_H

#include <memory>
#include <string>
#include <vector>

#include "gromacs/options/abstractoption.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/gmxassert.h"

#include "ivaluestore.h"

namespace gmx
{

class BooleanOptionInfo;
class BooleanOptionStorage;
class IntegerOptionInfo;
class IntegerOptionStorage;
class Int64OptionInfo;
class Int64OptionStorage;
class DoubleOptionInfo;
class DoubleOptionStorage;
class FloatOptionInfo;
class FloatOptionStorage;
class StringOptionInfo;
class StringOptionStorage;
class EnumOptionInfo;
class EnumOptionStorage;

//! \addtogroup module_options
//! \{

/*! \brief
 * Specifies an option that provides boolean values.
 *
 * Example:
 * \code
   bool  bPBC;
   using gmx::BooleanOption;
   options.addOption(BooleanOption("pbc").store(&bPBC));
 * \endcode
 *
 * Public methods in this class do not throw.
 *
 * \inpublicapi
 */
class BooleanOption : public OptionTemplate<bool, BooleanOption>
{
public:
    //! OptionInfo subclass corresponding to this option type.
    typedef BooleanOptionInfo InfoType;

    //! Initializes an option with the given name.
    explicit BooleanOption(const char* name) : MyBase(name) {}

private:
    //! Creates a BooleanOptionStorage object.
    AbstractOptionStorage* createStorage(const OptionManagerContainer& managers) const override;
};

/*! \brief
 * Specifies an option that provides integer values.
 *
 * Examples:
 * \code
   using gmx::IntegerOption;
   // Simple option
   int  rcut = 0;
   options.addOption(IntegerOption("rcut").store(&rcut));
   // Vector-valued option
   int  box[3] = {1, 1, 1};  // Default value
   options.addOption(IntegerOption("box").store(box).vector());
 * \endcode
 *
 * Public methods in this class do not throw.
 *
 * \inpublicapi
 */
class IntegerOption : public OptionTemplate<int, IntegerOption>
{
public:
    //! OptionInfo subclass corresponding to this option type.
    typedef IntegerOptionInfo InfoType;

    //! Initializes an option with the given name.
    explicit IntegerOption(const char* name) : MyBase(name) {}

    /*! \brief
     * Sets the option to return a vector value.
     *
     * A vector value returns a fixed number of values, the default being
     * three (can be changed with valueCount()).  However, it also accepts
     * a single value, in which case the value is used to fill the whole
     * vector.
     */
    MyClass& vector()
    {
        setVector();
        return me();
    }

private:
    //! Creates an IntegerOptionStorage object.
    AbstractOptionStorage* createStorage(const OptionManagerContainer& managers) const override;

    /*! \brief
     * Needed to initialize IntegerOptionStorage from this class without
     * otherwise unnecessary accessors.
     */
    friend class IntegerOptionStorage;
};

/*! \brief
 * Specifies an option that provides 64-bit integer values.
 *
 * Public methods in this class do not throw.
 *
 * \see IntegerOption
 *
 * \inpublicapi
 */
class Int64Option : public OptionTemplate<int64_t, Int64Option>
{
public:
    //! OptionInfo subclass corresponding to this option type.
    typedef Int64OptionInfo InfoType;

    //! Initializes an option with the given name.
    explicit Int64Option(const char* name) : MyBase(name) {}

private:
    //! Creates an Int64OptionStorage object.
    AbstractOptionStorage* createStorage(const OptionManagerContainer& managers) const override;

    /*! \brief
     * Needed to initialize Int64OptionStorage from this class without
     * otherwise unnecessary accessors.
     */
    friend class Int64OptionStorage;
};

/*! \brief
 * Specifies an option that provides floating-point (double) values.
 *
 * Public methods in this class do not throw.
 *
 * \inpublicapi
 */
class DoubleOption : public OptionTemplate<double, DoubleOption>
{
public:
    //! OptionInfo subclass corresponding to this option type.
    typedef DoubleOptionInfo InfoType;

    //! Initializes an option with the given name.
    explicit DoubleOption(const char* name) : MyBase(name), bTime_(false) {}

    //! \copydoc IntegerOption::vector()
    MyClass& vector()
    {
        setVector();
        return me();
    }
    /*! \brief
     * Marks this option as providing a time value whose unit can be changed.
     *
     * By itself, this option does nothing.  It marks the option as a time
     * value such that TimeUnitManager::scaleTimeOptions() can process it.
     * In typical cases, \Gromacs scales the time options just before
     * Options::finish() has been called, so the option value is only
     * available after all option values have been processed.
     * All values in the program are in ps (including any default value);
     * user-provided values are scaled according to the time unit set in
     * TimeUnitManager.
     */
    MyClass& timeValue()
    {
        bTime_ = true;
        return me();
    }

private:
    //! Creates a DoubleOptionStorage object.
    AbstractOptionStorage* createStorage(const OptionManagerContainer& managers) const override;

    bool bTime_;

    /*! \brief
     * Needed to initialize DoubleOptionStorage from this class without
     * otherwise unnecessary accessors.
     */
    friend class DoubleOptionStorage;
};

/*! \brief
 * Specifies an option that provides floating-point (float) values.
 *
 * Public methods in this class do not throw.
 *
 * \see DoubleOption
 *
 * \inpublicapi
 */
class FloatOption : public OptionTemplate<float, FloatOption>
{
public:
    //! OptionInfo subclass corresponding to this option type.
    typedef FloatOptionInfo InfoType;

    //! Initializes an option with the given name.
    explicit FloatOption(const char* name) : MyBase(name), bTime_(false) {}

    //! \copydoc IntegerOption::vector()
    MyClass& vector()
    {
        setVector();
        return me();
    }
    //! \copydoc DoubleOption::timeValue()
    MyClass& timeValue()
    {
        bTime_ = true;
        return me();
    }

private:
    //! Creates a FloatOptionStorage object.
    AbstractOptionStorage* createStorage(const OptionManagerContainer& managers) const override;

    bool bTime_;

    /*! \brief
     * Needed to initialize FloatOptionStorage from this class without
     * otherwise unnecessary accessors.
     */
    friend class FloatOptionStorage;
};

/*! \brief
 * Specifies an option that provides string values.
 *
 * Examples:
 * \code
   using gmx::StringOption;
   // Simple option
   std::string  str;
   options.addOption(StringOption("str").store(&str));
   // Option that only accepts predefined values
   const char * const  allowed[] = { "atom", "residue", "molecule" };
   std::string  str;
   options.addOption(StringOption("type").enumValue(allowed).store(&str));
 * \endcode
 *
 * Public methods in this class do not throw.
 *
 * \inpublicapi
 */
class StringOption : public OptionTemplate<std::string, StringOption>
{
public:
    //! OptionInfo subclass corresponding to this option type.
    typedef StringOptionInfo InfoType;

    //! Initializes an option with the given name.
    explicit StringOption(const char* name) :
        MyBase(name), enumValues_(nullptr), enumValuesCount_(0), defaultEnumIndex_(-1)
    {
    }

    /*! \brief
     * Sets the option to only accept one of a fixed set of strings.
     *
     * \param[in] values  Array of strings to accept.
     *
     * Also accepts prefixes of the strings; if a prefix matches more than
     * one of the possible strings, the shortest one is used (in a tie, the
     * first one is).
     *
     * The strings are copied once the option is created.
     */
    template<size_t count>
    MyClass& enumValue(const char* const (&values)[count])
    {
        GMX_ASSERT(enumValues_ == nullptr, "Multiple sets of enumerated values specified");
        enumValues_      = values;
        enumValuesCount_ = count;
        return me();
    }
    /*! \brief
     * Sets the option to only accept one of a fixed set of strings.
     *
     * \param[in] values  Array of strings to accept, with a NULL pointer
     *      following the last string.
     *
     * Works otherwise as the array version, but accepts a pointer to
     * an array of undetermined length.  The end of the array is indicated
     * by a NULL pointer in the array.
     *
     * \see enumValue()
     */
    MyClass& enumValueFromNullTerminatedArray(const char* const* values)
    {
        GMX_ASSERT(enumValues_ == nullptr, "Multiple sets of enumerated values specified");
        enumValues_      = values;
        enumValuesCount_ = -1;
        return me();
    }
    /*! \brief
     * Sets the default value using an index into the enumeration table.
     *
     * Cannot be specified without enumValue().
     */
    MyClass& defaultEnumIndex(int index)
    {
        GMX_ASSERT(index >= 0, "Invalid enumeration index");
        defaultEnumIndex_ = index;
        return me();
    }

private:
    //! Creates a StringOptionStorage object.
    AbstractOptionStorage* createStorage(const OptionManagerContainer& managers) const override;

    const char* const* enumValues_;
    int                enumValuesCount_;
    int                defaultEnumIndex_;

    /*! \brief
     * Needed to initialize StringOptionStorage from this class without
     * otherwise unnecessary accessors.
     */
    friend class StringOptionStorage;
};

//! \}

namespace internal
{

/*! \internal
 * \brief
 * Type-specific implementation for IOptionValueStore for an enum option.
 *
 * This class is instantiated for each enum type for which EnumOption is used,
 * and takes care of managing `int`-to-`enum` conversions.  Having this part in
 * the header allows the actual storage implementation to not be in the header,
 * which would require exposing all the internals through this one header...
 *
 * \ingroup module_options
 */
template<typename EnumType>
class EnumIndexStore : public IOptionValueStore<int>
{
public:
    //! Initializes the storage for the given actual enum variables.
    EnumIndexStore(EnumType* store, std::vector<EnumType>* storeVector) :
        store_(store), storeVector_(storeVector)
    {
        if (storeVector_ != nullptr)
        {
            for (EnumType value : *storeVector_)
            {
                intStore_.push_back(static_cast<int>(value));
            }
        }
        else if (store_ != nullptr)
        {
            // TODO: Copy more than one value if that would make sense.
            intStore_.push_back(static_cast<int>(store_[0]));
        }
    }

    int           valueCount() override { return gmx::ssize(intStore_); }
    ArrayRef<int> values() override { return intStore_; }
    void          clear() override
    {
        intStore_.clear();
        if (storeVector_ != nullptr)
        {
            storeVector_->clear();
        }
    }
    void reserve(size_t count) override
    {
        intStore_.reserve(intStore_.size() + count);
        if (storeVector_ != nullptr)
        {
            storeVector_->reserve(storeVector_->size() + count);
        }
    }
    void append(const int& value) override
    {
        const size_t count = intStore_.size();
        intStore_.push_back(value);
        if (store_ != nullptr)
        {
            store_[count] = static_cast<EnumType>(value);
        }
        if (storeVector_ != nullptr)
        {
            storeVector_->push_back(static_cast<EnumType>(value));
        }
    }

private:
    //! Stores the integer values for values().
    std::vector<int>       intStore_;
    EnumType*              store_;
    std::vector<EnumType>* storeVector_;
};

//! \cond internal
/*! \internal
 * \brief
 * Helper to create EnumOptionStorage instances.
 *
 * This function works as a proxy between EnumOption::createStorage() and the
 * EnumOptionStorage constructor, such that the latter does not need to be
 * exposed in the header.
 *
 * \ingroup module_options
 */
AbstractOptionStorage* createEnumOptionStorage(const AbstractOption& option,
                                               const char* const*    enumValues,
                                               int                   count,
                                               int                   defaultValue,
                                               int                   defaultValueIfSet,
                                               std::unique_ptr<IOptionValueStore<int>> store);
//! \endcond

} // namespace internal

//! \addtogroup module_options
//! \{

/*! \brief
* Specifies an option that accepts an EnumerationArray of string values and writes the
* selected index into an `enum` variable.
*
* \tparam EnumType  DataType of the variable that receives the values
*
* Examples:
* \code
  enum class MyEnum { Atom, Res, Mol, Count } : int;
  EnumerationArray<MyEnum, const char *> myEnumNames = { "atom", "residue", "molecule" };
  MyEnum       value = MyEnum::Atom; // default value
  options.addOption(EnumOption<MyEnum>("type").enumValue(myEnumNames).store(&value));
* \endcode
*
* storeCount() is not currently implemented for this option type, and
* providing multiple default values through an array passed to store() does
* not work consistently in all cases.
* In the current implementation, the values of the enum type should correspond
* to indices in the array passed to enumValue(), i.e., be consecutive
* starting from zero.  Only values corresponding to valid indices are accepted
* as parameters to, e.g., defaultValue().  However, other values can be used
* as the initial value of the variable (`value` in the above example), and
* those will be preserved if the option is not set.
*
* Public methods in this class do not throw.
*
* \inpublicapi
*/
template<typename EnumType>
class EnumOption : public OptionTemplate<EnumType, EnumOption<EnumType>>
{
public:
    //! OptionInfo subclass corresponding to this option type.
    typedef EnumOptionInfo InfoType;

    // This needs to be duplicated from OptionTemplate because this class
    // is a template.
    //! Short-hand for the base class.
    typedef OptionTemplate<EnumType, EnumOption<EnumType>> MyBase;

    //! Initializes an option with the given name.
    explicit EnumOption(const char* name) : MyBase(name), enumValues_(nullptr), enumValuesCount_(0)
    {
    }

    /*! \brief
     * Sets the option to only accept one of a fixed set of strings.
     *

     * \param[in] values  Array of strings to accept.
     *
     * Also accepts prefixes of the strings; if a prefix matches more than
     * one of the possible strings, the shortest one is used (in a tie, the
     * first one is).
     *
     * The strings are copied once the option is created.
     */
    EnumOption& enumValue(const EnumerationArray<EnumType, const char*>& values)
    {
        GMX_ASSERT(enumValues_ == nullptr, "Multiple sets of enumerated values specified");
        enumValues_      = values.data();
        enumValuesCount_ = values.size();
        return MyBase::me();
    }

private:
    //! Helper function to convert default values for storage initialization.
    static int convertToInt(const EnumType* defaultValue)
    {
        return defaultValue != nullptr ? static_cast<int>(*defaultValue) : -1;
    }

    //! Creates a EnumOptionStorage object.
    AbstractOptionStorage* createStorage(const OptionManagerContainer& /*managers*/) const override
    {
        // TODO: Implement storeCount() if necessary.
        return internal::createEnumOptionStorage(*this,
                                                 enumValues_,
                                                 enumValuesCount_,
                                                 convertToInt(MyBase::defaultValue()),
                                                 convertToInt(MyBase::defaultValueIfSet()),
                                                 std::make_unique<internal::EnumIndexStore<EnumType>>(
                                                         MyBase::store(), MyBase::storeVector()));
    }

    const char* const* enumValues_;
    int                enumValuesCount_;

    /*! \brief
     * Needed to initialize EnumOptionStorage from this class without
     * otherwise unnecessary accessors.
     */
    friend class EnumOptionStorage;
};

/*! \brief
* Specifies an option that accepts enumerated string values and writes the
* selected index into an `enum` variable.
*
* \tparam EnumType  Type of the variable that receives the values
*     (can also be `int`).
*
* Examples:
* \code
  enum MyEnum { eAtom, eRes, eMol };
  using gmx::LegacyEnumOption;
  const char * const  allowed[] = { "atom", "residue", "molecule" };
  MyEnum       value = eAtom; // default value
  options.addOption(LegacyEnumOption<MyEnum>("type").enumValue(allowed).store(&value));
* \endcode
*
* Works exactly as EnumOption.
*
* This is legacy support for pargsToOptions and can be removed when it
* is removed.  No new uses of it should be made.
*
* Public methods in this class do not throw.
*
* \inpublicapi
*/
template<typename EnumType>
class LegacyEnumOption : public OptionTemplate<EnumType, LegacyEnumOption<EnumType>>
{
public:
    //! OptionInfo subclass corresponding to this option type.
    typedef EnumOptionInfo InfoType;

    // This needs to be duplicated from OptionTemplate because this class
    // is a template.
    //! Short-hand for the base class.
    typedef OptionTemplate<EnumType, LegacyEnumOption<EnumType>> MyBase;

    //! Initializes an option with the given name.
    explicit LegacyEnumOption(const char* name) :
        MyBase(name), enumValues_(nullptr), enumValuesCount_(0)
    {
    }

    /*! \brief
     * Sets the option to only accept one of a fixed set of strings.
     *
     * \param[in] values  Array of strings to accept.
     *
     * Also accepts prefixes of the strings; if a prefix matches more than
     * one of the possible strings, the shortest one is used (in a tie, the
     * first one is).
     *
     * The strings are copied once the option is created.
     */
    template<size_t count>
    LegacyEnumOption& enumValue(const char* const (&values)[count])
    {
        GMX_ASSERT(enumValues_ == nullptr, "Multiple sets of enumerated values specified");
        enumValues_      = values;
        enumValuesCount_ = count;
        return MyBase::me();
    }
    /*! \brief
     * Sets the option to only accept one of a fixed set of strings.
     *
     * \param[in] values  Array of strings to accept, with a NULL pointer
     *      following the last string.
     *
     * Works otherwise as the array version, but accepts a pointer to
     * an array of undetermined length.  The end of the array is indicated
     * by a NULL pointer in the array.
     *
     * \see enumValue()
     */
    LegacyEnumOption& enumValueFromNullTerminatedArray(const char* const* values)
    {
        GMX_ASSERT(enumValues_ == nullptr, "Multiple sets of enumerated values specified");
        enumValues_      = values;
        enumValuesCount_ = -1;
        return MyBase::me();
    }

private:
    //! Helper function to convert default values for storage initialization.
    static int convertToInt(const EnumType* defaultValue)
    {
        return defaultValue != nullptr ? static_cast<int>(*defaultValue) : -1;
    }

    //! Creates a EnumOptionStorage object.
    AbstractOptionStorage* createStorage(const OptionManagerContainer& /*managers*/) const override
    {
        // TODO: Implement storeCount() if necessary.
        return internal::createEnumOptionStorage(*this,
                                                 enumValues_,
                                                 enumValuesCount_,
                                                 convertToInt(MyBase::defaultValue()),
                                                 convertToInt(MyBase::defaultValueIfSet()),
                                                 std::make_unique<internal::EnumIndexStore<EnumType>>(
                                                         MyBase::store(), MyBase::storeVector()));
    }

    const char* const* enumValues_;
    int                enumValuesCount_;

    /*! \brief
     * Needed to initialize EnumOptionStorage from this class without
     * otherwise unnecessary accessors.
     */
    friend class EnumOptionStorage;
};

/*! \brief
 * Wrapper class for accessing boolean option information.
 *
 * \inpublicapi
 */
class BooleanOptionInfo : public OptionInfo
{
public:
    //! Creates an option info object for the given option.
    explicit BooleanOptionInfo(BooleanOptionStorage* option);

    //! Returns the default value for this option.
    bool defaultValue() const;

private:
    const BooleanOptionStorage& option() const;
};

/*! \brief
 * Wrapper class for accessing integer option information.
 *
 * \inpublicapi
 */
class IntegerOptionInfo : public OptionInfo
{
public:
    //! Creates an option info object for the given option.
    explicit IntegerOptionInfo(IntegerOptionStorage* option);
};

/*! \brief
 * Wrapper class for accessing 64-bit integer option information.
 *
 * \inpublicapi
 */
class Int64OptionInfo : public OptionInfo
{
public:
    //! Creates an option info object for the given option.
    explicit Int64OptionInfo(Int64OptionStorage* option);
};

/*! \brief
 * Wrapper class for accessing floating-point option information.
 *
 * \inpublicapi
 */
class DoubleOptionInfo : public OptionInfo
{
public:
    //! Creates an option info object for the given option.
    explicit DoubleOptionInfo(DoubleOptionStorage* option);

    //! Whether the option specifies a time value.
    bool isTime() const;

    /*! \brief
     * Sets a scale factor for user-provided values.
     *
     * Any user-provided value is scaled by the provided factor.
     * Programmatically set default values are not scaled.
     * If called multiple times, later calls override the previously set
     * value.  In other words, the scaling is not cumulative.
     */
    void setScaleFactor(double factor);

private:
    DoubleOptionStorage&       option();
    const DoubleOptionStorage& option() const;
};

/*! \brief
 * Wrapper class for accessing floating-point option information.
 *
 * \inpublicapi
 */
class FloatOptionInfo : public OptionInfo
{
public:
    //! Creates an option info object for the given option.
    explicit FloatOptionInfo(FloatOptionStorage* option);

    //! Whether the option specifies a time value.
    bool isTime() const;

    //! \copydoc DoubleOptionInfo::setScaleFactor()
    void setScaleFactor(double factor);

private:
    FloatOptionStorage&       option();
    const FloatOptionStorage& option() const;
};

/*! \brief
 * Wrapper class for accessing string option information.
 *
 * \inpublicapi
 */
class StringOptionInfo : public OptionInfo
{
public:
    //! Creates an option info object for the given option.
    explicit StringOptionInfo(StringOptionStorage* option);

    /*! \brief
     * Whether this option accepts an enumerated set of values.
     *
     * Returns true if StringOption::enumValues() was used when creating
     * this option.
     */
    bool isEnumerated() const;
    /*! \brief
     * Returns the set of allowed values for this option.
     *
     * Returns an empty vector if isEnumerated() returns false.
     */
    const std::vector<std::string>& allowedValues() const;

private:
    const StringOptionStorage& option() const;
};

/*! \brief
 * Wrapper class for accessing enum option information.
 *
 * \inpublicapi
 */
class EnumOptionInfo : public OptionInfo
{
public:
    //! Creates an option info object for the given option.
    explicit EnumOptionInfo(EnumOptionStorage* option);

    /*! \brief
     * Returns the set of allowed values for this option.
     */
    const std::vector<std::string>& allowedValues() const;

private:
    const EnumOptionStorage& option() const;
};

/*! \typedef RealOption
 * \brief
 * Typedef for either DoubleOption or FloatOption, depending on precision.
 *
 * Generally, new would be better using DoubleOption, but this is provided for
 * cases where the output value needs to be of type `real` for some reason.
 */
/*! \typedef RealOptionInfo
 * \brief
 * Typedef for either DoubleOptionInfo or FloatOptionInfo, depending on precision.
 *
 * Generally, new would be better using DoubleOption, but this is provided for
 * cases where the output value needs to be of type `real` for some reason.
 */
#if GMX_DOUBLE
typedef DoubleOption     RealOption;
typedef DoubleOptionInfo RealOptionInfo;
#else
typedef FloatOption     RealOption;
typedef FloatOptionInfo RealOptionInfo;
#endif

//! \}

} // namespace gmx

#endif
