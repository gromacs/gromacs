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

#include <string>

#include "gromacs/options/abstractoption.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

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

/*! \addtogroup module_options
 * \{
 */

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
        explicit BooleanOption(const char *name) : MyBase(name) {}

    private:
        //! Creates a BooleanOptionStorage object.
        virtual AbstractOptionStorage *createStorage(
            const OptionManagerContainer &managers) const;
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
        explicit IntegerOption(const char *name) : MyBase(name) {}

        /*! \brief
         * Sets the option to return a vector value.
         *
         * A vector value returns a fixed number of values, the default being
         * three (can be changed with valueCount()).  However, it also accepts
         * a single value, in which case the value is used to fill the whole
         * vector.
         */
        MyClass &vector() { setVector(); return me(); }

    private:
        //! Creates an IntegerOptionStorage object.
        virtual AbstractOptionStorage *createStorage(
            const OptionManagerContainer &managers) const;

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
class Int64Option : public OptionTemplate<gmx_int64_t, Int64Option>
{
    public:
        //! OptionInfo subclass corresponding to this option type.
        typedef Int64OptionInfo InfoType;

        //! Initializes an option with the given name.
        explicit Int64Option(const char *name) : MyBase(name) {}

    private:
        //! Creates an Int64OptionStorage object.
        virtual AbstractOptionStorage *createStorage(
            const OptionManagerContainer &managers) const;

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
        explicit DoubleOption(const char *name) : MyBase(name), bTime_(false)
        {
        }

        //! \copydoc IntegerOption::vector()
        MyClass &vector() { setVector(); return me(); }
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
        MyClass &timeValue() { bTime_ = true; return me(); }

    private:
        //! Creates a DoubleOptionStorage object.
        virtual AbstractOptionStorage *createStorage(
            const OptionManagerContainer &managers) const;

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
        explicit FloatOption(const char *name) : MyBase(name), bTime_(false)
        {
        }

        //! \copydoc IntegerOption::vector()
        MyClass &vector() { setVector(); return me(); }
        //! \copydoc DoubleOption::timeValue()
        MyClass &timeValue() { bTime_ = true; return me(); }

    private:
        //! Creates a FloatOptionStorage object.
        virtual AbstractOptionStorage *createStorage(
            const OptionManagerContainer &managers) const;

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
   int          type;
   options.addOption(StringOption("type").enumValue(allowed).store(&str)
                        .storeEnumIndex(&type));
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
        explicit StringOption(const char *name)
            : MyBase(name), enumValues_(NULL), enumValuesCount_(0),
              defaultEnumIndex_(-1), enumIndexStore_(NULL)
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
         * It is not possible to provide multiple values for an option with
         * this property set, i.e., valueCount() and similar attributes cannot
         * be set.
         *
         * The strings are copied once the option is created.
         */
        template <size_t count>
        MyClass &enumValue(const char *const (&values)[count])
        {
            GMX_ASSERT(enumValues_ == NULL,
                       "Multiple sets of enumerated values specified");
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
        MyClass &enumValueFromNullTerminatedArray(const char *const *values)
        {
            GMX_ASSERT(enumValues_ == NULL,
                       "Multiple sets of enumerated values specified");
            enumValues_      = values;
            enumValuesCount_ = -1;
            return me();
        }
        /*! \brief
         * Sets the default value using an index into the enumeration table.
         *
         * Cannot be specified without enumValue().
         */
        MyClass &defaultEnumIndex(int index)
        {
            GMX_ASSERT(index >= 0, "Invalid enumeration index");
            defaultEnumIndex_ = index;
            return me();
        }
        /*! \brief
         * Stores the index of the selected value into the provided memory
         * location.
         *
         * The index (zero-based) of the selected value in the array \p values
         * provided to enumValues() is written into \p *store after the
         * option gets its value.  If the option has not been provided,
         * and there is no default value, -1 is stored.  If store(),
         * storeVector() or defaultEnumIndex() is not present, the value in
         * \p *store is kept as a default value, otherwise it is always
         * overwritten.
         *
         * Cannot be specified without enumValue().
         *
         * \todo
         * Implement this such that it is also possible to store the value
         * directly into a real enum type.
         */
        MyClass &storeEnumIndex(int *store)
        { enumIndexStore_ = store; return me(); }

    private:
        //! Creates a StringOptionStorage object.
        virtual AbstractOptionStorage *createStorage(
            const OptionManagerContainer &managers) const;

        const char *const      *enumValues_;
        int                     enumValuesCount_;
        int                     defaultEnumIndex_;
        int                    *enumIndexStore_;

        /*! \brief
         * Needed to initialize StringOptionStorage from this class without
         * otherwise unnecessary accessors.
         */
        friend class StringOptionStorage;
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
        explicit BooleanOptionInfo(BooleanOptionStorage *option);

        //! Returns the default value for this option.
        bool defaultValue() const;

    private:
        const BooleanOptionStorage &option() const;
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
        explicit IntegerOptionInfo(IntegerOptionStorage *option);
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
        explicit Int64OptionInfo(Int64OptionStorage *option);
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
        explicit DoubleOptionInfo(DoubleOptionStorage *option);

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
        DoubleOptionStorage &option();
        const DoubleOptionStorage &option() const;
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
        explicit FloatOptionInfo(FloatOptionStorage *option);

        //! Whether the option specifies a time value.
        bool isTime() const;

        //! \copydoc DoubleOptionInfo::setScaleFactor()
        void setScaleFactor(double factor);

    private:
        FloatOptionStorage &option();
        const FloatOptionStorage &option() const;
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
        explicit StringOptionInfo(StringOptionStorage *option);

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
        const std::vector<std::string> &allowedValues() const;

    private:
        StringOptionStorage &option();
        const StringOptionStorage &option() const;
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
#ifdef GMX_DOUBLE
typedef DoubleOption     RealOption;
typedef DoubleOptionInfo RealOptionInfo;
#else
typedef FloatOption      RealOption;
typedef FloatOptionInfo  RealOptionInfo;
#endif

/*!\}*/

} // namespace gmx

#endif
