/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Declares storage classes for basic option types.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_BASICOPTIONSTORAGE_H
#define GMX_OPTIONS_BASICOPTIONSTORAGE_H

#include <memory>
#include <string>
#include <vector>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/optionstoragetemplate.h"

namespace gmx
{

/*! \addtogroup module_options
 * \{
 */

/*! \internal \brief
 * Converts, validates, and stores boolean values.
 */
class BooleanOptionStorage : public OptionStorageTemplateSimple<bool>
{
    public:
        /*! \brief
         * Initializes the storage from option settings.
         *
         * \param[in] settings   Storage settings.
         */
        explicit BooleanOptionStorage(const BooleanOption &settings)
            : MyBase(settings), info_(this)
        {
        }

        virtual OptionInfo &optionInfo() { return info_; }
        virtual std::string typeString() const { return "bool"; }
        virtual std::string formatSingleValue(const bool &value) const;

        //! \copydoc BooleanOptionInfo::defaultValue()
        bool defaultValue() const { return valueCount() > 0 && values()[0]; }

    private:
        virtual void initConverter(ConverterType *converter);

        BooleanOptionInfo       info_;
};

/*! \internal \brief
 * Converts, validates, and stores integer values.
 */
class IntegerOptionStorage : public OptionStorageTemplateSimple<int>
{
    public:
        //! \copydoc BooleanOptionStorage::BooleanOptionStorage()
        explicit IntegerOptionStorage(const IntegerOption &settings)
            : MyBase(settings), info_(this)
        {
        }

        virtual OptionInfo &optionInfo() { return info_; }
        virtual std::string typeString() const
        { return isVector() ? "vector" : "int"; }
        virtual std::string formatSingleValue(const int &value) const;

    private:
        virtual void initConverter(ConverterType *converter);
        virtual void processSetValues(ValueList *values);

        IntegerOptionInfo       info_;
};

/*! \internal \brief
 * Converts, validates, and stores integer values.
 */
class Int64OptionStorage : public OptionStorageTemplateSimple<gmx_int64_t>
{
    public:
        //! \copydoc BooleanOptionStorage::BooleanOptionStorage()
        explicit Int64OptionStorage(const Int64Option &settings)
            : MyBase(settings), info_(this)
        {
        }

        virtual OptionInfo &optionInfo() { return info_; }
        virtual std::string typeString() const { return "int"; }
        virtual std::string formatSingleValue(const gmx_int64_t &value) const;

    private:
        virtual void initConverter(ConverterType *converter);

        Int64OptionInfo       info_;
};

/*! \internal \brief
 * Converts, validates, and stores floating-point (double) values.
 */
class DoubleOptionStorage : public OptionStorageTemplateSimple<double>
{
    public:
        //! \copydoc IntegerOptionStorage::IntegerOptionStorage()
        explicit DoubleOptionStorage(const DoubleOption &settings);

        virtual OptionInfo &optionInfo() { return info_; }
        virtual std::string typeString() const;
        virtual std::string formatSingleValue(const double &value) const;

        //! \copydoc DoubleOptionInfo::isTime()
        bool isTime() const { return bTime_; }
        //! \copydoc DoubleOptionInfo::setScaleFactor()
        void setScaleFactor(double factor);

    private:
        virtual void initConverter(ConverterType *converter);
        virtual double processValue(const double &value) const;
        virtual void processSetValues(ValueList *values);

        DoubleOptionInfo        info_;
        bool                    bTime_;
        double                  factor_;
};

/*! \internal \brief
 * Converts, validates, and stores floating-point (float) values.
 */
class FloatOptionStorage : public OptionStorageTemplateSimple<float>
{
    public:
        //! \copydoc IntegerOptionStorage::IntegerOptionStorage()
        explicit FloatOptionStorage(const FloatOption &settings);

        virtual OptionInfo &optionInfo() { return info_; }
        virtual std::string typeString() const;
        virtual std::string formatSingleValue(const float &value) const;

        //! \copydoc DoubleOptionStorage::isTime()
        bool isTime() const { return bTime_; }
        //! \copydoc DoubleOptionStorage::setScaleFactor()
        void setScaleFactor(double factor);

    private:
        virtual void initConverter(ConverterType *converter);
        virtual float processValue(const float &value) const;
        virtual void processSetValues(ValueList *values);

        FloatOptionInfo         info_;
        bool                    bTime_;
        double                  factor_;
};

/*! \internal \brief
 * Converts, validates, and stores string values.
 */
class StringOptionStorage : public OptionStorageTemplateSimple<std::string>
{
    public:
        //! \copydoc DoubleOptionStorage::DoubleOptionStorage()
        explicit StringOptionStorage(const StringOption &settings);

        virtual OptionInfo &optionInfo() { return info_; }
        virtual std::string typeString() const
        { return allowed_.empty() ? "string" : "enum"; }
        virtual std::string formatExtraDescription() const;
        virtual std::string formatSingleValue(const std::string &value) const;

        //! \copydoc StringOptionInfo::allowedValues()
        const ValueList &allowedValues() const { return allowed_; }

    private:
        virtual void initConverter(ConverterType *converter);
        virtual std::string processValue(const std::string &value) const;

        StringOptionInfo        info_;
        ValueList               allowed_;
};

/*! \internal \brief
 * Converts, validates, and stores enum values.
 */
class EnumOptionStorage : public OptionStorageTemplateSimple<int>
{
    public:
        /*! \brief
         * Initializes the storage from option settings.
         *
         * \param[in] settings      Basic storage settings.
         * \param[in] enumValues    Allowed values.
         * \param[in] count         Number of elements in \p enumValues,
         *     or -1 if \p enumValues is `NULL`-terminated.
         * \param[in] defaultValue  Default value, or -1 if no default.
         * \param[in] defaultValueIfSet  Default value if set, or -1 if none.
         * \param[in] store         Storage to convert the values to/from `int`.
         *
         * This constructor takes more parameters than other storage parameters
         * because the front-end option type is a template, and as such cannot
         * be passed here without exposing also this header as an installed
         * header.
         */
        EnumOptionStorage(const AbstractOption &settings,
                          const char *const *enumValues, int count,
                          int defaultValue, int defaultValueIfSet,
                          StorePointer store);

        virtual OptionInfo &optionInfo() { return info_; }
        virtual std::string typeString() const { return "enum"; }
        virtual std::string formatExtraDescription() const;
        virtual std::string formatSingleValue(const int &value) const;
        virtual Variant normalizeValue(const int &value) const;

        //! \copydoc EnumOptionInfo::allowedValues()
        const std::vector<std::string> &allowedValues() const { return allowed_; }

    private:
        virtual void initConverter(ConverterType *converter);

        EnumOptionInfo            info_;
        std::vector<std::string>  allowed_;
};

/*!\}*/

} // namespace gmx

#endif
