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
    explicit BooleanOptionStorage(const BooleanOption& settings) : MyBase(settings), info_(this) {}

    OptionInfo& optionInfo() override { return info_; }
    std::string typeString() const override { return "bool"; }
    std::string formatSingleValue(const bool& value) const override;

    //! \copydoc BooleanOptionInfo::defaultValue()
    bool defaultValue() const { return valueCount() > 0 && values()[0]; }

private:
    void initConverter(ConverterType* converter) override;

    BooleanOptionInfo info_;
};

/*! \internal \brief
 * Converts, validates, and stores integer values.
 */
class IntegerOptionStorage : public OptionStorageTemplateSimple<int>
{
public:
    //! \copydoc BooleanOptionStorage::BooleanOptionStorage()
    explicit IntegerOptionStorage(const IntegerOption& settings) : MyBase(settings), info_(this) {}

    OptionInfo& optionInfo() override { return info_; }
    std::string typeString() const override { return isVector() ? "vector" : "int"; }
    std::string formatSingleValue(const int& value) const override;

private:
    void initConverter(ConverterType* converter) override;
    void processSetValues(ValueList* values) override;

    IntegerOptionInfo info_;
};

/*! \internal \brief
 * Converts, validates, and stores integer values.
 */
class Int64OptionStorage : public OptionStorageTemplateSimple<int64_t>
{
public:
    //! \copydoc BooleanOptionStorage::BooleanOptionStorage()
    explicit Int64OptionStorage(const Int64Option& settings) : MyBase(settings), info_(this) {}

    OptionInfo& optionInfo() override { return info_; }
    std::string typeString() const override { return "int"; }
    std::string formatSingleValue(const int64_t& value) const override;

private:
    void initConverter(ConverterType* converter) override;

    Int64OptionInfo info_;
};

/*! \internal \brief
 * Converts, validates, and stores floating-point (double) values.
 */
class DoubleOptionStorage : public OptionStorageTemplateSimple<double>
{
public:
    //! \copydoc IntegerOptionStorage::IntegerOptionStorage()
    explicit DoubleOptionStorage(const DoubleOption& settings);

    OptionInfo& optionInfo() override { return info_; }
    std::string typeString() const override;
    std::string formatSingleValue(const double& value) const override;

    //! \copydoc DoubleOptionInfo::isTime()
    bool isTime() const { return bTime_; }
    //! \copydoc DoubleOptionInfo::setScaleFactor()
    void setScaleFactor(double factor);

private:
    void   initConverter(ConverterType* converter) override;
    double processValue(const double& value) const override;
    void   processSetValues(ValueList* values) override;

    DoubleOptionInfo info_;
    bool             bTime_;
    double           factor_;
};

/*! \internal \brief
 * Converts, validates, and stores floating-point (float) values.
 */
class FloatOptionStorage : public OptionStorageTemplateSimple<float>
{
public:
    //! \copydoc IntegerOptionStorage::IntegerOptionStorage()
    explicit FloatOptionStorage(const FloatOption& settings);

    OptionInfo& optionInfo() override { return info_; }
    std::string typeString() const override;
    std::string formatSingleValue(const float& value) const override;

    //! \copydoc DoubleOptionStorage::isTime()
    bool isTime() const { return bTime_; }
    //! \copydoc DoubleOptionStorage::setScaleFactor()
    void setScaleFactor(double factor);

private:
    void  initConverter(ConverterType* converter) override;
    float processValue(const float& value) const override;
    void  processSetValues(ValueList* values) override;

    FloatOptionInfo info_;
    bool            bTime_;
    double          factor_;
};

/*! \internal \brief
 * Converts, validates, and stores string values.
 */
class StringOptionStorage : public OptionStorageTemplateSimple<std::string>
{
public:
    //! \copydoc DoubleOptionStorage::DoubleOptionStorage()
    explicit StringOptionStorage(const StringOption& settings);

    OptionInfo& optionInfo() override { return info_; }
    std::string typeString() const override { return allowed_.empty() ? "string" : "enum"; }
    std::string formatExtraDescription() const override;
    std::string formatSingleValue(const std::string& value) const override;

    //! \copydoc StringOptionInfo::allowedValues()
    const ValueList& allowedValues() const { return allowed_; }

private:
    void        initConverter(ConverterType* converter) override;
    std::string processValue(const std::string& value) const override;

    StringOptionInfo info_;
    ValueList        allowed_;
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
    EnumOptionStorage(const AbstractOption& settings,
                      const char* const*    enumValues,
                      int                   count,
                      int                   defaultValue,
                      int                   defaultValueIfSet,
                      StorePointer          store);

    OptionInfo& optionInfo() override { return info_; }
    std::string typeString() const override { return "enum"; }
    std::string formatExtraDescription() const override;
    std::string formatSingleValue(const int& value) const override;
    Any         normalizeValue(const int& value) const override;

    //! \copydoc EnumOptionInfo::allowedValues()
    const std::vector<std::string>& allowedValues() const { return allowed_; }

private:
    void initConverter(ConverterType* converter) override;

    EnumOptionInfo           info_;
    std::vector<std::string> allowed_;
};

/*!\}*/

} // namespace gmx

#endif
