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
/*! \internal \file
 * \brief
 * Declares storage classes for basic option types.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_BASICOPTIONSTORAGE_H
#define GMX_OPTIONS_BASICOPTIONSTORAGE_H

#include <string>
#include <vector>

#include "basicoptions.h"
#include "optionstoragetemplate.h"

namespace gmx
{

class BooleanOption;
class IntegerOption;
class DoubleOption;
class StringOption;

/*! \addtogroup module_options
 * \{
 */

/*! \internal \brief
 * Converts, validates, and stores boolean values.
 */
class BooleanOptionStorage : public OptionStorageTemplate<bool>
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
        virtual const char *typeString() const { return "bool"; }
        virtual std::string formatSingleValue(const bool &value) const;

    private:
        virtual void convertValue(const std::string &value);

        BooleanOptionInfo       info_;
};

/*! \internal \brief
 * Converts, validates, and stores integer values.
 */
class IntegerOptionStorage : public OptionStorageTemplate<int>
{
    public:
        //! \copydoc BooleanOptionStorage::BooleanOptionStorage()
        explicit IntegerOptionStorage(const IntegerOption &settings)
            : MyBase(settings), info_(this)
        {
        }

        virtual OptionInfo &optionInfo() { return info_; }
        virtual const char *typeString() const
        { return isVector() ? "vector" : "int"; }
        virtual std::string formatSingleValue(const int &value) const;

    private:
        virtual void convertValue(const std::string &value);
        virtual void processSetValues(ValueList *values);

        IntegerOptionInfo       info_;
};

/*! \internal \brief
 * Converts, validates, and stores floating-point (double) values.
 */
class DoubleOptionStorage : public OptionStorageTemplate<double>
{
    public:
        //! \copydoc IntegerOptionStorage::IntegerOptionStorage()
        explicit DoubleOptionStorage(const DoubleOption &settings);

        virtual OptionInfo &optionInfo() { return info_; }
        virtual const char *typeString() const;
        virtual std::string formatSingleValue(const double &value) const;

        //! \copydoc DoubleOptionInfo::isTime()
        bool isTime() const { return bTime_; }
        //! \copydoc DoubleOptionInfo::setScaleFactor()
        void setScaleFactor(double factor);

    private:
        virtual void convertValue(const std::string &value);
        virtual void processSetValues(ValueList *values);
        virtual void processAll();

        DoubleOptionInfo        info_;
        bool                    bTime_;
        double                  factor_;
};

/*! \internal \brief
 * Converts, validates, and stores string values.
 */
class StringOptionStorage : public OptionStorageTemplate<std::string>
{
    public:
        //! \copydoc DoubleOptionStorage::DoubleOptionStorage()
        explicit StringOptionStorage(const StringOption &settings);

        virtual OptionInfo &optionInfo() { return info_; }
        virtual const char *typeString() const { return allowed_.empty() ? "string" : "enum"; }
        virtual std::string formatSingleValue(const std::string &value) const;

    private:
        virtual void convertValue(const std::string &value);
        virtual void refreshValues();

        StringOptionInfo        info_;
        ValueList               allowed_;
        int                    *enumIndexStore_;
};

/*!\}*/

} // namespace gmx

#endif
