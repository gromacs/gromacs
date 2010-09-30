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
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_BASICOPTIONSTORAGE_H
#define GMX_OPTIONS_BASICOPTIONSTORAGE_H

#include <string>
#include <vector>

#include "optionfiletype.h"
#include "optionstoragetemplate.h"

namespace gmx
{

class IntegerOption;
class DoubleOption;
class StringOption;
class FileNameOption;

/*! \addtogroup module_options
 * \{
 */

/*! \internal \brief
 * Converts, validates, and stores boolean values.
 */
class BooleanOptionStorage : public OptionStorageTemplate<bool>
{
    public:
        virtual const char *typeString() const { return "bool"; }
        virtual int appendValue(const std::string &value,
                                AbstractErrorReporter *errors);
        virtual int finishSet(int /*nvalues*/,
                              AbstractErrorReporter * /*errors*/)
        { return 0; }
        virtual std::string formatValue(int i) const;
};

/*! \internal \brief
 * Converts, validates, and stores integer values.
 */
class IntegerOptionStorage : public OptionStorageTemplate<int>
{
    public:
        IntegerOptionStorage();

        virtual const char *typeString() const
        { return hasFlag(efVector) ? "vector" : "int"; }
        virtual int appendValue(const std::string &value,
                                AbstractErrorReporter *errors);
        virtual int finishSet(int nvalues, AbstractErrorReporter *errors);
        virtual std::string formatValue(int i) const;
};

/*! \internal \brief
 * Converts, validates, and stores floating-point (double) values.
 */
class DoubleOptionStorage : public OptionStorageTemplate<double>
{
    public:
        DoubleOptionStorage();

        /*! \brief
         * Initializes the storage from option settings.
         *
         * \param[in] settings   Storage settings.
         * \param[in] options    Options object.
         * \retval 0 on success.
         */
        int init(const DoubleOption &settings, Options *options);

        virtual const char *typeString() const;
        virtual int appendValue(const std::string &value,
                                AbstractErrorReporter *errors);
        virtual int finishSet(int nvalues, AbstractErrorReporter *errors);
        virtual int finish(AbstractErrorReporter *errors);
        virtual std::string formatValue(int i) const;

    private:
        bool                    _bTime;
};

/*! \internal \brief
 * Converts, validates, and stores string values.
 */
class StringOptionStorage : public OptionStorageTemplate<std::string>
{
    public:
        StringOptionStorage();

        //! \copydoc DoubleOptionStorage::init()
        int init(const StringOption &settings, Options *options);

        virtual const char *typeString() const { return _allowed.empty() ? "string" : "enum"; }
        virtual int appendValue(const std::string &value,
                                AbstractErrorReporter *errors);
        virtual int finishSet(int /*nvalues*/,
                              AbstractErrorReporter * /*errors*/)
        { return 0; }
        virtual std::string formatValue(int i) const;

    private:
        ValueList               _allowed;
        int                    *_enumIndexStore;
};

/*! \internal \brief
 * Converts, validates, and stores file names.
 */
class FileNameOptionStorage : public OptionStorageTemplate<std::string>
{
    public:
        FileNameOptionStorage();

        //! \copydoc StringOptionStorage::init()
        int init(const FileNameOption &settings, Options *options);

        virtual const char *typeString() const { return "file"; }
        virtual int appendValue(const std::string &value,
                                AbstractErrorReporter *errors);
        virtual int finishSet(int /*nvalues*/,
                              AbstractErrorReporter * /*errors*/)
        { return 0; }
        virtual std::string formatValue(int i) const;

    private:
        OptionFileType          _filetype;
};

/*!\}*/

} // namespace gmx

#endif
