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
 * Declares option settings objects for basic option types.
 *
 * Together with options.h, this header forms the part of the public API
 * that most classes will use to provide options.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_BASICOPTIONS_H
#define GMX_OPTIONS_BASICOPTIONS_H

#include <cassert>

#include <string>

#include "abstractoption.h"
#include "optionfiletype.h"

namespace gmx
{

class IntegerOptionStorage;
class DoubleOptionStorage;
class StringOptionStorage;
class FileNameOptionStorage;
class Options;

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
 * \inpublicapi
 */
class BooleanOption : public OptionTemplate<bool, BooleanOption>
{
    public:
        //! Initializes an option with the given name.
        explicit BooleanOption(const char *name) : MyBase(name)
        {
            setFlag(efBoolean);
        }

    protected:
        virtual int createDefaultStorage(Options *options,
                                         OptionStorageInterface **storage) const;
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
 * \inpublicapi
 */
class IntegerOption : public OptionTemplate<int, IntegerOption>
{
    public:
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

    protected:
        virtual int createDefaultStorage(Options *options,
                                         OptionStorageInterface **storage) const;

        /*! \brief
         * Needed to initialize IntegerOptionStorage from this class without
         * otherwise unnecessary accessors.
         */
        friend class IntegerOptionStorage;
};

/*! \brief
 * Specifies an option that provides floating-point (double) values.
 *
 * \inpublicapi
 */
class DoubleOption : public OptionTemplate<double, DoubleOption>
{
    public:
        //! Initializes an option with the given name.
        explicit DoubleOption(const char *name) : MyBase(name), _bTime(false)
        {
        }

        //! \copydoc IntegerOption::vector()
        MyClass &vector() { setVector(); return me(); }
        /*! \brief
         * Sets the option to obey time conversion rules.
         */
        MyClass &timeValue() { _bTime = true; return me(); }

    private:
        virtual int createDefaultStorage(Options *options,
                                         OptionStorageInterface **storage) const;

        bool _bTime;

        /*! \brief
         * Needed to initialize DoubleOptionStorage from this class without
         * otherwise unnecessary accessors.
         */
        friend class DoubleOptionStorage;
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
const char * const  allowed[] = { "atom", "residue", "molecule", NULL };
std::string  str;
int          type;
options.addOption(StringOption("type").enumValue(allowed).store(&str)
                     .storeEnumIndex(&type));
 * \endcode
 *
 * \inpublicapi
 */
class StringOption : public OptionTemplate<std::string, StringOption>
{
    public:
        //! Initializes an option with the given name.
        explicit StringOption(const char *name)
            : MyBase(name), _enumValues(NULL), _defaultEnumIndex(-1),
              _enumIndexStore(NULL)
        {}

        /*! \brief
         * Sets the option to only accept one of a fixed set of strings.
         *
         * \param[in] values  Array of strings to accept, terminated with a
         *     NULL value.
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
        MyClass &enumValue(const char *const *values)
        { _enumValues = values; return me(); }
        /*! \brief
         * Sets the default value using an index into the enumeration table.
         *
         * Cannot be specified without enumValue().
         */
        MyClass &defaultEnumIndex(int index)
        { assert(index >= 0); _defaultEnumIndex = index; return me(); }
        /*! \brief
         * Stores the index of the selected value into the provided memory
         * location.
         *
         * The index (zero-based) of the selected value in the array \p values
         * provided to enumValues() is written into \p *store after the
         * option gets its value.  If the option has not been provided,
         * and there is no default value, -1 is stored.
         *
         * Cannot be specified without enumValue().
         */
        MyClass &storeEnumIndex(int *store)
        { _enumIndexStore = store; return me(); }

    protected:
        virtual int createDefaultStorage(Options *options,
                                         OptionStorageInterface **storage) const;
        virtual std::string createDescription() const;

    private:
        const char *const      *_enumValues;
        int                     _defaultEnumIndex;
        int                    *_enumIndexStore;

        /*! \brief
         * Needed to initialize StringOptionStorage from this class without
         * otherwise unnecessary accessors.
         */
        friend class StringOptionStorage;
};

/*! \brief
 * Specifies an option that provides file names.
 *
 * \inpublicapi
 */
class FileNameOption : public OptionTemplate<std::string, FileNameOption>
{
    public:
        //! Initializes an option with the given name.
        explicit FileNameOption(const char *name)
            : MyBase(name), _filetype(eftUnknown)
        {
            setFlag(efFile);
        }

        /*! \brief
         * Sets the type of the file this option accepts.
         *
         * This attribute must be provided.
         */
        MyClass &filetype(OptionFileType type)
        { _filetype = type; return me(); }
        //! Tells that the file provided by this option is used read-only.
        MyClass &readOnly()
        { setFlag(efFileRead); clearFlag(efFileWrite); return me(); }
        //! Tells that the file provided by this option is used write-only.
        MyClass &writeOnly()
        { setFlag(efFileWrite); clearFlag(efFileRead); return me(); }
        /*! \brief
         * Tells that the file provided by this option is used for reading and
         * writing.
         */
        MyClass &readwrite()
        { setFlag(efFileRead); setFlag(efFileWrite); return me(); }
        /*! \brief
         * Tells that the file will be looked up in library directories in
         * addition to working directory.
         */
        MyClass &libraryFile() { setFlag(efFileLibrary); return me(); }

    protected:
        virtual int createDefaultStorage(Options *options,
                                         OptionStorageInterface **storage) const;

    private:
        OptionFileType          _filetype;

        /*! \brief
         * Needed to initialize FileNameOptionStorage from this class without
         * otherwise unnecessary accessors.
         */
        friend class FileNameOptionStorage;
};

/*!\}*/

} // namespace gmx

#endif
