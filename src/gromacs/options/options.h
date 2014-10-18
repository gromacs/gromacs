/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2014, by the GROMACS development team, led by
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
 * Declares gmx::Options.
 *
 * Together with basicoptions.h, this header forms the part of the public
 * API that most classes will use to provide options.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_OPTIONS_H
#define GMX_OPTIONS_OPTIONS_H

#include <string>

#include "gromacs/options/abstractoption.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

template <typename T> class ConstArrayRef;

class AbstractOption;
class OptionsAssigner;
class OptionsIterator;

/*! \brief
 * Base class for option managers.
 *
 * This class is used as a marker for all classes that are used with
 * Options::addManager().  It doesn't provide any methods, but only supports
 * transporting these classes through the Options collection into the
 * individual option implementation classes.
 *
 * The virtual destructor is present to make this class polymorphic, such that
 * `dynamic_cast` can be used when retrieving a manager of a certain type for
 * the individual options.
 *
 * \inlibraryapi
 * \ingroup module_options
 */
class OptionManagerInterface
{
    protected:
        virtual ~OptionManagerInterface();
};

/*! \brief
 * Collection of options.
 *
 * This class provides a standard interface for implementing input options.
 * Standard usage is to write a method that creates an Options that is owned by
 * the object, populates it with supported options, and then returns it:
 * \code
   // <as class attributes>
   using gmx::Options;
   Options      options("common", "Common Options");
   std::string  arg1;
   int          arg2;

   // <populating>
   using gmx::StringOption;
   using gmx::IntegerOption;
   options.addOption(StringOption("arg1").store(&arg1));
   options.addOption(IntegerOption("arg2").store(&arg2));
   return &options;
   \endcode
 * The caller of that method can then use a parser implementation such as
 * CommandLineParser to provide values for the options.
 *
 * Header basicoptions.h provides declarations of several standard
 * option types for use with addOption().  Documentation of those classes
 * also give more examples of how to define options.
 *
 * In order to keep the public interface of this class simple and to reduce
 * build dependencies on objects that simply provide options, functionality
 * to assign values to options is provided by a separate OptionsAssigner class.
 * Similarly, functionality for looping over all options (e.g., for writing out
 * help) is provided by OptionsIterator.
 *
 * \inpublicapi
 * \ingroup module_options
 */
class Options
{
    public:
        /*! \brief
         * Initializes the name and title of an option collection.
         *
         * \param[in] name  Single-word name.
         * \param[in] title Descriptive title.
         *
         * Copies the input strings.
         */
        Options(const char *name, const char *title);
        ~Options();

        //! Returns the short name of the option collection.
        const std::string &name() const;
        //! Returns the title of the option collection.
        const std::string &title() const;
        //! Returns the full description of the option collection.
        const std::string &description() const;

        /*! \brief
         * Sets the full description of the option collection.
         *
         * \param[in] desc  String to set as the description.
         *
         * This overload is mainly useful if the description is very short.
         * Currently this is mostly the case in unit testing.
         */
        void setDescription(const std::string &desc);
        /*! \brief
         * Sets the full description of the option collection from string array.
         *
         * \param[in] descArray  String array to set as the description.
         *
         * All strings in `descArray` are concatenated to form the description.
         * Spaces are inserted between strings if they are missing.
         *
         * Example usage:
         * \code
           const char *const desc[] = {
               "This is the description",
               "for the options"
           };

           gmx::Options options(NULL, NULL);
           options.setDescription(desc);
           \endcode
         *
         * To use this overload, you must also include
         * `gromacs/utility/arrayref.h`.
         */
        void setDescription(const ConstArrayRef<const char *> &descArray);

        /*! \brief
         * Adds an option manager.
         *
         * \param    manager Manager to add.
         * \throws   std::bad_alloc if out of memory.
         *
         * Option managers are used by some types of options that require
         * interaction between different option instances (e.g., selection
         * options), or need to support globally set properties (e.g., a global
         * default file prefix).  Option objects can retrieve the pointer to
         * their manager when they are created, and the caller can alter the
         * behavior of the options through the manager.
         * See the individual managers for details.
         *
         * Caller is responsible for memory management of \p manager.
         * The Options object (and its contained options) only stores a
         * reference to the object.
         *
         * This method cannot be called after adding options or subsections.
         */
        void addManager(OptionManagerInterface *manager);

        /*! \brief
         * Adds an option collection as a subsection of this collection.
         *
         * \param[in] section Subsection to add.
         *
         * The name() field of \p section is used as the name of the
         * subsection.  If an attempt is made to add two different subsections
         * with the same name, this function asserts.
         *
         * \p section should not have any options added at the point this
         * method is called.
         *
         * Only a pointer to the provided object is stored.  The caller is
         * responsible that the object exists for the lifetime of the
         * collection.
         * It is not possible to add the same Options object as a subsection to
         * several different Options.
         * If an attempt is made, the function asserts.
         */
        void addSubSection(Options *section);
        /*! \brief
         * Adds a recognized option to the collection.
         *
         * \param[in] settings Option description.
         * \returns   OptionInfo object for the created option (never NULL).
         * \throws    APIError if invalid option settings are provided.
         *
         * This method provides the internal implementation, but in most cases
         * the templated method is called from user code.
         * See the templated method for more details.
         */
        OptionInfo *addOption(const AbstractOption &settings);
        /*! \brief
         * Adds a recognized option to the collection.
         *
         * \tparam    OptionType Type of the options description object.
         * \param[in] settings   Option description.
         * \returns   OptionInfo object for the created option (never NULL).
         * \throws    APIError if invalid option settings are provided.
         *
         * The return value is a pointer for more convenient use in callers:
         * often callers need to declare the variable that will hold the return
         * value in wider scope than would be achieved by declaring it at the
         * site where addOption() is called.
         * The returned pointer must not be freed.
         *
         * See \link Options class documentation \endlink for example usage.
         *
         * \libinternal
         * \p OptionType::InfoType must specify a type that derives from
         * OptionInfo and matches the type that is returned by
         * AbstractOptionStorage::optionInfo() for the storage object that
         * corresponds to \p OptionType.
         */
        template <class OptionType>
        typename OptionType::InfoType *addOption(const OptionType &settings)
        {
            OptionInfo *info
                = addOption(static_cast<const AbstractOption &>(settings));
            GMX_ASSERT(info->isType<typename OptionType::InfoType>(),
                       "Mismatching option info type declaration and implementation");
            return info->toType<typename OptionType::InfoType>();
        }

        //! Returns true if option \p name is set.
        bool isSet(const char *name) const;
        /*! \brief
         * Notifies the collection that all option values are assigned.
         *
         * \throws InvalidInputError if invalid user input is detected.
         *
         * This function should be called after no more option values are
         * to be assigned.  Values in storage variables are guaranteed to be
         * available only after this call, although in most cases, they are
         * available already during assignment.
         *
         * If invalid option values, e.g., missing required option, is detected
         * at this point, this function throws.  The thrown exception contains
         * information on all errors detected during the call.
         */
        void finish();

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;

        //! Needed to be able to extend the interface of this object.
        friend class OptionsAssigner;
        //! Needed to be able to extend the interface of this object.
        friend class OptionsIterator;
        //! Needed to be able to extend the interface of this object.
        friend class OptionsModifyingIterator;
};

} // namespace gmx

#endif
