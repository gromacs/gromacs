/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Declares gmx::IOptionsContainer.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_IOPTIONSCONTAINER_H
#define GMX_OPTIONS_IOPTIONSCONTAINER_H

#include "gromacs/options/abstractoption.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

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
class IOptionManager
{
    protected:
        virtual ~IOptionManager();
};

/*! \brief
 * Interface for adding input options.
 *
 * This interface provides methods to add new options.
 * Standard usage is for code to receive this interface and populate it with
 * supported options:
 * \code
   // <as class attributes>
   std::string  arg1;
   int          arg2;

   void MyClass::initOptions(gmx::IOptionsContainer *options)
   {
       options->addOption(gmx::StringOption("arg1").store(&arg1));
       options->addOption(gmx::IntegerOption("arg2").store(&arg2));
   }
   \endcode
 * The caller can collect options from multiple sources into a single container
 * (a gmx::Options), and use a parser implementation such as CommandLineParser
 * to provide values for the options.
 *
 * Header basicoptions.h provides declarations of several standard
 * option types for use with addOption().  Documentation of those classes
 * also give more examples of how to define options.
 *
 * \inpublicapi
 * \ingroup module_options
 */
class IOptionsContainer
{
    public:
        /*! \brief
         * Adds a recognized option.
         *
         * \param[in] settings Option description.
         * \returns   OptionInfo object for the created option (never NULL).
         * \throws    APIError if invalid option settings are provided.
         *
         * This method provides the internal implementation, but in most cases
         * the templated method is called from user code.
         * See the templated method for more details.
         */
        virtual OptionInfo *addOption(const AbstractOption &settings) = 0;
        /*! \brief
         * Adds a recognized option.
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
        virtual void addManager(IOptionManager *manager) = 0;

        /*! \brief
         * Creates a subgroup of options within the current options.
         *
         * To add options to the group, use the returned interface.
         *
         * Currently, this is only used to influence the order of options:
         * all options in a group appear before options in a group added after
         * it, no matter in which order the options are added to the groups.
         * In the future, the groups could also be used to influence the help
         * output.
         */
        virtual IOptionsContainer &addGroup() = 0;
    protected:
        // Disallow deletion through the interface.
        // (no need for the virtual, but some compilers warn otherwise)
        virtual ~IOptionsContainer();

};

} // namespace

#endif
