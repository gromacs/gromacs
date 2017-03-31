/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * \brief Declares wrappers for gmx::Options and such.
 * \ingroup module_python
 */
#ifndef PYGMX_OPTIONS_H
#define PYGMX_OPTIONS_H


#include "gmxpre.h"

#include <string>

#include "gromacs/options/options.h"

namespace gmx
{
class OptionsVisitor;

/*! \brief API client code from which to export Python bindings
 *
 * \internal
 * \ingroup module_python
 */
namespace pyapi
{
using std::shared_ptr;
using std::unique_ptr;
using std::make_shared;

/*! \brief Wraps an Options collection for exporting to Python.
 *
 * \internal \ingroup module_python
 */
class PyOptions
{
    public:
        /// Create an empty options container.
        PyOptions();
        /*! \brief Create an options container with our only known option.
         *
         * This is a temporary proof-of-concept shim until proper options handling
         * is implemented.
         * \param filename trajectory filename
         */
        PyOptions(std::string filename);
        /// clean up
        ~PyOptions();
        /// don't yet have a use for copies
        PyOptions(const PyOptions &)                  = delete;
        /// don't yet have a use for copies, though move would be handy
        const PyOptions &operator=(const PyOptions &) = delete;
        // Copy semantics seem likely to involve multiple pointers to the same object rather than copies of the options object, but we'll see...
        // gmx::Options objects have implementation members that look like they are not meant to be copied...

        /// Get a raw pointer to the member data.
        /*! \internal \ingroup module_python
         *
         * \return raw pointer to the gmx::Options object owned by *this.
         */
        gmx::Options* data();

        /// Provide a manager for OptionsVisitors
        /*! This method is misleadingly named, because the visitor provides its
         *  own iterator with which to traverse the options tree.
         *  Visitor may modify itself during traversal.
         * \param visitor object will receive the root gmx::OptionSectionInfo.
         */
        void view_traverse(gmx::OptionsVisitor &&visitor) const;
//    void modify_traverse(gmx::OptionsVisitor& visitor);

        /*! \brief Process input for registered Options.
         *
         * No way to access processing errors yet...
         * \return true if successful, else false.
         */
        bool parse();

    private:
        /// Wrapped gmx::Options object
        gmx::Options options_;
        /// Right now, the only option we address
        std::string  filename_;
};

/*! \brief Apply a OptionsVisitor that prints out the contents of the Options collection.
 *
 * \param pyoptions options tree to be inspected.
 * \internal \ingroup module_python
 */
void print_options(const PyOptions &pyoptions); // can't, really...

};                                              //end namespace pyapi
};                                              // end namespace gmx
#endif                                          // header guard
