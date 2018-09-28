/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
#ifndef GMXAPI_SYSTEM_H
#define GMXAPI_SYSTEM_H
/*! \file
 * \brief Declare container for molecular systems
 *
 * \ingroup gmxapi
 */
#include <memory>

namespace gmxapi
{

/*!
 * \brief Deprecated: A wrapper for gmx::Mdrunner
 *
 * Container for molecular model and simulation parameters.
 *
 * \internal
 * It was not intended as such, but gmxapi::System has ended up as basically a
 * wrapper for the gmx::Mdrunner simulator object and serves as the aspect of a
 * gmxapi Session that performs simulation work. As such, this class does not fit the
 * current gmxapi paradigm and will be removed, reworked, or renamed soon.
 *
 * # Protocol
 *
 * In this version, the System class is a non-functioning placeholder.
 *
 * \ingroup gmxapi
 */
class System final
{
    public:
        /*! \brief Private implementation class.
         *
         * System::Impl does not have a public interface and is only exposed in opaque pointers.
         */
        class Impl;

        /*! \brief No copy.
         *
         * The semantics of copying a System are ambiguous, so disallow implicit
         * copy.
         * \internal
         * Some sort of prototype or clone idiom could be useful, but
         * needs to explicitly identify any expensive operations.
         * \{
         */
        System(const System &)                         = delete;
        System              &operator=(const System &) = delete;
        /*! \} */

        /*! \brief Allow move.
         *
         * \{
         */
        System(System &&) noexcept;
        System &operator=(System &&) noexcept;
        /*! \} */

        /*!
         * \brief Create by taking ownership of an implementation object.
         *
         * \param implementation
         */
        explicit System(std::unique_ptr<Impl> implementation);

        /*! \cond internal
         *  Destructor defined later to allow unique_ptr members of partially-defined types.
         */
        ~System();
        /*! \endcond */

    private:
        /*!
         * \brief Opaque pointer to implementation.
         */
        std::unique_ptr<Impl> impl_;
};


/*! \brief Defines an MD workflow from a TPR file.
 *
 * The TPR file has sufficient information to fully specify an MD run, though
 * various parameters are implicit until the work is launched. The TPR filename
 * provided must refer to identical TPR files at the API client and at the
 * master rank of the execution host.
 *
 * \param filename Filesystem path of TPR file.
 * \returns gmxapi::System object with the specified workflow.
 * \ingroup gmxapi
 */
System fromTprFile(std::string filename);

}      // end namespace gmxapi

#endif // include guard
