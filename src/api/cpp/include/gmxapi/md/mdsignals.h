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

#ifndef GMXAPI_MDSIGNALS_H
#define GMXAPI_MDSIGNALS_H

/*!
 * \file
 * \brief Temporary infrastructure for signalling MD simulations.
 *
 * These interfaces are not considered to be part of the gmxapi spec, but will exist in 0.0.6 and
 * possibly 0.0.7 until more abstract data flow is available to MD plugin developers, at which point
 * any remaining functionality here will be moved to private implementation details.
 *
 * \ingroup gmxapi_md
 */

#include <memory>

namespace gmxapi
{

/*!
 * \brief Internal details of gmxapi MD functionality
 */
namespace md
{

/*!
 * \brief Symbolic signal slots for MD signalling.
 *
 * \see getMdrunnerSignal()
 */
enum class signals
{
    STOP
};

}   // end namespace md


class SessionResources;  // reference gmxapi/session/resources.h

/*!
 * \brief Proxy for signalling function objects.
 *
 * Objects of this type are simple callables that issue a specific signal.
 *
 * \todo This class and its implementations should be replaced with a std::function.
 *
 * \ingroup gmxapi_md
 */
class Signal
{
    public:
        //! \internal
        class SignalImpl;

        /*!
         * \brief Construct by taking ownership of an implementation object.
         *
         * \param signal
         */
        explicit Signal(std::unique_ptr<SignalImpl> signal);

        /*!
         * \brief Object is trivially moveable.
         *
         * \{
         */
        Signal(Signal &&) noexcept;
        Signal &operator=(Signal &&) noexcept;
        //! \}

        //! \cond
        ~Signal();
        //! \endcond

        /*!
         * \brief Signal interface is a function object that issues a signal when called.
         *
         * \todo replace with more concise named type based on std::function.
         */
        void operator()();

    private:
        //! Wrapped signaller.
        std::unique_ptr<SignalImpl> impl_;
};

/*!
 * \brief Get a function object that issues a signal to the currently active MD runner.
 *
 * For use during execution of extension code. E.g. within a MD simulation plugin
 * during an MD step.
 *
 * \param resources non-null pointer to the active Session resources.
 * \param signal type of signal the client would like to issue.
 * \return Callable function object handle to issue a stop signal.
 *
 * \throws gmxapi::NotImplementedError for unknown values of signal.
 * \throws gmxapi::UsageError for invalid resources argument.
 *
 * \ingroup gmxapi_md
 */
Signal getMdrunnerSignal(SessionResources *resources,
                         md::signals       signal);

}      // end namespace gmxapi

#endif //GMXAPI_MDSIGNALS_H
