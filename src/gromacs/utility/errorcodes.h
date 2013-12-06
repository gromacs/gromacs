/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013, by the GROMACS development team, led by
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
 * Declares error codes and related functions for fatal error handling.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_ERRORCODES_H
#define GMX_UTILITY_ERRORCODES_H

namespace gmx
{

/*! \addtogroup module_utility
 * \{
 */

/*! \brief
 * Possible error return codes from Gromacs functions.
 */
enum ErrorCode
{
    //! Zero for successful return.
    eeOK,

    //! Not enough memory to complete operation.
    eeOutOfMemory,
    //! Provided file could not be opened.
    eeFileNotFound,
    //! System I/O error.
    eeFileIO,
    //! Invalid user input (could not be understood).
    eeInvalidInput,
    //! Invalid user input (conflicting or unsupported settings).
    eeInconsistentInput,
    //! Simulation instability detected.
    eeInstability,

    /*! \name Error codes for buggy code
     *
     * Error codes below are for internal error checking; if triggered, they
     * should indicate a bug in the code.
     * \{
     */
    //! Requested feature not yet implemented.
    eeNotImplemented,
    //! Input value violates API specification.
    eeInvalidValue,
    //! Invalid routine called or wrong calling sequence detected.
    eeInvalidCall,
    //! Internal consistency check failed.
    eeInternalError,
    //! API specification was violated.
    eeAPIError,
    //! Range consistency check failed.
    eeRange,
    //! Communication consistency check failed.
    eeCommunication,
    //!\}

    //! Unknown error detected.
    eeUnknownError,
};

/*! \brief
 * Returns a short string description of an error code.
 *
 * \param[in] errorcode Error code to retrieve the string for.
 * \returns   A constant string corresponding to \p errorcode.
 *
 * If \p errorcode is not one of those defined for ::gmx::ErrorCode,
 * the string corresponding to ::eeUnknownError is returned.
 *
 * This function does not throw.
 */
const char *getErrorCodeString(int errorcode);

/*! \brief
 * Callback function pointer type for error handlers.
 *
 * \param[in] retcode Code of the error that has occurred.
 * \param[in] msg     More detailed description of the error.
 * \param[in] file    Name of the file where the error occurred.
 * \param[in] line    Line in \p file on which the error occurred.
 */
typedef void (*ErrorHandlerFunc)(int retcode, const char *msg,
                                 const char *file, int line);

/*! \brief
 * Sets callback function for handling errors.
 *
 * \param[in] handler New error handler function.
 * \returns   Old error handler function.
 *
 * The default error handler prints out the location and reason of the error to
 * stderr, and then calls std::abort().
 */
ErrorHandlerFunc setFatalErrorHandler(ErrorHandlerFunc handler);

/*! \brief
 * Macro for raising an error and returning from a function.
 *
 * The function should return \c int.
 * If it doesn't, use ::GMX_ERROR_NORET.
 */
#define GMX_ERROR(retcode, msg) \
    do { \
        int _rc_internal = (retcode); \
        ::gmx::internal::fatalError(_rc_internal, msg, __FILE__, __LINE__); \
        return _rc_internal; \
    } while (0)

/*! \brief
 * Macro for raising an error in a function that does not return \c int.
 *
 * \see GMX_ERROR
 */
#define GMX_ERROR_NORET(retcode, msg) \
        ::gmx::internal::fatalError(retcode, msg, __FILE__, __LINE__)

/*!\}*/

/*! \cond internal */
namespace internal
{

/*! \internal \brief
 * Raises a fatal error.
 *
 * \param[in] retcode Error code to raise.
 * \param[in] msg     More detailed description of the error.
 * \param[in] file    Name of the source file where the error occurred.
 * \param[in] line    Line in \p file on which the error occurred.
 *
 * Should not be called directly, but instead through ::GMX_ERROR or
 * ::GMX_ERROR_NORET.
 *
 * \ingroup module_utility
 */
void fatalError(int retcode, const char *msg, const char *file, int line);

}   // namespace internal
//! \endcond

} // namespace gmx

#endif
