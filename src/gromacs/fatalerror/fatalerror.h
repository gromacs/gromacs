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
 *
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
 * Declares common return codes and functions for fatal error handling.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_fatalerror
 */
#ifndef GMX_FATALERROR_FATALERROR_H
#define GMX_FATALERROR_FATALERROR_H

namespace gmx
{

/*! \addtopublicapi
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

    // Error codes below are for internal error checking; if triggered, they
    // should indicate a bug in the code.
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
 * stderr, and then calls abort().
 */
ErrorHandlerFunc setFatalErrorHandler(ErrorHandlerFunc handler);

/*! \brief
 * Raises a fatal error.
 *
 * \param[in] retcode Error code to raise.
 * \param[in] msg     More detailed description of the error.
 * \param[in] file    Name of the source file where the error occurred.
 * \param[in] line    Line in \p file on which the error occurred.
 */
void fatalError(int retcode, const char *msg, const char *file, int line);
/*! \brief
 * Raises an error with a formatted message.
 *
 * Use like this:
 * \code
::gmx::fatalErrorFormatted(::gmx::eeInvalidInput, GMX_ERRORLOC,
    "Invalid command-line argument: %s", argname);
 * \endcode
 *
 * \param[in] retcode Error code to raise.
 * \param[in] file    Name of the source file where the error occurred.
 * \param[in] line    Line in \p file on which the error occurred.
 * \param[in] fmt     printf format string.
 */
void fatalErrorFormatted(int retcode, const char *file, int line,
                         const char *fmt, ...);

/*! \brief
 * Helper macro for fatalErrorFormatted().
 *
 * \see fatalErrorFormatted()
 */
#define GMX_ERRORLOC   __FILE__, __LINE__

/*! \brief
 * Macro for raising an error and returning from a function.
 *
 * The function should return \c int ; if it doesn't, use GMX_ERROR_NORET.
 */
#define GMX_ERROR(retcode, msg) \
    do { \
        int _rc_internal = (retcode); \
        ::gmx::fatalError(_rc_internal, msg, __FILE__, __LINE__); \
        return _rc_internal; \
    } while (0)

/*! \brief
 * Macro for raising an error in a function that does not return \c int.
 *
 * \see GMX_ERROR
 */
#define GMX_ERROR_NORET(retcode, msg) \
        ::gmx::fatalError(retcode, msg, __FILE__, __LINE__)

/*! \def GMX_ERROR_DEBUG
 * \brief
 * Macro for raising a debugging error and returning from a function.
 *
 * If NDEBUG is defined, this macro only returns the given error code.
 * If it is not defined, behaves exactly like GMX_ERROR.
 *
 * The function should return \c int ; if it doesn't, use GMX_ERROR_DEBUG_NORET.
 */
/*! \def GMX_ERROR_DEBUG_NORET
 * \brief
 * Macro for raising a debugging error in a function that does not return \c int.
 *
 * If NDEBUG is defined, this macro does nothing.  If it is not defined,
 * behaves exactly like GMX_ERROR_NORET.
 */
#ifdef NDEBUG
#define GMX_ERROR_DEBUG(retcode, msg) \
        return (retcode)
#define GMX_ERROR_DEBUG_NORET(retcode, msg)
#else
#define GMX_ERROR_DEBUG(retcode, msg) \
    do { \
        int _rc_internal = (retcode); \
        ::gmx::fatalError(_rc_internal, msg, __FILE__, __LINE__); \
        return _rc_internal; \
    } while (0)
#define GMX_ERROR_DEBUG_NORET(retcode, msg) \
        ::gmx::fatalError(retcode, msg, __FILE__, __LINE__)
#endif

/*!\}*/

} // namespace gmx

#endif
