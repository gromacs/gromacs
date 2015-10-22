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

    \brief
     The header file for the verbosity management
     Manage the amount of information that reaches the user and
     how it reaches the user.

     Four verbosity levels adjust the amount of information that
     is provided (inspired by MEAD "blab" levels 0, 1, 2, 3).

     Current GROMACS policy is using C-style output via (f)printf.
     If this policy should change, undefining the preprocessor macro
     GMX_OSTREAMS_USE_CSTDIO should be enough to switch to using
     std::cout/cerr. The StreamAdapter class enables the rest of the
     code to use C++ ostream syntax, while the output can internally
     be redirected to different sinks using either of the following:

     -C++ ostreams std::cout, std::cerr
     -C-style handles stdout, stderr
     -C++ output file streams
     -no effect stream for silencing

    \author
     R. Thomas Ullmann <tullman@gwdg.de>

    \date Mar 2015

    \copyright
     GROMACS licence

    \inlibraryapi \file
    \ingroup module_verbosity_management

 */

#ifndef GMX_UTILITY_VERBOSITY_MANAGEMENT_H
#define GMX_UTILITY_VERBOSITY_MANAGEMENT_H

//! use ostreams that internally employ fprintf to stdout/stderr
#define GMX_OSTREAMS_USE_CSTDIO
#ifdef GMX_OSTREAMS_USE_CSTDIO
 #include <cstdio>
#endif

#include "cstdio_ostream.h"

namespace gmx
{

/*! \ingroup module_verbosity_management

    \brief
     Verbosity Management inspired by MEAD
     but avoiding potentially dangerous globally accessible macro definitions
     of dereferenced pointers to ostreams and the platform dependent /dev/null
 */
namespace verbosity_management
{

/* -------------------------------------------------------------------------------------

   Use these ostream handles like std::cout/cerr in the rest of the code, e.g.:
   debug << "Debugging Info" << endl;

   Internally, fprintf will be used for GROMACS. If the policy should change or for use
   elsewhere: Switching to std::cout/std::cerr is easily done by commenting out the
   definition of the preprocessor macro GMX_OSTREAMS_USE_CSTDIO above.
 */

//! \ingroup module_verbosity_management
//! \brief stream for error messages
extern StreamAdapter error;
//! \ingroup module_verbosity_management
//! \brief stream for debugging information
extern StreamAdapter debug;
//! \ingroup module_verbosity_management
//! \brief high detail stream, for highly detailed progress information/intermediate results
extern StreamAdapter blab3;
//! \ingroup module_verbosity_management
//! \brief medium detail stream, for moderatly detailed progress information/intermediate results
extern StreamAdapter blab2;
//! \ingroup module_verbosity_management
//! \brief low    detail stream, for brief progress information/intermediate results
extern StreamAdapter blab1;
//! \ingroup module_verbosity_management
//! \brief no detail stream, for essential information only
extern StreamAdapter blab0;

//! -------------------------------------------------------------------------------------

/*! \brief set the verbosity level, that is, how important a piece of information has to reach the user

      \param[in] blab  verbosity level 0 <= blab <= 3
                     >= 3:  print every detail on the simulation available except for
                            purely technical stuff like which function was entered
                        2:  print more additional information like intermediate results of a free
                            energy calculation with potentially many iteration cycles
                        1:  print some additional information like independent intermediate
                            contributions to a final result
                     <= 0:  only provide the final results, otherwise keep silence
 */
void setBlabLevel(const unsigned int blab);

/*! \brief  enable or disable debugging output

      \param[in] flag  true or false
                       true: provide information on the program flow, e.g., which function was entered
                             for which class under which conditions etc.
                       false: no debug output, blab level determines remaining verbosity
 */
void setDebugOut(const bool flag);

//! enable debug output
void setDebugOut();

/*! \brief set the verbosity level -- how important a piece of information has to be to reach the user

      \param[in] flag  true or false
                       true: provide information on the program flow, e.g., which function was entered
                             for which class under which conditions etc.
                       false: no debug output, blab level determines remaining verbosity
 */
void setErrorOut(const bool flag);

//! enable error output
void setErrorOut();

/*! accepts MEAD-like strings for the verbosity level and the debug flag blab3, blab2, blab1, blab0, debug

      \param[in] s   "debug":           high/blab3 below plus information on the program flow, e.g.,
                                        which function was entered under which conditions etc.
                     "blab3"/"high":    print every detail on the simulation available except for
                                        purely technical stuff like which function was entered
                     "blab2"/"medium":  print more additional information like intermediate results of a free
                                        energy calculation with potentially many iteration cycles
                     "blab1"/"low":     print some additional information like independent intermediate
                                        contributions to a final result
                     everything else:   only provide the final results, otherwise keep silence
 */
void setVerbosity(const std::string s);

//! initialize the verbosity management with default settings
void setDefaults();

//! ----------------------------------------------------------------------------------------
//! redirect module streams normally writing to stderr to file
//! ----------------------------------------------------------------------------------------

/*! \brief (un)trigger redirection of error output, normally printed to stderr, to file

      \param[in] flag  true or false
                       true : redirect output streams of this module that are normally
                              printed to stderr to file
                       false: redirect these ouput streams to stderr again
      \param[in] filename    filename for writing the output
 */
void setRedirectErrOutToFile(const bool flag, const std::string &filename);

/*! \brief (un)trigger redirection of error output, normally printed to stderr, to file

      \param[in] flag  true or false
                       true : redirect output streams of this module that are normally
                              printed to stderr to file
                       false: redirect these ouput streams to stderr again
 */
void setRedirectErrOutToFile(const bool flag);

/*! \brief (un)trigger redirection of error output, normally printed to stderr, to file

      \param[in] filename    filename for writing the output
 */
void setRedirectErrOutToFile(const std::string &filename);

//! ----------------------------------------------------------------------------------------
//! redirect module streams normally writing to stdout to file
//! ----------------------------------------------------------------------------------------

/*! \brief (un)trigger redirection of regular output, normally printed to stdout, to file

      \param[in] flag  true or false
                       true : redirect output streams of this module that are normally
                              printed to stdout to file
                       false: redirect these ouput streams to stdout again
      \param[in] filename    filename for writing the output
 */
void setRedirectStdOutToFile(const bool flag, const std::string &filename);

/*! \brief (un)trigger redirection of error output, normally printed to stdout, to file

      \param[in] flag  true or false
                       true : redirect output streams of this module that are normally
                              printed to stdout to file
                       false: redirect these ouput streams to stdout again
 */
void setRedirectStdOutToFile(const bool flag);

/*! \brief (un)trigger redirection of error output, normally printed to stdout, to file

      \param[in] filename    filename for writing the output
 */
void setRedirectStdOutToFile(const std::string &filename);

//! ----------------------------------------------------------------------------------------

/*! initialize the lambda site module with custom settings

      \param[in] print_dbg  flag triggering error output
      \param[in] print_err  flag triggering debug output
      \param[in] verbosity  string indicating the amount of verbosity
                            (possible choices documented for set_verbosity()
 */
void init(const bool print_dbg, const bool print_err, const std::string &verbosity);

// ---------------------------------------------------------------------------------------


}      // end namespace verbosity_management

}      // end namespace gmx

#endif // end header
