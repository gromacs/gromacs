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
     The implementation file for the verbosity management

    \author
     R. Thomas Ullmann <tullman@gwdg.de>

    \date Mar 2015

    \copyright
     GROMACS licence

    \libinternal \file
    \ingroup module_utility
 */

#include "gmxpre.h"

#include "verbosity_management.h"

// smart pointer with exclusive ownership of the managed object
#include "external/boost/boost/scoped_ptr.hpp"
// gromacs exceptions
#include "gromacs/utility/exceptions.h"

namespace gmx
{

namespace verbosity_management
{

//! This unnamed namespace is meant to clarify that the enclosed objects and
//! pointers should not be used directly in the rest of the code.
namespace
{

//! enables redirecting regular output to a file
boost::scoped_ptr<std::ofstream> module_outfstr_pt;
//! enables redirecting error output to a file
boost::scoped_ptr<std::ofstream> module_errfstr_pt;

#ifdef GMX_OSTREAMS_USE_CSTDIO
//! flushed whenever std::endl is inserted into the stream
CstdioStreamBuffer<20>  module_errbuf(stderr);
//! flushed whenever std::endl is inserted into the stream
CstdioStreamBuffer<255> module_outbuf(stdout);
std::ostream            module_errstr(&module_errbuf);
std::ostream            module_outstr(&module_outbuf);
//! output to stderr
std::ostream* const     cerrpt = &module_errstr;
//! output to stdout
std::ostream* const     coutpt = &module_outstr;
#else
//! output to stderr
std::ostream* const cerrpt = &std::cerr;
//! output to stdout
std::ostream* const coutpt = &std::cout;
#endif

//! targets of regular and error output in general (plus cnullpt outside the unnamed namespace, below)
//! a pointer to an ostream for redirecting error output as needed
std::ostream* module_errpt = cerrpt;
//! a pointer to an ostream for redirecting regular output as needed
std::ostream* module_outpt = coutpt;
}

//! the StreamAdapters are the entities to be used like std::ostreams in the rest of the code
//! for example: debug << "This is a debugging note" << endl;

//! all module parts access the debug output stream through this object
StreamAdapter debug(&cnullpt);
//! all module parts access the debug output stream through this object
StreamAdapter error(&module_errpt);
//! all module parts access the high   verbosity output stream through this object
StreamAdapter blab3(&cnullpt);
//! all module parts access the medium verbosity output stream through this object
StreamAdapter blab2(&cnullpt);
//! all module parts access the low    verbosity output stream through this object
StreamAdapter blab1(&cnullpt);
//! all module parts access the no     verbosity output stream through this object
StreamAdapter blab0(&cnullpt);


// ===========================================================================================
// free functions for settings that affect all members of the gmx::verbosity_management namespace ...
// ===========================================================================================

//! set the verbosity level -- how important a piece of information has to be to be printed/written to file
void setBlabLevel(const unsigned int blab)
{
    if      (blab >= 3)
    {
        // output every detail on the simulation available except for purely technical stuff like which function was entered
        blab3.setStream(&module_outpt);
        blab2.setStream(&module_outpt);
        blab1.setStream(&module_outpt);
    }
    else if (blab == 2)
    {
        // output more information
        blab3.setStream(&cnullpt);
        blab2.setStream(&module_outpt);
        blab1.setStream(&module_outpt);
    }
    else if (blab == 1)
    {
        // output some additional information
        blab3.setStream(&cnullpt);
        blab2.setStream(&cnullpt);
        blab1.setStream(&module_outpt);
    }
    else
    {
        // otherwise all output stays redirect to the no-effect ostream
        // and the computation proceeds silently
        blab3.setStream(&cnullpt);
        blab2.setStream(&cnullpt);
        blab1.setStream(&cnullpt);
    }
    // essential output is never silenced
    blab0.setStream(&module_outpt);
}

//! enable or disable debug output
void setDebugOut(const bool flag)
{
    if (flag)
    {
        // ouput to stdout or a file, by whichever means
        debug.setStream(&module_outpt);
    }
    else
    {
        // no ouput
        debug.setStream(&cnullpt);
    }
}

//! enable debug output
void setDebugOut()
{
    setDebugOut(true);
}

//! enable or disable error output
void setErrorOut(const bool flag)
{
    if (flag)
    {
        // ouput to stderr or a file, by whichever means
        error.setStream(&module_errpt);
    }
    else
    {
        // no ouput
        error.setStream(&cnullpt);
    }
}

//! enable error output
void setErrorOut()
{
    setErrorOut(true);
}

//! \brief set the verbosity level -- how important a piece of information has to be to reach the user
void setVerbosity(const std::string s)
{
    if      (s == "debug"  || s == "dbg")
    {
        // high/blab3 below plus information on the program flow, e.g., which function was entered under which conditions etc.
        setBlabLevel(3);
        setDebugOut();
    }
    else if (s == "high"   || s == "blab3")
    {
        // print every detail on the simulation available except for purely technical stuff like which function was entered
        setBlabLevel(3);
    }
    else if (s == "medium" || s == "blab2")
    {
        // print some more progress information, e.g., on the progress of substeps
        setBlabLevel(2);
    }
    else if (s == "low"    || s == "blab1")
    {
        // print some additional information, e.g., on the progress of multi-step calculations
        setBlabLevel(1);
    }
    else
    {
        // only print the final results, otherwise keep silence
        setBlabLevel(0);
    }
}

//! initialize the lambda site module with custom settings
void initModule(const bool print_dbg, const bool print_err, const std::string &verbosity)
{
    //! activate debug output if requested and set the requested level of verbosity
    if (print_dbg)
    {
        setDebugOut();
    }
    //! activate debug output if requested and set the requested level of verbosity
    if (print_err)
    {
        setErrorOut();
    }
    setVerbosity(verbosity);
}


//! initialize the lambda site module with default settings
void setDefaults()
{
    // defaults
    // no debug output
    const bool        print_dbg = false;
    // print error output
    const bool        print_err = true;
    // do not print any extra information
    const std::string verbosity("low");
    initModule(print_dbg, print_err, verbosity);
}

//! redirect module output normally written to stdout to a file or reverse
void setRedirectStdOutToFile(const bool flag, const std::string &filename)
{
    if (flag)
    {
        boost::scoped_ptr<std::ofstream> tmp_ptr(new std::ofstream(filename.c_str(), std::ofstream::app));
        if (!tmp_ptr->good())
        {
            std::fprintf(stderr, "Error in gmx::verbosity_management::setRedirectStdOutToFile(): failed to open file \"%s\".", filename.c_str());
            throw InternalError("Failed to open file for writing");
        }
        // the previously managed ofstream will be deleted by tmp_ptr's destructor when going out of scope
        module_outfstr_pt.swap(tmp_ptr);
        module_outpt = module_outfstr_pt.get();
    }
    else
    {
        module_outfstr_pt.reset();
        module_outpt = coutpt;
    }

}
void setRedirectStdOutToFile(const std::string &filename) { return setRedirectStdOutToFile(true, filename); }
void setRedirectStdOutToFile(const bool flag) { return setRedirectStdOutToFile(flag, std::string("LambdaSiteOut.txt")); }

//! redirect module output normally written to stderr to a file or reverse
void setRedirectErrOutToFile(const bool flag, const std::string &filename)
{
    if (flag)
    {
        boost::scoped_ptr<std::ofstream> tmp_ptr(new std::ofstream(filename.c_str(), std::ofstream::app));
        if (!tmp_ptr->good())
        {
            std::fprintf(stderr, "Error in gmx::verbosity_management::setRedirectStdErrToFile(): failed to open file \"%s\".", filename.c_str());
            throw InternalError("Failed to open file for writing");
        }
        // the previously managed ofstream will be deleted by tmp_ptr's destructor when going out of scope
        module_errfstr_pt.swap(tmp_ptr);
        module_errpt = module_errfstr_pt.get();
    }
    else
    {
        module_errfstr_pt.reset();
        module_errpt = cerrpt;
    }

}

void setRedirectErrOutToFile(const std::string &filename) { return setRedirectErrOutToFile(true, filename); }
void setRedirectErrOutToFile(const bool flag) { return setRedirectErrOutToFile(flag, std::string("LambdaSiteErr.txt")); }

// ===========================================================================================
// ... free functions for settings that affect all members of the gmx::verbosity_management namespace
// ===========================================================================================


} // end namespace verbosity_management

} // end namespace gmx
