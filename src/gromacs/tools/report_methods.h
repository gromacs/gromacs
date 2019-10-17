/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
#ifndef GMX_TOOLS_REPORT_METHODS_H
#define GMX_TOOLS_REPORT_METHODS_H

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/textwriter.h"

namespace gmx
{

class ReportMethodsInfo
{
public:
    static const char                       name[];
    static const char                       shortDescription[];
    static ICommandLineOptionsModulePointer create();
};

// Helper functions of the class

/*! \brief
 * Write appropiate Header to output stream.
 *
 * \param[in] writer TextWriter object for writing information.
 * \param[in] text String with the header before writing.
 * \param[in] section String with section text for header.
 * \param[in] writeFormattedText If we need to format the text for LaTeX output or not
 */
void writeHeader(TextWriter* writer, const std::string& text, const std::string& section, bool writeFormattedText);

/*! \brief
 * Write information about the molecules in the system.
 *
 * This method should write all possible information about
 * the molecular composition of the system.
 *
 * \param[in] writer TextWriter object for writing information.
 * \param[in] top Local topology used to derive the information to write out.
 * \param[in] writeFormattedText Decide if we want formatted text output or not.
 *
 */
void writeSystemInformation(TextWriter* writer, const gmx_mtop_t& top, bool writeFormattedText);

/*! \brief
 * Write information about system parameters.
 *
 * This method writes the basic information for the system parameters
 * and simulation settings as reported in the \p ir.
 *
 * \param[in] writer TextWriter object for writing information.
 * \param[in] ir Reference to inputrec of the run input.
 * \param[in] writeFormattedText Decide if we want formatted text output or not.
 */
void writeParameterInformation(TextWriter* writer, const t_inputrec& ir, bool writeFormattedText);

/*! \brief
 * Wrapper for writing out information.
 *
 * This function is actual called from within the run method
 * to write the information to the terminal or to file.
 * New write out methods should be added to it instead of adding them in run.
 *
 * \param[in] outputStream The filestream used to write the information to.
 * \param[in] ir Reference to inputrec of the run input.
 * \param[in] top Local topology used to derive the information to write out.
 * \param[in] writeFormattedText Decide if we want formatted text output or not.
 * \param[in] notStdout Bool to see if we can close the file after writing or not in case of stdout.
 */
void writeInformation(TextOutputFile*   outputStream,
                      const t_inputrec& ir,
                      const gmx_mtop_t& top,
                      bool              writeFormattedText,
                      bool              notStdout);

} // namespace gmx

#endif
