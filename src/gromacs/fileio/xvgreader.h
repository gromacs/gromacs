/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014,2015, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Defines data structures for containing data read from an .xvg file,
 * and functions for doing reading.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 */
#ifndef GMX_FILEIO_XVGREADER_H
#define GMX_FILEIO_XVGREADER_H

#include <string>
#include <vector>

#include "gromacs/utility/textstream.h"

namespace gmx
{

//! Contains data corresponding to a column of an .xvg file
struct XvgColumn
{
    /*! \brief The contents of the .xvg legend field (or an empty string) */
    std::string         legend_;
    /*! \brief Column of data */
    std::vector<double> rows_;
};

//! Contains data corresponding to the contents of an .xvg file
struct XvgTable
{
    /*! \brief Table of data, in column-major format.
     *
     * All columns have equal numbers of rows. */
    std::vector<XvgColumn> columns_;
    //! The contents of the .xvg subtitle field (or an empty string)
    std::string            subtitle_;
    //! The label of the x axis (or an empty string)
    std::string            xAxisLabel_;
    //! The label of the y axis (or an empty string)
    std::string            yAxisLabel_;
};

/*! \brief Read the named .xvg file to return the data it contains
 *
 * \param[in]  filename  Name of .xvg filename to read
 *
 * All data from the .xvg file is returned, e.g. we don't remove the
 * first column because it contains abscissa values, because we can't
 * assume the .xvg contains abscissa-ordinate columns. That's the
 * caller's problem.
 *
 * Does not throw. Calls readXvgTable(const char*), but exits with a
 * fatal error when any exception is thrown.
 *
 * No new GROMACS code should call this function.
 *
 * \todo Replace calls to at least read_xvg() and read_xvg_legend()
 * with calls to this function. */
XvgTable readXvgTableFromLegacyCode(const char *filename) noexcept;
/*! \brief Read the named .xvg file to return the data it contains.
 *
 * \param[in]  filename  Name of .xvg filename to read
 *
 * All data from the .xvg file is returned, e.g. we don't remove the
 * first column because it contains abscissa values, because we can't
 * assume the .xvg contains abscissa-ordinate columns. That's the
 * caller's problem.
 *
 * \throws std::bad_alloc     when out of memeory
 * \throws InvalidInputError  when the file cannot be parsed correctly
 * \throws FileIOError        when an I/O error occurs (containing errno diagnostics)
 */
XvgTable readXvgTable(const char *filename);
/*! \brief Read a stream containing an .xvg file to return the data it
 * contains
 *
 * \param[in]  stream    Open stream for reading lines (corresponds to \c filename)
 *
 * All data from the .xvg file is returned, e.g. we don't remove the
 * first column because it contains abscissa values, because we can't
 * assume the .xvg contains abscissa-ordinate columns. That's the
 * caller's problem.
 *
 * This function might be called when an open stream is already
 * available, e.g. when using an in-memory stream during testing. */
XvgTable readXvgTable(TextInputStream *stream);

} // namespace

#endif
