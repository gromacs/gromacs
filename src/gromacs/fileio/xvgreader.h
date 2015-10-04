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
/*! \internal \file
 * \brief
 * Defines data structures for containing data read from an .xvg file,
 * and functions for doing reading.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#ifndef GMX_FILEIO_XVGREADER_H
#define GMX_FILEIO_XVGREADER_H

#include <string>
#include <vector>

#include "gromacs/utility/textreader.h"

namespace gmx
{

/*! \internal
 * \brief Contains data corresponding to a column of an .xvg file */
struct XvgColumn
{
    /*! \brief The contents of the .xvg legend field (or an empty string) */
    std::string         legend_;
    /*! \brief Column of data */
    std::vector<double> rows_;
};

/*! \internal
 * \brief Contains data corresponding to a column of an .xvg graph file */
struct XvgTable
{
    //! The contents of the .xvg title field (or an empty string)
    std::string title_;
    //! The contents of the .xvg subtitle field (or an empty string)
    std::string subtitle_;
    /*! \brief The label of the "y axis" (or an empty string)
     *
     * Note that the "x axis" is labelled in the ordinate_ member */
    std::string yAxisLabel_;
    //! The "x axis" data, including label
    XvgColumn   ordinate_;
    /*! \brief Table of columns of "y axis" abscissa data, in column-major format.
     *
     * All columns have equal numbers of rows. */
    std::vector<XvgColumn> columns_;
};

/*! \brief Read the named .xvg file to return the data it contains
 *
 * \param[in]  filename  Name of .xvg filename to read
 *
 * \see readXvgTable(const char *) for requirements on the contents of
 * the file
 *
 * Does not throw. Calls readXvgTable(const char*), but exits with a
 * fatal error when any exception is thrown.
 *
 * No new GROMACS code should call this function.
 *
 * \todo Replace calls to at least read_xvg() and read_xvg_legend()
 * with calls to this function or readXvgTable() (and probably fix
 * associated bugs). */
XvgTable readXvgTableFromLegacyCode(const char *filename);
/*! \brief Read the named .xvg file to return the data it contains.
 *
 * \param[in]  filename  Name of .xvg filename to read
 *
 * All data from the .xvg file is returned, where the first column is
 * assumed to be ordinate values. Data is required to be numerical,
 * and to exactly fill a regular array (i.e same number of rows for
 * each column and vice versa). Fields can be separated, preceded or
 * trailed by any amount or kind of whitespace. Completely blank lines
 * are not permitted. It is permitted for no legend strings to be
 * present, but if any are present, then no legend string may be
 * duplicated. Axis labels, and graph title and subtitle are read, if
 * present.
 *
 * \throws std::bad_alloc     when out of memeory
 * \throws InvalidInputError  when the file cannot be parsed correctly
 * \throws FileIOError        when an I/O error occurs (containing errno diagnostics)
 */
XvgTable readXvgTable(const char *filename);
/*! \brief Read a stream containing an .xvg file to return the data it
 * contains
 *
 * \param[in]  reader    Valid object for reading lines
 *
 * \see readXvgTable(const char *) for requirements on the contents of
 * the stream
 *
 * \throws std::bad_alloc     when out of memeory
 * \throws InvalidInputError  when the file cannot be parsed correctly
 * \throws FileIOError        when an I/O error occurs (containing errno diagnostics)
 *
 * This function might be called when an open stream is already
 * available, e.g. during testing, when a StringInputStream can be
 * used to configure a TextReader. */
XvgTable readXvgTable(TextReader *reader);

} // namespace

#endif
