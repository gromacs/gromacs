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
 * Declares functions for file handling.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_FILE_H
#define GMX_UTILITY_FILE_H

#include <cstdio>

#include <string>

#include "common.h"

namespace gmx
{

/*! \brief
 * Basic file object.
 *
 * This class provides basic file I/O functionality and uses exceptions
 * (FileIOError) for error reporting.
 *
 * \inpublicapi
 * \ingroup module_utility
 */
class File
{
    public:
        /*! \brief
         * Creates a file object and opens a file.
         *
         * \param[in] filename  Path of the file to open.
         * \param[in] mode      Mode to open the file in (for fopen()).
         * \throws    FileIOError on any I/O error.
         *
         * \see open(const char *, const char *)
         */
        File(const char *filename, const char *mode);
        //! \copydoc File(const char *, const char *)
        File(const std::string &filename, const char *mode);
        /*! \brief
         * Destroys the file object.
         *
         * If the file is still open, it is closed.
         * Any error conditions will be ignored.
         */
        ~File();

        /*! \brief
         * Opens a file.
         *
         * \param[in] filename  Path of the file to open.
         * \param[in] mode      Mode to open the file in (for fopen()).
         * \throws    FileIOError on any I/O error.
         *
         * The file object must not be open.
         */
        void open(const char *filename, const char *mode);
        //! \copydoc open(const char *, const char *)
        void open(const std::string &filename, const char *mode);
        /*! \brief
         * Closes the file object.
         *
         * \throws  FileIOError on any I/O error.
         *
         * The file must be open.
         */
        void close();

        /*! \brief
         * Returns a file handle for interfacing with C functions.
         *
         * The file must be open.
         * Does not throw.
         */
        FILE *handle();

        /*! \brief
         * Reads contents of a file to a std::string.
         *
         * \param[in] filename  File to read.
         * \returns   The contents of \p filename.
         * \throws    std::bad_alloc if out of memory.
         * \throws    FileIOError on any I/O error.
         */
        static std::string readToString(const char *filename);
        //! \copydoc readToString(const char *)
        static std::string readToString(const std::string &filename);

    private:
        FILE                   *fp_;

        GMX_DISALLOW_COPY_AND_ASSIGN(File);
};

} // namespace gmx

#endif
