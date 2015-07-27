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
/*! \libinternal \file
 * \brief
 * Declares implementations for textstream.h interfaces for input/output to
 * in-memory strings.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_STRINGSTREAM_H
#define GMX_UTILITY_STRINGSTREAM_H

#include <string>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/textstream.h"

namespace gmx
{

/*! \libinternal \brief
 * Text output stream implementation for writing to an in-memory string.
 *
 * Implementations for the TextOutputStream methods throw std::bad_alloc if
 * reallocation of the string fails.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class StringOutputStream : public TextOutputStream
{
    public:
        //! Returns the text written to the stream so far.
        const std::string &toString() const { return str_; }

        // From TextOutputStream
        virtual void write(const char *text);
        virtual void close();

    private:
        std::string str_;
};

/*! \libinternal \brief
 * Helper class to convert static string data to a stream
 */
class StringInputStream : public TextInputStream
{
    public:
        /*! \brief
         * Constructor that initializes local variables including iterator.
         *
         * \param[in] input Pointer to vector of strings to be served by the stream.
         */
        StringInputStream(ConstArrayRef<const char *> input)
        {
            for (ConstArrayRef<const char *>::iterator i = input.begin(); (i < input.end()); ++i)
            {
                input_.append(*i);
                input_.append("\n");
            }
            pos_   = 0;
        }
        virtual ~StringInputStream() {}

        /*! \brief
         * Reads a line (with newline included) from the stream.
         *
         * \param[out] line    String to receive the line.
         * \returns    `false` if nothing was read because the stream ended.
         *
         * On error or when `false` is returned, \p line will be empty.
         */
        virtual bool readLine(std::string *line)
        {
            if (pos_ == input_.size())
            {
                return false;
            }
            else
            {
                size_t newpos = input_.find("\n", pos_);
                if (newpos == std::string::npos)
                {
                    newpos = input_.size();
                }
                else
                {
                    // To include the newline as well!
                    newpos += 1;
                }
                line->assign(input_.substr(pos_, newpos-pos_));
                pos_ = newpos;
                return true;
            }
        }
        /*! \brief
         * Closes the stream (does nothing in this case).
         */
        virtual void close() {};
    private:
        std::string input_;
        size_t      pos_;
};

} // namespace gmx

#endif
