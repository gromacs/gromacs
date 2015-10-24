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
/*! \internal
    \file

    \brief
     ostream-related classes that internally use cstdio

     The CstdioStreamBuffer class allows ostreams to set up ostreams to
     internally use fprintf() from cstdio

     Class CstdioOstream behaves like a std::ostream with the addition of
     a bool operator for testing availability of the stream for output.

    \author
     R. Thomas Ullmann <tullman@gwdg.de>

    \date Mar 2015

    \copyright
     GROMACS licence

    \ingroup module_testutils
 */

#ifndef GMX_TESTUTILS_CSTDIO_STREAMBUF_H
#define GMX_TESTUTILS_CSTDIO_STREAMBUF_H

#include <cstdio>

#include <algorithm>
#include <ostream>
#include <streambuf>
#include <vector>

namespace gmx
{

namespace test
{

//! \cond internal

//! \internal
//! \class CstdioStreamBuffer
//!
//! \brief Buffer for std::ostream that internally uses C style fprintf
//!
//! Buffer for std::ostream that internally uses C style fprintf for
//! compliance with the current GROMACS convention of not using std::cout/cerr.
//! Can also be used as substitute for the buffer in std::cout/cerr
//! to make them use fprintf internally.
//!
//! \author R. Thomas Ullmann <tullman@gwdg.de>
//!
//! \copyright GROMACS license
//!
//! \date Oct 2015
class CstdioStreamBuffer : public std::streambuf
{
    public:

        //! default constructor, sets up a no-effect streambuf
        CstdioStreamBuffer() : fp_(nullptr), bufsize_(0), buffer_(bufsize_ + 2, '\0')
        {
            // put area pointers to work with buffer_
            setAreaPointers();
        }

        //! \brief constructor with custom FILE handle
        //! \param[in]   size      number of characters the buffer will be able to store
        //! \param[in]   fpin      FILE handle to be used as target of the stream output
        explicit CstdioStreamBuffer(size_t size, FILE *fpin) : fp_(fpin), bufsize_(size), buffer_(bufsize_ + 2, '\0')
        {
            // put area pointers to work with buffer_
            setAreaPointers();
        }

        //! the destructor, flush buffer to make sure that all information was printed/written
        ~CstdioStreamBuffer()
        {
            flush();
        }

    protected:

        //! \brief put area pointers to work with buffer
        void setAreaPointers()
        {
            // setp(pbase, epptr)
            // pbase points to the first buffer array element, epptr one past the last writable array element
            // the extra two elements buffer_[bufsize_], buffer_[bufsize_+1] at the end will be used in this
            // variant to store the last character that was inserted into the stream and led to overflow()
            // and the terminating '\0' for C-style fprintf
            setp(&(buffer_.data()[0]), &(buffer_.data()[bufsize_]));
        }

        //! write buffer contents to file handle fp_
        bool printBuffer()
        {
            if (fp_ != nullptr)
            {
                // upon success, fprintf returns the number of written characters
                // otherwise it returns a result < 0
                return (fprintf(fp_, "%s", (char*)buffer_.data()) >= (int)0 ? true : false);
            }
            else
            {
                return false;
            }
        }

        //! custom overload for the overflow function
        //! \param[in]   ch   last inserted character
        virtual int_type overflow(int_type ch)
        {
            if (fp_ != nullptr)
            {
                // write the last inserted character and the terminator
                // '\0' to the array elements pointed to by pptr and pptr+1.
                // This is always possible since we reserved two array
                // elements at the end of the buffer.
                // EOF is not written
                if (ch != traits_type::eof())
                {
                    *pptr() = ch;
                    pbump((std::vector<char_type>::size_type) 1);
                }
                *pptr() = '\0';
                pbump((std::vector<char_type>::size_type) 1);
                // write data and clear buffer array
                if (flush())
                {
                    return ch;
                }
            }
            // should never get here, unless fp_ == nullptr
            return traits_type::eof();
        }
        //! synchronize the contents in the buffer with those
        //! of the associated character sequence, that is, fp_
        virtual int sync()
        {
            return flush() ? 0 : -1;
        }
    private:
        //! let StreamAdapter access private members
        friend class StreamAdapter;

        //! output the contents of buffer to the associated character
        //! sequence and empty buffer
        bool flush()
        {
            // output buffer contents to the destination FILE handle
            const bool result = printBuffer();

            // clear the buffer, irregardless of the success of printBuffer()
            std::fill(buffer_.begin(), buffer_.end(), '\0');
            // reset the the put pointer to pbase
            const std::ptrdiff_t n = pptr() - pbase();
            pbump(-n);

            return result;
        }

        //! file handle as target for writing the contents of the buffer
        FILE* const            fp_;
        //! number of characters that can be stored in the buffer
        const size_t           bufsize_;
        //! \brief the char buffer itself
        //!
        //! the extra two elements buffer[bufsize_], buffer[buf_size+1] at the end
        //! are needed to store the last character that was inserted into the stream
        //! and led to overflow() and the terminating '\0' for C-style fprintf
        std::vector<char_type> buffer_;
};

//! \internal
//! \class TestableOstream
//!
//! \brief like std::ostream, with the addition of a bool operator
//!        for testing availability of the stream
//!
//!        TestableOstream debug(pointerToStreamBuffer);
//!        if (debug) { doSomething(); }
//!
//! \author R. Thomas Ullmann <tullman@gwdg.de>
//!
//! \copyright GROMACS license
//!
//! \date Oct 2015
//!
//! \tparam   SIZE   buffer size, that is, number of characters that can be stored/buffered
class TestableOstream : public std::ostream
{
    public:
        //! default constructor
        TestableOstream()
            : std::ostream(nullptr)
        { }

        //! constructor with streambuf
        //! \param[in]   bufptr   pointer to the output stream buffer
        explicit TestableOstream(std::streambuf *bufptr)
            : std::ostream(bufptr)
        { }

        //! default destructor
        ~TestableOstream()
        { }

        //! bool operator  meant for testing whether the stream is active/has an attached streambuf
        operator bool() const
        {
            if (std::ostream::rdbuf() == nullptr)
            {
                return false;
            }
            else
            {
                return std::ostream::good();
            }
        }
};

//! \endcond

}      // end namespace test

}      // end namespace gmx

#endif // end header
