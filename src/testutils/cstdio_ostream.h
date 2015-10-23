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

     Current GROMACS policy is using C-style output via (f)printf
     to avoid interference of output via (f)printf and std::cout/cerr.
     The StreamAdapter class enables the rest of the code to use
     C++ ostream syntax, while the output can internally be redirected
     to different sinks using either of the following:

     - C-style handles stdout, stderr
     - C++ ostreams std::cout, std::cerr
     - C++ output file streams
     - no effect stream for silencing

    \author
     R. Thomas Ullmann <tullman@gwdg.de>

    \date Mar 2015

    \copyright
     GROMACS licence

    \ingroup module_testutils
 */

#ifndef GMX_TESTUTILS_CSTDIO_OSTREAM_H
#define GMX_TESTUTILS_CSTDIO_OSTREAM_H

#include <cstdio>

#include <fstream>
#include <iomanip>
#include <ios>
#include <ostream>
#include <vector>

namespace gmx
{

namespace test
{

//! \cond internal

//! \internal
//! \brief pointer to a no-effect ostream needed for
//! a) the default constructor of class StreamAdapter and
//! b) for redirecting unwanted output as needed
extern std::ostream* cnullpt;

//! \internal
//! \class StreamAdapter
//!
//! \brief adapter class to redirect output to different streams
//!
//! Enables stream output without using std::cout/cerr. If wanted,
//! output can easily be redirected to std::cout/cerr also.
//! Employed exclusively within unit tests with Google test for
//! debug/error output, which can be redirected either to the
//! no-effect ostream or to a custom ostream that internally uses
//! fprintf for screen or file ouput and thus avoids interference
//! with (f)printf. Debugging information can be easily turned on
//! and off without reverting to preprocessor macros.
//!
//! \author R. Thomas Ullmann <tullman@gwdg.de>
//!
//! \copyright GROMACS license
//!
//! \date Mar 2015
class StreamAdapter : public std::ostream
{
    public:
        //! default constructor
        StreamAdapter() : std::ostream(), ost_(&cnullpt) { }
        /*! \brief constructor with target stream

            \param[in]    ost_ptr   pointer to a pointer to an ostream used for stream output
         */
        explicit StreamAdapter(std::ostream **ost_ptr) : std::ostream(), ost_(ost_ptr) { }
        //! default destructor
        ~StreamAdapter() { }

        /*! \brief set the target for stream  output

            \param[in]    onew   pointer to a pointer to an ostream used for stream output
         */
        void setStream(std::ostream** onew)
        {
            ost_ = onew;
        }
        //! return a handle to the target ostream
        std::ostream &getStream()
        {
            return **ost_;
        }

        //! conversion operator for beeing able to determine whether a stream adapter has an effect/leads to actual output
        operator bool () const { return *ost_ != cnullpt; }

        //! comparison operator < necessary for sorting (only here for completeness, likely not needed)
        //!
        //! \param[in]   a   StreamAdapter on the right-hand side of the (in)equality
        bool operator<  (const StreamAdapter &a) const { return *ost_ < *a.ost_; }
        //! comparison operator > necessary for sorting (only here for completeness, likely not needed)
        //!
        //! \param[in]   a   StreamAdapter on the right-hand side of the (in)equality
        bool operator>  (const StreamAdapter &a) const { return *ost_ > *a.ost_; }
        //! comparison operator for beeing able to determine whether two stream adapters have the same target stream
        //!
        //! \param[in]   a   StreamAdapter on the right-hand side of the (in)equality
        bool operator== (const StreamAdapter &a) const { return *ost_ == *a.ost_; }
        //! comparison operator for beeing able to determine whether two stream adapters have distinct target streams
        //!
        //! \param[in]   a   StreamAdapter on the right-hand side of the (in)equality
        bool operator!= (const StreamAdapter &a) const { return *ost_ != *a.ost_; }
    private:
        //! the target ostream
        std::ostream** ost_;
};

//! \internal
//! \ingroup module_testutils
//! \brief    ostream operator forwarding function for inserting an object
//!           into an ostream encapsulated by the StreamAdapter class
//!
//! \tparam   T   object type to be inserted into the ostream
//!
//! \param[in]   str   the insertion target for obj
//! \param[in]   obj   the object to be inserted into ost
template<typename T>
std::ostream &operator<<(StreamAdapter &str, const T &obj)
{
    str.getStream() << obj;
    return str.getStream();
}


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
//! \date Mar 2015
//!
//! \tparam   SIZE   buffer size, that is, number of characters that can be stored/buffered
template<std::size_t SIZE>
class CstdioStreamBuffer : public std::streambuf
{
    public:

        //! adapt the member function setp from the base class
        using std::streambuf::setp;
        //! adapt the data type used for storing characters from the base class
        using std::streambuf::char_type;
        //! adapt the data type used for representing integers from the base class
        using std::streambuf::int_type;
        //! adapt the traits_type from teh base class
        using std::streambuf::traits_type;

        //! default constructor
        CstdioStreamBuffer() : fp_(nullptr), buffer_(SIZE + 2, '\0')
        {
            // put area pointers to work with buffer_
            setAreaPointers();
        }

        //! \brief constructor with custom FILE handle
        //! \param[in]    fpin   FILE handle to be used as target of the stream output
        explicit CstdioStreamBuffer(FILE *fpin) : fp_(fpin), buffer_(SIZE + 2, '\0')
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
            //  setp(pbase, epptr)
            //  pbase points to the first buffer array element, epptr one past the last writable array element
            //  the extra two elements buffer_[SIZE], buffer_[SIZE+1] at the end will be used in this variant
            //  to store the last character that was inserted into the stream and led to overflow() and the
            //  terminating '\0' for C-style fprintf
            setp(&(buffer_.data()[0]), &(buffer_.data()[SIZE]));
        }

        //! write buffer contents to file handle fp_
        bool printBuffer()
        {
            // upon success, fprintf returns the number of written characters
            // otherwise it returns a result < 0
            return (fprintf(fp_, "%s", (char*)buffer_.data()) >= (int)0 ? true : false);
        }

        //! custom overload for the overflow function
        virtual int_type overflow(int_type ch)
        {
            if (fp_ != nullptr)
            {
                // write the last inserted character and the terminator
                // '\0' to the array elements pointed to by pptr and pptr+1.
                // This is always possible since we reserved two array
                // elements at the end of the buffer.
                *pptr() = ch;
                pbump((std::vector<char_type>::size_type) 1);
                *pptr() = '\0';
                pbump((std::vector<char_type>::size_type) 1);
                // write data and clear buffer array
                if (flush())
                {
                    return ch;
                }
            }
            // should never get here
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
            for (std::vector<char_type>::size_type i = 0; i <= (std::vector<char_type>::size_type)(SIZE + 1); ++i)
            {
                buffer_[i] = '\0';
            }
            // reset the the put pointer to pbase
            const std::ptrdiff_t n = pptr() - pbase();
            pbump(-n);

            return result;
        }

        //! file handle as target for writing the contents of the buffer
        FILE*                  fp_;
        //! \brief the char buffer itself
        //!
        //! the extra two elements buffer[SIZE], buffer[SIZE+1] at the end are needed to store
        //! the last character that was inserted into the stream and led to overflow() and the
        //! terminating '\0' for C-style fprintf
        std::vector<char_type> buffer_;
};

//! \endcond

}      // end namespace test

}      // end namespace gmx

#endif // end header
