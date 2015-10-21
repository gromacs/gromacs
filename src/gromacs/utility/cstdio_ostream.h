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
     ostream classes that internally use cstdio

     Current GROMACS policy is using C-style output via (f)printf
     to avoid interference of output via (f)printf and std::cout/cerr.
     The StreamAdapter class enables the rest of the code to use
     C++ ostream syntax, while the output can internally be redirected
     to different sinks using either of the following:

     -C-style handles stdout, stderr
     -C++ ostreams std::cout, std::cerr
     -C++ output file streams
     -no effect stream for silencing

    \author
     R. Thomas Ullmann <tullman@gwdg.de>

    \date Mar 2015

    \copyright
     GROMACS licence

    \libinternal \file
    \ingroup module_utility
 */

#ifndef GMX_UTILITY_ERRORCODES_H
#define GMX_UTILITY_ERRORCODES_H

#include <fstream>
#include <ios>
#include <iomanip>
#include <ostream>
#include <vector>

namespace gmx
{

// enclose cnullpt in the verbosity_management subnamespace to avoid potential name clashes with outside code
namespace verbosity_management
{

//! pointer to a no-effect ostream needed for
//! a) the default constructor of class StreamAdapter and
//! b) for redirecting unwanted output as needed
extern std::ostream* cnullpt;

}

/*! \class StreamAdapter

    \brief adapter class to redirect output to different streams

    default: no-effect ostream
    within GROMACS: debug/error output and output for different verbosity
                    levels redirected either to the no-effect ostream
                    or to a custom ostream that internally uses fprintf
                    for screen or file ouput.
                    In this way, the convenience of stream output for
                    whole objects can still be used without using
                    std::cout/cerr. Output can easily be redirected to
                    std::cout/cerr again if desired by setting the ostream
                    pointer accordingly.

    \author R. Thomas Ullmann <tullman@gwdg.de>

    \copyright GROMACS license

    \date Mar 2015

    \ingroup module_utility
    \inlibraryapi
 */
class StreamAdapter : public std::ostream
{
    public:
        //! default constructor
        StreamAdapter() : ost_(&verbosity_management::cnullpt)
        {

        }
        /*! \brief constructor with target stream

            \param[in]    ost_ptr   pointer to a pointer to an ostream used for stream output
         */
        explicit StreamAdapter(std::ostream **ost_ptr) : ost_(ost_ptr)
        {

        }
        //! default destructor
        ~StreamAdapter()
        {

        }
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
        operator bool () const { return *ost_ != verbosity_management::cnullpt; }

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

//! ostream operator forwarding function for inserting an object into
//! an ostream encapsulated by the StreamAdapter class
//!
//! \tparam   T   object type to be inserted into the ostream
//!
//! \param[in]   str   the insertion target for obj
//! \param[in]   obj   the object to be inserted into ost
template<typename T>
inline std::ostream &operator<<(StreamAdapter &str, const T &obj)
{
    str.getStream() << obj;
    return str.getStream();
}

/*! \class CstdioStreamBuffer

    \brief Buffer for std::ostream that internally uses C style fprintf

    Buffer for std::ostream that internally uses C style fprintf for
    compliance with the current GROMACS convention of not using std::cout/cerr
    while retaining the clarity of the ostream notation for output of
    entire objects.

    \author R. Thomas Ullmann <tullman@gwdg.de>

    \copyright GROMACS license

    \date Mar 2015

    \tparam   SIZE   buffer size, that is, number of characters that can be stored/buffered

    \ingroup module_utility
    \inlibraryapi
 */
template<std::size_t SIZE>
class CstdioStreamBuffer : public std::streambuf
{
    public:

        //! use the member functions of the base class
        using Base = std::streambuf;
        //! adapt the data type used for storing characters from the base class
        using char_type = typename Base::char_type;
        //! adapt the data type used for representing integers from the base class
        using int_type = typename Base::int_type;

        //! default constructor
        CstdioStreamBuffer() : fp_(nullptr), buffer_(SIZE + 2, '\0') // value-initialize buffer
        {
            // put area pointers to work with "buffer"
            setAreaPointers();
        }
        //! \brief constructor with custom FILE handle
        //! \param[in]    fpin   FILE handle to be used as target of the stream output
        explicit CstdioStreamBuffer(FILE *fpin) : fp_(fpin), buffer_(SIZE + 2, '\0') // value-initialize buffer
        {
            // put area pointers to work with "buffer"
            setAreaPointers();
        }
        //! the destructor, flush buffer to make sure that all information was printed/written
        ~CstdioStreamBuffer()
        {
            flush();
        }

    protected:
        //! \brief put area pointers to work with "buffer"
        void setAreaPointers()
        {
            //  setp(pbase, epptr)
            //  pbase points to the first buffer array element, epptr one past the last writable array element
            //  the extra two elements buffer[SIZE], buffer[SIZE+1] at the end will be used in this variant
            //  to store the last character that was inserted into the stream and led to overflow() and the
            //  terminating '\0' for C-style fprintf
            Base::setp(&(buffer_.data()[0]), &(buffer_.data()[SIZE])); // set std::basic_streambuf
        }

        //! write buffer contents to file handle fp
        bool printBuffer()
        {
            // upon success, fprintf returns the number of written characters
            // otherwise it returns a result < 0
            return (fprintf(fp_, "%s", buffer_.data()) >= (int)0 ? true : false);
        }

        //! custom overload for the overflow function
        virtual int_type overflow(int_type ch)
        {
            if (fp_ != nullptr)
            {
                // write the last inserted character and the terminator
                // '\0' to the array elements pointed to by pptr and pptr+1.
                // This is always, possible since we reserved two array
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
            //! should never get here
            return traits_type::eof();
        }
        //! synchronize the contents in the buffer with those
        //! of the associated character sequence, that is, fp_
        virtual int sync()
        {
            return flush() ? 0 : -1;
        }
    private:
        friend class StreamAdapter;
        //! output the contents of buffer to the associated character
        //! sequence and empty buffer
        bool flush()
        {
            //! output buffer contents to the destination FILE handle
            const bool result = printBuffer();

            //! clear the buffer, irregardless of the success of printBuffer()
            for (std::vector<char_type>::size_type i = 0; i <= (std::vector<char_type>::size_type)(SIZE + 1); ++i)
            {
                buffer_[i] = '\0';
            }
            //! reset the the put pointer to pbase
            const std::ptrdiff_t n = pptr() - pbase();
            pbump(-n);

            return result;
        }

        //! file handle as target for writing the contents of the buffer
        FILE*                   fp_;
        //! \brief the char buffer itself
        //!
        //! the extra two elements buffer[SIZE], buffer[SIZE+1] at the end are needed to store
        //! the last character that was inserted into the stream and led to overflow() and the
        //! terminating '\0' for C-style fprintf
        std::vector<char_type> buffer_;
};


}      // end namespace gmx

#endif // end header
