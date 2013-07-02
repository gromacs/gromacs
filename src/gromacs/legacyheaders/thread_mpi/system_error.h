/*
   This source code file is part of thread_mpi.
   Written by Sander Pronk, Erik Lindahl, and possibly others.

   Copyright (c) 2009, Sander Pronk, Erik Lindahl.
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   1) Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
   2) Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
   3) Neither the name of the copyright holders nor the
   names of its contributors may be used to endorse or promote products
   derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY US ''AS IS'' AND ANY
   EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL WE BE LIABLE FOR ANY
   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   If you want to redistribute modifications, please consider that
   scientific software is very special. Version control is crucial -
   bugs must be traceable. We will be happy to consider code for
   inclusion in the official distribution, but derived work should not
   be called official thread_mpi. Details are found in the README & COPYING
   files.
 */

/** \file
 * \brief A C++11 compatible system_error class for reporting exceptions
 *
 * This header contains class definitions for system_error.
 */

#ifndef TMPI_SYSTEM_ERROR_H_
#define TMPI_SYSTEM_ERROR_H_

#include <stdexcept>

#include "visibility.h"

#ifdef __cplusplus

namespace tMPI
{
/*! \brief Subset of the C++11 system_error class

   Only contains the errno-based constructor. */
class system_error : public std::runtime_error
{
    public:
        //! Type to represent error codes within this class.
        typedef int error_code;

        //system_error(error_code ec, const std::string& what_arg);
        //system_error(error_code ec, const char* what_arg);
        /*! \brief Constuctor that takes an system error number */
        system_error(error_code ec);

        /*! \brief Returns the error code */
        const error_code &code() const
        {
            return ec_;
        }
    private:
        error_code ec_;
};
}

#endif /* __cplusplus */

#endif /* TMPI_SYSTEM_ERROR_H_ */
