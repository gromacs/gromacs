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

#ifndef TMPI_VISIBILITY_H_
#define TMPI_VISIBILITY_H_

/** \file
 *
 * \brief Visibility macros
 *
 * These macros enable dynamic library visibility support. Either set the
 * 'TMPI_USE_VISIBILITY', or set the TMPI_EXPORT to the right export
 * statement, when incorporating thread_mpi into a larger project.
 *
 * All exported functions and classes will be tagged by the visibility macro.
 *
 * \sa http://gcc.gnu.org/wiki/Visibility
 *
 */

#ifdef TMPI_USE_VISIBILITY  /* off by default */

/* only set if the macro hasn't been set elsewhere */
#ifndef TMPI_EXPORT

/* gcc-like */
#if defined(__GNUC__)

#define TMPI_EXPORT __attribute__((__visibility__("default")))

#elif defined _WIN32 || defined __CYGWIN__ || defined WINDOWS

#ifdef TMPI_EXPORTS
#define TMPI_EXPORT __declspec(dllexport)
#else
#define TMPI_EXPORT __declspec(dllimport)
#endif

#else /* no viable visibility */

#define TMPI_EXPORT

#endif /* compiler check */

#endif /* TMPI_EXPORT */

#else  /* TMPI_USE_VISIBILITY */

#define TMPI_EXPORT

#endif /* TMPI_USE_VISIBILITY */

#endif /* TMPI_VISIBILITY_H_ */
