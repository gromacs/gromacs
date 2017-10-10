/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Declares and defines architecture booleans to minimize preprocessed code
 *
 * \ingroup module_hardware
 */
#ifndef GMX_ARCHITECTURE_H
#define GMX_ARCHITECTURE_H

namespace gmx
{

//! Enum for GROMACS CPU hardware detection support
enum class Architecture
{
    Unknown, //! Not one of the cases below
    X86,     //! X86
    Arm,     //! ARM
    PowerPC  //! IBM PowerPC
};

//! Whether the compilation is targeting 32-bit x86.
# if (defined __i386__ || defined __i386 || defined _X86_ || defined _M_IX86)
#define GMX_IS_X86_32 1
#else
#define GMX_IS_X86_32 0
#endif

//! Whether the compilation is targeting 64-bit x86.
#if (defined __x86_64__ || defined __x86_64 || defined __amd64__ || defined __amd64 || defined _M_X64 || defined _M_AMD64)
#define GMX_IS_X86_64 1
#else
#define GMX_IS_X86_64 0
#endif

//! Constant that tells what the architecture is
static constexpr Architecture c_architecture =
#if GMX_IS_X86_32 || GMX_IS_X86_64
    Architecture::X86;
#elif defined __arm__ || defined __arm || defined _M_ARM || defined __aarch64_
    Architecture::Arm;
#elif defined __powerpc__ || defined __ppc__ || defined __PPC__
    Architecture::PowerPC;
#else
    Architecture::Unknown;
#endif

}      // namespace gmx

#endif // GMX_ARCHITECTURE_H
