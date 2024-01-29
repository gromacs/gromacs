/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
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
    Unknown,    //!< Not one of the cases below
    X86,        //!< X86
    Arm,        //!< ARM
    PowerPC,    //!< IBM PowerPC
    RiscV32,    //!< 32-bit RISC-V
    RiscV64,    //!< 64-bit RISC-V
    Loongarch64 //!< 64-bit Loongarch
};

//! Whether the compilation is targeting 32-bit x86.
#if (defined __i386__ || defined __i386 || defined _X86_ || defined _M_IX86)
#    define GMX_IS_X86_32 1
#else
#    define GMX_IS_X86_32 0
#endif

//! Whether the compilation is targeting 64-bit x86.
#if (defined __x86_64__ || defined __x86_64 || defined __amd64__ || defined __amd64 \
     || defined _M_X64 || defined _M_AMD64)
#    define GMX_IS_X86_64 1
#else
#    define GMX_IS_X86_64 0
#endif

//! Constant that tells what the architecture is
static constexpr Architecture c_architecture =
#if GMX_IS_X86_32 || GMX_IS_X86_64
        Architecture::X86;
#elif defined __arm__ || defined __arm || defined _M_ARM || defined __aarch64__
        Architecture::Arm;
#elif defined __powerpc__ || defined __ppc__ || defined __PPC__
        Architecture::PowerPC;
#elif defined __riscv && defined __riscv_xlen && (__riscv_xlen == 32)
        Architecture::RiscV32;
#elif defined __riscv && defined __riscv_xlen && (__riscv_xlen == 64)
        Architecture::RiscV64;
#elif defined __loongarch__ && defined __loongarch64
        Architecture::Loongarch64;
#else
        Architecture::Unknown;
#endif

} // namespace gmx

#endif // GMX_ARCHITECTURE_H
