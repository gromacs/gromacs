/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015 by the GROMACS development team, led by
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
    \ingroup module_data_structures
    \inpublicapi

   \brief
   This header chooses between C++11 and C++98 variants of the symbol table.
   t_symbtab. Class t_symtab stores std::strings non-redundantly (thread-safe)
   and trimmed to the chars actually used. Global replacement is possible,
   e.g., for renaming all residues of a certain type in GROMACS.

   \author
    R. Thomas Ullmann <thomas.ullmann@mpibpc.mpg.de>

   \copyright
    GROMACS license

   \date Sep 2015
 */

#ifndef gmx_symtab_h
#define gmx_symtab_h

// NVCC supports C++11 beginning with version 7.0, but most
// likely we will anyway not need the shared_strings on the GPU
#ifdef  __NVCC__
 #ifdef __CUDA_API_VERSION
  #if __CUDA_API_VERSION < 7000
   #define old_nvcc_symtab_cpp_h
  #endif
 #else
  #define old_nvcc_symtab_cpp_h
 #endif
#endif

#include <string>
#include <cstring>
#include <map>
#include <iostream>
#include <ios>
#include <stdexcept>

// check for C++11 standard availability once
#if __cplusplus < 201103L || defined old_nvcc_symtab_cpp_h
 #include "gromacs/data_structures/symtab_cpp98.h"
#else
 #include "gromacs/data_structures/symtab_cpp11.h"
#endif

#endif // end header
