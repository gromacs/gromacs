/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 *  \brief Declares common utilities to use with rocfft.
 *
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *  \ingroup module_fft
 */

#ifndef GMX_FFT_ROCFFT_COMMON_UTILS_H
#define GMX_FFT_ROCFFT_COMMON_UTILS_H

#include <rocfft/rocfft.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{
//! Model the kinds of 3D FFT implemented
enum class FftDirection : int
{
    RealToComplex,
    ComplexToReal,
    Count,
};

//! Helper for consistent error handling
void handleRocFftError(rocfft_status result, const std::string& msg);

//! Helper for consistent error handling
void handleRocFftError(rocfft_status result, const std::string& direction, const std::string& msg);

//! Provides RAII-style initialization of rocFFT library
class RocfftInitializer
{
public:
    RocfftInitializer();
    ~RocfftInitializer();
};

//! All the persistent data for planning an executing a 3D FFT
struct RocfftPlan
{
    //! High level information about the plan
    rocfft_plan plan = nullptr;
    //! Execution details (working buffer, HIP stream to use, etc)
    rocfft_execution_info info = nullptr;
    //! Persistent work buffer (left unallocated if not needed)
    void* workBuffer = nullptr;
    //! Destructor
    ~RocfftPlan();
    //! Default constructor
    RocfftPlan() = default;
    //! Full constructor
    RocfftPlan(rocfft_plan inPlan, rocfft_execution_info inInfo, void* inWorkBuffer) :
        plan(inPlan), info(inInfo), workBuffer(inWorkBuffer)
    {
    }
    //! Move construct is allowed
    RocfftPlan(RocfftPlan&& other) = default;
    //! Move assign is allowed
    RocfftPlan& operator=(RocfftPlan&& other) = default;
    //! Copy construct is not allowed
    RocfftPlan(const RocfftPlan& other) = delete;
    //! Copy assign is not allowed
    RocfftPlan& operator=(const RocfftPlan& other) = delete;
};

//! Helper struct to reduce repetitive code setting up a 3D FFT plan
struct PlanSetupData
{
    //! Format of the input array (real or hermitian)
    rocfft_array_type arrayType;
    //! Strides through the input array for the three dimensions
    std::array<size_t, DIM> strides;
    //! Total size of the input array (including padding)
    size_t totalSize;
};

//! Compute the stride through the real 1D array
inline std::array<size_t, DIM> makeRealStrides(ivec realGridSizePadded)
{
    return { 1, size_t(realGridSizePadded[ZZ]), size_t(realGridSizePadded[ZZ] * realGridSizePadded[YY]) };
};

//! Compute the stride through the complex 1D array
inline std::array<size_t, DIM> makeComplexStrides(ivec complexGridSizePadded)
{
    return { 1,
             size_t(complexGridSizePadded[ZZ]),
             size_t(complexGridSizePadded[ZZ] * complexGridSizePadded[YY]) };
}

//! Compute total grid size
inline size_t computeTotalSize(ivec gridSize)
{
    return size_t(gridSize[XX] * gridSize[YY] * gridSize[ZZ]);
}

} // namespace gmx

#endif
