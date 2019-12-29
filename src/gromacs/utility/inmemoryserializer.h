/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019, by the GROMACS development team, led by
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
 * Declares gmx::ISerializer implementation for in-memory serialization.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_INMEMORYSERIALIZER_H
#define GMX_UTILITY_INMEMORYSERIALIZER_H

#include <cstddef>

#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/iserializer.h"

namespace gmx
{

//! Specify endian swapping behavoir.
//
// The host-dependent choices avoid the calling file having to
// depend on config.h.
//
enum class EndianSwapBehavior : int
{
    DoNotSwap,                //!< Don't touch anything
    Swap,                     //!< User-enforced swapping
    SwapIfHostIsBigEndian,    //!< Only swap if machine we execute on is big-endian
    SwapIfHostIsLittleEndian, //!< Only swap if machine we execute on is little-endian
    Count                     //!< Number of possible behaviors
};

class InMemorySerializer : public ISerializer
{
public:
    explicit InMemorySerializer(EndianSwapBehavior endianSwapBehavior = EndianSwapBehavior::DoNotSwap);
    ~InMemorySerializer() override;

    std::vector<char> finishAndGetBuffer();

    // From ISerializer
    bool reading() const override { return false; }
    void doBool(bool* value) override;
    void doUChar(unsigned char* value) override;
    void doChar(char* value) override;
    void doUShort(unsigned short* value) override;
    void doInt(int* value) override;
    void doInt32(int32_t* value) override;
    void doInt64(int64_t* value) override;
    void doFloat(float* value) override;
    void doDouble(double* value) override;
    void doReal(real* value) override;
    void doIvec(ivec* value) override;
    void doRvec(rvec* value) override;
    void doString(std::string* value) override;
    void doOpaque(char* data, std::size_t size) override;

private:
    class Impl;

    PrivateImplPointer<Impl> impl_;
};

class InMemoryDeserializer : public ISerializer
{
public:
    InMemoryDeserializer(ArrayRef<const char> buffer,
                         bool                 sourceIsDouble,
                         EndianSwapBehavior   endianSwapBehavior = EndianSwapBehavior::DoNotSwap);
    ~InMemoryDeserializer() override;

    //! Get if the source data was written in double precsion
    bool sourceIsDouble() const;

    // From ISerializer
    bool reading() const override { return true; }
    void doBool(bool* value) override;
    void doUChar(unsigned char* value) override;
    void doChar(char* value) override;
    void doUShort(unsigned short* value) override;
    void doInt(int* value) override;
    void doInt32(int32_t* value) override;
    void doInt64(int64_t* value) override;
    void doFloat(float* value) override;
    void doDouble(double* value) override;
    void doReal(real* value) override;
    void doIvec(ivec* value) override;
    void doRvec(rvec* value) override;
    void doString(std::string* value) override;
    void doOpaque(char* data, std::size_t size) override;

private:
    class Impl;

    PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
