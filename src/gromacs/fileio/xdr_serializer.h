/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
#ifndef GMX_FILEIO_XDR_SERIALIZER_H
#define GMX_FILEIO_XDR_SERIALIZER_H

#include "config.h"

#include <filesystem>
#include <string>

#include "gromacs/serialization/iserializer.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/vectypes.h"

#if GMX_INTERNAL_XDR
#    include "rpc_xdr/xdr.h"
#else
#    include <rpc/rpc.h>
#    include <rpc/xdr.h>
#endif

struct t_fileio;

namespace gmx
{

/*! \internal
 * \brief Serializer to read/write XDR data.
 *
 * Most methods on this class throw FileIOError when the
 * serialization failed.
 */
class XdrSerializer : public ISerializer
{
public:
    XdrSerializer(const std::filesystem::path& filename, const char* mode);
    ~XdrSerializer() override;

    //! If file is open in reading mode.
    bool reading() const override;
    /*! \brief Set when file should read in double-precision mode
     *
     * Often this is described in the GROMACS header of an XDR file,
     * so when reading from a build of GROMACS that does not match the
     * one that wrote the XDR file, the serializer needs to be able to
     * be told how to expect "real" values to look on disk. */
    void setDoublePrecision(bool doublePrecision);
    //! Handle bool I/O.
    void doBool(bool* value) override;
    //! Handle unsigned char I/O.
    void doUChar(unsigned char* value) override;
    //! Handle char I/O.
    void doChar(char* value) override;
    //! Handle unsigned short I/O.
    void doUShort(unsigned short* value) override;
    //! Handle default integer I/O.
    void doInt(int* value) override;
    //! Handle int32 I/O.
    void doInt32(int32_t* value) override;
    //! Handle int64 I/O.
    void doInt64(int64_t* value) override;
    //! Handle single precision float I/O.
    void doFloat(float* value) override;
    //! Handle double precision float I/O.
    void doDouble(double* value) override;
    //! Handle GROMACS floating point number I/O.
    void doReal(real* value) override;
    //! Handle I/O of integer vector of size DIM.
    void doIvec(ivec* value) override;
    //! Handle I/O of integer vector of size DIM.
    void doIvec(IVec* value) override;
    //! Handle I/O of GROMACS real vector of size DIM.
    void doRvec(rvec* value) override;
    //! Handle I/O of GROMACS real vector of size DIM.
    void doRvec(RVec* value) override;
    //! Handle I/O if string.
    void doString(std::string* value) override;
    //! Handle opaque data.
    void doOpaque(char* data, std::size_t size) override;
    //! Special case for handling I/O of a vector of characters.
    void doCharArray(char* values, int elements) override;
    //! Special case for handling I/O of a vector of unsigned characters.
    void doUCharArray(unsigned char* values, int elements) override;
    //! Special case for handling I/O of an ArrayRef of RVec.
    void doRvecArray(ArrayRef<RVec> values) override;

    //! Temporary getter while XDR serialization still uses legacy interfaces
    XDR* xdr();

private:
    //! File handle opened during construction
    FILE* fp_ = nullptr;
    /*! \brief XDR I/O handle using external file pointer.
     *
     * Used when interacting with legacy routines not yet ported
     * to use this serializer. */
    XDR xdr_;
    //! Whether the file is reading or writing
    bool reading_ = false;
    //! Whether the build that wrote (or is writing) this file was double precision
    bool doublePrecision_ = false;
};

} // namespace gmx

#endif
