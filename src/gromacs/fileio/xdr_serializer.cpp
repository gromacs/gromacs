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
#include "gmxpre.h"

#include "xdr_serializer.h"

#include "gromacs/fileio/xdrf.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

XdrSerializer::XdrSerializer(const std::filesystem::path& filename, const char* mode) :
    reading_(mode[0] == 'r'), doublePrecision_(GMX_DOUBLE)
{
    std::string binaryMode(mode);
    if (binaryMode.find('b') == std::string::npos)
    {
        binaryMode += 'b'; // Always open XDR in binary mode
    }
    // TODO this always makes a backup, which is not clearly what is
    // wanted in all cases.
    fp_ = gmx_ffopen(filename, binaryMode.c_str());
    if (!fp_)
    {
        GMX_THROW(FileIOError("Failed to open XDR serializer to file " + filename.string()));
    }
    const enum xdr_op xdrMode = reading_ ? XDR_DECODE : XDR_ENCODE;
    // Note there is no documented failure mode for this call
    xdrstdio_create(&xdr_, fp_, xdrMode);
}

XdrSerializer::~XdrSerializer()
{
    xdr_destroy(&xdr_);
    gmx_ffclose(fp_);
}

bool XdrSerializer::reading() const
{
    return reading_;
}

void XdrSerializer::setDoublePrecision(const bool doublePrecision)
{
    doublePrecision_ = doublePrecision;
}

void XdrSerializer::doBool(bool* value)
{
    // XDR can't do bool natively, so convert to integer
    int intValue = static_cast<int>(*value);
    doInt(&intValue);
    *value = static_cast<bool>(intValue);
}

void XdrSerializer::doUChar(unsigned char* value)
{
    if (xdr_u_char(&xdr_, value) == 0)
    {
        GMX_THROW(FileIOError("Failed to serialize unsigned char value to XDR"));
    }
}

void XdrSerializer::doChar(char* value)
{
    if (xdr_char(&xdr_, value) == 0)
    {
        GMX_THROW(FileIOError("Failed to serialize char value to XDR"));
    }
}

void XdrSerializer::doUShort(unsigned short* value)
{
    if (xdr_u_short(&xdr_, value) == 0)
    {
        GMX_THROW(FileIOError("Failed to serialize unsigned short value to XDR"));
    }
}

void XdrSerializer::doInt(int* value)
{
    if (xdr_int(&xdr_, value) == 0)
    {
        GMX_THROW(FileIOError("Failed to serialize int value to XDR"));
    }
}

void XdrSerializer::doInt32(int32_t* value)
{
    if (xdr_int32(&xdr_, value) == 0)
    {
        GMX_THROW(FileIOError("Failed to serialize int32 value to XDR"));
    }
}

void XdrSerializer::doInt64(int64_t* value)
{
    if (xdr_int64(&xdr_, value) == 0)
    {
        GMX_THROW(FileIOError("Failed to serialize int64 value to XDR"));
    }
}

void XdrSerializer::doFloat(float* value)
{
    if (xdr_float(&xdr_, value) == 0)
    {
        GMX_THROW(FileIOError("Failed to serialize float value to XDR"));
    }
}

void XdrSerializer::doDouble(double* value)
{
    if (xdr_double(&xdr_, value) == 0)
    {
        GMX_THROW(FileIOError("Failed to serialize double value to XDR"));
    }
}

void XdrSerializer::doReal(real* value)
{
    if (doublePrecision_)
    {
        double temp = static_cast<double>(*value);
        if (xdr_double(&xdr_, &temp) == 0)
        {
            GMX_THROW(FileIOError("Failed to serialize real(double) value to XDR"));
        }
        *value = static_cast<real>(temp);
    }
    else
    {
        float temp = static_cast<float>(*value);
        if (xdr_float(&xdr_, &temp) == 0)
        {
            GMX_THROW(FileIOError("Failed to serialize real(float) value to XDR"));
        }
        *value = static_cast<real>(temp);
    }
}

void XdrSerializer::doIvec(ivec* value)
{
    doInt(&(*value)[XX]);
    doInt(&(*value)[YY]);
    doInt(&(*value)[ZZ]);
}

void XdrSerializer::doIvec(IVec* value)
{
    doInt(&(*value)[XX]);
    doInt(&(*value)[YY]);
    doInt(&(*value)[ZZ]);
}

void XdrSerializer::doRvec(rvec* value)
{
    // RVec (and rvec) tend to be default-initialized (e.g. with
    // std::vector<RVec>::resize(n) which means aggregate
    // initialization, i.e no initialization. Such unintialized
    // memory can be filled here when reading, but we must take care
    // not to try to convert such values to the other precision, which
    // can trigger a floating-point exception.
    if (doublePrecision_)
    {
        dvec    temp;
        double* d;
        if constexpr (std::is_same_v<rvec, dvec>)
        {
            // No need for conversions
            d = reinterpret_cast<double*>(&(*value)[0]);
        }
        else
        {
            if (!reading())
            {
                // Convert precision when writing
                temp[XX] = (*value)[XX];
                temp[YY] = (*value)[YY];
                temp[ZZ] = (*value)[ZZ];
            }
            d = &(temp[0]);
        }
        if (xdr_vector(&xdr_,
                       reinterpret_cast<char*>(d),
                       DIM,
                       static_cast<unsigned int>(sizeof(double)),
                       reinterpret_cast<xdrproc_t>(xdr_double))
            == 0)
        {
            GMX_THROW(FileIOError("Failed to serialize rvec(double) value to XDR"));
        }
        if constexpr (!std::is_same_v<rvec, dvec>)
        {
            if (reading())
            {
                // Convert precision when reading
                (*value)[XX] = static_cast<real>(temp[XX]);
                (*value)[YY] = static_cast<real>(temp[YY]);
                (*value)[ZZ] = static_cast<real>(temp[ZZ]);
            }
        }
    }
    else
    {
        float  temp[DIM];
        float* f;
        if constexpr (std::is_same_v<rvec, dvec>)
        {
            if (!reading())
            {
                temp[XX] = (*value)[XX];
                temp[YY] = (*value)[YY];
                temp[ZZ] = (*value)[ZZ];
            }
            f = &(temp[0]);
        }
        else
        {
            // No need for conversions
            f = reinterpret_cast<float*>(&((*value)[0]));
        }
        if (xdr_vector(&xdr_,
                       reinterpret_cast<char*>(f),
                       DIM,
                       static_cast<unsigned int>(sizeof(float)),
                       reinterpret_cast<xdrproc_t>(xdr_float))
            == 0)
        {
            GMX_THROW(FileIOError("Failed to serialize rvec(float) value to XDR"));
        }
        if constexpr (std::is_same_v<rvec, dvec>)
        {
            if (reading())
            {
                (*value)[XX] = static_cast<real>(temp[XX]);
                (*value)[YY] = static_cast<real>(temp[YY]);
                (*value)[ZZ] = static_cast<real>(temp[ZZ]);
            }
        }
    }
}

void XdrSerializer::doRvec(RVec* value)
{
    doRvec(reinterpret_cast<rvec*>(value));
}

void XdrSerializer::doCharArray(char* values, int elements)
{
    GMX_RELEASE_ASSERT(elements < static_cast<int>(std::numeric_limits<int>::max()),
                       "The XDR interface cannot handle array lengths > 2^31");
    if (xdr_vector(&xdr_,
                   values,
                   static_cast<int>(elements),
                   static_cast<unsigned int>(sizeof(char)),
                   reinterpret_cast<xdrproc_t>(xdr_char))
        == 0)
    {
        GMX_THROW(FileIOError("Failed to serialize array of chars to XDR"));
    }
}

void XdrSerializer::doUCharArray(unsigned char* values, int elements)
{
    GMX_RELEASE_ASSERT(elements < static_cast<int>(std::numeric_limits<int>::max()),
                       "The XDR interface cannot handle array lengths > 2^31");
    if (xdr_vector(&xdr_,
                   // xdr_vector expects a char* regardless
                   reinterpret_cast<char*>(values),
                   static_cast<int>(elements),
                   static_cast<unsigned int>(sizeof(unsigned char)),
                   reinterpret_cast<xdrproc_t>(xdr_u_char))
        == 0)
    {
        GMX_THROW(FileIOError("Failed to serialize array of unsigned chars to XDR"));
    }
}

void XdrSerializer::doRvecArray(ArrayRef<RVec> values)
{
    // This would probably be more efficient as a single call to
    // xdr_vector (perhaps preceded by a length value?), but we have
    // many files already written that followed this approach and we
    // have to be able to read them.
    for (RVec& value : values)
    {
        doRvec(&value);
    }
}

void XdrSerializer::doString(std::string* value)
{
    // Historically GROMACS serialized both the terminating null
    // character and a string length (which included that null
    // character). In particular, an empty string was serialized with
    // a length of 1 and a single null character. So that XDR files
    // written by GROMACS can still be read, this is perpetuated.
    const int stringLength            = static_cast<int>(value->size());
    int       stringLengthToSerialize = stringLength + 1;
    if (xdr_int(&xdr_, &stringLengthToSerialize) == 0)
    {
        GMX_THROW(FileIOError("Failed to serialize string length to XDR"));
    }
    if (reading())
    {
        value->resize(stringLengthToSerialize - 1);
    }
    // This is guaranteed by the C++ standard to have a terminating
    // null character, i.e we can do both reading and writing and
    // access that null character, so long as we replace it when
    // writing.
    char* charPtr = value->data();
    if (xdr_string(&xdr_, &charPtr, stringLengthToSerialize) == 0)
    {
        GMX_THROW(FileIOError("Failed to serialize string contents to XDR"));
    }
}

void XdrSerializer::doOpaque(char* data, std::size_t size)
{
    // GROMACS 2020 and later wants to embed the entire TPR body as a
    // single opaque object. The default XDR interface only uses
    // integers for the size field, which limits the size of such an
    // object. So we need to make sure that we break opaque data into
    // suitably small chunks, so that TPR files larger than 2GB will
    // work.
    //
    // In general not all chunks will have the same size (and we do not
    // want to pad them), so we calculate the chunk size as:
    // - The max value of a signed integer + 1
    // - Subtract 4 (the XDR object size) to get a size within the
    //   range of the signed int.
    const std::size_t maxChunk = static_cast<std::size_t>(std::numeric_limits<int>::max()) + 1 - 4;

    size_t offset = 0;
    bool_t result;
    for (result = 1; result > 0 && size > 0;)
    {
        const std::size_t thisChunk = std::min(maxChunk, size);
        result                      = xdr_opaque(&xdr_, data + offset, thisChunk);
        offset += thisChunk;
        size -= thisChunk;
    }
    if (result == 0)
    {
        GMX_THROW(FileIOError("Failed to serialize string contents to XDR"));
    }
}

XDR* XdrSerializer::xdr()
{
    return &xdr_;
}

} // namespace gmx
