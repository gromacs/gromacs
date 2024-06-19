/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 *
 * \brief
 * Implements mrc/ccp4-file format handling.
 *
 * \author Christian Blau <blau@kth.se>
 *
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "mrcdensitymap.h"

#include <cstdio>

#include <algorithm>
#include <iterator>
#include <vector>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/gmxfio_xdr.h"
#include "gromacs/fileio/mrcdensitymapheader.h"
#include "gromacs/math/coordinatetransformation.h"
#include "gromacs/mdspan/extents.h"
#include "gromacs/mdspan/layouts.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/inmemoryserializer.h"
#include "gromacs/utility/iserializer.h"
#include "gromacs/utility/stringutil.h"

#include "mrcserializer.h"

namespace gmx
{

namespace
{

/*! \brief Read file into memory as vector of chars
 *
 * \param[in] filename of the file to be read
 * \returns the file contents as a vector
 * \throws FileIOError if file not found
 * \throws FileIOError if reading was not successful
 */
std::vector<char> readCharBufferFromFile(const std::filesystem::path& filename)
{
    if (!gmx_fexist(filename))
    {
        GMX_THROW(FileIOError(gmx::formatString("Error while reading '%s' - file not found.",
                                                filename.string().c_str())));
    }
    t_fileio* mrcFile = gmx_fio_open(filename, "r");

    // Determine file size
    gmx_fseek(gmx_fio_getfp(mrcFile), 0, SEEK_END);
    gmx_off_t fileSize = gmx_fio_ftell(mrcFile);
    gmx_fseek(gmx_fio_getfp(mrcFile), 0, SEEK_SET);
    // Read whole file into buffer the size of the file
    std::vector<char> fileContentBuffer(fileSize);
    size_t            readSize = fread(
            fileContentBuffer.data(), sizeof(char), fileContentBuffer.size(), gmx_fio_getfp(mrcFile));
    gmx_fio_close(mrcFile);

    if (fileContentBuffer.size() != readSize)
    {
        GMX_THROW(FileIOError(gmx::formatString(
                "Error while reading '%s' - file size and read buffer size do not match.",
                filename.string().c_str())));
    }

    return fileContentBuffer;
}

} // namespace

/********************************************************************
 * MrcDensityMapOfFloatReader::Impl
 */

/*! \internal \brief
 * Private implementation class for MrcDensityMapOfFloatReader.
 */
class MrcDensityMapOfFloatReader::Impl
{
public:
    //! Build the map reader from a serializer.
    explicit Impl(ISerializer* serializer);
    ~Impl() {}
    //! The header of the read mrc file
    MrcDensityMapHeader header_;
    //! The data of the mrc file
    std::vector<float> data_;
};

MrcDensityMapOfFloatReader::Impl::Impl(ISerializer* serializer)
{
    if (!serializer->reading())
    {
        GMX_THROW(InternalError("Cannot use writing serializer to read."));
    }

    header_             = deserializeMrcDensityMapHeader(serializer);
    const auto dataSize = numberOfExpectedDataItems(header_);
    data_.resize(dataSize);
    for (auto& value : data_)
    {
        serializer->doFloat(&value);
    }
}

/********************************************************************
 * MrcDensityMapOfFloatReader
 */

MrcDensityMapOfFloatReader::MrcDensityMapOfFloatReader(ISerializer* serializer) :
    impl_(new Impl(serializer))
{
}

ArrayRef<const float> MrcDensityMapOfFloatReader::constView() const
{
    return impl_->data_;
}

const MrcDensityMapHeader& MrcDensityMapOfFloatReader::header() const
{
    return impl_->header_;
}

MrcDensityMapOfFloatReader::~MrcDensityMapOfFloatReader() {}

/********************************************************************
 * MrcDensityMapOfFloatFromFileReader::Impl
 */


class MrcDensityMapOfFloatFromFileReader::Impl
{
public:
    explicit Impl(const std::filesystem::path& fileName);
    ~Impl() = default;
    const MrcDensityMapOfFloatReader& reader() const;

private:
    const std::vector<char>                     buffer_;
    std::unique_ptr<InMemoryDeserializer>       serializer_;
    std::unique_ptr<MrcDensityMapOfFloatReader> reader_;
};

MrcDensityMapOfFloatFromFileReader::Impl::Impl(const std::filesystem::path& filename) :
    buffer_(readCharBufferFromFile(filename)),
    serializer_(std::make_unique<InMemoryDeserializer>(buffer_, false)),
    reader_(std::make_unique<MrcDensityMapOfFloatReader>(serializer_.get()))
{
    if (!mrcHeaderIsSane(reader_->header()))
    {
        serializer_ = std::make_unique<InMemoryDeserializer>(buffer_, false, EndianSwapBehavior::Swap);
        reader_     = std::make_unique<MrcDensityMapOfFloatReader>(serializer_.get());
        if (!mrcHeaderIsSane(reader_->header()))
        {
            GMX_THROW(FileIOError(gmx::formatString(
                    "Header of '%s' fails sanity check for little- as well as big-endian reading.",
                    filename.string().c_str())));
        }
    }

    layout_right::mapping<dynamicExtents3D> map(getDynamicExtents3D(reader_->header()));
    if (map.required_span_size() != reader_->constView().ssize())
    {
        GMX_THROW(
                FileIOError(gmx::formatString("File header density extent information of '%s'"
                                              " does not match density data size",
                                              filename.string().c_str())));
    }
}

const MrcDensityMapOfFloatReader& MrcDensityMapOfFloatFromFileReader::Impl::reader() const
{
    return *reader_;
}

/********************************************************************
 * MrcDensityMapOfFloatFromFileReader
 */

MrcDensityMapOfFloatFromFileReader::MrcDensityMapOfFloatFromFileReader(const std::filesystem::path& filename) :
    impl_(new Impl(filename))
{
}

MrcDensityMapOfFloatFromFileReader::~MrcDensityMapOfFloatFromFileReader() = default;

TranslateAndScale MrcDensityMapOfFloatFromFileReader::transformationToDensityLattice() const
{
    return getCoordinateTransformationToLattice(impl_->reader().header());
}

MultiDimArray<std::vector<float>, dynamicExtents3D> MrcDensityMapOfFloatFromFileReader::densityDataCopy() const
{
    MultiDimArray<std::vector<float>, dynamicExtents3D> result(
            getDynamicExtents3D(impl_->reader().header()));
    std::copy(std::begin(impl_->reader().constView()),
              std::end(impl_->reader().constView()),
              begin(result.asView()));
    return result;
}

/********************************************************************
 * MrcDensityMapOfFloatWriter::Impl
 */

/*! \internal \brief
 * Private implementation class for MrcDensityMapOfFloatWriter.
 */
class MrcDensityMapOfFloatWriter::Impl
{
public:
    //! Construct mrc file writer by providing header and data to be written.
    Impl(const MrcDensityMapHeader& header, ArrayRef<const float> data);
    ~Impl() {}
    //! Write the header and data from the writer to a given serialier
    void write(ISerializer* serializer) const;
    //! The mrc density map header data
    const MrcDensityMapHeader header_;
    //! The density data
    const ArrayRef<const float> data_;
};

MrcDensityMapOfFloatWriter::Impl::Impl(const MrcDensityMapHeader& header, ArrayRef<const float> data) :
    header_(header), data_(data)
{
}

void MrcDensityMapOfFloatWriter::Impl::write(ISerializer* serializer) const
{
    if (serializer->reading())
    {
        GMX_THROW(InternalError("Cannot use reading serializer to write."));
    }

    serializeMrcDensityMapHeader(serializer, header_);

    if (numberOfExpectedDataItems(header_) != data_.size())
    {
        GMX_THROW(InternalError("Mrc data size does not match header information."));
    }

    for (float value : data_)
    {
        serializer->doFloat(&value);
    }
}

/********************************************************************
 * MrcDensityMapOfFloatWriter
 */

MrcDensityMapOfFloatWriter::MrcDensityMapOfFloatWriter(const MrcDensityMapHeader& header,
                                                       ArrayRef<const float>      data) :
    impl_(new Impl(header, data))
{
}

void MrcDensityMapOfFloatWriter::write(ISerializer* serializer) const
{
    impl_->write(serializer);
}

MrcDensityMapOfFloatWriter::~MrcDensityMapOfFloatWriter() {}

} // namespace gmx
