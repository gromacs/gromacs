/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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

#include <algorithm>
#include <vector>

#include "gromacs/fileio/mrcdensitymapheader.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/iserializer.h"

#include "mrcserializer.h"

namespace gmx
{

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
        explicit Impl(ISerializer *serializer);
        ~Impl() {}
        //! The header of the read mrc file
        MrcDensityMapHeader header_;
        //! The data of the mrc file
        std::vector<float>  data_;
};

MrcDensityMapOfFloatReader::Impl::Impl(ISerializer *serializer)
{
    if (!serializer->reading())
    {
        GMX_THROW(InternalError("Cannot use writing serializer to read."));
    }

    header_ = deserializeMrcDensityMapHeader(serializer);
    const auto dataSize = numberOfExpectedDataItems(header_);
    data_.resize(dataSize);
    for (auto &value : data_)
    {
        serializer->doFloat(&value);
    }
}

/********************************************************************
 * MrcDensityMapOfFloatReader
 */

MrcDensityMapOfFloatReader::MrcDensityMapOfFloatReader(ISerializer *serializer) : impl_(new Impl(serializer))
{
}

ArrayRef<const float> MrcDensityMapOfFloatReader::data() const
{
    return impl_->data_;
}

const MrcDensityMapHeader &MrcDensityMapOfFloatReader::header() const
{
    return impl_->header_;
}

MrcDensityMapOfFloatReader::~MrcDensityMapOfFloatReader()
{
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
        Impl(const MrcDensityMapHeader &header, ArrayRef<const float> data);
        ~Impl() {}
        //! Write the header and data from the writer to a given serialier
        void write(ISerializer *serializer) const;
        //! The mrc density map header data
        const MrcDensityMapHeader   header_;
        //! The density data
        const ArrayRef<const float> data_;
};

MrcDensityMapOfFloatWriter::Impl::Impl(const MrcDensityMapHeader &header, ArrayRef<const float> data) : header_(header), data_(data)
{
}

void MrcDensityMapOfFloatWriter::Impl::write(ISerializer *serializer) const
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

MrcDensityMapOfFloatWriter::MrcDensityMapOfFloatWriter(const MrcDensityMapHeader &header, ArrayRef<const float> data) : impl_(new Impl(header, data))
{
}

void MrcDensityMapOfFloatWriter::write(ISerializer *serializer) const
{
    impl_->write(serializer);
}

MrcDensityMapOfFloatWriter::~MrcDensityMapOfFloatWriter()
{
}

} // namespace gmx
