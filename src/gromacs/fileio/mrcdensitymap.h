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
/*! \libinternal \file
 * \brief
 * Declars mrc/ccp4-file format handling.
 *
 * \author Christian Blau <blau@kth.se>
 *
 * \inlibraryapi
 * \ingroup module_fileio
 */
#include <string>
#include <vector>

#include "gromacs/math/multidimarray.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

struct MrcDensityMapHeader;
class ISerializer;
class TranslateAndScale;

/*! \libinternal \brief Read an mrc/ccp4 file that contains float values.
 */
class MrcDensityMapOfFloatReader
{
public:
    /*! \brief Construct from directly de-serializing data into the object.
     * \throws InternalError if serializer is not reading
     * \throws InternalError if header is inconsistent
     * \throws if serializer throws error upon failed reading
     * \param[in] serializer Serializer to read the object data from
     */
    explicit MrcDensityMapOfFloatReader(ISerializer* serializer);

    ~MrcDensityMapOfFloatReader();

    //! Return a view on the data of the density grid
    ArrayRef<const float> constView() const;
    //! Return the header
    const MrcDensityMapHeader& header() const;

private:
    class Impl;
    PrivateImplPointer<Impl> impl_;
};

/*! \libinternal \brief Read an mrc density map from a given file.
 *
 * Higher level class than MrcDensityMapOfFloatReader that takes a file name
 * upon construction and returns coordinate transformation into the density
 * lattice as well as the density data.
 *
 * Attempts reading with swapped endianess if header is not sane.
 *
 * Performs basic sanity checks on header information and data size.
 *
 * \note File reading is completed during construction. When the constructor
 *       completes succesfully, transformation to density lattice and density
 *       data are valid, irrespective of the state of the read file.
 */
class MrcDensityMapOfFloatFromFileReader
{
public:
    MrcDensityMapOfFloatFromFileReader();

    /*! \brief Read from filename.
     * \throws FileIOError if file does not exist
     * \throws FileIOError if read in buffer size does not match file size
     * \throws FileIOError if header information does not match density
     *                     data size
     */
    explicit MrcDensityMapOfFloatFromFileReader(const std::string& filename);

    ~MrcDensityMapOfFloatFromFileReader();

    //! Return the coordinate transformation into the density
    TranslateAndScale transformationToDensityLattice() const;

    //! Return a copy of the density data
    MultiDimArray<std::vector<float>, dynamicExtents3D> densityDataCopy() const;

private:
    class Impl;
    PrivateImplPointer<Impl> impl_;
};

/*! \libinternal \brief Write an mrc/ccp4 file that contains float values.
 */
class MrcDensityMapOfFloatWriter
{
public:
    /*! \brief Construct by setting the data and the header.
     *
     * \throws if the header data description does not match the provided data
     *
     * \param[in] header mrc density map header
     * \param[in] data the density map data
     */
    MrcDensityMapOfFloatWriter(const MrcDensityMapHeader& header, ArrayRef<const float> data);

    ~MrcDensityMapOfFloatWriter();

    //! Serialize the mrc density data.
    void write(ISerializer* serializer) const;

private:
    class Impl;
    PrivateImplPointer<Impl> impl_;
};

} // namespace gmx
