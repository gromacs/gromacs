/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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

/*! \brief Declarations of H5md compression routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FILEIO_H5MD_COMPRESSION_H
#define GMX_FILEIO_H5MD_COMPRESSION_H

#include <hdf5.h>

#include <optional>

namespace gmx
{
//! \brief Registered id for the SZ3 filter in HDF5
constexpr int c_sz3FilterId = 32024;

//! \brief Data set compression options
enum class H5mdCompression : int
{
    Uncompressed,      //!< No compression
    LosslessNoShuffle, //!< Lossless, only gzip (deflate)
    LosslessShuffle,   //!< Lossless, byte shuffle followed by gzip (deflate)
    LossySz3BioMD,     //!< Lossy, SZ3 BioMD algorithm with an absolute error
    Count,             //!< Number of compression schemes
};

//! \brief Return string description for input \p compression scheme.
const char* enumValueToString(H5mdCompression compression);

/*! \brief Return whether input \p dataSet uses SZ3 compression.
 *
 * Can be called regardless of whether GROMACS is compiled with SZ3 support
 * or not. Should be used by non-SZ3 builds to check whether we are
 * opening a data set which we cannot interact with, i.e. if SZ3
 * compression is used for it.
 *
 * \throws H5mdError if the creation property list could not be retrieved
 * from \p dataSet.
 */
bool isCompressedWithSz3(hid_t dataSet);

/*! \brief Return the absolute error for SZ3 compression in \p dataSet.
 *
 * If a lossy SZ3 compression scheme is applied to the given \p dataSet
 * the absolute error is retrieved and returned. Otherwise, including
 * if no filter is applied, std::nullopt is returned.
 *
 * \throws H5mdError if there was an error getting the SZ3 filter
 * configuration.
 *
 * \throws NotImplementedError if GROMACS is not built with SZ3 support.
 */
std::optional<double> getSz3AbsoluteError(hid_t dataSet);

/*! \brief Register the active SZ3 filter and return if successful.
 *
 * This replaces, for the current process, an SZ3 filter that HDF5
 * may have loaded from HDF5_PLUGIN_PATH with the filter that GROMACS
 * is built against.
 *
 * \throws NotImplementedError if GROMACS is not built with SZ3 support.
 */
bool registerSz3Filter();


/*! \brief Set HDF5 filters for input \p compression scheme to \p creationPropertyList.
 *
 * \note Filters must be applied before any checksum filters are
 * (e.g. Fletcher-32).
 *
 * \param[in] creationPropertyList Property list for data set creation.
 * \param[in] compression          Compression scheme.
 * \param[in] absoluteError        Maximum absolute error for lossy compression schemes
 *                                 (unused for non-lossy schemes).
 *
 * \throws H5mdError if \p compression is a lossy compression scheme
 * and \p absoluteError is not positive, or if the compression filter
 * could not be applied to \p creationPropertyList.
 *
 * \throws NotImplementedError if GROMACS is not built with SZ3 support.
 */
void setCompressionFilters(hid_t creationPropertyList, H5mdCompression compression, double absoluteError = 0.0);

} // namespace gmx

#endif // GMX_FILEIO_H5MD_COMPRESSION_H
