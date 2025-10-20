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

/*! \brief Declarations of H5md fixed data set class.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FILEIO_H5MD_FIXEDDATASET_H
#define GMX_FILEIO_H5MD_FIXEDDATASET_H

#include <hdf5.h>

#include "gromacs/fileio/h5md/h5md_datasetbase.h"

namespace gmx
{

template<typename ValueType>
class ArrayRef;
template<typename ValueType>
class BasicVector;

/*! \brief Class which provides an interface for reading and writing data to data sets of fixed size.
 *
 * \tparam ValueType Type of values stored in data set (not limited to primitives).
 */
template<typename ValueType>
class H5mdFixedDataSet : private H5mdDataSetBase<ValueType>
{
public:
    using Base = H5mdDataSetBase<ValueType>;

    //! \brief Constructor to manage a given \p dataSet.
    H5mdFixedDataSet(H5mdDataSetBase<ValueType>&& dataSet);

    //! \brief Constructor to open an existing data set (called \p name, in \p container).
    H5mdFixedDataSet(hid_t container, const char* name);

    /*! \brief Return the data set dimensions of the templated ValueType.
     *
     * \note For data sets of e.g. BasicVector<T> = T[3] this refers to the actual dimensions
     * in terms of BasicVector<T>, not the base primitive T. So while data set with shape
     * BasicVector<float>[30][50] has primitive dimensions float[30][50][3], the dimensions
     * returned here are [30][50]. The returned dimensions are thus consistent with the size
     * of input buffers ArrayRef<T> used for the I/O methods \c readData and \c writeData.
     */
    const DataSetDims& dims() const;

    /*! \brief Return the total number of values, of ValueType, stored in the data set.
     */
    hsize_t numValues() const;

    /*! \brief Read all data from the data set into the given buffer \p data.
     *
     * The size of the input \p data buffer must match that of the data set dims(). If the data set
     * is multidimensional one must cast the data into a 1-dimensional form of row-major order.
     *
     * \param[out] data Buffer to read data into from the data set.
     *
     * \throws gmx::FileIOError if the buffer size does not match the data set.
     */
    void readData(ArrayRef<ValueType> data) const;

    /*! \brief Write to all values in the data set from the given buffer \p data.
     *
     * The size of the input \p data buffer must match that of the data set dims(). If the data set
     * is multidimensional one must cast the data into a 1-dimensional form of row-major order.
     *
     * \param[in] data Buffer with values to write into the data set.
     *
     * \throws gmx::FileIOError if the buffer size does not match the data set.
     */
    void writeData(ArrayRef<const ValueType> data) const;

private:
    //! \brief Data set dimensions for the templated ValueType.
    DataSetDims dims_;

    //! \brief The total number of values in the data set.
    hsize_t numValues_;
};

extern template class H5mdFixedDataSet<int32_t>;

extern template class H5mdFixedDataSet<int64_t>;

extern template class H5mdFixedDataSet<float>;

extern template class H5mdFixedDataSet<double>;

extern template class H5mdFixedDataSet<BasicVector<float>>;

extern template class H5mdFixedDataSet<BasicVector<double>>;

} // namespace gmx

#endif // GMX_FILEIO_H5MD_FIXEDDATASET_H
