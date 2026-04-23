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

/*! \brief Declarations of H5md data set base class.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FILEIO_H5MD_DATASETBASE_H
#define GMX_FILEIO_H5MD_DATASETBASE_H

#include <hdf5.h>

#include <string>
#include <vector>

#include "gromacs/fileio/h5md/h5md_guard.h"

namespace gmx
{

template<typename ValueType>
class BasicVector;
template<typename ValueType>
class H5mdDataSetBuilder;

/*! \brief Dimensions of a HDF5 data set.
 */
using DataSetDims = std::vector<hsize_t>;

/*! \brief Class which manages a HDF5 data set handle and associated data.
 *
 * Can be made either by opening an existing data set via its name
 * within a container via a public constructor (typically when reading),
 *  or via \c H5mdDataSetBuilder<ValueType> to make a new data set
 *  (typically) when writing). The latter is implemented via a
 * private constructor.
 *
 *  \tparam ValueType Type of values stored in data set (not limited to base primitives).
 */
template<typename ValueType>
class H5mdDataSetBase
{
public:
    // Only make constructors with direct hid_t handles public to our builder classes.
    friend class H5mdDataSetBuilder<ValueType>;

    //! \brief Constructor to open an existing data set (called \p name, in \p container).
    H5mdDataSetBase(hid_t container, const char* name);

    //! \brief Destructor.
    ~H5mdDataSetBase() noexcept;
    //! \brief Copy constructor explicitly deleted.
    H5mdDataSetBase(const H5mdDataSetBase&) noexcept = delete;
    //! \brief Copy assignment explicitly deleted.
    H5mdDataSetBase& operator=(const H5mdDataSetBase&) = delete;
    //! \brief Move constructor.
    H5mdDataSetBase(H5mdDataSetBase&&) noexcept;
    //! \brief Move assignment explicitly deleted (not supported by sg::scope_guard).
    H5mdDataSetBase& operator=(H5mdDataSetBase&&) noexcept = delete;

    /*! \brief Return the HDF5 handle of the data set.
     *
     *  \warning Must not be closed by the caller.
     *
     *  \deprecated This handle will not be exposed once I/O is implemented as class methods.
     */
    hid_t id() const;

    /*! \brief Return the HDF5 handle of the data type for the managed data set.
     *
     *  \warning Must not be closed by the caller.
     *
     *  \deprecated This handle will not be exposed once I/O is implemented as class methods.
     */
    hid_t dataType() const;

    /*! \brief Return the HDF5 handle of the native data type for the managed data set.
     *
     *  \warning Must not be closed by the caller.
     *
     *  \deprecated This handle will not be exposed once I/O is implemented as class methods.
     */
    hid_t nativeDataType() const;

    /*! \brief Return the dimensions of the data set.
     *
     *  \throws H5mdError if the dimensions cannot be read.
     */
    DataSetDims dims() const;

private:
    /*! \brief Constructor to manage a data set with given \p dataSetHandle.
     *
     * \note The responsibility of closing this handle is taken over by this class:
     * it should not be closed by the caller.
     */
    explicit H5mdDataSetBase(hid_t dataSetHandle);

    //!< Handle to the managed data set.
    const hid_t dataSet_;

    /*! \brief Scope guard which closes the data set handle when the destructor is called.
     *
     * \note Must be declared directly below and thus initialized just after the guarded
     * handle. If not the guard may not be active if a subsequent initialization fails,
     * resulting in a memory leak.
     */
    H5mdGuard dataSetGuard_ = sg::make_scope_guard(H5mdCloser(dataSet_, H5Dclose));

    //!< Handle to the data type of the managed data set used when writing the set.
    const hid_t dataType_;

    /*! \brief Scope guard which closes the data type handle when the destructor is called.
     *
     * \note Must be declared directly below and thus initialized just after the guarded
     * handle. If not the guard may not be active if a subsequent initialization fails,
     * resulting in a memory leak.
     */
    H5mdGuard dataTypeGuard_ = sg::make_scope_guard(H5mdCloser(dataType_, H5Tclose));

    /*! \brief Handle to the immutable native data type of the managed data set.
     *
     * This may be different from dataType_ when opening a data set from a file which was created
     * on another machine. When reading data from a data set this is used along with dataType_ to
     * convert the data into the correct type of the current compiler target (i.e. native).
     */
    const hid_t nativeDataType_;

    /*! \brief Scope guard which closes the native data type handle when the destructor is called.
     *
     * \note Must be declared directly below and thus initialized just after the guarded
     * handle. If not the guard may not be active if a subsequent initialization fails,
     * resulting in a memory leak.
     */
    H5mdGuard nativeDataTypeGuard_ = sg::make_scope_guard(H5mdCloser(nativeDataType_, H5Tclose));

    //!< Number of dimensions of data set.
    const hsize_t numDims_;
};

extern template class H5mdDataSetBase<int32_t>;

extern template class H5mdDataSetBase<int64_t>;

extern template class H5mdDataSetBase<float>;

extern template class H5mdDataSetBase<double>;

extern template class H5mdDataSetBase<BasicVector<float>>;

extern template class H5mdDataSetBase<BasicVector<double>>;

extern template class H5mdDataSetBase<std::string>;

} // namespace gmx

#endif // GMX_FILEIO_H5MD_DATASETBASE_H
