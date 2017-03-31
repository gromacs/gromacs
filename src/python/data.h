/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
/// \cond
/*! \internal \file
 * \brief Declares trajectory data structures for export.
 *
 * \ingroup module_python
 */
#ifndef PYGMX_DATA_H
#define PYGMX_DATA_H

#include <memory>
#include <stdexcept>

#include "gromacs/options/options.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/runner.h"

namespace gmx
{

namespace pyapi
{
using std::shared_ptr;
using std::unique_ptr;
using std::make_shared;

/// Wrapper for flat data structure
/*! This is a temporary shim to experiment with how to manage multidimensional
 * data of arbitrary size in ways that are friendly to old code, new code, and STL.
 * \internal \ingroup module_python
 */
template<typename Scalar, size_t D>
class TrajDataArray
{
    public:
        /*
           // Create empty container
           TrajDataArray() :
            data_{}
           {};
         */
        /// An empty constructor call doesn't make sense
        TrajDataArray() = delete;
        /// We don't have good semantics for copy construction yet
        TrajDataArray(const TrajDataArray &)                  = delete;
        /// Disalloy copy assignment
        const TrajDataArray &operator=(const TrajDataArray &) = delete;

        /// Allocate space for an NxD array
        /*! \internal \ingroup module_python
         * \param N number of rows of width D to allocate
         */
        TrajDataArray(size_t N) :
            data_(N*D),
            N_(N)
        {
        };

        /*! \brief Copy from raw data pointer.
         *
         * The caller must tell this constructor the number of rows of D columns.
         * Since TrajDataArray has no influence over the lifetime of the source
         * data, we have to copy it.
         * \param data_src rvec* or similar raw pointer to unraveled array of
         * dimension known to the caller.
         * \param N number of rows in the array
         */
        TrajDataArray(Scalar* data_src, size_t N) :
            data_(N*D),
            N_(N)
        {
            std::copy(data_src, data_src + N*D, data_.begin());
        };

        // Move constructor
        //TrajDataArray(TrajDataArray&& )

        /// Simple destructor.
        ~TrajDataArray() { };

        /// Get width of data
        /*! \return number of columns (dimensionality for arrays of vectors)
         */
        size_t dim() const { return D; };
        /// Get number of elements
        /*! \return number of rows (number of elements for arrays of vectors)
         */
        size_t N() const { return N_; };

        /*! \brief Get raw pointer to managed data
         *
         * \return the pointer to the underlying contiguous data storage
         * This is very bad if the caller tries to delete the result. This is a
         * temporary workaround for the way we are using the python buffer protocol
         * in this iteration.
         */
        Scalar* data() { return data_.data(); };

        /*! \brief indexing operator
         *
         * \param i (first) index in the two-dimensional array
         * \return a copy (not a reference) of row i as a std::vector of length D
         * TODO: Not sure why this isn't a std::array. span or array view could be
         * even better. But first some thought would need to go into the expected
         * readable/writable behavior. Note that the return value *will be copied*
         * when returned to Python unless shared via a shared_ptr or transfered via
         * unique_ptr to a class that Python is aware of (no casting necessary).
         */
        std::vector<Scalar> operator[](const size_t i) const
        {
            if (i < N_)
            {
                return std::vector<Scalar>(&data_[i*D], &data_[i*D] + D);
            }
            else
            {
                throw std::out_of_range("bad index value to Trajectory data");
            }
        }

    private:
        /// Flattened array of data
        std::vector<Scalar> data_;
        /// Actual dimensions are N_ x D
        const size_t        N_;
};

template class TrajDataArray<real, 3>;
/// Provide an easier named alias.
using Data3 = TrajDataArray<real, 3>;

/*! \brief scoped enum to help identify requested fields and types.
 *
 * The enum value can be used as a light-weight way to select requested
 * trajectory data fields.
 *
 * The enum type can be used to specialize templates or for overload
 * resolution. Other enum types can be added at the cost of additional
 * such specializations. To be useful, each enum type so used should
 * be associated with data of a consistent structure.
 * \internal \ingroup module_python
 */
enum class trjvectorfield
{
    POSITION,
    VELOCITY,
    FORCE,
    ACCELERATION,
    MOMENTUM
};

/// Set up a compile-time mapping of fields to types
template<class T> struct t_field {
    /// specializations of t_field<fieldname> can specify an expected data type.
    typedef void type;
};

/// trjvectorfield have type Data3
template<> struct t_field<trjvectorfield> {
    /// Data fields enumerated in trjvectorfield are Nx3 arrays.
    typedef Data3 type;
};


/*! \brief API data object handle
 *
 * Proxy object to maintain a binding to a Gromacs data object.
 * Objects of this type are returned by getters of other API objects,
 * which construct handles by providing a means to keep data alive
 * and to retrieve it.
 *
 * Must be available for various returned data types for proper type
 * checking. Concrete classes implemented according to circumstance.
 * \internal \ingroup module_python
 */
template<class T>
class PyDataHandle
{
    public:
        /// TODO: Empty handle doesn't really make sense...
        PyDataHandle() = default;
        /// Handles should probably not be copied, but moving is probably fine.
        PyDataHandle(const PyDataHandle &)            = delete;
        /// Handles should probably not be copied, but moving is probably fine.
        PyDataHandle &operator=(const PyDataHandle &) = delete;

        virtual ~PyDataHandle() = default;

        // TODO: Derived classes must implement a read-only fetch.
        //virtual const shared_ptr<T> fetch_data() const = 0;
        // TODO: Derived classes do not need to provide a non-const fetch,
        // but if they do they should not do so for read-only policy.
        /*! \brief Retrieve data to be exported outside of the API.
         *
         * To operate on Gromacs data outside of the API may require
         * data copies or communication that is not required within
         * library calls, so exporting data to external calling code
         * may incur performance penalties.
         * \return shared ownership of raw array data.
         */
        virtual shared_ptr<T> fetch_data() { return nullptr; };
};

template class PyDataHandle<Data3>;
/// Provide an easy name for a common specialization.
using Data3Handle = PyDataHandle<Data3>;


/*! \brief Implement a handle. Trivial for local data.
 *
 * \internal \ingroup module_python
 */
class LocalTrajDataHandle : public Data3Handle
{
    public:
        /// A handle is meaningless on its own. No empty handles.
        LocalTrajDataHandle() = delete;
        /// Construct handle from a shareable data object.
        explicit LocalTrajDataHandle(shared_ptr<Data3> data);
        /// clean up
        virtual ~LocalTrajDataHandle();

        /// Make data available to caller
        /*! In this case, it is trivial.
         * \return shared_ptr to data already available locally
         */
        shared_ptr<Data3> fetch_data();
    private:
        /// Keep the data alive and accessible. Trivial in this case.
        shared_ptr<Data3> data_;
};

};     //end namespace pyapi
};     // end namespace gmx
#endif // header guard
/// \endcond
