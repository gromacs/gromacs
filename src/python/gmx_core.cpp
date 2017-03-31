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
/*! \internal \file
 * \brief Exports Python bindings for gmx.core module.
 *
 * \ingroup module_python
 */

#include "gmxpre.h"

#include "gmx_core.h"

// TODO: Tell doxygen to suggest quotes instead of angle brackets for headers.
// TODO: remove extraneous headers
#include <sstream>

#include "pybind11/pybind11.h"

#include "gromacs/options/filenameoptionmanager.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/runner.h"
#include "gromacs/trajectoryanalysis/modules/caching.h"
#include "python/data.h"

#include "options.h"
#include "runner.h"
#include "trajectory.h"

using std::shared_ptr;
using std::unique_ptr;
using std::make_shared;

namespace gmx
{
namespace pyapi
{


} // end namespace pyapi
} // end namespace gmx

// Export Python module.

namespace py = pybind11;

const char* const name = "core";                     ///< used to set __name__
// pybind11 uses const char* objects for docstrings. C++ raw literals can be used.
const char* const docstring = "Gromacs core module"; ///< used to set __doc__

/*! \brief Export gmx.core Python module in shared object file.
 *
 * One goal of these bindings is to declare a buffer type suitable for numpy Nx3 array output.
 * If we want to pass access but not ownership to Python, we need to make
 * sure we can allow a C++ shared pointer ref count to be increased.
 * The buffer protocol requires that the exporter (this code) keeps the
 * memory valid for the exported view until all consumers are done and
 * the PyBuffer_Release(buffer *view) is issued. I'm not sure, but I assume
 * pybind11 manages that for us by holding a shared_ptr to this. However, there
 * seem to be subtleties I overlooked in some testing, so this warrants
 * further investigation.
 *
 * To be safe (and in general) we can use the keep_alive<>() call policy
 * and return_value_policy annotations.
 * The keep_alive<Nurse, Patient>() call policy keeps Patient alive as long
 * as Nurse is alive. Indices for Nurse or Patient are 0 for the return value
 * of the annotated function
 * and higher numbers for arguments. For member functions, index 1 is used
 * for the this* object.
 *
 * The pybind11 documentation notes "For functions returning smart pointers, it is not necessary to specify a return value policy."
 * and
 * "It is even possible to completely avoid copy operations with Python expressions like np.array(matrix_instance, copy = False)"
 * \internal \ingroup module_python
 */
PYBIND11_PLUGIN(core) {
    using namespace gmx::pyapi;

    // Instantiate the module
    py::module m(name, docstring);

    // Export options class
    py::class_< PyOptions, std::shared_ptr<PyOptions> >(m, "Options")
        .def(
            py::init<const std::string>(),
            py::arg("filename")
            )
    ;

    // Export data classes


    // Export buffer class that exposes trajectory data arrays
    py::class_< Data3,
                std::shared_ptr<Data3 >
                >(m, "TrajData3", py::buffer_protocol())
    // A buffer interface exported to Python.
    // I'm not sure that def_buffer implies return_value_policy::reference_internal,
    // which implies keep_alive<0,1> to make sure that C++ will keep the
    // object alive as long as the return value is alive, def_buffer does
    // not take return_value_policy arguments. It would be nice if we could
    // instead use the support for Eigen / numpy compatibility.
    // Note that the bindings for std::vector use reference_internal for __getitem__
    // and keep_alive<0,1> for __iter__ but nothing for the def_buffer interface to
    // std::vector. However, it copies data from the py::buffer in the vector_buffer __init__.
    // We should probably use Eigen for data on the C++ side and export directly
    // to numpy arrays instead.
        .def_buffer(
            [](Data3 &data) -> py::buffer_info
            {
                return py::buffer_info(
                        data.data(),                                /* Pointer to buffer */
                        sizeof(real),                               /* Size of one scalar */
                        py::format_descriptor<real>::format(),      /* Python struct-style format descriptor */
                        2,                                          /* Number of dimensions */
                        { data.N(), data.dim() },                   /* Python buffer dimensions */
                        { sizeof(real) * data.dim(), sizeof(real) } /* Strides (in bytes) for each index in C++ */
                        );
            }
            )
    // Accept buffers from Python.
    // We should whether and how to safely perform a no-copy construction from a
    // buffer when TrajDataArray can have multiple references on the C++ side.
    // If the Python buffer views are all closed and there are no more
    // Python references to the object, then any remaining C++ references
    // to the object will have their data become invalid if the buffer object
    // is allowed to be released. I need to look into this more, but if pybind11 doesn't
    // already do it, we can keep the source of the buffer alive by using
    // keep_alive<1,3> in the init method to keep the py::buffer argument
    // alive as long as the TrajDataArray. return_value_policy::reference
    // prevents Python from taking ownership (C++ manages the lifetime).
    // That may interfere with our ability to release the global interpreter
    // lock. Instead, we should probably only attempt non-temporary objects
    // with non-copy construction on memory that is already owned by the API and
    // not the Python interpreter or retain handles to Python objects (providing
    // a buffer interface) that are only exposed in controlled and limited scope.
    // If we want to set
    // data in TrajDataArray objects with minimal copies, we allocate our own
    // memory and implement
    // element access methods to write directly from Python to the managed
    // array.
    // TO DO: only accept dense numpy arrays with array_t arguments.
    // py::array_t<real, py::array::c_style | py::array::forcecast> data
        .def("__init__",
             [](Data3 &data, py::buffer b)
             {
                 /* Request a buffer descriptor from Python */
                 py::buffer_info info = b.request();

                 /* Some sanity checks ... */
                 if (info.format != py::format_descriptor<real>::format())
                 {
                     throw std::runtime_error("Incompatible format: expected a array of type real!");
                 }
                 ;
                 if (info.ndim != 2 || info.shape[0] != 3)
                 {
                     throw std::runtime_error("Incompatible buffer dimension!");
                 }
                 ;

                 // Construct in place
                 // It is important that the reference count of the buffer object b
                 // should be incremented to prevent Python garbage collection from
                 // deallocating the memory in the TrajDataArray object. I assume
                 // pybind11 takes care of that.
                 new (&data)Data3(static_cast<real *>(info.ptr), info.shape[0]);
             },
             py::keep_alive<1, 3>() // keep py::buffer b alive while *this is alive
             )
    // Inspect...
        .def("__repr__", [](const Data3 &t)
             {
                 std::stringstream repr;
                 //std::string repr{"tah dah!"};//
                 repr << t.N() << "x" << t.dim() << " array of trajectory data of type 'real'\n";
                 // Ugh...
                 for (size_t i = 0; i < t.N(); ++i)
                 {
                     repr << t[i][0] << "\t" << t[i][1] << "\t" << t[i][2] << "\n";
                 }
                 return repr.str();
             }
             )
        .def_property_readonly("N", &Data3::N, "number of elements")
    /* Needs better TrajDataArray
       .def("__iter__", [](const Data3 &s)
         {
             return py::make_iterator(s.begin(), s.end());
         },
         py::keep_alive<0, 1>() // Essential: keep object alive while iterator exists
        )
     */
    // Use generator or make a list instead...
    //.def("__getitem__", [](const Data3 &s, py::slice slice) -> Data3* {})
    ;

    // Export an API object to serve as a light-weight handle.
    py::class_< Data3Handle > (m, "GmxData3Base");
    py::class_< LocalTrajDataHandle, Data3Handle > (m, "GmxData3")
        .def("extract", &LocalTrajDataHandle::fetch_data, "Extract API object to Python interpreter");
    ;
    // Export trajectory frame class
    py::class_< PyTrajectoryFrame, shared_ptr<PyTrajectoryFrame> > (m, "Frame")
    // For some reason this doesn't work right if PyTrajectoryFrame::x returns
    // a unique_ptr instead of a shared_ptr.
    // Note that since TrajDataArray is a shim for the unmanaged arrays in t_trxframe
    // and since this result is not cached (yet), a TrajDataArray is
    // constructed and destroyed each time the property is accessed.
        .def_property_readonly("position",
                               [](PyTrajectoryFrame &f) -> std::unique_ptr<Data3Handle>
                               {
                                   return f.get_read_handle<trjvectorfield>(trjvectorfield::POSITION);
                               },
                               "API handle to positions");

    //
    // Export module classes
    //

    // Base class so Python and C++ can understand each other better.
    py::class_< gmx::TrajectoryAnalysisModule,
                shared_ptr<gmx::TrajectoryAnalysisModule>
                >(m, "TafModuleAbstractBase");
    // Caching dummy module.
    // Default holder is std::unique_ptr, but we allow multiple handles to module.
    py::class_< gmx::trajectoryanalysis::CachingTafModule,
                shared_ptr<gmx::trajectoryanalysis::CachingTafModule>,
                gmx::TrajectoryAnalysisModule
                >(m, "CachingTafModule")
        .def(py::init())
    //.def("select", py:make_iterator...)
        .def("frame",
             [](const gmx::trajectoryanalysis::CachingTafModule &cache) -> shared_ptr<PyTrajectoryFrame>
             {
                 return std::make_shared<PyTrajectoryFrame>(cache.frame());
             },
             "Retrieve cached trajectory frame handle."
             );

    // Export runner class
    py::class_< PyRunner > runner(m, "TafRunner");
    // We shouldn't need a keep_alive<>() for the module we're attaching since
    // pybind knows how to handle shared_ptr and this object does not need to
    // survive after the associated module is done with it, but we need to
    // reconsider when the usage model changes for chained or more general modules.
    runner.def(py::init< shared_ptr<gmx::TrajectoryAnalysisModule> >())
        .def("initialize", &PyRunner::initialize, "handle options")
        .def("next", &PyRunner::next, "Advance the current frame one step.");


/*
    py::class_< gmx::Options >(m, "Options")
        .def(py::init<>()); // Need to figure out options passing...


    py::class_< gmx::Selection >(m, "Selection");
 */

    return m.ptr();
}
