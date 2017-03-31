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

#include "gmx_core.h"

// TODO: Tell doxygen to suggest quotes instead of angle brackets for headers.
#include "gromacs/trajectoryanalysis/runner.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/modules/caching.h"
#include "gromacs/options/options.h"
#include "gromacs/options/filenameoptionmanager.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/options/optionsassigner.h"
#include "gromacs/options/optionsvisitor.h"
#include "gromacs/trajectory/trajectoryframe.h"

//#include <iostream>
#include <sstream>

#include "pybind11/pybind11.h"

using std::shared_ptr;
using std::unique_ptr;
using std::make_shared;

namespace gmx
{
namespace pyapi
{

PyTrajectoryFrame::PyTrajectoryFrame(std::shared_ptr<t_trxframe> frame) :
    frame_ {frame}
{
}

PyTrajectoryFrame::PyTrajectoryFrame(const t_trxframe &frame) :
    frame_ {gmx::trajectory::trxframe_copy(frame)}
{
}

shared_ptr< TrajDataArray<real, 3> > PyTrajectoryFrame::x()
{
    return make_shared< TrajDataArray<real, 3> >(frame_->x[0], frame_->natoms);
    //return std::unique_ptr< TrajDataArray<real, 3> >(new TrajDataArray<real, 3> (frame_->x[0], frame_->natoms));
}

PyOptions::PyOptions() :
    options_ {}
{
    //print_options(*this);
}

PyOptions::PyOptions(std::string filename) :
    options_ {},
filename_ {
    filename
}
{}

PyOptions::~PyOptions()
{}

gmx::Options* PyOptions::data()
{
    // Return a pointer to the Options object
    return &options_;
}


bool PyOptions::parse()
{
    /// Helper class for RAII wrapping of gmx::OptionsAssigner
    class Assigner
    {
        public:
            Assigner(gmx::Options &options) : assigner_ {&options}
            {
                assigner_.start();
            };
            ~Assigner()
            {
                assigner_.finish();
            };

            bool startOption(const std::string &name)
            {
                assigner_.startOption(name.c_str());
                return true;
            };

            bool addSingleValue(const std::string &value)
            {
                try
                {
                    assigner_.appendValue(value.c_str());
                }
                catch (GromacsException &ex)
                {
                    assigner_.finishOption();
                    return false;
                }
                assigner_.finishOption();
                return true;
            };

        private:
            gmx::OptionsAssigner assigner_;
    };
    // use gmx::OptionsAssigner to set parsed options, but refer to
    // CommandLineParser().impl_->assigner for example helper functions.
    // Note that toOptionName just identifies command-line flags and strips
    // them into the form of an option name. We can do this with py dict...
    // E.g.
    //    assigner.start();
    //    // loop over args
    //    const char* const optionName = toOptionName(arg);
    //    // note, arg sets optionName = nullptr implies append to previous option
    //    try
    //        is_valid_name = assigner.tryStartOption(optionName);
    //        // handle invalid options...
    //    catch (UserInputError &ex)
    //        // badly formed option argument
    //    // If instead we are appending:
    //    assigner.appendValue(arg) // throws GromacsException
    //    // conclude single or multivalue options. Multivalue options receive
    //    // additional checking and get another chance to throw.
    //    assigner.finishOption() // throws UserInputError
    //    // end loop
    //    assigner.finish();
    //
    // In the longer term, the options object has a parser available after
    // interacting with the module runner, and can then
    // provide useful user help. It probably doesn't make sense to allow kwargs
    // in a generic and configurable Python class, but maybe the various
    // components can be independently configured or accessed statically to
    // build an argparse-like object. I suppose that's what the help generator
    // does...

    // Scope the lifetime of assigner
    {
        Assigner assigner {
            options_
        };
        // TrajectoryRunnerCommon names the filename option "f"

        const std::string name {
            "f"
        };

        try     // assign a new option
        {
            if (assigner.startOption(name.c_str()))
            {
                // assign (all of) the value(s) for name
                if (!assigner.addSingleValue(filename_.c_str()))
                {
                    throw(InvalidInputError("bad option value"));
                }
                // else value successfully appended
            }
            else
            {
                throw(InvalidInputError("bad option name"));
            };
        }
        catch (InvalidInputError &ex)
        {
            // InvalidInputError is thrown if option is not recognized or inappropriate (e.g. specified more than once)
            throw(ex);
        }
    }
    options_.finish();
    //print_options(*this);
    return true;
}

void PyOptions::view_traverse(gmx::OptionsVisitor &&visitor) const
{
    visitor.visitSection(options_.rootSection());
}

void print_options(const PyOptions &pyoptions)
{
    /// Provide a way to explore the PyOptions object
    /*! gmx::OptionsVisitor defines an interface for visitor classes.
     *  gmx::OptionsIterator provides a means to traverse an Options collection
     *  as a non-member from arbitrary calling code, rather than as a member function of the collection,
     *  which would be more like the Visitor pattern I'm used to.
     *  gmx::OptionsIterator decorates an Options object to provide
     *  acceptSection(OptionsVisitor*) and acceptOptions(OptionsVisitor*)
     *  so that a visitor object should have its visit methods
     *  called directly. E.g. Visitor().visitSection(options.rootSection()) calls
     *  OptionsIteratory(options.rootSection()).acceptSections(visitor) and again
     *  for acceptOptions(visitor). It is not documented whether OptionsIterator.acceptSections(visitor) is made recursive through the Visitor's implementation of visitSection.
     */
    // TODO: IOptionsContainer can only throw APIError and has no way to tell
    // the caller that an option is already defined. However, it is implemented
    // in gmx::Options, so the calling code can get there in a roundabout way.
    // TODO: There is no const version of gmx::OptionsVisitor
    class Visitor : public gmx::OptionsVisitor
    {
        virtual void visitSection(const gmx::OptionSectionInfo &section)
        {
            // note hierarchy...
            // const std::string name = section.name()
            gmx::OptionsIterator iterator(section);
            iterator.acceptSections(this);
            iterator.acceptOptions(this);
        }
        virtual void visitOption(const OptionInfo &option)
        {
            // Do something...
            const std::string name   = option.name();
            const bool        is_set = option.isSet();

            // Can't see values? OptionInfo::option() returns a AbstractOptionStorage& but is a protected function.
            // There does not appear to be a way to get the OptionType (template
            // parameter) settings object used in addOption. OptionType is
            // derived from AbstractOption, I think. Unless there is a default
            // value, there is no option value until the Assigner runs, which
            // operates on a full gmx::Options object. Then the value is owned
            // by the caller of IOptionsContainer.addOption()
            // Where does the reference in AbstractOption.store(T*) end up?
            // Options have a T* store_ that points to the storage defined
            // in the object passed to addOption().
            // OptionSectionImpl::Group::addOptionImpl(AbstractOption& settings) calls
            // settings.createStorage() to get a OptionSectionImpl::AbstractOptionStoragePointer object.
            // When the Assigner gets to appendValue, there is ultimately a
            // commitValues() template method that calls OptionStorageTemplate<T>::store_->append(value). Here, store_ was
            // initialized in the constructor and is a
            // std::unique_ptr<IOptionValueStore<T> >
            // IOptionValueStore<T> is a wrapper that is implemented by various
            // basic ValueStore types to
            // provide some methods and a wrapped ArrayRef<T>.
            // Assigner::Impl uses Options::impl_->rootSection_ to get a
            // internal::OptionSectionImpl which provides findOption() to get
            // an AbstractOptionStorage*. AbstractOptionsStorage->appendValue()
            // results in the values actually being set through virtual functions
            // in the inaccessible OptionStorage object.
            // There appears to be no way to get from either an Options or
            // OptionInfo object to the OptionStorageTemplate object that can
            // see the actual storage. I would need to implement an alternative
            // IOptionsContainer or cause the modules to provide some interface.
            // Also, the parsed option values appear temporarily in the machinery
            // but I think are gone after assignment completes.
            //
            // In short, if I want to see the values being passed now, I can
            // look at the raw memory for the storage destinations, the values
            // being handled by the Assigner, or add printf debugging into the
            // core options handling code...
        }
    };
    pyoptions.view_traverse(Visitor());
}


PyRunner::PyRunner(shared_ptr<gmx::TrajectoryAnalysisModule> module) :
    runner_(),
    module_(module)
{
    runner_.add_module(module_);
}

PyRunner::~PyRunner() {}

bool PyRunner::next()
{
    return runner_.next();
}

void PyRunner::initialize(PyOptions &options)
{
    //gmx::FileNameOptionManager filename_option_manager;
    //options.data()->addManager(&filename_option_manager);
    //print_options(options);
    runner_.register_options(*options.data());
    // parse options...
    if (!options.parse())
    {
        throw(InvalidInputError("could not parse"));
    }
    //print_options(options);
    options.data()->finish();
    //print_options(options);
    runner_.initialize(*options.data());
    ////print_options(options);
}

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

    // Export runner class
    py::class_< PyRunner > runner(m, "TafRunner");
    // We shouldn't need a keep_alive<>() for the module we're attaching since
    // pybind knows how to handle shared_ptr and this object does not need to
    // survive after the associated module is done with it, but we need to
    // reconsider when the usage model changes for chained or more general modules.
    runner.def(py::init< shared_ptr<gmx::TrajectoryAnalysisModule> >())
        .def("initialize", &PyRunner::initialize, "handle options")
        .def("next", &PyRunner::next, "Advance the current frame one step.");

    // Export module classes
    py::class_< gmx::TrajectoryAnalysisModule,
                shared_ptr<gmx::TrajectoryAnalysisModule>
                >(m, "TafModuleAbstractBase");

    // Export trajectory frame class
    py::class_< PyTrajectoryFrame, shared_ptr<PyTrajectoryFrame> > (m, "Frame")
    // For some reason this doesn't work right if PyTrajectoryFrame::x returns
    // a unique_ptr instead of a shared_ptr.
    // Note that since TrajDataArray is a shim for the unmanaged arrays in t_trxframe
    // and since this result is not cached, a TrajDataArray is
    // constructed and destroyed each time the property is accessed.
        .def_property_readonly("x", &PyTrajectoryFrame::x, "positions");

    // Export the caching dummy module.
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
             "Retrieve cached trajectory frame."
             );

    // Export options class
    py::class_< PyOptions, std::shared_ptr<PyOptions> >(m, "Options")
        .def(
            py::init<const std::string>(),
            py::arg("filename")
            )
    ;

    // Export buffer class that exposes trajectory data arrays
    py::class_< TrajDataArray<real, 3>,
                std::shared_ptr<TrajDataArray<real, 3> >
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
            [](TrajDataArray<real, 3> &data) -> py::buffer_info
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
             [](TrajDataArray<real, 3> &data, py::buffer b)
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
                 new (&data)TrajDataArray<real, 3>(static_cast<real *>(info.ptr), info.shape[0]);
             },
             py::keep_alive<1, 3>() // keep py::buffer b alive while *this is alive
             )
    // Inspect...
        .def("__repr__", [](const TrajDataArray<real, 3> &t)
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
        .def_property_readonly("N", &TrajDataArray<real, 3>::N, "number of elements")
    /* Needs better TrajDataArray
       .def("__iter__", [](const TrajDataArray<real, 3> &s)
         {
             return py::make_iterator(s.begin(), s.end());
         },
         py::keep_alive<0, 1>() // Essential: keep object alive while iterator exists
        )
     */
    // Use generator or make a list instead...
    //.def("__getitem__", [](const TrajDataArray<real, 3> &s, py::slice slice) -> TrajDataArray<real, 3>* {})
    ;
/*
    py::class_< gmx::Options >(m, "Options")
        .def(py::init<>()); // Need to figure out options passing...


    py::class_< gmx::Selection >(m, "Selection");
 */

    return m.ptr();
}
