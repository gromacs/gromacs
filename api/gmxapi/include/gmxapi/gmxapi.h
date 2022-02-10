/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
/*! \file
 * \brief Header for public GROMACS C++ API
 *
 * API clients include this header. It is intended to provide a minimal set of
 * declarations to allow independently implemented API clients to make compatible
 * references to gmxapi objects. Each client will still need to include additional
 * headers and to link against gmxapi, but clients don't need to be completely ABI
 * or API compatible with each other. Clients should use the gmxapi versioning utilities
 * to check for compatibility before accessing members of an object passed by
 * another client.
 *
 * \ingroup gmxapi
 */
/*! \mainpage
 * The `gmxapi` library allows GROMACS extension and client code to interact with
 * GROMACS internals without compile-time dependence on the GROMACS source code.
 * It is sufficient to use the headers described here and to link against `libgmxapi`
 * CMake helpers make it
 * easy to compile your own code against an installed copy of GROMACS.
 *
 * Using the `Gromacs::gmxapi` CMake target, the headers described here can be
 * included from `gmxapi/...` and symbols in the `::gmxapi` C++ namespace can be
 * resolved in the library by link time.
 *
 * For API stability and unambiguous versioning, headers in the top-level `gmxapi/`
 * directory should not refer to the `gmx` namespace used by the core GROMACS library,
 * except for gromacsfwd.h, which consolidates forward declarations in the `gmx` namespace
 * required by the extension APIs in subdirectories like `gmxapi/md`.
 * \todo gromacsfwd.h probably shouldn't be in the top level either...
 * To reduce dependencies, headers should not depend on headers at a deeper level
 * than themselves, where versioning and compatibility guarantees are weaker
 * (possibly dependent on GROMACS versions) or feature-dependent.
 *
 * Refer to the <a href="modules.html">modules</a> section for a hierarchical
 * overview of the API documentation.
 */
/*!
 * \defgroup gmxapi gmxapi
 *
 * \brief Provide external access to an installed GROMACS binary through the gmxapi library.
 *
 * API client code primarily operates on Proxy objects, which can be copied,
 * serialized, and configured independently of execution context, data location,
 * or parallelism.
 *
 * When data needs to be accessed or moved, smart Handle objects are used.
 * DataHandle objects explicitly describe ownership and access characteristics,
 * and so are not copyable, but can have ownership transferred. DataHandle objects
 * should be held briefly, since the implementation can depend on the current
 * execution context and can affect performance.
 *
 * Computation Modules and Runners have Proxy objects whose lifetime is independent
 * of execution context, and Implementation objects whose lifetimes are necessarily
 * bounded by the execution Context in which they run.
 *
 * Proxy and Implementation classes use Interfaces defined by the API, which can also be
 * subclassed in client code to extend the library functionality.
 *
 * The Implementation class supporting a given proxy object is likely to change between
 * successive handle requests, using a "state" behavioral design pattern.
 *
 * API classes participating in behavioral protocols beyond RAII construction and
 * destruction use implementation inheritance or separate manager classes to
 * ensure proper behavior by classes
 * implementing the relevant Interface. E.g. if the API specifies that an
 * Interface requires calls to object->initialize() to be followed by calls to
 * object->deinitialize(), the API library provides mechanisms to guarantee
 * that clients of the interface execute an appropriate state machine.
 *
 * The runner for a workflow cannot be fully instantiated until the execution context
 * is active, but the execution context cannot be configured without some knowledge
 * of the work to be performed. Thus, initialization requires several stages in
 * which the Context is configured with knowledge of the work to be performed,
 * a Session is launched and can provide resources to the computational
 * components at run time.
 *
 * For example:
 *
 * A client may create a MDProxy to pass to a RunnerProxy to pass to a Context, optionally
 * keeping references to the proxies or not.
 * The client asks the Context for a Session.
 * Context asks RunnerProxy for a RunnerBuilder, which asks MDProxy for an MD Engine builder.
 * Each builder (and product) gets a handle to the Context.
 * MDEngineBuilder produces and returns a handle to a new MDEngine.
 * RunnerBuilder produces a Runner with a handle to the MDEngine and returns runner handle to session.
 * Session makes API calls to runner handle until it no longer needs it and releases it.
 * Releasing the handle can provide notification to the runner.
 * When no references remain (as via handles), the runner is destroyed (same for MDEngine).
 * In the case of complex shutdown, the Context holds the last references to the runners and manages shutdown.
 *
 * For maximal compatibility with other libgmxapi clients (such as third-party
 * Python modules), client code should use the wrappers and protocols in the
 * gmxapi.h header. Note that there is a separate CMake target to build the full
 * developer documentation for gmxapi: `gmxapi-cppdocs-dev`.
 *
 * # Possible Example
 *
 * Configure a system and work specification from a TPR file, configure locally-detected computing resources,
 * and run the workflow.
 *
 *      std::string filename;
 *      // ...
 *      auto system = gmxapi::fromTprFile(filename);
 *      auto session = system.launch();
 *      auto status = session.run();
 *      return status.success();
 *
 * Load a TPR file, extract some data, modify the workflow, then run.
 *
 *
 *      auto system = gmxapi::fromTprFile(filename);
 *      // Load some custom code:
 *      auto potential = myplugin::MyPotential::newSpec();
 *      // ...
 *      auto atoms = system.atoms();
 *      {
 *          // Acquire read handle to global data in the current thread and/or MPI rank.
 *          const auto scopedPositionsHandle = gmxapi::extractLocally(atoms.positions, gmxapi::readOnly());
 *          for (auto&& position : scopedPositionsHandle)
 *          {
 *              // do something
 *          }
 *          // Release the read handle in the current code block (a shared reference may still be held elsewhere)
 *      }
 *      system.addPotential(potential);
 *      auto session = system.launch();
 *      auto status = session.run();
 *      return status.success();
 *
 * \internal
 * To extend the API may require GROMACS library
 * headers and possibly linking against `libgromacs`.
 * Where possible, though, gmxapi public interfaces never require the instantiation
 * of a gmx object, but instead take them as parameters. The (private) API implementation
 * code then uses binding protocols defined in the library API to connect the
 * GROMACS library objects (or interfaces), generally in an exchange such as the
 * following.
 *
 * gmxapi client code:
 *
 *     MyExtension local_object;
 *     gmxapi::SomeThinWrapper wrapper_object = gmxapi::some_function();
 *     local_object.bind(wrapper_object);
 *     gmxapi::some_other_function(wrapper_object);
 *
 * `MyExtension` implementation:
 *
 *     class MyExtension : public gmxapi::SomeExtensionInterface
 *     {
 *          void bind(gmxapi::SomeThinWrapper wrapper)
 *          {
 *              std::shared_ptr<gmxapi::SomeWrappedInterface> api_object = wrapper.getSpec();
 *              api_object->register(this);
 *          };
 *          //...
 *     };
 *
 * gmxapi protocol:
 *
 *     void gmxapi::SomeWrappedInterface::register(gmxapi::SomeExtensionInterface* client_object)
 *     {
 *          gmx::InterestingCoreLibraryClass* core_object = this->non_public_interface_function();
 *          gmx::SomeLibraryClass* library_object = client_object->required_protocol_function();
 *          core_object->do_something(library_object);
 *     };
 *
 *
 * Refer to GMXAPI developer docs for the protocols that map gmxapi interfaces to
 * GROMACS library interfaces.
 * Refer to the GROMACS developer
 * documentation for details on library interfaces forward-declared in this header.
 *
 * For binding Molecular Dynamics modules into a simulation, the thin wrapper class is gmxapi::MDHolder.
 * It provides the MD work specification that can provide and interface that an IRestraintPotential
 * can bind to.
 *
 * The gmxpy Python module provides a method gmx.core.MD().add_force(force) that immediately calls
 * force.bind(mdholder), where mdholder is an object with Python bindings to gmxapi::MDHolder. The
 * Python interface for a restraint must have a bind() method that takes the thin wrapper. The
 * C++ extension implementing the restraint unpacks the wrapper and provides a gmxapi::MDModule to
 * the gmxapi::MDWorkSpec. The MDModule provided must be able to produce a gmx::IRestraintPotential
 * when called upon later during simulation launch.
 *
 * external Restraint implementation:
 *
 *      class MyRestraintModule : public gmxapi::MDModule
 *      {
 *          // not a specified interface at the C++ level. Shown for Python module implementation.
 *          void bind(gmxapi::MDHolder wrapper)
 *          {
 *              auto workSpec = wrapper->getSpec();
 *              shared_ptr<gmxapi::MDModule> module = this->getModule();
 *              workSpec->addModule(module);
 *          }
 *      };
 *
 * When the simulation is launched, the work specification is passed to the runner, which binds to
 * the restraint module as follows. The execution session is launched with access to the work specification.
 *
 * Possible client code
 *
 *      gmxapi::System system;
 *      // System can provide an interface with which
 *      //...
 *      auto session = gmxapi::Session::create(system);
 *      session.run();
 *
 * Possible API implementation
 *
 *      static std::unique_ptr<gmxapi::Session> gmxapi::Session::create(shared_ptr<gmxapi::System> system)
 *      {
 *          auto spec = system->getWorkSpec();
 *          for (auto&& module : spec->getModules())
 *          {
 *              if (providesRestraint(module))
 *              {
 *                  // gmxapi restraint management protocol
 *                  runner->setRestraint(module);
 *              }
 *          }
 *
 *      }
 *
 *      // gmxapi restraint management protocol
 *      void gmxapi::Session::setRestraint(std::shared_ptr<gmxapi::MDModule> module)
 *      {
 *          auto runner = impl_->getRunner();
 *          auto restraint = module->getRestraint();
 *          runner->addPullPotential(restraint, module->name());
 *      }
 *
 */

#ifndef GMXAPI_H
#define GMXAPI_H

#include <memory>
#include <string>

/*! \brief Contains the external C++ Gromacs API.
 *
 * High-level interfaces for client code is provided in the gmxapi namespace.
 *
 * \ingroup gmxapi
 */
namespace gmxapi
{

// Forward declarations for other gmxapi classes.
class MDWorkSpec;
class MDModule;

/*!
 * \brief API name for MDHolder struct.
 *
 * Any object that includes this header (translation unit) will have the const cstring embedded as
 * defined here. The char array is not a global symbol.
 *
 * \todo Consider if this is in the right place.
 *
 * This header defines the API, but we may want to assure ABI compatibility. On the
 * other hand, we may want to embed the API version in the structure itself and leave
 * the higher-level name more general.
 */
static constexpr const char MDHolder_Name[] = "__GMXAPI_MDHolder_v1__";

/*!
 * \brief Minimal holder for exchanging MD specifications.
 *
 * The interface is minimal in the hopes of backwards and forwards compatibility
 * across build systems. This type specification can be embedded in any client code
 * so that arbitrary client code has a robust way to exchange data, each depending
 * only on the gmxapi library and not on each other.
 *
 * Objects of this type are intended to serve as wrappers used briefly to establish
 * shared ownership of a MDWorkSpec between binary objects.
 *
 * A sufficiently simple struct can be defined for maximum forward/backward compatibility
 * and given a name that can be used to uniquely identify Python capsules or other data
 * members that provide a pointer to such an object for a highly compatible interface.
 *
 * \todo Consider whether/how we might use a wrapper like this to enable a C API.
 * We may be able to replace this public type specification with an opaque declaration
 * and a set of free functions, but the change would need to be tested in the
 * context of the Python bindings and use cases at https://github.com/kassonlab/gmxapi
 *
 * Example:
 *
 *     py::class_< ::gmxapi::MDModule, std::shared_ptr<::gmxapi::MDModule> >
 *         gmxapi_mdmodule(m, "MDModule", py::module_local());
 *     gmxapi_mdmodule.def(
 *       "bind",
 *       [](std::shared_ptr<TestModule> self, py::object object){
 *           if (PyCapsule_IsValid(object.ptr(), gmxapi::MDHolder::api_name))
 *           {
 *               auto holder = (gmxapi::MDHolder*) PyCapsule_GetPointer(
 *                   object.ptr(),
 *                   gmxapi::MDHolder::api_name);
 *               auto spec = holder->getSpec();
 *               std::cout << self->name() << " received " << holder->name();
 *               std::cout << " containing spec of size ";
 *               std::cout << spec->getModules().size();
 *               std::cout << std::endl;
 *               spec->addModule(self);
 *           }
 *           else
 *           {
 *               throw gmxapi::ProtocolError("MDModule bind method requires properly named PyCapsule input.");
 *           }
 *       }
 *     );
 *
 *     py::class_<System, std::shared_ptr<System> > system(m, "MDSystem");
 *     system.def(py::init(), "A blank system object is possible, but not useful. Use a helper function.");
 *     system.def("launch",
 *                 [](System* system, std::shared_ptr<Context> context)
 *                 {
 *                     return system->launch(context);
 *                 },
 *                 "Launch the configured workflow in the provided context.");
 *
 *     system.def(
 *         "add_mdmodule",
 *         [](System* system, py::object force_object){
 *             // If force_object has a bind method, give it a PyCapsule with a pointer
 *             // to our C++ object.
 *             if (py::hasattr(force_object, "bind"))
 *             {
 *                 auto spec = system->getSpec();
 *                 auto holder = new gmxapi::MDHolder(spec);
 *                 holder->name_ = "pygmx holder";
 *                 auto deleter = [](PyObject *o) {
 *                     if (PyCapsule_IsValid(o, gmxapi::MDHolder_Name))
 *                     {
 *                         auto holder_ptr = (gmxapi::MDHolder *) PyCapsule_GetPointer(o, gmxapi::MDHolder_Name);
 *                         delete holder_ptr;
 *                     };
 *                 };
 *                 auto capsule = py::capsule(holder,
 *                                            gmxapi::MDHolder_Name,
 *                                            deleter);
 *                 py::object bind = force_object.attr("bind");
 *                 // py::capsule does not have bindings and does not implicitly convert to py::object
 *                 py::object obj = capsule;
 *                 bind(obj);
 *                 std::cout << "Work specification now has " << spec->getModules().size() << " modules." << std::endl;
 *             }
 *             else
 *             {
 *                 // Note: Exception behavior is likely to change.
 *                 // Ref: https://github.com/kassonlab/gmxapi/issues/125
 *                 throw PyExc_RuntimeError;
 *             }
 *         },
 *         "Set a restraint potential for the system.");
 *
 * \ingroup gmxapi_md
 */
class MDHolder
{
public:
    /*!
     * \brief Declare the schema for extra checks before casting.
     *
     * MDHolder exists for cases where API objects can only be passed by
     * casting pointers. api_name is suitable as a descriptor of the schema
     * of the object pointed to. See, for instance, usage of the PyCapsule_IsValid
     * Python C API function.
     */
    static const char* api_name;

    MDHolder();

    /*!
     * \brief For convenience and logging, give the object an identifying string.
     *
     * \param name
     */
    explicit MDHolder(std::string name);

    /*!
     * \brief Wrap a Molecular Dynamics work specification.
     *
     * The container allows portable specification of MD work to be performed.
     * It is used when setting up and then launching the simulation.
     *
     * \param spec references a container with interfaces for client and library APIs
     *
     * Example:
     *
     *     # With `system` as a gmxapi::System object
     *     auto spec = system->getSpec();
     *     auto holder = std::make_unique<gmxapi::MDHolder>(spec);
     *
     * A PyCapsule object with the name given by gmxapi::MDHolder_Name is assumed to
     * contain a pointer to an MDHolder and to have an appropriate deleter attached.
     *
     * Example:
     *
     *     auto deleter = [](PyObject *o) {
     *         if (PyCapsule_IsValid(o, gmxapi::MDHolder_Name))
     *         {
     *             auto holder_ptr = (gmxapi::MDHolder *) PyCapsule_GetPointer(o,
     * gmxapi::MDHolder_Name); delete holder_ptr;
     *         };
     *     };
     *     # With pybind11 PyCapsule bindings:
     *     auto capsule = py::capsule(holder,
     *                                gmxapi::MDHolder_Name,
     *                                deleter);
     *
     * The gmxapi Python package gives modules a chance to associate themselves with a
     * gmxapi::System object by passing such a PyCapsule to its `bind` method, if implemented.
     *
     * Such a bind method could be implemented as follows. Assume object.ptr() returns a
     * `PyObject*`
     *
     * Example:
     *
     *    PyObject* capsule = object.ptr();
     *    if (PyCapsule_IsValid(capsule, gmxapi::MDHolder::api_name))
     *    {
     *        auto holder = static_cast<gmxapi::MDHolder*>(PyCapsule_GetPointer(capsule,
     *            gmxapi::MDHolder::api_name));
     *        auto workSpec = holder->getSpec();
     *        workSpec->addModule(module);
     *    }
     *    else
     *    {
     *        throw gmxapi::ProtocolError("bind method requires a python capsule as input");
     *    }
     */
    explicit MDHolder(std::shared_ptr<MDWorkSpec> spec);

    /*!
     * \brief Get client-provided name.
     *
     * \return Name as string.
     */
    std::string name() const;

    /*! \brief Instance name.
     */
    std::string name_{};

    /// \{
    /*!
     * \brief Get the wrapped work specification
     * \return shared ownership of the api object.
     */
    std::shared_ptr<MDWorkSpec> getSpec();
    /*!
     * \brief Get the wrapped work specification
     * \return smart pointer to const object
     */
    std::shared_ptr<const MDWorkSpec> getSpec() const;
    /// \}
private:
    /*! \cond internal
     * \brief private implementation class
     */
    class Impl;
    /// \brief opaque pointer to implementation
    std::shared_ptr<Impl> impl_{ nullptr };
    /*! \endcond */
};

/*!
 * \brief Label the types recognized by gmxapi.
 *
 * Provide an enumeration to aid in translating data between languages, APIs,
 * and storage formats.
 *
 * \todo The spec should explicitly map these to types in APIs already used.
 * e.g. MPI, Python, numpy, GROMACS, JSON, etc.
 * \todo Actually check the size of the types.
 *
 * \see https://gitlab.com/gromacs/gromacs/-/issues/2993 for discussion.
 */
enum class GmxapiType
{
    NULLTYPE, //! Reserved
    MAP,      //! Mapping of key name (string) to a value of some MdParamType
    BOOL,     //! Boolean logical type
    INT64,    //! 64-bit integer type
    FLOAT64,  //! 64-bit float type
    STRING,   //! string with metadata
    NDARRAY,  //! multi-dimensional array with metadata
};
} // end namespace gmxapi

#endif // header guard
