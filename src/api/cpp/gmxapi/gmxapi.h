/*! \file
 * \brief Header for public Gromacs C++ API
 *
 * API clients include this header.
 *
 * \defgroup gmxapi Gromacs external C++ interface.
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
 * of the work to be performed. Thus, the generic runner API object participates in initializing
 * the execution context, and a particular context implementation can provide
 * specialized factory extensions to instantiate and initialize the concrete runner
 * implementation, which exists at the library API level and is not directly exposed to
 * the external API.
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
 *
 */

#ifndef GMXAPI_H
#define GMXAPI_H

#include <memory>

#include "gmxapi/exceptions.h"
#include "gmxapi/version.h"

/*! \brief Contains the external C++ Gromacs API.
 *
 * High-level interfaces for client code is provided in the gmxapi namespace.
 * The interface is semantically the same as the accompanying Python module,
 * though syntactically distinct to be idiomatically consistent with the language.
 *
 * A lower-level interface to implement or extend the API is in gmxapi::core.
 *
 * \ingroup gmxapi
 */
namespace gmxapi
{

/*! \brief Container for results of API operations.
 *
 * \internal
 * I'm leaning towards saying this should not be derivable, but that it
 * may contain one or more additional objects, perhaps including exceptions
 * or chains of status / exceptions. Maybe it is a stack. Maybe all
 * API objects should have a Status member that can accumulate Status
 * objects of child objects/operations.
 */
class Status
{
    public:
        /*!
         * \brief Default constructor.
         */
        Status();
        /*!
         * \brief Copy constructor
         * \param status
         */
        Status(const Status &status);
        /*!
         * \brief Move constructor.
         * \param status
         */
        Status(Status &&status) noexcept;
        /*!
         * \brief Copy assignment operator.
         * \param status
         * \return reference to lhs.
         */
        Status &operator=(const Status &status);
        /*!
         * \brief Move assignment operator.
         * \param status
         * \return reference to lhs
         */
        Status &operator=(Status &&status) noexcept;
        /*!
         * \brief Converting assignment operator.
         * \param success
         * \return reference to lhs
         */
        Status &operator=(const bool &success);

        /*!
         * \brief Converting constructor.
         * \param success
         */
        explicit Status(const bool &success);

        /*!
         * \brief non-virtual destructor
         *
         * Do not inherit from this class.
         */
        ~Status();
        /*
         * \brief Check success status.
         *
         * \return true if the operation described was successful.
         */
        bool success() const;
    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
};

}

#endif // header guard
