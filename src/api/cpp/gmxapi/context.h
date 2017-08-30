#ifndef GMXAPI_CONTEXT_H
#define GMXAPI_CONTEXT_H
/*! \file
 * \brief Declares classes representing execution context.
 *
 * \ingroup gmxapi
 */


#include <memory>

namespace gmxapi
{

class Status;
class IRunner;
class MDInput;

/// Execution context.
/*!
 * A proxy can be configured with information needed to initialize a runtime
 * environment capable of executing a work load, independently of defining the
 * work.
 * The actual execution
 * environment is not necessarily instantiated / configured until the work is
 * performed.
 * Thus, construction of a Context object does not necessarily imply
 * initialization of compute resources, but any active compute resources are
 * appropriately deinitialized when the object is destroyed. However, to allow
 * opportunities for exception handling, resources should be deinitialized when
 * and as work is completed by the Runner.
 *
 * Ultimately, the class in this header file should just define an interface.
 * Implementations for different execution environments will be provided by the
 * library or extension code and documented separately,
 * and it should be straight-forward to
 * write extension code for other execution environments.
 *
 * In the first draft, it is difficult to separate object definition from
 * initialization.
 * \ingroup gmxapi
 */
class Context
{
    public:
        Context();
        ~Context();

        // Disallow copy
        Context(const Context &)            = delete;
        Context &operator=(const Context &) = delete;

        // Allow move
        Context(Context &&)            = default;
        Context &operator=(Context &&) = default;

        /*! \brief Initialize execution context.
         *
         * After the call, the execution context will either be fully configured
         * and running, or the status will describe why the work cannot be executed.
         */
//        Status initialize();
        /*! \brief Deinitialize execution context.
         *
         * The context should be deinitialized after work completes. Status will
         * describe errors and exceptions that occurred during attempts to shutdown
         * and free resources.
         */
//        Status deinitialize();

        /// Returns true while context is initialized / executing.
//        bool isInitialized() const;

        /// Bind to a runner.
//        Status setRunner(std::shared_ptr<IRunner> runner);

        /// Get a handle to the Runner.
//        std::shared_ptr<IRunner> runner();

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
};

std::unique_ptr<Context> defaultContext();

}      // end namespace gmxapi

#endif // header guard
