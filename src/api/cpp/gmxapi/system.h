#ifndef GMXAPI_SYSTEM_H
#define GMXAPI_SYSTEM_H
/*! \file
 * \brief Declare container for molecular systems
 *
 * \ingroup gmxapi
 */
#include <memory>

#include "gmxapi/gmxapi.h"

namespace gmxapi
{

class Atoms;
class MDEngine;
class IMDRunner;

/// Container for molecular model and simulation parameters.
/*!
 * \cond
 * A system instance is sort of a container of builders, and a Context is sort of a factory. together they allow a simulation to be constructed and initialized with the appropriate implementations of runner, integrator, and data objects.
 * \endcond
 * \ingroup gmxapi
 */
class System
{
    public:
        /// A blank system object is possible, but not yet useful.
        System();
        /// No copy.
        /*! The semantics of copying a System are ambiguous, so disallow implicit
         * copy. Some sort of prototype or clone idiom is probably useful, but
         * needs to explicitly identify any expensive operations.
         */
        System(const System &)            = delete;
        /// No copy.
        System &operator=(const System &) = delete;

        /// Allow move.
        System(System &&) noexcept;
        /// Allow move.
        System &operator=(System &&) noexcept;

        /// \cond internal
        /// Destructor defined later to allow unique_ptr members of partially-defined types.
        ~System();
        /// \endcond

        /// The mechanics of building a System object are deferred to the builder(s).
        /// The Builder class(es) do not yet have a public interface.
        class Builder;

//        /// Get a handle to system atoms.
//        std::unique_ptr<Atoms> atoms();

        /// Get a handle to bound MD engine.
        std::shared_ptr<MDEngine> md();

        /// Set the MD engine
        void md(std::shared_ptr<MDEngine> md);

        /// Get a handle to bound runner.
        std::shared_ptr<IMDRunner> runner();

        /// Set the runner.
        void runner(std::shared_ptr<IMDRunner> runner);

//        /// Invoke an appropriate runner, if possible.
//        /// Equivalent to system->runner()->initialize(defaultContext())->run();
//        Status run();

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
};


/// Process a TPR file into a ready-to-run system.
/*!  Manages the mdrun module, options processing, input record, and initialization
 * to produce a ready-to-run MD simulation.
 *
 * Reads the input record to identify the appropriate MDEngine and
 * extract the appropriate System members. Returns System with bound
 * runner, mdengine, and atoms.
 *
 *
 * \param filename Filesystem path of TPR file.
 * \returns gmxapi::System with bound objects and parameters specified in TPR.
 * \ingroup gmxapi
 */
std::unique_ptr<gmxapi::System> fromTprFile(std::string filename);

}      // end namespace gmxapi

#endif // include guard
