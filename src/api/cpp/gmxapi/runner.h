#ifndef GMXAPI_RUNNER_H
#define GMXAPI_RUNNER_H

/*! \file
 * \brief Public C++ API for Gromacs module runners.
 *
 * \ingroup gmxapi
 */

#include <memory>
#include <string>

#include "gmxapi/gmxapi.h"
#include "gmxapi/system.h"

namespace gmxapi
{
class MDProxy;
class Context;

class MDBuilder;

/*!
 * \brief Interface provided by runners for MD tasks.
 *
 * A runner implements this interface in order to bind to an MD task. A caller passes an IMDRunner pointer to the
 * bind() method of an MDProxy object. The MDProxy object then provides the runner with a builder for an MD task by
 * calling IMDRunner::registerMDBuilder(std::unique_ptr<MDBuilder>).
 *
 * The caller of bind() should guarantee the lifetime of IMDRunner through the subsequent call to registerMDBuilder().
 * The call to bind() should be made when the caller is ready to construct an MD task, such that the state of the
 * MDProxy at the time of the call is appropriate for the MD object to be built. The likely use case is for the call
 * to be made during the builder of a runner that is about to execute.
 *
 * Todo: Check or enforce assumptions: the base class for the builder can guarantee that registerMDBuilder()
 * will be called before destruction.
 *
 * Example:
 * \code
 * class MyRunner : public IMDRunner
 * {
 * private:
 *      std::shared_ptr<MDProxy> md_;
 *      std::unique_ptr<ModuleBuilder> mdbuilder_;
 * public:
 *      IMDRunner* actualRunner() {return this;};
 *      virtual void registerMDBuilder(std::unique_ptr<ModuleBuilder> builder) override
 *      {
 *          mdbuilder_ = std::move(builder);
 *      };
 *
 *      void run()
 *      {
 *          md_->bind(actualRunner);
 *          auto actualMd = mdbuilder_.build();
 *          gmx::MDRunner runner{};
 *          runner.setMD(actualMD);
 *          runner.run();
 *      };
 * };
 *
 * class MyMDProxy : public MDState
 * {
 * public:
 *      std::unique_ptr<ModuleBuilder> builder();
 *      virtual void bind(IMDRunner* runner)
 *      {
 *          runner->registerMDBuilder(builder());
 *      };
 * };
 *
 * \endcode
 */
class IMDRunner
{
    public:
        virtual ~IMDRunner() = default;

        virtual void registerMDBuilder(std::unique_ptr<MDBuilder> builder) = 0;
        virtual Status run() = 0;
        // This is an overly complicated and awkward pattern. Probably replace with
        // templated implementation in base class using curiously-recurring-template idiom
        //virtual Status run(const EndCondition& condition) = 0;

        /*!
         * \brief Activate a GROMACS runner in the given execution context.
         *
         * Launches a fully configured workflow.
         * Changes the state of the runner from uninitialized to initialized.
         * \param context Execution environment and resources.
         * \return Handle to a launched and runnable workflow.
         */
        virtual std::shared_ptr<IMDRunner> initialize(std::shared_ptr<Context> context) = 0;
};

/*!
 * \brief Get a builder for a concrete runner.
 *
 * A class provides this interface to allow a computational object to be created and launched. The
 * object returned provides a build() method from which to get a runnable object. Other interface
 * features TBD.
 *
 * In general, the class providing this interface will bind to a concrete task before returning from
 * build(). \see registerMDBuiler()
 */
class IMDRunnerBuilder
{
    public:
        virtual ~IMDRunnerBuilder();

        /// Build a runner. Return a handle to something that can be run.
        virtual std::shared_ptr<IMDRunner> build() = 0;
};

/// \cond internal
/*! \brief Runner for trivial graphs.
 *
 * The requirements of running only a single Gromacs tool once are substantially
 * different than those of a runner for a chain or data graph. This interface
 * allows classes to offer only a simple implementation.
 *
 * Implementations for object states will depend on execution context and, possibly,
 * on the module to be run.
 * \ingroup gmxapi
 */
class RunnerProxy : public IMDRunner, public std::enable_shared_from_this<RunnerProxy>
{
    public:
        RunnerProxy();

        // Disallow copy to avoid ambiguity of state ownership
        RunnerProxy(const RunnerProxy &)            = delete;
        RunnerProxy &operator=(const RunnerProxy &) = delete;

        // Allow move
        RunnerProxy(RunnerProxy &&)            = default;
        RunnerProxy &operator=(RunnerProxy &&) = default;

        explicit RunnerProxy(std::shared_ptr<MDProxy> md);

        ~RunnerProxy() override = default;

        Status run() override;

        void registerMDBuilder(std::unique_ptr<MDBuilder> builder) override;

        std::shared_ptr <IMDRunner> initialize(std::shared_ptr<Context> context) override;

        void setState(std::shared_ptr<IMDRunner> state);

    private:
        /// bound task, if any
        std::shared_ptr<MDProxy> module_;

        /// Implementation object we are proxying for
        std::shared_ptr<IMDRunner> instanceState_;
};
/// \endcond



}; // namespace gmxapi


#endif // header guard
