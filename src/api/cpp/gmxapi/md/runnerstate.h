//
// Created by Eric Irrgang on 8/16/17.
//

#ifndef GROMACS_RUNNERPROXY_H
#define GROMACS_RUNNERPROXY_H

#include "gmxapi/runner.h"
#include <memory>

// Declaring classes from other namespaces is iffy and turns compile errors into linking errors
// but reduces coupling. If we go this route, TODO: consider consolidating incomplete gmx types in a single header.
//namespace gmx
//{
//    class TpxState;
//}
// For now, I prefer the compiler type-checking of the coupled code.
#include "gromacs/mdtypes/TpxState.h"

namespace gmxapi
{

class EmptyMDRunnerState : public IMDRunner
{
    public:
        EmptyMDRunnerState() = default;

        Status run() override;

        std::shared_ptr<IMDRunner> initialize(std::shared_ptr<Context> context) override;

        void registerMDBuilder(std::unique_ptr<MDBuilder> builder) override;

//        std::shared_ptr<IMDRunnerBuilder> builder() override;
};

/*!
 * \brief An MDRunner that has not yet started.
 *
 * Accumulates configuration information that can be used to launch a gmx::Mdrunner.
 */
class UninitializedMDRunnerState : public IMDRunner
{
    public:
        ~UninitializedMDRunnerState() override;

        // Disallow copy
        UninitializedMDRunnerState(const UninitializedMDRunnerState &)            = delete;
        UninitializedMDRunnerState &operator=(const UninitializedMDRunnerState &) = delete;

        // Allow move
        UninitializedMDRunnerState(UninitializedMDRunnerState &&) noexcept            = default;
        UninitializedMDRunnerState &operator=(UninitializedMDRunnerState &&) noexcept = default;

        Status run() override;

        std::shared_ptr<IMDRunner> initialize(std::shared_ptr<Context> context) override;

        void registerMDBuilder(std::unique_ptr<MDBuilder> builder) override;

        class Builder
        {
            public:
                Builder();
                ~Builder();
                Builder(const Builder &)                = delete;
                Builder(Builder &&) noexcept            = default;
                Builder &operator=(const Builder &)     = delete;
                Builder &operator=(Builder &&) noexcept = default;

                Builder &mdEngine(std::shared_ptr<MDEngine> md);
                Builder &tpxState(std::shared_ptr<gmx::TpxState> input);

                std::unique_ptr<UninitializedMDRunnerState> build();
            private:
                std::unique_ptr<UninitializedMDRunnerState> runner_;
        };

    private:
        UninitializedMDRunnerState();
        /// Private implementation class
        class Impl;
        /// pointer to implementation
        std::unique_ptr<Impl> impl_;
};

/*!
 * \brief Handle to an active gmx::Mdrunner
 */
class RunningMDRunnerState : public IMDRunner
{
    public:
        ~RunningMDRunnerState() override;

        Status run() override;

        std::shared_ptr<IMDRunner> initialize(std::shared_ptr<Context> context) override;

        void registerMDBuilder(std::unique_ptr<MDBuilder> builder) override;

        // Todo: we can just template some of this: class Builder : public RunnerBuilder<RunningMDRunnerState>
        class Builder
        {
            public:
                Builder();
                ~Builder();
                Builder(const Builder &)                = delete;
                Builder(Builder &&) noexcept            = default;
                Builder &operator=(const Builder &)     = delete;
                Builder &operator=(Builder &&) noexcept = default;

                Builder &tpxState(std::shared_ptr<gmx::TpxState> input);

                std::unique_ptr<RunningMDRunnerState> build();
            private:
                std::unique_ptr<RunningMDRunnerState> runner_;
                std::shared_ptr<gmx::TpxState>        tpxState_;
        };
    private:
        RunningMDRunnerState();
        /// Private implementation class
        class Impl;
        /// pointer to implementation
        std::unique_ptr<Impl> impl_;
};


}      // end namespace gmxapi
#endif //GROMACS_RUNNERPROXY_H
