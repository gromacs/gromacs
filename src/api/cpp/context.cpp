/*! \file
 * \brief Implementation details of gmxapi::Context
 *
 */

#include <memory>

#include "gmxapi/context.h"
#include "gmxapi/gmxapi.h"
#include "gmxapi/runner.h"

#include "gromacs/compat/make_unique.h"

namespace gmxapi
{

/*! \brief Temporary shim until proper messaging.
 *
 */
class warn
{
    public:
        /*! \brief Create a warning message.
         *
         * \param message must outlive warn instance...
         */
        warn(const char* message) :
            message_ {message}
        { };
        /// pointer to string managed somewhere else.
        const char* message_;
};

class Context::Impl
{
    public:
        Impl() = default;
//        bool                       initialized_;
//        std::shared_ptr<IRunner>    runner_;
//        std::unique_ptr<t_commrec> commRec_;
};

//Context::Impl::Impl() :
//    initialized_ {false},
//runner_ {
//    nullptr
//},
//commRec_ {
//    nullptr }
//{}

Context::Context() :
    impl_ {gmx::compat::make_unique<Context::Impl>()}
{}

Context::~Context() = default;

//Status Context::initialize()
//{
//    impl_->initialized_ = true;
//    impl_->commRec_.reset(init_commrec());
// Multi-sim not yet supported.
// i.e. mdrun.cpp: init_multisystem(...)
// Restart not yet supported
// mdrun.cpp:     handleRestart(...);
// Log file not yet supported

// GPU not supported yet
// runner.cpp: mdrunner()->gmx_detect_hardware()

//MPI not yet supported.

//    return Status(true);
//}

//Status Context::deinitialize()
//{
//    if (impl_->commRec_ != nullptr)
//    {
//        // TODO: This is NOT a proper MPI shutdown! It just mimics gmx_mdrun()...
//        if (impl_->commRec_->mpb != nullptr)
//        {
//            delete impl_->commRec_->mpb;
//        }
//        impl_->commRec_.reset();
//    }
//    impl_->initialized_ = false;
//    return Status(true);
//}

//bool Context::isInitialized() const
//{
//    return impl_->initialized_;
//}

//Status Context::setRunner(std::shared_ptr<IRunner> runner)
//{
//    return Status(runner != nullptr);
//}

//std::shared_ptr<IRunner> Context::runner()
//{
//    return impl_->runner_;
//}

//Context::~Context()
//{
//    if (isInitialized())
//    {
//        warn("Context was not deinitialized before destruction.\n");
//        try
//        {
//            deinitialize();
//        }
//        catch (const std::exception &)
//        {
//            warn("Exceptions occurred during Context destruction.");
//        }
//    }
//}

// std::unique_ptr<Context> Context::fromInputRec()
// {
//
//     context->impl_->initialized_ = true;
//
//     context->impl_->commRec_.reset(init_commrec());
//
//     return gmx::compat::make_unique<Context>(std::move(context));
// }

std::unique_ptr<Context> defaultContext()
{
    auto context = gmx::compat::make_unique<Context>();
    return context;
}

} // end namespace gmxapi
