#include "gmxapi/gmxapi.h"
#include "gromacs/compat/make_unique.h"

namespace gmxapi
{

/*! \cond internal
 * \brief Implementation class for Status objects.
 */
class Status::Impl
{
    public:
        /*!
         * \brief Default construct as unsuccessful status.
         */
        Impl() : success_ {false}
        {};
        /*!
         * \brief Construct with success for true input.
         * \param success let Boolean true == success.
         */
        explicit Impl(const bool &success) :
            success_ {success}
        {};
        /*!
         * \brief Query success status
         * \return true if successful
         */
        bool success() const
        {
            return success_;
        };
    private:
        bool success_;
};
/// \endcond

Status::Status() :
    impl_ {gmx::compat::make_unique<Status::Impl>()}
{}

Status::Status(const Status &status)
{
    impl_.reset(new Impl {status.success()});
}

Status &Status::operator=(const Status &status)
{
    this->impl_.reset(new Impl {status.success()});
    return *this;
}

Status &Status::operator=(Status &&status) noexcept
{
    this->impl_ = std::move(status.impl_);
    return *this;
}

Status &Status::operator=(const bool &success)
{
    this->impl_.reset(new Impl {success});
    return *this;
}

Status::Status(Status &&status) noexcept
{
    this->impl_ = std::move(status.impl_);
}

Status::Status(const bool &success) :
    impl_ {gmx::compat::make_unique<Status::Impl>(success)}
{}

bool Status::success() const
{
    return impl_->success();
}

// Destructor must be defined after Impl to use unique_ptr<Impl>
Status::~Status() = default;

} // end namespace gmxapi
