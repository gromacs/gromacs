#ifndef GMXAPI_EXCEPTIONS_H
#define GMXAPI_EXCEPTIONS_H
/*! \file
 * \brief Declare exception classes for external API.
 *
 * \ingroup gmxapi
 */

#include <exception>
#include <string>

namespace gmxapi
{

/*! \brief Base exception for gmxapi library.
 *
 * Exceptions thrown in the gmxapi namespace are descended from gmxapi::Exception
 * or there is a bug.
 *
 * \ingroup gmxapi
 */
class Exception : public std::exception
{
    public:
        Exception()                             = default;
        ~Exception() override                   = default;
        Exception(const Exception &)            = default;
        Exception &operator=(const Exception &) = default;

        Exception(Exception &&) noexcept            = default;
        Exception &operator=(Exception &&) noexcept = default;

        const char* what() const noexcept override
        {
            return "Gromacs API error";
        };
};

template<class E>
class BasicException : public Exception
{
    private:
        std::string what_;
    public:
        BasicException() : BasicException{std::string()}
        {};

        explicit BasicException(std::string &&message) noexcept :
            what_ {std::move(message)}
        {};

        explicit BasicException(const char* message)
        {
            what_ = message;
        }

        const char* what() const noexcept override
        {
            return what_.c_str();
        }
};

/*! \brief Behavioral protocol violated.
 *
 * Indicates that a behavioral protocol specified in the API is not being followed. The class
 * throwing this exception expects certain methods to be called in a certain order.
 *
 * If this exception is encountered in client code, the API is being misused or there is a bug.
 * Generally, required behaviors should be implemented in templates or base classes rather than
 * exposing and requiring complete implementation of the protocol in client code.
 */
class ProtocolError : public BasicException<ProtocolError>
{
    public:
        using BasicException<ProtocolError>::BasicException;
};

class NotImplementedError : public BasicException<NotImplementedError>
{
    public:
        using BasicException<NotImplementedError>::BasicException;
};

}      // end namespace gmxapi

#endif // header guard
