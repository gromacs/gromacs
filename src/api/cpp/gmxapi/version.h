#ifndef GMXAPI_VERSION_H
#define GMXAPI_VERSION_H
/*! \file
 * \brief Implement versioning API for C++ external Gromacs interface.
 *  Versioning follows semantic versioning scheme in which the major version
 *  specifies the API compatibility level, and minor version indicates additional
 *  features that may or may not be ABI compatible.
 *
 *  Defines a class Version in the gmxapi namespace providing static methods
 *  for use by API client code at compile time and run time.
 *
 *  Todo: provide versioning for headers and library so clients can do abstract comparison of build versus runtime.
 *  \ingroup gmxapi
 */

#include <string>

namespace gmxapi
{

/// Provide API version information for client code.
/// \ingroup gmxapi
class Version
{
    public:
        /// Query gmxapi major version.
        /// \returns major version number
        static unsigned int major();
        /// Query gmxapi minor version.
        /// \returns minor version number
        static unsigned int minor();
        /// Query gmxapi patch level.
        /// \returns patch level number
        static unsigned int patch();
        /// Get formatted release string.
        /// Format is major.minor.patch
        /// \returns release as string
        static std::string release();
        /// Check features availability
        /// \returns `true` if the named feature is available.
        /// Features introduced after release 1.0.0 may be named in the documentation
        /// to improve readability of client code and simplify development. Prefer
        /// this mechanism when checking for features still under development or to
        /// distinguish between interface levels of a specific feature.
        /// \param featurename Feature name described in the feature's documentation.
        static bool has_feature(std::string featurename);
        /// Check for sufficiently high API version number.
        /// \returns `true` if gmxapi version is the same or greater than the argument(s).
        /// \param major gmxapi major version number.
        /// \param minor gmxapi minor version number (optional).
        /// \param patch patch level of the api library (optional).
        static bool is_at_least(unsigned int major, unsigned int minor = 0, unsigned int patch = 0);
};

}      // namespace gmxapi

#endif // version.h include guard
