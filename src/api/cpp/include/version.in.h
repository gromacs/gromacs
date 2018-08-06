/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
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
 * Goals: do better than compiler error, build time linker errors, and run time linker errors.
 *
 *  Todo: provide versioning for headers and library so clients can do abstract comparison of build versus runtime.
 *
 *  Todo: Provide better boilerplate (at least) for client self-test of API version compatibility at build time.
 *
 *  Todo: Add compatibility test/warning at module load when client was compiled against a different libgmxapi.
 *  \ingroup gmxapi
 */

// The version numbers for gmxapi are encoded in the repository solely in the `src/api/CMakeLists.txt` file.
// During CMake configuration, the `src/api/cpp/gmxapi/version.h` file is created so that the built library can
// report this version through the `gmxapi::Version` interface. Client code should include the installed
// `gmxapi/version.h` header in order to embed the constants `gmxapi::GMXAPI_MAJOR`, `gmxapi::GMXAPI_MINOR`,
// and `gmxapi::GMXAPI_PATCH` so that API compatibility checks can be performed at runtime.
//
// When a new software release is tagged, the next commit on the development branch should increment the patch level
// to distinguish development builds from the tagged release. As incompatibilities are introduced
// in feature branches, minor or major version number should be incremented as appropriate. At this time,
// client code has no indication of whether the version presented in a development build of gmxapi is an
// officially specified API revision or is subject to change. Developers coding against development branches
// should keep this in mind. If this becomes problematic, please offer your suggestions or propose a revision
// to the `gmxapi::Version` API.

#include <string>

namespace gmxapi
{

using version_t = unsigned short int;

// Todo: it may be preferable for CMake to get the version from the header instead of the other way around.
// It would be nice to be able to pull the headers straight from the repository...
static constexpr const version_t GMXAPI_MAJOR     = @GMXAPI_MAJOR@;
static constexpr const version_t GMXAPI_MINOR     = @GMXAPI_MINOR@;
static constexpr const version_t GMXAPI_PATCH     = @GMXAPI_PATCH@;
static const char                GMXAPI_RELEASE[] = "@GMXAPI_RELEASE@";

/*!
 * \brief Provide API library version information for client code.
 *
 * Allow client code to query the currently loaded gmxapi library object to find the built version. Provide helpers
 * to compare against the features for which the client was written and the headers against which it was compiled.
 *
 * \ingroup gmxapi
 */
class Version
{
    public:
        /// Query gmxapi major version.
        /// \returns major version number
        static version_t major();
        /// Query gmxapi minor version.
        /// \returns minor version number
        static version_t minor();
        /// Query gmxapi patch level.
        /// \returns patch level number
        static version_t patch();
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
        static bool has_feature(const std::string &featurename);
        /// Check for sufficiently high API version number.
        /// \returns `true` if gmxapi library version is the same or greater than the argument(s).
        /// \param major gmxapi major version number.
        /// \param minor gmxapi minor version number (optional).
        /// \param patch patch level of the api library (optional).
        static bool is_at_least(version_t major, version_t minor = 0, version_t patch = 0);
};

}      // namespace gmxapi

#endif // version.h include guard
