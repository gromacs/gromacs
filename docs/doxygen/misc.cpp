/*!
\internal \file
\brief
Doxygen documentation file for group declarations etc.

\author Teemu Murtola <teemu.murtola@gmail.com>
*/

/*!
\defgroup group_publicapi Public API
\brief
Classes and other symbols that are publicly accessible from user code.
*/

/*!
\defgroup group_libraryapi Library API
\brief
Classes and other symbols that are publicly accessible within the \Gromacs
library.

\see group_publicapi
*/

/*!
\defgroup group_utilitymodules Utility Modules
\brief
Modules with generic utility functions.
*/

/*!
\defgroup group_analysismodules Analysis Modules
\brief
Modules used in analysis tools.

A separate page describes the responsibilities of these modules:
\ref page_analysisframework
*/

/*! \libinternal
\defgroup group_mdrun Modules for simulation functionality
\brief
Modules used in running simulations with mdrun
*/

/*!
\namespace gmx
\brief
Generic \Gromacs namespace.

\inpublicapi
*/

/*!
\internal \namespace gmx::internal
\brief
Internal \Gromacs namespace.

This namespace is used to contain some implementation-specific functions and
classes.  These are not meant for direct user access, but typically reside
in public headers because of implementation reasons.
*/

/*!
\libinternal
\namespace gmx::compat
\brief Compatibility aliases for standard library features.

Provide consistent naming for standard library features that must be back-ported
on some platforms.
gmx::compat::some_feature may map to back-ported code or to a feature provided by
the STL available on a given build platform, but by including the compatibility
header and using the gmx::compat namespace, forward and backward compatible code
is cleaner and clearer. In the future, when a feature is determined to be provided
by the system on all supported platforms, references to gmx::compat::some_feature
can be replaced, e.g. with std::some_feature, and gmx::compat::some_feature
deprecated.

Since compatibility headers are likely to be included by other
headers, there is a risk of ambiguity if code in the gmx namespace refers to an
unqualified name in the std namespace. To reduce ambiguity, symbol names from
gmx::compat should not be imported into scopes that are shared between multiple
translation units (e.g. via `using` statements in header files).

\ingroup group_compatibility
*/

/*!
\libinternal
\defgroup group_compatibility C++ standard library compatibility helpers.
\brief Provide uniform interface to selected C++ standard library features.

For some features not available on all platforms supported by \Gromacs,
provide back-ports or mappings to available standard library implementations
as appropriate.
*/

/*!
\file share/template/template.cpp
\brief Template code for writing analysis programs.

See \ref page_analysistemplate for more information.
 */

/*!
\example template.cpp
\brief Template code for writing analysis programs.

See \ref page_analysistemplate for more information.
 */
