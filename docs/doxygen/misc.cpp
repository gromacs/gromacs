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
\file share/template/template.cpp
\brief Template code for writing analysis programs.

See \ref page_analysistemplate for more information.
 */

/*!
\example template.cpp
\brief Template code for writing analysis programs.

See \ref page_analysistemplate for more information.
 */
