/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
Classes and other symbols that are publicly accessible within the Gromacs
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

/*!
\namespace gmx
\brief
Generic Gromacs namespace.

\inpublicapi
*/

/*!
\internal \namespace gmx::internal
\brief
Internal Gromacs namespace.

This namespace is used to contain some implementation-specific functions and
classes.  These are not meant for direct user access, but typically reside
in public headers because of implementation reasons.
*/

/*!
\dir src
\brief
Main source code directory.
 */

/*!
\dir src/gromacs
\brief
Source code directory for building the `libgromacs` library.
 */

/*!
\libinternal
\dir src/programs
\brief
Source code directory for building executables.
 */

/*!
\dir share
\brief
Directory that contains installed data files.
 */

/*!
\dir share/template
\brief
Template code for writing analysis programs.
 */

/*!
\file share/template/template.cpp
\brief
Template code for writing analysis programs.

See \ref page_analysistemplate for more information.
 */

/*!
\example template.cpp
\brief
Template code for writing analysis programs.

See \ref page_analysistemplate for more information.
 */
