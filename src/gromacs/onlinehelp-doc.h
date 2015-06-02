/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
/*! \libinternal \defgroup module_onlinehelp Help Formatting for Online Help (onlinehelp)
 * \ingroup group_utilitymodules
 * \brief
 * Provides functionality for formatting help text for console and reStructuredText.
 *
 * This module provides helper functions and classes for formatting text to the
 * console and as reStructuredText through a single interface.  The main
 * components of the module are:
 *  - gmx::HelpWriterContext provides a single interface that can produce both
 *    output formats from the same input strings and API calls.  Whenever
 *    possible, the output format should be abstracted using this interface,
 *    but in some cases code still writes out raw reStructuredText.
 *  - rstparser.h provides the functionality to parse reStructuredText such that
 *    it can be rewrapped for console output.
 *  - helpformat.h provides some general text-processing classes, currently
 *    focused on producing aligned tables for console output.
 *  - helptopicinterface.h, helptopic.h, and helpmanager.h provide classes for
 *    managing a hierarchy of help topics and printing out help from this
 *    hierarchy.
 *
 * The formatting syntax for strings accepted by this module is described in
 * \ref page_onlinehelp.  The module is currently exposed outside \Gromacs only
 * through this formatting syntax, not any API calls.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 */
/*! \internal \file
 * \brief
 * Dummy header for \ref module_onlinehelp documentation.
 */
