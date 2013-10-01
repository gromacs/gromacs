/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Declares utilities for reading xml files using libxml2.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_waxsdebye
 */
#ifndef _xml_util_h
#define _xml_util_h

#include <libxml/parser.h>
#include <libxml/tree.h>

/*! \brief
 * Look for string in list of strings
 *
 * \param[in] name  The string to look for
 * \param[in] nr    The length of the list of strings
 * \param[in] names The list of strings
 * \return The index in the list or -1 if not found
 */
int find_elem(char *name, int nr, const char *names[]);

/*! \brief
 * Add integer attribute in xml tree
 *
 * \param[in] ptr   The xml tree
 * \param[in] name  The name of the attribute
 * \param[in] val   The value of the attribute
 */
void add_xml_int(xmlNodePtr ptr, const char *name, int val);

/*! \brief
 * Add double attribute in xml tree
 *
 * \param[in] ptr   The xml tree
 * \param[in] name  The name of the attribute
 * \param[in] val   The value of the attribute
 */
void add_xml_double(xmlNodePtr ptr, const char *name, double val);

/*! \brief
 * Add string attribute in xml tree
 *
 * \param[in] ptr   The xml tree
 * \param[in] name  The name of the attribute
 * \param[in] val   The value of the attribute
 */
void add_xml_char(xmlNodePtr ptr, const char *name, const char *val);

/*! \brief
 * Add branch to xml tree
 *
 * \param[in] parent The xml tree
 * \param[in] type   The name of the branch
 * \return The pointer to the branch
 */
xmlNodePtr add_xml_child(xmlNodePtr parent, const char *type);

/*! \brief
 * Add branch with character constant to xml tree
 *
 * \param[in] parent The xml tree
 * \param[in] type   The name of the branch
 * \param[in] val    The content of the character constant
 * \return The pointer to the branch
 */
xmlNodePtr add_xml_child_char(xmlNodePtr parent, const char *type, const char *val);

/*! \brief
 * Add branch with double constant to xml tree
 *
 * \param[in] parent The xml tree
 * \param[in] type   The name of the branch
 * \param[in] val    The double constant
 * \return The pointer to the branch
 */
xmlNodePtr add_xml_child_double(xmlNodePtr parent, const char *type, double val);

/*! \brief
 * Add branch with integer constant to xml tree
 *
 * \param[in] parent The xml tree
 * \param[in] type   The name of the branch
 * \param[in] val    The integer constant
 * \return The pointer to the branch
 */
xmlNodePtr add_xml_child_int(xmlNodePtr parent, const char *type, int val);

/*! \brief
 * Add top comment to xml tree
 *
 * \param[in] doc     The xml document
 * \param[in] prev    The xml tree
 * \param[in] comment The comment to add
 * \return The pointer to the comment branch
 */
xmlNodePtr add_xml_comment(xmlDocPtr doc,
                           xmlNodePtr prev, char *comment);

#endif
