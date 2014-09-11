/*
 * This source file is part of the Alexandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef XML_UTIL_H
#define XML_UTIL_H

#include <stdio.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "gromacs/utility/fatalerror.h"

/* This source code file is part of the Alexandria project */

#ifdef __cplusplus
extern "C" {
#endif

extern int find_elem(char *name, int nr, const char *names[]);

extern void add_xml_int(xmlNodePtr ptr, const char *name, int val);

extern void add_xml_double(xmlNodePtr ptr, const char *name, double val);

extern void add_xml_char(xmlNodePtr ptr, const char *name, const char *val);

extern xmlNodePtr add_xml_child(xmlNodePtr parent, const char *type);

extern xmlNodePtr add_xml_child_val(xmlNodePtr parent, const char *type, const char *value);

extern xmlNodePtr add_xml_comment(xmlDocPtr doc,
                                  xmlNodePtr prev, char *comment);

#ifdef __cplusplus
}
#endif

#endif
