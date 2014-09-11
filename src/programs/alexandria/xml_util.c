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
#include <string.h>
#include "gromacs/utility/fatalerror.h"
#include "xml_util.h"

int find_elem(char *name, int nr, const char *names[])
{
    int i;

    for (i = 0; (i < nr); i++)
    {
        if (strcmp(name, names[i]) == 0)
        {
            return i;
        }
    }

    return -1;
}

void add_xml_int(xmlNodePtr ptr, const char *name, int val)
{
    xmlChar buf[32];

    sprintf((char *)buf, "%d", val);
    if (xmlSetProp(ptr, (xmlChar *)name, buf) == 0)
    {
        gmx_fatal(FARGS, "Setting", (char *)name);
    }
}

void add_xml_double(xmlNodePtr ptr, const char *name, double val)
{
    xmlChar buf[32];

    sprintf((char *)buf, "%g", val);
    if (xmlSetProp(ptr, (xmlChar *)name, buf) == 0)
    {
        gmx_fatal(FARGS, "Setting", (char *)name);
    }
}

void add_xml_char(xmlNodePtr ptr, const char *name, const char *val)
{
    if (xmlSetProp(ptr, (xmlChar *)name, (xmlChar *)val) == 0)
    {
        gmx_fatal(FARGS, "Setting", (char *)name);
    }
}

xmlNodePtr add_xml_child(xmlNodePtr parent, const char *type)
{
    xmlNodePtr child;

    if ((child = xmlNewChild(parent, NULL, (xmlChar *)type, NULL)) == NULL)
    {
        gmx_fatal(FARGS, "Creating element", (char *)type);
    }

    return child;
}

xmlNodePtr add_xml_child_val(xmlNodePtr parent, const char *type, const char *value)
{
    xmlNodePtr child;

    if ((child = xmlNewChild(parent, NULL, (xmlChar *)type, (xmlChar *)value)) == NULL)
    {
        gmx_fatal(FARGS, "Creating element", (char *)type);
    }

    return child;
}

xmlNodePtr add_xml_comment(xmlDocPtr doc,
                           xmlNodePtr prev, char *comment)
{
    xmlNodePtr comm, ptr;

    if ((comm = xmlNewComment((xmlChar *)comment)) == NULL)
    {
        gmx_fatal(FARGS, "Creating doc comment element", "");
    }
    ptr = prev;
    while (ptr->next != NULL)
    {
        ptr = ptr->next;
    }
    ptr->next    = comm;
    comm->prev   = ptr;
    comm->doc    = doc;

    return comm;
}
