/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
#include <string.h>
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "xml_util.h"

int find_elem(char *name, int nr, const char *names[])
{
    int i;

    for (i = 0; (i < nr); i++)
    {
        if (gmx_strcasecmp(name, names[i]) == 0)
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

    if ((child = xmlNewChild(parent, /*parent->ns*/ NULL, (xmlChar *)type, NULL)) == NULL)
    {
        gmx_fatal(FARGS, "Creating element", (char *)type);
    }

    return child;
}

xmlNodePtr add_xml_child_char(xmlNodePtr parent, const char *type, const char *val)
{
    xmlNodePtr child;

    if ((child = xmlNewChild(parent, NULL, (xmlChar *)type, (xmlChar *)val)) == NULL)
    {
        gmx_fatal(FARGS, "Creating element", (char *)type);
    }

    return child;
}

xmlNodePtr add_xml_child_double(xmlNodePtr parent, const char *type, double val)
{
    char buf[16];

    sprintf(buf, "%12g", val);
    return add_xml_child_char(parent, type, buf);
}

xmlNodePtr add_xml_child_int(xmlNodePtr parent, const char *type, int val)
{
    char buf[16];

    sprintf(buf, "%12d", val);
    return add_xml_child_char(parent, type, buf);
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
