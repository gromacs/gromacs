/*
 * $Id: xml_util.h,v 1.3 2009/03/15 21:22:10 spoel Exp $
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */

#ifndef _xml_util_h
#define _xml_util_h


#include <stdio.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "gromacs/utility/fatalerror.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int find_elem(char *name, int nr, const char *names[]);

extern void add_xml_int(xmlNodePtr ptr, const char *name, int val);

extern void add_xml_double(xmlNodePtr ptr, const char *name, double val);

extern void add_xml_char(xmlNodePtr ptr, const char *name, const char *val);

extern xmlNodePtr add_xml_child(xmlNodePtr parent, const char *type);

extern xmlNodePtr add_xml_child_char(xmlNodePtr parent, const char *type, const char *val);

extern xmlNodePtr add_xml_child_double(xmlNodePtr parent, const char *type, double val);

extern xmlNodePtr add_xml_child_int(xmlNodePtr parent, const char *type, int val);

extern xmlNodePtr add_xml_comment(xmlDocPtr doc,
                                  xmlNodePtr prev, char *comment);

#ifdef __cplusplus
}
#endif

#endif
