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
