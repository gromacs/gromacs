/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.5
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gmx_fatal.h"


#ifdef HAVE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "gmx_xml.h"
	
void add_xml_int(xmlNodePtr ptr,xmlChar *name,int val)
{
  xmlChar buf[32];
  
  sprintf((char *)buf,"%d",val);
  if (xmlSetProp(ptr,name,buf) == 0)
    gmx_fatal(FARGS,"Setting %s",(char *)name);
}

void add_xml_double(xmlNodePtr ptr,xmlChar *name,double val)
{
  xmlChar buf[32];
  
  sprintf((char *)buf,"%g",val);
  if (xmlSetProp(ptr,name,buf) == 0)
    gmx_fatal(FARGS,"Setting %s",(char *)name);
}

void add_xml_char(xmlNodePtr ptr,xmlChar *name,char *val)
{
  if (xmlSetProp(ptr,name,(xmlChar *)val) == 0)
    gmx_fatal(FARGS,"Setting %s",(char *)name);
}

xmlNodePtr add_xml_child(xmlNodePtr parent,char *type)
{
  xmlNodePtr child;
  
  if ((child = xmlNewChild(parent,NULL,(xmlChar *)type,NULL)) == NULL)
    gmx_fatal(FARGS,"Creating element %s",(char *)type);
  
  return child;
}

xmlNodePtr add_xml_comment(xmlDocPtr doc,
			   xmlNodePtr prev,xmlChar *comment)
{
  xmlNodePtr comm,ptr;
  
  if ((comm = xmlNewComment(comment)) == NULL)
    gmx_fatal(FARGS,"Creating doc comment element");
  ptr = prev;
  while (ptr->next != NULL)
    ptr=ptr->next;
  ptr->next    = comm;
  comm->prev   = ptr;
  comm->doc    = doc;
  
  return comm;
}

#endif
