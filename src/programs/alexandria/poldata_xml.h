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
#ifndef POLDATA_XML_H
#define POLDATA_XML_H

#include "gromacs/topology/atomprop.h"

#include "poldata.h"
#include "xml_util.h"


//extern int xmlDoValidityCheckingDefaultValue;
namespace alexandria
{
//using namespace alexandria;
/* This source code file is part of the Alexandria project */
class PoldataXml
{
    public:

        PoldataXml(){}

        static void write(const std::string fn, Poldata * pd,
                          gmx_bool bCompress);

        static Poldata * read(const char *fn, gmx_atomprop_t aps);

    private:
        static void sp(int n, char buf[], int maxindent);

        static void processAttr(FILE *fp, xmlAttrPtr attr, int elem,
                                int indent, Poldata  * pd);

        static void processTree(FILE *fp, xmlNodePtr tree, int indent,
                                Poldata * pd, gmx_atomprop_t aps);

        static void addXmlPoldata(xmlNodePtr parent, Poldata * pd);


};
}
#endif
