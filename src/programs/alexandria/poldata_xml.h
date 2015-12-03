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
