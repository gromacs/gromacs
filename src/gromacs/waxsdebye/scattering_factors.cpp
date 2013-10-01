/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#include <ios>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <string.h>
#include <libxml/tree.h>
//#include "maths.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "xml_util.h"
#include "scattering_factors.h"

extern int xmlDoValidityCheckingDefaultValue;

using namespace std;

namespace gmx
{

//! Enumerate to distinguish the different input types
enum {
    exmlGROMACS,
    exmlSCATTERINGFACTORS,
    exmlSCATTERINGFACTOR,
    exmlDISPLACEDSOLVENT,
    exmlFORCEFIELD,
    exmlREFERENCE,
    exmlSOURCE,
    exmlRESIDUE,
    exmlATOM,
    exmlIDENTIFIER,
    exmlTYPE,
    exmlUNIT,
    exmlQ0, exmlQ1,
    exmlP0, exmlP1, exmlP2, exmlP3,
    exmlA0, exmlA1, exmlA2, exmlA3, exmlA4, exmlA5,
    exmlB0, exmlB1, exmlB2, exmlB3, exmlB4, exmlB5,
    exmlC,
    exmlNR
};

//! Strings with corresponding names
static const char *exml_names[exmlNR] =
{
    "Gromacs",
    "ScatteringFactors",
    "ScatteringFactor",
    "DisplacedSolvent",
    "ForceField",
    "Reference",
    "Source",
    "Residue",
    "Atom",
    "Identifier",
    "Type",
    "Unit",
    "q0", "q1",
    "p0", "p1", "p2", "p3",
    "a0", "a1", "a2", "a3", "a4", "a5",
    "b0", "b1", "b2", "b3", "b4", "b5",
    "c"
};

/**************************** Class definitions *****************************/

/*! \brief
 * Heir of ScatteringFactor, for storing parameters as polynomial+fourier series.
 */
class FourierScatteringFactor : public ScatteringFactor
{
    private:
        //! Range of scattering angles q
        double q0_, q1_;

        //! Parameters describing the Fourier model of the scattering
        std::vector<double> p_, a_, b_;
    public:
        /*! \brief
         * Constructor
         *
         * \param[in] residue        The residue
         * \param[in] atom           The atom
         * \param[in] type           The type indicator
         */
        FourierScatteringFactor(const char *residue,
                                const char *atom,
                                const char *type);

        //! Destructor
        virtual ~FourierScatteringFactor() {}

        /*! \brief
         * Function to compute the scattering for a given q
         * \param[in] q scattering vector
         * \return scattering factor
         */
        virtual double computeScatteringFactor(double q) const;

        /*! \brief
         * Function to compute the scattering for a given q
         * \param[in] theta scattering angle
         * \param[in] lambda wave length
         * \return scattering factor
         */
        virtual double computeScatteringFactor(double theta, double lambda) const;

        /*! \brief
         * Function to load data from an array retrieved from an xml structure
         * \param[in] tree The libxml2 data structure
         */
        virtual void loadScatteringFactor(xmlNodePtr tree);

        /*! \brief
         * Function to store data in xml structure
         * \param[in] parent xml data structure
         */
        virtual void storeScatteringFactor(xmlNodePtr parent);
};

/*! \brief
 * Heir of ScatteringFactor, for storing Cromer-Mann parameters.
 */
class CromerMannScatteringFactor : public ScatteringFactor
{
    private:
        //! Parameters that describe the polynomial in the Cromer Mann description
        std::vector<double> a_, b_;
        //! Parameter that describe the polynomial in the Cromer Mann description
        double              c_;
    public:
        /*! \brief
         * Constructor
         *
         * \param[in] residue        The residue
         * \param[in] atom           The atom
         * \param[in] type           The type indicator
         */
        CromerMannScatteringFactor(const char *residue,
                                   const char *atom,
                                   const char *type);

        //! Destructor
        virtual ~CromerMannScatteringFactor() {}

        /*! \brief
         * Function to compute the scattering for a given q
         * \param[in] q scattering vector
         * \return scattering factor
         */
        virtual double computeScatteringFactor(double q) const;

        /*! \brief
         * Function to compute the scattering for a given q
         * \param[in] theta scattering angle
         * \param[in] lambda wave length
         * \return scattering factor
         */
        virtual double computeScatteringFactor(double theta, double lambda) const;

        /*! \brief
         * Function to load data from an array retrieved from an xml structure
         * \param[in] tree the libxml2 data structure
         */
        virtual void loadScatteringFactor(xmlNodePtr tree);

        /*! \brief
         * Function to store data in xml structure
         * \param[in] parent xml data structure
         */
        virtual void storeScatteringFactor(xmlNodePtr parent);
};

/**************************** XML help routines *****************************/

/************************************
 * I/O help routines using libxml2
 ************************************/

#define NN(x) (NULL != (x))

//! Strings corresponding to XML elements
/*static const char *xmltypes[] = {
    NULL,
    "XML_ELEMENT_NODE",
    "XML_ATTRIBUTE_NODE",
    "XML_TEXT_NODE",
    "XML_CDATA_SECTION_NODE",
    "XML_ENTITY_REF_NODE",
    "XML_ENTITY_NODE",
    "XML_PI_NODE",
    "XML_COMMENT_NODE",
    "XML_DOCUMENT_NODE",
    "XML_DOCUMENT_TYPE_NODE",
    "XML_DOCUMENT_FRAG_NODE",
    "XML_NOTATION_NODE",
    "XML_HTML_DOCUMENT_NODE",
    "XML_DTD_NODE",
    "XML_ELEMENT_DECL",
    "XML_ATTRIBUTE_DECL",
    "XML_ENTITY_DECL",
    "XML_NAMESPACE_DECL",
    "XML_XINCLUDE_START",
    "XML_XINCLUDE_END",
    "XML_DOCB_DOCUMENT_NODE",
   };
   #define NXMLTYPES (sizeof(xmltypes)/sizeof(xmltypes[0]))
 */
//! Auxiliary routine converting boolean to char *
static char *bool2char(bool bbb)
{
    if (bbb)
    {
        return strdup("true");
    }
    else
    {
        return strdup("false");
    }
}

//! Return true if the strings are equal (case-insensitive)
static bool str_equal(const char *a, const char *b)
{
    std::string sa(a), sb(b);
    std::transform(sa.begin(), sa.end(), sa.begin(), ::toupper);
    std::transform(sb.begin(), sb.end(), sb.begin(), ::toupper);
    return (sa.compare(sb) == 0);
}

//! Auxiliary routine to convert char * to boolean
static bool char2bool(const char *esf)
{
    if ((str_equal(esf, "TRUE")) ||
        (str_equal(esf, "YES")) ||
        (str_equal(esf, "ON")))
    {
        return true;
    }
    else if ((str_equal(esf, "FALSE")) ||
             (str_equal(esf, "NO")) ||
             (str_equal(esf, "OFF")))
    {
        return false;
    }
    return false;
}

//! Ugly global for CromerMann scattering factor type
static const char *cm = "CromerMann";

//! Ugly global for Fourier scattering factor type
static const char *fr = "Fourier";

//! Convert char * to enumerated
static esfType char2esfType(const char *esf)
{
    if (str_equal(esf, cm))
    {
        return esfCromerMann;
    }
    else if (str_equal(esf, fr))
    {
        return esfFourier;
    }
    else
    {
        char errbuf[100];
        sprintf(errbuf, "Invalid type %s", esf);
        GMX_THROW(InvalidInputError(errbuf));
    }
}

//! Convert enumerated to char *
static const char *esfType2char(esfType esf)
{
    if (esfCromerMann == esf)
    {
        return cm;
    }
    else if (esfFourier == esf)
    {
        return fr;
    }
    else
    {
        char errbuf[100];
        sprintf(errbuf, "Invalid type %d", (int) esf);
        GMX_THROW(InvalidInputError(errbuf));
    }
}

//! Fill buffer with spaces for debug output
static void fill_spaces(int n, char buf[], int maxindent)
{
    int i;

    /* Don't indent more than maxindent characters */
    for (i = 0; (i < std::min(n, maxindent-1)); i++)
    {
        buf[i] = ' ';
    }
    buf[i] = '\0';
}

//! Fetch the attributes for this element
static void get_attributes(FILE *fp, bool bZero, int indent, xmlAttrPtr attr, char *xbuf[])
{
    char  buf[100];
    char *attrname, *attrval;
    int   i, kkk;

    if (bZero)
    {
        for (i = 0; (i < exmlNR); i++)
        {
            xbuf[i] = NULL;
        }
    }

    while (NULL != attr)
    {
        attrname = (char *)attr->name;
        attrval  = (char *)attr->children->content;

        if ((kkk = find_elem(attrname, exmlNR, exml_names)) != -1)
        {
            if (attrval != NULL)
            {
                xbuf[kkk] = strdup(attrval);
            }
        }
        if (NULL != fp)
        {
            fill_spaces(indent, buf, 99);
            fprintf(fp, "%sProperty: '%s' Value: '%s'\n",
                    buf, attrname, attrval);
        }
        attr = attr->next;
    }
}

//! Process a whole XML tree using DOM (Document Object Model)
static void process_tree(FILE                  *fp,
                         xmlNodePtr             tree,
                         int                    indent,
                         ScatteringFactorTable &sf)
{
    while (tree != NULL)
    {
        if (fp)
        {
            fprintf(fp, "Node type %d encountered\n", (int)tree->type);
        }

        if (tree->type == XML_ELEMENT_NODE)
        {
            int   elem;

            if (-1 != (elem = find_elem((char *)tree->name, exmlNR, exml_names)))
            {
                if (fp)
                {
                    char  mybuf[100];
                    fill_spaces(indent, mybuf, 99);
                    fprintf(fp, "%sElement node name %s\n", mybuf, (char *)tree->name);
                }
                char *xbuf[exmlNR];

                get_attributes(fp, true, indent, tree->properties, xbuf);
                /* Done processing attributes for this element. Let's see if we still need
                 * to interpret them.
                 */
                for (int k = 0; (k < exmlNR); k++)
                {
                    if (NN(xbuf[k]))
                    {
                        //printf("KOKO xbuf[%d] (%s) = %s\n", k, exml_names[k], xbuf[k]);
                    }
                }
                switch (elem)
                {
                    case exmlSCATTERINGFACTORS:
                        if (NN(xbuf[exmlTYPE]) &&
                            NN(xbuf[exmlFORCEFIELD]) &&
                            NN(xbuf[exmlDISPLACEDSOLVENT]))
                        {
                            sf.setForceField(xbuf[exmlFORCEFIELD]);
                            sf.setSfType(char2esfType(xbuf[exmlTYPE]));
                            sf.setDisplacedSolvent(char2bool(xbuf[exmlDISPLACEDSOLVENT]));
                            if (NN(xbuf[exmlREFERENCE]))
                            {
                                sf.setReference(xbuf[exmlREFERENCE]);
                            }
                        }
                        break;
                    case exmlSCATTERINGFACTOR:
                        if (NN(xbuf[exmlRESIDUE]) &&
                            NN(xbuf[exmlATOM]) &&
                            NN(xbuf[exmlIDENTIFIER]))
                        {
                            switch (sf.getSfType())
                            {
                                case esfFourier:
                                {
                                    ScatteringFactorPointer
                                        temp(new FourierScatteringFactor(xbuf[exmlRESIDUE],
                                                                         xbuf[exmlATOM],
                                                                         xbuf[exmlIDENTIFIER]));
                                    temp->loadScatteringFactor(tree);
                                    sf.addSFP(move(temp));
                                    break;
                                }
                                case esfCromerMann:
                                {
                                    ScatteringFactorPointer
                                        temp(new CromerMannScatteringFactor(xbuf[exmlRESIDUE],
                                                                            xbuf[exmlATOM],
                                                                            xbuf[exmlIDENTIFIER]));
                                    temp->loadScatteringFactor(tree);
                                    sf.addSFP(move(temp));
                                    break;
                                }
                                default:
                                    GMX_THROW(InvalidInputError("Wrong type"));
                                    break;
                            }
                        }
                        break;
                    default:
                        if (NULL != debug)
                        {
                            fprintf(debug, "Unknown combination of attributes:\n");
                            for (int i = 0; (i < exmlNR); i++)
                            {
                                if (xbuf[i] != NULL)
                                {
                                    fprintf(debug, "%s = %s\n", exml_names[i], xbuf[i]);
                                }
                            }
                        }
                }

                /* Clean up */
                for (int i = 0; (i < exmlNR); i++)
                {
                    if (xbuf[i] != NULL)
                    {
                        sfree(xbuf[i]);
                        xbuf[i] = NULL;
                    }
                }
                if (tree->children)
                {
                    process_tree(fp, tree->children, indent+2, sf);
                }
            }
        }
        tree = tree->next;
    }
}

/**************************** Class implementations *****************************/

FourierScatteringFactor::FourierScatteringFactor(const char *residue,
                                                 const char *atom,
                                                 const char *type)
{
    residue_.assign(residue);
    atom_.assign(atom);
    type_ = atoi(type);

    q0_ = 0;
    q1_ = 0;
    p_.resize(4);
    a_.resize(6);
    b_.resize(6);
    b_[0] = 0;
}

double FourierScatteringFactor::computeScatteringFactor(double theta, double lambda) const
{
    double q;

    GMX_RELEASE_ASSERT((lambda <= 0),
                       "FourierScatteringFactor::computeScatteringFactor called with lambda <= 0");

    q = 4*M_PI*sin(theta)/lambda;
    return computeScatteringFactor(q);
}

double FourierScatteringFactor::computeScatteringFactor(double q) const
{
    unsigned int i;
    double       s = a_[0];

    // add the polynomial
    s += p_[0]*q*q*q + p_[1]*q*q + p_[2]*q;

    // add the fourier series
    for (i = 1; i < a_.size(); i++)
    {
        s += a_[i]*cos(2*M_PI*(i+1)*(q-q0_)/(q1_-q0_));
        s += b_[i]*sin(2*M_PI*(i+1)*(q-q0_)/(q1_-q0_));
    }
    return s;
}

/*! \brief
 * Super-handy routine turns ptr into a double anf throws in case of problems
 * \param[in] ptr  ptr to a string containing a double
 * \param[in] name description of the variable
 * \param[in] unit the physical unit of the data
 * \return   the double
 */
static double my_atof(char *ptr, const char *name, const char *unit)
{
    char errbuf[STRLEN];
    if (NULL == ptr)
    {
        sprintf(errbuf, "No value for variable %s", name);
        GMX_THROW(InvalidInputError(errbuf));
    }
    double d;
    if (sscanf(ptr, "%20lf", &d) != 1)
    {
        sprintf(errbuf, "No double value for variable %s in field %s", name, ptr);
        GMX_THROW(InvalidInputError(errbuf));
    }

    return convert2gmx(d, string2unit(unit));
}

void FourierScatteringFactor::loadScatteringFactor(xmlNodePtr parent)
{
    int        node;
    xmlNodePtr tree;
    char      *xbuf[exmlNR];

    tree = parent->children;
    while (NULL != tree)
    {
        // We need properties and content of the children to be non NULL
        if ((NULL != tree->name) &&
            (NULL != tree->properties) &&
            (NULL != tree->children) &&
            (NULL != tree->children->content) &&
            ((node = find_elem((char *)tree->name, exmlNR, exml_names)) != -1))
        {
            if (NULL != debug)
            {
                fprintf(debug, "Found node %s ", exml_names[node]);
                fprintf(debug, "tree->name = %s tree->children->content %s\n",
                        (char *) tree->name, (char *)tree->children->content);
            }
            if (NULL != tree->properties)
            {
                get_attributes(debug, true, 0, tree->properties, xbuf);
            }
            if (NN(xbuf[exmlUNIT]))
            {
                char  *ptr   = strdup((char *)tree->children->content);
                double value = my_atof(ptr, exml_names[node], xbuf[exmlUNIT]);
                switch (node)
                {
                    case exmlQ0:
                        q0_ = value;
                        break;
                    case exmlQ1:
                        q1_ = value;
                        break;
                    case exmlP0:
                    case exmlP1:
                    case exmlP2:
                    case exmlP3:
                        p_[node-exmlP0] = value;
                        break;
                    case exmlA0:
                    case exmlA1:
                    case exmlA2:
                    case exmlA3:
                    case exmlA4:
                    case exmlA5:
                        a_[node-exmlA0] = value;
                        break;
                    case exmlB0:
                    case exmlB1:
                    case exmlB2:
                    case exmlB3:
                    case exmlB4:
                    case exmlB5:
                        b_[node-exmlB0] = value;
                        break;
                    default:
                        printf("Unexpected child %s\n", exml_names[node]);
                }

                sfree(ptr);
                sfree(xbuf[exmlUNIT]);
                xbuf[exmlUNIT] = NULL;
            }
        }
        tree = tree->next;
    }

}

void FourierScatteringFactor::storeScatteringFactor(xmlNodePtr parent)
{
    xmlNodePtr child = add_xml_child(parent, exml_names[exmlSCATTERINGFACTOR]);
    add_xml_char(child, exml_names[exmlRESIDUE], residue().c_str());
    add_xml_char(child, exml_names[exmlATOM], atom().c_str());
    add_xml_int(child, exml_names[exmlIDENTIFIER], type());

    add_xml_char(add_xml_child_double(child, exml_names[exmlQ0],
                                      gmx2convert(q0_, eg2c1_Angstrom)),
                 exml_names[exmlUNIT], unit2string(eg2c1_Angstrom));
    add_xml_char(add_xml_child_double(child, exml_names[exmlQ1],
                                      gmx2convert(q1_, eg2c1_Angstrom)),
                 exml_names[exmlUNIT], unit2string(eg2c1_Angstrom));

    for (unsigned int i = 1; (i < p_.size()); i++)
    {
        add_xml_char(add_xml_child_double(child, exml_names[i+exmlP0], p_[i]),
                     exml_names[exmlUNIT], unit2string(eg2cElectron));
    }

    for (unsigned int i = 0; (i < a_.size()); i++)
    {
        add_xml_char(add_xml_child_double(child, exml_names[i+exmlA0], a_[i]),
                     exml_names[exmlUNIT], unit2string(eg2cElectron));
    }

    for (unsigned int i = 1; (i < b_.size()); i++)
    {
        add_xml_char(add_xml_child_double(child, exml_names[i+exmlB0], b_[i]),
                     exml_names[exmlUNIT], unit2string(eg2cElectron));
    }
}

CromerMannScatteringFactor::CromerMannScatteringFactor(const char *residue,
                                                       const char *atom,
                                                       const char *type)
{
    residue_.assign(residue);
    atom_.assign(atom);
    type_ = atoi(type);
    c_    = 0;

    a_.resize(4);
    b_.resize(4);
}

double CromerMannScatteringFactor::computeScatteringFactor(double theta, double lambda) const
{
    double q;

    GMX_RELEASE_ASSERT((lambda <= 0),
                       "CromerMannScatteringFactor::computeScatteringFactor called with lambda <= 0");

    // added 4pi here for q -AB
    q = 4*M_PI*sin(theta)/lambda;

    return computeScatteringFactor(q);
}


double CromerMannScatteringFactor::computeScatteringFactor(double q) const
{
    int    i;
    double cm = c_;
    double q4 = q/(4*M_PI);

    for (i = 0; (i < 4); i++)
    {
        cm += a_[i]*exp(-b_[i]*q4*q4);
    }
    return cm;
}

void CromerMannScatteringFactor::loadScatteringFactor(xmlNodePtr tree)
{
    int   node;
    char *xbuf[exmlNR];

    printf("CromerMannScatteringFactor::loadScatteringFactor\n");
    while (NULL != tree)
    {
        if (((node = find_elem((char *)tree->name, exmlNR, exml_names)) != -1) &&
            (NULL != tree->children) &&
            (NULL != tree->children->content) &&
            (NULL != tree->children->properties))
        {
            get_attributes(debug, true, 0, tree->children->properties, xbuf);
            if (NN(xbuf[exmlUNIT]))
            {
                char  *ptr   = strdup((char *)tree->children->content);
                double value = my_atof(ptr, exml_names[node], xbuf[exmlUNIT]);
                switch (node)
                {
                    case exmlA0:
                    case exmlA1:
                    case exmlA2:
                    case exmlA3:
                        a_[node-exmlA0] = value;
                        break;
                    case exmlB0:
                    case exmlB1:
                    case exmlB2:
                    case exmlB3:
                        b_[node-exmlB0] = value;
                        break;
                    case exmlC:
                        c_   = value;
                        break;
                    default:
                        printf("Unexpected child %s reading CromerMann\n", exml_names[node]);
                }

                sfree(ptr);
                sfree(xbuf[exmlUNIT]);
            }
        }
        tree = tree->next;
    }
}

void CromerMannScatteringFactor::storeScatteringFactor(xmlNodePtr parent)
{
    xmlNodePtr child = add_xml_child(parent, exml_names[exmlSCATTERINGFACTOR]);

    add_xml_char(child, exml_names[exmlRESIDUE], residue().c_str());
    add_xml_char(child, exml_names[exmlATOM], atom().c_str());
    add_xml_int(child, exml_names[exmlIDENTIFIER], type());

    for (unsigned int i = 0; (i < a_.size()); i++)
    {
        add_xml_char(add_xml_child_double(child, exml_names[i+exmlA0], a_[i]),
                     exml_names[exmlUNIT], unit2string(eg2cElectron));
    }

    for (unsigned int i = 0; (i < b_.size()); i++)
    {
        add_xml_char(add_xml_child_double(child, exml_names[i+exmlB0], b_[i]),
                     exml_names[exmlUNIT], unit2string(eg2cElectron));
    }

    add_xml_char(add_xml_child_double(child, exml_names[exmlC], c_),
                 exml_names[exmlUNIT], unit2string(eg2cElectron));
}

//! Compare the types in two ScatteringFactors
static bool ScatteringFactorComp(const ScatteringFactorPointer &a, const ScatteringFactorPointer &b)
{
    return (a->type() < b->type());
}

ScatteringFactorTable::ScatteringFactorTable()
{
    setSfType(esfFourier);
    setDisplacedSolvent(false);
}

bool ScatteringFactorTable::read(const char *datafile)
{
    xmlDocPtr  doc;
    char      *fn2 = NULL;

    if (NULL != datafile)
    {
        fn2 = (char *)gmxlibfn(datafile);
    }
    if (NULL == fn2)
    {
        fn2 = (char *)gmxlibfn("scatteringfactor.xml");
    }
    if (NULL != debug)
    {
        fprintf(debug, "Opening library file %s\n", fn2);
    }
    xmlDoValidityCheckingDefaultValue = 0;
    if ((doc = xmlParseFile(fn2)) == NULL)
    {
        fprintf(stderr, "\nError reading XML file %s. Run a syntax checker such as nsgmls.\n",
                fn2);
    }
    else
    {
        process_tree(debug, doc->children, 0, *this);
        xmlFreeDoc(doc);
    }
    sfree(fn2);

    sort(cm_.begin(), cm_.end(), ScatteringFactorComp);

    printf("There are %d entries in the scatteringfactortable\n",
           (int) cm_.size());
    return true;
}

bool ScatteringFactorTable::write(const char *datafile)
{
    if (NULL != debug)
    {
        fprintf(debug, "Going to dump ScatteringFactorTable to %s\n", datafile);
    }
    xmlDocPtr      doc;
    xmlNodePtr     gmx_ptr, sfac_ptr;
    xmlNsPtr       gmx_ns;
    const xmlChar *gmx                = (const xmlChar *) exml_names[exmlGROMACS];
    //const xmlChar *scatteringfactors  = (const xmlChar *) exml_names[exmlSCATTERINGFACTORS];
    const xmlChar *href               = (const xmlChar *) "http://www.gromacs.org/schemas";
    const xmlChar *prefix             = (const xmlChar *) "gmx";

    if (NULL == (doc = xmlNewDoc((xmlChar *)"1.0")))
    {
        GMX_THROW(FileIOError("Creating XML document"));
    }
    if (NULL == (gmx_ptr = xmlNewDocNode(doc, NULL, gmx, NULL)))
    {
        GMX_THROW(APIError("Creating root element"));
    }
    if (NULL == (gmx_ns = xmlNewNs(gmx_ptr, href, prefix)))
    {
        GMX_THROW(APIError("Creating namespace"));
    }

    xmlDocSetRootElement(doc, gmx_ptr);
    sfac_ptr = add_xml_child(gmx_ptr, exml_names[exmlSCATTERINGFACTORS]);

    /* Add molecule definitions */
    add_xml_char(sfac_ptr, exml_names[exmlSOURCE], "photons");
    add_xml_char(sfac_ptr, exml_names[exmlTYPE],
                 esfType2char(getSfType()));
    add_xml_char(sfac_ptr, exml_names[exmlFORCEFIELD],
                 getForceField().c_str());
    add_xml_char(sfac_ptr, exml_names[exmlDISPLACEDSOLVENT],
                 bool2char(getDisplacedSolvent()));
    add_xml_char(sfac_ptr, exml_names[exmlREFERENCE],
                 getReference().c_str());
    for (std::vector<ScatteringFactorPointer>::iterator sfp = beginSFP(); (sfp < endSFP()); ++sfp)
    {
        (*sfp)->storeScatteringFactor(sfac_ptr);
    }

    xmlSetDocCompressMode(doc, 0);
    xmlIndentTreeOutput = 1;
    if (xmlSaveFormatFileEnc(datafile, doc, "UTF-8", 2) == 0)
    {
        GMX_THROW(FileIOError("Saving scatteringfactor file"));
    }
    xmlFreeDoc(doc);

    return true;
}

double ScatteringFactorTable::computeScatteringFactor(int type, double q) const
{
    unsigned int i = 0;
    // cm_ is a vector of scattering_factor:s.
    for (i = 0; (i < cm_.size()); i++)
    {
        // class scattering_factor has a function type() that returns
        // it's private variable type_:
        if (cm_[i]->type() == type)
        {
            // class scattering_factor has a method calc(),
            // so call that when the right entry is found.
            return cm_[i]->computeScatteringFactor(q);
        }
    }
    return -1;
}

double ScatteringFactorTable::computeScatteringFactor(const std::string &residue,
                                                      const std::string &atom,
                                                      double             q) const
{
    unsigned int i = 0;
    for (i = 0; (i < cm_.size()); i++)
    {
        if ((cm_[i]->atom()    == atom) &&
            (cm_[i]->residue() == residue))
        {
            return cm_[i]->computeScatteringFactor(q);
        }
    }
    return -1;
}

int ScatteringFactorTable::maxType() const
{
    unsigned int i         = 0;
    int          maxnumber = cm_[0]->type();
    for (i = 0; (i < cm_.size()); i++)
    {
        maxnumber = std::max(maxnumber, cm_[i]->type());
    }
    return maxnumber;
}

}
