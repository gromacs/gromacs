/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#include "maths.h"
#include "smalloc.h"
#include "xml_util.h"
#include "scattering_factors.h"
#include "gromacs/fileio/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"

extern int xmlDoValidityCheckingDefaultValue;

using namespace std;

namespace gmx
{

//! Enumerate to distinguish the different input types
enum {
    exmlSFACTORS,
    exmlSFACTOR,
    exmlDISPLACEDSOLVENT,
    exmlFORCEFIELD,
    exmlREFERENCE,
    exmlRESIDUE,
    exmlATOM,
    exmlTYPE,
    exmlUNIT,
    exmlAA,
    exmlQ0, exmlQRANGE,
    exmlP0, exmlP1, exmlP2,
    exmlA0, exmlA1, exmlA2, exmlA3, exmlA4,
    exmlB0, exmlB1, exmlB2, exmlB3, exmlB4,
    exmlC,
    exmlNR
};

//! Strings with corresponding names
static const char *exml_names[exmlNR] =
{
    "sfactors",
    "sfactor",
    "displaced_solvent",
    "force_field",
    "reference",
    "residue",
    "atom",
    "type",
    "unit",
    "aa",
    "q0", "qrange",
    "p0", "p1", "p2",
    "a0", "a1", "a2", "a3", "a4",
    "b0", "b1", "b2", "b3", "b4",
    "c"
};

/*! \brief
 * Heir of ScatteringFactor, for storing parameters as polynomial+fourier series.
 */
class FourierSfactor : public ScatteringFactor
{
    private:
        //! Parameters describing the Fourier descrption of the scattering
        double aa_, q0_, qrange_;

        //! More parameters
        std::vector<double> p_, a_, b_;
    public:
        /*! \brief
         * Constructor
         *
         * \param[in] residue        The residue
         * \param[in] atom           The atom
         * \param[in] type           The type indicator
         * \param[in] unit           The unit of the values
         */
        FourierSfactor(const char *residue,
                       const char *atom,
                       const char *type,
                       const char *unit);

        //! Destructor
        virtual ~FourierSfactor() {}

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
         * \param[in] xbuf array of strings
         */
        virtual void loadScatteringFactor(char *xbuf[]);

        /*! \brief
         * Function to store data in xml structure
         * \param[in] parent xml data structure
         */
        virtual void storeScatteringFactor(xmlNodePtr parent);
};

FourierSfactor::FourierSfactor(const char *residue,
                               const char *atom,
                               const char *type,
                               const char *unit)
{
    residue_.assign(residue);
    atom_.assign(atom);
    if (NULL == unit)
    {
        unit = "1/A";
    }
    unit_.assign(unit);
    
    type_ = atoi(type);

    aa_     = 0;
    q0_     = 0;
    qrange_ = 0;
    p_.resize(3);
    a_.resize(5);
    b_.resize(5);
}

double FourierSfactor::computeScatteringFactor(double theta, double lambda) const
{
    double q;

    GMX_RELEASE_ASSERT((lambda <= 0),
                       "FourierSfactor::computeScatteringFactor called with lambda <= 0");

    q = 4*M_PI*sin(theta)/lambda; 
    return computeScatteringFactor(q);
}

double FourierSfactor::computeScatteringFactor(double q) const
{
    unsigned int i;
    double       s = aa_;

    // q must be in inverse Angstroms
    q *= 0.1;

    // add the polynomial
    s += p_[0]*q*q*q + p_[1]*q*q + p_[2]*q;

    // add the fourier series
    for (i = 0; i < a_.size(); i++)
    {
        s += a_[i]*cos(2*M_PI*(i+1)*(q-q0_)/qrange_);
        s += b_[i]*sin(2*M_PI*(i+1)*(q-q0_)/qrange_);
    }
    return s;
}

/*! \brief
 * Super-handy routine turns ptr into a double anf throws in case of problems
 * \param[in] ptr  ptr to a string containing a double
 * \param[in] name description of the variabl
 * \return   the double
 */
static double my_atof(char *ptr, const char *name)
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
    return d;
}

void FourierSfactor::loadScatteringFactor(char *xbuf[])
{
    q0_     = my_atof(xbuf[exmlQ0], exml_names[exmlQ0]);
    qrange_ = my_atof(xbuf[exmlQRANGE], exml_names[exmlQRANGE]);
    p_[0]   = my_atof(xbuf[exmlP0], exml_names[exmlP0]);
    p_[1]   = my_atof(xbuf[exmlP1], exml_names[exmlP1]);
    p_[2]   = my_atof(xbuf[exmlP2], exml_names[exmlP2]);
    aa_     = my_atof(xbuf[exmlAA], exml_names[exmlAA]);
    a_[0]   = my_atof(xbuf[exmlA0], exml_names[exmlA0]);
    a_[1]   = my_atof(xbuf[exmlA1], exml_names[exmlA1]);
    a_[2]   = my_atof(xbuf[exmlA2], exml_names[exmlA2]);
    a_[3]   = my_atof(xbuf[exmlA3], exml_names[exmlA3]);
    a_[4]   = my_atof(xbuf[exmlA4], exml_names[exmlA4]);
    b_[0]   = my_atof(xbuf[exmlB0], exml_names[exmlB0]);
    b_[1]   = my_atof(xbuf[exmlB1], exml_names[exmlB1]);
    b_[2]   = my_atof(xbuf[exmlB2], exml_names[exmlB2]);
    b_[3]   = my_atof(xbuf[exmlB3], exml_names[exmlB3]);
    b_[4]   = my_atof(xbuf[exmlB4], exml_names[exmlB4]);
}

void FourierSfactor::storeScatteringFactor(xmlNodePtr parent)
{
    xmlNodePtr child = add_xml_child(parent, exml_names[exmlSFACTOR]);
    add_xml_char(child, exml_names[exmlRESIDUE], residue().c_str());
    add_xml_char(child, exml_names[exmlATOM], atom().c_str());
    add_xml_char(child, exml_names[exmlUNIT], unit().c_str());
    add_xml_int(child, exml_names[exmlTYPE], type());
    add_xml_child_double(child, exml_names[exmlAA], aa_);
    add_xml_child_double(child, exml_names[exmlQ0], q0_);
    add_xml_child_double(child, exml_names[exmlQRANGE], qrange_);
    add_xml_child_double(child, exml_names[exmlP0], p_[0]);
    add_xml_child_double(child, exml_names[exmlP1], p_[1]);
    add_xml_child_double(child, exml_names[exmlP2], p_[2]);
    add_xml_child_double(child, exml_names[exmlA0], a_[0]);
    add_xml_child_double(child, exml_names[exmlA1], a_[1]);
    add_xml_child_double(child, exml_names[exmlA2], a_[2]);
    add_xml_child_double(child, exml_names[exmlA3], a_[3]);
    add_xml_child_double(child, exml_names[exmlA4], a_[4]);
    add_xml_child_double(child, exml_names[exmlB0], b_[0]);
    add_xml_child_double(child, exml_names[exmlB1], b_[1]);
    add_xml_child_double(child, exml_names[exmlB2], b_[2]);
    add_xml_child_double(child, exml_names[exmlB3], b_[3]);
    add_xml_child_double(child, exml_names[exmlB4], b_[4]);
}

/*! \brief
 * Heir of ScatteringFactor, for storing Cromer-Mann parameters.
 */
class CromerMannSfactor : public ScatteringFactor
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
         * \param[in] unit           The unit of the values
         */
        CromerMannSfactor(const char *residue,
                          const char *atom,
                          const char *type,
                          const char *unit);

        //! Destructor
        virtual ~CromerMannSfactor() {}

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
         * \param[in] xbuf array of strings
         */
        virtual void loadScatteringFactor(char *xbuf[]);

        /*! \brief
         * Function to store data in xml structure
         * \param[in] parent xml data structure
         */
        virtual void storeScatteringFactor(xmlNodePtr parent);
};


CromerMannSfactor::CromerMannSfactor(const char *residue,
                                     const char *atom,
                                     const char *type,
                                     const char *unit)
{
    residue_.assign(residue);
    atom_.assign(atom);
    if (NULL == unit)
    {
        unit = "1/A";
    }
    unit_.assign(unit);
    type_ = atoi(type);
    c_    = 0;
    
    a_.resize(4);
    b_.resize(4);
}

double CromerMannSfactor::computeScatteringFactor(double theta, double lambda) const
{
    double q;

    GMX_RELEASE_ASSERT((lambda <= 0),
                       "CromerMannSfactor::computeScatteringFactor called with lambda <= 0");

    q = 4*M_PI*sin(theta)/lambda; // added 4pi here for q -AB
    return computeScatteringFactor(q);
}


double CromerMannSfactor::computeScatteringFactor(double q) const
{
    int    i;
    double cm = c_;
    double q4 = q/(4*M_PI);

    // q must be in inverse Angstroms
    q4 *= 0.1;

    for (i = 0; (i < 4); i++)
    {
        cm += a_[i]*exp(-b_[i]*q4*q4);
    }
    return cm;
}

void CromerMannSfactor::loadScatteringFactor(char *xbuf[])
{
    a_[0]   = my_atof(xbuf[exmlA0], exml_names[exmlA0]);
    a_[1]   = my_atof(xbuf[exmlA1], exml_names[exmlA1]);
    a_[2]   = my_atof(xbuf[exmlA2], exml_names[exmlA2]);
    a_[3]   = my_atof(xbuf[exmlA3], exml_names[exmlA3]);
    b_[0]   = my_atof(xbuf[exmlB0], exml_names[exmlB0]);
    b_[1]   = my_atof(xbuf[exmlB1], exml_names[exmlB1]);
    b_[2]   = my_atof(xbuf[exmlB2], exml_names[exmlB2]);
    b_[3]   = my_atof(xbuf[exmlB3], exml_names[exmlB3]);
    c_      = my_atof(xbuf[exmlC], exml_names[exmlC]);
}

void CromerMannSfactor::storeScatteringFactor(xmlNodePtr parent)
{
    xmlNodePtr child = add_xml_child(parent, exml_names[exmlSFACTOR]);
    add_xml_char(child, exml_names[exmlRESIDUE], residue().c_str());
    add_xml_char(child, exml_names[exmlATOM], atom().c_str());
    add_xml_char(child, exml_names[exmlUNIT], atom().c_str());
    add_xml_int(child, exml_names[exmlTYPE], type());
    add_xml_child_double(child, exml_names[exmlA0], a_[0]);
    add_xml_child_double(child, exml_names[exmlA1], a_[1]);
    add_xml_child_double(child, exml_names[exmlA2], a_[2]);
    add_xml_child_double(child, exml_names[exmlA3], a_[3]);
    add_xml_child_double(child, exml_names[exmlB0], b_[0]);
    add_xml_child_double(child, exml_names[exmlB1], b_[1]);
    add_xml_child_double(child, exml_names[exmlB2], b_[2]);
    add_xml_child_double(child, exml_names[exmlB3], b_[3]);
    add_xml_child_double(child, exml_names[exmlC], c_);
}

/************************************
 * I/O help routines using libxml2
 ************************************/

#define NN(x) (NULL != (x))

//! Strings corresponding to XML elements
static const char *xmltypes[] = {
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
    fprintf(stderr, "Got esf = %d\n", (int) esf);
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

//! Process XML children for elementes that have such
static void process_children(xmlNodePtr tree, char *xbuf[])
{
    int node;

    while (NULL != tree)
    {
        if (((node = find_elem((char *)tree->name, exmlNR, exml_names)) != -1) &&
            (NULL != tree->children) &&
            (NULL != tree->children->content))
        {
            if (NULL == xbuf[node])
            {
                xbuf[node] = strdup((char *)tree->children->content);
            }
        }
        tree = tree->next;
    }
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
    int   elem;
    char  buf[100];
    int   i;
    char *xbuf[exmlNR];

    while (tree != NULL)
    {
        if (fp)
        {
            fprintf(fp, "Node type %d encountered\n", (int)tree->type);
        }

        if (tree->type == XML_ELEMENT_NODE)
        {
            if (-1 != (elem = find_elem((char *)tree->name, exmlNR, exml_names)))
            {
                if (fp)
                {
                    fill_spaces(indent, buf, 99);
                    fprintf(fp, "%sElement node name %s\n", buf, (char *)tree->name);
                }
                get_attributes(fp, true, indent, tree->properties, xbuf);
                /* Done processing attributes for this element. Let's see if we still need
                 * to interpret them.
                 */
                switch (elem)
                {
                    case exmlSFACTORS:
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
                    case exmlSFACTOR:
                        if (NN(xbuf[exmlRESIDUE]) &&
                            NN(xbuf[exmlATOM]) &&
                            NN(xbuf[exmlTYPE]))
                        {
                            process_children(tree->children, xbuf);
                            switch (sf.getSfType())
                            {
                                case esfFourier:
                                {
                                    ScatteringFactorPointer
                                        temp(new FourierSfactor(xbuf[exmlRESIDUE],
                                                                xbuf[exmlATOM],
                                                                xbuf[exmlTYPE],
                                                                xbuf[exmlUNIT]));
                                    temp->loadScatteringFactor(xbuf);
                                    sf.cm_.push_back(move(temp));
                                    break;
                                }
                                case esfCromerMann:
                                {
                                    ScatteringFactorPointer
                                        temp(new CromerMannSfactor(xbuf[exmlRESIDUE],
                                                                   xbuf[exmlATOM],
                                                                   xbuf[exmlTYPE],
                                                                   xbuf[exmlUNIT]));
                                    temp->loadScatteringFactor(xbuf);
                                    sf.cm_.push_back(move(temp));
                                    break;
                                }
                                default:
                                    break;
                            }
                        }
                        break;
                    default:
                        if (NULL != debug)
                        {
                            fprintf(debug, "Unknown combination of attributes:\n");
                            for (i = 0; (i < exmlNR); i++)
                            {
                                if (xbuf[i] != NULL)
                                {
                                    fprintf(debug, "%s = %s\n", exml_names[i], xbuf[i]);
                                }
                            }
                        }
                }

                /* Clean up */
                for (i = 0; (i < exmlNR); i++)
                {
                    if (xbuf[i] != NULL)
                    {
                        sfree(xbuf[i]);
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

//! Compare the types in two ScatteringFactors
static bool sfactor_comp(const ScatteringFactorPointer &a, const ScatteringFactorPointer &b)
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
    if (1)
    {
        xmlDocPtr  doc;
        char      *fn2 = NULL;

        if (NULL != datafile)
        {
            fn2 = (char *)gmxlibfn(datafile);
        }
        if (NULL == fn2)
        {
            fn2 = (char *)gmxlibfn("sfactor.xml");
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
            sfree(fn2);
        }
        else
        {
            process_tree(debug, doc->children, 0, *this);
            xmlFreeDoc(doc);
        }
        sfree(fn2);
    }
    sort(cm_.begin(), cm_.end(), sfactor_comp);

    return true;
}

bool ScatteringFactorTable::write(const char *datafile)
{
    if (NULL != debug)
    {
        fprintf(debug, "Going to dump ScatteringFactorTable to %s\n", datafile);
    }
    xmlDocPtr   doc;
    xmlDtdPtr   dtd;
    xmlNodePtr  myroot;
    xmlChar    *libdtdname, *dtdname, *gmx;

    gmx        = (xmlChar *) "sfactors";
    dtdname    = (xmlChar *) "sfactors.dtd";
    libdtdname = dtdname;

    if ((doc = xmlNewDoc((xmlChar *)"1.0")) == NULL)
    {
        GMX_THROW(FileIOError("Creating XML document"));
    }
    if ((dtd = xmlCreateIntSubset(doc, dtdname, libdtdname, dtdname)) == NULL)
    {
        GMX_THROW(FileIOError("Creating XML DTD"));
    }
    if ((myroot = xmlNewDocNode(doc, NULL, gmx, NULL)) == NULL)
    {
        GMX_THROW(APIError("Creating root element"));
    }
    dtd->next    = myroot;
    myroot->prev = (xmlNodePtr) dtd;

    /* Add molecule definitions */
    add_xml_char(myroot, exml_names[exmlTYPE],
                 esfType2char(getSfType()));
    add_xml_char(myroot, exml_names[exmlFORCEFIELD],
                 getForceField().c_str());
    add_xml_char(myroot, exml_names[exmlDISPLACEDSOLVENT],
                 bool2char(getDisplacedSolvent()));
    add_xml_char(myroot, exml_names[exmlREFERENCE],
                 getReference().c_str());
    for (std::vector<ScatteringFactorPointer>::iterator sfp = beginSFP(); (sfp < endSFP()); sfp++)
    {
        (*sfp)->storeScatteringFactor(myroot);
    }

    xmlSetDocCompressMode(doc, 0);
    xmlIndentTreeOutput = 1;
    if (xmlSaveFormatFileEnc(datafile, doc, "UTF-8", 2) == 0)
    {
        GMX_THROW(FileIOError("Saving sfactor file"));
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
