/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include "poldata_xml.h"

#include <cstdlib>
#include <cstring>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#include "poldata.h"
#include "poldata-low.h"
#include "xml_util.h"

extern int xmlDoValidityCheckingDefaultValue;
namespace alexandria
{

#define NN(x) (0 != (x.size()))

const char *xmltypes[] = {
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
    "XML_XINCLUDE_END"
};
#define NXMLTYPES sizeof(xmltypes)/sizeof(xmltypes[0])

enum {
    exmlGENTOP, exmlREFERENCE,
    exmlATOMTYPES, exmlATOMTYPE,
    exmlGT_FORCEFIELD, exmlPOLAR_UNIT, exmlCOMB_RULE, exmlNEXCL,
    exmlFUDGEQQ, exmlFUDGELJ,
    exmlPOLTYPES, exmlPOLTYPE, exmlPTYPE,
    exmlELEM, exmlNAME, exmlDESC,
    exmlATYPE, exmlMILLER, exmlVALENCE, exmlBOSQUE,
    exmlBTYPE,
    exmlNEIGHBORS, exmlAROMATIC,
    exmlGEOMETRY, exmlNUMBONDS, exmlPOLARIZABILITY, exmlSIGPOL, exmlVDWPARAMS, exmlEREF,
    exmlFUNCTION,
    exmlATOM1, exmlATOM2, exmlATOM3, exmlATOM4,
    exmlSIGMA, exmlBONDORDER, exmlPARAMS,
    exmlLENGTH, exmlLENGTH_UNIT, exmlNTRAIN,
    exmlANGLE, exmlANGLE_UNIT,
    exmlGT_BONDS, exmlGT_BOND,
    exmlGT_ANGLES, exmlGT_ANGLE,
    exmlGT_DIHEDRALS, exmlGT_DIHEDRAL,
    exmlGT_IMPROPERS, exmlGT_IMPROPER,
    exmlBSATOMS, exmlBSATOM,
    exmlMILATOMS, exmlTAU_UNIT, exmlAHP_UNIT,
    exmlMILATOM, exmlMILNAME, exmlALEXANDRIA_EQUIV,
    exmlATOMNUMBER, exmlTAU_AHC, exmlALPHA_AHP,
    exmlSYMMETRIC_CHARGES, exmlSYM_CHARGE,
    exmlCENTRAL, exmlATTACHED, exmlNUMATTACH,
    exmlEEMPROPS, exmlEEMPROP, exmlMODEL, exmlJ0, exmlCHI0, exmlZETA, exmlROW,
    exmlEEMPROP_REF, exmlEPREF, exmlCHARGES,
    exmlNR
};

const char * exml_names[exmlNR] = {
    "gentop", "reference",
    "atomtypes", "atomtype", "forcefield", "polarizability_unit", "combination_rule", "nexclusions",
    "fudgeQQ", "fudgeLJ",
    "poltypes", "poltype", "ptype",
    "elem", "name", "description",
    "atype", "miller", "valence", "bosque",
    "btype",
    "neighbors", "aromatic",
    "geometry", "numbonds", "polarizability", "sigma_pol", "vdwparams", "ref_enthalpy",
    "function",
    "atom1", "atom2", "atom3", "atom4",
    "sigma", "bondorder", "params",
    "length", "length_unit", "ntrain",
    "angle", "angle_unit",
    "gt_bonds", "gt_bond",
    "gt_angles", "gt_angle",
    "gt_dihedrals", "gt_dihedral",
    "gt_impropers", "gt_improper",
    "bsatoms", "bsatom",
    "milatoms", "tau_ahc_unit", "alpha_ahp_unit",
    "milatom", "milname", "alexandria_equiv",
    "atomnumber", "tau_ahc", "alpha_ahp",
    "symmetric_charges", "sym_charge",
    "central", "attached", "numattach",
    "eemprops", "eemprop", "model", "jaa", "chi", "zeta", "row",
    "eemprop_ref", "epref", "charge"
};

static void sp(int n, char buf[], int maxindent)
{
    int i;
    if (n >= maxindent)
    {
        n = maxindent-1;
    }

    /* Don't indent more than maxindent characters */
    for (i = 0; (i < n); i++)
    {
        buf[i] = ' ';
    }
    buf[i] = '\0';
}

static double my_atof(const char *str)
{
    char   *ptr = NULL;
    double  d   = strtod(str, &ptr);
    GMX_RELEASE_ASSERT(ptr == NULL || strcmp(ptr, str) != 0, "Could not read double precision number");
    return d;
}

static void processAttr(FILE *fp, xmlAttrPtr attr, int elem,
                        int indent, Poldata &pd)
{
    std::string attrname, attrval;
    char        buf[100];
    int         i, kkk;
    std::string xbuf[exmlNR];

    for (i = 0; (i < exmlNR); i++)
    {
        xbuf[i] = "";
    }
    while (attr != NULL)
    {
        attrname.assign((char *)attr->name);
        attrval.assign((char *)attr->children->content);

#define atest(s) ((strcasecmp(attrname, s) == 0) && (attrval != NULL))
        kkk = find_elem((char *)attrname.c_str(), exmlNR, exml_names);
        if (-1 != kkk)
        {
            if (attrval.size() != 0)
            {
                xbuf[kkk] = attrval;
            }

            if (NULL != fp)
            {
                sp(indent, buf, 99);
                fprintf(fp, "%sProperty: '%s' Value: '%s'\n", buf, attrname.c_str(), attrval.c_str());
            }
        }
        else
        {
            fprintf(stderr, "Ignoring invalid attribute %s\n", attrname.c_str());
        }
        attr = attr->next;
#undef atest
    }
    /* Done processing attributes for this element. Let's see if we still need
     *  to interpret them.
     */
    switch (elem)
    {
        case exmlPOLTYPES:
            if (NN(xbuf[exmlPOLAR_UNIT]) && NN(xbuf[exmlREFERENCE]))
            {
                pd.setPolarUnit(xbuf[exmlPOLAR_UNIT]);
                pd.setPolarRef(xbuf[exmlREFERENCE]);
            }
            break;
        case exmlATOMTYPES:
            if (NN(xbuf[exmlGT_FORCEFIELD]))
            {
                pd.setForceField(xbuf[exmlGT_FORCEFIELD]);
            }
            if (NN(xbuf[exmlFUNCTION]))
            {
                pd.setVdwFunction(xbuf[exmlFUNCTION]);
            }
            if (NN(xbuf[exmlCOMB_RULE]))
            {
                pd.setCombinationRule(xbuf[exmlCOMB_RULE]);
            }
            if (NN(xbuf[exmlNEXCL]))
            {
                pd.setNexcl(atoi(xbuf[exmlNEXCL].c_str()));
            }
            if (NN(xbuf[exmlFUDGEQQ]))
            {
                pd.setFudgeQQ(atoi(xbuf[exmlFUDGEQQ].c_str()));
            }
            if (NN(xbuf[exmlFUDGELJ]))
            {
                pd.setFudgeLJ(atoi(xbuf[exmlFUDGELJ].c_str()));
            }
            break;
        case exmlBSATOMS:
            if (NN(xbuf[exmlPOLAR_UNIT]) &&
                NN(xbuf[exmlREFERENCE]))
            {
                pd.setBosqueFlags(xbuf[exmlPOLAR_UNIT], xbuf[exmlREFERENCE]);
            }
            break;
        case exmlGT_BONDS:
            if (NN(xbuf[exmlFUNCTION]) && NN(xbuf[exmlLENGTH_UNIT]))
            {
                GtBonds gtbs (xbuf[exmlFUNCTION], xbuf[exmlLENGTH_UNIT]);
                pd.addGtBonds(gtbs);
            }
            break;
        case exmlGT_ANGLES:
            if (NN(xbuf[exmlFUNCTION]) && NN(xbuf[exmlANGLE_UNIT]))
            {
                GtAngles gtas (xbuf[exmlFUNCTION], xbuf[exmlANGLE_UNIT]);
                pd.addGtAngles(gtas);
            }
            break;
        case exmlGT_DIHEDRALS:
            if (NN(xbuf[exmlFUNCTION]) && NN(xbuf[exmlANGLE_UNIT]))
            {
                GtDihedrals gtds (xbuf[exmlFUNCTION], xbuf[exmlANGLE_UNIT]);
                pd.addGtDihedrals(gtds);
            }
            break;
        case exmlMILATOMS:
            if (NN(xbuf[exmlTAU_UNIT]) &&
                NN(xbuf[exmlAHP_UNIT]) &&
                NN(xbuf[exmlREFERENCE]))
            {
                pd.setMillerFlags(xbuf[exmlTAU_UNIT], xbuf[exmlAHP_UNIT],
                                  xbuf[exmlREFERENCE]);
            }
            break;
        case exmlPOLTYPE:
            if (NN(xbuf[exmlPTYPE]) && NN(xbuf[exmlMILLER]) && NN(xbuf[exmlBOSQUE]) &&
                NN(xbuf[exmlPOLARIZABILITY]) && NN(xbuf[exmlSIGPOL]))
            {
                pd.addPtype(xbuf[exmlPTYPE],
                            xbuf[exmlMILLER],
                            xbuf[exmlBOSQUE],
                            my_atof(xbuf[exmlPOLARIZABILITY].c_str()),
                            my_atof(xbuf[exmlSIGPOL].c_str()));
            }
            break;
        case exmlATOMTYPE:
            if (NN(xbuf[exmlELEM]) &&
                NN(xbuf[exmlATYPE]) &&
                NN(xbuf[exmlBTYPE]) &&
                NN(xbuf[exmlPTYPE]) &&
                NN(xbuf[exmlVDWPARAMS]) &&
                NN(xbuf[exmlEREF]))
            {
                pd.addAtype(xbuf[exmlELEM],
                            xbuf[exmlDESC],
                            xbuf[exmlATYPE],
                            xbuf[exmlPTYPE],
                            xbuf[exmlBTYPE],
                            xbuf[exmlVDWPARAMS],
                            xbuf[exmlEREF]);
            }
            break;
        case exmlMILATOM:
            if (NN(xbuf[exmlMILNAME]) && NN(xbuf[exmlATOMNUMBER]) &&
                NN(xbuf[exmlTAU_AHC]) && NN(xbuf[exmlALPHA_AHP]))
            {
                pd.addMiller(xbuf[exmlMILNAME],
                             atoi(xbuf[exmlATOMNUMBER].c_str()),
                             my_atof(xbuf[exmlTAU_AHC].c_str()),
                             my_atof(xbuf[exmlALPHA_AHP].c_str()),
                             xbuf[exmlALEXANDRIA_EQUIV]);
            }
            break;
        case exmlBSATOM:
            if (NN(xbuf[exmlELEM]) && NN(xbuf[exmlPOLARIZABILITY]))
            {
                pd.addBosque( xbuf[exmlELEM], my_atof(xbuf[exmlPOLARIZABILITY].c_str()));
            }
            break;
        case exmlGT_BOND:
            if (NN(xbuf[exmlATOM1])     && NN(xbuf[exmlATOM2]) &&
                NN(xbuf[exmlLENGTH])    && NN(xbuf[exmlSIGMA]) && 
		NN(xbuf[exmlBONDORDER]) && NN(xbuf[exmlNTRAIN]))
            {
                auto& gtbs = pd.getLastGtBonds();
                gtbs.addBond(xbuf[exmlATOM1], xbuf[exmlATOM2],
                             my_atof(xbuf[exmlLENGTH].c_str()),
                             my_atof(xbuf[exmlSIGMA].c_str()),
                             atoi(xbuf[exmlNTRAIN].c_str()),
                             my_atof(xbuf[exmlBONDORDER].c_str()),
                             xbuf[exmlPARAMS]);
            }
            break;
        case exmlGT_ANGLE:
            if (NN(xbuf[exmlATOM1]) && NN(xbuf[exmlATOM2]) &&
                NN(xbuf[exmlATOM3]) && NN(xbuf[exmlANGLE]) && 
		NN(xbuf[exmlSIGMA]) && NN(xbuf[exmlNTRAIN]))
            {
                auto& gtas = pd.getLastGtAngles();
                gtas.addAngle(xbuf[exmlATOM1], xbuf[exmlATOM2], xbuf[exmlATOM3], 
			      my_atof(xbuf[exmlANGLE].c_str()),
                              my_atof(xbuf[exmlSIGMA].c_str()), 
			      atoi(xbuf[exmlNTRAIN].c_str()),
                              xbuf[exmlPARAMS].c_str());
            }
            break;
        case exmlGT_DIHEDRAL:
            if (NN(xbuf[exmlATOM1]) && NN(xbuf[exmlATOM2]) &&
                NN(xbuf[exmlATOM3]) && NN(xbuf[exmlATOM4]) &&
                NN(xbuf[exmlANGLE]) && NN(xbuf[exmlSIGMA]) &&
                NN(xbuf[exmlNTRAIN]))
            {
                auto& gtds = pd.getLastGtDihedrals();
                gtds.addDihedral(xbuf[exmlATOM1], xbuf[exmlATOM2],
                                 xbuf[exmlATOM3], xbuf[exmlATOM4],
                                 my_atof(xbuf[exmlANGLE].c_str()), 
				 my_atof(xbuf[exmlSIGMA].c_str()),
                                 atoi(xbuf[exmlNTRAIN].c_str()),
                                 xbuf[exmlPARAMS]);
            }
            break;
        case exmlSYM_CHARGE:
            if (NN(xbuf[exmlCENTRAL]) && NN(xbuf[exmlATTACHED]) &&
                NN(xbuf[exmlNUMATTACH]))
            {
                pd.addSymcharges(xbuf[exmlCENTRAL],
                                 xbuf[exmlATTACHED],
                                 atoi(xbuf[exmlNUMATTACH].c_str()));
            }
            break;
        case exmlEEMPROP:
            if (NN(xbuf[exmlMODEL]) && NN(xbuf[exmlNAME]) &&
                NN(xbuf[exmlCHI0])  && NN(xbuf[exmlJ0]) &&
                NN(xbuf[exmlZETA])  && NN(xbuf[exmlCHARGES]) &&
                NN(xbuf[exmlROW]))
            {
                Eemprops eep(name2eemtype(xbuf[exmlMODEL]),
                             xbuf[exmlNAME],
                             xbuf[exmlROW],
                             xbuf[exmlZETA],
                             xbuf[exmlCHARGES],
                             my_atof(xbuf[exmlJ0].c_str()),
                             my_atof(xbuf[exmlCHI0].c_str()) );
                pd.addEemprops(eep);
            }
            break;
        case exmlEEMPROP_REF:
            if (NN(xbuf[exmlMODEL]) && NN(xbuf[exmlEPREF]))
            {
                pd.setEpref(name2eemtype(xbuf[exmlMODEL]), xbuf[exmlEPREF]);
            }
            break;
        default:
            if (NULL != debug)
            {
                fprintf(debug, "Unknown combination of attributes:\n");
                for (i = 0; (i < exmlNR); i++)
                {
                    if (xbuf[i].size() != 0)
                    {
                        fprintf(debug, "%s = %s\n", exml_names[i], xbuf[i].c_str());
                    }
                }
            }
    }

    /* Clean up */
    /*for (i = 0; (i < exmlNR); i++)
       {
       if (xbuf[i].size() != 0)
        {
            sfree(xbuf[i]);
        }
       }*/

}

static void processTree(FILE *fp, xmlNodePtr tree, int indent,
                        Poldata &pd, gmx_atomprop_t aps)
{
    int           elem;
    char          buf[100];

    while (tree != NULL)
    {
        if (fp)
        {
            if ((tree->type > 0) && (tree->type < NXMLTYPES))
            {
                fprintf(fp, "Node type %s encountered with name %s\n",
                        xmltypes[tree->type], (char *)tree->name);
            }
            else
            {
                fprintf(fp, "Node type %d encountered\n", tree->type);
            }
        }

        if (tree->type == XML_ELEMENT_NODE)
        {
            elem = find_elem((char *)tree->name, exmlNR, exml_names);
            if (fp)
            {
                sp(indent, buf, 99);
                fprintf(fp, "%sElement node name %s\n", buf, (char *)tree->name);
            }
            if (-1 != elem)
            {
                if (elem != exmlGENTOP)
                {
                    processAttr(fp, tree->properties, elem, indent+2, pd);
                }

                if (tree->children)
                {
                    processTree(fp, tree->children, indent+2, pd, aps);
                }
            }
        }
        tree = tree->next;
    }
}

void readPoldata(const std::string &fileName,
                 Poldata           &pd,
                 gmx_atomprop_t     aps)
{
    xmlDocPtr   doc;
    std::string fn2;

    if (fileName.size() > 0)
    {
        const char *f = low_gmxlibfn(fileName.c_str(), TRUE, FALSE);
        if (NULL != f)
        {
            fn2.assign(f);
        }
    }
    if (!gmx_fexist(fn2.c_str()))
    {
        fn2 = gmxlibfn("alexandria.ff/gentop.dat");
    }

    if (NULL != debug)
    {
        fprintf(debug, "Opening library file %s\n", fn2.c_str());
    }
    xmlDoValidityCheckingDefaultValue = 0;
    doc = xmlParseFile(fn2.c_str());
    if (doc == NULL)
    {
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "Error reading XML file %s. Run a syntax checker such as nsgmls.",
                 fn2.c_str());
        GMX_THROW(gmx::FileIOError(buf));
    }

    pd.setFilename(fn2);
    processTree(debug, doc->children, 0, pd, aps);

    xmlFreeDoc(doc);

    if (NULL != debug)
    {
        writePoldata("pdout.dat", pd, false);
    }
}

static void addXmlPoldata(xmlNodePtr parent, const Poldata &pd)
{
    xmlNodePtr                    child, grandchild;
    int                           nexcl;
    std::string                   geometry, name,
                                  acentral, attached, tau_unit, ahp_unit,
                                  epref, desc, params;
    std::string  neighbors, zeta, qstr, rowstr;
    double       fudgeQQ, fudgeLJ;
    std::string  tmp, func, blu;
    child = add_xml_child(parent, exml_names[exmlATOMTYPES]);
    tmp   = pd.getForceField();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names[exmlGT_FORCEFIELD], tmp.c_str());
    }
    tmp = pd.getVdwFunction();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names[exmlFUNCTION], tmp.c_str());
    }
    tmp = pd.getCombinationRule();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names[exmlCOMB_RULE], tmp.c_str());
    }
    nexcl = pd.getNexcl();
    add_xml_int(child, exml_names[exmlNEXCL], nexcl);
    fudgeQQ = pd.getFudgeQQ();
    add_xml_double(child, exml_names[exmlFUDGEQQ], fudgeQQ);
    fudgeLJ = pd.getFudgeLJ();
    add_xml_double(child, exml_names[exmlFUDGELJ], fudgeLJ);
    {
        for (FfatypeConstIterator aType = pd.getAtypeBegin();
             aType != pd.getAtypeEnd(); aType++)
        {
            grandchild = add_xml_child(child, exml_names[exmlATOMTYPE]);
            add_xml_char(grandchild, exml_names[exmlELEM], aType->getElem().c_str());
            add_xml_char(grandchild, exml_names[exmlDESC], aType->getDesc().c_str());
            add_xml_char(grandchild, exml_names[exmlATYPE], aType->getType().c_str());
            add_xml_char(grandchild, exml_names[exmlPTYPE], aType->getPtype().c_str());
            add_xml_char(grandchild, exml_names[exmlBTYPE], aType->getBtype().c_str());
            add_xml_char(grandchild, exml_names[exmlVDWPARAMS], aType->getVdwparams().c_str());
            add_xml_char(grandchild, exml_names[exmlEREF], aType->getRefEnthalpy().c_str());
        }
    }
    child = add_xml_child(parent, exml_names[exmlPOLTYPES]);
    tmp   = pd.getPolarUnit();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names[exmlPOLAR_UNIT], tmp.c_str());
    }
    tmp = pd.getPolarRef();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names[exmlREFERENCE], tmp.c_str());
    }
    {
        for (PtypeConstIterator pType = pd.getPtypeBegin();
             pType != pd.getPtypeEnd(); pType++)
        {
            grandchild = add_xml_child(child, exml_names[exmlPOLTYPE]);
            add_xml_char(grandchild, exml_names[exmlPTYPE], pType->getType().c_str());
            add_xml_char(grandchild, exml_names[exmlMILLER], pType->getMiller().c_str());
            add_xml_char(grandchild, exml_names[exmlBOSQUE], pType->getBosque().c_str());
            add_xml_double(grandchild, exml_names[exmlPOLARIZABILITY], pType->getPolarizability());
            add_xml_double(grandchild, exml_names[exmlSIGPOL], pType->getSigPol());
        }
    }

    for (auto gtbs = pd.getBondsBegin(); gtbs != pd.getBondsEnd(); gtbs++)
    {
        child = add_xml_child(parent, exml_names[exmlGT_BONDS]);
        blu   = gtbs->getLengthUnit();
        if (blu.size() != 0)
        {
            add_xml_char(child, exml_names[exmlLENGTH_UNIT], blu.c_str());
        }
        func = gtbs->getBondFunction();
        if (func.size() != 0)
        {
            add_xml_char(child, exml_names[exmlFUNCTION], func.c_str());
        }
        for (auto gtb = gtbs->getBondBegin(); gtb != gtbs->getBondEnd(); gtb++)
        {
            grandchild = add_xml_child(child, exml_names[exmlGT_BOND]);
            add_xml_char(grandchild, exml_names[exmlATOM1], gtb->getAtom1().c_str());
            add_xml_char(grandchild, exml_names[exmlATOM2], gtb->getAtom2().c_str());
            add_xml_double(grandchild, exml_names[exmlLENGTH], gtb->getLength());
            add_xml_double(grandchild, exml_names[exmlSIGMA], gtb->getSigma());
            add_xml_int(grandchild, exml_names[exmlNTRAIN], gtb->getNtrain());
            add_xml_double(grandchild, exml_names[exmlBONDORDER], gtb->getBondorder());
            add_xml_char(grandchild, exml_names[exmlPARAMS], gtb->getParams().c_str());
        }
    }

    for (auto gtas = pd.getAnglesBegin(); gtas != pd.getAnglesEnd(); gtas++)
    {
        child = add_xml_child(parent, exml_names[exmlGT_ANGLES]);
        blu   = gtas->getAngleUnit();
        if (blu.size() != 0)
        {
            add_xml_char(child, exml_names[exmlANGLE_UNIT], blu.c_str());
        }
        func = gtas->getAngleFunction();
        if (func.size() != 0)
        {
            add_xml_char(child, exml_names[exmlFUNCTION], func.c_str());
        }
        for (auto gta = gtas->getAngleBegin(); gta != gtas->getAngleEnd(); gta++)
        {
            grandchild = add_xml_child(child, exml_names[exmlGT_ANGLE]);
            add_xml_char(grandchild, exml_names[exmlATOM1], gta->getAtom1().c_str());
            add_xml_char(grandchild, exml_names[exmlATOM2], gta->getAtom2().c_str());
            add_xml_char(grandchild, exml_names[exmlATOM3], gta->getAtom3().c_str());
            add_xml_double(grandchild, exml_names[exmlANGLE], gta->getAngle());
            add_xml_double(grandchild, exml_names[exmlSIGMA], gta->getSigma());
            add_xml_int(grandchild, exml_names[exmlNTRAIN], gta->getNtrain());
            add_xml_char(grandchild, exml_names[exmlPARAMS], gta->getParams().c_str());
        }
    }

    for (auto gtds = pd.getDihedralsBegin(); gtds != pd.getDihedralsEnd(); gtds++)
    {
        child = add_xml_child(parent, exml_names[exmlGT_DIHEDRALS]);
        blu   = gtds->getDihedralUnit();
        if (blu.size() != 0)
        {
            add_xml_char(child, exml_names[exmlANGLE_UNIT], blu.c_str());
        }
        func = gtds->getDihedralFunction();
        if (func.size() != 0)
        {
            add_xml_char(child, exml_names[exmlFUNCTION], func.c_str());
        }
        for (auto gtd = gtds->getDihedralBegin(); gtd != gtds->getDihedralEnd(); gtd++)
        {
            grandchild = add_xml_child(child, exml_names[exmlGT_DIHEDRAL]);
            add_xml_char(grandchild, exml_names[exmlATOM1], gtd->getAtom1().c_str());
            add_xml_char(grandchild, exml_names[exmlATOM2], gtd->getAtom2().c_str());
            add_xml_char(grandchild, exml_names[exmlATOM3], gtd->getAtom3().c_str());
            add_xml_char(grandchild, exml_names[exmlATOM4], gtd->getAtom4().c_str());
            add_xml_double(grandchild, exml_names[exmlANGLE], gtd->getDihedral());
            add_xml_double(grandchild, exml_names[exmlSIGMA], gtd->getSigma());
            add_xml_int(grandchild, exml_names[exmlNTRAIN], gtd->getNtrain());
            add_xml_char(grandchild, exml_names[exmlPARAMS], gtd->getParams().c_str());
        }
    }
    child = add_xml_child(parent, exml_names[exmlBSATOMS]);
    std::string ref;
    pd.getBosqueFlags(tmp, ref);
    add_xml_char(child, exml_names[exmlPOLAR_UNIT], tmp.c_str());
    add_xml_char(child, exml_names[exmlREFERENCE], ref.c_str());

    for (BosqueConstIterator bosque = pd.getBosqueBegin();
         bosque != pd.getBosqueEnd(); bosque++)
    {
        grandchild = add_xml_child(child, exml_names[exmlBSATOM]);
        add_xml_char(grandchild, exml_names[exmlELEM], bosque->getBosque().c_str());
        add_xml_double(grandchild, exml_names[exmlPOLARIZABILITY], bosque->getPolarizability());
    }
    child = add_xml_child(parent, exml_names[exmlMILATOMS]);
    std::string milref;
    pd.getMillerFlags(tau_unit, ahp_unit, milref);
    add_xml_char(child, exml_names[exmlTAU_UNIT], tau_unit.c_str());
    add_xml_char(child, exml_names[exmlAHP_UNIT], ahp_unit.c_str());
    add_xml_char(child, exml_names[exmlREFERENCE], milref.c_str());
    for (MillerConstIterator miller = pd.getMillerBegin();
         miller != pd.getMillerEnd(); miller++)
    {
        grandchild = add_xml_child(child, exml_names[exmlMILATOM]);
        add_xml_char(grandchild, exml_names[exmlMILNAME], miller->getMiller().c_str());
        add_xml_int(grandchild, exml_names[exmlATOMNUMBER], miller->getAtomnumber());
        add_xml_double(grandchild, exml_names[exmlTAU_AHC], miller->getTauAhc());
        add_xml_double(grandchild, exml_names[exmlALPHA_AHP], miller->getAlphaAhp());
        const std::string ae = miller->getAlexandriaEquiv();
        if (ae.size() > 0)
        {
            add_xml_char(grandchild, exml_names[exmlALEXANDRIA_EQUIV], ae.c_str());
        }
    }

    child = add_xml_child(parent, exml_names[exmlSYMMETRIC_CHARGES]);
    for (SymchargesConstIterator symcharges = pd.getSymchargesBegin();
         symcharges != pd.getSymchargesEnd(); symcharges++)
    {
        grandchild = add_xml_child(child, exml_names[exmlSYM_CHARGE]);
        add_xml_char(grandchild, exml_names[exmlCENTRAL], symcharges->getCentral().c_str());
        add_xml_char(grandchild, exml_names[exmlATTACHED], symcharges->getAttached().c_str());
        add_xml_int(grandchild, exml_names[exmlNUMATTACH], symcharges->getNumattach());
    }

    child = add_xml_child(parent, exml_names[exmlEEMPROPS]);

    for (EempropsConstIterator eep = pd.BeginEemprops();
         eep != pd.EndEemprops(); eep++)
    {
        ChargeDistributionModel model = eep->getEqdModel();

        grandchild = add_xml_child(child, exml_names[exmlEEMPROP]);
        add_xml_char(grandchild, exml_names[exmlMODEL], getEemtypeName(model));
        add_xml_char(grandchild, exml_names[exmlNAME], eep->getName());
        add_xml_double(grandchild, exml_names[exmlJ0], eep->getJ0());
        add_xml_double(grandchild, exml_names[exmlCHI0], eep->getChi0());
        add_xml_char(grandchild, exml_names[exmlZETA], eep->getZetastr());
        add_xml_char(grandchild, exml_names[exmlCHARGES], eep->getQstr());
        add_xml_char(grandchild, exml_names[exmlROW], eep->getRowstr());
    }
    for (auto eep = pd.epRefBegin(); eep < pd.epRefEnd(); ++eep)
    {
        grandchild = add_xml_child(child, exml_names[exmlEEMPROP_REF]);
        add_xml_char(grandchild, exml_names[exmlMODEL], getEemtypeName(eep->getEqdModel()));
        add_xml_char(grandchild, exml_names[exmlEPREF], eep->getEpref());
    }
}

void writePoldata(const std::string &fileName,
                  const Poldata     &pd,
                  bool               compress)
{
    xmlDocPtr   doc;
    xmlDtdPtr   dtd;
    xmlNodePtr  myroot;
    xmlChar    *libdtdname, *dtdname, *gmx;

    gmx        = (xmlChar *) "gentop";
    dtdname    = (xmlChar *) "gentop.dtd";
    libdtdname = dtdname;

    if ((doc = xmlNewDoc((xmlChar *)"1.0")) == NULL)
    {
        gmx_fatal(FARGS, "Creating XML document", "");
    }

    if ((dtd = xmlCreateIntSubset(doc, dtdname, libdtdname, dtdname)) == NULL)
    {
        gmx_fatal(FARGS, "Creating XML DTD", "");
    }

    if ((myroot = xmlNewDocNode(doc, NULL, gmx, NULL)) == NULL)
    {
        gmx_fatal(FARGS, "Creating root element", "");
    }
    dtd->next    = myroot;
    myroot->prev = (xmlNodePtr) dtd;

    /* Add molecule definitions */
    addXmlPoldata(myroot, pd);

    xmlSetDocCompressMode(doc, compress ? 1 : 0);
    xmlIndentTreeOutput = 1;
    if (xmlSaveFormatFileEnc(fileName.c_str(), doc, "ISO-8859-1", 2) == 0)
    {
        gmx_fatal(FARGS, "Saving file", fileName.c_str());
    }
    xmlFreeDoc(doc);
}

}
