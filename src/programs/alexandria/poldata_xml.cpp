/*
 * This source file is part of the Aleandria project.
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
#include "gmxpre.h"
#include <stdlib.h>
#include <string.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "gromacs/utility/futil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/utility/smalloc.h"
#include "poldata_xml.h"
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
#define NXMLTYPES asize(xmltypes)

enum {
    exmlGENTOP, exmlREFERENCE,
    exmlATOMTYPES, exmlATOMTYPE,
    exmlGT_FORCEFIELD, exmlPOLAR_UNIT, exmlCOMB_RULE, exmlNEXCL,
    exmlFUDGEQQ, exmlFUDGELJ,
    exmlPOLTYPES, exmlPOLTYPE, exmlPTYPE,
    exmlBONDING_RULES, exmlBONDING_RULE, exmlBTYPE,
    exmlELEM, exmlNAME, exmlDESC,
    exmlATYPE, exmlMILLER, exmlVALENCE, exmlBOSQUE,
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
    "bonding_rules", "bonding_rule", "btype",
    "elem", "name", "description",
    "atype", "miller", "valence", "bosque",
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

void PoldataXml::sp(int n, char buf[], int maxindent)
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

void PoldataXml::processAttr(FILE *fp, xmlAttrPtr attr, int elem,
                         int indent, Poldata * pd)
{
  std::string attrname, attrval;
    char  buf[100];
    int   i, kkk;
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
                pd->setPolarUnit(xbuf[exmlPOLAR_UNIT]);
                pd->setPolarRef(xbuf[exmlREFERENCE]);
            }
            break;
        case exmlATOMTYPES:
            if (NN(xbuf[exmlGT_FORCEFIELD]))
            {
                pd->setForceField(xbuf[exmlGT_FORCEFIELD]);
            }
            if (NN(xbuf[exmlFUNCTION]))
            {
                pd->setVdwFunction(xbuf[exmlFUNCTION]);
            }
            if (NN(xbuf[exmlCOMB_RULE]))
            {
                pd->setCombinationRule(xbuf[exmlCOMB_RULE]);
            }
            if (NN(xbuf[exmlNEXCL]))
            {
                pd->setNexcl(atoi(xbuf[exmlNEXCL].c_str()));
            }
            if (NN(xbuf[exmlFUDGEQQ]))
            {
                pd->setFudgeQQ(atoi(xbuf[exmlFUDGEQQ].c_str()));
            }
            if (NN(xbuf[exmlFUDGELJ]))
            {
                pd->setFudgeLJ(atoi(xbuf[exmlFUDGELJ].c_str()));
            }
            break;
        case exmlBSATOMS:
            if (NN(xbuf[exmlPOLAR_UNIT]))
            {
                pd->setBosqueUnit(xbuf[exmlPOLAR_UNIT]);
            }
            break;
        case exmlGT_DIHEDRALS:
            if (NN(xbuf[exmlFUNCTION]))
            {
                pd->setDihedralFunction( egdPDIHS, xbuf[exmlFUNCTION]);
            }
            break;
        case exmlGT_IMPROPERS:
            if (NN(xbuf[exmlFUNCTION]))
            {
                pd->setDihedralFunction(egdIDIHS, xbuf[exmlFUNCTION]);
            }
            break;
        case exmlGT_ANGLES:
            if (NN(xbuf[exmlANGLE_UNIT]))
            {
                pd->setAngleUnit(xbuf[exmlANGLE_UNIT]);
            }
            if (NN(xbuf[exmlFUNCTION]))
            {
                pd->setAngleFunction(xbuf[exmlFUNCTION]);
            }
            break;
        case exmlGT_BONDS:
            if (NN(xbuf[exmlLENGTH_UNIT]))
            {
                pd->setLengthUnit(xbuf[exmlLENGTH_UNIT]);
            }
            if (NN(xbuf[exmlFUNCTION]))
            {
                pd->setBondFunction(xbuf[exmlFUNCTION]);
            }
            break;
        case exmlMILATOMS:
            if (NN(xbuf[exmlTAU_UNIT]) && NN(xbuf[exmlAHP_UNIT]))
            {
                pd->setMillerUnits(xbuf[exmlTAU_UNIT], xbuf[exmlAHP_UNIT]);
            }
            break;
        case exmlPOLTYPE:
            if (NN(xbuf[exmlPTYPE]) && NN(xbuf[exmlMILLER]) && NN(xbuf[exmlBOSQUE]) &&
                NN(xbuf[exmlPOLARIZABILITY]) && NN(xbuf[exmlSIGPOL]))
            {
                pd->addPtype(xbuf[exmlPTYPE],
                                      xbuf[exmlMILLER],
                                      xbuf[exmlBOSQUE],
                                      atof(xbuf[exmlPOLARIZABILITY].c_str()),
                                      atof(xbuf[exmlSIGPOL].c_str()));
            }
            break;
        case exmlATOMTYPE:
            if (NN(xbuf[exmlELEM]) &&
                NN(xbuf[exmlATYPE]) &&
                NN(xbuf[exmlPTYPE]) &&
                NN(xbuf[exmlBTYPE]) &&
                NN(xbuf[exmlVDWPARAMS]) &&
                NN(xbuf[exmlEREF]))
            {
                pd->addAtype(xbuf[exmlELEM],
                                      xbuf[exmlDESC],
                                      xbuf[exmlATYPE],
                                      xbuf[exmlPTYPE],
                                      xbuf[exmlBTYPE],
                                      xbuf[exmlVDWPARAMS],
                                      atof(xbuf[exmlEREF].c_str()));
            }
            break;
        case exmlBONDING_RULE:
            if (NN(xbuf[exmlGEOMETRY]) &&
                NN(xbuf[exmlNUMBONDS]) &&
                NN(xbuf[exmlNEIGHBORS]) &&
                NN(xbuf[exmlATYPE]) &&
                NN(xbuf[exmlAROMATIC]) &&
                NN(xbuf[exmlVALENCE]) &&
                NN(xbuf[exmlNAME]))
            {
                pd->addBondingRule(xbuf[exmlNAME], xbuf[exmlATYPE],
                                             xbuf[exmlGEOMETRY],
                                             atoi(xbuf[exmlNUMBONDS].c_str()),
                                             atof(xbuf[exmlVALENCE].c_str()),
                                             atoi(xbuf[exmlAROMATIC].c_str()),
                                             xbuf[exmlNEIGHBORS]);
            }
            break;
        case exmlMILATOM:
            if (NN(xbuf[exmlMILNAME]) && NN(xbuf[exmlATOMNUMBER]) &&
                NN(xbuf[exmlTAU_AHC]) && NN(xbuf[exmlALPHA_AHP]))
            {
                pd->addMiller(xbuf[exmlMILNAME], atoi(xbuf[exmlATOMNUMBER].c_str()),
                                       atof(xbuf[exmlTAU_AHC].c_str()), atof(xbuf[exmlALPHA_AHP].c_str()));
            }
            break;
        case exmlBSATOM:
            if (NN(xbuf[exmlELEM]) && NN(xbuf[exmlPOLARIZABILITY]))
            {
                pd->addBosque( xbuf[exmlELEM], atof(xbuf[exmlPOLARIZABILITY].c_str()));
            }
            break;
        case exmlGT_BOND:
            if (NN(xbuf[exmlATOM1]) && NN(xbuf[exmlATOM2]) &&
                NN(xbuf[exmlLENGTH]) && NN(xbuf[exmlSIGMA]) && NN(xbuf[exmlBONDORDER]) &&
                NN(xbuf[exmlPARAMS]) && NN(xbuf[exmlNTRAIN]))
            {
                pd->addBond( xbuf[exmlATOM1], xbuf[exmlATOM2],
                                     atof(xbuf[exmlLENGTH].c_str()),
                                     atof(xbuf[exmlSIGMA].c_str()),
                                     atoi(xbuf[exmlNTRAIN].c_str()),
                                     atof(xbuf[exmlBONDORDER].c_str()),
                                     xbuf[exmlPARAMS]);
            }
            break;
        case exmlGT_ANGLE:
            if (NN(xbuf[exmlATOM1]) && NN(xbuf[exmlATOM2]) &&
                NN(xbuf[exmlATOM3]) && NN(xbuf[exmlANGLE]) && NN(xbuf[exmlSIGMA]) &&
                NN(xbuf[exmlPARAMS]) && NN(xbuf[exmlNTRAIN]))
            {
                pd->addAngle( xbuf[exmlATOM1], xbuf[exmlATOM2],
                                      xbuf[exmlATOM3], atof(xbuf[exmlANGLE].c_str()),
                                      atof(xbuf[exmlSIGMA].c_str()), atoi(xbuf[exmlNTRAIN].c_str()),
                                      xbuf[exmlPARAMS].c_str());
            }
            break;
        case exmlGT_DIHEDRAL:
            if (NN(xbuf[exmlATOM1]) && NN(xbuf[exmlATOM2]) &&
                NN(xbuf[exmlATOM3]) && NN(xbuf[exmlATOM4]) &&
                NN(xbuf[exmlANGLE]) && NN(xbuf[exmlSIGMA]) &&
                NN(xbuf[exmlPARAMS]) && NN(xbuf[exmlNTRAIN]))
            {
                pd->addDihedral( egdPDIHS,
                                         xbuf[exmlATOM1], xbuf[exmlATOM2],
                                         xbuf[exmlATOM3], xbuf[exmlATOM4],
                                         atof(xbuf[exmlANGLE].c_str()), atof(xbuf[exmlSIGMA].c_str()),
                                         atoi(xbuf[exmlNTRAIN].c_str()),
                                         xbuf[exmlPARAMS]);
            }
            break;
        case exmlGT_IMPROPER:
            if (NN(xbuf[exmlATOM1]) && NN(xbuf[exmlATOM2]) &&
                NN(xbuf[exmlATOM3]) && NN(xbuf[exmlATOM4]) &&
                NN(xbuf[exmlANGLE]) && NN(xbuf[exmlSIGMA]) &&
                NN(xbuf[exmlPARAMS]) && NN(xbuf[exmlNTRAIN]))
            {
                pd->addDihedral( egdIDIHS,
                                         xbuf[exmlATOM1], xbuf[exmlATOM2],
                                         xbuf[exmlATOM3], xbuf[exmlATOM4],
                                         atof(xbuf[exmlANGLE].c_str()), atof(xbuf[exmlSIGMA].c_str()),
				 atoi(xbuf[exmlNTRAIN].c_str()),
                                         xbuf[exmlPARAMS]);
            }
            break;
        case exmlSYM_CHARGE:
            if (NN(xbuf[exmlCENTRAL]) && NN(xbuf[exmlATTACHED]) &&
                NN(xbuf[exmlNUMATTACH]))
            {
                pd->addSymcharges( xbuf[exmlCENTRAL],
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
	      pd->setEemprops(pd->name2eemtype(xbuf[exmlMODEL]), xbuf[exmlNAME].c_str(),
                                         atof(xbuf[exmlJ0].c_str()), atof(xbuf[exmlCHI0].c_str()),
                                         xbuf[exmlZETA], xbuf[exmlCHARGES], xbuf[exmlROW]);
            }
            break;
        case exmlEEMPROP_REF:
            if (NN(xbuf[exmlMODEL]) && NN(xbuf[exmlEPREF]))
            {
                pd->setEpref(pd->name2eemtype(xbuf[exmlMODEL]), xbuf[exmlEPREF]);
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

void PoldataXml::processTree(FILE *fp, xmlNodePtr tree, int indent,
                         Poldata * pd, gmx_atomprop_t aps)
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

  Poldata * PoldataXml::read(const char *fn, gmx_atomprop_t aps)
{
    xmlDocPtr     doc;
    Poldata * pd;
    std::string         fn2;

    if (NULL != fn)
    {
      fn2 = gmxlibfn(fn);
    }
    if (0 == fn2.size())
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
        fprintf(stderr, "\nError reading XML file %s. Run a syntax checker such as nsgmls.\n",
                fn2.c_str());
        //sfree(fn2);
        return NULL;
    }

    pd = new Poldata();
    pd->setFilename( fn2);
    processTree(debug, doc->children, 0, pd, aps);

    xmlFreeDoc(doc);

    if (NULL != debug)
    {
        write("pdout.dat", pd, 0);
    }
    //sfree(fn2);

    return pd;
}

void PoldataXml::addXmlPoldata(xmlNodePtr parent, Poldata * pd)
{
    xmlNodePtr              child, grandchild;
    int                     i, atomnumber, numbonds, nexcl,
                            numattach, bAromatic, ntrain;
    ChargeDistributionModel model;
    std::string                   elem, geometry, name, atype, vdwparams, btype,
    atom1, atom2, atom3, atom4, central, attached, tau_unit, ahp_unit,
    epref, desc, params;
    std::string  neighbors, zeta, qstr, rowstr;
   double polarizability, sig_pol, length, tau_ahc, alpha_ahp, angle, J0, chi0,
        bondorder, sigma, fudgeQQ, fudgeLJ, valence, ref_enthalpy;
    std::string tmp, func,blu;
    child = add_xml_child(parent, exml_names[exmlATOMTYPES]);
    tmp   = pd->getForceField();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names[exmlGT_FORCEFIELD], tmp.c_str());
    }
    tmp = pd->getVdwFunction();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names[exmlFUNCTION], tmp.c_str());
    }
    tmp = pd->getCombinationRule();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names[exmlCOMB_RULE], tmp.c_str());
    }
    nexcl = pd->getNexcl();
    add_xml_int(child, exml_names[exmlNEXCL], nexcl);
    fudgeQQ = pd->getFudgeQQ();
    add_xml_double(child, exml_names[exmlFUDGEQQ], fudgeQQ);
    fudgeLJ = pd->getFudgeLJ();
    add_xml_double(child, exml_names[exmlFUDGELJ], fudgeLJ);
    {
      std::string ptype;
        while (1 == pd->getAtype( &elem, &desc, &atype, &ptype, &btype,
                                          &vdwparams, &ref_enthalpy))
        {
            grandchild = add_xml_child(child, exml_names[exmlATOMTYPE]);
            add_xml_char(grandchild, exml_names[exmlELEM], elem.c_str());
            add_xml_char(grandchild, exml_names[exmlDESC], desc.c_str());
            add_xml_char(grandchild, exml_names[exmlATYPE], atype.c_str());
            add_xml_char(grandchild, exml_names[exmlPTYPE], ptype.c_str());
            add_xml_char(grandchild, exml_names[exmlBTYPE], btype.c_str());
            add_xml_char(grandchild, exml_names[exmlVDWPARAMS], vdwparams.c_str());
            add_xml_double(grandchild, exml_names[exmlEREF], ref_enthalpy);
	    /* sfree(elem);
            sfree(desc);
            sfree(atype);
            sfree(ptype);
            sfree(btype);
            sfree(vdwparams);*/
        }
    }
    child = add_xml_child(parent, exml_names[exmlPOLTYPES]);
    tmp   = pd->getPolarUnit();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names[exmlPOLAR_UNIT], tmp.c_str());
    }
    tmp = pd->getPolarRef();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names[exmlREFERENCE], tmp.c_str());
    }
    {
      std::string miller, bosque, ptype;
        while (1 == pd->getPtype( &ptype, &miller, &bosque, &polarizability, &sig_pol))
        {
            grandchild = add_xml_child(child, exml_names[exmlPOLTYPE]);
            add_xml_char(grandchild, exml_names[exmlPTYPE], ptype.c_str());
            add_xml_char(grandchild, exml_names[exmlMILLER], miller.c_str());
            add_xml_char(grandchild, exml_names[exmlBOSQUE], bosque.c_str());
            add_xml_double(grandchild, exml_names[exmlPOLARIZABILITY], polarizability);
            add_xml_double(grandchild, exml_names[exmlSIGPOL], sig_pol);
            /*sfree(ptype);
            sfree(miller);
            sfree(bosque);*/
        }
    }

    child = add_xml_child(parent, exml_names[exmlBONDING_RULES]);
    while (1 == pd->getBondingRule( &name, &atype, &geometry,
                                             &numbonds, &valence, &bAromatic, &neighbors))
    {
        grandchild = add_xml_child(child, exml_names[exmlBONDING_RULE]);
        add_xml_char(grandchild, exml_names[exmlNAME], name.c_str());
        add_xml_char(grandchild, exml_names[exmlATYPE], atype.c_str());
        add_xml_char(grandchild, exml_names[exmlGEOMETRY], geometry.c_str());
        add_xml_int(grandchild, exml_names[exmlNUMBONDS], numbonds);
        add_xml_double(grandchild, exml_names[exmlVALENCE], valence);
        add_xml_int(grandchild, exml_names[exmlAROMATIC], bAromatic);
        add_xml_char(grandchild, exml_names[exmlNEIGHBORS], neighbors.c_str());
        /*sfree(name);
        sfree(atype);
        sfree(geometry);
        sfree(neighbors);*/
    }

    child = add_xml_child(parent, exml_names[exmlGT_BONDS]);
    blu = pd->getLengthUnit();
    if (blu.size() != 0)
      {
        add_xml_char(child, exml_names[exmlLENGTH_UNIT], blu.c_str());
      }
    func = pd->getBondFunction();
    if (func.size() != 0)
      {
        add_xml_char(child, exml_names[exmlFUNCTION], func.c_str());
      }
    while (pd->getBond( &atom1, &atom2, &length, &sigma,
			&ntrain, &bondorder, &params) > 0)
    {
        grandchild = add_xml_child(child, exml_names[exmlGT_BOND]);
        add_xml_char(grandchild, exml_names[exmlATOM1], atom1.c_str());
        add_xml_char(grandchild, exml_names[exmlATOM2], atom2.c_str());
        add_xml_double(grandchild, exml_names[exmlLENGTH], length);
        add_xml_double(grandchild, exml_names[exmlSIGMA], sigma);
        add_xml_int(grandchild, exml_names[exmlNTRAIN], ntrain);
        add_xml_double(grandchild, exml_names[exmlBONDORDER], bondorder);
        add_xml_char(grandchild, exml_names[exmlPARAMS], params.c_str());
        /*sfree(atom1);
        sfree(atom2);
        sfree(params);*/
    }

    child = add_xml_child(parent, exml_names[exmlGT_ANGLES]);
    blu = pd->getAngleUnit();
    if (blu.size() != 0)
    {
        add_xml_char(child, exml_names[exmlANGLE_UNIT], blu.c_str());
    }
    func = pd->getAngleFunction();
    if (func.size() != 0)
    {
        add_xml_char(child, exml_names[exmlFUNCTION], func.c_str());
    }
    while (pd->getAngle( &atom1, &atom2, &atom3, &angle, &sigma,
                                 &ntrain, &params) > 0)
    {
        grandchild = add_xml_child(child, exml_names[exmlGT_ANGLE]);
        add_xml_char(grandchild, exml_names[exmlATOM1], atom1.c_str());
        add_xml_char(grandchild, exml_names[exmlATOM2], atom2.c_str());
        add_xml_char(grandchild, exml_names[exmlATOM3], atom3.c_str());
        add_xml_double(grandchild, exml_names[exmlANGLE], angle);
        add_xml_double(grandchild, exml_names[exmlSIGMA], sigma);
        add_xml_int(grandchild, exml_names[exmlNTRAIN], ntrain);
        add_xml_char(grandchild, exml_names[exmlPARAMS], params.c_str());
	/* sfree(atom1);
        sfree(atom2);
        sfree(atom3);
        sfree(params);*/
    }

    for (i = 0; (i < egdNR); i++)
    {
        int exs[egdNR] = { exmlGT_DIHEDRALS, exmlGT_IMPROPERS };
        int ex[egdNR]  = { exmlGT_DIHEDRAL, exmlGT_IMPROPER };

        child = add_xml_child(parent, exml_names[exs[i]]);
	blu = pd->getAngleUnit();
        if (blu.size() != 0)
        {
            add_xml_char(child, exml_names[exmlANGLE_UNIT], blu.c_str());
        }
	func = pd->getDihedralFunction( i);
        if (func.size() != 0)
        {
            add_xml_char(child, exml_names[exmlFUNCTION], func.c_str());
        }
        while (pd->getDihedral( i, &atom1, &atom2, &atom3, &atom4,
                                        &angle, &sigma, &ntrain, &params) > 0)
        {
            grandchild = add_xml_child(child, exml_names[ex[i]]);
            add_xml_char(grandchild, exml_names[exmlATOM1], atom1.c_str());
            add_xml_char(grandchild, exml_names[exmlATOM2], atom2.c_str());
            add_xml_char(grandchild, exml_names[exmlATOM3], atom3.c_str());
            add_xml_char(grandchild, exml_names[exmlATOM4], atom4.c_str());
            add_xml_double(grandchild, exml_names[exmlANGLE], angle);
            add_xml_double(grandchild, exml_names[exmlSIGMA], sigma);
            add_xml_int(grandchild, exml_names[exmlNTRAIN], ntrain);
            add_xml_char(grandchild, exml_names[exmlPARAMS], params.c_str());
            /*sfree(atom1);
            sfree(atom2);
            sfree(atom3);
            sfree(atom4);
            sfree(params);*/
        }
    }
    child = add_xml_child(parent, exml_names[exmlBSATOMS]);
    tmp = pd->getBosqueUnit();
    if (tmp.size() != 0)
    {
        add_xml_char(child, exml_names[exmlPOLAR_UNIT], tmp.c_str());
    }

    while (1 == pd->getBosque( &name, &polarizability))
    {
        grandchild = add_xml_child(child, exml_names[exmlBSATOM]);
        add_xml_char(grandchild, exml_names[exmlELEM], name.c_str());
        add_xml_double(grandchild, exml_names[exmlPOLARIZABILITY], polarizability);
    }
    child = add_xml_child(parent, exml_names[exmlMILATOMS]);
    pd->getMillerUnits( &tau_unit, &ahp_unit);
    if (tau_unit.size() != 0)
    {
        add_xml_char(child, exml_names[exmlTAU_UNIT], tau_unit.c_str());
        //sfree(tau_unit);
    }
    if (ahp_unit.size() != 0)
    {
        add_xml_char(child, exml_names[exmlAHP_UNIT], ahp_unit.c_str());
        //sfree(ahp_unit);
    }

    while (1 == pd->getMiller( &name, &atomnumber, &tau_ahc, &alpha_ahp))
    {
        grandchild = add_xml_child(child, exml_names[exmlMILATOM]);
        add_xml_char(grandchild, exml_names[exmlMILNAME], name.c_str());
        add_xml_int(grandchild, exml_names[exmlATOMNUMBER], atomnumber);
        add_xml_double(grandchild, exml_names[exmlTAU_AHC], tau_ahc);
        add_xml_double(grandchild, exml_names[exmlALPHA_AHP], alpha_ahp);
        //sfree(name);
    }

    child = add_xml_child(parent, exml_names[exmlSYMMETRIC_CHARGES]);

    while (pd->getSymcharges( &central, &attached, &numattach) == 1)
    {
        grandchild = add_xml_child(child, exml_names[exmlSYM_CHARGE]);
        add_xml_char(grandchild, exml_names[exmlCENTRAL], central.c_str());
        add_xml_char(grandchild, exml_names[exmlATTACHED], attached.c_str());
        add_xml_int(grandchild, exml_names[exmlNUMATTACH], numattach);
        /*sfree(central);
	  sfree(attached);*******/
    }

    child = add_xml_child(parent, exml_names[exmlEEMPROPS]);
    while (pd->getEemprops( &model, &name, &J0, &chi0, &zeta, &qstr, &rowstr) == 1)
    {
        grandchild = add_xml_child(child, exml_names[exmlEEMPROP]);
        add_xml_char(grandchild, exml_names[exmlMODEL],
                     pd->getEemtypeName(model).c_str());
        add_xml_char(grandchild, exml_names[exmlNAME], name.c_str());
        add_xml_double(grandchild, exml_names[exmlJ0], J0);
        add_xml_double(grandchild, exml_names[exmlCHI0], chi0);
        add_xml_char(grandchild, exml_names[exmlZETA], zeta.c_str());
        add_xml_char(grandchild, exml_names[exmlCHARGES], qstr.c_str());
        add_xml_char(grandchild, exml_names[exmlROW], rowstr.c_str());
	/* sfree(zeta);
        sfree(qstr);
        sfree(rowstr);*/
    }
    while (pd->listEpref( &model, &epref) == 1)
    {
        grandchild = add_xml_child(child, exml_names[exmlEEMPROP_REF]);
        add_xml_char(grandchild, exml_names[exmlMODEL], pd->getEemtypeName(model).c_str());
        add_xml_char(grandchild, exml_names[exmlEPREF], epref.c_str());
        //sfree(epref);
    }
}

  void PoldataXml::write(const std::string fn, Poldata * pd,
                       gmx_bool bCompress)
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

    xmlSetDocCompressMode(doc, (int)bCompress);
    xmlIndentTreeOutput = 1;
    if (xmlSaveFormatFileEnc(fn.c_str(), doc, "ISO-8859-1", 2) == 0)
    {
      gmx_fatal(FARGS, "Saving file", fn.c_str());
    }
    xmlFreeDoc(doc);
}
}
