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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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

#define NN(x) (NULL != (x))

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
    exmlGEOMETRY, exmlNUMBONDS, exmlPOLARIZABILITY, exmlSIGPOL, exmlVDWPARAMS,
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

static const char *exml_names[exmlNR] = {
    "gentop", "reference",
    "atomtypes", "atomtype", "forcefield", "polarizability_unit", "combination_rule", "nexclusions",
    "fudgeQQ", "fudgeLJ",
    "poltypes", "poltype", "ptype",
    "bonding_rules", "bonding_rule", "btype",
    "elem", "name", "description",
    "atype", "miller", "valence", "bosque",
    "neighbors", "aromatic",
    "geometry", "numbonds", "polarizability", "sigma_pol", "vdwparams",
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

static void process_attr(FILE *fp, xmlAttrPtr attr, int elem,
                         int indent, gmx_poldata_t pd)
{
    char *attrname, *attrval;
    char  buf[100];
    int   i, kkk;
    char *xbuf[exmlNR];

    for (i = 0; (i < exmlNR); i++)
    {
        xbuf[i] = NULL;
    }
    while (attr != NULL)
    {
        attrname = (char *)attr->name;
        attrval  = (char *)attr->children->content;

#define atest(s) ((strcasecmp(attrname, s) == 0) && (attrval != NULL))
        kkk = find_elem(attrname, exmlNR, exml_names);
        if (-1 != kkk)
        {
            if (attrval != NULL)
            {
                xbuf[kkk] = strdup(attrval);
            }

            if (NULL != fp)
            {
                sp(indent, buf, 99);
                fprintf(fp, "%sProperty: '%s' Value: '%s'\n", buf, attrname, attrval);
            }
        }
        else
        {
            fprintf(stderr, "Ignoring invalid attribute %s\n", attrname);
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
                gmx_poldata_set_polar_unit(pd, xbuf[exmlPOLAR_UNIT]);
                gmx_poldata_set_polar_ref(pd, xbuf[exmlREFERENCE]);
            }
            break;
        case exmlATOMTYPES:
            if (NN(xbuf[exmlGT_FORCEFIELD]))
            {
                gmx_poldata_set_force_field(pd, xbuf[exmlGT_FORCEFIELD]);
            }
            if (NN(xbuf[exmlFUNCTION]))
            {
                gmx_poldata_set_vdw_function(pd, xbuf[exmlFUNCTION]);
            }
            if (NN(xbuf[exmlCOMB_RULE]))
            {
                gmx_poldata_set_combination_rule(pd, xbuf[exmlCOMB_RULE]);
            }
            if (NN(xbuf[exmlNEXCL]))
            {
                gmx_poldata_set_nexcl(pd, atoi(xbuf[exmlNEXCL]));
            }
            if (NN(xbuf[exmlFUDGEQQ]))
            {
                gmx_poldata_set_fudgeQQ(pd, atoi(xbuf[exmlFUDGEQQ]));
            }
            if (NN(xbuf[exmlFUDGELJ]))
            {
                gmx_poldata_set_fudgeLJ(pd, atoi(xbuf[exmlFUDGELJ]));
            }
            break;
        case exmlBSATOMS:
            if (NN(xbuf[exmlPOLAR_UNIT]))
            {
                gmx_poldata_set_bosque_unit(pd, xbuf[exmlPOLAR_UNIT]);
            }
            break;
        case exmlGT_DIHEDRALS:
            if (NN(xbuf[exmlFUNCTION]))
            {
                gmx_poldata_set_dihedral_function(pd, egdPDIHS, xbuf[exmlFUNCTION]);
            }
            break;
        case exmlGT_IMPROPERS:
            if (NN(xbuf[exmlFUNCTION]))
            {
                gmx_poldata_set_dihedral_function(pd, egdIDIHS, xbuf[exmlFUNCTION]);
            }
            break;
        case exmlGT_ANGLES:
            if (NN(xbuf[exmlANGLE_UNIT]))
            {
                gmx_poldata_set_angle_unit(pd, xbuf[exmlANGLE_UNIT]);
            }
            if (NN(xbuf[exmlFUNCTION]))
            {
                gmx_poldata_set_angle_function(pd, xbuf[exmlFUNCTION]);
            }
            break;
        case exmlGT_BONDS:
            if (NN(xbuf[exmlLENGTH_UNIT]))
            {
                gmx_poldata_set_length_unit(pd, xbuf[exmlLENGTH_UNIT]);
            }
            if (NN(xbuf[exmlFUNCTION]))
            {
                gmx_poldata_set_bond_function(pd, xbuf[exmlFUNCTION]);
            }
            break;
        case exmlMILATOMS:
            if (NN(xbuf[exmlTAU_UNIT]) && NN(xbuf[exmlAHP_UNIT]))
            {
                gmx_poldata_set_miller_units(pd, xbuf[exmlTAU_UNIT], xbuf[exmlAHP_UNIT]);
            }
            break;
        case exmlPOLTYPE:
            if (NN(xbuf[exmlPTYPE]) && NN(xbuf[exmlMILLER]) && NN(xbuf[exmlBOSQUE]) &&
                NN(xbuf[exmlPOLARIZABILITY]) && NN(xbuf[exmlSIGPOL]))
            {
                gmx_poldata_add_ptype(pd,
                                      xbuf[exmlPTYPE],
                                      xbuf[exmlMILLER],
                                      xbuf[exmlBOSQUE],
                                      atof(xbuf[exmlPOLARIZABILITY]),
                                      atof(xbuf[exmlSIGPOL]));
            }
            break;
        case exmlATOMTYPE:
            if (NN(xbuf[exmlELEM]) &&
                NN(xbuf[exmlATYPE]) &&
                NN(xbuf[exmlPTYPE]) &&
                NN(xbuf[exmlBTYPE]) &&
                NN(xbuf[exmlVDWPARAMS]))
            {
                gmx_poldata_add_atype(pd, xbuf[exmlELEM],
                                      xbuf[exmlDESC] ? xbuf[exmlDESC] : (char *) "",
                                      xbuf[exmlATYPE],
                                      xbuf[exmlPTYPE],
                                      xbuf[exmlBTYPE],
                                      xbuf[exmlVDWPARAMS]);
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
                gmx_poldata_add_bonding_rule(pd, xbuf[exmlNAME], xbuf[exmlATYPE],
                                             xbuf[exmlGEOMETRY],
                                             atoi(xbuf[exmlNUMBONDS]),
                                             atof(xbuf[exmlVALENCE]),
                                             atoi(xbuf[exmlAROMATIC]),
                                             xbuf[exmlNEIGHBORS]);
            }
            break;
        case exmlMILATOM:
            if (NN(xbuf[exmlMILNAME]) && NN(xbuf[exmlATOMNUMBER]) &&
                NN(xbuf[exmlTAU_AHC]) && NN(xbuf[exmlALPHA_AHP]))
            {
                gmx_poldata_add_miller(pd, xbuf[exmlMILNAME], atoi(xbuf[exmlATOMNUMBER]),
                                       atof(xbuf[exmlTAU_AHC]), atof(xbuf[exmlALPHA_AHP]));
            }
            break;
        case exmlBSATOM:
            if (NN(xbuf[exmlELEM]) && NN(xbuf[exmlPOLARIZABILITY]))
            {
                gmx_poldata_add_bosque(pd, xbuf[exmlELEM], atof(xbuf[exmlPOLARIZABILITY]));
            }
            break;
        case exmlGT_BOND:
            if (NN(xbuf[exmlATOM1]) && NN(xbuf[exmlATOM2]) &&
                NN(xbuf[exmlLENGTH]) && NN(xbuf[exmlSIGMA]) && NN(xbuf[exmlBONDORDER]) &&
                NN(xbuf[exmlPARAMS]) && NN(xbuf[exmlNTRAIN]))
            {
                gmx_poldata_add_bond(pd, xbuf[exmlATOM1], xbuf[exmlATOM2],
                                     atof(xbuf[exmlLENGTH]),
                                     atof(xbuf[exmlSIGMA]),
                                     atoi(xbuf[exmlNTRAIN]),
                                     atof(xbuf[exmlBONDORDER]),
                                     xbuf[exmlPARAMS]);
            }
            break;
        case exmlGT_ANGLE:
            if (NN(xbuf[exmlATOM1]) && NN(xbuf[exmlATOM2]) &&
                NN(xbuf[exmlATOM3]) && NN(xbuf[exmlANGLE]) && NN(xbuf[exmlSIGMA]) &&
                NN(xbuf[exmlPARAMS]) && NN(xbuf[exmlNTRAIN]))
            {
                gmx_poldata_add_angle(pd, xbuf[exmlATOM1], xbuf[exmlATOM2],
                                      xbuf[exmlATOM3], atof(xbuf[exmlANGLE]),
                                      atof(xbuf[exmlSIGMA]), atoi(xbuf[exmlNTRAIN]),
                                      xbuf[exmlPARAMS]);
            }
            break;
        case exmlGT_DIHEDRAL:
            if (NN(xbuf[exmlATOM1]) && NN(xbuf[exmlATOM2]) &&
                NN(xbuf[exmlATOM3]) && NN(xbuf[exmlATOM4]) &&
                NN(xbuf[exmlANGLE]) && NN(xbuf[exmlSIGMA]) &&
                NN(xbuf[exmlPARAMS]) && NN(xbuf[exmlNTRAIN]))
            {
                gmx_poldata_add_dihedral(pd, egdPDIHS,
                                         xbuf[exmlATOM1], xbuf[exmlATOM2],
                                         xbuf[exmlATOM3], xbuf[exmlATOM4],
                                         atof(xbuf[exmlANGLE]), atof(xbuf[exmlSIGMA]),
                                         atoi(xbuf[exmlNTRAIN]),
                                         xbuf[exmlPARAMS]);
            }
            break;
        case exmlGT_IMPROPER:
            if (NN(xbuf[exmlATOM1]) && NN(xbuf[exmlATOM2]) &&
                NN(xbuf[exmlATOM3]) && NN(xbuf[exmlATOM4]) &&
                NN(xbuf[exmlANGLE]) && NN(xbuf[exmlSIGMA]) &&
                NN(xbuf[exmlPARAMS]) && NN(xbuf[exmlNTRAIN]))
            {
                gmx_poldata_add_dihedral(pd, egdIDIHS,
                                         xbuf[exmlATOM1], xbuf[exmlATOM2],
                                         xbuf[exmlATOM3], xbuf[exmlATOM4],
                                         atof(xbuf[exmlANGLE]), atof(xbuf[exmlSIGMA]),
                                         atoi(xbuf[exmlNTRAIN]),
                                         xbuf[exmlPARAMS]);
            }
            break;
        case exmlSYM_CHARGE:
            if (NN(xbuf[exmlCENTRAL]) && NN(xbuf[exmlATTACHED]) &&
                NN(xbuf[exmlNUMATTACH]))
            {
                gmx_poldata_add_symcharges(pd, xbuf[exmlCENTRAL],
                                           xbuf[exmlATTACHED],
                                           atoi(xbuf[exmlNUMATTACH]));
            }
            break;
        case exmlEEMPROP:
            if (NN(xbuf[exmlMODEL]) && NN(xbuf[exmlNAME]) &&
                NN(xbuf[exmlCHI0])  && NN(xbuf[exmlJ0]) &&
                NN(xbuf[exmlZETA])  && NN(xbuf[exmlCHARGES]) &&
                NN(xbuf[exmlROW]))
            {
                gmx_poldata_set_eemprops(pd, name2eemtype(xbuf[exmlMODEL]), xbuf[exmlNAME],
                                         atof(xbuf[exmlJ0]), atof(xbuf[exmlCHI0]),
                                         xbuf[exmlZETA], xbuf[exmlCHARGES], xbuf[exmlROW]);
            }
            break;
        case exmlEEMPROP_REF:
            if (NN(xbuf[exmlMODEL]) && NN(xbuf[exmlEPREF]))
            {
                gmx_poldata_set_epref(pd, name2eemtype(xbuf[exmlMODEL]), xbuf[exmlEPREF]);
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

}

static void process_tree(FILE *fp, xmlNodePtr tree, int indent,
                         gmx_poldata_t pd, gmx_atomprop_t aps)
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
                    process_attr(fp, tree->properties, elem, indent+2, pd);
                }

                if (tree->children)
                {
                    process_tree(fp, tree->children, indent+2, pd, aps);
                }
            }
        }
        tree = tree->next;
    }
}

gmx_poldata_t gmx_poldata_read(const char *fn, gmx_atomprop_t aps)
{
    xmlDocPtr     doc;
    gmx_poldata_t pd;
    char         *fn2 = NULL;

    if (NULL != fn)
    {
        fn2 = (char *)gmxlibfn(fn);
    }
    if (NULL == fn2)
    {
        fn2 = (char *)gmxlibfn("alexandria.ff/gentop.dat");
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
        return NULL;
    }
    pd = gmx_poldata_init();
    gmx_poldata_set_filename(pd, fn2);
    process_tree(debug, doc->children, 0, pd, aps);

    xmlFreeDoc(doc);

    if (NULL != debug)
    {
        gmx_poldata_write("pdout.dat", pd, 0);
    }
    sfree(fn2);

    return pd;
}

static void add_xml_poldata(xmlNodePtr parent, gmx_poldata_t pd)
{
    xmlNodePtr            child, grandchild;
    int                   i, atomnumber, numbonds, nexcl,
                          numattach, bAromatic, ntrain;
    ChargeGenerationModel model;
    char                 *elem, *geometry, *name, *atype, *vdwparams, *blu, *btype,
    *atom1, *atom2, *atom3, *atom4, *tmp, *central, *attached, *tau_unit, *ahp_unit,
    *epref, *desc, *params, *func;
    char  *neighbors, *zeta, *qstr, *rowstr;
    double polarizability, sig_pol, length, tau_ahc, alpha_ahp, angle, J0, chi0,
           bondorder, sigma, fudgeQQ, fudgeLJ, valence;

    child = add_xml_child(parent, exml_names[exmlATOMTYPES]);
    tmp   = gmx_poldata_get_force_field(pd);
    if (NULL != tmp)
    {
        add_xml_char(child, exml_names[exmlGT_FORCEFIELD], tmp);
    }
    tmp = gmx_poldata_get_vdw_function(pd);
    if (NULL != tmp)
    {
        add_xml_char(child, exml_names[exmlFUNCTION], tmp);
    }
    tmp = gmx_poldata_get_combination_rule(pd);
    if (NULL != tmp)
    {
        add_xml_char(child, exml_names[exmlCOMB_RULE], tmp);
    }
    nexcl = gmx_poldata_get_nexcl(pd);
    add_xml_int(child, exml_names[exmlNEXCL], nexcl);
    fudgeQQ = gmx_poldata_get_fudgeQQ(pd);
    add_xml_double(child, exml_names[exmlFUDGEQQ], fudgeQQ);
    fudgeLJ = gmx_poldata_get_fudgeLJ(pd);
    add_xml_double(child, exml_names[exmlFUDGELJ], fudgeLJ);
    {
        char *ptype;
        while (1 == gmx_poldata_get_atype(pd, &elem, &desc, &atype, &ptype, &btype,
                                          &vdwparams))
        {
            grandchild = add_xml_child(child, exml_names[exmlATOMTYPE]);
            add_xml_char(grandchild, exml_names[exmlELEM], elem);
            add_xml_char(grandchild, exml_names[exmlDESC], desc);
            add_xml_char(grandchild, exml_names[exmlATYPE], atype);
            add_xml_char(grandchild, exml_names[exmlPTYPE], ptype);
            add_xml_char(grandchild, exml_names[exmlBTYPE], btype);
            add_xml_char(grandchild, exml_names[exmlVDWPARAMS], vdwparams);
            sfree(elem);
            sfree(desc);
            sfree(atype);
            sfree(ptype);
            sfree(btype);
            sfree(vdwparams);
        }
    }
    child = add_xml_child(parent, exml_names[exmlPOLTYPES]);
    tmp   = gmx_poldata_get_polar_unit(pd);
    if (NULL != tmp)
    {
        add_xml_char(child, exml_names[exmlPOLAR_UNIT], tmp);
    }
    tmp = gmx_poldata_get_polar_ref(pd);
    if (NULL != tmp)
    {
        add_xml_char(child, exml_names[exmlREFERENCE], tmp);
    }
    {
        char *miller, *bosque, *ptype;
        while (1 == gmx_poldata_get_ptype(pd, &ptype, &miller, &bosque, &polarizability, &sig_pol))
        {
            grandchild = add_xml_child(child, exml_names[exmlPOLTYPE]);
            add_xml_char(grandchild, exml_names[exmlPTYPE], ptype);
            add_xml_char(grandchild, exml_names[exmlMILLER], miller);
            add_xml_char(grandchild, exml_names[exmlBOSQUE], bosque);
            add_xml_double(grandchild, exml_names[exmlPOLARIZABILITY], polarizability);
            add_xml_double(grandchild, exml_names[exmlSIGPOL], sig_pol);
            sfree(ptype);
            sfree(miller);
            sfree(bosque);
        }
    }

    child = add_xml_child(parent, exml_names[exmlBONDING_RULES]);
    while (1 == gmx_poldata_get_bonding_rule(pd, &name, &atype, &geometry,
                                             &numbonds, &valence, &bAromatic, &neighbors))
    {
        grandchild = add_xml_child(child, exml_names[exmlBONDING_RULE]);
        add_xml_char(grandchild, exml_names[exmlNAME], name);
        add_xml_char(grandchild, exml_names[exmlATYPE], atype);
        add_xml_char(grandchild, exml_names[exmlGEOMETRY], geometry);
        add_xml_int(grandchild, exml_names[exmlNUMBONDS], numbonds);
        add_xml_double(grandchild, exml_names[exmlVALENCE], valence);
        add_xml_int(grandchild, exml_names[exmlAROMATIC], bAromatic);
        add_xml_char(grandchild, exml_names[exmlNEIGHBORS], neighbors);
        sfree(name);
        sfree(atype);
        sfree(geometry);
        sfree(neighbors);
    }

    child = add_xml_child(parent, exml_names[exmlGT_BONDS]);
    if ((blu = gmx_poldata_get_length_unit(pd)) != NULL)
    {
        add_xml_char(child, exml_names[exmlLENGTH_UNIT], blu);
    }
    if ((func = gmx_poldata_get_bond_function(pd)) != NULL)
    {
        add_xml_char(child, exml_names[exmlFUNCTION], func);
    }
    while (gmx_poldata_get_bond(pd, &atom1, &atom2, &length, &sigma,
                                &ntrain, &bondorder, &params) > 0)
    {
        grandchild = add_xml_child(child, exml_names[exmlGT_BOND]);
        add_xml_char(grandchild, exml_names[exmlATOM1], atom1);
        add_xml_char(grandchild, exml_names[exmlATOM2], atom2);
        add_xml_double(grandchild, exml_names[exmlLENGTH], length);
        add_xml_double(grandchild, exml_names[exmlSIGMA], sigma);
        add_xml_int(grandchild, exml_names[exmlNTRAIN], ntrain);
        add_xml_double(grandchild, exml_names[exmlBONDORDER], bondorder);
        add_xml_char(grandchild, exml_names[exmlPARAMS], params);
        sfree(atom1);
        sfree(atom2);
        sfree(params);
    }

    child = add_xml_child(parent, exml_names[exmlGT_ANGLES]);
    if ((blu = gmx_poldata_get_angle_unit(pd)) != NULL)
    {
        add_xml_char(child, exml_names[exmlANGLE_UNIT], blu);
    }
    if ((func = gmx_poldata_get_angle_function(pd)) != NULL)
    {
        add_xml_char(child, exml_names[exmlFUNCTION], func);
    }
    while (gmx_poldata_get_angle(pd, &atom1, &atom2, &atom3, &angle, &sigma,
                                 &ntrain, &params) > 0)
    {
        grandchild = add_xml_child(child, exml_names[exmlGT_ANGLE]);
        add_xml_char(grandchild, exml_names[exmlATOM1], atom1);
        add_xml_char(grandchild, exml_names[exmlATOM2], atom2);
        add_xml_char(grandchild, exml_names[exmlATOM3], atom3);
        add_xml_double(grandchild, exml_names[exmlANGLE], angle);
        add_xml_double(grandchild, exml_names[exmlSIGMA], sigma);
        add_xml_int(grandchild, exml_names[exmlNTRAIN], ntrain);
        add_xml_char(grandchild, exml_names[exmlPARAMS], params);
        sfree(atom1);
        sfree(atom2);
        sfree(atom3);
        sfree(params);
    }

    for (i = 0; (i < egdNR); i++)
    {
        int exs[egdNR] = { exmlGT_DIHEDRALS, exmlGT_IMPROPERS };
        int ex[egdNR]  = { exmlGT_DIHEDRAL, exmlGT_IMPROPER };

        child = add_xml_child(parent, exml_names[exs[i]]);
        if ((blu = gmx_poldata_get_angle_unit(pd)) != NULL)
        {
            add_xml_char(child, exml_names[exmlANGLE_UNIT], blu);
        }
        if ((func = gmx_poldata_get_dihedral_function(pd, i)) != NULL)
        {
            add_xml_char(child, exml_names[exmlFUNCTION], func);
        }
        while (gmx_poldata_get_dihedral(pd, i, &atom1, &atom2, &atom3, &atom4,
                                        &angle, &sigma, &ntrain, &params) > 0)
        {
            grandchild = add_xml_child(child, exml_names[ex[i]]);
            add_xml_char(grandchild, exml_names[exmlATOM1], atom1);
            add_xml_char(grandchild, exml_names[exmlATOM2], atom2);
            add_xml_char(grandchild, exml_names[exmlATOM3], atom3);
            add_xml_char(grandchild, exml_names[exmlATOM4], atom4);
            add_xml_double(grandchild, exml_names[exmlANGLE], angle);
            add_xml_double(grandchild, exml_names[exmlSIGMA], sigma);
            add_xml_int(grandchild, exml_names[exmlNTRAIN], ntrain);
            add_xml_char(grandchild, exml_names[exmlPARAMS], params);
            sfree(atom1);
            sfree(atom2);
            sfree(atom3);
            sfree(atom4);
            sfree(params);
        }
    }
    child = add_xml_child(parent, exml_names[exmlBSATOMS]);
    if ((tmp = gmx_poldata_get_bosque_unit(pd)) != NULL)
    {
        add_xml_char(child, exml_names[exmlPOLAR_UNIT], tmp);
    }

    while (1 == gmx_poldata_get_bosque(pd, &name, &polarizability))
    {
        grandchild = add_xml_child(child, exml_names[exmlBSATOM]);
        add_xml_char(grandchild, exml_names[exmlELEM], name);
        add_xml_double(grandchild, exml_names[exmlPOLARIZABILITY], polarizability);
    }
    child = add_xml_child(parent, exml_names[exmlMILATOMS]);
    gmx_poldata_get_miller_units(pd, &tau_unit, &ahp_unit);
    if (tau_unit)
    {
        add_xml_char(child, exml_names[exmlTAU_UNIT], tau_unit);
        sfree(tau_unit);
    }
    if (ahp_unit)
    {
        add_xml_char(child, exml_names[exmlAHP_UNIT], ahp_unit);
        sfree(ahp_unit);
    }

    while (1 == gmx_poldata_get_miller(pd, &name, &atomnumber, &tau_ahc, &alpha_ahp))
    {
        grandchild = add_xml_child(child, exml_names[exmlMILATOM]);
        add_xml_char(grandchild, exml_names[exmlMILNAME], name);
        add_xml_int(grandchild, exml_names[exmlATOMNUMBER], atomnumber);
        add_xml_double(grandchild, exml_names[exmlTAU_AHC], tau_ahc);
        add_xml_double(grandchild, exml_names[exmlALPHA_AHP], alpha_ahp);
        sfree(name);
    }

    child = add_xml_child(parent, exml_names[exmlSYMMETRIC_CHARGES]);

    while (gmx_poldata_get_symcharges(pd, &central, &attached, &numattach) == 1)
    {
        grandchild = add_xml_child(child, exml_names[exmlSYM_CHARGE]);
        add_xml_char(grandchild, exml_names[exmlCENTRAL], central);
        add_xml_char(grandchild, exml_names[exmlATTACHED], attached);
        add_xml_int(grandchild, exml_names[exmlNUMATTACH], numattach);
        sfree(central);
        sfree(attached);
    }

    child = add_xml_child(parent, exml_names[exmlEEMPROPS]);
    while (gmx_poldata_get_eemprops(pd, &model, &name, &J0, &chi0, &zeta, &qstr, &rowstr) == 1)
    {
        grandchild = add_xml_child(child, exml_names[exmlEEMPROP]);
        add_xml_char(grandchild, exml_names[exmlMODEL],
                     get_eemtype_name(model));
        add_xml_char(grandchild, exml_names[exmlNAME], name);
        add_xml_double(grandchild, exml_names[exmlJ0], J0);
        add_xml_double(grandchild, exml_names[exmlCHI0], chi0);
        add_xml_char(grandchild, exml_names[exmlZETA], zeta);
        add_xml_char(grandchild, exml_names[exmlCHARGES], qstr);
        add_xml_char(grandchild, exml_names[exmlROW], rowstr);
        sfree(zeta);
        sfree(qstr);
        sfree(rowstr);
    }
    while (gmx_poldata_list_epref(pd, &model, &epref) == 1)
    {
        grandchild = add_xml_child(child, exml_names[exmlEEMPROP_REF]);
        add_xml_char(grandchild, exml_names[exmlMODEL], get_eemtype_name(model));
        add_xml_char(grandchild, exml_names[exmlEPREF], epref);
        sfree(epref);
    }
}

void gmx_poldata_write(const char *fn, gmx_poldata_t pd,
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
    add_xml_poldata(myroot, pd);

    xmlSetDocCompressMode(doc, (int)bCompress);
    xmlIndentTreeOutput = 1;
    if (xmlSaveFormatFileEnc(fn, doc, "ISO-8859-1", 2) == 0)
    {
        gmx_fatal(FARGS, "Saving file", fn);
    }
    xmlFreeDoc(doc);
}
