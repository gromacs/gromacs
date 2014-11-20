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
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/cstringutil.h"
#include "xml_util.h"
#include "molprop.h"
#include "molprop_xml.h"
#include "stringutil.h"

extern int xmlDoValidityCheckingDefaultValue;

static gmx_bool NN(char *x)
{
    return (NULL != x);
//    return ((NULL != (x)) && (strlen(x) > 0));
}

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
    exmlMOLECULES, exmlMOLECULE, exmlFORMULA, exmlMOLNAME, exmlMASS,
    exmlMOLINFO, exmlIUPAC, exmlCAS, exmlCID, exmlINCHI,
    exmlMULTIPLICITY, exmlCHARGE,
    exmlCATEGORY, exmlCATNAME, exmlEXPERIMENT,
    exmlPOLARIZABILITY, exmlENERGY, exmlDIPOLE, exmlQUADRUPOLE, exmlPOTENTIAL,
    exmlNAME, exmlAVERAGE, exmlERROR,
    exmlMETHOD, exmlREFERENCE, exmlTYPE, exmlSOURCE,
    exmlBOND, exmlAI, exmlAJ, exmlBONDORDER,
    exmlCOMPOSITION, exmlCOMPNAME, exmlCATOM, exmlC_NAME, exmlC_NUMBER,
    exmlCALCULATION, exmlPROGRAM, exmlBASISSET, exmlCONFORMATION, exmlDATAFILE,
    exmlUNIT, exmlATOM, exmlATOMID, exmlOBTYPE, exmlX_UNIT, exmlV_UNIT, exmlESPID,
    exmlX, exmlY, exmlZ, exmlV, exmlXX, exmlYY, exmlZZ,
    exmlXY, exmlXZ, exmlYZ, exmlQ,
    exmlNR
};

static const char *exml_names[exmlNR] = {
    "molecules", "molecule", "formula", "molname", "mass",
    "molinfo", "iupac", "cas", "cid", "inchi",
    "multiplicity", "charge",
    "category", "catname", "experiment",
    "polarizability", "energy", "dipole", "quadrupole", "potential",
    "name", "average", "error",
    "method", "reference", "type", "source",
    "bond", "ai", "aj", "bondorder",
    "composition", "compname", "catom", "cname", "cnumber",
    "calculation", "program", "basisset", "conformation", "datafile",
    "unit", "atom", "atomid", "obtype", "coord_unit", "potential_unit", "espid",
    "x", "y", "z", "V", "xx", "yy", "zz", "xy", "xz", "yz", "q"
};

static void add_xml_string(xmlNodePtr ptr, const char *name, std::string val)
{
    if (xmlSetProp(ptr, (xmlChar *)name, (xmlChar *)val.c_str()) == 0)
    {
        gmx_fatal(FARGS, "Setting", (char *)name);
    }
}

static char *sp(int n, char buf[], int maxindent)
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

    return buf;
}

double my_atof(char *ptr)
{
    if (NULL != ptr)
    {
        return atof(ptr);
    }
    else
    {
        return 0;
    }
}

static void get_attributes(FILE *fp, gmx_bool bZero, int indent, xmlAttrPtr attr, char *xbuf[])
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

    while (attr != NULL)
    {
        attrname = (char *)attr->name;
        attrval  = (char *)attr->children->content;

#define atest(s) ((strcasecmp(attrname, s) == 0) && (attrval != NULL))
        if ((kkk = find_elem(attrname, exmlNR, exml_names)) != -1)
        {
            if (attrval != NULL)
            {
                xbuf[kkk] = strdup(attrval);
            }
        }
        if (fp)
        {
            fprintf(fp, "%sProperty: '%s' Value: '%s'\n", sp(indent, buf, 99),
                    attrname, attrval);
        }
        attr = attr->next;
#undef atest
    }
}

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

static void mp_process_tree(FILE *fp, xmlNodePtr tree,
                            int indent,
                            std::vector<alexandria::MolProp> &molprops,
                            gmx_bool *bExperiment)
{
    xmlNodePtr                   tc;
    char                         buf[100];
    alexandria::MolProp         *mpt;
    alexandria::CalcAtomIterator atom_it;
    gmx_bool                     bCompIt = FALSE;
    int                          i;
    char                        *xbuf[exmlNR];
    int                          node, elem = -1;

    while (tree != NULL)
    {
        if (fp)
        {
            if ((tree->type > 0) && ((unsigned)tree->type < NXMLTYPES))
            {
                fprintf(fp, "Node type %s encountered with name %s\n",
                        xmltypes[tree->type], (char *)tree->name);
            }
            else
            {
                fprintf(fp, "Node type %d encountered\n", tree->type);
            }
        }
        if (molprops.size() > 0)
        {
            mpt = &(molprops.back());
        }
        else
        {
            mpt = NULL;
        }
        switch (tree->type)
        {
            case XML_TEXT_NODE:
            {
                //   fprintf(stderr,"Text node %s encountered\n",(char *)tree->name);
                break;
            }
            case XML_ELEMENT_NODE:
            {
                if ((elem = find_elem((char *)tree->name, exmlNR, exml_names)) != -1)
                {
                    if (fp)
                    {
                        fprintf(fp, "%sElement node name %s\n", sp(indent, buf, 99),
                                (char *)tree->name);
                    }
                    get_attributes(fp, TRUE, indent, tree->properties, xbuf);

                    switch (elem)
                    {
                        case exmlMOLECULES:
                            break;
                        case exmlMOLECULE:
                        {
                            alexandria::MolProp mp;
                            if (NN(xbuf[exmlFORMULA]))
                            {
                                mp.SetFormula(xbuf[exmlFORMULA]);
                            }
                            if (NN(xbuf[exmlMOLNAME]))
                            {
                                mp.SetMolname(xbuf[exmlMOLNAME]);
                            }
                            if (NN(xbuf[exmlMASS]))
                            {
                                mp.SetMass(atof(xbuf[exmlMASS]));
                            }
                            if (NN(xbuf[exmlCHARGE]))
                            {
                                mp.SetCharge(atoi(xbuf[exmlCHARGE]));
                            }
                            if (NN(xbuf[exmlMULTIPLICITY]))
                            {
                                mp.SetMultiplicity(atoi(xbuf[exmlMULTIPLICITY]));
                            }
                            molprops.push_back(mp);
                            mpt = &(molprops.back());
                        }
                        break;
                        /* The items below are handled when treating attributes */
                        case exmlMOLINFO:
                            if (NN(xbuf[exmlIUPAC]))
                            {
                                mpt->SetIupac(xbuf[exmlIUPAC]);
                            }
                            if (NN(xbuf[exmlCAS]))
                            {
                                mpt->SetCas(xbuf[exmlCAS]);
                            }
                            if (NN(xbuf[exmlCID]))
                            {
                                mpt->SetCid(xbuf[exmlCID]);
                            }
                            if (NN(xbuf[exmlINCHI]))
                            {
                                mpt->SetInchi(xbuf[exmlINCHI]);
                            }
                            break;
                        case exmlCATEGORY:
                            if (NN(xbuf[exmlCATNAME]))
                            {
                                mpt->AddCategory(xbuf[exmlCATNAME]);
                            }
                            break;
                        case exmlPOLARIZABILITY:
                            process_children(tree->children, xbuf);
                            if (NN(xbuf[exmlTYPE])  && NN(xbuf[exmlUNIT]) &&
                                ((NN(xbuf[exmlAVERAGE]) && NN(xbuf[exmlERROR])) ||
                                 (NN(xbuf[exmlXX]) && NN(xbuf[exmlYY]) && NN(xbuf[exmlZZ]))))
                            {
                                alexandria::MolecularPolarizability mdp(xbuf[exmlTYPE], xbuf[exmlUNIT],
                                                                        NN(xbuf[exmlXX]) ? my_atof(xbuf[exmlXX]) : 0,
                                                                        NN(xbuf[exmlYY]) ? my_atof(xbuf[exmlYY]) : 0,
                                                                        NN(xbuf[exmlZZ]) ? my_atof(xbuf[exmlZZ]) : 0,
                                                                        NN(xbuf[exmlXY]) ? my_atof(xbuf[exmlXY]) : 0,
                                                                        NN(xbuf[exmlXZ]) ? my_atof(xbuf[exmlXZ]) : 0,
                                                                        NN(xbuf[exmlYZ]) ? my_atof(xbuf[exmlYZ]) : 0,
                                                                        my_atof(xbuf[exmlAVERAGE]), my_atof(xbuf[exmlERROR]));
                                if (*bExperiment)
                                {
                                    mpt->LastExperiment()->AddPolar(mdp);
                                }
                                else
                                {
                                    mpt->LastCalculation()->AddPolar(mdp);
                                }
                            }
                            break;
                        case exmlPOTENTIAL:
                            process_children(tree->children, xbuf);
                            if (NN(xbuf[exmlX_UNIT]) && NN(xbuf[exmlV_UNIT]) &&
                                NN(xbuf[exmlESPID]) &&
                                NN(xbuf[exmlX]) && NN(xbuf[exmlY]) &&
                                NN(xbuf[exmlZ]) && NN(xbuf[exmlV]))
                            {
                                alexandria::ElectrostaticPotential ep(xbuf[exmlX_UNIT], xbuf[exmlV_UNIT],
                                                                      atoi(xbuf[exmlESPID]),
                                                                      my_atof(xbuf[exmlX]), my_atof(xbuf[exmlY]),
                                                                      my_atof(xbuf[exmlZ]), my_atof(xbuf[exmlV]));
                                mpt->LastCalculation()->AddPotential(ep);
                            }
                            break;
                        case exmlDIPOLE:
                            process_children(tree->children, xbuf);
                            if (NN(xbuf[exmlTYPE]) && NN(xbuf[exmlUNIT]) &&
                                NN(xbuf[exmlAVERAGE]) && NN(xbuf[exmlERROR]))
                            {
                                alexandria::MolecularDipole mdp(xbuf[exmlTYPE], xbuf[exmlUNIT],
                                                                  NN(xbuf[exmlX]) ? my_atof(xbuf[exmlX]) : 0,
                                                                  NN(xbuf[exmlY]) ? my_atof(xbuf[exmlY]) : 0,
                                                                  NN(xbuf[exmlZ]) ? my_atof(xbuf[exmlZ]) : 0,
                                                                  my_atof(xbuf[exmlAVERAGE]), my_atof(xbuf[exmlERROR]));

                                if (*bExperiment)
                                {
                                    mpt->LastExperiment()->AddDipole(mdp);
                                }
                                else
                                {
                                    mpt->LastCalculation()->AddDipole(mdp);
                                }
                            }
                            break;
                        case exmlQUADRUPOLE:
                            process_children(tree->children, xbuf);
                            if (NN(xbuf[exmlTYPE]) && NN(xbuf[exmlUNIT]) &&
                                NN(xbuf[exmlXX]) && NN(xbuf[exmlYY]) && NN(xbuf[exmlZZ]) &&
                                NN(xbuf[exmlXY]) && NN(xbuf[exmlXZ]) && NN(xbuf[exmlYZ]))
                            {
                                alexandria::MolecularQuadrupole mq(xbuf[exmlTYPE], xbuf[exmlUNIT],
                                                                   my_atof(xbuf[exmlXX]), my_atof(xbuf[exmlYY]),
                                                                   my_atof(xbuf[exmlZZ]), my_atof(xbuf[exmlXY]),
                                                                   my_atof(xbuf[exmlXZ]), my_atof(xbuf[exmlYZ]));
                                if (*bExperiment)
                                {
                                    mpt->LastExperiment()->AddQuadrupole(mq);
                                }
                                else
                                {
                                    mpt->LastCalculation()->AddQuadrupole(mq);
                                }
                            }
                            break;
                        case exmlBOND:
                            process_children(tree->children, xbuf);
                            if (NN(xbuf[exmlAI]) && NN(xbuf[exmlAJ]) &&
                                NN(xbuf[exmlBONDORDER]))
                            {
                                alexandria::Bond b(atoi(xbuf[exmlAI]), atoi(xbuf[exmlAJ]),
                                                   atoi(xbuf[exmlBONDORDER]));
                                mpt->AddBond(b);
                            }
                            break;
                        case exmlENERGY:
                            process_children(tree, xbuf);
                            if (NN(xbuf[exmlTYPE]) && NN(xbuf[exmlUNIT]) &&
                                NN(xbuf[exmlENERGY]))
                            {
                                alexandria::MolecularEnergy me(xbuf[exmlTYPE], xbuf[exmlUNIT],
                                                               my_atof(xbuf[exmlENERGY]),
                                                               xbuf[exmlERROR] ? my_atof(xbuf[exmlERROR]) : 0.0);
                                if (*bExperiment)
                                {
                                    mpt->LastExperiment()->AddEnergy(me);
                                }
                                else
                                {
                                    mpt->LastCalculation()->AddEnergy(me);
                                }
                            }
                            break;

                        case exmlCOMPOSITION:
                            if (NN(xbuf[exmlCOMPNAME]))
                            {
                                alexandria::MolecularComposition mc(xbuf[exmlCOMPNAME]);
                                mpt->AddComposition(mc);
                                bCompIt = TRUE;
                            }
                            break;
                        case exmlCATOM:
                            if (NN(xbuf[exmlC_NAME]) && NN(xbuf[exmlC_NUMBER]) && bCompIt)
                            {
                                alexandria::AtomNum an(xbuf[exmlC_NAME], atoi(xbuf[exmlC_NUMBER]));
                                mpt->LastMolecularComposition()->AddAtom(an);
                            }
                            break;
                        case exmlCALCULATION:
                            if (NN(xbuf[exmlPROGRAM]) && NN(xbuf[exmlMETHOD]) &&
                                NN(xbuf[exmlBASISSET]) && NN(xbuf[exmlREFERENCE]) &&
                                NN(xbuf[exmlCONFORMATION]) && NN(xbuf[exmlDATAFILE]))
                            {
                                alexandria::Calculation mycalc(xbuf[exmlPROGRAM], xbuf[exmlMETHOD],
                                                               xbuf[exmlBASISSET], xbuf[exmlREFERENCE],
                                                               xbuf[exmlCONFORMATION], xbuf[exmlDATAFILE]);
                                mpt->AddCalculation(mycalc);
                                *bExperiment = FALSE;
                            }
                            else
                            {
                                gmx_fatal(FARGS, "Trying to add calculation with program %s, method %s, basisset %s and reference %s conformation %s datafile %s",
                                          xbuf[exmlPROGRAM], xbuf[exmlMETHOD],
                                          xbuf[exmlBASISSET], xbuf[exmlREFERENCE],
                                          xbuf[exmlCONFORMATION], xbuf[exmlDATAFILE]);
                            }
                            break;
                        case exmlATOM:
                            if (NN(xbuf[exmlNAME]) && NN(xbuf[exmlOBTYPE]) && NN(xbuf[exmlATOMID]))
                            {
                                alexandria::CalcAtom ca(xbuf[exmlNAME], xbuf[exmlOBTYPE], atoi(xbuf[exmlATOMID]));
                                sfree(xbuf[exmlNAME]);   xbuf[exmlNAME]   = NULL;
                                sfree(xbuf[exmlOBTYPE]); xbuf[exmlOBTYPE] = NULL;
                                sfree(xbuf[exmlATOMID]); xbuf[exmlATOMID] = NULL;
                                for (tc = tree->children; (NULL != tc); tc = tc->next)
                                {
                                    get_attributes(fp, FALSE, indent, tc->properties, xbuf);
                                    if (((node = find_elem((char *)tc->name, exmlNR, exml_names)) != -1) &&
                                        (NULL != tc->children) &&
                                        (NULL != tc->children->content))
                                    {
                                        xbuf[node] = strdup((char *)tc->children->content);
                                    }

                                    if (NN(xbuf[exmlX]) && NN(xbuf[exmlY]) && NN(xbuf[exmlZ])
                                        && NN(xbuf[exmlUNIT]))
                                    {
                                        ca.SetUnit(xbuf[exmlUNIT]);
                                        ca.SetCoords(my_atof(xbuf[exmlX]), my_atof(xbuf[exmlY]), my_atof(xbuf[exmlZ]));
                                        sfree(xbuf[exmlX]); xbuf[exmlX]       = NULL;
                                        sfree(xbuf[exmlY]); xbuf[exmlY]       = NULL;
                                        sfree(xbuf[exmlZ]); xbuf[exmlZ]       = NULL;
                                        sfree(xbuf[exmlUNIT]); xbuf[exmlUNIT] = NULL;
                                    }
                                    if (NN(xbuf[exmlQ]) && NN(xbuf[exmlTYPE]) && NN(xbuf[exmlUNIT]))
                                    {
                                        alexandria::AtomicCharge aq(xbuf[exmlTYPE], xbuf[exmlUNIT], my_atof(xbuf[exmlQ]));
                                        ca.AddCharge(aq);
                                        sfree(xbuf[exmlQ]);    xbuf[exmlQ]    = NULL;
                                        sfree(xbuf[exmlUNIT]); xbuf[exmlUNIT] = NULL;
                                        sfree(xbuf[exmlTYPE]); xbuf[exmlTYPE] = NULL;
                                    }
                                }
                                /* Now finally add the atom */
                                mpt->LastCalculation()->AddAtom(ca);
                            }
                            break;

                        case exmlEXPERIMENT:
                            if (NN(xbuf[exmlREFERENCE]) && NN(xbuf[exmlCONFORMATION]))
                            {
                                alexandria::Experiment myexp(xbuf[exmlREFERENCE], xbuf[exmlCONFORMATION]);
                                mpt->AddExperiment(myexp);
                                *bExperiment = TRUE;
                            }
                            else
                            {
                                gmx_fatal(FARGS, "Experimental data without reference");
                            }
                            break;
                        default:
                            break;
                    }
                }
                for (i = 0; (i < exmlNR); i++)
                {
                    if (NN(xbuf[i]))
                    {
                        sfree(xbuf[i]);
                    }
                }
                if (tree->children)
                {
                    if ((elem = find_elem((char *)tree->name, exmlNR, exml_names)) != -1)
                    {
                        mp_process_tree(fp, tree->children, indent+2,
                                        molprops, bExperiment);
                    }
                }
                break;
            }
            default:
                break;
        }
        tree = tree->next;
    }
}

void MolPropRead(const char *fn, std::vector<alexandria::MolProp> &mpt)
{
    xmlDocPtr     doc;
    const char   *db          = "alexandria.ff/molprops.dat";
    gmx_bool      bExperiment = FALSE;

    xmlDoValidityCheckingDefaultValue = 0;
    if (NULL == fn)
    {
        fn = gmxlibfn(db);
    }
    if (NULL != debug)
    {
        fprintf(debug, "Opening %s\n", fn);
    }
    if ((doc = xmlParseFile(fn)) == NULL)
    {
        fprintf(stderr, "Reading XML file %s. Run a syntax checker such as nsgmls.",
                fn);
        exit(1);
    }

    mp_process_tree(NULL, doc->children, 0,
                    mpt, &bExperiment);

    xmlFreeDoc(doc);
}

static void add_exper_properties(xmlNodePtr              exp,
                                 alexandria::Experiment &exper)
{
    xmlNodePtr  child;
    const char *ptr;
    double      value, error, x, y, z, xx, yy, zz, xy, xz, yz;

    for (alexandria::MolecularEnergyIterator me_it = exper.BeginEnergy();
         (me_it < exper.EndEnergy()); me_it++)
    {
        me_it->Get(&value, &error);

        ptr   = gmx_ftoa(value).c_str();
        child = add_xml_child_val(exp, exml_names[exmlENERGY], ptr);
        add_xml_string(child, exml_names[exmlTYPE], me_it->GetType());
        add_xml_string(child, exml_names[exmlUNIT], me_it->GetUnit());
    }

    for (alexandria::MolecularDipoleIterator mdp_it = exper.BeginDipole();
         (mdp_it < exper.EndDipole()); mdp_it++)
    {
        mdp_it->Get(&x, &y, &z, &value, &error);

        child = add_xml_child(exp, exml_names[exmlDIPOLE]);
        add_xml_string(child, exml_names[exmlTYPE], mdp_it->GetType());
        add_xml_string(child, exml_names[exmlUNIT], mdp_it->GetUnit());
        ptr = gmx_ftoa(value).c_str();
        add_xml_child_val(child, exml_names[exmlAVERAGE], ptr);
        ptr = gmx_ftoa(error).c_str();
        add_xml_child_val(child, exml_names[exmlERROR], ptr);
        if ((x != 0) || (y != 0) || (z != 0))
        {
            ptr = gmx_ftoa(x).c_str();
            add_xml_child_val(child, exml_names[exmlX], ptr);
            ptr = gmx_ftoa(y).c_str();
            add_xml_child_val(child, exml_names[exmlY], ptr);
            ptr = gmx_ftoa(z).c_str();
            add_xml_child_val(child, exml_names[exmlZ], ptr);
        }
    }

    for (alexandria::MolecularPolarizabilityIterator mdp_it = exper.BeginPolar(); (mdp_it < exper.EndPolar()); mdp_it++)
    {
        mdp_it->Get(&xx, &yy, &zz, &xy, &xz, &yz, &value, &error);

        child = add_xml_child(exp, exml_names[exmlPOLARIZABILITY]);
        add_xml_string(child, exml_names[exmlTYPE], mdp_it->GetType());
        add_xml_string(child, exml_names[exmlUNIT], mdp_it->GetUnit());
        ptr = gmx_ftoa(value).c_str();
        add_xml_child_val(child, exml_names[exmlAVERAGE], ptr);
        ptr = gmx_ftoa(error).c_str();
        add_xml_child_val(child, exml_names[exmlERROR], ptr);
        if ((xx != 0) || (yy != 0) || (zz != 0))
        {
            ptr = gmx_ftoa(xx).c_str();
            add_xml_child_val(child, exml_names[exmlXX], ptr);
            ptr = gmx_ftoa(yy).c_str();
            add_xml_child_val(child, exml_names[exmlYY], ptr);
            ptr = gmx_ftoa(zz).c_str();
            add_xml_child_val(child, exml_names[exmlZZ], ptr);
        }
        if ((xy != 0) || (xz != 0) || (yz != 0))
        {
            ptr = gmx_ftoa(xy).c_str();
            add_xml_child_val(child, exml_names[exmlXY], ptr);
            ptr = gmx_ftoa(xz).c_str();
            add_xml_child_val(child, exml_names[exmlXZ], ptr);
            ptr = gmx_ftoa(yz).c_str();
            add_xml_child_val(child, exml_names[exmlYZ], ptr);
        }
    }

    for (alexandria::MolecularQuadrupoleIterator mq_it = exper.BeginQuadrupole();
         (mq_it < exper.EndQuadrupole()); mq_it++)
    {
        mq_it->Get(&xx, &yy, &zz, &xy, &xz, &yz);

        child = add_xml_child(exp, exml_names[exmlQUADRUPOLE]);
        add_xml_string(child, exml_names[exmlTYPE], mq_it->GetType());
        add_xml_string(child, exml_names[exmlUNIT], mq_it->GetUnit());
        ptr = gmx_ftoa(xx).c_str();
        add_xml_child_val(child, exml_names[exmlXX], ptr);
        ptr = gmx_ftoa(yy).c_str();
        add_xml_child_val(child, exml_names[exmlYY], ptr);
        ptr = gmx_ftoa(zz).c_str();
        add_xml_child_val(child, exml_names[exmlZZ], ptr);
        ptr = gmx_ftoa(xy).c_str();
        add_xml_child_val(child, exml_names[exmlXY], ptr);
        ptr = gmx_ftoa(xz).c_str();
        add_xml_child_val(child, exml_names[exmlXZ], ptr);
        ptr = gmx_ftoa(yz).c_str();
        add_xml_child_val(child, exml_names[exmlYZ], ptr);
    }
}

static void add_calc_properties(xmlNodePtr               exp,
                                alexandria::Calculation &calc)
{
    for (alexandria::ElectrostaticPotentialIterator ep_it = calc.BeginPotential();
         (ep_it < calc.EndPotential()); ep_it++)
    {
        char  *x_unit, *v_unit;
        double x, y, z, V;
        int    espid;

        ep_it->Get(&x_unit, &v_unit, &espid, &x, &y, &z, &V);

        xmlNodePtr child = add_xml_child(exp, exml_names[exmlPOTENTIAL]);
        add_xml_char(child, exml_names[exmlX_UNIT], x_unit);
        add_xml_char(child, exml_names[exmlV_UNIT], v_unit);
        add_xml_int(child, exml_names[exmlESPID], espid);
        if ((x != 0) || (y != 0) || (z != 0) || (V != 0))
        {
            const char *ptr = gmx_ftoa(x).c_str();
            add_xml_child_val(child, exml_names[exmlX], ptr);
            ptr = gmx_ftoa(y).c_str();
            add_xml_child_val(child, exml_names[exmlY], ptr);
            ptr = gmx_ftoa(z).c_str();
            add_xml_child_val(child, exml_names[exmlZ], ptr);
            ptr = gmx_ftoa(V).c_str();
            add_xml_child_val(child, exml_names[exmlV], ptr);
        }
        sfree(x_unit);
        sfree(v_unit);
    }
}

static void add_xml_molprop(xmlNodePtr                                 parent,
                            std::vector<alexandria::MolProp>::iterator mp_it)
{
    xmlNodePtr ptr = add_xml_child(parent, exml_names[exmlMOLECULE]);
    add_xml_string(ptr, exml_names[exmlMOLNAME], mp_it->GetMolname());
    add_xml_string(ptr, exml_names[exmlFORMULA], mp_it->GetFormula());
    add_xml_double(ptr, exml_names[exmlMASS], mp_it->GetMass());
    add_xml_double(ptr, exml_names[exmlCHARGE], mp_it->GetCharge());
    add_xml_double(ptr, exml_names[exmlMULTIPLICITY], mp_it->GetMultiplicity());

    xmlNodePtr child = add_xml_child(ptr, exml_names[exmlMOLINFO]);
    add_xml_string(child, exml_names[exmlIUPAC], mp_it->GetIupac());
    add_xml_string(child, exml_names[exmlCAS], mp_it->GetCas());
    add_xml_string(child, exml_names[exmlCID], mp_it->GetCid());
    add_xml_string(child, exml_names[exmlINCHI], mp_it->GetInchi());

    for (alexandria::BondIterator b_it = mp_it->BeginBond();
         (b_it < mp_it->EndBond()); b_it++)
    {
        xmlNodePtr child = add_xml_child(ptr, exml_names[exmlBOND]);
        add_xml_int(child, exml_names[exmlAI], b_it->GetAi());
        add_xml_int(child, exml_names[exmlAJ], b_it->GetAj());
        add_xml_int(child, exml_names[exmlBONDORDER], b_it->GetBondOrder());
    }

    for (alexandria::ExperimentIterator e_it = mp_it->BeginExperiment();
         (e_it < mp_it->EndExperiment()); e_it++)
    {
        xmlNodePtr child = add_xml_child(ptr, exml_names[exmlEXPERIMENT]);
        add_xml_string(child, exml_names[exmlREFERENCE], e_it->GetReference());
        add_xml_string(child, exml_names[exmlCONFORMATION], e_it->GetConformation());

        add_exper_properties(child, *e_it);
    }

    for (alexandria::CalculationIterator c_it = mp_it->BeginCalculation();
         (c_it < mp_it->EndCalculation()); c_it++)
    {
        xmlNodePtr child = add_xml_child(ptr, exml_names[exmlCALCULATION]);
        add_xml_string(child, exml_names[exmlPROGRAM], c_it->GetProgram());
        add_xml_string(child, exml_names[exmlMETHOD], c_it->GetMethod());
        add_xml_string(child, exml_names[exmlBASISSET], c_it->GetBasisset());
        add_xml_string(child, exml_names[exmlREFERENCE], c_it->GetReference());
        add_xml_string(child, exml_names[exmlCONFORMATION], c_it->GetConformation());
        add_xml_string(child, exml_names[exmlDATAFILE], c_it->GetDatafile());

        add_exper_properties(child, *c_it);
        add_calc_properties(child, *c_it);

        for (alexandria::CalcAtomIterator ca_it = c_it->BeginAtom();
             (ca_it < c_it->EndAtom()); ca_it++)
        {
            xmlNodePtr grandchild = add_xml_child(child, exml_names[exmlATOM]);
            add_xml_string(grandchild, exml_names[exmlNAME], ca_it->GetName());
            add_xml_string(grandchild, exml_names[exmlOBTYPE], ca_it->GetObtype());
            add_xml_int(grandchild, exml_names[exmlATOMID], ca_it->GetAtomid());

            double x, y, z;
            ca_it->GetCoords(&x, &y, &z);

            const char *p    = gmx_ftoa(x).c_str();
            xmlNodePtr  baby = add_xml_child_val(grandchild, exml_names[exmlX], p);
            add_xml_string(baby, exml_names[exmlUNIT], ca_it->GetUnit());
            p    = gmx_ftoa(y).c_str();
            baby = add_xml_child_val(grandchild, exml_names[exmlY], p);
            add_xml_string(baby, exml_names[exmlUNIT], ca_it->GetUnit());
            p    = gmx_ftoa(z).c_str();
            baby = add_xml_child_val(grandchild, exml_names[exmlZ], p);
            add_xml_string(baby, exml_names[exmlUNIT], ca_it->GetUnit());

            for (alexandria::AtomicChargeIterator q_it = ca_it->BeginQ();
                 (q_it < ca_it->EndQ()); q_it++)
            {
                p = gmx_ftoa(q_it->GetQ()).c_str();
                xmlNodePtr atomptr = add_xml_child_val(grandchild, exml_names[exmlQ], p);
                add_xml_string(atomptr, exml_names[exmlTYPE], q_it->GetType());
                add_xml_string(atomptr, exml_names[exmlUNIT], q_it->GetUnit());
            }
        }
    }
    for (std::vector<std::string>::iterator s_it = mp_it->BeginCategory();
         (s_it < mp_it->EndCategory()); s_it++)
    {
        xmlNodePtr child = add_xml_child(ptr, exml_names[exmlCATEGORY]);
        add_xml_string(child, exml_names[exmlCATNAME], *s_it);
    }

    for (alexandria::MolecularCompositionIterator mc_it = mp_it->BeginMolecularComposition();
         (mc_it < mp_it->EndMolecularComposition()); mc_it++)
    {
        xmlNodePtr child = add_xml_child(ptr, exml_names[exmlCOMPOSITION]);
        add_xml_string(child, exml_names[exmlCOMPNAME], mc_it->GetCompName());
        for (alexandria::AtomNumIterator an_it = mc_it->BeginAtomNum();
             (an_it < mc_it->EndAtomNum()); an_it++)
        {
            xmlNodePtr grandchild = add_xml_child(child, exml_names[exmlCATOM]);
            add_xml_string(grandchild, exml_names[exmlC_NAME], an_it->GetAtom());
            add_xml_int(grandchild, exml_names[exmlC_NUMBER], an_it->GetNumber());
        }
    }
}

void MolPropWrite(const char *fn, std::vector<alexandria::MolProp> mpt, gmx_bool bCompress)
{
    xmlDocPtr                   doc;
    xmlDtdPtr                   dtd;
    xmlNodePtr                  myroot;
    xmlChar                    *libdtdname, *dtdname, *gmx;
    alexandria::MolPropIterator mp_it;

    gmx        = (xmlChar *) "molecules";
    dtdname    = (xmlChar *) "molprops.dtd";
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
    for (mp_it = mpt.begin(); (mp_it < mpt.end()); mp_it++)
    {
        if (NULL != debug)
        {
            fprintf(debug, "Adding %d/%d %s\n", (int)(mp_it - mpt.begin()),
                    (int)mpt.size(), mp_it->GetMolname().c_str());
        }
        //mp_it->Stats();
        add_xml_molprop(myroot, mp_it);
    }
    xmlSetDocCompressMode(doc, (int)bCompress);
    xmlIndentTreeOutput = 1;
    if (xmlSaveFormatFileEnc(fn, doc, "ISO-8859-1", 1) == 0)
    {
        gmx_fatal(FARGS, "Saving file", fn);
    }
    xmlFreeDoc(doc);
}
