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

#include "molprop_xml.h"

#include <stdlib.h>
#include <string.h>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"

#include "molprop.h"
#include "stringutil.h"
#include "xml_util.h"


static const char *job_name[alexandria::JOB_NR] = 
{
    "Opt", "Pop", "POLAR", "G2", "G3", 
    "G4", "CBSQB3", "W1U", "W1BD", "SP", "unknown"
};

const char *alexandria::jobType2string(alexandria::jobType jType)

{
    if (jType < alexandria::JOB_NR)
    {
        return job_name[jType];
    }
    return nullptr;
}

alexandria::jobType alexandria::string2jobType(const std::string &str)
{
    int i;

    for (i = 0; (i < alexandria::JOB_NR); i++)
    {
        if (str.compare(job_name[i]) == 0)
        {
            return static_cast<alexandria::jobType>(i);
        }
    }
    return alexandria::JOB_NR;
}

static bool NN(const std::string &s)
{
    return s.size() > 0;
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
#define NXMLTYPES sizeof(xmltypes)/sizeof(xmltypes[0])

enum {
    exmlMOLECULES, exmlMOLECULE, exmlFORMULA, exmlMOLNAME, exmlMASS,
    exmlMOLINFO, exmlIUPAC, exmlCAS, exmlCID, exmlINCHI,
    exmlMULTIPLICITY, exmlCHARGE,
    exmlCATEGORY, exmlCATNAME, exmlEXPERIMENT,
    exmlPOLARIZABILITY, exmlENERGY, exmlDIPOLE, exmlQUADRUPOLE, exmlPOTENTIAL,
    exmlNAME, exmlAVERAGE, exmlERROR, exmlTEMPERATURE, exmlPHASE,
    exmlMETHOD, exmlREFERENCE, exmlTYPE, exmlSOURCE,
    exmlBOND, exmlAI, exmlAJ, exmlBONDORDER,
    exmlCOMPOSITION, exmlCOMPNAME, exmlCATOM, exmlC_NAME, exmlC_NUMBER,
    exmlDATASOURCE, exmlPROGRAM, exmlBASISSET, exmlJOBTYPE, exmlCONFORMATION, exmlDATAFILE,
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
    "name", "average", "error", "temperature", "phase",
    "method", "reference", "type", "source",
    "bond", "ai", "aj", "bondorder",
    "composition", "compname", "catom", "cname", "cnumber",
    "datasource", "program", "basisset", "jobtype", "conformation", "datafile",
    "unit", "atom", "atomid", "obtype", "coord_unit", "potential_unit", "espid",
    "x", "y", "z", "V", "xx", "yy", "zz", "xy", "xz", "yz", "q"
};

static void add_xml_string(xmlNodePtr ptr, const char *name, std::string val)
{
    if (xmlSetProp(ptr, xmlCharStrdup(name), xmlCharStrdup(val.c_str())) == 0)
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

static double my_atof(const std::string &ptr)
{
    if (ptr.size() > 0)
    {
        return atof(ptr.c_str());
    }
    else
    {
        return 0;
    }
}

static void get_attributes(FILE *fp, gmx_bool bZero, int indent, xmlAttrPtr attr,
                           std::vector<std::string> &xbuf)
{
    if (bZero)
    {
        for (auto &x : xbuf)
        {
            x.clear();
        }
    }

    while (attr != NULL)
    {
        char *attrname = (char *)attr->name;
        char *attrval  = (char *)attr->children->content;

        int   kkk;
        if ((kkk = find_elem(attrname, exmlNR, exml_names)) != -1)
        {
            if (attrval != NULL)
            {
                xbuf[kkk].assign(attrval);
            }
        }
        if (fp)
        {
            char  buf[100];

            fprintf(fp, "%sProperty: '%s' Value: '%s'\n",
                    sp(indent, buf, sizeof(buf)-1),
                    attrname, attrval);
        }
        attr = attr->next;
    }
}

static void process_children(xmlNodePtr tree, std::vector<std::string> &xbuf)
{
    int node;

    while (NULL != tree)
    {
        if (((node = find_elem((char *)tree->name, exmlNR, exml_names)) != -1) &&
            (NULL != tree->children) &&
            (NULL != tree->children->content))
        {
            if (xbuf[node].size() == 0)
            {
                xbuf[node].assign(reinterpret_cast<char *>(tree->children->content));
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
    std::vector<std::string>     xbuf;
    int                          node, elem = -1;
    std::string                  xxx;

    xxx.clear();
    xbuf.resize(exmlNR, xxx);
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
            mpt = nullptr;
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

                    alexandria::Experiment *last = nullptr;
                    if (nullptr != mpt)
                    {
                        last = mpt->LastExperiment();
                    }
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
                                mp.SetMass(my_atof(xbuf[exmlMASS]));
                            }
                            if (NN(xbuf[exmlCHARGE]))
                            {
                                mp.SetCharge(atoi(xbuf[exmlCHARGE].c_str()));
                            }
                            if (NN(xbuf[exmlMULTIPLICITY]))
                            {
                                mp.SetMultiplicity(atoi(xbuf[exmlMULTIPLICITY].c_str()));
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
                            if ((nullptr != last) &&
                                NN(xbuf[exmlTYPE])  && NN(xbuf[exmlUNIT]) &&
                                NN(xbuf[exmlTEMPERATURE]) &&
                                ((NN(xbuf[exmlAVERAGE]) && NN(xbuf[exmlERROR])) ||
                                 (NN(xbuf[exmlXX]) && NN(xbuf[exmlYY]) && NN(xbuf[exmlZZ]))))
                            {
                                alexandria::MolecularPolarizability mdp(xbuf[exmlTYPE], xbuf[exmlUNIT],
                                                                        my_atof(xbuf[exmlTEMPERATURE]),
                                                                        my_atof(xbuf[exmlXX]),
                                                                        my_atof(xbuf[exmlYY]),
                                                                        my_atof(xbuf[exmlZZ]),
                                                                        my_atof(xbuf[exmlXY]),
                                                                        my_atof(xbuf[exmlXZ]),
                                                                        my_atof(xbuf[exmlYZ]),
                                                                        my_atof(xbuf[exmlAVERAGE]),
                                                                        my_atof(xbuf[exmlERROR]));
                                last->AddPolar(mdp);
                            }
                            break;
                        case exmlPOTENTIAL:
                            process_children(tree->children, xbuf);
                            if ((nullptr != last) &&
                                NN(xbuf[exmlX_UNIT]) && NN(xbuf[exmlV_UNIT]) &&
                                NN(xbuf[exmlESPID]) &&
                                NN(xbuf[exmlX]) && NN(xbuf[exmlY]) &&
                                NN(xbuf[exmlZ]) && NN(xbuf[exmlV]))
                            {
                                alexandria::ElectrostaticPotential ep(xbuf[exmlX_UNIT], xbuf[exmlV_UNIT],
                                                                      atoi(xbuf[exmlESPID].c_str()),
                                                                      my_atof(xbuf[exmlX]), my_atof(xbuf[exmlY]),
                                                                      my_atof(xbuf[exmlZ]), my_atof(xbuf[exmlV]));
                                last->AddPotential(ep);
                            }
                            break;
                        case exmlDIPOLE:
                            process_children(tree->children, xbuf);
                            if ((nullptr != last) &&
                                NN(xbuf[exmlTYPE]) && NN(xbuf[exmlUNIT]) &&
                                NN(xbuf[exmlAVERAGE]) && NN(xbuf[exmlERROR]) &&
                                NN(xbuf[exmlTEMPERATURE]))
                            {
                                alexandria::MolecularDipole mdp(xbuf[exmlTYPE], xbuf[exmlUNIT],
                                                                my_atof(xbuf[exmlTEMPERATURE]),
                                                                my_atof(xbuf[exmlX]),
                                                                my_atof(xbuf[exmlY]),
                                                                my_atof(xbuf[exmlZ]),
                                                                my_atof(xbuf[exmlAVERAGE]),
                                                                my_atof(xbuf[exmlERROR]));

                                last->AddDipole(mdp);
                            }
                            break;
                        case exmlQUADRUPOLE:
                            process_children(tree->children, xbuf);
                            if ((nullptr != last) &&
                                NN(xbuf[exmlTYPE]) && NN(xbuf[exmlUNIT]) &&
                                NN(xbuf[exmlTEMPERATURE]) &&
                                NN(xbuf[exmlXX]) && NN(xbuf[exmlYY]) && NN(xbuf[exmlZZ]) &&
                                NN(xbuf[exmlXY]) && NN(xbuf[exmlXZ]) && NN(xbuf[exmlYZ]))
                            {
                                alexandria::MolecularQuadrupole mq(xbuf[exmlTYPE], xbuf[exmlUNIT],
                                                                   my_atof(xbuf[exmlTEMPERATURE]),
                                                                   my_atof(xbuf[exmlXX]), my_atof(xbuf[exmlYY]),
                                                                   my_atof(xbuf[exmlZZ]), my_atof(xbuf[exmlXY]),
                                                                   my_atof(xbuf[exmlXZ]), my_atof(xbuf[exmlYZ]));
                                last->AddQuadrupole(mq);
                            }
                            break;
                        case exmlBOND:
                            process_children(tree->children, xbuf);
                            if (NN(xbuf[exmlAI]) && NN(xbuf[exmlAJ]) &&
                                NN(xbuf[exmlBONDORDER]))
                            {
                                alexandria::Bond b(atoi(xbuf[exmlAI].c_str()), atoi(xbuf[exmlAJ].c_str()),
                                                   atoi(xbuf[exmlBONDORDER].c_str()));
                                mpt->AddBond(b);
                            }
                            break;
                        case exmlENERGY:
                            process_children(tree, xbuf);
                            if ((nullptr != last) &&
                                NN(xbuf[exmlTYPE]) && NN(xbuf[exmlUNIT]) &&
                                NN(xbuf[exmlENERGY]) && NN(xbuf[exmlTEMPERATURE]) &&
                                NN(xbuf[exmlPHASE]))
                            {
                                alexandria::MolecularEnergy me(xbuf[exmlTYPE],
                                                               xbuf[exmlUNIT],
                                                               my_atof(xbuf[exmlTEMPERATURE]),
                                                               string2phase(xbuf[exmlPHASE]),
                                                               my_atof(xbuf[exmlENERGY]),
                                                               my_atof(xbuf[exmlERROR]));
                                last->AddEnergy(me);
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
                                alexandria::AtomNum               an(xbuf[exmlC_NAME], atoi(xbuf[exmlC_NUMBER].c_str()));
                                alexandria::MolecularComposition *l = mpt->LastMolecularComposition();
                                if (nullptr != l)
                                {
                                    l->AddAtom(an);
                                }
                            }
                            break;
                        case exmlATOM:
                            if ((nullptr != last) &&
                                NN(xbuf[exmlNAME]) && NN(xbuf[exmlOBTYPE]) && NN(xbuf[exmlATOMID]))
                            {
                                alexandria::CalcAtom ca(xbuf[exmlNAME], xbuf[exmlOBTYPE],
                                                        atoi(xbuf[exmlATOMID].c_str()));
                                xbuf[exmlNAME].clear();
                                xbuf[exmlOBTYPE].clear();
                                xbuf[exmlATOMID].clear();
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
                                        xbuf[exmlX].clear();
                                        xbuf[exmlY].clear();
                                        xbuf[exmlZ].clear();
                                        xbuf[exmlUNIT].clear();
                                    }
                                    if (NN(xbuf[exmlQ]) && NN(xbuf[exmlTYPE]) &&
                                        NN(xbuf[exmlTEMPERATURE]) &&
                                        NN(xbuf[exmlUNIT]))
                                    {
                                        alexandria::AtomicCharge aq(xbuf[exmlTYPE], xbuf[exmlUNIT],
                                                                    my_atof(xbuf[exmlTEMPERATURE]),
                                                                    my_atof(xbuf[exmlQ]));
                                        ca.AddCharge(aq);
                                        xbuf[exmlQ].clear();
                                        xbuf[exmlUNIT].clear();
                                        xbuf[exmlTYPE].clear();
                                    }
                                }
                                /* Now finally add the atom */
                                last->AddAtom(ca);
                            }
                            break;

                        case exmlEXPERIMENT:
                            if (NN(xbuf[exmlDATASOURCE]))
                            {
                                alexandria::DataSource ds = alexandria::dataSourceFromName(xbuf[exmlDATASOURCE]);

                                if (ds == alexandria::dsTheory &&
                                    NN(xbuf[exmlPROGRAM]) && NN(xbuf[exmlMETHOD]) &&
                                    NN(xbuf[exmlBASISSET]) && NN(xbuf[exmlREFERENCE]) &&
                                    NN(xbuf[exmlCONFORMATION]) && NN(xbuf[exmlDATAFILE]) &&
                                    NN(xbuf[exmlJOBTYPE]))
                                {
                                    alexandria::Experiment mycalc(xbuf[exmlPROGRAM], xbuf[exmlMETHOD],
                                                                  xbuf[exmlBASISSET], xbuf[exmlREFERENCE],
                                                                  xbuf[exmlCONFORMATION], xbuf[exmlDATAFILE],
                                                                  alexandria::string2jobType(xbuf[exmlJOBTYPE]));
                                    mpt->AddExperiment(mycalc);
                                }
                                else if (ds == alexandria::dsExperiment)
                                {
                                    if (NN(xbuf[exmlREFERENCE]))
                                    {
                                        const char *unknown = "unknown";
                                        alexandria::Experiment myexp(xbuf[exmlREFERENCE], 
                                                                     NN(xbuf[exmlCONFORMATION]) ? xbuf[exmlCONFORMATION] : unknown);
                                        mpt->AddExperiment(myexp);
                                    }
                                    else
                                    {
                                        gmx_fatal(FARGS, "Experimental data without reference");
                                    }
                                }
                                break;
                            }
                        default:
                            break;
                    }
                }
                for (auto &i : xbuf)
                {
                    if (NN(i))
                    {
                        i.clear();
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
        me_it->get(&value, &error);

        ptr   = gmx_ftoa(value).c_str();
        child = add_xml_child_val(exp, exml_names[exmlENERGY], ptr);
        add_xml_string(child, exml_names[exmlTYPE], me_it->getType());
        add_xml_string(child, exml_names[exmlUNIT], me_it->getUnit());
        add_xml_double(child, exml_names[exmlTEMPERATURE], me_it->getTemperature());
        add_xml_string(child, exml_names[exmlPHASE], phase2string(me_it->getPhase()));
    }

    for (alexandria::MolecularDipoleIterator mdp_it = exper.BeginDipole();
         (mdp_it < exper.EndDipole()); mdp_it++)
    {
        mdp_it->get(&x, &y, &z, &value, &error);

        child = add_xml_child(exp, exml_names[exmlDIPOLE]);
        add_xml_string(child, exml_names[exmlTYPE], mdp_it->getType());
        add_xml_string(child, exml_names[exmlUNIT], mdp_it->getUnit());
        add_xml_double(child, exml_names[exmlTEMPERATURE], mdp_it->getTemperature());
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
        mdp_it->get(&xx, &yy, &zz, &xy, &xz, &yz, &value, &error);

        child = add_xml_child(exp, exml_names[exmlPOLARIZABILITY]);
        add_xml_string(child, exml_names[exmlTYPE], mdp_it->getType());
        add_xml_string(child, exml_names[exmlUNIT], mdp_it->getUnit());
        add_xml_double(child, exml_names[exmlTEMPERATURE], mdp_it->getTemperature());
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
        mq_it->get(&xx, &yy, &zz, &xy, &xz, &yz);

        child = add_xml_child(exp, exml_names[exmlQUADRUPOLE]);
        add_xml_string(child, exml_names[exmlTYPE], mq_it->getType());
        add_xml_string(child, exml_names[exmlUNIT], mq_it->getUnit());
        add_xml_double(child, exml_names[exmlTEMPERATURE], mq_it->getTemperature());
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

static void add_calc_properties(xmlNodePtr              exp,
                                alexandria::Experiment &calc)
{
    for (alexandria::ElectrostaticPotentialIterator ep_it = calc.BeginPotential();
         (ep_it < calc.EndPotential()); ep_it++)
    {
        char  *x_unit, *v_unit;
        double x, y, z, V;
        int    espid;

        ep_it->get(&x_unit, &v_unit, &espid, &x, &y, &z, &V);

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
        free(x_unit);
        free(v_unit);
    }
}

static void add_xml_molprop(xmlNodePtr                                 parent,
                            std::vector<alexandria::MolProp>::iterator mp_it)
{
    xmlNodePtr ptr = add_xml_child(parent, exml_names[exmlMOLECULE]);
    add_xml_string(ptr, exml_names[exmlMOLNAME], mp_it->getMolname());
    add_xml_string(ptr, exml_names[exmlFORMULA], mp_it->formula());
    add_xml_double(ptr, exml_names[exmlMASS], mp_it->getMass());
    add_xml_double(ptr, exml_names[exmlCHARGE], mp_it->getCharge());
    add_xml_double(ptr, exml_names[exmlMULTIPLICITY], mp_it->getMultiplicity());

    xmlNodePtr child = add_xml_child(ptr, exml_names[exmlMOLINFO]);
    add_xml_string(child, exml_names[exmlIUPAC], mp_it->getIupac());
    add_xml_string(child, exml_names[exmlCAS], mp_it->getCas());
    add_xml_string(child, exml_names[exmlCID], mp_it->getCid());
    add_xml_string(child, exml_names[exmlINCHI], mp_it->getInchi());

    for (alexandria::BondIterator b_it = mp_it->BeginBond();
         (b_it < mp_it->EndBond()); b_it++)
    {
        xmlNodePtr child = add_xml_child(ptr, exml_names[exmlBOND]);
        add_xml_int(child, exml_names[exmlAI], b_it->getAi());
        add_xml_int(child, exml_names[exmlAJ], b_it->getAj());
        add_xml_int(child, exml_names[exmlBONDORDER], b_it->getBondOrder());
    }

    for (alexandria::ExperimentIterator e_it = mp_it->BeginExperiment();
         (e_it < mp_it->EndExperiment()); e_it++)
    {
        xmlNodePtr             child = add_xml_child(ptr, exml_names[exmlEXPERIMENT]);
        alexandria::DataSource ds    = e_it->dataSource();
        add_xml_string(child, exml_names[exmlDATASOURCE], dataSourceName(ds));
        add_xml_string(child, exml_names[exmlREFERENCE], e_it->getReference());
        add_xml_string(child, exml_names[exmlCONFORMATION], e_it->getConformation());
        if (alexandria::dsTheory == ds)
        {
            add_xml_string(child, exml_names[exmlPROGRAM], e_it->getProgram());
            add_xml_string(child, exml_names[exmlMETHOD], e_it->getMethod());
            add_xml_string(child, exml_names[exmlBASISSET], e_it->getBasisset());
            add_xml_string(child, exml_names[exmlJOBTYPE], jobType2string(e_it->getJobtype()));
            add_xml_string(child, exml_names[exmlDATAFILE], e_it->getDatafile());
        }

        add_exper_properties(child, *e_it);
        add_calc_properties(child, *e_it);

        for (alexandria::CalcAtomIterator ca_it = e_it->BeginAtom();
             (ca_it < e_it->EndAtom()); ca_it++)
        {
            xmlNodePtr grandchild = add_xml_child(child, exml_names[exmlATOM]);
            add_xml_string(grandchild, exml_names[exmlNAME], ca_it->getName());
            add_xml_string(grandchild, exml_names[exmlOBTYPE], ca_it->getObtype());
            add_xml_int(grandchild, exml_names[exmlATOMID], ca_it->getAtomid());

            double x, y, z;
            ca_it->getCoords(&x, &y, &z);

            const char *p    = gmx_ftoa(x).c_str();
            xmlNodePtr  baby = add_xml_child_val(grandchild, exml_names[exmlX], p);
            add_xml_string(baby, exml_names[exmlUNIT], ca_it->getUnit());
            p    = gmx_ftoa(y).c_str();
            baby = add_xml_child_val(grandchild, exml_names[exmlY], p);
            add_xml_string(baby, exml_names[exmlUNIT], ca_it->getUnit());
            p    = gmx_ftoa(z).c_str();
            baby = add_xml_child_val(grandchild, exml_names[exmlZ], p);
            add_xml_string(baby, exml_names[exmlUNIT], ca_it->getUnit());

            for (alexandria::AtomicChargeIterator q_it = ca_it->BeginQ();
                 (q_it < ca_it->EndQ()); q_it++)
            {
                p = gmx_ftoa(q_it->getQ()).c_str();
                xmlNodePtr atomptr = add_xml_child_val(grandchild, exml_names[exmlQ], p);
                add_xml_string(atomptr, exml_names[exmlTYPE], q_it->getType());
                add_xml_string(atomptr, exml_names[exmlUNIT], q_it->getUnit());
                add_xml_double(atomptr, exml_names[exmlTEMPERATURE], q_it->getTemperature());
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
        add_xml_string(child, exml_names[exmlCOMPNAME], mc_it->getCompName());
        for (alexandria::AtomNumIterator an_it = mc_it->BeginAtomNum();
             (an_it < mc_it->EndAtomNum()); an_it++)
        {
            xmlNodePtr grandchild = add_xml_child(child, exml_names[exmlCATOM]);
            add_xml_string(grandchild, exml_names[exmlC_NAME], an_it->getAtom());
            add_xml_int(grandchild, exml_names[exmlC_NUMBER], an_it->getNumber());
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
                    (int)mpt.size(), mp_it->getMolname().c_str());
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
