/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2019
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
 *             Paul J. van Maaren,
 *             David van der Spoel (Project leader)
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 */

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "poldata_xml.h"

#include <cstdlib>
#include <cstring>
#include <map>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"

#include "poldata.h"
#include "poldata_low.h"
#include "xml_util.h"

extern int xmlDoValidityCheckingDefaultValue;
namespace alexandria
{

const char *xmltypes[] = {
    nullptr,
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
    exmlGENTOP             = 0,
    exmlREFERENCE          = 1,
    exmlATOMTYPES          = 2,
    exmlATOMTYPE           = 3,
    exmlCHARGEMODEL        = 4,
    exmlPOLAR_UNIT         = 5,
    exmlCOMB_RULE          = 6,
    exmlNEXCL              = 7,
    exmlVERSION            = 8,
    exmlPOLTYPES           = 10,
    exmlPOLTYPE            = 11,
    exmlPTYPE              = 12,
    exmlELEM               = 13,
    exmlNAME               = 14,
    exmlDESC               = 15,
    exmlATYPE              = 16,
    exmlMILLER             = 17,
    exmlVALENCE            = 18,
    exmlBOSQUE             = 19,
    exmlBTYPE              = 20,
    exmlZTYPE              = 21,
    exmlNEIGHBORS          = 22,
    exmlAROMATIC           = 23,
    exmlGEOMETRY           = 24,
    exmlNUMBONDS           = 25,
    exmlPOLARIZABILITY     = 26,
    exmlSIGPOL             = 27,
    exmlVDWPARAMS          = 28,
    exmlEREF               = 29,
    exmlFUNCTION           = 30,
    exmlINTERACTION        = 31,
    exmlATOM1              = 32,
    exmlATOM2              = 33,
    exmlATOM3              = 34,
    exmlATOM4              = 35,
    exmlSIGMA              = 36,
    exmlBONDORDER          = 37,
    exmlPARAMS             = 38,
    exmlREFVALUE           = 39,
    exmlUNIT               = 40,
    exmlNTRAIN             = 41,
    exmlGT_VSITES          = 42,
    exmlGT_VSITE           = 43,
    exmlGT_BONDS           = 44,
    exmlGT_BOND            = 45,
    exmlGT_ANGLES          = 46,
    exmlGT_ANGLE           = 47,
    exmlGT_DIHEDRALS       = 48,
    exmlGT_DIHEDRAL        = 49,
    exmlGT_IMPROPERS       = 50,
    exmlGT_IMPROPER        = 51,
    exmlBSATOMS            = 52,
    exmlBSATOM             = 53,
    exmlMILATOMS           = 54,
    exmlTAU_UNIT           = 55,
    exmlAHP_UNIT           = 56,
    exmlMILATOM            = 57,
    exmlMILNAME            = 58,
    exmlALEXANDRIA_EQUIV   = 59,
    exmlATOMNUMBER         = 60,
    exmlTAU_AHC            = 61,
    exmlALPHA_AHP          = 62,
    exmlSYMMETRIC_CHARGES  = 63,
    exmlSYM_CHARGE         = 64,
    exmlCENTRAL            = 65,
    exmlATTACHED           = 66,
    exmlNUMATTACH          = 67,
    exmlEEMPROPS           = 68,
    exmlEEMPROP            = 69,
    exmlMODEL              = 70,
    exmlJ0                 = 71,
    exmlJ0_SIGMA           = 72,
    exmlCHI0               = 73,
    exmlCHI0_SIGMA         = 74,
    exmlZETA               = 75,
    exmlZETA_SIGMA         = 76,
    exmlROW                = 77,
    exmlCHARGES            = 80,
    exmlANGLE_UNIT         = 81,
    exmlLENGTH_UNIT        = 82,
    exmlDISTANCE           = 83,
    exmlNCONTROLATOMS      = 84,
    exmlNUMBER             = 85,
    exmlVTYPE              = 86,
    exmlANGLE              = 87,
    exmlFIXED              = 88,
    exmlNR                 = 89
};

std::map<const std::string, int> xmlxxx =
{
    { "gentop",                 exmlGENTOP           },
    { "reference",              exmlREFERENCE        },
    { "atomtypes",              exmlATOMTYPES        },
    { "atomtype",               exmlATOMTYPE         },
    { "chargemodel",            exmlCHARGEMODEL      },
    { "version",                exmlVERSION          },
    { "polarizability_unit",    exmlPOLAR_UNIT       },
    { "combination_rule",       exmlCOMB_RULE        },
    { "nexclusions",            exmlNEXCL            },
    { "poltypes",               exmlPOLTYPES         },
    { "poltype",                exmlPOLTYPE          },
    { "ptype",                  exmlPTYPE            },
    { "elem",                   exmlELEM             },
    { "name",                   exmlNAME             },
    { "description",            exmlDESC             },
    { "atype",                  exmlATYPE            },
    { "miller",                 exmlMILLER           },
    { "valence",                exmlVALENCE          },
    { "bosque",                 exmlBOSQUE           },
    { "btype",                  exmlBTYPE            },
    { "ztype",                  exmlZTYPE            },
    { "neighbors",              exmlNEIGHBORS        },
    { "aromatic",               exmlAROMATIC         },
    { "geometry",               exmlGEOMETRY         },
    { "numbonds",               exmlNUMBONDS         },
    { "polarizability",         exmlPOLARIZABILITY   },
    { "sigma_pol",              exmlSIGPOL           },
    { "fixed",                  exmlFIXED            },
    { "vdwparams",              exmlVDWPARAMS        },
    { "ref_enthalpy",           exmlEREF             },
    { "function",               exmlFUNCTION         },
    { "interaction",            exmlINTERACTION      },
    { "atom1",                  exmlATOM1            },
    { "atom2",                  exmlATOM2            },
    { "atom3",                  exmlATOM3            },
    { "atom4",                  exmlATOM4            },
    { "sigma",                  exmlSIGMA            },
    { "bondorder",              exmlBONDORDER        },
    { "params",                 exmlPARAMS           },
    { "refValue",               exmlREFVALUE         },
    { "unit",                   exmlUNIT             },
    { "ntrain",                 exmlNTRAIN           },
    { "gt_vsites",              exmlGT_VSITES        },
    { "gt_vsite",               exmlGT_VSITE         },
    { "gt_bonds",               exmlGT_BONDS         },
    { "gt_bond",                exmlGT_BOND          },
    { "gt_angles",              exmlGT_ANGLES        },
    { "gt_angle",               exmlGT_ANGLE         },
    { "gt_dihedrals",           exmlGT_DIHEDRALS     },
    { "gt_dihedral",            exmlGT_DIHEDRAL      },
    { "gt_impropers",           exmlGT_IMPROPERS     },
    { "gt_improper",            exmlGT_IMPROPER      },
    { "bsatoms",                exmlBSATOMS          },
    { "bsatom",                 exmlBSATOM           },
    { "milatoms",               exmlMILATOMS         },
    { "tau_ahc_unit",           exmlTAU_UNIT         },
    { "alpha_ahp_unit",         exmlAHP_UNIT         },
    { "milatom",                exmlMILATOM          },
    { "milname",                exmlMILNAME          },
    { "alexandria_equiv",       exmlALEXANDRIA_EQUIV },
    { "atomnumber",             exmlATOMNUMBER       },
    { "tau_ahc",                exmlTAU_AHC          },
    { "alpha_ahp",              exmlALPHA_AHP        },
    { "symmetric_charges",      exmlSYMMETRIC_CHARGES},
    { "sym_charge",             exmlSYM_CHARGE       },
    { "central",                exmlCENTRAL          },
    { "attached",               exmlATTACHED         },
    { "numattach",              exmlNUMATTACH        },
    { "eemprops",               exmlEEMPROPS         },
    { "eemprop",                exmlEEMPROP          },
    { "model",                  exmlMODEL            },
    { "jaa",                    exmlJ0               },
    { "jaa_sigma",              exmlJ0_SIGMA         },
    { "chi",                    exmlCHI0             },
    { "chi_sigma",              exmlCHI0_SIGMA       },
    { "zeta",                   exmlZETA             },
    { "zeta_sigma",             exmlZETA_SIGMA       },
    { "row",                    exmlROW              },
    { "charge",                 exmlCHARGES          },
    { "angle_unit",             exmlANGLE_UNIT       },
    { "length_unit",            exmlLENGTH_UNIT      },
    { "distance",               exmlDISTANCE         },
    { "ncontrolatoms",          exmlNCONTROLATOMS    },
    { "number",                 exmlNUMBER           },
    { "vtype",                  exmlVTYPE            },
    { "angle",                  exmlANGLE            }
};

std::map<int, const std::string> rmap;

static const char *exml_names(int xml)
{
    if (rmap.empty())
    {
        for (auto iter = xmlxxx.begin(); iter != xmlxxx.end(); ++iter)
        {
            rmap.insert({iter->second, iter->first});
        }
    }
    auto iter = rmap.find(xml);
    if (iter != rmap.end())
    {
        return iter->second.c_str();
    }
    return nullptr;
}

static bool NN(const std::string &x)
{
    return (0 != x.size());
}

static bool NNobligatory(const std::string xbuf[], int xml)
{
    if (0 == xbuf[xml].size())
    {
        gmx_fatal(FARGS, "Missing required variable '%s'", exml_names(xml));
    }
    return true;
}

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

static double xbuf_atof(const std::string xbuf[],
                        int               xbuf_index)
{
    return my_atof(xbuf[xbuf_index].c_str(), rmap[xbuf_index].c_str());
}

static void processAttr(FILE *fp, xmlAttrPtr attr, int elem,
                        int indent, Poldata &pd)
{
    std::string attrname, attrval;
    char        buf[100];
    int         i;
    std::string xbuf[exmlNR];

    for (i = 0; (i < exmlNR); i++)
    {
        xbuf[i] = "";
    }
    while (attr != nullptr)
    {
        attrname.assign((char *)attr->name);
        attrval.assign((char *)attr->children->content);

#define atest(s) ((strcasecmp(attrname, s) == 0) && (attrval != nullptr))
        auto iter = xmlxxx.find(attrname);
        if (iter != xmlxxx.end())
        {
            int kkk = iter->second;
            if (attrval.size() != 0)
            {
                xbuf[kkk] = attrval;
            }

            if (nullptr != fp)
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
    bool fixed = false;
    if (NN(xbuf[exmlFIXED]))
    {
        fixed = (xbuf[exmlFIXED] == "true" ||
                 xbuf[exmlFIXED] == "yes");
    }
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
            if (NNobligatory(xbuf, exmlCHARGEMODEL))
            {
                pd.setChargeModel(xbuf[exmlCHARGEMODEL]);
            }
            if (NNobligatory(xbuf, exmlVERSION))
            {
                pd.setVersion(xbuf[exmlVERSION]);
            }
            if (NNobligatory(xbuf, exmlFUNCTION))
            {
                pd.setVdwFunction(xbuf[exmlFUNCTION]);
            }
            if (NNobligatory(xbuf, exmlCOMB_RULE))
            {
                pd.setCombinationRule(xbuf[exmlCOMB_RULE]);
            }
            if (NNobligatory(xbuf, exmlNEXCL))
            {
                pd.setNexcl(atoi(xbuf[exmlNEXCL].c_str()));
            }
            break;
        case exmlBSATOMS:
            if (NN(xbuf[exmlPOLAR_UNIT]) &&
                NN(xbuf[exmlREFERENCE]))
            {
                pd.setBosqueFlags(xbuf[exmlPOLAR_UNIT], xbuf[exmlREFERENCE]);
            }
            break;
        case exmlGT_VSITES:
            if (NN(xbuf[exmlANGLE_UNIT]) && NN(xbuf[exmlLENGTH_UNIT]))
            {
                pd.setVsite_angle_unit(xbuf[exmlANGLE_UNIT]);
                pd.setVsite_length_unit(xbuf[exmlLENGTH_UNIT]);
            }
            break;
        case exmlGT_BONDS:
            if (NN(xbuf[exmlINTERACTION]) &&
                NN(xbuf[exmlFUNCTION])    &&
                NN(xbuf[exmlUNIT]))
            {
                ListedForces bonds (xbuf[exmlINTERACTION],
                                    xbuf[exmlFUNCTION],
                                    xbuf[exmlUNIT]);
                pd.addForces(bonds);
            }
            break;
        case exmlGT_ANGLES:
            if (NN(xbuf[exmlINTERACTION]) &&
                NN(xbuf[exmlFUNCTION])    &&
                NN(xbuf[exmlUNIT]))
            {
                ListedForces angles (xbuf[exmlINTERACTION],
                                     xbuf[exmlFUNCTION],
                                     xbuf[exmlUNIT]);
                pd.addForces(angles);
            }
            break;
        case exmlGT_DIHEDRALS:
            if (NN(xbuf[exmlINTERACTION]) &&
                NN(xbuf[exmlFUNCTION])    &&
                NN(xbuf[exmlUNIT]))
            {
                ListedForces dihedrals (xbuf[exmlINTERACTION],
                                        xbuf[exmlFUNCTION],
                                        xbuf[exmlUNIT]);
                pd.addForces(dihedrals);
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
                            xbuf_atof(xbuf, exmlPOLARIZABILITY),
                            xbuf_atof(xbuf, exmlSIGPOL));
            }
            break;
        case exmlATOMTYPE:
            if (NN(xbuf[exmlELEM]) &&
                NN(xbuf[exmlATYPE]) &&
                NN(xbuf[exmlBTYPE]) &&
                NN(xbuf[exmlPTYPE]) &&
                NN(xbuf[exmlZTYPE]) &&
                NN(xbuf[exmlVDWPARAMS]) &&
                NN(xbuf[exmlEREF]))
            {
                pd.addAtype(xbuf[exmlELEM],
                            xbuf[exmlDESC],
                            xbuf[exmlATYPE],
                            xbuf[exmlPTYPE],
                            xbuf[exmlBTYPE],
                            xbuf[exmlZTYPE],
                            fixed,
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
                             xbuf_atof(xbuf, exmlTAU_AHC),
                             xbuf_atof(xbuf, exmlALPHA_AHP),
                             xbuf[exmlALEXANDRIA_EQUIV]);
            }
            break;
        case exmlBSATOM:
            if (NN(xbuf[exmlELEM]) && NN(xbuf[exmlPOLARIZABILITY]))
            {
                pd.addBosque(xbuf[exmlELEM],
                             xbuf_atof(xbuf, exmlPOLARIZABILITY));
            }
            break;
        case exmlGT_VSITE:
            if (NN(xbuf[exmlATYPE])  && NN(xbuf[exmlVTYPE])    &&
                NN(xbuf[exmlNUMBER]) && NN(xbuf[exmlDISTANCE]) &&
                NN(xbuf[exmlANGLE])  && NN(xbuf[exmlNCONTROLATOMS]))
            {
                pd.addVsite(xbuf[exmlATYPE],
                            xbuf[exmlVTYPE],
                            atoi(xbuf[exmlNUMBER].c_str()),
                            atof(xbuf[exmlDISTANCE].c_str()),
                            atof(xbuf[exmlANGLE].c_str()),
                            atoi(xbuf[exmlNCONTROLATOMS].c_str()));
            }
            break;
        case exmlGT_BOND:
            if (NN(xbuf[exmlATOM1])    && NN(xbuf[exmlATOM2]) &&
                NN(xbuf[exmlREFVALUE]) && NN(xbuf[exmlSIGMA]) &&
                NN(xbuf[exmlNTRAIN]))
            {

                const std::vector<std::string> &atoms = {
                    xbuf[exmlATOM1].c_str(),
                    xbuf[exmlATOM2].c_str()
                };

                auto &force = pd.lastForces();

                force.addForce(atoms, xbuf[exmlPARAMS].c_str(),
                               fixed,
                               xbuf_atof(xbuf, exmlREFVALUE),
                               xbuf_atof(xbuf, exmlSIGMA),
                               atoi(xbuf[exmlNTRAIN].c_str()),
                               atoi(xbuf[exmlBONDORDER].c_str()));
            }
            break;
        case exmlGT_ANGLE:
            if (NN(xbuf[exmlATOM1]) && NN(xbuf[exmlATOM2])    &&
                NN(xbuf[exmlATOM3]) && NN(xbuf[exmlREFVALUE]) &&
                NN(xbuf[exmlSIGMA]) && NN(xbuf[exmlNTRAIN]))
            {
                const std::vector<std::string> &atoms = {
                    xbuf[exmlATOM1].c_str(),
                    xbuf[exmlATOM2].c_str(),
                    xbuf[exmlATOM3].c_str()
                };

                auto &force  = pd.lastForces();

                force.addForce(atoms, xbuf[exmlPARAMS].c_str(),
                               fixed,
                               xbuf_atof(xbuf, exmlREFVALUE),
                               xbuf_atof(xbuf, exmlSIGMA),
                               atoi(xbuf[exmlNTRAIN].c_str()));
            }
            break;
        case exmlGT_DIHEDRAL:
            if (NN(xbuf[exmlATOM1])    && NN(xbuf[exmlATOM2]) &&
                NN(xbuf[exmlATOM3])    && NN(xbuf[exmlATOM4]) &&
                NN(xbuf[exmlREFVALUE]) && NN(xbuf[exmlSIGMA]) &&
                NN(xbuf[exmlNTRAIN]))
            {
                const std::vector<std::string> &atoms = {
                    xbuf[exmlATOM1].c_str(),
                    xbuf[exmlATOM2].c_str(),
                    xbuf[exmlATOM3].c_str(),
                    xbuf[exmlATOM4].c_str()
                };

                auto &force = pd.lastForces();

                force.addForce(atoms, xbuf[exmlPARAMS].c_str(),
                               fixed,
                               xbuf_atof(xbuf, exmlREFVALUE),
                               xbuf_atof(xbuf, exmlSIGMA),
                               atoi(xbuf[exmlNTRAIN].c_str()));
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
        case exmlEEMPROPS:
            if (NN(xbuf[exmlREFERENCE]))
            {
                pd.setEpref(xbuf[exmlREFERENCE]);
            }
            break;
        case exmlEEMPROP:
            if (NN(xbuf[exmlNAME])       &&
                NN(xbuf[exmlJ0])      && NN(xbuf[exmlJ0_SIGMA])   &&
                NN(xbuf[exmlCHI0])    && NN(xbuf[exmlCHI0_SIGMA]) &&
                NN(xbuf[exmlZETA])    && NN(xbuf[exmlZETA_SIGMA]) &&
                NN(xbuf[exmlCHARGES]) && NN(xbuf[exmlROW]))
            {
                Eemprops eep(xbuf[exmlNAME],
                             xbuf[exmlROW],
                             xbuf[exmlZETA],
                             xbuf[exmlZETA_SIGMA],
                             xbuf[exmlCHARGES],
                             xbuf_atof(xbuf, exmlJ0),
                             xbuf_atof(xbuf, exmlJ0_SIGMA),
                             xbuf_atof(xbuf, exmlCHI0),
                             xbuf_atof(xbuf, exmlCHI0_SIGMA));
                pd.addEemprops(std::move(eep));
            }
            break;
        default:
            if (nullptr != debug)
            {
                fprintf(debug, "Unknown combination of attributes:\n");
                for (i = 0; (i < exmlNR); i++)
                {
                    if (xbuf[i].size() != 0)
                    {
                        fprintf(debug, "%s = %s\n", exml_names(i), xbuf[i].c_str());
                    }
                }
            }
    }
}

static void processTree(FILE *fp, xmlNodePtr tree, int indent,
                        Poldata &pd, gmx_atomprop_t aps)
{
    char buf[100];

    while (tree != nullptr)
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
            if (fp)
            {
                sp(indent, buf, 99);
                fprintf(fp, "%sElement node name %s\n", buf, (char *)tree->name);
            }
            auto iter = xmlxxx.find((const char *)tree->name);
            if (iter != xmlxxx.end())
            {
                int elem = iter->second;
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
    std::string fn, fn2;

    fn = fileName;
    if (fn.empty())
    {
        fn.assign("ACM-g_2020.dat");
    }
    fn2 = gmx::findLibraryFile(fn, true, false);
    if (fn2.empty())
    {
        fn = "alexandria.ff/" + fn;
        fn2 = gmx::findLibraryFile(fn, true, false);
        if (fn2.empty())
        {
            fn = "top/" + fn;
            fn2 = gmx::findLibraryFile(fn, true, false);
        }
    }
    if (fn2.empty())
    {
        gmx_fatal(FARGS, "Could not find %s\n", fn.c_str());
    }
    if (nullptr != debug)
    {
        fprintf(debug, "Opening library file %s\n", fn2.c_str());
    }
    printf("Opening library file %s\n", fn2.c_str());
    xmlDoValidityCheckingDefaultValue = 0;
    doc = xmlParseFile(fn2.c_str());
    if (doc == nullptr)
    {
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "Error reading XML file '%s'. Run a syntax checker such as nsgmls.",
                 fn2.c_str());
        GMX_THROW(gmx::FileIOError(buf));
    }

    pd.setFilename(fn2);
    processTree(debug, doc->children, 0, pd, aps);

    xmlFreeDoc(doc);

    // Generate maps
    pd.makeMappings();
    pd.checkConsistency(debug);
    if (nullptr != debug)
    {
        writePoldata("pdout.dat", &pd, false);
    }
}

static void addXmlPoldata(xmlNodePtr parent, const Poldata *pd)
{
    xmlNodePtr   child, grandchild;
    int          nexcl;
    std::string  geometry, name,
                 acentral, attached, tau_unit, ahp_unit,
                 epref, desc, params;
    std::string  neighbors, zeta, qstr, rowstr;
    std::string  tmp, func, blu;

    child = add_xml_child(parent, exml_names(exmlATOMTYPES));
    add_xml_char(child, exml_names(exmlCHARGEMODEL),
                 getEemtypeName(pd->getChargeModel()));
    tmp   = pd->getVersion();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names(exmlVERSION), tmp.c_str());
    }
    tmp = pd->getVdwFunction();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names(exmlFUNCTION), tmp.c_str());
    }
    tmp = pd->getCombinationRule();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names(exmlCOMB_RULE), tmp.c_str());
    }
    nexcl = pd->getNexcl();
    add_xml_int(child, exml_names(exmlNEXCL), nexcl);

    for (auto aType = pd->getAtypeBegin();
         aType != pd->getAtypeEnd(); aType++)
    {
        grandchild = add_xml_child(child, exml_names(exmlATOMTYPE));
        add_xml_char(grandchild, exml_names(exmlELEM), aType->getElem().c_str());
        add_xml_char(grandchild, exml_names(exmlDESC), aType->getDesc().c_str());
        add_xml_char(grandchild, exml_names(exmlATYPE), aType->getType().c_str());
        add_xml_char(grandchild, exml_names(exmlPTYPE), aType->getPtype().c_str());
        add_xml_char(grandchild, exml_names(exmlBTYPE), aType->getBtype().c_str());
        add_xml_char(grandchild, exml_names(exmlZTYPE), aType->getZtype().c_str());
        add_xml_char(grandchild, exml_names(exmlVDWPARAMS), aType->getVdwparams().c_str());
        if (aType->fixed())
        {
            add_xml_char(grandchild, exml_names(exmlFIXED), "true");
        }
        add_xml_char(grandchild, exml_names(exmlEREF), aType->getRefEnthalpy().c_str());
    }
    child = add_xml_child(parent, exml_names(exmlPOLTYPES));
    tmp   = pd->getPolarUnit();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names(exmlPOLAR_UNIT), tmp.c_str());
    }
    tmp = pd->getPolarRef();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names(exmlREFERENCE), tmp.c_str());
    }
    for (auto pType = pd->getPtypeBegin();
         pType != pd->getPtypeEnd(); pType++)
    {
        grandchild = add_xml_child(child, exml_names(exmlPOLTYPE));
        add_xml_char(grandchild, exml_names(exmlPTYPE), pType->getType().c_str());
        add_xml_char(grandchild, exml_names(exmlMILLER), pType->getMiller().c_str());
        add_xml_char(grandchild, exml_names(exmlBOSQUE), pType->getBosque().c_str());
        add_xml_double(grandchild, exml_names(exmlPOLARIZABILITY), pType->getPolarizability());
        add_xml_double(grandchild, exml_names(exmlSIGPOL), pType->getSigPol());
    }
    tmp   = pd->getVsite_angle_unit();
    if (0 != tmp.size())
    {
        child = add_xml_child(parent, exml_names(exmlGT_VSITES));
        add_xml_char(child, exml_names(exmlANGLE_UNIT), tmp.c_str());
    }
    tmp   = pd->getVsite_length_unit();
    if (0 != tmp.size())
    {
        add_xml_char(child, exml_names(exmlLENGTH_UNIT), tmp.c_str());
    }
    for (auto vsite = pd->getVsiteBegin(); vsite != pd->getVsiteEnd(); vsite++)
    {
        grandchild = add_xml_child(child, exml_names(exmlGT_VSITE));
        add_xml_char(grandchild, exml_names(exmlATYPE), vsite->atype().c_str());
        add_xml_char(grandchild, exml_names(exmlVTYPE), vsiteType2string(vsite->type()));
        add_xml_int(grandchild, exml_names(exmlNUMBER), vsite->nvsite());
        add_xml_double(grandchild, exml_names(exmlDISTANCE), vsite->distance());
        add_xml_double(grandchild, exml_names(exmlANGLE), vsite->angle());
        add_xml_int(grandchild, exml_names(exmlNCONTROLATOMS), vsite->ncontrolatoms());
    }
    for (auto fs = pd->forcesBegin(); fs != pd->forcesEnd(); fs++)
    {
        if (eitBONDS == fs->iType())
        {
            child = add_xml_child(parent, exml_names(exmlGT_BONDS));
            blu   = fs->unit();
            if (blu.size() != 0)
            {
                add_xml_char(child, exml_names(exmlUNIT), blu.c_str());
            }

            add_xml_char(child, exml_names(exmlINTERACTION), iType2string(fs->iType()));

            func = fs->function();
            if (func.size() != 0)
            {
                add_xml_char(child, exml_names(exmlFUNCTION), func.c_str());
            }
            for (auto f = fs->forceBegin(); f != fs->forceEnd(); f++)
            {
                const std::vector<std::string> atoms = f->atoms();

                grandchild = add_xml_child(child, exml_names(exmlGT_BOND));
                add_xml_char(grandchild, exml_names(exmlATOM1), atoms[0].c_str());
                add_xml_char(grandchild, exml_names(exmlATOM2), atoms[1].c_str());
                add_xml_int(grandchild, exml_names(exmlBONDORDER), f->bondOrder());
                add_xml_double(grandchild, exml_names(exmlREFVALUE), f->refValue());
                add_xml_double(grandchild, exml_names(exmlSIGMA), f->sigma());
                add_xml_int(grandchild, exml_names(exmlNTRAIN), f->ntrain());
                if (f->fixed())
                {
                    add_xml_char(grandchild, exml_names(exmlFIXED), "true");
                }
                add_xml_char(grandchild, exml_names(exmlPARAMS), f->params().c_str());
            }
        }
        else if (eitANGLES == fs->iType() ||
                 eitLINEAR_ANGLES == fs->iType())
        {
            child = add_xml_child(parent, exml_names(exmlGT_ANGLES));
            blu   = fs->unit();
            if (blu.size() != 0)
            {
                add_xml_char(child, exml_names(exmlUNIT), blu.c_str());
            }

            add_xml_char(child, exml_names(exmlINTERACTION), iType2string(fs->iType()));

            func = fs->function();
            if (func.size() != 0)
            {
                add_xml_char(child, exml_names(exmlFUNCTION), func.c_str());
            }
            for (auto f = fs->forceBegin(); f != fs->forceEnd(); f++)
            {
                const std::vector<std::string> atoms = f->atoms();

                grandchild = add_xml_child(child, exml_names(exmlGT_ANGLE));
                add_xml_char(grandchild, exml_names(exmlATOM1), atoms[0].c_str());
                add_xml_char(grandchild, exml_names(exmlATOM2), atoms[1].c_str());
                add_xml_char(grandchild, exml_names(exmlATOM3), atoms[2].c_str());
                add_xml_double(grandchild, exml_names(exmlREFVALUE), f->refValue());
                add_xml_double(grandchild, exml_names(exmlSIGMA), f->sigma());
                add_xml_int(grandchild, exml_names(exmlNTRAIN), f->ntrain());
                if (f->fixed())
                {
                    add_xml_char(grandchild, exml_names(exmlFIXED), "true");
                }
                add_xml_char(grandchild, exml_names(exmlPARAMS), f->params().c_str());
            }
        }
        else if (eitPROPER_DIHEDRALS   == fs->iType() ||
                 eitIMPROPER_DIHEDRALS == fs->iType())
        {
            child = add_xml_child(parent, exml_names(exmlGT_DIHEDRALS));
            blu   = fs->unit();
            if (blu.size() != 0)
            {
                add_xml_char(child, exml_names(exmlUNIT), blu.c_str());
            }

            add_xml_char(child, exml_names(exmlINTERACTION), iType2string(fs->iType()));

            func = fs->function();
            if (func.size() != 0)
            {
                add_xml_char(child, exml_names(exmlFUNCTION), func.c_str());
            }

            for (auto f = fs->forceBegin(); f != fs->forceEnd(); f++)
            {
                const std::vector<std::string> atoms = f->atoms();

                grandchild = add_xml_child(child, exml_names(exmlGT_DIHEDRAL));
                add_xml_char(grandchild, exml_names(exmlATOM1), atoms[0].c_str());
                add_xml_char(grandchild, exml_names(exmlATOM2), atoms[1].c_str());
                add_xml_char(grandchild, exml_names(exmlATOM3), atoms[2].c_str());
                add_xml_char(grandchild, exml_names(exmlATOM4), atoms[3].c_str());
                add_xml_double(grandchild, exml_names(exmlREFVALUE), f->refValue());
                add_xml_double(grandchild, exml_names(exmlSIGMA), f->sigma());
                add_xml_int(grandchild, exml_names(exmlNTRAIN), f->ntrain());
                if (f->fixed())
                {
                    add_xml_char(grandchild, exml_names(exmlFIXED), "true");
                }
                add_xml_char(grandchild, exml_names(exmlPARAMS), f->params().c_str());
            }
        }
    }

    child = add_xml_child(parent, exml_names(exmlBSATOMS));
    std::string ref;
    pd->getBosqueFlags(tmp, ref);
    add_xml_char(child, exml_names(exmlPOLAR_UNIT), tmp.c_str());
    add_xml_char(child, exml_names(exmlREFERENCE), ref.c_str());

    for (auto bosque = pd->getBosqueBegin();
         bosque != pd->getBosqueEnd(); bosque++)
    {
        grandchild = add_xml_child(child, exml_names(exmlBSATOM));
        add_xml_char(grandchild, exml_names(exmlELEM), bosque->getBosque().c_str());
        add_xml_double(grandchild, exml_names(exmlPOLARIZABILITY), bosque->getPolarizability());
    }
    child = add_xml_child(parent, exml_names(exmlMILATOMS));
    std::string milref;
    pd->getMillerFlags(tau_unit, ahp_unit, milref);
    add_xml_char(child, exml_names(exmlTAU_UNIT), tau_unit.c_str());
    add_xml_char(child, exml_names(exmlAHP_UNIT), ahp_unit.c_str());
    add_xml_char(child, exml_names(exmlREFERENCE), milref.c_str());
    for (auto miller = pd->getMillerBegin();
         miller != pd->getMillerEnd(); miller++)
    {
        grandchild = add_xml_child(child, exml_names(exmlMILATOM));
        add_xml_char(grandchild, exml_names(exmlMILNAME), miller->getMiller().c_str());
        add_xml_int(grandchild, exml_names(exmlATOMNUMBER), miller->getAtomnumber());
        add_xml_double(grandchild, exml_names(exmlTAU_AHC), miller->getTauAhc());
        add_xml_double(grandchild, exml_names(exmlALPHA_AHP), miller->getAlphaAhp());
        const std::string ae = miller->getAlexandriaEquiv();
        if (ae.size() > 0)
        {
            add_xml_char(grandchild, exml_names(exmlALEXANDRIA_EQUIV), ae.c_str());
        }
    }

    child = add_xml_child(parent, exml_names(exmlSYMMETRIC_CHARGES));
    for (auto symcharges = pd->getSymchargesBegin();
         symcharges != pd->getSymchargesEnd(); symcharges++)
    {
        grandchild = add_xml_child(child, exml_names(exmlSYM_CHARGE));
        add_xml_char(grandchild, exml_names(exmlCENTRAL), symcharges->getCentral().c_str());
        add_xml_char(grandchild, exml_names(exmlATTACHED), symcharges->getAttached().c_str());
        add_xml_int(grandchild, exml_names(exmlNUMATTACH), symcharges->getNumattach());
    }

    child = add_xml_child(parent, exml_names(exmlEEMPROPS));
    add_xml_char(child, exml_names(exmlREFERENCE), pd->getEpref().c_str());
    for (auto eep = pd->BeginEemprops();
         eep != pd->EndEemprops(); eep++)
    {
        grandchild = add_xml_child(child, exml_names(exmlEEMPROP));
        add_xml_char(grandchild, exml_names(exmlNAME), eep->getName());
        add_xml_double(grandchild, exml_names(exmlJ0), eep->getJ0());
        add_xml_double(grandchild, exml_names(exmlJ0_SIGMA), eep->getJ0_sigma());
        add_xml_double(grandchild, exml_names(exmlCHI0), eep->getChi0());
        add_xml_double(grandchild, exml_names(exmlCHI0_SIGMA), eep->getChi0_sigma());
        add_xml_char(grandchild, exml_names(exmlZETA), eep->getZetastr());
        add_xml_char(grandchild, exml_names(exmlZETA_SIGMA), eep->getZeta_sigma());
        add_xml_char(grandchild, exml_names(exmlCHARGES), eep->getQstr());
        add_xml_char(grandchild, exml_names(exmlROW), eep->getRowstr());
    }
}

void writePoldata(const std::string &fileName,
                  const Poldata     *pd,
                  bool               compress)
{
    xmlDocPtr   doc;
    xmlDtdPtr   dtd;
    xmlNodePtr  myroot;
    xmlChar    *libdtdname, *dtdname, *gmx;

    gmx        = (xmlChar *) "gentop";
    dtdname    = (xmlChar *) "gentop.dtd";
    libdtdname = dtdname;

    if ((doc = xmlNewDoc((xmlChar *)"1.0")) == nullptr)
    {
        gmx_fatal(FARGS, "Creating XML document %s", fileName.c_str());
    }

    if ((dtd = xmlCreateIntSubset(doc, dtdname, libdtdname, dtdname)) == nullptr)
    {
        gmx_fatal(FARGS, "Creating XML DTD in %s", fileName.c_str());
    }

    if ((myroot = xmlNewDocNode(doc, nullptr, gmx, nullptr)) == nullptr)
    {
        gmx_fatal(FARGS, "Creating root element for %s", fileName.c_str());
    }
    dtd->next    = myroot;
    myroot->prev = (xmlNodePtr) dtd;

    /* Add molecule definitions */
    addXmlPoldata(myroot, pd);

    xmlSetDocCompressMode(doc, compress ? 1 : 0);
    xmlIndentTreeOutput = 1;
    if (xmlSaveFormatFileEnc(fileName.c_str(), doc, "ISO-8859-1", 2) == 0)
    {
        gmx_fatal(FARGS, "Saving file %s", fileName.c_str());
    }
    xmlFreeDoc(doc);
}

}
