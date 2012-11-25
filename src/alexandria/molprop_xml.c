/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: molprop_xml.c,v 1.23 2009/06/01 06:13:18 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "gmx_fatal.h"
#include "futil.h"
#include "macros.h"
#include "grompp.h"
#include "smalloc.h"
#include "string2.h"
#include "xml_util.h"
#include "molprop.h"
#include "molprop_xml.h"

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
    exmlNAME, exmlVALUE, exmlERROR, 
    exmlMETHOD, exmlREFERENCE, exmlTYPE, exmlSOURCE,
    exmlBOND, exmlAI, exmlAJ, exmlBONDORDER,
    exmlCOMPOSITION, exmlCOMPNAME, exmlCATOM, exmlC_NAME, exmlC_NUMBER,
    exmlCALCULATION, exmlPROGRAM, exmlBASISSET, exmlCONFORMATION,
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
    "calculation", "program", "basisset", "conformation",
    "unit", "atom", "atomid", "obtype", "coord_unit", "potential_unit", "espid",
    "x", "y", "z", "V", "xx", "yy", "zz", "xy", "xz", "yz", "q"
};

typedef struct {
    int nmolprop;
    gmx_molprop_t *mpt;
} t_xmlrec;

static char *sp(int n, char buf[], int maxindent)
{
    int i;
    if (n>=maxindent)
    {
        n=maxindent-1;
    }
    
    /* Don't indent more than maxindent characters */
    for (i=0; (i<n); i++)
    {
        buf[i] = ' ';
    }
    buf[i] = '\0';
  
    return buf;
}

double my_atof(char *ptr)
{
    if (NULL != ptr) 
        return atof(ptr);
    else
        return 0;
}

static void get_attributes(FILE *fp,gmx_bool bZero,int indent,xmlAttrPtr attr,char *xbuf[])
{
    char buf[100];
    char *attrname,*attrval;
    int  i,kkk;    
    
    if (bZero) 
    {
        for(i=0; (i<exmlNR); i++)
        {
            xbuf[i] = NULL;
        }
    }
    
    while (attr != NULL) 
    {
        attrname = (char *)attr->name;
        attrval  = (char *)attr->children->content;
        
#define atest(s) ((strcasecmp(attrname,s) == 0) && (attrval != NULL))
        if ((kkk = find_elem(attrname,exmlNR,exml_names)) != -1)
        {
            if (attrval != NULL)
                xbuf[kkk] = strdup(attrval);
        }
        if (fp)
            fprintf(fp,"%sProperty: '%s' Value: '%s'\n",sp(indent,buf,99),
                    attrname,attrval);
        attr = attr->next;
#undef atest
    }
}

static void process_children(xmlNodePtr tree,char *xbuf[])
{
    int node;
    
    while (NULL != tree) 
    {
        if (((node = find_elem((char *)tree->name,exmlNR,exml_names)) != -1) &&
            (NULL != tree->children) &&
            (NULL != tree->children->content))
        {
            xbuf[node] = strdup((char *)tree->children->content);
        }
        tree = tree->next;
    }
}

static void mp_process_tree(FILE *fp,xmlNodePtr tree,int parent,
                            int indent,t_xmlrec *xml,int expref)
{
    xmlNodePtr tc;
    char buf[100];
    gmx_molprop_t mpt;
    int  i,kkk,atomref;
    char *xbuf[exmlNR];
    int  node,elem=-1;
  
    while (tree != NULL) 
    {
        if (fp) 
        {
            if ((tree->type > 0) && (tree->type < NXMLTYPES))
                fprintf(fp,"Node type %s encountered with name %s\n",
                        xmltypes[tree->type],(char *)tree->name);
            else
                fprintf(fp,"Node type %d encountered\n",tree->type);
        }
        if (xml->nmolprop > 0)
            mpt = xml->mpt[xml->nmolprop-1];
        else
            mpt = NULL;
            
        switch (tree->type)
        {
        case XML_TEXT_NODE:
        {
            //   fprintf(stderr,"Text node %s encountered\n",(char *)tree->name);
            break;
        }
        case XML_ELEMENT_NODE:
        {
            if ((elem = find_elem((char *)tree->name,exmlNR,exml_names)) != -1)
            {
                if (fp)
                {
                    fprintf(fp,"%sElement node name %s\n",sp(indent,buf,99),
                            (char *)tree->name);
                }
                get_attributes(fp,TRUE,indent,tree->properties,xbuf);
                
                switch (elem) 
                {
                case exmlMOLECULES:
                    break;
                case exmlMOLECULE:
                    xml->nmolprop++;
                    srenew(xml->mpt,xml->nmolprop);
                    mpt = xml->mpt[xml->nmolprop-1] = gmx_molprop_init();
                    if (NN(xbuf[exmlFORMULA]))
                        gmx_molprop_set_formula(mpt,xbuf[exmlFORMULA]);
                    if (NN(xbuf[exmlMOLNAME]))
                        gmx_molprop_set_molname(mpt,xbuf[exmlMOLNAME]);
                    if (NN(xbuf[exmlMASS]))
                        gmx_molprop_set_mass(mpt,atof(xbuf[exmlMASS]));
                    if (NN(xbuf[exmlCHARGE]))
                        gmx_molprop_set_charge(mpt,atoi(xbuf[exmlCHARGE]));
                    if (NN(xbuf[exmlMULTIPLICITY]))
                        gmx_molprop_set_multiplicity(mpt,atoi(xbuf[exmlMULTIPLICITY]));
                    break;
                    /* The items below are handled when treating attributes */
                case exmlMOLINFO:
                    if (NN(xbuf[exmlIUPAC]))
                        gmx_molprop_set_iupac(mpt,xbuf[exmlIUPAC]);
                    if (NN(xbuf[exmlCAS]))
                        gmx_molprop_set_cas(mpt,xbuf[exmlCAS]);
                    if (NN(xbuf[exmlCID]))
                        gmx_molprop_set_cid(mpt,xbuf[exmlCID]);
                    if (NN(xbuf[exmlINCHI]))
                        gmx_molprop_set_inchi(mpt,xbuf[exmlINCHI]);
                    break;
                case exmlCATEGORY:
                    if (NN(xbuf[exmlCATNAME]))
                        gmx_molprop_add_category(mpt,xbuf[exmlCATNAME]);
                    break;
                case exmlPOLARIZABILITY:
                    process_children(tree->children,xbuf);
                    if (NN(xbuf[exmlTYPE])  && NN(xbuf[exmlUNIT]) &&
                        NN(xbuf[exmlVALUE]) && NN(xbuf[exmlERROR]))
                        gmx_molprop_add_polar(mpt,expref,xbuf[exmlTYPE],xbuf[exmlUNIT],
                                              my_atof(xbuf[exmlXX]),my_atof(xbuf[exmlYY]),
                                              my_atof(xbuf[exmlZZ]),
                                              my_atof(xbuf[exmlVALUE]),my_atof(xbuf[exmlERROR]));
                    break;
                case exmlPOTENTIAL:
                    process_children(tree->children,xbuf);
                    if (NN(xbuf[exmlX_UNIT]) && NN(xbuf[exmlV_UNIT]) &&
                        NN(xbuf[exmlESPID]) &&
                        NN(xbuf[exmlX]) && NN(xbuf[exmlY]) &&
                        NN(xbuf[exmlZ]) && NN(xbuf[exmlV]))
                        gmx_molprop_add_potential(mpt,expref,xbuf[exmlX_UNIT],
                                                  xbuf[exmlV_UNIT],
                                                  atoi(xbuf[exmlESPID]),
                                                  my_atof(xbuf[exmlX]),my_atof(xbuf[exmlY]),
                                                  my_atof(xbuf[exmlZ]),my_atof(xbuf[exmlV]));
                    break;
                case exmlDIPOLE: 
                    process_children(tree->children,xbuf);
                    if (NN(xbuf[exmlTYPE]) && NN(xbuf[exmlUNIT]) &&
                        NN(xbuf[exmlVALUE]) && NN(xbuf[exmlERROR]))
                        gmx_molprop_add_dipole(mpt,expref,xbuf[exmlTYPE],xbuf[exmlUNIT],
                                               my_atof(xbuf[exmlX]),my_atof(xbuf[exmlY]),
                                               my_atof(xbuf[exmlZ]),
                                               my_atof(xbuf[exmlVALUE]),my_atof(xbuf[exmlERROR]));
                    break;
                case exmlQUADRUPOLE: 
                    process_children(tree->children,xbuf);
                    if (NN(xbuf[exmlTYPE]) && NN(xbuf[exmlUNIT]) &&
                        NN(xbuf[exmlXX]) && NN(xbuf[exmlYY]) && NN(xbuf[exmlZZ]))
                         gmx_molprop_add_quadrupole(mpt,expref,xbuf[exmlTYPE],xbuf[exmlUNIT],
                                                    my_atof(xbuf[exmlXX]),my_atof(xbuf[exmlYY]),
                                                    my_atof(xbuf[exmlZZ]),my_atof(xbuf[exmlXY]),
                                                    my_atof(xbuf[exmlXZ]),my_atof(xbuf[exmlYZ]));
                    break;
                case exmlBOND: 
                    process_children(tree->children,xbuf);
                    if (NN(xbuf[exmlAI]) && NN(xbuf[exmlAJ]) &&
                        NN(xbuf[exmlBONDORDER]))
                        gmx_molprop_add_bond(mpt,atoi(xbuf[exmlAI]),
                                             atoi(xbuf[exmlAJ]),
                                             atoi(xbuf[exmlBONDORDER]));
                    break;
                case exmlENERGY:
                    process_children(tree,xbuf);
                    if (NN(xbuf[exmlTYPE]) && NN(xbuf[exmlUNIT]) &&
                        NN(xbuf[exmlENERGY]))
                        gmx_molprop_add_energy(mpt,expref,xbuf[exmlTYPE],xbuf[exmlUNIT],
                                               my_atof(xbuf[exmlENERGY]),
                                               xbuf[exmlERROR] ? my_atof(xbuf[exmlERROR]) : 0.0);
                    break;
                    
                case exmlCOMPOSITION: 
                    if (NN(xbuf[exmlCOMPNAME]))
                        gmx_molprop_add_composition(mpt,xbuf[exmlCOMPNAME]);
                    break;
                case exmlCATOM:
                    if (NN(xbuf[exmlC_NAME]) && NN(xbuf[exmlC_NUMBER]))
                        gmx_molprop_add_composition_atom(mpt,NULL,xbuf[exmlC_NAME],atoi(xbuf[exmlC_NUMBER]));
                    break;
                case exmlCALCULATION:
                    if (NN(xbuf[exmlPROGRAM]) && NN(xbuf[exmlMETHOD]) && 
                        NN(xbuf[exmlBASISSET]) && NN(xbuf[exmlREFERENCE]) &&
                        NN(xbuf[exmlCONFORMATION]))
                        gmx_molprop_add_calculation(mpt,xbuf[exmlPROGRAM],xbuf[exmlMETHOD],
                                                    xbuf[exmlBASISSET],xbuf[exmlREFERENCE],
                                                    xbuf[exmlCONFORMATION],&expref);
                    else 
                        gmx_fatal(FARGS,"Trying to add calculation with program %s, method %s, basisset %s and reference %s and conformation %s",
                                  xbuf[exmlPROGRAM],xbuf[exmlMETHOD],
                                  xbuf[exmlBASISSET],xbuf[exmlREFERENCE],
                                  xbuf[exmlCONFORMATION]);
                    break;
                case exmlATOM:
                    if (NN(xbuf[exmlNAME]) && NN(xbuf[exmlOBTYPE]) && NN(xbuf[exmlATOMID]))
                    {
                        gmx_molprop_calc_add_atom(mpt,expref,xbuf[exmlNAME],xbuf[exmlOBTYPE],
                                                  atoi(xbuf[exmlATOMID]),&atomref);
                        for(tc = tree->children; (NULL != tc); tc = tc->next)
                        {
                            get_attributes(fp,FALSE,indent,tc->properties,xbuf);
                            if (((node = find_elem((char *)tc->name,exmlNR,exml_names)) != -1) &&
                                (NULL != tc->children) &&
                                (NULL != tc->children->content))
                            {
                                xbuf[node] = strdup((char *)tc->children->content);
                            }
                            
                            if (NN(xbuf[exmlX]) && NN(xbuf[exmlY]) && NN(xbuf[exmlZ])
                                && NN(xbuf[exmlUNIT])) 
                            {
                                gmx_molprop_calc_set_atomcoords(mpt,expref,atomref,
                                                                xbuf[exmlUNIT],
                                                                my_atof(xbuf[exmlX]),
                                                                my_atof(xbuf[exmlY]),
                                                                my_atof(xbuf[exmlZ]));
                                sfree(xbuf[exmlX]); xbuf[exmlX] = NULL;
                                sfree(xbuf[exmlY]); xbuf[exmlY] = NULL;
                                sfree(xbuf[exmlZ]); xbuf[exmlZ] = NULL;
                                sfree(xbuf[exmlUNIT]); xbuf[exmlUNIT] = NULL;
                            }
                            if (NN(xbuf[exmlQ]) && NN(xbuf[exmlTYPE]) && NN(xbuf[exmlUNIT]))
                            {
                                gmx_molprop_calc_set_atomcharge(mpt,expref,atomref,
                                                                xbuf[exmlTYPE],
                                                                xbuf[exmlUNIT],
                                                                my_atof(xbuf[exmlQ]));
                                sfree(xbuf[exmlQ]);    xbuf[exmlQ]    = NULL;
                                sfree(xbuf[exmlUNIT]); xbuf[exmlUNIT] = NULL;
                                sfree(xbuf[exmlTYPE]); xbuf[exmlTYPE] = NULL;
                            }
                        }
                    }
                    break;
                    
                case exmlEXPERIMENT:
                    if (NN(xbuf[exmlREFERENCE]) && NN(xbuf[exmlCONFORMATION]))
                    {
                        gmx_molprop_add_experiment(mpt,xbuf[exmlREFERENCE],
                                                   xbuf[exmlCONFORMATION],&expref);
                    }
                    else
                        gmx_fatal(FARGS,"Experimental data without reference");
                    break;
                default:
                    break;
                }
            }
            for(i=0; (i<exmlNR); i++)
                if (NN(xbuf[i]))
                    sfree(xbuf[i]);
            if (tree->children) {
                if ((elem = find_elem((char *)tree->name,exmlNR,exml_names)) != -1)
                {
                    mp_process_tree(fp,tree->children,elem,indent+2,xml,expref);
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

gmx_molprop_t *gmx_molprops_read(const char *fn,int *nmolprop)
{
    xmlDocPtr     doc;
    int           i,npd;
    t_xmlrec      *xml;
    const char *db="alexandria.ff/molprops.dat";
  
    xmlDoValidityCheckingDefaultValue = 0;
    if (NULL == fn) 
        fn = gmxlibfn(db);
    if (NULL != debug)
        fprintf(debug,"Opening %s\n",fn);
    if ((doc = xmlParseFile(fn)) == NULL) 
    {
        fprintf(stderr,"Reading XML file %s. Run a syntax checker such as nsgmls.",
                fn);
        exit(1);
    }

    snew(xml,1);
    mp_process_tree(NULL,doc->children,0,0,xml,-12345);

    xmlFreeDoc(doc);
  
    *nmolprop = xml->nmolprop;
  
    return xml->mpt;
}

static void add_properties(xmlNodePtr exp,gmx_molprop_t mpt,int ref)
{
    xmlNodePtr child;
    char *type,*unit,*x_unit,*v_unit,*ptr;
    int espid;
    double value,error,x,y,z,V,xx,yy,zz,xy,xz,yz;
    
    while (gmx_molprop_get_energy(mpt,ref,&type,&unit,&value,&error) == 1) 
    {
        ptr = gmx_ftoa(value);
        child = add_xml_child_val(exp,exml_names[exmlENERGY],ptr);
        sfree(ptr);
        add_xml_char(child,exml_names[exmlTYPE],type);
        add_xml_char(child,exml_names[exmlUNIT],unit);
        sfree(type);
        sfree(unit);
    }
    while (gmx_molprop_get_potential(mpt,ref,&x_unit,&v_unit,
                                     &espid,&x,&y,&z,&V) == 1) 
    {
        child = add_xml_child(exp,exml_names[exmlPOTENTIAL]);
        add_xml_char(child,exml_names[exmlX_UNIT],x_unit);
        add_xml_char(child,exml_names[exmlV_UNIT],v_unit);
        add_xml_int(child,exml_names[exmlESPID],espid);
        if ((x != 0) || (y != 0) || (z != 0) || (V != 0)) 
        {
            ptr = gmx_ftoa(x);
            add_xml_child_val(child,exml_names[exmlX],ptr);
            sfree(ptr);
            ptr = gmx_ftoa(y);
            add_xml_child_val(child,exml_names[exmlY],ptr);
            sfree(ptr);
            ptr = gmx_ftoa(z);
            add_xml_child_val(child,exml_names[exmlZ],ptr);
            sfree(ptr);
            ptr = gmx_ftoa(V);
            add_xml_child_val(child,exml_names[exmlV],ptr);
            sfree(ptr);
        }
        sfree(x_unit);
        sfree(v_unit);
    }
    while (gmx_molprop_get_dipole(mpt,ref,&type,&unit,&x,&y,&z,&value,&error) == 1) 
    {
        child = add_xml_child(exp,exml_names[exmlDIPOLE]);
        add_xml_char(child,exml_names[exmlTYPE],type);
        add_xml_char(child,exml_names[exmlUNIT],unit);
        ptr = gmx_ftoa(value);
        add_xml_child_val(child,exml_names[exmlVALUE],ptr);
        sfree(ptr);
        ptr = gmx_ftoa(error);
        add_xml_child_val(child,exml_names[exmlERROR],ptr);
        sfree(ptr);
        if ((x != 0) || (y != 0) || (z != 0)) 
        {
            ptr = gmx_ftoa(x);
            add_xml_child_val(child,exml_names[exmlX],ptr);
            sfree(ptr);
            ptr = gmx_ftoa(y);
            add_xml_child_val(child,exml_names[exmlY],ptr);
            sfree(ptr);
            ptr = gmx_ftoa(z);
            add_xml_child_val(child,exml_names[exmlZ],ptr);
            sfree(ptr);
        }   
        sfree(type);
        sfree(unit);
    }
    while (gmx_molprop_get_polar(mpt,ref,&type,&unit,&x,&y,&z,&value,&error) == 1) 
    {
        child = add_xml_child(exp,exml_names[exmlPOLARIZABILITY]);
        add_xml_char(child,exml_names[exmlTYPE],type);
        add_xml_char(child,exml_names[exmlUNIT],unit);
        ptr = gmx_ftoa(value);
        add_xml_child_val(child,exml_names[exmlVALUE],ptr);
        sfree(ptr);
        ptr = gmx_ftoa(error);
        add_xml_child_val(child,exml_names[exmlERROR],ptr);
        sfree(ptr);
        if ((x != 0) || (y != 0) || (z != 0)) 
        {
            ptr = gmx_ftoa(x);
            add_xml_child_val(child,exml_names[exmlX],ptr);
            sfree(ptr);
            ptr = gmx_ftoa(y);
            add_xml_child_val(child,exml_names[exmlY],ptr);
            sfree(ptr);
            ptr = gmx_ftoa(z);
            add_xml_child_val(child,exml_names[exmlZ],ptr);
            sfree(ptr);
        }
        sfree(type);
        sfree(unit);
    }
    while (gmx_molprop_get_quadrupole(mpt,ref,&type,&unit,&xx,&yy,
                                      &zz,&xy,&xz,&yz) == 1) 
    {
        child = add_xml_child(exp,exml_names[exmlQUADRUPOLE]);
        add_xml_char(child,exml_names[exmlTYPE],type);
        add_xml_char(child,exml_names[exmlUNIT],unit);
        ptr = gmx_ftoa(xx);
        add_xml_child_val(child,exml_names[exmlXX],ptr);
        sfree(ptr);
        ptr = gmx_ftoa(yy);
        add_xml_child_val(child,exml_names[exmlYY],ptr);
        sfree(ptr);
        ptr = gmx_ftoa(zz);
        add_xml_child_val(child,exml_names[exmlZZ],ptr);
        sfree(ptr);
        ptr = gmx_ftoa(xy);
        add_xml_child_val(child,exml_names[exmlXY],ptr);
        sfree(ptr);
        ptr = gmx_ftoa(xz);
        add_xml_child_val(child,exml_names[exmlXZ],ptr);
        sfree(ptr);
        ptr = gmx_ftoa(yz);
        add_xml_child_val(child,exml_names[exmlYZ],ptr);
        sfree(ptr);
        sfree(type);
        sfree(unit);
    }
}

static void add_xml_molprop(xmlNodePtr parent,gmx_molprop_t mpt)
{
    xmlNodePtr ptr,child,grandchild,comp,atomptr,baby;
    int    i,ai,aj,bondorder,eMP,atomid,cnumber,expref,calcref,atomref;
    char   *program,*basisset,*method,*name,*type,*unit,*reference;
    char   *catom,*atomname,*obtype,*coords,*conformation,*p;
    const  char   *iupac,*cas,*cid,*inchi,*category,*composition;
    double q,x,y,z,xx,yy,zz,xy,xz,yz,aver,error;
  
    gmx_molprop_reset_experiment(mpt);
    gmx_molprop_reset_calculation(mpt);
    ptr = add_xml_child(parent,exml_names[exmlMOLECULE]);
    add_xml_char(ptr,exml_names[exmlMOLNAME],gmx_molprop_get_molname(mpt));
    add_xml_char(ptr,exml_names[exmlFORMULA],gmx_molprop_get_formula(mpt));
    add_xml_double(ptr,exml_names[exmlMASS],gmx_molprop_get_mass(mpt));
    add_xml_double(ptr,exml_names[exmlCHARGE],gmx_molprop_get_charge(mpt));
    add_xml_double(ptr,exml_names[exmlMULTIPLICITY],gmx_molprop_get_multiplicity(mpt));
  
    iupac = gmx_molprop_get_iupac(mpt);  
    cas   = gmx_molprop_get_cas(mpt);  
    cid   = gmx_molprop_get_cid(mpt);  
    inchi = gmx_molprop_get_inchi(mpt);  
    child = add_xml_child(ptr,exml_names[exmlMOLINFO]);
    add_xml_char(child,exml_names[exmlIUPAC],iupac);
    add_xml_char(child,exml_names[exmlCAS],cas);
    add_xml_char(child,exml_names[exmlCID],cid);
    add_xml_char(child,exml_names[exmlINCHI],inchi);
    while (gmx_molprop_get_bond(mpt,&ai,&aj,&bondorder) == 1)
    {
        child = add_xml_child(ptr,exml_names[exmlBOND]);
        add_xml_int(child,exml_names[exmlAI],ai);
        add_xml_int(child,exml_names[exmlAJ],aj);
        add_xml_int(child,exml_names[exmlBONDORDER],bondorder);
    }
    while (gmx_molprop_get_experiment(mpt,&reference,&conformation,&expref) == 1) 
    {
        child = add_xml_child(ptr,exml_names[exmlEXPERIMENT]);
        add_xml_char(child,exml_names[exmlREFERENCE],reference);
        add_xml_char(child,exml_names[exmlCONFORMATION],conformation);
        
        add_properties(child,mpt,expref);
        
        sfree(reference);
    }
    gmx_molprop_reset_experiment(mpt);
    while (gmx_molprop_get_calculation(mpt,&program,&method,&basisset,
                                       &reference,&conformation,&calcref) == 1) 
    {
        child = add_xml_child(ptr,exml_names[exmlCALCULATION]);
        add_xml_char(child,exml_names[exmlPROGRAM],program);
        add_xml_char(child,exml_names[exmlMETHOD],method);
        add_xml_char(child,exml_names[exmlBASISSET],basisset);
        add_xml_char(child,exml_names[exmlREFERENCE],reference);
        add_xml_char(child,exml_names[exmlCONFORMATION],conformation);

        add_properties(child,mpt,calcref);
        
        while (gmx_molprop_calc_get_atom(mpt,calcref,&atomname,&obtype,&atomid,&atomref) == 1) 
        {
            grandchild = add_xml_child(child,exml_names[exmlATOM]);
            add_xml_char(grandchild,exml_names[exmlNAME],atomname);
            add_xml_char(grandchild,exml_names[exmlOBTYPE],obtype);
            add_xml_int(grandchild,exml_names[exmlATOMID],atomid);
            
            if (gmx_molprop_calc_get_atomcoords(mpt,calcref,atomref,&unit,&x,&y,&z) == 1) 
            {
                p = gmx_ftoa(x);
                baby = add_xml_child_val(grandchild,exml_names[exmlX],p);
                sfree(p);
                add_xml_char(baby,exml_names[exmlUNIT],unit);
                p = gmx_ftoa(y);
                baby = add_xml_child_val(grandchild,exml_names[exmlY],p);
                sfree(p);
                add_xml_char(baby,exml_names[exmlUNIT],unit);
                p = gmx_ftoa(z);
                baby = add_xml_child_val(grandchild,exml_names[exmlZ],p);
                sfree(p);
                add_xml_char(baby,exml_names[exmlUNIT],unit);
                sfree(unit);
            }
            while (gmx_molprop_calc_get_atomcharge(mpt,calcref,atomref,
                                                   &type,&unit,&q) == 1) 
            {
                p = gmx_ftoa(q);
                atomptr = add_xml_child_val(grandchild,exml_names[exmlQ],p);
                sfree(p);
                add_xml_char(atomptr,exml_names[exmlTYPE],type);
                add_xml_char(atomptr,exml_names[exmlUNIT],unit);
                sfree(type);
                sfree(unit);
            }
            sfree(atomname);
        }
        sfree(program);
        sfree(method);
        sfree(basisset);
        sfree(reference);
    }
    gmx_molprop_reset_calculation(mpt);
    while ((category = gmx_molprop_get_category(mpt)) != NULL) 
    {
        child = add_xml_child(ptr,exml_names[exmlCATEGORY]);
        add_xml_char(child,exml_names[exmlCATNAME],category);
    }
	
    while ((composition = gmx_molprop_get_composition(mpt)) != NULL) 
    {
        child = add_xml_child(ptr,exml_names[exmlCOMPOSITION]);
        add_xml_char(child,exml_names[exmlCOMPNAME],composition);
        while ((gmx_molprop_get_composition_atom(mpt,composition,&catom,&cnumber)) == 1) 
        {
            grandchild = add_xml_child(child,exml_names[exmlCATOM]);
            add_xml_char(grandchild,exml_names[exmlC_NAME],catom);
            add_xml_int(grandchild,exml_names[exmlC_NUMBER],cnumber);
        }
    }
  
}

void gmx_molprops_write(const char *fn,int nmolprop,gmx_molprop_t mpt[],
                        gmx_bool bCompress)
{
    xmlDocPtr  doc;
    xmlDtdPtr  dtd;
    xmlNodePtr myroot;
    int        i,nmt;
    xmlChar    *libdtdname,*dtdname,*gmx;
  
    gmx        = (xmlChar *) "molecules";
    dtdname    = (xmlChar *) "molprops.dtd";
    libdtdname = dtdname;
  
    if ((doc = xmlNewDoc((xmlChar *)"1.0")) == NULL)
        gmx_fatal(FARGS,"Creating XML document","");
    
    if ((dtd = xmlCreateIntSubset(doc,dtdname,libdtdname,dtdname)) == NULL)
        gmx_fatal(FARGS,"Creating XML DTD","");
    
    if ((myroot = xmlNewDocNode(doc,NULL,gmx,NULL)) == NULL)
        gmx_fatal(FARGS,"Creating root element","");
    dtd->next = myroot;
    myroot->prev = (xmlNodePtr) dtd;
    
    /* Add molecule definitions */
    for(i=0; (i<nmolprop); i++)
    {
        printf("Adding %d/%d %s\n",i+1,nmolprop,gmx_molprop_get_molname(mpt[i]));
        add_xml_molprop(myroot,mpt[i]);
    }
    xmlSetDocCompressMode(doc,(int)bCompress);
    xmlIndentTreeOutput = 1;
    if (xmlSaveFormatFileEnc(fn,doc,"ISO-8859-1",1) == 0)
        gmx_fatal(FARGS,"Saving file",fn);
    xmlFreeDoc(doc);
}
