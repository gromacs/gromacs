#include <stdio.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "typedefs.h"
#include "fatal.h"
#include "string.h"
#include "smalloc.h"
#include "names.h"
#include "assert.h"
#include "symtab.h"
#include "macros.h"
#include "symtab.h"
#include "xmlio.h"

typedef struct {
  int        nx,nv,nf,ntop,nbox,ninprec;
  int        step,natoms;
  real       t,lambda;
  t_inputrec *ir;
  rvec       *box;
  rvec       *x,*v,*f;
  t_topology *top;
} t_xmlrec;

typedef struct {
  char *name,*value;
} t_masstype;

static char *xmltypes[] = { 
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
	
extern int xmlDoValidityCheckingDefaultValue;
	
enum { 
  exmlGROMACS, 
  exmlINPUTREC,    exmlOUTPUT,      exmlCOUPLING,    exmlCUTOFF,
  exmlPMEPARM,     exmlTCOUPLING,   exmlPCOUPLING,   exmlTCPARM,      
  exmlREFPRES,     exmlCOMPRESS,    exmlRVEC,
  exmlSYSTEM,      exmlCOMPOSITION, exmlMOLECULE,    exmlATOM,
  exmlTOPOLOGY,    exmlBONDS,       exmlANGLES,      exmlDIHEDRALS,
  exmlFORCEFIELD,  exmlMASSTYPE, 
  exmlCOORDINATES, exmlVELOCITIES, exmlFORCES,     exmlNR };
  
static char *exml_names[] = {
  "gromacs",
  /* Inputrec stuff */
  "inputrec",    "output",     "coupling", "cutoff",   "pmeparm",
  "tcoupling",   "pcoupling",  "tcparm",   "refpres",  "compress", "rvec",
  /* System description */
  "system",      "composition", "molecule","atom",     
  /* Topology description */
  "topology",
  "bonds",       "angles",      "dihedrals",
  /* Force field */
  "forcefield",  "masstype",    
  /* Coordinates etc. */
  "coordinates", "velocities", "forces" 
};
	
static int find_elem(char *name,int nr,char *names[])
{
  int i;
  
  for(i=0; (i<nr); i++)
    if (strcmp(name,names[i]) == 0) 
      break;
  if (i == nr)
    fatal_error(0,"Unknown element name %s",name);
    
  return i;
}
	
static char *sp(int n)
{
  static char buf[80];
  int i;
  
  /* Don't indent more than 80 characters */
  n = n % 80;
  for(i=0; (i<n); i++)
    buf[i] = ' ';
  buf[i] = '\0';
  
  return buf;
}

static void process_attr(FILE *fp,xmlAttrPtr attr,int elem,
			 int indent,t_xmlrec *xml)
{
  char *attrname,*attrval;
  
  while (attr != NULL) {
    attrname = (char *)attr->name;
    attrval  = (char *)attr->children->content;
    
#define atest(s) ((strcmp(attrname,s) == 0) && (attrval != NULL))
    switch (elem) {
    case exmlGROMACS:
      if (atest("title"))
	xml->top->name = put_symtab(&xml->top->symtab,attrval);
      break;
    case exmlINPUTREC:
      if (atest("algorithm"))
	xml->ir->eI = find_elem(attrval,eiNR,ei_names);
      break;
    case exmlOUTPUT:      
    case exmlCOUPLING:    
    case exmlCUTOFF:
    case exmlPMEPARM:     
    case exmlTCOUPLING:  
    case exmlPCOUPLING:  
    case exmlTCPARM:      
    case exmlREFPRES:    
    case exmlCOMPRESS:
    case exmlRVEC:
    case exmlSYSTEM:   
    case exmlCOMPOSITION:
    case exmlMOLECULE:
    case exmlATOM:
    case exmlTOPOLOGY:  
    case exmlBONDS:  
    case exmlANGLES:  
    case exmlDIHEDRALS:
    case exmlFORCEFIELD: 
    case exmlMASSTYPE: 
    case exmlCOORDINATES:
    case exmlVELOCITIES: 
    case exmlFORCES: 
    default:
      if (fp)
	fprintf(fp,"%sProperty: '%s' Value: '%s'\n",sp(indent),
		attrname,attrval);
    }
    attr = attr->next;
#undef atest
  }
}

static void process_tree(FILE *fp,xmlNodePtr tree,int indent,t_xmlrec *xml)
{
  int elem;
  
  while (tree != NULL) {
    switch (tree->type) {
    case XML_ELEMENT_NODE:
      elem = find_elem((char *)tree->name,exmlNR,exml_names);
      if (fp)
	fprintf(fp,"%sElement node name %s\n",sp(indent),(char *)tree->name);
      
      process_attr(fp,tree->properties,elem,indent+2,xml);
      
      if (tree->children)
	process_tree(fp,tree->children,indent+2,xml);
      break;
    case XML_COMMENT_NODE:
      if (fp)
	fprintf(fp,"Comment node encountered\n");
      break;
    case XML_PI_NODE:
    case XML_TEXT_NODE:
      /* Silently ignore these for now */
      break;
    default:
      if (fp) {
	if ((tree->type > 0) && (tree->type < NXMLTYPES))
	  fprintf(fp,"Node type %s encountered with name %s\n",
		  xmltypes[tree->type],(char *)tree->name);
	else
	  fprintf(fp,"Node type %d encountered\n",tree->type);
      }
    }
    tree = tree->next;
  }
}
	
void read_xml(char *fn,int *step,real *t,real *lambda,
	      t_inputrec *ir,rvec *box,int *natoms,
	      rvec **x,rvec **v,rvec **f,t_topology *top)
{
  xmlDocPtr  doc;
  t_xmlrec   *xml;
  
  xmlDoValidityCheckingDefaultValue = 1;
  assert(asize(exml_names) == exmlNR);
  if ((doc = xmlParseFile(fn)) == NULL)
    fatal_error(0,"Reading XML file %s. Run a syntax checker such as nsgmls.",
		fn);

  snew(xml,1);
  xml->ir     = ir;
  xml->box    = box;
  xml->top    = top;
  process_tree(debug,doc->children,0,xml);
  
  xmlFreeDoc(doc);
  sfree(xml);
}

static void add_xml_inputrec(xmlDocPtr doc,t_inputrec *ir)
{
  
}

static void add_xml_molecule(xmlDocPtr doc,t_atoms *atoms,
			     int nmt,t_masstype mt[])
{
  
}

static void add_xml_idef(xmlDocPtr doc,t_idef *idef)
{
  
}

static void add_xml_rvecs(xmlDocPtr doc,char *name,int natoms,rvec *xvf)
{
  
}

static t_masstype *mk_masstype(int nmol,t_atoms atoms[],int *nmt)
{
  int        i,j,k,nm;
  t_masstype *mt=NULL;
  char       buf[12];
  nm = 0;
  for(i=0; (i<nmol); i++) {
    for(j=0; (j<atoms[i].nr); j++) {
      for(k=0; (k<nm); k++)
	if (strcmp(*atoms[i].atomname[j],mt[k].name) == 0)
	  break;
      if (k == nm) {
	srenew(mt,nm+1);
	mt[nm].name = strdup(*atoms[i].atomname[j]);
	sprintf(buf,"%.5f",atoms[i].atom[j].m);
	mt[nm].value = strdup(buf);
	nm++;
      }
    }
  }
  *nmt = nm;
  return mt;
}

static void add_xml_masstype(xmlDocPtr doc,int nmt,t_masstype mt[])
{
  
}

void write_xml(char *fn,char *title,t_inputrec *ir,rvec *box,
	       int natoms,rvec *x,rvec *v,rvec *f,
	       int nmol,t_atoms atoms[],t_idef *idef)
{
  xmlDocPtr  doc;
  t_masstype *mt;
  int        i,nmt;
  
  if ((doc = xmlNewDoc("1.0")) == NULL)
    fatal_error(0,"Creating document");
  if ((doc->children = xmlNewDocNode(doc,NULL,"gromacs",NULL)) == NULL)
    fatal_error(0,"Creating root element");
    
  /* Add DTD declaration and stylesheet information here */
  
  /* Title of the system */
  if (title)
    if (xmlSetProp(doc->children,"title",title) == 0)
      fatal_error(0,"Setting title %s",title);
  
  /* Add inputrec */
  if (ir) 
    add_xml_inputrec(doc,ir);
  
  /* Generate masstypes */
  mt = mk_masstype(nmol,atoms,&nmt);
  
  /* Add molecule definitions */
  for(i=0; (i<nmol); i++)
    add_xml_molecule(doc,&(atoms[i]),nmt,mt);

  /* Add force field */  
  if (idef)
    add_xml_idef(doc,idef);
    
  /* Coordinates */
  if (x)
    add_xml_rvecs(doc,"coordinates",natoms,x);
  
  /* Velocities */
  if (v)
    add_xml_rvecs(doc,"velocities",natoms,v);
  
  /* Forces */  
  if (f)
    add_xml_rvecs(doc,"forces",natoms,f);
  
  if (xmlSaveFile(fn,doc) == 0)
    fatal_error(0,"Saving file %s",fn);
  xmlFreeDoc(doc);
}
