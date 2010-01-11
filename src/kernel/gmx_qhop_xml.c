#include <string.h>
#include "gmx_fatal.h"
#include "smalloc.h"
#include "macros.h"
#include "futil.h"
#include "types/gmx_qhop_types.h"
#include "gmx_qhop_parm.h"
#include "gmx_qhop_db.h"

#ifdef HAVE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>
	
extern int xmlDoValidityCheckingDefaultValue;
	
#define NN(x) (NULL != (x))
#define CPYLEN 128 /* This is the default length used for strncpy.
		    * I avoid using strcpy becase of potential buffer-overflows
		    * Note that the careless use of CPYLEN isn't safe either, but
		    * as the code maturates, CPYLEN may be replaced by explicit numbers.
		    * - Erik Marklund */
#define RNMLEN 6 /* Ditto for resnames. This can probabaly be a bit shorter,
		  * but let's keep it for now. */
#define ANMLEN 6 /* Ditto for atomnames. */

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
  exmlQHOPS,
  exmlRESBLOCKS,
  exmlRESIDUE_TYPE, exmlFILE, exmlRES, exmlPROTON, exmlPRODUCT,
  exmlQHOP, exmlDONOR, exmlACCEPTOR,
  exmlPARAM, exmlNAME, exmlVALUE, 
  exmlUNIT, exmlDON_ATOM, exmlACC_ATOM,
  exmlNR 
};
  
static const char *exml_names[exmlNR] = {
  "qhops",
  "resblocks",
  "residue_type", "file", "res", "proton", "product",
  "qhop", "donor", "acceptor", "parameter", 
  "name", "value", "unit", "don_atom", "acc_atom"
};

static int find_elem(char *name,int nr,const char *names[])
{
  int i;
  
  for(i=0; (i<nr); i++)
    if (strcmp(name,names[i]) == 0) 
      break;
  if (i == nr)
    gmx_fatal(FARGS,"Unknown element name %s",name);
    
  return i;
}

void add_xml_int(xmlNodePtr ptr,const char *name,int val)
{
  xmlChar buf[32];
  
  sprintf((char *)buf,"%d",val);
  if (xmlSetProp(ptr,(xmlChar *)name,buf) == 0)
    gmx_fatal(FARGS,"Setting",(char *)name);
}

void add_xml_double(xmlNodePtr ptr,const char *name,double val)
{
  xmlChar buf[32];
  
  sprintf((char *)buf,"%g",val);
  if (xmlSetProp(ptr,(xmlChar *)name,buf) == 0)
    gmx_fatal(FARGS,"Setting",(char *)name);
}

void add_xml_char(xmlNodePtr ptr,const char *name,char *val)
{
  if (xmlSetProp(ptr,(xmlChar *)name,(xmlChar *)val) == 0)
    gmx_fatal(FARGS,"Setting",(char *)name);
}

xmlNodePtr add_xml_child(xmlNodePtr parent,const char *type)
{
  xmlNodePtr child;
  
  if ((child = xmlNewChild(parent,NULL,(xmlChar *)type,NULL)) == NULL)
    gmx_fatal(FARGS,"Creating element",(char *)type);
  
  return child;
}

xmlNodePtr add_xml_comment(xmlDocPtr doc,
			   xmlNodePtr prev,char *comment)
{
  xmlNodePtr comm,ptr;
  
  if ((comm = xmlNewComment((xmlChar *)comment)) == NULL)
    gmx_fatal(FARGS,"Creating doc comment element","");
  ptr = prev;
  while (ptr->next != NULL)
    ptr=ptr->next;
  ptr->next    = comm;
  comm->prev   = ptr;
  comm->doc    = doc;
  
  return comm;
}

static char *sp(int n, char buf[], int maxindent)
{
  int i;
  if(n>=maxindent)
    n=maxindent-1;
  
  /* Don't indent more than maxindent characters */
  for(i=0; (i<n); i++)
    buf[i] = ' ';
  buf[i] = '\0';
  
  return buf;
}

/* Finds the first matching string in an array of strings. Returns -1 if not found. */
static int qhop_find_name(char **ss, const char *name, int n)
{
  int i;
  for (i=0; i<n; i++)
    if (strcmp(ss[i], name) == 0)
      return i;

  return -1;/* not found */
}


/* add_to_record should be replaced with strdup!*/

/* Copies the text in src to a new element at the end of rec and,
 * if update_index == eUPDATE_INDEX, it also updates n.
 * len is the maximal length of the data being copied. */
/* static void add_to_record(const char *src, char **rec, int *n, int len, index_switch update_index) */
/* { */
/*   snew(*rec, len+1); /\* Allocate name buffer *\/ */
/*   strncpy(*rec, src, len); */
/*   (*rec)[len] = '\0';  /\* just in case *\/ */
/*   if (update_index == eUPDATE_INDEX) */
/*     (*n)++; */
/* } */

/* Adds a file with name *f to a qhop_interactions */
static void qi_add_file(qhop_interactions_t qi, char *f)
{
  if NN(f)
    {
      if (qhop_find_name(qi->files, f, qi->nf) < 0) /* New file? */
	{
	  printf("Adding file %s\n", f);
	  /* Copy the filename to the array in qi.*/
	  srenew(qi->files, qi->nf+1);
	  qi->files[qi->nf++] = strdup(f);
	  /*add_to_record(f, qi->files, &(qi->nf), CPYLEN, eUPDATE_INDEX);*/
	}
    }
  else
    gmx_fatal(FARGS, "No filename found.");
}

static void rb_add_restype(qhop_resblocks_t rb, currentRes *ri, char *name)
{
  
  if NN(name)
    {
      printf("Adding restype %s\n", name);
      /* Is this a new resblock? If not, make a new one. */
      if ((ri->rt = qhop_find_name(rb->restype, name, ri->r)) < 0)
	{
	  srenew(rb->restype, rb->nrestypes+1);       /* Make new resblock */
	  srenew(rb->nres, rb->nrestypes+1);          /* Make new resindex */
	  rb->nres[rb->nrestypes] = 0;
	  srenew(rb->res, rb->nrestypes+1);  /* Make new resarray for this resblock */	    
	  ri->rt = rb->nrestypes;            /* Make this resblock the current one */
	  ri->r  = -1;                       /* just to be sure... */
	  ri->da = -1;
	  rb->restype[rb->nrestypes++] = strdup(name);
	  /*add_to_record(name, &(rb->restype[rb->nrestypes]), &(rb->nrestypes), RNMLEN, eUPDATE_INDEX);*/
	}
    }
  else
    gmx_fatal(FARGS, "No name found for this resblock.");
}

static void rb_add_res(qhop_resblocks_t rb, currentRes *ri, char *name)
{
  bool found;
  int i;
  if NN(name)
    {
      printf("Adding res %s\n", name);
      /* Have we seen this res before for this resblock? */
      found = FALSE;
      for (i=0; i < rb->nres[ri->rt] && !found; i++)
	if (strcmp(rb->res[ri->rt][i].name, name) == 0)
	  found = TRUE;

      if (!found)
	{
	  ri->r = rb->nres[ri->rt];
	  if (ri->r != 0)
	    srenew(rb->res[ri->rt], (rb->nres[ri->rt])+1);
	  else
	    snew(rb->res[ri->rt], 1);
	  rb->res[ri->rt][ri->r].na = 0;
	  rb->res[ri->rt][ri->r].nd = 0;
	  rb->res[ri->rt][rb->nres[ri->rt]++].name = strdup(name);
	  /*add_to_record(name, &(rb->res[ri->rt][ri->r].name),
	    &(rb->nres[ri->rt]), RNMLEN, eUPDATE_INDEX);*/
	}
      else
	{
	  ri->r = i;
	}
    }
  else
    gmx_fatal(FARGS, "No name found for this residue.");

}

static void rb_add_proton(qhop_res **res, currentRes *ri, char *name)
{
  qhop_reactant *r;
  int nr;
 
  r = (ri->Acc == exmlACCEPTOR) ?
    &(res[ri->rt][ri->r].acc[ri->da]) :
    &(res[ri->rt][ri->r].don[ri->da]);
      
  if NN(name)
    {
      printf("Adding proton %s\n", name);

      if (qhop_find_name(r->H, name, r->nH) < 0) /* New proton? */
	{
	  if (r->nH != 0)
	    srenew(r->H, r->nH + 1);
	  else
	    snew(r->H, 1);

	  r->H[r->nH++] = strdup(name);
	  /*add_to_record(name,
	    rb->res[ri->rt][ri->r].acc[ri->da].H,
	    &(rb->res[ri->rt][ri->r].acc[ri->da].nH),
	    ANMLEN, eUPDATE_INDEX);*/
	}
    }
}

/* adds a qhop_reactant with name and product to a qhop_res 'res'.
 * If acc_or_don == exmlACCEPTOR then this is an acceptor, otherwise a donor. */
static void rb_add_name_product(qhop_res_t res, const char *name, const char *product, currentRes *ri)
{
  int i;
  bool found;
  
  if (NN(name) && NN(product))
    {
      printf("Adding name %s and product %s \n", name, product);
            
      found = FALSE;
      for(i=0; i < (ri->Acc == exmlACCEPTOR ? res->na : res->nd) && !found; i++)
	if (strcmp(ri->Acc == exmlACCEPTOR ?
		   res->acc[i].name, name : res->don[i].name, name) == 0)
	  {
	    found = TRUE;
	    ri->da = i;
	  }
      if (!found)
	{
	  if (ri->Acc == exmlACCEPTOR)
	    {
	      if (res->na != 0)
		srenew(res->acc, (res->na)+1);
	      else
		snew(res->acc, 1);
	      
	      res->acc[res->na].name    = strdup(name);
	      res->acc[res->na].product = strdup(product);
	      res->acc[res->na].nH = 0;
	      ri->da = (res->na)++;
	    }
	  else /* A donor */
	    {
	      if (res->nd != 0)
		srenew(res->don, res->nd+1);
	      else
		snew(res->don, 1);
	      
	      res->don[res->nd].name    = strdup(name);
	      res->don[res->nd].product = strdup(product);
	      res->don[res->nd].nH = 0;
	      ri->da = (res->nd)++;
	    }
	}
      else
	/* Otherwise there's nothing to add. */
	printf("This %s has been found before: %s.\n",
	       ri->Acc==exmlACCEPTOR ? "acceptor" : "donor", name);
    }
  else
    gmx_fatal(FARGS, "%s without %s%s%s. Can't have that!",
	      ri->Acc==exmlACCEPTOR   ? "Acceptor" : "Donor",
	      NN(name)                ? "": "name",
	      NN(name) || NN(product) ? "":" and ",
	      NN(product)             ? "": "product");
}

/* rt: residue type/resblock
 * r:  res in resblock. */
static void qhop_process_attr(FILE *fp,xmlAttrPtr attr,int parent,
			      int elem,int indent, xmlrec_t xml, currentRes *ri)
{
  char *attrname,*attrval;
  char buf[100];
  int  i,kkk,eprop;
  bool found;
  char *xbuf[exmlNR];
  
  for(i=0; (i<exmlNR); i++)
    xbuf[i] = NULL;
  while (attr != NULL) {
    attrname = (char *)attr->name;
    attrval  = (char *)attr->children->content;
    
#define atest(s) ((strcasecmp(attrname,s) == 0) && (attrval != NULL))
    kkk = find_elem(attrname,exmlNR,exml_names);
    if (attrval != NULL)
      xbuf[kkk] = strdup(attrval);
      
    if (fp)
      fprintf(fp,"%sProperty: '%s' Value: '%s'\n",sp(indent,buf,99),
	      attrname,attrval);
    attr = attr->next;
#undef atest
  }
  
  switch (elem) {
  case exmlQHOP:
    if (NN(xbuf[exmlDONOR]) && NN(xbuf[exmlACCEPTOR])) {
      /*      gmx_qhop_set_donor(xml->qht,xbuf[exmlDONOR]);
	      gmx_qhop_set_acceptor(xml->qht,xbuf[exmlACCEPTOR]);*/
      printf("Adding qhop: ");
      qhop_set_donor(xml->gqh[xml->nqh-1],xbuf[exmlDONOR]);
      qhop_set_acceptor(xml->gqh[xml->nqh-1],xbuf[exmlACCEPTOR]);
    }
    break;
  case exmlPARAM:
    if (NN(xbuf[exmlNAME]) && NN(xbuf[exmlUNIT]) && NN(xbuf[exmlVALUE])) {
      /*      gmx_qhop_add_param(xml->qht,xbuf[exmlNAME],xbuf[exmlVALUE],
			 xbuf[exmlUNIT]);*/
      qhop_add_param(xml->gqh[xml->nqh-1],xbuf[exmlNAME],xbuf[exmlVALUE],
			 xbuf[exmlUNIT]);
    }
    break;
    /* Here's where some resblocks stuff needs to be implemented. */
  case exmlRESIDUE_TYPE: /* 'name' */
    rb_add_restype(xml->rb, ri, xbuf[exmlNAME]);
    break;
  case exmlRES: /* 'name' */
    rb_add_res(xml->rb, ri, xbuf[exmlNAME]);
    break;

  case exmlACCEPTOR: /* 'name', 'product' */
    rb_add_name_product(&(xml->rb->res[ri->rt][ri->r]), xbuf[exmlNAME], xbuf[exmlPRODUCT], ri);
    break;
  case exmlDONOR: /* 'name', 'product' */
    rb_add_name_product(&(xml->rb->res[ri->rt][ri->r]), xbuf[exmlNAME], xbuf[exmlPRODUCT], ri);
    break;

  default:
    break;
  }
  

  for(i=0; (i<exmlNR); i++)
    if (NN(xbuf[i]))
      sfree(xbuf[i]);
}

static void qhop_process_element(FILE *fp, xmlNodePtr tree, int parent,
				 int indent, xmlrec_t xml, currentRes *ri)
{
  int elem, i;
  char buf[100];
  
  elem = find_elem((char *)tree->name,exmlNR,exml_names);
  if (fp)
    fprintf(fp,"%sElement node name %s\n",sp(indent,buf,99),
	    (char *)tree->name);
  switch(elem)
    {
    case exmlQHOP:
      xml->nqh++;
      srenew(xml->gqh,xml->nqh);
      xml->gqh[xml->nqh-1] = qhop_init();
      break;

    case exmlRESBLOCKS:
      if (!xml->rb)
	snew(xml->rb,1);
      break;

    case exmlRESIDUE_TYPE:
      /* A new resblock will be created if this is indeed a new one.
       * That all happens in qhop_process_attr below. *rt will also be set accordingly.
       * Thus, nothing need to be done here.*/
      break;

    case exmlFILE:
      qi_add_file(xml->qi, (char*) (tree->children->content));
      break;

    case exmlRES:
      /* The identity of the resblock was determined at
       * the last encounter of an exmlRESIDUE_TYPE and the
       * subsequent call to qhop_process_attr.
       * r will be set in the next call to qhop_process_attr.
       * Thus, nothing to be done here */
      break;

      /* ri->da will be set when any of the following two are encountered. */
    case exmlACCEPTOR:
      /* A new acceptor will be added in the next call
       * to qhop_process_attr. */
      ri->Acc = exmlACCEPTOR;
      break;

    case exmlDONOR:
      /* A new donor will be added in the next call
       * to qhop_process_attr. */
      ri->Acc = exmlDONOR;
      break;

    case exmlPROTON:
      /* Now add this proton to the current acceptor/donor. */
      rb_add_proton(xml->rb->res, ri, (char*)(tree->children->content));
      break;

    case exmlDON_ATOM:
      /* This is a child of qhop. */
      qhop_set_don_atom(xml->gqh[xml->nqh-1], (char*)(tree->children->content));
      break;

    case exmlACC_ATOM:
      /* This is a child of qhop. */
      qhop_set_acc_atom(xml->gqh[xml->nqh-1], (char*)(tree->children->content));
      break;

    default:
      break;
    }
  /*  if (elem != exmlQHOPS)*/
  if NN(tree->properties)
    qhop_process_attr(fp,tree->properties,parent,
		      elem,indent+2,/*xml->gqh[xml->nqh-1],*/
		      xml, ri);
}
 
static void qhop_process_tree(FILE *fp, xmlNodePtr tree, int parent,
			      int indent, xmlrec_t xml, currentRes *ri)
{
  char buf[100];
  int  elem;
  
  while (tree != NULL) {
    if (fp) {
      if ((tree->type > 0) && (tree->type < NXMLTYPES))
	fprintf(fp,"Node type %s encountered with name %s\n",
		xmltypes[tree->type],(char *)tree->name);
      else
	fprintf(fp,"Node type %d encountered\n",tree->type);
    }
    
    switch (tree->type) {
    case XML_ELEMENT_NODE:
      qhop_process_element(fp, tree, parent, indent+2, xml, ri)/*, rb, &rt, &r, &da, qi)*/;
      
      if (tree->children) {
	elem = find_elem((char *)tree->name,exmlNR,exml_names);
	qhop_process_tree(fp,tree->children,elem,indent+2,xml,ri);
      }
      break;
    default:
      break;
    }
    tree = tree->next;
  }
}

void qhops_read(char *fn, qhop_db_t qdb)
{
  xmlDocPtr     doc;
  int           i,npd;
  xmlrec        *xml;
  const char *db="qhops.dat";
  bool fna=FALSE;
  currentRes    ri;
  
  ri.Acc = FALSE;
  ri.rt  = -1;
  ri.r   = -1;
  ri.da  = -1;

  xmlDoValidityCheckingDefaultValue = 0;
  if (NULL == fn) 
    {
      fn = (char *)libfn(db);
      fna=TRUE;
    }
  if ((doc = xmlParseFile(fn)) == NULL)
    {
      fprintf(stderr,"Reading XML file %s. Run a syntax checker such as nsgmls.", fn);
      exit(1);
    }

  snew(xml,1);
  snew(xml->qi, 1);
  xml->qi->nf = 0;

  qhop_process_tree(NULL,doc->children,0,0,xml,&ri);
  
  xmlFreeDoc(doc);
  if (fna)
    sfree(fn);

  /*snew(qdb,1);*/

  

  qdb->ngqh = xml->nqh;
  qdb->gqh = xml->gqh;
  qdb->qi = xml->qi;
  qdb->rb = xml->rb;
  
  /* *nqhop = xml->nqh; */

}

static void add_xml_qhop(xmlNodePtr parent,qhop_t qht)
{
  xmlNodePtr ptr,child,grandchild,comp;
  char   *name,*type,*value,*unit;
  
  ptr = add_xml_child(parent,exml_names[exmlQHOP]);
  add_xml_char(ptr,exml_names[exmlDONOR],qhop_get_donor(qht));
  add_xml_char(ptr,exml_names[exmlACCEPTOR],qhop_get_acceptor(qht));
  
  while (qhop_get_param(qht,&name,&value,&unit) == 1) {
    child = add_xml_child(ptr,exml_names[exmlPARAM]);
    add_xml_char(child,exml_names[exmlNAME],name);
    add_xml_char(child,exml_names[exmlVALUE],value);
    add_xml_char(child,exml_names[exmlUNIT],unit);
    sfree(name);
    sfree(value);
    sfree(unit);
  }
}

void qhops_write(char *fn,int nqhop,qhop_t qht[])
{
  xmlDocPtr  doc;
  xmlDtdPtr  dtd;
  xmlNodePtr myroot;
  int        i,nmt;
  xmlChar    *libdtdname,*dtdname,*gmx;
  
  gmx        = (xmlChar *) "qhops";
  dtdname    = (xmlChar *) "qhops.dtd";
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
  for(i=0; (i<nqhop); i++)
    add_xml_qhop(myroot,qht[i]);

  xmlSetDocCompressMode(doc,0);
  xmlIndentTreeOutput = 1;
  if (xmlSaveFormatFileEnc(fn,doc,"ISO-8859-1",2) == 0)
    gmx_fatal(FARGS,"Saving file",fn);
  xmlFreeDoc(doc);
}



#else

qhop_t qhops_read(char *fn,int *nqhop)
{
  gmx_fatal(FARGS,"You need to configure the software with --with-xml for function qhops_read to work");
  return NULL;
}

void qhops_write(char *fn,int nqhop,qhop_t qht)
{
  gmx_fatal(FARGS,"You need to configure the software with --with-xml for function qhops_write to work");
}

#endif
