/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_matio_c = "$Id$";

#include <ctype.h>
#include "sysstuff.h"
#include "futil.h"
#include "string2.h"
#include "macros.h"
#include "smalloc.h"
#include "fatal.h"
#include "matio.h"
#include "statutil.h"

#define round(a) (int)(a+0.5)

static char mapper[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!@#$%^&*()-_=+{}|;:',<.>/?";
#define NMAP strlen(mapper)

bool matelmt_cmp(t_xpmelmt e1, t_xpmelmt e2) 
{ 
  return (e1.c1 == e2.c1) && (e1.c2 == e2.c2);
}
    
t_matelmt searchcmap(int n,t_mapping map[],t_xpmelmt c)
{
  int i;
  
  for(i=0; (i<n); i++)
    if (matelmt_cmp(map[i].code, c))
      return i;
  
  return -1;
}

int getcmap(FILE *in,char *fn,t_mapping **map)
{
  int       i,n;
  char      line[STRLEN];
  char      code,desc[STRLEN];
  double    r,g,b;
  t_mapping *m;
  
  if (fgets2(line,STRLEN-1,in) == NULL)
    fatal_error(0,"Not enough lines in colormap file %s"
		"(just wanted to read number of entries)",fn);
  sscanf(line,"%d",&n);
  snew(m,n);
  for(i=0; (i<n); i++) {
    if (fgets2(line,STRLEN-1,in) == NULL)
      fatal_error(0,"Not enough lines in colormap file %s"
		  "(should be %d, found only %d)",fn,n+1,i);
    sscanf(line,"%c%s%lf%lf%lf",&code,desc,&r,&g,&b);
    m[i].code.c1=code;
    m[i].code.c2=0;
    m[i].desc=strdup(desc);
    m[i].rgb.r=r;
    m[i].rgb.g=g;
    m[i].rgb.b=b;
  }
  *map=m;
  
  return n;
}

int readcmap(char *fn,t_mapping **map)
{
  FILE      *in;
  int       n;
  
  in=libopen(fn);
  n=getcmap(in,fn,map);
  fclose(in);
  
  return n;
}

void printcmap(FILE *out,int n,t_mapping map[])
{
  int i;
  
  fprintf(out,"%d\n",n);
  for(i=0; (i<n); i++)
    fprintf(out,"%c%c  %20s  %10g  %10g  %10g\n",
	    map[i].code.c1?map[i].code.c1:' ',
	    map[i].code.c2?map[i].code.c2:' ',
	    map[i].desc,map[i].rgb.r,map[i].rgb.g,map[i].rgb.b);
}

void writecmap(char *fn,int n,t_mapping map[])
{
  FILE *out;
  
  out=ffopen(fn,"w");
  printcmap(out,n,map);
  fclose(out);
}

void do_wmap(FILE *out,int i0,int imax,
	     int nlevels,t_rgb rlo,t_rgb rhi,real lo,real hi)
{
  int  i,nlo;
  real r,g,b;
  
  for(i=0; (i<imax); i++) {
    nlo=nlevels-i;
    r=(nlo*rlo.r+i*rhi.r)/nlevels;
    g=(nlo*rlo.g+i*rhi.g)/nlevels;
    b=(nlo*rlo.b+i*rhi.b)/nlevels;
    fprintf(out,"%c %10.3g %10g  %10g  %10g\n",
	    mapper[i+i0],(nlo*lo+i*hi)/nlevels,r,g,b);
  }
}

char *fgetline(char **line, FILE *in)
{
  static char line0[2048];
  char *fg;

  fg=fgets(line0,2048,in);
  *line=line0;
  trim(*line);

  return fg;
}

void skipstr(char **line)
{
  ltrim(*line);
  while((*line[0] != ' ') && (*line[0] != '\0'))
    (*line)++;
}

char *line2string(char **line)
{
  int i;
  
  if (*line != NULL) {
    while (((*line)[0] != '\"' ) && ( (*line)[0] != '\0' ))
      (*line)++;
		       
    if ((*line)[0] != '\"')
      return NULL;
    (*line)++;
  
      i=0;
    while (( (*line)[i] != '\"' ) && ( (*line)[i] != '\0' ))
      i++;
    
    if ((*line)[i] != '\"')
      *line=NULL;
    else
      (*line)[i] = 0;
  }
    
  return *line;
}

void parsestring(char *line,char *label, char *string)
{
  if (strstr(line,label)) {
    if (strstr(line,label) < strchr(line,'\"')) {
      line2string(&line);
      strcpy(string,line);
    }
  }
}

void read_xpm_entry(FILE *in,t_matrix *mm)
{
  t_mapping *map;
  char *line=NULL,*str,buf[256];
  int i,m,col_len,nch;
  unsigned int r,g,b;
  double u;
  char *fg;
  bool bGetOnWithIt;
  t_xpmelmt c;
  
  mm->title[0]=0;
  mm->legend[0]=0;
  mm->label_x[0]=0;
  mm->label_y[0]=0;
  mm->matrix=NULL;
  mm->axis_x=NULL;
  mm->axis_y=NULL;
  mm->bDiscrete=FALSE;

  while (fgetline(&line,in) && (strncmp(line,"static",6) != 0)) {
    parsestring(line,"title",(mm->title));
    parsestring(line,"legend",(mm->legend));
    parsestring(line,"x-label",(mm->label_x));
    parsestring(line,"y-label",(mm->label_y));
    parsestring(line,"type",buf);
  }
  if (buf[0] && (strcasecmp(buf,"Discrete")==0))
    mm->bDiscrete=TRUE;
   
  if (debug)
    printf("%s %s %s %s\n",mm->title,mm->legend,mm->label_x,mm->label_y);

  if  (strncmp(line,"static",6) != 0)
    fatal_error(0,"Invalid XPixMap\n");
  /* Read sizes */
  bGetOnWithIt=FALSE;
  while (!bGetOnWithIt && fgetline(&line,in)) {
    while (( line[0] != '\"' ) && ( line[0] != '\0' ))
      line++;

    if  ( line[0] == '\"' ) {
      line2string(&line);
      sscanf(line,"%d %d %d %d",&(mm->nx),&(mm->ny),&(mm->nmap),&nch);
      if (nch>2)
	fatal_error(0,"Sorry can only read xpm's with at most 2 caracters per pixel\n");
      bGetOnWithIt=TRUE;
    }
  }
  if (debug)
    printf("mm->nx %d mm->ny %d mm->nmap %d nch %d\n",
	   mm->nx,mm->ny,mm->nmap,nch);
  
  /* Read color map */
  snew(map,mm->nmap);
  m=0;
  while ((m < mm->nmap) && fgetline(&line,in)) {
    line=strchr(line,'\"');

    if  ( line == NULL ) 
      fatal_error(0,"Not enough color map entries");
    else {
      line++;
      /* Read xpm color map entry */
      str=strchr(line+2,'#')+1;
      col_len=0;
      while (isxdigit(str[col_len]))
	col_len++;
      if (col_len==6) {
	sscanf(line,"%s %*c #%2x%2x%2x",
	       buf,&r,&g,&b);
	map[m].code.c1=buf[0];
	if (nch==1)
	  map[m].code.c2=0;
	else
	  map[m].code.c2=buf[1];
	map[m].rgb.r=r/255.0;
	map[m].rgb.g=g/255.0;
	map[m].rgb.b=b/255.0;
      }
      else if (col_len==12) {
	sscanf(line,"%s %*c #%4x%4x%4x",
	       buf,&r,&g,&b);
	map[m].code.c1=buf[0];
	if (nch==1)
	  map[m].code.c2=0;
	else
	  map[m].code.c2=buf[1];
	map[m].rgb.r=r/65535.0;
	map[m].rgb.g=g/65535.0;
	map[m].rgb.b=b/65535.0;
      }
      else
	fatal_error(0,"Unsupported or invalid colormap in X PixMap\n");
      line=strchr(line,'\"');
      line++;
      line2string(&line);
      map[m].desc = strdup(line);
      m++;
    }
  }
  mm->map = map;
  if (debug)
    for(m=0;m<mm->nmap;m++) 
      printf("%c %f %f %f %s\n",map[m].code.c1,map[m].rgb.r,map[m].rgb.g,
	     map[m].rgb.b,map[m].desc);

  /* Read axes, if the are any */ 
  bGetOnWithIt=FALSE;
  do {
    if (strstr(line,"x-axis")) {
      line=strstr(line,"x-axis");
      skipstr(&line);
      if (debug)
	printf("%s\n",line);
      snew(mm->axis_x,mm->nx);
      for(i=0;i<mm->nx;i++) {
	sscanf(line,"%lf",&u);
	mm->axis_x[i] = u;
	skipstr(&line);
      }
    }
    else if (strstr(line,"y-axis")) {
      line=strstr(line,"y-axis");
      skipstr(&line);
      if (debug)
	printf("%s\n",line);
      snew(mm->axis_y,mm->ny);
      for(i=0;i<mm->ny;i++) {
	sscanf(line,"%lf",&u);
	   mm->axis_y[i] = u;
	   skipstr(&line);
      }
    }
  } while ((line[0] != '\"') && fgetline(&line,in));

  if (debug)
    for(i=0;i<mm->nx;i++)
      printf("%f %f\n",mm->axis_x[i],mm->axis_y[i]);

  /* Read matrix */
  snew(mm->matrix,mm->nx);
  for(i=0; i<mm->nx; i++)
    snew(mm->matrix[i],mm->ny);
  m=mm->ny-1;
  do {
    while ((line[0] != '\"') && (line[0] != '\0'))
      line++;
    if (line[0] != '\"')
      fatal_error(0,"Not enough caracters in row %d of the matrix\n",m+1);
    else {
      line++;
      for(i=0; i<mm->nx; i++) {
	c.c1=line[nch*i];
	if (nch==1)
	  c.c2=0;
	else
	  c.c2=line[nch*i+1];
	mm->matrix[i][m]=searchcmap(mm->nmap,mm->map,c);
	}
      m--;
    }
  } while ((m>=0) && fgetline(&line,in));
  if (m>=0)
    fatal_error(0,"Not enough rows in the matrix\n");
}

int read_xpm_matrix(char *fnm,t_matrix **matrix)
{
  FILE *in;
  char *line;
  int nmat;

  in=ffopen(fnm,"r");
  
  nmat=0;
  while (fgetline(&line,in)) {
    if (strstr(line,"/* XPM */")) {
      srenew(*matrix,nmat+1);
      read_xpm_entry(in,&(*matrix)[nmat]);
      nmat++;
    }
  }
  fclose(in);

  if (nmat==0)
    fatal_error(0,"Invalid XPixMap\n");

  return nmat;
}

bool matrix2real(t_matrix *matrix, real ***mat)
{
  t_mapping *map;
  real *rmap;
  int i,j,nmap;
  
  nmap=matrix->nmap;
  map=matrix->map;
  snew(rmap,nmap);
  
  for(i=0; i<nmap; i++)
    if ((map[i].desc==NULL) || (sscanf(map[i].desc,"%g",&(rmap[i]))!=1)) {
      fprintf(stderr,"Could not convert matrix to reals,\n"
	      "color map entry %d has a non-real description: \"%s\"\n",
	      i,map[i].desc);
      return FALSE;
    }
  
  snew(*mat,matrix->nx);
  for(i=0; i<matrix->nx; i++) {
    snew((*mat)[i],matrix->ny);
    for(j=0; j<matrix->ny; j++) 
      (*mat)[i][j]=rmap[matrix->matrix[i][j]];
  }
  
  fprintf(stderr,"Converted a matrix with %d levels to reals\n",nmap);

  return TRUE;
}

void write_xpm_header(FILE *out,
		      char *title,char *legend,char *label_x,char *label_y,
		      bool bDiscrete)
{
  fprintf(out,  "/* XPM */\n");
  fprintf(out,  "/* This matrix is generated by %s. */\n",Program());
  fprintf(out,  "/* title:   \"%s\" */\n",title);
  fprintf(out,  "/* legend:  \"%s\" */\n",legend);
  fprintf(out,  "/* x-label: \"%s\" */\n",label_x);
  fprintf(out,  "/* y-label: \"%s\" */\n",label_y);
  if (bDiscrete) 
    fprintf(out,"/* type:    \"Discrete\" */\n");
  else
    fprintf(out,"/* type:    \"Continuous\" */\n");
}

void write_xpm_map3(FILE *out,int n_x,int n_y,int *nlevels,
		    real lo,real mid,real hi,
		    t_rgb rlo,t_rgb rmid,t_rgb rhi)
{
  int    i,nlo,nmid;
  real   r,g,b,clevels;

  if (*nlevels > NMAP*NMAP) {
    fprintf(stderr,"Warning, too many levels (%d) in matrix, using %u only\n",
	    *nlevels,NMAP*NMAP);
    *nlevels=NMAP*NMAP;
  }
  else if (*nlevels < 2) {
    fprintf(stderr,"Warning, too few levels (%d) in matrix, using 2 instead\n",
	    *nlevels);
    *nlevels=2;
  }   
  if (!((mid > lo) && (mid < hi)))
    fatal_error(0,"Lo: %f, Mid: %f, Hi: %f\n",lo,mid,hi);

  fprintf(out,"static char * gv_xpm[] = {\n");
  fprintf(out,"\"%d %d   %d %u\",\n",n_x,n_y,*nlevels,1 + (*nlevels-1)/NMAP);

  nmid = ((mid-lo)/(hi-lo))*(*nlevels-1);
  if (*nlevels > (2 * nmid))
    clevels = *nlevels-1-nmid;
  else
    clevels = nmid;
  for(i=0; (i<nmid); i++) {
    nlo=nmid-i;
    r=rmid.r+(nlo*(rlo.r-rmid.r)/clevels);
    g=rmid.g+(nlo*(rlo.g-rmid.g)/clevels);
    b=rmid.b+(nlo*(rlo.b-rmid.b)/clevels);
    fprintf(out,"\"%c%c c #%02X%02X%02X \" /* \"%.3g\" */,\n",
	    mapper[i % NMAP],
	    (*nlevels <= NMAP) ? ' ' : mapper[i/NMAP],
	    (unsigned int)round(255*r),
	    (unsigned int)round(255*g),
	    (unsigned int)round(255*b),
	    (nlo*lo+i*mid)/nmid);
     }
  for(i=0; (i<(*nlevels-nmid)); i++) {
    nlo=(*nlevels-1-nmid)-i;
    r=rmid.r+(i*(rhi.r-rmid.r)/clevels);
    g=rmid.g+(i*(rhi.g-rmid.g)/clevels);
    b=rmid.b+(i*(rhi.b-rmid.b)/clevels);
    fprintf(out,"\"%c%c c #%02X%02X%02X \" /* \"%.3g\" */,\n",
	    mapper[(i+nmid) % NMAP],
	    (*nlevels <= NMAP) ? ' ' : mapper[(i+nmid)/NMAP],
	    (unsigned int)round(255*r),
	    (unsigned int)round(255*g),
	    (unsigned int)round(255*b),
	    (nlo*mid+i*hi)/(*nlevels-1-nmid));
  }
}


void write_xpm_map(FILE *out,int n_x, int n_y,int *nlevels,real lo,real hi,
		   t_rgb rlo,t_rgb rhi)
{
  int    i,nlo;
  real   invlevel,r,g,b;

  if (*nlevels > NMAP*NMAP) {
    fprintf(stderr,"Warning, too many levels (%d) in matrix, using %u only\n",
	    *nlevels,NMAP*NMAP);
    *nlevels=NMAP*NMAP;
  }
  else if (*nlevels < 2) {
    fprintf(stderr,"Warning, too few levels (%d) in matrix, using 2 instead\n",*nlevels);
    *nlevels=2;
  }

  fprintf(out,"static char * gv_xpm[] = {\n");
  fprintf(out,"\"%d %d   %d %d\",\n",n_x,n_y,*nlevels,1);

  invlevel=1.0/(*nlevels-1);
  for(i=0; (i<*nlevels); i++) {
    nlo=*nlevels-1-i;
    r=(nlo*rlo.r+i*rhi.r)*invlevel;
    g=(nlo*rlo.g+i*rhi.g)*invlevel;
    b=(nlo*rlo.b+i*rhi.b)*invlevel;
    fprintf(out,"\"%c%c c #%02X%02X%02X \" /* \"%.3g\" */,\n",
	    mapper[i % NMAP],(*nlevels <= NMAP) ? ' ' : mapper[i/NMAP],
	    (unsigned int)round(255*r),
	    (unsigned int)round(255*g),
	    (unsigned int)round(255*b),
	    (nlo*lo+i*hi)*invlevel);
  }
}

void write_xpm_axis(FILE *out, char *axis, int n, real *label)
{
  int i;

  if (label) {
    fprintf(out,"/* %s-axis:  ",axis);
    for(i=0;i<n;i++) {
      fprintf(out,"%g ",label[i]);
    }
    fprintf(out,"*/\n");
  }
}

void write_xpm_data(FILE *out, int n_x, int n_y, real **matrix, 
		    real lo, real hi, int nlevels)
{
  int i,j,c;
  real invlevel;

  invlevel=(nlevels-1)/(hi-lo);
  for(j=n_y-1; (j>=0); j--) {
    fprintf(out,"\"");
    for(i=0; (i<n_x); i++) {
      c=round((matrix[i][j]-lo)*invlevel);
      if (c<0) c=0;
      if (c>=nlevels) c=nlevels-1;
      if (nlevels <= NMAP)
	fprintf(out,"%c",mapper[c]);
      else
	fprintf(out,"%c%c",mapper[c % NMAP],mapper[c / NMAP]);
    }
    if (j > 0)
      fprintf(out,"\",\n");
    else
      fprintf(out,"\"\n");
  }
}

void write_xpm_m(FILE *out, t_matrix m)
{
  /* Writes a t_matrix struct to .xpm file */ 
     
  int       i,j;
  bool      bOneChar;
  t_xpmelmt c;

  bOneChar=(m.map[0].code.c2 == 0);
  write_xpm_header(out,m.title,m.legend,m.label_x,m.label_y,
		   m.bDiscrete);
  fprintf(out,"static char * gv_xpm[] = {\n");
  fprintf(out,"\"%d %d   %d %d\",\n",m.nx,m.ny,m.nmap,bOneChar ? 1 : 2);
  for(i=0; (i<m.nmap); i++) 
    fprintf(out,"\"%c%c c #%02X%02X%02X \" /* \"%s\" */,\n",
	    m.map[i].code.c1,
	    bOneChar ? ' ' : m.map[i].code.c2,
	    (unsigned int)round(m.map[i].rgb.r*255),
	    (unsigned int)round(m.map[i].rgb.g*255),
	    (unsigned int)round(m.map[i].rgb.b*255),m.map[i].desc);
  write_xpm_axis(out,"x",m.nx,m.axis_x);
  write_xpm_axis(out,"y",m.ny,m.axis_y);
  for(j=m.ny-1; (j>=0); j--) {
    fprintf(out,"\"");
    if (bOneChar)
      for(i=0; (i<m.nx); i++)
	fprintf(out,"%c",m.map[m.matrix[i][j]].code.c1);
    else
      for(i=0; (i<m.nx); i++) {
	c=m.map[m.matrix[i][j]].code;
	fprintf(out,"%c%c",c.c1,c.c2);
      }
    if (j > 0)
      fprintf(out,"\",\n");
    else
      fprintf(out,"\"\n");
  }
}

void write_xpm3(FILE *out,
		char *title,char *legend,char *label_x,char *label_y,
		int n_x,int n_y,real axis_x[],real axis_y[],
		real *matrix[],real lo,real mid,real hi,
		t_rgb rlo,t_rgb rmid,t_rgb rhi,int *nlevels)
{
  /* See write_xpm.
   * Writes a colormap varying as rlo -> rmid -> rhi.
   */

  if (hi <= lo) 
    fatal_error(0,"hi (%g) <= lo (%g)",hi,lo);

  write_xpm_header(out,title,legend,label_x,label_y,FALSE);
  write_xpm_map3(out,n_x,n_y,nlevels,lo,mid,hi,rlo,rmid,rhi);
  write_xpm_axis(out,"x",n_x,axis_x);
  write_xpm_axis(out,"y",n_y,axis_y);
  write_xpm_data(out,n_x,n_y,matrix,lo,hi,*nlevels);
}

void write_xpm(FILE *out,
	       char *title,char *legend,char *label_x,char *label_y,
	       int n_x,int n_y,real axis_x[],real axis_y[],
	       real *matrix[],real lo,real hi,
	       t_rgb rlo,t_rgb rhi,int *nlevels)
{
  /* out        xpm file
   * title      matrix title
   * legend     label for the continuous legend
   * label_x    label for the x-axis
   * label_y    label for the y-axis
   * n_x, n_y   size of the matrix
   * axis_x[]   the x-ticklabels
   * axis_y[]   the y-ticklables
   * *matrix[]  element x,y is matrix[x][y]
   * lo         output lower than lo is set to lo
   * hi         output higher than hi is set to hi
   * rlo        rgb value for level lo
   * rhi        rgb value for level hi
   * nlevels    number of levels for the output minus one
   */

  if (hi <= lo) 
    fatal_error(0,"hi (%g) <= lo (%g)",hi,lo);

  write_xpm_header(out,title,legend,label_x,label_y,FALSE);
  write_xpm_map(out,n_x,n_y,nlevels,lo,hi,rlo,rhi);
  write_xpm_axis(out,"x",n_x,axis_x);
  write_xpm_axis(out,"y",n_y,axis_y);
  write_xpm_data(out,n_x,n_y,matrix,lo,hi,*nlevels);
}

