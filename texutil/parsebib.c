#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void ltrim (char *str)
{
  char *tr;
  int c;

  if (!str)
    return;

  tr = strdup (str);
  c  = 0;
  while ((tr[c] == ' ') || (tr[c] == '\t'))
    c++;

  strcpy (str,tr+c);
  free (tr);
}

void rtrim (char *str)
{
  int nul;

  if (!str)
    return;

  nul = strlen(str)-1;
  while ((nul > 0) && ((str[nul] == ' ') || (str[nul] == '\t')) ) {
    str[nul] = '\0';
    nul--;
  }
}

void trim (char *str)
{
  ltrim (str);
  rtrim (str);
}

char *fgets2(char *line, int n, FILE *stream)
/* This routine reads a string from stream of max length n
 * and zero terminated, without newlines
 * line should be long enough (>= n)
 */
{
  char *c;
  if (fgets(line,n,stream)==NULL) return NULL;
  if ((c=strchr(line,'\n'))!=NULL) *c=0;
  return line;
}

int main(int argc,char *argv[])
{
  FILE *in,*out;
  char fn[234];
  char buf[1024];
  char ref[1024];
  int  i,j;
  char *aap,*acco;
  char *header;
  
  header=argv[1];
  for(j=2; (j<argc); j++) {
    sprintf(fn,"%s.bib",argv[j]);
    in=fopen(fn,"r");
    sprintf(fn,"%s.tex",argv[j]);
    out=fopen(fn,"w");
    
    fprintf(out,"\\bibliographystyle{plain}\n");
    /*fprintf(out,"\\nocite{TitlesOn}\n");*/
    /* Let's assume we have not more than one entry per line */
    while (fgets2(buf,1023,in) != NULL) {
      if ((aap=strchr(buf,'@')) != NULL) {
	if ((acco=strchr(aap,'{')) != NULL) {
	  acco++;
	for(i=0; (*acco != '\0') && (*acco != ','); i++,acco++)
	  ref[i]=*acco;
	  ref[i]='\0';
	  trim(ref);
	  if (strlen(ref) > 0)
	    fprintf(out,"\\nocite{%s}\n",ref);
	}
      }
    }
    fprintf(out,"\\centerline{\\Huge\\bf %s}\n",header);
    fprintf(out,"\\bibliography{%s}",argv[j]);
    fclose(in);
    fclose(out);
  }
}
