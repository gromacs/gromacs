#include <stdio.h>
#include <ctype.h>
#include "strdb.h"
#include "copyrite.h"
#include "smalloc.h"

void add_quote(char *q)
{
  FILE *fp;
  int  i,n;
  char **str = NULL;
  char *db   = "gurgle.dat";
  
  n = get_strings(db,&str);
  srenew(str,n+1);
  snew(str[n],strlen(q)+1);
  for(i=0; (i<strlen(q)); i++)
    str[n][i] = ~q[i];
  str[n][i] = '\0';
  n++;
  fp = fopen(db,"w");
  fprintf(fp,"%d\n",n);
  for(i=0; (i<n); i++) 
    fprintf(fp,"%s\n",str[i]);
  fclose(fp);
}

int main(int argc,char *argv[])
{
  int  i;
  char c;
  
  for(i=1; (i<argc); i++) {
    do {
      fprintf(stderr,"Add quote '%s' (y/n)? ",argv[i]);
      c = toupper(fgetc(stdin));
    } while ((c != 'Y') && (c != 'N'));
    if (c == 'Y') {
      add_quote(argv[i]);
    }
  }
  thanx(stdout);
  
  return 0;
}
