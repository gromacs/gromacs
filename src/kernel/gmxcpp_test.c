#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include "gmxcpp.h"

static void die(char *s)
{
  fprintf(stderr,"%s\n",s);
  exit(1);
}

int main(int argc,char *argv[])
{
  int status,done;
#define BUFLEN 256
  char buf[BUFLEN];
  void *handle;
  
  if (argc < 2) {
    printf("Usage %s filename\n",argv[0]);
    exit(1);
  }
  cpp_debug_off();
  status = cpp_open_file(argv[1],&handle);
  if (status != 0) 
    die(cpp_error(&handle,status));
  do {
    status = cpp_read_line(&handle,BUFLEN,buf);
    done = (status == eCPP_EOF);
    if (!done) {
      if (status != eCPP_OK)
	die(cpp_error(&handle,status));
      else 
	printf("%s\n",buf);
    }
  } while (!done);
  status = cpp_close_file(&handle);
  if (status != eCPP_OK) 
    die(cpp_error(&handle,status));
    
  return 0;
}
