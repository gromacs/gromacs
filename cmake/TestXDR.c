#include<rpc/rpc.h>
#include<rpc/xdr.h>

int
main()
{
  /* This should only compile, not run, so set xd to NULL */
  XDR *xd = NULL;
  float f; 
  xdr_float(xd,&f);
  return 0;
}    
