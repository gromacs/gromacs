#include<rpc/rpc.h>
#include<rpc/xdr.h>

int
main()
{
  XDR *xd; 
  float f; 
  xdr_float(xd,&f);
}    
