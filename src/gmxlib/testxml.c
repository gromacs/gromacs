#include "typedefs.h"
#include "macros.h"
#include "smalloc.h"
#include "xmlio.h"
#include "statutil.h"
#include "tpxio.h"

int main(int argc,char *argv[])
{
  int        step,natoms;
  real       t,lambda;
  t_inputrec ir;
  t_topology top;
  matrix     box;
  rvec       *x,*v,*f;
  t_filenm fnm[] = {
    { efTPX, NULL, NULL, ffREAD  },
    { efXML, "-r", NULL, ffREAD  },
    { efXML, "-o", NULL, ffWRITE }
  };  
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,0,NULL,0,NULL,0,NULL);
  
  init_top(&top);
  if (opt2bSet("-r",NFILE,fnm))
    read_xml(opt2fn("-r",NFILE,fnm),&step,&t,&lambda,&ir,
	     box,&natoms,&x,&v,&f,&top);
  else {
    t_tpxheader tpx;
    
    read_tpxheader(ftp2fn(efTPX,NFILE,fnm),&tpx);
    snew(x,tpx.natoms);
    snew(v,tpx.natoms);
    f = NULL;
    read_tpx(ftp2fn(efTPX,NFILE,fnm),&step,&t,&lambda,&ir,
	     box,&natoms,x,v,f,&top);
  }
  /*write_xml(opt2fn("-o",NFILE,fnm),step,t,lambda,&ir,box,natoms,x,v,f,&top);*/
  
  return 0;
}



