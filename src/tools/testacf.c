#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "xvgr.h"
#include "gstat.h"
#include "copyrite.h"
#include "macros.h"
#include "random.h"
#include "smalloc.h"

int main(int argc,char *argv[])
{
  FILE *fp;
  static char *desc[] = {
    "testac tests the functioning of the GROMACS acf routines"
  };
  static int nframes = 1024;
  static int datatp  = 0;
  static real a=0.02*M_PI;
  t_pargs pa[] = {
    { "-np", FALSE, etINT, &nframes,
      "Number of data points" },
    { "-dtp",FALSE, etINT, &datatp,
      "Which data: 0=all 0.0, 1=all 1.0, 2=cos(a t), 3=random, 4=cos(a t)+random, 5=sin(a t)/(a t)" }
  };
  static char *str[] = {
    "all 0.0", 
    "all 1.0",
    "cos(a t)",
    "random", 
    "cos(a t)+random",
    "sin(a t)/(a t)"
  };
  t_filenm fnm[] = {
    { efXVG, "-d", "acf-data", ffWRITE },
    { efXVG, "-c", "acf-corr", ffWRITE },
    { efXVG, "-comb", "acf-comb.xvg", ffWRITE }
  };
#define NFILE asize(fnm)
  int     npargs,i,nlag;
  int     seed=1993;
  real    *data,*data2,x;
  t_pargs *ppa;
  
  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  ppa    = add_acf_pargs(&npargs,pa);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,TRUE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL);
  snew(data,nframes);
  snew(data2,nframes);
  
  fp = xvgropen(opt2fn("-d",NFILE,fnm),"testac","x","y");
  for(i=0; (i<nframes); i++) {
    x = a*i;
    switch (datatp) {
    case 1:
      data[i] = 1;
      break;
    case 2:
      data[i] = cos(x);
      break;
    case 3:
      data[i] = 2*rando(&seed)-1.0;
      break;
    case 4:
      data[i] = cos(x)+2*rando(&seed)-1.0;
      break;
    case 5:
      if (i==0)
	data[i] = 1;
      else
	data[i] = sin(x)/(x);
    default:
      /* Data remains 0.0 */
      break;
    }
    fprintf(fp,"%10g  %10g\n",x,data[i]);
    data2[i] = data[i];
  }
  fclose(fp);
  
  do_autocorr(opt2fn("-c",NFILE,fnm),str[datatp],
	      nframes,1,&data,a,eacNormal,TRUE,NULL,NULL);
	      
  nlag = get_acflag();
  fp = xvgropen(opt2fn("-comb",NFILE,fnm),"testac","x","y");
  for(i=0; (i<nlag); i++) {
    fprintf(fp,"%10g  %10g  %10g\n",a*i,data2[i],data[i]);
  }
  ffclose(fp);

  xvgr_file(opt2fn("-comb",NFILE,fnm),"-nxy");
    
  thanx(stdout);

  return 0;
}
