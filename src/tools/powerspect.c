#include <smalloc.h>
#include <gmx_fft.h>
#include <gmx_fatal.h>
#include <futil.h>
#include "interf.h"
#include "powerspect.h"

void addtoavgenergy(t_complex *list, real *result, int size, int tsteps){
  int i;
	for (i=0;i<size;i++){
		result[i]+=cabs2(list[i])/tsteps;
	}

} 

 
void powerspectavg(real ***intftab, int tsteps, int xbins, int ybins, char **outfiles){
	//Fourier plans and output;
  gmx_fft_t  fftp;
 t_complex *ftspect1; /* Spatial FFT of interface for each time frame and interface ftint[time,xycoord][0], ftintf[time,xycoord][1] for interface 1 and 2 				respectively */
  t_complex *ftspect2;
  real *pspectavg1; //power -spectrum 1st interface
  real *pspectavg2; // -------------- 2nd interface
  real *temp;
  FILE *datfile1,*datfile2; //data-files with interface data
  int n; //time index
  int fy=ybins/2+1; // number of (symmetric) fourier y elements; 
  int rfl=xbins*fy; //length of real - DFT == Symmetric 2D matrix
  int status;
 
/*Prepare data structures for FFT, with time averaging of power spectrum*/
if ( (status=gmx_fft_init_2d_real(&fftp,xbins,ybins,GMX_FFT_FLAG_NONE) )!=0)
  {
        free(fftp);
        gmx_fatal(status,__FILE__,__LINE__,"Error allocating FFT");
    }

//Initialize arrays
snew(ftspect1,rfl);
snew(ftspect2,rfl);
snew(temp,rfl);
snew(pspectavg1,rfl);
snew(pspectavg2,rfl);

//Fouriertransform directly (no normalization or anything)
//NB! Check carefully indexes here

	for (n=0;n<tsteps;n++){
		gmx_fft_2d_real(fftp,GMX_FFT_REAL_TO_COMPLEX,intftab[0][n],ftspect1);
		gmx_fft_2d_real(fftp,GMX_FFT_REAL_TO_COMPLEX,intftab[1][n],ftspect2);
		
		//Add to average for interface 1 here
		addtoavgenergy(ftspect1,pspectavg1,rfl,tsteps);
		//Add to average for interface 2 here
		addtoavgenergy(ftspect2,pspectavg2,rfl,tsteps);
	}
	//Print out average energy-spectrum to outfiles[0] and outfiles[1];

datfile1 = ffopen(outfiles[0],"w");
datfile2 = ffopen(outfiles[1],"w");

//Filling dat files with spectral data
fprintf(datfile1,"%s\n","kx\t ky\t\tPower(kx,ky)");
fprintf(datfile2,"%s\n","kx\t ky\t\tPower(kx,ky)");
        for(n=0;n<rfl;n++){
		fprintf(datfile1,"%d\t%d\t %8.6f\n",(n / fy),(n % fy),pspectavg1[n]);
		fprintf(datfile2,"%d\t%d\t %8.6f\n",(n /fy),(n % fy),pspectavg2[n]);
	}  
ffclose(datfile1); 
ffclose(datfile2);

free(ftspect1);
free(ftspect2);

}

void powerspectavg_intf(t_interf ***if1, t_interf ***if2, int t, int xb, int yb, char **fnms)
{
	real ***surf;

	int xy=xb*yb;
	int i,n;
	
	snew(surf,2);
	snew(surf[0],t);
	snew(surf[1],t);
	for (n=0;n<t;n++){
		snew(surf[0][n],xy);
		snew(surf[1][n],xy);
		for(i=0;i<xy;i++){
			surf[0][n][i]=if1[n][i]->Z;
			surf[1][n][i]=if2[n][i]->Z;
		}
	}
        powerspectavg(surf, t,xb,yb,fnms);
}
	
