#include <stdio.h>
#include <stdlib.h>

typedef double matrix[3][3];
typedef float  fmatrix[3][3];

int
main(int argc, char **argv)
{
    int i,j,k,bdouble;
    int nx,ny,nz;
    int d,dd;
    float   fbox[3][3];
    double  box[3][3];
    
    matrix *localpressure;
    matrix *slice;
    fmatrix tmpmat;
    FILE *fp;
    matrix average;

    fp=fopen(argv[1],"r");
    
    
    fread(&bdouble,sizeof(int),1,fp);
    printf("Reading %s precision local pressure file...\n", bdouble ? "double" : "single");

    if(bdouble)
    {
        fread(box,sizeof(double),9,fp);    
    }
    else
    {
        fread(fbox,sizeof(float),9,fp);
        for(i=0;i<3;i++)
            for(j=0;j<3;j++)
                box[i][j]=fbox[i][j];
    }

    fread(&nx,sizeof(int),1,fp);
    fread(&ny,sizeof(int),1,fp);
    fread(&nz,sizeof(int),1,fp);
    
    localpressure=malloc(sizeof(matrix)*nx*ny*nz);
    slice = malloc(sizeof(matrix)*nz);
    
    if(bdouble)
    {
        fread(localpressure,sizeof(matrix),nx*ny*nz,fp);
    }
    else
    {
        for(k=0;k<nx*ny*nz;k++)
        {
            fread(tmpmat,sizeof(tmpmat),1,fp);
            for(i=0;i<3;i++)
                for(j=0;j<3;j++)
                    localpressure[k][i][j]=tmpmat[i][j];
        }
    }
    
    fclose(fp);
    
    for(d=0;d<3;d++)
      for(dd=0;dd<3;dd++)
          average[d][dd]=0;

    for(k=0;k<nz;k++)
    {
        for(d=0;d<3;d++)
            for(dd=0;dd<3;dd++)
                slice[k][d][dd]=0;
        
        for(i=0;i<nx;i++)
            for(j=0;j<ny;j++)
                for(d=0;d<3;d++)
                    for(dd=0;dd<3;dd++)
                        slice[k][d][dd] += localpressure[i*ny*nz+j*nz+k][d][dd]/(nx*ny);

        for(d=0;d<3;d++)
	  for(dd=0;dd<3;dd++)
	    average[d][dd] += slice[k][d][dd]/nz;

        printf("P[z%2d]: %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n",
               k,
               slice[k][0][0],slice[k][0][1],slice[k][0][2],
               slice[k][1][0],slice[k][1][1],slice[k][1][2],               
               slice[k][2][0],slice[k][2][1],slice[k][2][2]);
    }
    printf("Average:  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n",
	   average[0][0],average[0][1],average[0][2],
	   average[1][0],average[1][1],average[1][2],
	   average[2][0],average[2][1],average[2][2]);

    return 0;
}
