#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sz.h"
#include "rw.h"

float absEB = 1E-4;

int main(int argc, char * argv[])
{
    size_t r5=0,r4=0,r3=0,r2=0,r1=0;
    char oriFilePath[640];
    char *cfgFile;
    
    if(argc < 3)
    {
	printf("Test case: testfloat_CompDecomp_subblock [config_file] [srcFilePath] [dimension sizes...]\n");
	printf("Example: testfloat_CompDecomp_subblock sz.config testfloat_8_8_128.dat 8 8 128\n");
	exit(0);
    }
   
    cfgFile=argv[1];
    sprintf(oriFilePath, "%s", argv[2]);
    if(argc>=4)
	r1 = atoi(argv[3]); //8
    if(argc>=5)
	r2 = atoi(argv[4]); //8
    if(argc>=6)
	r3 = atoi(argv[5]); //128
    if(argc>=7)
        r4 = atoi(argv[6]);
    if(argc>=8)
        r5 = atoi(argv[7]);
   
    printf("cfgFile=%s\n", cfgFile); 

    int status = SZ_Init(cfgFile);
    if(status == SZ_NSCS)
	exit(0);
   
    size_t nbEle;
    float *oriData = NULL, *decompData = NULL;

    oriData = readFloatData(oriFilePath, &nbEle, &status);
    if(status != SZ_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }
  
    size_t outSize; 
    unsigned char *bytes = (unsigned char *)malloc(nbEle*sizeof(float));

    /* Compress a subblock of the original data */
    SZ_compress_args3(SZ_FLOAT, oriData, bytes, &outSize, ABS, absEB, 0,
    		r5, r4, r3, r2, r1, 0, 0, 0, 0, 0, r5/2, r4/2, r3/2, r2/2, r1/64);
    printf ("Subblock data's compression is done.\n");

    /* Decompress the subblock */
    if (r2 == 0)
    	decompData = SZ_decompress(SZ_FLOAT, bytes, outSize, 0, 0, 0, 0, r1/64+1);
    else
    if (r3 == 0)
    	decompData = SZ_decompress(SZ_FLOAT, bytes, outSize, 0, 0, 0, r2/2+1, r1/64+1);
    else
    if (r4 == 0)
    	decompData = SZ_decompress(SZ_FLOAT, bytes, outSize, 0, 0, r3/2+1, r2/2+1, r1/64+1);
    else
    if (r5 == 0)
    	decompData = SZ_decompress(SZ_FLOAT, bytes, outSize, 0, r4/2+1, r3/2+1, r2/2+1, r1/64+1);
    else
		printf("Error: doesn't support 5 dimensions for now.\n");

    printf ("Subblock data's decompression is done.\n");

    float maxDiff = 0;

    size_t i1, i2, i3, i4, i5;
    size_t index1 = 0, index2 = 0;
    for (i5 = 0; i5 <= r5/2; i5++)
        for (i4 = 0; i4 <= r4/2; i4++)
            for (i3 = 0; i3 <= r3/2; i3++)
                for (i2 = 0; i2 <= r2/2; i2++)
                    for (i1 = 0; i1 <= r1/64; i1++)
                    {
                    	index1 = i5*(r4*r3*r2*r1)+i4*(r3*r2*r1)+i3*(r2*r1)+i2*r1+i1;
                    	float data1 = oriData[index1];
                    	float data2 = decompData[index2++];
                    	float diff = fabs(data1-data2);
                    	if (diff > maxDiff)
                    		maxDiff = diff;
                    }

    if (maxDiff <= absEB)
    {
    	printf ("Maximum Absolute Error is %f\n", maxDiff);
    	printf ("Absolute Error bound is %f\n", absEB);
    	printf ("Test passed.\n");
    }

    free(bytes); 
    free(oriData);
    free(decompData);

    SZ_Finalize();
    
    return 0;
}
