/**
 *  @file test_compress_ts.c
 *  @author Sheng Di
 *  @date May, 2018
 *  @brief This is an example of using compression interface
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include "sz.h"
#include "rw.h"

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;


void cost_start()
{
	totalCost = 0;
        gettimeofday(&costStart, NULL);
}

void cost_end()
{
        double elapsed;
        struct timeval costEnd;
        gettimeofday(&costEnd, NULL);
        elapsed = ((costEnd.tv_sec*1000000+costEnd.tv_usec)-(costStart.tv_sec*1000000+costStart.tv_usec))/1000000.0;
        totalCost += elapsed;
}


int main(int argc, char * argv[])
{
    int i = 0;
    size_t r5=0,r4=0,r3=0,r2=0,r1=0;
    char cmprFilePath[640], outputDir[640], outputFilePath[600];
    int status = 0;
    
    if(argc < 3)
    {
		printf("Test case: testdouble_decompress_ts [srcDir] [dimension sizes...]\n");
		printf("Example: testdouble_decompress_ts /home/sdi/Data/Hurricane-ISA/consecutive-steps 500 500 100\n");
		exit(0);
    }
  
    sprintf(outputDir, "%s", argv[1]);
    if(strcmp(outputDir, "sz.config")==0)
    {
    	printf("Error: wrong input\n");
	printf("Test case: testdouble_decompress_ts [srcDir] [dimension sizes...]\n");
	exit(0);
    } 
    if(argc>=3)
		r1 = atoi(argv[2]); //8
    if(argc>=4)
		r2 = atoi(argv[3]); //8
    if(argc>=5)
		r3 = atoi(argv[4]); //128
    if(argc>=6)
        r4 = atoi(argv[5]);
    if(argc>=7)
        r5 = atoi(argv[6]);
      
    char oriFilePath[600];
    size_t byteLen = 0;
    size_t dataLength = computeDataLength(r5,r4,r3,r2,r1);
    double *data = (double*)malloc(sizeof(double)*dataLength);
    float *data_out = (float*)malloc(sizeof(float)*dataLength);
    SZ_registerVar(1, "dump", SZ_DOUBLE, data, REL, 0, 0.001, 0, r5, r4, r3, r2, r1);

    if(status != SZ_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }
  
    int j = 0;
    for(i=1;i<20;i++)
        {
                printf("simulation time step %d\n", i);
                sprintf(cmprFilePath, "%s/QCLOUDf%02d-double.bin.dat.sz2", outputDir, i);
                printf("compressed data file: %s\n", cmprFilePath);
                unsigned char *bytes = readByteData(cmprFilePath, &byteLen, &status);
                cost_start();
                SZ_decompress_ts(bytes, byteLen);
                cost_end();
                printf("timecost=%f\n",totalCost);
                sprintf(outputFilePath, "%s/QCLOUDf%02d-double.bin.dat.sz2.out", outputDir, i);
                printf("writing decompressed data to %s\n", outputFilePath);
		for(j=0;j<dataLength;j++)
			data_out[j] = (float)data[j];
                writeFloatData_inBytes(data_out, dataLength, outputFilePath, &status);
                free(bytes);
        }


    printf("done\n");
    free(data);
    SZ_Finalize();
    
    return 0;
}
