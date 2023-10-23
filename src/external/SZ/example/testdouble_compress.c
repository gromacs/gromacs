/**
 *  @file test_compress.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief This is an example of using compression interface
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <sys/time.h>
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
    size_t r5=0,r4=0,r3=0,r2=0,r1=0;
    char oriFilePath[640], outputFilePath[650];
    char *cfgFile;
    
    if(argc < 3)
    {
	printf("Test case: testdouble_compress [config_file] [srcFilePath] [dimension sizes...]\n");
	printf("Example: testdouble_compress sz.config testdouble_8_8_128.dat 8 8 128\n");
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
 
    int status = 0;
    printf("cfgFile=%s\n", cfgFile); 
    status = SZ_Init(cfgFile);

    /*sz_params sz;
    memset(&sz, 0, sizeof(sz_params));
    sz.sol_ID = SZ;
    sz.sampleDistance = 100;
    sz.quantization_intervals = 0;
    sz.max_quant_intervals = 65536;
    sz.predThreshold = 0.98;
    sz.szMode = SZ_BEST_COMPRESSION;
    sz.losslessCompressor = ZSTD_COMPRESSOR;
    sz.gzipMode = 1;
    sz.errorBoundMode = REL;
    sz.absErrBound = 1E-6;
    sz.relBoundRatio = 1E-5;

    SZ_Init_Params(&sz);*/

    sprintf(outputFilePath, "%s.sz", oriFilePath);
   
    size_t nbEle = 0;
    double *data = readDoubleData(oriFilePath, &nbEle, &status);
    if(status!=SZ_SCES)
    {
	printf("Error: file %s cannot be read!\n", oriFilePath);
	exit(0);
    }
   
    size_t outSize;
    cost_start(); 
    unsigned char *bytes = SZ_compress(SZ_DOUBLE, data, &outSize, r5, r4, r3, r2, r1);
    //char *bytes = (char *)malloc(nbEle*sizeof(double)); //
    //char* bytes = SZ_compress_args(SZ_DOUBLE, data, &outSize, ABS, 1E-12, 0.000001, 0.01, r5, r4, r3, r2, r1);
    cost_end();
    printf("timecost=%f\n",totalCost);

    writeByteData(bytes, outSize, outputFilePath, &status);
    if(status!=SZ_SCES)
    {
	printf("Error: file %s cannot be written!\n", outputFilePath);
	free(data);
	exit(0);
    }
    free(data);
    free(bytes);
    printf("done\n");
    
    SZ_Finalize();
    
    return 0;
}
