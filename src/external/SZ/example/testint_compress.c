/**
 *  @file test_compress.c
 *  @author Sheng Di
 *  @date Aug, 2017
 *  @brief This is an example of using compression interface
 *  (C) 2017 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
    int dataType = SZ_INT32;
    int status; 
    
    if(argc < 4)
    {
		printf("Test case: testint_compress [datatype(-i8/-i16/-i32/-i64/-ui8/-ui16/-ui32/-ui64)] [config_file] [data_file]\n");
		printf("Example: testint_compress -i32 sz.config testdata/x86/testint32_8x8x8.dat 8 8 8\n");
		exit(0);
    }
   
	if(strcmp(argv[1], "-i8")==0)
		dataType = SZ_INT8;
	else if(strcmp(argv[1], "-i16")==0)
		dataType = SZ_INT16;
	else if(strcmp(argv[1], "-i32")==0)
		dataType = SZ_INT32;
	else if(strcmp(argv[1], "-i64")==0)
		dataType = SZ_INT64;	
	else if(strcmp(argv[1], "-ui8")==0)
		dataType = SZ_UINT8;
	else if(strcmp(argv[1], "-ui16")==0)
		dataType = SZ_UINT16;
	else if(strcmp(argv[1], "-ui32")==0)
		dataType = SZ_UINT32;
	else if(strcmp(argv[1], "-ui64")==0)
		dataType = SZ_UINT64;			
	else
	{
		printf("Error: missing/unrecoganized data type: %s. \n", argv[1]);
		printf("Test case: testint_compress [datatype(-i8/-i16/-i32/-i64)] [config_file] [data_file]\n");
		printf("Example: testint_compress -i32 sz.config testdata/x86/testint32_8x8x8.dat 8 8 8\n");
		exit(0);		
	}
    cfgFile=argv[2];

    sprintf(oriFilePath, "%s", argv[3]);
    if(argc>=5)
		r1 = atoi(argv[4]); //8
    if(argc>=6)
		r2 = atoi(argv[5]); //8
    if(argc>=7)
		r3 = atoi(argv[6]); //8
    if(argc>=8)
        r4 = atoi(argv[7]);
    if(argc>=9)
    {
	   r5 = atoi(argv[8]);
	}
	printf("cfgFile=%s\n", cfgFile); 
    status = SZ_Init(cfgFile);
    if(status == SZ_NSCS)
		exit(0);
    sprintf(outputFilePath, "%s.sz", oriFilePath);
   
    size_t nbEle, outSize; 
    unsigned char *bytes = NULL;
	if(dataType==SZ_INT8)
	{
		int8_t *data = (int8_t *) readByteData(oriFilePath, &nbEle, &status);
		if(status != SZ_SCES)
		{
			printf("Error: data file %s cannot be read!\n", oriFilePath);
			exit(0);
		}

		cost_start();
		bytes = SZ_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
		cost_end();
		writeByteData(bytes, outSize, outputFilePath, &status);
		if(status != SZ_SCES)
		{
			printf("Error: data file %s cannot be written!\n", outputFilePath);
			exit(0);
		}

		free(data);				
	}
	else if(dataType==SZ_INT16)
	{
		int16_t *data = readInt16Data(oriFilePath, &nbEle, &status);
		if(status != SZ_SCES)
		{
			printf("Error: data file %s cannot be read!\n", oriFilePath);
			exit(0);
		}

		cost_start();
		if(confparams_cpr->sol_ID==SZ)
			bytes = SZ_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
		else if(confparams_cpr->sol_ID==SZ_Transpose)
		{
			int status = 0;
			bytes = SZ_compress_customize("SZ_Transpose", NULL, dataType, data, r5, r4, r3, r2, r1, &outSize, &status);
		}
		cost_end();
		writeByteData(bytes, outSize, outputFilePath, &status);
		if(status != SZ_SCES)
		{
			printf("Error: data file %s cannot be written!\n", outputFilePath);
			exit(0);
		}

		free(data);				
	}
	else if(dataType==SZ_INT32)
	{
		int32_t *data = readInt32Data(oriFilePath, &nbEle, &status);
		if(status != SZ_SCES)
		{
			printf("Error: data file %s cannot be read!\n", oriFilePath);
			exit(0);
		}

		cost_start();
		bytes = SZ_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
		cost_end();
		writeByteData(bytes, outSize, outputFilePath, &status);
		if(status != SZ_SCES)
		{
			printf("Error: data file %s cannot be written!\n", outputFilePath);
			exit(0);
		}

		free(data);		
	}
	else if(dataType==SZ_INT64)
	{
		int64_t *data = readInt64Data(oriFilePath, &nbEle, &status);
		if(status != SZ_SCES)
		{
			printf("Error: data file %s cannot be read!\n", oriFilePath);
			exit(0);
		}

		cost_start();
		bytes = SZ_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
		cost_end();
		writeByteData(bytes, outSize, outputFilePath, &status);
		if(status != SZ_SCES)
		{
			printf("Error: data file %s cannot be written!\n", outputFilePath);
			exit(0);
		}

		free(data);				
	}
	else if(dataType==SZ_UINT8)
    {
		uint8_t *data = readByteData(oriFilePath, &nbEle, &status);
		if(status != SZ_SCES)
		{
			printf("Error: data file %s cannot be read!\n", oriFilePath);
			exit(0);
		}

		cost_start();
		bytes = SZ_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
		cost_end();
		writeByteData(bytes, outSize, outputFilePath, &status);
		if(status != SZ_SCES)
		{
			printf("Error: data file %s cannot be written!\n", outputFilePath);
			exit(0);
		}

		free(data);				
	}
	else if(dataType==SZ_UINT16)
	{
		uint16_t *data = readUInt16Data(oriFilePath, &nbEle, &status);
		if(status != SZ_SCES)
		{
			printf("Error: data file %s cannot be read!\n", oriFilePath);
			exit(0);
		}

		cost_start();
		if(confparams_cpr->sol_ID==SZ)
			bytes = SZ_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
		else if(confparams_cpr->sol_ID==SZ_Transpose)
		{
			int status = 0;
			bytes = SZ_compress_customize("SZ_Transpose", NULL, dataType, data, r5, r4, r3, r2, r1, &outSize, &status);
		}
		cost_end();
		writeByteData(bytes, outSize, outputFilePath, &status);
		if(status != SZ_SCES)
		{
			printf("Error: data file %s cannot be written!\n", outputFilePath);
			exit(0);
		}

		free(data);				
	}
	else if(dataType==SZ_UINT32)
	{
		uint32_t *data = readUInt32Data(oriFilePath, &nbEle, &status);
		if(status != SZ_SCES)
		{
			printf("Error: data file %s cannot be read!\n", oriFilePath);
			exit(0);
		}

		cost_start();
		bytes = SZ_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
		cost_end();
		writeByteData(bytes, outSize, outputFilePath, &status);
		if(status != SZ_SCES)
		{
			printf("Error: data file %s cannot be written!\n", outputFilePath);
			exit(0);
		}

		free(data);		
	}
	else if(dataType==SZ_UINT64)
	{
		uint64_t *data = readUInt64Data(oriFilePath, &nbEle, &status);
		if(status != SZ_SCES)
		{
			printf("Error: data file %s cannot be read!\n", oriFilePath);
			exit(0);
		}

		cost_start();
		bytes = SZ_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
		cost_end();
		writeByteData(bytes, outSize, outputFilePath, &status);
		if(status != SZ_SCES)
		{
			printf("Error: data file %s cannot be written!\n", outputFilePath);
			exit(0);
		}

		free(data);				
	}
	
	free(bytes); 
	    
    printf("timecost=%f, output compressed file: %s\n",totalCost, outputFilePath);     
    SZ_Finalize();
    printf("done\n");
    
    return 0;
}
