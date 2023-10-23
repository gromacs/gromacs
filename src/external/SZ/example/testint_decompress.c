/**
 *  @file test_decompress.c
 *  @author Sheng Di
 *  @date Aug, 2017
 *  @brief This is an example of using Decompression interface.
 *  (C) 2017 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sz.h"
#include "rw.h"

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;

void assessDeCompressionData(int dataType, char* zipFilePath, void* decompressedData, size_t nbEle);

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
    size_t nbEle;
    char zipFilePath[640], outputFilePath[650];
    if(argc < 2)
    {
		printf("Test case: testint_decompress [datatype(-i8/-i16/-i32/-i64/-ui8/-ui16/-ui32/-ui64)] [srcFilePath] [dimension sizes...]\n");
		printf("Example: testint_decompress -i32 testdata/x86/testint32_8x8x8.dat.sz 8 8 8\n");
		exit(0);
	}	
   
    int dataType = SZ_INT32;
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
		printf("Test case: testint_decompress [datatype(-i8/-i16/-i32/-i64)] [data_file]\n");
		printf("Example: testint_decompress -i32 testdata/x86/testint32_8x8x8.dat.sz 8 8 8\n");
		exit(0);		
	}    
    
    sprintf(zipFilePath, "%s", argv[2]);
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
    
    if(r2==0)
		nbEle = r1;
    else if(r3==0)
		nbEle = r1*r2;
    else if(r4==0) 
		nbEle = r1*r2*r3;
    else if(r5==0)
		nbEle = r1*r2*r3*r4;
    else
		nbEle = r1*r2*r3*r4*r5;

    sprintf(outputFilePath, "%s.out", zipFilePath);
    
    size_t byteLength; 
    int status;
    unsigned char *bytes = readByteData(zipFilePath, &byteLength, &status);
    if(status!=SZ_SCES)
    {
        printf("Error: %s cannot be read!\n", zipFilePath);
        exit(0);
    }
  
    //printf("r1=%d,r2=%d,r3=%d,r4=%d,r5=%d\n", r1,r2,r3,r4,r5);
 
	if(dataType == SZ_INT8)
	{
		cost_start();
		uint8_t *data = SZ_decompress(SZ_INT8, bytes, byteLength, r5, r4, r3, r2, r1);
		cost_end();
		free(bytes); 
		
		if(status!=SZ_SCES)
		{
			printf("Error: %s cannot be written!\n", outputFilePath);
			exit(0);
		}
				
		writeByteData(data, nbEle, outputFilePath, &status);	
		assessDeCompressionData(dataType, zipFilePath, data, nbEle);
		free(data);			
	}
	else if(dataType == SZ_INT16)
	{
		cost_start();
		int16_t *data = NULL;
		int16_t *data_ = SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		if(confparams_dec->sol_ID==SZ_Transpose)
			data = detransposeData(data_, dataType, r5, r4, r3, r2, r1);
		else //confparams_dec->sol_ID==SZ
			data = data_;		
		cost_end();
		free(bytes); 
		
		if(status!=SZ_SCES)
		{
			printf("Error: %s cannot be written!\n", outputFilePath);
			exit(0);
		}
				
		writeShortData_inBytes(data, nbEle, outputFilePath, &status);	
		assessDeCompressionData(dataType, zipFilePath, data, nbEle);	
		free(data);					
	}
	else if(dataType == SZ_INT32)
	{
		cost_start();
		int32_t *data = SZ_decompress(SZ_INT32, bytes, byteLength, r5, r4, r3, r2, r1);
		cost_end();
		free(bytes); 
		
		if(status!=SZ_SCES)
		{
			printf("Error: %s cannot be written!\n", outputFilePath);
			exit(0);
		}
				
		writeIntData_inBytes(data, nbEle, outputFilePath, &status);	
		assessDeCompressionData(dataType, zipFilePath, data, nbEle);			
		free(data);			
	}
	else if(dataType == SZ_INT64)
	{
		cost_start();
		int64_t *data = SZ_decompress(SZ_INT64, bytes, byteLength, r5, r4, r3, r2, r1);
		cost_end();
		free(bytes); 
		
		if(status!=SZ_SCES)
		{
			printf("Error: %s cannot be written!\n", outputFilePath);
			exit(0);
		}
				
		writeLongData_inBytes(data, nbEle, outputFilePath, &status);	
		assessDeCompressionData(dataType, zipFilePath, data, nbEle);	
		free(data);			
	}
	else if(dataType == SZ_UINT8)
	{
		cost_start();
		uint8_t *data = SZ_decompress(SZ_UINT8, bytes, byteLength, r5, r4, r3, r2, r1);
		cost_end();
		free(bytes); 
		
		if(status!=SZ_SCES)
		{
			printf("Error: %s cannot be written!\n", outputFilePath);
			exit(0);
		}
				
		writeByteData(data, nbEle, outputFilePath, &status);	
		assessDeCompressionData(dataType, zipFilePath, data, nbEle);		
		free(data);		
	}
	else if(dataType == SZ_UINT16)
	{
		cost_start();
		uint16_t *data = NULL;
		uint16_t *data_ = SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		if(confparams_dec->sol_ID==SZ_Transpose)
			data = detransposeData(data_, dataType, r5, r4, r3, r2, r1);
		else //confparams_dec->sol_ID==SZ
			data = data_;
		cost_end();
		free(bytes); 
		
		if(status!=SZ_SCES)
		{
			printf("Error: %s cannot be written!\n", outputFilePath);
			exit(0);
		}
				
		writeUShortData_inBytes(data, nbEle, outputFilePath, &status);	
		assessDeCompressionData(dataType, zipFilePath, data, nbEle);
		free(data);	
	}
	else if(dataType == SZ_UINT32)
	{
		cost_start();
		uint32_t *data = SZ_decompress(SZ_UINT32, bytes, byteLength, r5, r4, r3, r2, r1);
		cost_end();
		free(bytes); 
		
		if(status!=SZ_SCES)
		{
			printf("Error: %s cannot be written!\n", outputFilePath);
			exit(0);
		}
				
		writeUIntData_inBytes(data, nbEle, outputFilePath, &status);	
		assessDeCompressionData(dataType, zipFilePath, data, nbEle);				
		free(data);
	}
	else if(dataType == SZ_UINT64)
	{
		cost_start();
		uint64_t *data = SZ_decompress(SZ_UINT64, bytes, byteLength, r5, r4, r3, r2, r1);
		cost_end();
		free(bytes); 
		
		if(status!=SZ_SCES)
		{
			printf("Error: %s cannot be written!\n", outputFilePath);
			exit(0);
		}
				
		writeULongData_inBytes(data, nbEle, outputFilePath, &status);	
		assessDeCompressionData(dataType, zipFilePath, data, nbEle);	
		free(data);
	}	
	
    printf("timecost=%f\n",totalCost); 
    printf("done\n");
    
    SZ_Finalize();
   
    return 0;
}

/**
 * Assess the compression error..
 * 
 * */
void assessDeCompressionData(int dataType, char* zipFilePath, void* decompressedData, size_t nbEle)
{
	size_t i, totalNbEle;
	int status;
    char oriFilePath[640];
    strcpy(oriFilePath, zipFilePath);
    oriFilePath[strlen(zipFilePath)-3] = '\0';
	int64_t *data = (int64_t*)malloc(sizeof(int64_t)*nbEle);//decompressed data
	int64_t *ori_data = (int64_t*)malloc(sizeof(int64_t)*nbEle); //original data
	
	if(dataType==SZ_INT8)
    {
		uint8_t *oData = readByteData(oriFilePath, &totalNbEle, &status);
		if(status!=SZ_SCES)
		{
			printf("Error: %s cannot be read!\n", oriFilePath);
			exit(0);
		}
		int8_t* data_ = (int8_t*)decompressedData;    		
		for(i=0;i<nbEle;i++)
		{	
			ori_data[i] = (int8_t)oData[i];	
			data[i] = data_[i];
			//printf("data[%d]=%d %d\n", i, ori_data[i], data[i]);
		}
		
	}
    else if(dataType==SZ_INT16)
    {
		int16_t *oData = readInt16Data(oriFilePath, &totalNbEle, &status);
		if(status!=SZ_SCES)
		{
			printf("Error: %s cannot be read!\n", oriFilePath);
			exit(0);
		}    	
		int16_t* data_ = (int16_t*)decompressedData; 	
		for(i=0;i<nbEle;i++)
		{
			ori_data[i] = oData[i];
			data[i] = data_[i];
		}
		free(oData);
	}
    else if(dataType==SZ_INT32)
    {
		int32_t *oData = readInt32Data(oriFilePath, &totalNbEle, &status);
		if(status!=SZ_SCES)
		{
			printf("Error: %s cannot be read!\n", oriFilePath);
			exit(0);
		}
		int32_t* data_ = (int32_t*)decompressedData;
		for(i=0;i<nbEle;i++)
		{
			ori_data[i] = oData[i];
			data[i] = data_[i];
		}
		free(oData);
	}
    else if(dataType==SZ_INT64)
    {
		free(ori_data);
		int64_t *oData = readInt64Data(oriFilePath, &totalNbEle, &status);
		ori_data = oData;
		int64_t* data_ = (int64_t*)decompressedData;
		for(i=0;i<nbEle;i++)
			data[i] = data_[i];
	}
    else if(dataType==SZ_UINT8)
    {
		uint8_t *oData = readByteData(oriFilePath, &totalNbEle, &status);
		if(status!=SZ_SCES)
		{
			printf("Error: %s cannot be read!\n", oriFilePath);
			exit(0);
		}    		
		uint8_t* data_ = (uint8_t*)decompressedData;		
		for(i=0;i<nbEle;i++)
		{
			ori_data[i] = oData[i];
			data[i] = data_[i];
		}
		free(oData);			
	}
    else if(dataType==SZ_UINT16)
    {
		uint16_t *oData = readUInt16Data(oriFilePath, &totalNbEle, &status);
		if(status!=SZ_SCES)
		{
			printf("Error: %s cannot be read!\n", oriFilePath);
			exit(0);
		}
		uint16_t* data_ = (uint16_t*)decompressedData;
		for(i=0;i<nbEle;i++)
		{
			ori_data[i] = oData[i];
			data[i] = data_[i];
		}
		free(oData);		
	}
    else if(dataType==SZ_UINT32)
    {
		uint32_t *oData = readUInt32Data(oriFilePath, &totalNbEle, &status);
		if(status!=SZ_SCES)
		{
			printf("Error: %s cannot be read!\n", oriFilePath);
			exit(0);
		}
		uint32_t* data_ = (uint32_t*)decompressedData;
		for(i=0;i<nbEle;i++)
		{
			ori_data[i] = oData[i];
			data[i] = data_[i];	
		}
		free(oData);		
	}
    else if(dataType==SZ_UINT64)
    {
		uint64_t *oData = readUInt64Data(oriFilePath, &totalNbEle, &status);
		if(status!=SZ_SCES)
		{
			printf("Error: %s cannot be read!\n", oriFilePath);
			exit(0);
		}  
		uint64_t* data_ = (uint64_t*)decompressedData;
		for(i=0;i<nbEle;i++)
		{
			ori_data[i] = (int64_t)oData[i];			
			data[i] = data_[i];
		}
		free(oData);		
	}	
    
    long Max = 0, Min = 0, diffMax = 0;
    Max = ori_data[0];
    Min = ori_data[0];
    diffMax = llabs(data[0] - ori_data[0]);
    double sum1 = 0, sum2 = 0;
    for (i = 0; i < nbEle; i++)
    {
        sum1 += ori_data[i];
		sum2 += data[i];
    }
    double mean1 = sum1/nbEle;
    double mean2 = sum2/nbEle;

    double sum3 = 0, sum4 = 0;
    double sum = 0, prodSum = 0, relerr = 0;
   
    double maxpw_relerr = 0; 
    for (i = 0; i < nbEle; i++)
    {
        if (Max < ori_data[i]) Max = ori_data[i];
        if (Min > ori_data[i]) Min = ori_data[i];
        
        float err = llabs(data[i] - ori_data[i]);
        //printf("%d: %f, ori=%d, dec=%d\n", i, err, ori_data[i], data[i]);
		if(ori_data[i]!=0)
		{
			relerr = err/ori_data[i];
			if(maxpw_relerr<relerr)
				maxpw_relerr = relerr;
		}
		/*if(relerr>0.001)
		{
			printf("%d %d: err=%.20G ori=%.20G dec=%.20G\n", k, i, err, ori_data[i], data[i]);
			break;
		}*/
		if (diffMax < err)
			diffMax = err;
        prodSum += (ori_data[i]-mean1)*(data[i]-mean2);
        sum3 += (ori_data[i] - mean1)*(ori_data[i]-mean1);
        sum4 += (data[i] - mean2)*(data[i]-mean2);
		sum += err*err;	
    }
    double std1 = sqrt(sum3/nbEle);
    double std2 = sqrt(sum4/nbEle);
    double ee = prodSum/nbEle;
    double acEff = ee/std1/std2;
 
    double mse = sum/nbEle;
    double range = (long)Max - Min;
    double psnr = 20*log10(range)-10*log10(mse);
    double nrmse = sqrt(mse)/range;

    printf ("Min=%ld, Max=%ld, range=%f\n", Min, Max, range);
    printf ("Max absolute error = %ld\n", diffMax);
    printf ("Max relative error = %f\n", ((float)diffMax)/(Max-Min));
    printf ("Max pw relative error = %f\n", maxpw_relerr);
    printf ("PSNR = %f, NRMSE= %.20G\n", psnr,nrmse);
    printf ("acEff=%f\n", acEff);	
    
    free(ori_data);
    free(data);
}
