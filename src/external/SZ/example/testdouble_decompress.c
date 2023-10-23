/**
 *  @file test_decompress.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief This is an example of using Decompression interface.
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
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
    size_t nbEle, totalNbEle;
    char zipFilePath[640], outputFilePath[650];
    
    if(argc < 2)
    {
		printf("Test case: testdouble_decompress [srcFilePath] [dimension sizes...]\n");
		printf("Example: testdouble_decompress testdouble_8_8_128.dat.sz 8 8 128\n");
		exit(0);
	}	
    sprintf(zipFilePath, "%s", argv[1]);
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
        printf("Error: %s cannot be READ!\n", zipFilePath);
        exit(0);
    }

    cost_start();    
    double *data = SZ_decompress(SZ_DOUBLE, bytes, byteLength, r5, r4, r3, r2, r1);
    cost_end();    
    printf("timecost=%f\n",totalCost);
   
    free(bytes); 
    //int i=0;
    //for(;i<8192;i++)
    //	printf("i=%d, data=%f\n",i,data[i]);
    writeDoubleData_inBytes(data, nbEle, outputFilePath, &status);
    if(status!=SZ_SCES)
    {
        printf("Error: %s cannot be written!\n", outputFilePath);
        exit(0);
    }

    
    printf("done\n");
    
    //SZ_Finalize();
    

    char oriFilePath[650];
    strcpy(oriFilePath, zipFilePath);
    oriFilePath[strlen(zipFilePath)-3] = '\0';
    double *ori_data = readDoubleData(oriFilePath, &totalNbEle, &status);
    if(status!=SZ_SCES)
    {
        printf("Error: %s cannot be read!\n", oriFilePath);
        exit(0);
    }

    size_t i;
    double Max, Min, diffMax, err, maxpw_relerr = 0, relerr;
    Max = ori_data[0];
    Min = ori_data[0];
    diffMax = fabs(data[0] - ori_data[0]);

    for (i = 0; i < nbEle; i++)
    {
    	if (Max < ori_data[i]) Max = ori_data[i];
    	if (Min > ori_data[i]) Min = ori_data[i];
		err = fabs(data[i] - ori_data[i]);
    	if (diffMax < err)
    		diffMax = err;
        if(ori_data[i]!=0)
        {
                relerr = err/fabs(ori_data[i]);
		/*if(relerr>0.00001)
		{
			printf("error:i=%d, ori_data=%f, dec_data=%f\n",i, ori_data[i], data[i]);
			exit(0);
		}*/
                if(maxpw_relerr<relerr)
                        maxpw_relerr = relerr;
        }


	/*if(fabs(data[i] - ori_data[i]) > 1E-1)
	{
		printf("error: i=%d, %.20G, %.20G\n",i,ori_data[i], data[i]);
		exit(0);
	}*/
    }

    printf ("Max absolute error = %.20G\n", diffMax);
    printf ("Max relative error = %.20G\n", diffMax/(Max-Min));
    printf ("Max pw_relative err = %.20G\n", maxpw_relerr);

    free(ori_data);
    free(data);
    return 0;
}
