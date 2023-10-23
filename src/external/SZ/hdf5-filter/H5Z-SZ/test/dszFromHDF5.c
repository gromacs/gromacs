/**
 *  @file dszFromHDF5.c
 *  @author Sheng Di
 *  @date July, 2017
 *  @brief This is an example of using decompression interface (HDF5)
 *  (C) 2017 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>
#include "hdf5.h"
#include "sz.h"
#include "H5Z_SZ.h"

#define DATASET "testdata_compressed"
#define MAX_CHUNK_SIZE 4294967295 //2^32-1

int main(int argc, char * argv[])
{
	int dimSize = 0;
	size_t r5=0,r4=0,r3=0,r2=0,r1=0,nbEle = 0;
	char hdf5FilePath[640], outputFilePath[640];
	hid_t file, dset, dcpl, space_id, dtype; /*Handles*/
	H5Z_filter_t filter_id = 0;
	herr_t status;
	H5T_class_t type_class;
	H5T_sign_t dsign;
	H5T_order_t dorder;

	htri_t avail;
	char filter_name[80];
	unsigned int flags = 0;
	size_t nelmts = 0, dsize;
	unsigned int values_out[7] = {0,0,0,0,0,0,0}; //at most 7 parameters 

	if(argc < 2)
	{
		printf("Test case: dszFromHDF5 [hdf5FilePath]\n");
		printf("Example 1: dszFromHDF5 testdata/x86/testfloat_8_8_128.dat.sz.hdf5\n");
		printf("Example 2: dszFromHDF5 testdata/x86/testint32_8x8x8.dat.sz.hdf5\n");
		exit(0);
	}

	sprintf(hdf5FilePath, "%s", argv[1]);
	sprintf(outputFilePath, "%s.out.h5", hdf5FilePath);

	/*Open the hdf5 file with SZ-compressed data*/
    file = H5Fopen(hdf5FilePath, H5F_ACC_RDONLY, H5P_DEFAULT);
    dset = H5Dopen(file, DATASET, H5P_DEFAULT);
    
    /*Retrieve dataset creation property list.*/
    dcpl = H5Dget_create_plist(dset);
	
    /*Check that filter is not registered with the library yet*/
	avail = H5Zfilter_avail(H5Z_FILTER_SZ);
	if(!avail)
		printf("sz filter is not yet available after the H5Pget_filter call.\n");
	else
		printf("sz filter is available.\n");
	
	space_id = H5Dget_space(dset);	
	nbEle = H5Sget_simple_extent_npoints(space_id);
	
	if((dtype = H5Dget_type(dset)) < 0)
		printf("Error: H5Dget_type(dset) < 0\n");

	/*Read the data using the default properties.*/
	printf("....Reading SZ compressed data .....................\n");

	if((type_class = H5Tget_class(dtype)) < 0)
	{
		printf("Error: H5Tget_class<0\n");
		exit(0);
	}	
	if (0 == (dsize = H5Tget_size(dtype)))
	{
		printf("Error: H5Tget_size==0\n");
		exit(0);		
	}
		
	if((dorder = H5Tget_order(dtype)) < 0)
		printf("Error: H5Tget_order<0\n");

	switch (type_class)
	{
	case H5T_FLOAT:
		if (H5Tequal(dtype, H5T_IEEE_F32BE) == 1 || H5Tequal(dtype, H5T_IEEE_F32LE) == 1
		|| H5Tequal(dtype, H5T_NATIVE_FLOAT) == 1) 
		{
			printf("data type: float\n");
			float* data = (float*)malloc(sizeof(float)*nbEle);		
			if(dorder==H5T_ORDER_LE)		
				status = H5Dread(dset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			else //H5T_ORDER_BE
				status = H5Dread(dset, H5T_IEEE_F32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			/*Print the first 20 data values to check the correctness.*/	
			int i;
			printf("reconstructed data = ");
			for(i=0;i<20;i++)
				printf("%f ", data[i]);	
			printf("\n");		
			free(data);		
        }
		else //64bit: double 
		{
			printf("data type: double\n");
			double* data = (double*)malloc(sizeof(double)*nbEle);
			if(dorder==H5T_ORDER_LE)
				status = H5Dread(dset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			else
				status = H5Dread(dset, H5T_IEEE_F64BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			/*Print the first 10 data values to check the correctness.*/	
			int i;
			printf("reconstructed data = ");
			for(i=0;i<20;i++)
				printf("%f ", data[i]);	
			printf("\n");	
			free(data);						
		}
		break;
	case H5T_INTEGER:
		if (0 > (dsign = H5Tget_sign(dtype)))
		{
			printf("Error in calling H5Tget_sign(type_id)....\n");
			exit(0);
		}
		if(dsign == H5T_SGN_NONE) //unsigned
		{
			if(dsize==1)
			{
				printf("data type: unsigned char\n");
				unsigned char* data = (unsigned char*)malloc(sizeof(unsigned char)*nbEle);		
				if(dorder==H5T_ORDER_LE)	
					status = H5Dread(dset, H5T_STD_U8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
				else
					status = H5Dread(dset, H5T_STD_U8BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);		
				int i;
				printf("reconstructed data = ");
				for(i=0;i<20;i++)
					printf("%d ", data[i]);	
				printf("\n");	
				free(data);								
			}
			else if(dsize==2)
			{
				printf("data type: unsigned short\n");
				unsigned short* data = (unsigned short*)malloc(sizeof(unsigned short)*nbEle);		
				if(dorder==H5T_ORDER_LE)	
					status = H5Dread(dset, H5T_STD_U16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
				else
					status = H5Dread(dset, H5T_STD_U16BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);	
				int i;
				printf("reconstructed data = ");
				for(i=0;i<20;i++)
					printf("%d ", data[i]);	
				printf("\n");	
				free(data);									
			}
			else if(dsize==4)
			{
				printf("data type: unsigned int\n");
				unsigned int* data = (unsigned int*)malloc(sizeof(unsigned int)*nbEle);		
				if(dorder==H5T_ORDER_LE)	
					status = H5Dread(dset, H5T_STD_U32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
				else
					status = H5Dread(dset, H5T_STD_U32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);		
				int i;
				printf("reconstructed data = ");
				for(i=0;i<20;i++)
					printf("%d ", data[i]);	
				printf("\n");	
				free(data);								
			}
			else if(dsize==8)
			{
				printf("data type: unsigned long\n");
				unsigned long* data = (unsigned long*)malloc(sizeof(unsigned long)*nbEle);		
				if(dorder==H5T_ORDER_LE)	
					status = H5Dread(dset, H5T_STD_U64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
				else
					status = H5Dread(dset, H5T_STD_U64BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);	
				int i;
				printf("reconstructed data = ");
				for(i=0;i<20;i++)
					printf("%ld ", data[i]);	
				printf("\n");	
				free(data);									
			}
		}
		else
		{
			if(dsize==1)
			{
				printf("data type: char\n");
				char *data = (char*)malloc(sizeof(char)*nbEle);
				if(dorder==H5T_ORDER_LE)	
					status = H5Dread(dset, H5T_STD_I8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
				else
					status = H5Dread(dset, H5T_STD_I8BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
				int i;
				printf("reconstructed data = ");
				for(i=0;i<20;i++)
					printf("%d ", data[i]);	
				printf("\n");	
				free(data);										
			}
			else if(dsize==2)
			{
				printf("data type: short\n");
				short *data = (short*)malloc(sizeof(short)*nbEle);
				if(dorder==H5T_ORDER_LE)	
					status = H5Dread(dset, H5T_STD_I16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
				else
					status = H5Dread(dset, H5T_STD_I16BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);		
				int i;
				printf("reconstructed data = ");
				for(i=0;i<20;i++)
					printf("%d ", data[i]);	
				printf("\n");	
				free(data);
			}
			else if(dsize==4)
			{
				printf("data type: int\n");
				int *data = (int*)malloc(sizeof(int)*nbEle);
				if(dorder==H5T_ORDER_LE)	
					status = H5Dread(dset, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
				else
					status = H5Dread(dset, H5T_STD_I32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);		
				int i;
				printf("reconstructed data = ");
				for(i=0;i<20;i++)
					printf("%d ", data[i]);	
				printf("\n");	
				free(data);								
			}
			else if(dsize==8)
			{
				printf("data type: long\n");
				long *data = (long*)malloc(sizeof(long)*nbEle);
				if(dorder==H5T_ORDER_LE)	
					status = H5Dread(dset, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
				else
					status = H5Dread(dset, H5T_STD_I64BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
				int i;
				printf("reconstructed data = ");
				for(i=0;i<20;i++)
					printf("%ld ", data[i]);	
				printf("\n");	
				free(data);									
			}			
		}		
		
		break;
	default: 
		printf("Error: H5Z-SZ supports only float, double or integers.\n");
		exit(0);
	}
	
	status = H5Pclose(dcpl);
	status = H5Dclose(dset);
	status = H5Fclose(file);
	
	return 0;
}
