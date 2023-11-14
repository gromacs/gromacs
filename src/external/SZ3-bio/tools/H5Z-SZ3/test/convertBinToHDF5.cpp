#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <cstdint>
#include "hdf5.h"


int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
        int dimension;
        if(r1==0)
        {
                dimension = 0;
        }
        else if(r2==0)
        {
                dimension = 1;
        }
        else if(r3==0)
        {
                dimension = 2;
        }
        else if(r4==0)
        {
                dimension = 3;
        }
        else if(r5==0)
        {
                dimension = 4;
        }
        else
        {
                dimension = 5;
        }
        return dimension;
}


size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
        size_t dataLength;
        if(r1==0)
        {
                dataLength = 0;
        }
        else if(r2==0)
        {
                dataLength = r1;
        }
        else if(r3==0)
        {
                dataLength = r1*r2;
        }
        else if(r4==0)
        {
                dataLength = r1*r2*r3;
        }
        else if(r5==0)
        {
                dataLength = r1*r2*r3*r4;
        }
        else
        {
                dataLength = r1*r2*r3*r4*r5;
        }
        return dataLength;
}


int main(int argc, char* argv[]) {

	hid_t       file_id, dataset_id, dataspace_id;  /* identifiers */
	hsize_t     dims[3];
	herr_t      status;
	size_t r1 = 0, r2 = 0, r3 = 0, r4 = 0, r5 = 0;

	char varName[100], infile[100], outfile[100], database[100];
	char datatype[100];

	if(argc<5)
	{
		printf("Usage: convertBinToHDF5 [datatype] [varName] [infile] [r1, r2, r3, ....]\n");
		printf("[datatype] options: -f float32\t-d float64\t-i32 int32\t-i64 int64\n");
		exit(0);
	}



	strcpy(datatype, argv[1]);
	strcpy(varName, argv[2]);
	strcpy(infile, argv[3]);

        if(argc>=5)
                r1 = atoi(argv[4]); //8
        if(argc>=6)
                r2 = atoi(argv[5]); //8
        if(argc>=7)
                r3 = atoi(argv[6]); //128
        if(argc>=8)
                r4 = atoi(argv[7]);
        if(argc>=9)
                r5 = atoi(argv[8]);

	int dim = computeDimension(r5, r4, r3, r2, r1);
	size_t nbEle = computeDataLength(r5, r4, r3, r2, r1);

	snprintf(outfile, 100, "%s.h5", infile);

	/* Create a new file using default properties. */
	file_id = H5Fcreate(outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
                                                                               
        /* Create the data space for the dataset. */
        if(dim==1)
        	dims[0] = r1;
        else if(dim==2)
        {
        	dims[0] = r2;
        	dims[1] = r1;
        }
        else if(dim==3)
        {
        	dims[0] = r3;
        	dims[1] = r2;
        	dims[2] = r1;
        }
        else
        {
        	printf("Error: wrong dimension\n");
        	exit(0);
        }
        dataspace_id = H5Screate_simple(dim, dims, NULL);
                                                                               
        /* Create the dataset. */
        snprintf(database, 100, "/%s", varName);

	if(strcmp(datatype, "-f") == 0){

		FILE *f;
		f = fopen(infile, "rb");
		float *data = (float*)malloc(nbEle*sizeof(float));
		fread(data, sizeof(float), nbEle, f);
		fclose(f);

                dataset_id = H5Dcreate2(file_id, database, H5T_IEEE_F32LE, dataspace_id, 
                					  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                                                                                                  
                status = H5Dwrite(dataset_id, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                				 data);

		hsize_t di_ms[H5S_MAX_RANK];
		int ndims;
		ndims = H5Sget_simple_extent_dims(dataspace_id, di_ms, 0);
		printf("NDIMS: %i\n", ndims);

	}
	else if(strcmp(datatype, "-d") == 0){
		FILE *f;                                          		
                f = fopen(infile, "rb");
                double *data = (double*)malloc(nbEle*sizeof(double));
                fread(data, sizeof(double), nbEle, f);
	        fclose(f);

		dataset_id = H5Dcreate2(file_id, database, H5T_IEEE_F64LE, dataspace_id, 
                					  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                                                                                                  
                status = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                				 data);

	}
	else if(strcmp(datatype, "-i32") == 0){
        	FILE *f;                                          		
                f = fopen(infile, "rb");
                int *data = (int*)malloc(nbEle*sizeof(int));
                fread(data, sizeof(int), nbEle, f);
                fclose(f);

		dataset_id = H5Dcreate2(file_id, database, H5T_STD_I32LE, dataspace_id, 
							  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                                                                                                  
		status = H5Dwrite(dataset_id, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                				 data);
       
       }
	else if(strcmp(datatype, "-i64") == 0){
        	FILE *f;                                          		
                f = fopen(infile, "rb");
                int64_t *data = (int64_t*)malloc(nbEle*sizeof(int64_t));
                fread(data, sizeof(int64_t), nbEle, f);
                fclose(f);
        
		dataset_id = H5Dcreate2(file_id, database, H5T_STD_I64LE, dataspace_id, 
	        					  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	                                                                                          
	        status = H5Dwrite(dataset_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
	        				 data);
	
	}
	else{
		printf("Invalid datatype, use ./convertBinToHDF5 to see usage");
		exit(0);
	}


	/* End access to the dataset and release resources used by it. */
	status = H5Dclose(dataset_id);

	/* Terminate access to the data space. */ 
	status = H5Sclose(dataspace_id);

	/* Close the file. */
	status = H5Fclose(file_id);
}
