#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
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

	if(argc<3)
	{
		printf("Usage: convertBinToHDF5 [varName] [infile] [r1, r2, r3, ....]\n");
		exit(0);
	}

	strcpy(varName, argv[1]);
	strcpy(infile, argv[2]);
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

	int dim = computeDimension(r5, r4, r3, r2, r1);
	size_t nbEle = computeDataLength(r5, r4, r3, r2, r1);

	sprintf(outfile, "%s.h5", infile);
	
	FILE *f;
	f = fopen(infile, "rb");
	float *data = (float*)malloc(nbEle*sizeof(float));
	fread(data, sizeof(float), nbEle, f);
	fclose(f);

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
	sprintf(database, "/%s", varName);
	dataset_id = H5Dcreate2(file_id, database, H5T_IEEE_F32LE, dataspace_id, 
						  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(dataset_id, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
					 data);

	/* End access to the dataset and release resources used by it. */
	status = H5Dclose(dataset_id);

	/* Terminate access to the data space. */ 
	status = H5Sclose(dataspace_id);

	/* Close the file. */
	status = H5Fclose(file_id);
}
