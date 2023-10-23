#include <stdio.h>
#include <stdlib.h>
#include "NetCDFReader.h"
#include "rw.h"

void usage()
{
	printf("Options:\n");
	printf("* data type:\n");
	printf("	-f: single precision (float type)\n");
	printf("	-d: double precision (double type)\n");
	printf("* input data file:\n");
	printf("	-i <NetCDF data file> : NetCDF data file\n");
	printf("	-n <dataset path> : dataset path in NetCDF file\n");
	printf("* output data file:\n");
	printf("	-o <output binary file> : binary file in little_endian type\n");
	printf("* dimensions: \n");
	printf("	-1 <nx> : dimension for 1D data such as data[nx]\n");
	printf("	-2 <nx> <ny> : dimensions for 2D data such as data[ny][nx]\n");
	printf("	-3 <nx> <ny> <nz> : dimensions for 3D data such as data[nz][ny][nx] \n");
	printf("	-4 <nx> <ny> <nz> <np>: dimensions for 4D data such as data[np][nz][ny][nx] \n");
	printf("* examples: \n");
	printf("	test_NetCDFReader -f -i testdata/testdata.nc -n data -2 6 12 -o test.dat\n");
	exit(0);
}

int main(int argc, char* argv[])
{
	int dataType;
	char* inPath = NULL;
	char* outPath = NULL;
	char *dataset = NULL;

	size_t i = 0;
	size_t nbEle;

	size_t r5 = 0;
	size_t r4 = 0;
	size_t r3 = 0;
	size_t r2 = 0; 
	size_t r1 = 0;

	if(argc==1)
	  usage();

	for(i=1;i<argc;i++)
	{
		if (argv[i][0] != '-' || argv[i][2])
		  usage();
		switch (argv[i][1])
		{
			case 'o': 
				if (++i == argc)
				  usage();
				outPath = argv[i];	
				break;
			case 'f': 
				dataType = 0;
				break;
			case 'd':
				dataType = 1;
				break;
			case 'i':
				if (++i == argc)
				  usage();
				inPath = argv[i];		
				break;
			case '1': 
				if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1)
				  usage();
				break;
			case '2':
				if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1 ||
							++i == argc || sscanf(argv[i], "%zu", &r2) != 1)
				  usage();
				break;
			case '3':
				if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1 ||
							++i == argc || sscanf(argv[i], "%zu", &r2) != 1 ||
							++i == argc || sscanf(argv[i], "%zu", &r3) != 1)
				  usage();		
				break;
			case '4':
				if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1 ||
							++i == argc || sscanf(argv[i], "%zu", &r2) != 1 ||
							++i == argc || sscanf(argv[i], "%zu", &r3) != 1 ||
							++i == argc || sscanf(argv[i], "%zu", &r4) != 1)
				  usage();		
				break;
			case 'n':
				if (++i == argc)
				  usage();
				dataset = argv[i];
				break;
			default: 
				usage();
				break;
		}
	}

	if ((r1==0) && (r2==0) && (r3==0) && (r4==0) && (r5==0))
	{
		printf ("Error: please specify dimensions.\n");
		printf("-1 <nx> : dimension for 1D data such as data[nx]\n");
		printf("-2 <nx> <ny> : dimensions for 2D data such as data[ny][nx]\n");
		printf("-3 <nx> <ny> <nz> : dimensions for 3D data such as data[nz][ny][nx] \n");
		printf("-4 <nx> <ny> <nz> <np>: dimensions for 4D data such as data[np][nz][ny][nx] \n");
		exit(0);
	}

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

	int status = 0;
	if (dataType == 0)
	{
		float *data = (float*)malloc(nbEle*sizeof(float));
		netcdfReader((void*)data, inPath, dataset, dataType);

		for (i = 0; i < 10; i++)
		  printf ("%f\n", ((float*)data)[i]);

		writeFloatData_inBytes(data, nbEle, outPath, &status);
	}
	else
	  if (dataType == 1)
	  {
		  double *data = (double*)malloc(nbEle*sizeof(double));
		  netcdfReader((void*)data, inPath, dataset, dataType);

		  for (i = 0; i < 10; i++)
			printf ("%lf\n", ((double*)data)[i]);
		  writeDoubleData_inBytes(data, nbEle, outPath, &status);
	  }

	return 0;
}

