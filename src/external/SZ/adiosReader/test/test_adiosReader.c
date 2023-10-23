#include <stdio.h>
#include <stdlib.h>
#include "adiosReader.h"

void usage()
{
	printf("Options:\n");
	printf("* input data file:\n");
	printf("	-i <ADIOS bp data file> : ADIOS bp data file\n");
	printf("* dimensions: \n");
	printf("	-1 <nx> : dimension for 1D data such as data[nx]\n");
	printf("	-2 <nx> <ny> : dimensions for 2D data such as data[ny][nx]\n");
	printf("	-3 <nx> <ny> <nz> : dimensions for 3D data such as data[nz][ny][nx] \n");
	printf("	-4 <nx> <ny> <nz> <np>: dimensions for 4D data such as data[np][nz][ny][nx] \n");
	printf("	-4 <nx> <ny> <nz> <np> <nq>: dimensions for 5D data such as data[nq][np][nz][ny][nx] \n");
	printf("* examples: \n");
	printf("	test_adiosReader -f -i testdata/ADIOS2ADIOS1WriteADIOS1Read2D2x4Test.bp -2 4 2\n");
	exit(0);
}


int main(int argc, char **argv)
{	
	char* inPath = NULL;
	size_t i = 0;

	size_t r5 = 0;
	size_t r4 = 0;
	size_t r3 = 0;
	size_t r2 = 0; 
	size_t r1 = 0;

	int8_t *I8;
	int16_t *I16;
	int32_t *I32;
	int64_t *I64;
	uint8_t *U8;
	uint16_t *U16;
	uint32_t *U32;
	uint64_t *U64;
	float *R32;
	double *R64;

	if (argc == 1)
	  usage();

	for(i = 1;i < argc; i++)
	{
		if (argv[i][0] != '-' || argv[i][2])
		  usage();
		switch (argv[i][1])
		{
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
			case '5':
				if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1 ||
							++i == argc || sscanf(argv[i], "%zu", &r2) != 1 ||
							++i == argc || sscanf(argv[i], "%zu", &r3) != 1 ||
							++i == argc || sscanf(argv[i], "%zu", &r4) != 1 ||
							++i == argc || sscanf(argv[i], "%zu", &r5) != 1)
				  usage();
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
	  adiosReader_1D (inPath, r1, &I8, &I16, &I32, &I64, &U8, &U16, &U32, &U64, &R32, &R64);
	else if(r3==0)
	  adiosReader_2D (inPath, r1, r2,  &I8, &I16, &I32, &I64, &U8, &U16, &U32, &U64, &R32, &R64);
	else if(r4==0)
	  adiosReader_3D (inPath, r1, r2, r3, &I8, &I16, &I32, &I64, &U8, &U16, &U32, &U64, &R32, &R64);
	else if(r5==0)
	  adiosReader_4D (inPath, r1, r2, r3, r4, &I8, &I16, &I32, &I64, &U8, &U16, &U32, &U64, &R32, &R64);
	else
	  adiosReader_5D (inPath, r1, r2, r3, r4, r5, &I8, &I16, &I32, &I64, &U8, &U16, &U32, &U64, &R32, &R64);


	// Check
	for (i = 0; i < 5; i++)
	  printf ("I8 = %d, I16 = %d, I32 = %d, I64 = %lld, U8 = %u, U16 = %u, U32 = %u, U64 = %llu, R32 = %f, R64 = %lf\n", 
				  I8[i], I16[i], I32[i], I64[i], U8[i], U16[i], U32[i], U64[i], R32[i], R64[i]);


	free(I8);
	free(I16);
	free(I32);
	free(I64);
	free(U8);
	free(U16);
	free(U32);
	free(U64);
	free(R32);
	free(R64);



	return 0;
}
