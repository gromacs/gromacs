#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sz.h> 

void usage()
{
	printf("Usage: print_h5repack_args <options>\n");
	printf("Options:\n");
	printf("	-M <error bound mode> : 10 options as follows. \n");
	printf("		ABS (absolute error bound)\n");
	printf("		REL (value range based error bound, so a.k.a., VR_REL)\n");
	printf("		PSNR (peak signal-to-noise ratio)\n");
	printf("		PW_REL (point-wise relative error bound)\n");
	printf("	-A <absolute error bound>: specifying absolute error bound\n");
	printf("	-R <value_range based relative error bound>: specifying relative error bound\n");
	printf("	-P <point-wise relative error bound>: specifying point-wise relative error bound\n");
	printf("	-S <PSNR>: specifying PSNR\n");
	printf("* examples: \n");
	printf("	print_h5repack_args -M ABS -A 1E-3 (output: -f UD=32017,0,9,0,1062232653,3539053052,0,0,0,0,0,0)\n");
	printf("	print_h5repack_args -M REL -R 1E-4 (output: -f UD=32017,0,9,1,0,0,1058682594,3944497965,0,0,0,0)\n");
	exit(0);
}


int main(int argc, char* argv[])
{
	char* errBoundMode = NULL;
	char* absErrorBound = NULL;
	char* relErrorBound = NULL;
	char* pwrErrorBound = NULL;
	char* psnr_ = NULL;
	
	if(argc==1)
		usage();
	
	int i = 0;
	
	for(i=1;i<argc;i++)
	{
		if (argv[i][0] != '-' || argv[i][2])
			usage();
		switch (argv[i][1])
		{
		case 'h':
			usage();
			exit(0);
			break;
		case 'M':
			if (++i == argc)
				usage();
			errBoundMode = argv[i];
			break;
		case 'A':
			if (++i == argc)
				usage();
			absErrorBound = argv[i];
			break;
		case 'R':
			if (++i == argc)
				usage();
			relErrorBound = argv[i];
			break;
		case 'P':
			if (++i == argc)
				usage();
			pwrErrorBound = argv[i];
			break;
		case 'S': 
			if (++i == argc)
				usage();
			psnr_ = argv[i];
			break;
		default: 
			usage();
			break;
		}
	}
	unsigned char b[8];
	unsigned int cd_values[9];
	if(strcmp(errBoundMode, "ABS")==0)
	{
		cd_values[0] = ABS;
		doubleToBytes(b, atof(absErrorBound));
		cd_values[1] = bytesToInt32_bigEndian(b);
		cd_values[2] = bytesToInt32_bigEndian(b+4);
		cd_values[3] = 0;
		cd_values[4] = 0;
		cd_values[5] = 0;
		cd_values[6] = 0;
		cd_values[7] = 0;
		cd_values[8] = 0;
	}
	else if(strcmp(errBoundMode, "REL")==0)
	{
		cd_values[0] = REL;
		cd_values[1] = 0;
		cd_values[2] = 0;
		doubleToBytes(b, atof(relErrorBound));
		cd_values[3] = bytesToInt32_bigEndian(b);
		cd_values[4] = bytesToInt32_bigEndian(b+4);
		cd_values[5] = 0;
		cd_values[6] = 0;
		cd_values[7] = 0;
		cd_values[8] = 0;		
	}
	else if(strcmp(errBoundMode, "PW_REL")==0)
	{
		cd_values[0] = PW_REL;
		cd_values[1] = 0;
		cd_values[2] = 0;
		cd_values[3] = 0;
		cd_values[4] = 0;
		doubleToBytes(b, atof(pwrErrorBound));
		cd_values[5] = bytesToInt32_bigEndian(b);
		cd_values[6] = bytesToInt32_bigEndian(b+4);		
		cd_values[7] = 0;
		cd_values[8] = 0;
	}
	else if(strcmp(errBoundMode, "PSNR")==0)
	{
		cd_values[0] = PSNR;
		cd_values[1] = 0;
		cd_values[2] = 0;
		cd_values[3] = 0;
		cd_values[4] = 0;
		cd_values[5] = 0;
		cd_values[6] = 0;
		doubleToBytes(b, atof(psnr_));
		cd_values[7] = bytesToInt32_bigEndian(b);
		cd_values[8] = bytesToInt32_bigEndian(b+4);
	}
	else
	{
		printf("Error: wrong errBoundMode setting.\n");
		exit(0);
	}
	printf("-f UD=32017,0,9");

	for(i=0;i<9;i++)
		printf(",%u",cd_values[i]);
	printf("\n");
}
