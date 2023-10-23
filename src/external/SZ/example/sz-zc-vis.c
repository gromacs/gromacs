#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sz.h"
#include "rw.h"
#include "zc.h"

void usage()
{
	printf("Usage: sz <options>\n");
	printf("Options:\n");
	printf("* Z-checker parameters:\n");
	printf("	-v <variable name>: variable name\n");
	printf("	-y <solution name>: solution name\n");
	printf("	-C <zc configuration file> : configuraiton file zc.config\n");	
	printf("* operation type:\n");
	printf("	-h: print the help information\n");
	printf("* data type:\n");
	printf("	-f: single precision (float type)\n");
	printf("	-d: double precision (double type)\n");
	printf("* configuration file: \n");
	printf("	-c <sz configuration file> : configuration file sz.config\n");
	printf("* error control: (the error control parameters here will overwrite the setting in sz.config)\n");
	printf("	-M <error bound mode> : 10 options as follows. \n");
	printf("		ABS (absolute error bound)\n");
	printf("		REL (value range based error bound)\n");
	printf("		ABS_AND_REL (using min{ABS, REL})\n");
	printf("		ABS_OR_REL (using max{ABS, REL})\n");
	printf("		PSNR (peak signal-to-noise ratio)\n");
	printf("		PW_REL (point-wise relative error bound)\n");
	printf("	-A <absolute error bound>: specifying absolute error bound\n");
	printf("	-R <value_range based relative error bound>: specifying relative error bound\n");
	printf("	-P <point-wise relative error bound>: specifying point-wise relative error bound\n");
	printf("	-S <PSNR>: specifying PSNR\n");
	printf("* input data file:\n");
	printf("	-i <original data file> : original data file\n");
	printf("	-s <compressed data file> : compressed data file in decompression\n");
	printf("* output type of decompressed file: \n");
	printf("	-b (by default) : decompressed file stored in binary format\n");
	printf("	-t : decompreadded file stored in text format\n");
	printf("	-T : pre-processing with Tucker Tensor Decomposition\n");
	printf("* dimensions: \n");
	printf("	-1 <nx> : dimension for 1D data such as data[nx]\n");
	printf("	-2 <nx> <ny> : dimensions for 2D data such as data[ny][nx]\n");
	printf("	-3 <nx> <ny> <nz> : dimensions for 3D data such as data[nz][ny][nx] \n");
	printf("	-4 <nx> <ny> <nz> <np>: dimensions for 4D data such as data[np][nz][ny][nx] \n");
	printf("* print compression results: \n");
	printf("	-a : print compression results such as distortions\n");
	printf("* examples: \n");
	printf("	sz_zc -f -i ~/Data/Hurricane-ISA/CLOUDf48.bin.dat -3 500 500 100 -v CLOUDf48 -y \"SZ(1E-3):CLOUDf48\" -C zc.config -M REL -R 1E-3\n");
	exit(0);
}


int main(int argc, char* argv[])
{
	int binaryOutput = 1;
	int printCmpResults = 0;
	int isCompression = -1000; //1 : compression ; 0: decompression
	int printMeta = 0;
	int dataType = 0; //0: single precision ; 1: double precision
	int tucker = 0; //0: without tucker tensor decomposition preprocessing; 1: with tucker tensor decomposition
	char* inPath = NULL;
	char* cmpPath = NULL;
	char* conPath = NULL;
	char* zcConPath = NULL;
	char* decPath = NULL;
	
	char* varName = NULL;
	char* solName = NULL;
	
	char* errBoundMode = NULL;
	char* absErrorBound = NULL;
	char* relErrorBound = NULL;
	char* pwrErrorBound = NULL;
	char* psnr_ = NULL;
	
	size_t r5 = 0;
	size_t r4 = 0;
	size_t r3 = 0;
	size_t r2 = 0; 
	size_t r1 = 0;
	
	size_t i = 0;
	int status;
	size_t nbEle;
	if(argc==1)
		usage();
	
	for(i=1;i<argc;i++)
	{
		if (argv[i][0] != '-' || argv[i][2])
			usage();
		switch (argv[i][1])
		{
		case 'h':
			usage();
			exit(0);
		case 'b': 
			binaryOutput = 1;
			break;
		case 't': 
			binaryOutput = 0;
			break;
		case 'a':
			printCmpResults = 1;
			break;
		case 'p':
			printMeta = 1; //print metadata
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
		case 's':
			if (++i == argc)
				usage();
			cmpPath = argv[i];
			break;
		case 'c':
			if (++i == argc)
				usage();
			conPath = argv[i];
			break;
		case 'C':
			if (++i == argc)
				usage();
			zcConPath = argv[i];
			break;
		case 'v': 
			if (++i == argc)
				usage();
			varName = argv[i];
			break;	
		case 'y': 
			if (++i == argc)
				usage();
			solName = argv[i];
			break;						
		case 'T':
			tucker = 1;
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

	if((inPath==NULL) & (cmpPath == NULL))
	{
		printf("Error: you need to specify either a raw binary data file or a compressed data file as input\n");
		usage();
		exit(0);
	}

	if(printMeta == 0)
	{
		if ((r1==0) && (r2==0) && (r3==0) && (r4==0) && (r5==0))
		{
			printf ("Error: please specify dimensions.\n");
			printf("-1 <nx> : dimension for 1D data such as data[nx]\n");
			printf("-2 <nx> <ny> : dimensions for 2D data such as data[ny][nx]\n");
			printf("-3 <nx> <ny> <nz> : dimensions for 3D data such as data[nz][ny][nx] \n");
			printf("-4 <nx> <ny> <nz> <np>: dimensions for 4D data such as data[np][nz][ny][nx] \n");
			exit(0);
		}		
	}
	else
	{
		if(cmpPath == NULL && isCompression != 1) //if no compression file is provided and this is not a compression operation
		{
			printf("Error: -p can only be used when providing a compressed data file or in the compression step\n");
			printf("Solution: use -s to specify a compressed data file or use -c and -i to generate a compressed file\n");
			usage();
			exit(0);
		}
	}
	
	//Initialization (only for compression because decompression doesn't need the initialization)
	
	SZ_Init(conPath);	
	ZC_Init(zcConPath);

	plot_dec_data = 1;

	int errorBoundMode = 0;
	if(strcmp(errBoundMode, "ABS")==0)
		errorBoundMode = ABS;
	else if(strcmp(errBoundMode, "REL")==0)
		errorBoundMode = REL;
	else if(strcmp(errBoundMode, "ABS_AND_REL")==0)
		errorBoundMode = ABS_AND_REL;
	else if(strcmp(errBoundMode, "ABS_OR_REL")==0)
		errorBoundMode = ABS_OR_REL;
	else if(strcmp(errBoundMode, "PSNR")==0)
		errorBoundMode = PSNR;
	else if(strcmp(errBoundMode, "PW_REL")==0)
		errorBoundMode = PW_REL;
	else
	{
		printf("Error: wrong error bound mode setting by using the option '-M'\n");
		usage();
		exit(0);
	}
	confparams_cpr->errorBoundMode = errorBoundMode;
	
	char outputFilePath[256];	
	unsigned char *bytes = NULL; //the binary data read from "compressed data file"

	if(absErrorBound != NULL)
		confparams_cpr->absErrBound = atof(absErrorBound);
	
	if(relErrorBound != NULL)
		confparams_cpr->relBoundRatio = atof(relErrorBound);

	if(pwrErrorBound != NULL)
		confparams_cpr->pw_relBoundRatio = atof(pwrErrorBound);

	if(psnr_ != NULL)
		confparams_cpr->psnr = atof(psnr_);
	
	size_t outSize;	
	if(dataType == 0) //single precision
	{
		float *ori_data = readFloatData(inPath, &nbEle, &status);
		if(status!=SZ_SCES)
		{
			printf("Error: cannot read the input file: %s\n", inPath);
			exit(0);
		}

		//compression
		ZC_DataProperty* dataProperty = ZC_startCmpr(varName, ZC_FLOAT, ori_data, r5, r4, r3, r2, r1);
		bytes = SZ_compress(SZ_FLOAT, ori_data, &outSize, r5, r4, r3, r2, r1);
		ZC_CompareData* compareResult = ZC_endCmpr(dataProperty, solName, outSize);
			
		//decompression
		ZC_startDec();
		float *dec_data = SZ_decompress(SZ_FLOAT, bytes, outSize, r5, r4, r3, r2, r1);
		ZC_endDec(compareResult, dec_data);	
		
		freeDataProperty(dataProperty);
		freeCompareResult(compareResult);
		free(bytes);
		
		if(inPath==NULL)
		{
			printf("Error: Since you add -a option (analysis), please specify the original data path by -i <path>.\n");
			exit(0);
		}
		//compute the distortion / compression errors...
		if(status!=SZ_SCES)
		{
			printf("Error: %s cannot be read!\n", inPath);
			exit(0);
		}

		size_t i = 0;
		float Max = 0, Min = 0, diffMax = 0;
		Max = ori_data[0];
		Min = ori_data[0];
		diffMax = fabs(dec_data[0] - ori_data[0]);
		double sum1 = 0, sum2 = 0;
		for (i = 0; i < nbEle; i++)
		{
			sum1 += ori_data[i];
			sum2 += dec_data[i];
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
			
			float err = fabs(dec_data[i] - ori_data[i]);
			if(ori_data[i]!=0)
			{
				relerr = err/fabs(ori_data[i]);
				//if(relerr>0.1)
				//{printf("i=%zu, %.30G\n", i, relerr); break;}
				if(maxpw_relerr<relerr)
					maxpw_relerr = relerr;
			}

			if (diffMax < err)
				diffMax = err;
			prodSum += (ori_data[i]-mean1)*(dec_data[i]-mean2);
			sum3 += (ori_data[i] - mean1)*(ori_data[i]-mean1);
			sum4 += (dec_data[i] - mean2)*(dec_data[i]-mean2);
			sum += err*err;	
		}
		double std1 = sqrt(sum3/nbEle);
		double std2 = sqrt(sum4/nbEle);
		double ee = prodSum/nbEle;
		double acEff = ee/std1/std2;

		double mse = sum/nbEle;
		double range = Max - Min;
		double psnr = 20*log10(range)-10*log10(mse);
		double nrmse = sqrt(mse)/range;
		double compressionRatio = 1.0*nbEle*sizeof(float)/outSize;

		printf ("Min=%.20G, Max=%.20G, range=%.20G\n", Min, Max, range);
		printf ("Max absolute error = %.10f\n", diffMax);
		printf ("Max relative error = %f\n", diffMax/(Max-Min));
		printf ("Max pw relative error = %f\n", maxpw_relerr);
		printf ("PSNR = %f, NRMSE= %.20G\n", psnr,nrmse);
		printf ("acEff=%f\n", acEff);	
		printf ("compressionRatio=%f\n", compressionRatio);
		
		free(ori_data);
		free(dec_data);
	}
	else //dataType == 1: double precision
	{
		double *ori_data = readDoubleData(inPath, &nbEle, &status);	
		if(status!=SZ_SCES)
		{
			printf("Error: cannot read the input file: %s\n", inPath);
			exit(0);
		}
		
		ZC_DataProperty* dataProperty = ZC_startCmpr(varName, ZC_DOUBLE, ori_data, r5, r4, r3, r2, r1);
		bytes = SZ_compress(SZ_DOUBLE, ori_data, &outSize, r5, r4, r3, r2, r1);
		ZC_CompareData* compareResult = ZC_endCmpr(dataProperty, solName, outSize);
	
		//decompression
		ZC_startDec();
		double *dec_data = SZ_decompress(SZ_DOUBLE, bytes, outSize, r5, r4, r3, r2, r1);
		ZC_endDec(compareResult, dec_data);	
		
		freeDataProperty(dataProperty);
		freeCompareResult(compareResult);
		free(bytes);
		
		if(inPath==NULL)
		{
			printf("Error: Since you add -a option (analysis), please specify the original data path by -i <path>.\n");
			exit(0);
		}

		if(status!=SZ_SCES)
		{
			printf("Error: %s cannot be read!\n", inPath);
			exit(0);
		}

		size_t i = 0;
		double Max = 0, Min = 0, diffMax = 0;
		Max = ori_data[0];
		Min = ori_data[0];
		diffMax = dec_data[0]>ori_data[0]?dec_data[0]-ori_data[0]:ori_data[0]-dec_data[0];

		//diffMax = fabs(data[0] - ori_data[0]);
		double sum1 = 0, sum2 = 0;

		for (i = 0; i < nbEle; i++)
		{
			sum1 += ori_data[i];
			sum2 += dec_data[i];
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

			float err = fabs(dec_data[i] - ori_data[i]);
			if(ori_data[i]!=0)
			{
				relerr = err/fabs(ori_data[i]);
				if(maxpw_relerr<relerr)
				  maxpw_relerr = relerr;
			}

			if (diffMax < err)
			  diffMax = err;
			prodSum += (ori_data[i]-mean1)*(dec_data[i]-mean2);
			sum3 += (ori_data[i] - mean1)*(ori_data[i]-mean1);
			sum4 += (dec_data[i] - mean2)*(dec_data[i]-mean2);
			sum += err*err;	
		}
		double std1 = sqrt(sum3/nbEle);
		double std2 = sqrt(sum4/nbEle);
		double ee = prodSum/nbEle;
		double acEff = ee/std1/std2;

		double mse = sum/nbEle;
		double range = Max - Min;
		double psnr = 20*log10(range)-10*log10(mse);
		double nrmse = sqrt(mse)/range;

		double compressionRatio = 1.0*nbEle*sizeof(double)/outSize;

		printf ("Min = %.20G, Max = %.20G, range = %.20G\n", Min, Max, range);
		printf ("Max absolute error = %.10f\n", diffMax);
		printf ("Max relative error = %f\n", diffMax/(Max-Min));
		printf ("Max pw relative error = %f\n", maxpw_relerr);
		printf ("PSNR = %f, NRMSE = %.20G\n", psnr,nrmse);
		printf ("acEff = %f\n", acEff);
		printf ("compressionRatio = %f\n", compressionRatio);
		
		free(ori_data);
		free(dec_data);
	}
	
	SZ_Finalize();
	ZC_Finalize();
}
