#include <unistd.h>
#include <fcntl.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sz.h>
#include <math.h>

int main()
{
  //get a default set of parameters
  SZ_Init(NULL);
  //override the specifics we care about

  const size_t FULL = 500*500*100;
  const size_t STRIDE = 500*500*20;
  int hurricane_file = open("CLOUDf48.bin.f32", O_RDONLY);
  if(hurricane_file == -1) {
    perror("failed to open file: ");
    exit(1);
  }
  uint8_t* data = mmap(
    NULL,
    (500*500*100*sizeof(float)),
    PROT_READ,
    MAP_SHARED,
    hurricane_file,
    0
  );
  uint8_t* outdata = calloc(500*500*100, sizeof(float));
  if(data == MAP_FAILED) {
    perror("failed to map file: ");
    exit(1);
  }

   int stride_bytes = STRIDE*sizeof(float);
//#pragma omp parallel for default(none) shared(data,outdata,STRIDE)
//#pragma omp parallel for shared(data,outdata, stride_bytes)
#pragma omp parallel for shared(data,outdata)
  for (int i = 0; i < 5; ++i) {
    float* data_part = ((float*)data) + (i*STRIDE);
    size_t outSize = 0;

    //unsigned char* bytes = SZ_compress(SZ_FLOAT, data_part, &outsize, 0, 0, 20, 500, 500);
    
    unsigned char* bytes = SZ_compress_args(SZ_FLOAT, data_part, &outSize, REL, 0, 1E-4*(i+1), 0, 0, 0, 20, 500, 500);
    float* data_out = SZ_decompress(SZ_FLOAT, bytes, outSize, 0, 0, 20, 500, 500);

//================================================evaluation==================
    size_t nbEle = STRIDE;
    size_t j = 0;
    float Max = 0, Min = 0, diffMax = 0;
    float* ori_data = data_part;
    float* data2 = data_out;
    Max = ori_data[0];
    Min = ori_data[0];
    diffMax = fabs(data2[0] - ori_data[0]);
    double sum1 = 0, sum2 = 0;
    for (j = 0; j < nbEle; j++)
    {
        sum1 += ori_data[j];
	sum2 += data2[j];
    }
    double mean1 = sum1/nbEle;
    double mean2 = sum2/nbEle;

    double sum3 = 0, sum4 = 0;
    double sum = 0, prodSum = 0, relerr = 0;
   
    double maxpw_relerr = 0; 
    for (j = 0; j < nbEle; j++)
    {
        if (Max < ori_data[j]) Max = ori_data[j];
        if (Min > ori_data[j]) Min = ori_data[j];
        
        float err = fabs(data2[j] - ori_data[j]);
	if(ori_data[j]!=0)
	{
		if(fabs(ori_data[j])>1)
			relerr = err/ori_data[j];
		else
			relerr = err;
		if(maxpw_relerr<relerr)
			maxpw_relerr = relerr;
        }

	if (diffMax < err)
		diffMax = err;
        prodSum += (ori_data[j]-mean1)*(data2[j]-mean2);
        sum3 += (ori_data[j] - mean1)*(ori_data[j]-mean1);
        sum4 += (data2[j] - mean2)*(data2[j]-mean2);
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
     
    printf ("Min=%.20G, Max=%.20G, range=%.20G\n", Min, Max, range);
    printf ("Max absolute error = %.10f\n", diffMax);
    printf ("Max relative error = %f\n", diffMax/(Max-Min));
    printf ("Max pw relative error = %f\n", maxpw_relerr);
    printf ("PSNR = %f, NRMSE= %.20G\n", psnr,nrmse);
    printf ("acEff=%f\n", acEff);


    memcpy(outdata+(i*STRIDE)*sizeof(float), data_out, 20*500*500*sizeof(float));

    free(bytes);
    free(data_out);
  }

  float* in = (float*)data;
  float* out = (float*)outdata;
  double max = fabs(in[0]-out[0]);
  for (size_t i = 0; i < FULL; ++i) {
    max = fmax(fabs(in[i] - out[i]), max);
  }
  printf("max_error=%f\n", max);

  free(outdata);
  munmap(data, 100*500*500*sizeof(float));
  close(hurricane_file);
  SZ_Finalize();
  
  return 0;
}
