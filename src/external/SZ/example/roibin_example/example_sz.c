#include <stdio.h>
#include <assert.h>
#include "sz.h"

size_t determine_n_peaks() {
  assert(0 && "remove this asssert after you have the code to determine the number of peaks");
  return 0;
}

uint8_t* load_calib(const char* path, size_t n_dims, size_t dims[]) {
  size_t total_dims = 1;
  for (size_t i = 0; i < n_dims; ++i) {
    total_dims *= dims[i];
  }
  uint8_t* calib_panel_ptr = malloc(total_dims * sizeof(uint8_t));

  //TODO fill in the calibration panel
  // In reality this is known from HDF5 and conversions have to be done for some formats
  // we expect each index to be a 8bit 0/1 mask which is 1 when a pixel should be IGNORED
  // it should be the same size and shape as the input
  (void)path;
  assert(0 && "remove this assert after you have the load_calib function");

  return calib_panel_ptr;
}

uint16_t* load_location_info(const char* path, size_t n_peaks) {
  uint16_t* location_data_ptr = malloc(n_peaks*sizeof(uint16_t));

  //TODO fill in the location_data panel
  // In reality this is known from HDF5 and conversions have to be done for some formats
  // we expect each index to be a 16bit unsigned integer that indicates one part of of a location of a peak

  //TODO ignore the path parameter, until it is filled in sanely
  assert(0 && "remove this assert after you have the load_location_info function");
  (void)path;

  return location_data_ptr;
}

float* load_input_data(const char* path, size_t n_dims, size_t dims[]) {
  size_t total_dims = 1;
  for (size_t i = 0; i < n_dims; ++i) {
    total_dims *= dims[i];
  }
  float* input_data = malloc(total_dims*sizeof(float));

  //TODO fill in the input data
  // In reality this is known from HDF5 and conversions have to be done for some formats
  // we expect this to be the actual float 32bit data with 3 or 4 dimensions
  // if it is 3d, we expect it to be panels, rows, cols from slowest to fastest incrementing index; we assume a single event for 3d
  // if it is 4d, we expect it to be events, panels, rows, cols from slowest to fastest incrementing index
  (void) path;
  assert(0 && "remove this assert after you have the load_input_data function");

  return input_data;
}


int main(int argc, char *argv[])
{
  SZ_Init(NULL);
  exafelSZ_params params;

  //TODO fill in code that sets the input dimensions
  const size_t n_dims = 4;
  size_t dims[4] = {0,0,0,0};
  assert(0 && "remove this asssert after you have the code to determine input data dimensions");

  const size_t n_peaks = determine_n_peaks();

  float* input = load_input_data("input.f32", n_dims, dims);
  params.binSize = 2;
  params.numPeaks = n_peaks;
  params.szDim = 3;
  params.peakSize = 3; //must be odd
  params.tolerance = 1e-6;
  params.calibPanel = load_calib("calib.u8", n_dims, dims);
  params.peaksCols = load_location_info("cols.u16", n_peaks);
  params.peaksRows = load_location_info("rows.u16", n_peaks);
  params.peaksSegs = load_location_info("segs.u16", n_peaks);
  assert(params.peaksRows != NULL && "remote this assert fill in the code to determine the number of peaks");
  assert(params.peaksCols != NULL && "remove this assert after you have the code to determine the number of peaks");
  assert(params.peaksSegs != NULL && "remove this assert after you have the code to determine the number of peaks");
  assert(params.calibPanel != NULL && "remove this assert after you have the code to determine the number of peaks");
  

  size_t output_size = 0;
  int status = 0;
  unsigned char* compressed = SZ_compress_customize("ExaFEL", &params, SZ_FLOAT, input, 0, dims[0], dims[1], dims[2], dims[3], &output_size, &status);
  if(status) {
    fprintf(stderr, "something went wrong during compression\n");
    exit(1);
  }
  float * output = SZ_decompress_customize("ExaFEL", &params, SZ_FLOAT, compressed, output_size, 0, dims[0], dims[1], dims[2], dims[3], &status);
  if(status) {
    fprintf(stderr, "something went wrong during compression\n");
    exit(1);
  }

  free(params.peaksSegs);
  free(params.peaksCols);
  free(params.peaksRows);
  free(params.calibPanel);
  free(input);
  free(compressed);
  free(output);

  return 0;
}
