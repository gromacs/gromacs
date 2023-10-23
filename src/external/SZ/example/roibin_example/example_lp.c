#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "libpressio.h"
#include "libpressio_ext/io/pressio_io.h"

size_t determine_n_peaks() {
  assert(0 && "remove this asssert after you have the code to determine the number of peaks");
  return 0;
}

struct pressio_data* load_calib(const char* path, size_t n_dims, size_t dims[]) {
  struct pressio_data* calib_panel = pressio_data_new_owning(pressio_uint8_dtype, n_dims, dims);
  uint8_t* calib_panel_ptr = pressio_data_ptr(calib_panel, NULL);

  //TODO fill in the calibration panel
  // In reality this is known from HDF5 and conversions have to be done for some formats
  // we expect each index to be a 8bit 0/1 mask which is 1 when a pixel should be IGNORED
  // it should be the same size and shape as the input
  (void)path;
  assert(0 && "remove this assert after you have the load_calib function");

  return calib_panel;
}

struct pressio_data* load_location_info(const char* path, size_t n_peaks) {
  struct pressio_data* location_data = pressio_data_new_owning(pressio_uint16_dtype, 1, &n_peaks);
  uint16_t* location_data_ptr = pressio_data_ptr(location_data, NULL);

  //TODO fill in the location_data panel
  // In reality this is known from HDF5 and conversions have to be done for some formats
  // we expect each index to be a 16bit unsigned integer that indicates one part of of a location of a peak

  //TODO ignore the path parameter, until it is filled in sanely
  assert(0 && "remove this assert after you have the load_location_info function");
  (void)path;

  return location_data;
}

struct pressio_data* load_input_data(const char* path, size_t n_dims, size_t dims[]) {
  struct pressio_data* input_data = pressio_data_new_owning(pressio_float_dtype, n_dims, dims);
  float* input_data_ptr = pressio_data_ptr(input_data, NULL);

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
  struct pressio* library =  pressio_instance();
  assert(library != NULL && "something is gravely wrong with your installation of libpressio");
  struct pressio_compressor* compressor = pressio_get_compressor(library, "sz");
  assert(library != NULL && "you build libpressio without sz");

  //TODO fill in code that sets the input dimensions
  const size_t n_dims = 4;
  size_t dims[4] = {0,0,0,0};
  assert(0 && "remove this asssert after you have the code to determine input data dimensions");

  struct pressio_data* input = load_input_data("input.f32", n_dims, dims);
  struct pressio_data* calib_panel = load_calib("calib.u8", n_dims, dims);
  assert(calib_panel != NULL && "remove this assert after you have the code to determine calib_panel");

  const size_t n_peaks = determine_n_peaks();
  struct pressio_data* peak_rows = load_location_info("rows.u16", n_peaks);
  struct pressio_data* peak_cols = load_location_info("cols.u16", n_peaks);
  struct pressio_data* peak_segs = load_location_info("segs.u16", n_peaks);
  assert(peak_rows != NULL && "remote this assert fill in the code to determine the number of peaks");
  assert(peak_cols != NULL && "remove this assert after you have the code to determine the number of peaks");
  assert(peak_segs != NULL && "remove this assert after you have the code to determine the number of peaks");



  struct pressio_options* options = pressio_options_new();
  assert(options != NULL && "failed to create options; something is gravely wrong");
  pressio_options_set_string(options, "sz:app", "ExaFEL");
  pressio_options_set_uinteger(options, "sz:exafel:bin_size", 2);
  pressio_options_set_uinteger(options, "sz:exafel:num_peaks", n_peaks);
  pressio_options_set_uinteger(options, "sz:exafel:sz_dim", 3);
  pressio_options_set_uinteger(options, "sz:exafel:peak_size", 3); //must be odd
  pressio_options_set_double(options, "sz:exafel:tolerance", 1e-6);
  pressio_options_set_data(options, "sz:exafel:calib_panel", calib_panel);
  pressio_options_set_data(options, "sz:exafel:peak_cols", peak_rows);
  pressio_options_set_data(options, "sz:exafel:peak_rows", peak_cols);
  pressio_options_set_data(options, "sz:exafel:peak_segs", peak_segs);

  if(pressio_compressor_check_options(compressor, options) != 0) {
    fprintf(stderr, "%s\n", pressio_compressor_error_msg(compressor));
    exit(pressio_compressor_error_code(compressor));
  }
  if(pressio_compressor_set_options(compressor, options) != 0) {
    fprintf(stderr, "%s\n", pressio_compressor_error_msg(compressor));
    exit(pressio_compressor_error_code(compressor));
  }
  pressio_data_free(calib_panel);
  pressio_data_free(peak_segs);
  pressio_data_free(peak_rows);
  pressio_data_free(peak_cols);
  pressio_options_free(options);

  struct pressio_data* compressed = pressio_data_new_empty(pressio_byte_dtype, 0, NULL);
  if(pressio_compressor_compress(compressor, input, compressed) != 0) {
    fprintf(stderr, "%s\n", pressio_compressor_error_msg(compressor));
    exit(pressio_compressor_error_code(compressor));
  }

  struct pressio_data* output = pressio_data_new_clone(input);
  if(pressio_compressor_decompress(compressor, compressed, output) != 0) {
    fprintf(stderr, "%s\n", pressio_compressor_error_msg(compressor));
    exit(pressio_compressor_error_code(compressor));
  }


  pressio_data_free(input);
  pressio_data_free(compressed);
  pressio_data_free(output);
  pressio_compressor_release(compressor);
  pressio_release(library);
  
  return 0;
}
