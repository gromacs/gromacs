#include <cuda.h>
#include <cufft.h>
//#include <cufftXt.h>

void cuft_foo() {
  cufftHandle plan;
  cufftCreate(&plan);

  int rank = 3; // 3d
  int n[] = { 8, 8, 8 }; // dimensions
  int inembed[] = { 8, 8, 8 }; // input storage dimensions
  int istride = 4; // distance between two successive input elements 
  // in the least significant (i.e., innermost) dimension
  int idist = 8 * 8 * 8; // distance between the first element
  // of two consecutive signals in a batch of the input data

  int onembed[] = { 8, 8, 8 };
  int ostride = 8;
  int odist = 8 * 8 * 8;

  cufftType type = CUFFT_R2C;
  int batch = 4; // batch size
  cufftPlanMany(&plan, rank, n,
		inembed, istride, idist,
		onembed, ostride, odist,
		type, batch);
}
