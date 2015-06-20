#include <cuda.h>

int main(void)
{
  int *d;
  cudaMalloc(&d, sizeof(int));
}
