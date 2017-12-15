__global__ void kernel (void) {}

int
main()
{
  kernel<<<1,1>>>();
  return 0;
}
