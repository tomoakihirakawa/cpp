
#include "fundamental.hpp"
#include "OPENCL.hpp"
////////////////////////////////////
int main()
{

  const int n = 32*512;       // size of vectors
  const int k = 10;          // number of loop iterations
  
  OPENCL_kernel kernel("./kernel_code.h", "add_looped_1");

  std::vector<double> A(n,2.), B(n,3.), C(n,0);

  kernel.makeBuffer(0,A);
  kernel.makeBuffer(1,B);
  kernel.makeBuffer(2,C);  
  
  kernel.enqueueWriteBuffer(0,A);
  kernel.enqueueWriteBuffer(1,B);
  
  kernel.setArg(0, kernel.buffers[0]);
  kernel.setArg(1, kernel.buffers[1]);
  kernel.setArg(2, kernel.buffers[2]);
  kernel.setArg(3, n);
  kernel.setArg(4, k);
  
  kernel.enqueueNDRangeKernel(512,32);
  
  kernel.enqueueReadBuffer(2, C);
  std::cout << C << std::endl;
    
  return 0;
}
