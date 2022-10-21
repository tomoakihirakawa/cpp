// test.cpp
#include <stdlib.h>
#include <iostream>
#include <fstream>
//#include <string>
#include <vector>
// #include <functional>
// #include <algorithm>
// #include <numeric>
//#include <cmath>
// #include <chrono>
// #include <sys/stat.h>
// #include <map>
// #include <sstream>
// //#define debug_fundamental

#define full_debug_fundamental

// #include "scw.h"

// #include "GNUPLOT3.h"
#include "fundamental.hpp"

#include "OPENCL.h"

// #include "interp2parm.h"
// #include "EarthScience.h"

////////////////////////////////////
#define pi 3.1415926535897932385
using namespace std;

int main()
{

  const int n = 32*512;       // size of vectors
  const int k = 10;          // number of loop iterations
  

  OPENCL_kernel kernel("./kernel_code", "add_looped_1");

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
  //cout << C << endl;
    
  return 0;
}
