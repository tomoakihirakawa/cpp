#include <iostream>
#include <ctime>
#include <vector>
#include "cl.hpp"
#include "fundamental.hpp"

#define NUM_GLOBAL_WITEMS 1024

#if defined(cl_khr_fp64) // Khronos extension available?
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
// OS X is supported double precision and ther is no cl_khr_fp64, therefore no need to enable it
#define DOUBLE_SUPPORT_AVAILABLE
#elif defined(cl_amd_fp64) // AMD extension available?
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#define DOUBLE_SUPPORT_AVAILABLE
#endif

void compareResults(double CPUtime, double GPUtime, int trial)
{
  double time_ratio = (CPUtime / GPUtime);
  std::cout << "VERSION " << trial << " -----------" << std::endl;
  std::cout << "CPU time: " << CPUtime << std::endl;
  std::cout << "GPU time: " << GPUtime << std::endl;
  std::cout << "GPU is ";
  if (time_ratio > 1)
    std::cout << time_ratio << " times faster!" << std::endl;
  else
    std::cout << (1 / time_ratio) << " times slower :(" << std::endl;
}

double timeAddVectorsCPU(const int n, const int k)
{
  // adds two vectors of size n, k times, returns total duration
  std::clock_t start;
  double duration;
  double A[n], B[n], C[n];
  for (auto i = 0; i < n; i++)
  {
    A[i] = i;
    B[i] = n - i;
    C[i] = 0;
  }

  start = std::clock();
#ifdef _OPENMP
  std::cout << "Parallelized by OpenMP" << std::endl;
#pragma omp parallel for
#endif
  for (auto ii = 0; ii < 10; ii++)
    for (auto j = 0; j < n; j++)
      for (auto i = 0; i < k; i++)
        C[j] = A[j] * B[j];

  duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
  return duration;
}

void warmup(cl::Context &context,
            cl::CommandQueue &queue,
            cl::Kernel &add,
            const double *A,
            const double *B,
            const int n)
{
  double C[n];
  // allocate space
  cl::Buffer buffer_A(context, CL_MEM_READ_WRITE, sizeof(int) * n);
  cl::Buffer buffer_B(context, CL_MEM_READ_WRITE, sizeof(int) * n);
  cl::Buffer buffer_C(context, CL_MEM_READ_WRITE, sizeof(int) * n);

  // push write commands to queue
  queue.enqueueWriteBuffer(buffer_A, CL_TRUE, 0, sizeof(int) * n, &A);
  queue.enqueueWriteBuffer(buffer_B, CL_TRUE, 0, sizeof(int) * n, &B);

  // RUN ZE KERNEL
  add.setArg(1, buffer_B);
  add.setArg(0, buffer_A);
  add.setArg(2, buffer_C);

  for (auto i = 0; i < 5; i++)
  {
    std::cout << i << std::endl;
    queue.enqueueNDRangeKernel(add, cl::NullRange, cl::NDRange(NUM_GLOBAL_WITEMS), cl::NDRange(32));
  }
  queue.enqueueReadBuffer(buffer_C, CL_TRUE, 0, sizeof(int) * n, C);
  queue.finish();
}
//////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{

  bool verbose;
  if (argc == 1 || std::strcmp(argv[1], "0") == 0)
    verbose = true;
  else
    verbose = false;

  const int n = 10 * 32 * 512; // size of vectors
  const int k = 1000;          // number of loop iterations

  // get all platforms (drivers), e.g. NVIDIA
  std::vector<cl::Platform> all_platforms;
  cl::Platform::get(&all_platforms);

  //////////////////////
  for (const auto &pfm : all_platforms)
  {
    std::cout << std::setw(40) << "CL_PLATFORM_EXTENSIONS: " << std::setw(30) << pfm.getInfo<CL_PLATFORM_EXTENSIONS>() << std::endl;
    std::cout << std::setw(40) << "CL_PLATFORM_NAME: " << std::setw(30) << pfm.getInfo<CL_PLATFORM_NAME>() << std::endl;
    std::cout << std::setw(40) << "CL_PLATFORM_PROFILE: " << std::setw(30) << pfm.getInfo<CL_PLATFORM_PROFILE>() << std::endl;
    std::cout << std::setw(40) << "CL_PLATFORM_VENDOR: " << std::setw(30) << pfm.getInfo<CL_PLATFORM_VENDOR>() << std::endl;
    std::cout << std::setw(40) << "CL_PLATFORM_VERSIONN: " << std::setw(30) << pfm.getInfo<CL_PLATFORM_VERSION>() << std::endl;
  }

  cl::Platform default_platform{all_platforms[0]};
  // get default device of the default platform
  std::vector<cl::Device> all_devices;
  default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);

  for (const auto &dev : all_devices)
  {
    std::cout << std::setw(40) << "CL_DEVICE_NAME: " << std::setw(30) << dev.getInfo<CL_DEVICE_NAME>() << std::endl;
    std::cout << std::setw(40) << "CL_DEVICE_MAX_CLOCK_FREQUENCY: " << std::setw(30) << dev.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << std::endl;
    std::cout << std::setw(40) << "CL_DEVICE_MAX_WORK_ITEM_SIZES: " << std::setw(30) << dev.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>() << std::endl;
    std::cout << std::setw(40) << "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS: " << std::setw(30) << dev.getInfo<CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>() << std::endl;
    std::cout << std::setw(40) << "CL_DEVICE_MAX_WRITE_IMAGE_ARGS: " << std::setw(30) << dev.getInfo<CL_DEVICE_MAX_WRITE_IMAGE_ARGS>() << std::endl;
  }

  if (all_devices.size() == 0)
  {
    std::cout << std::setw(30) << " No devices found. Check OpenCL installation!\n";
    exit(1);
  }

  // use device[1] because that's a GPU; device[0] is the CPU
  cl::Device default_device = all_devices[1];

  std::cout << default_device.getInfo<CL_DEVICE_NAME>() << std::endl;
  ////////////////////////////////

  // std::cout<< "Using device: "<<default_device.getInfo<CL_DEVICE_NAME>()<<"\n";
  cl::Context context({default_device});
  cl::Program::Sources sources;

  // calculates for each element; C = A + B
  std::string kernel_code = "./kernel_code";
  sources.push_back({kernel_code.c_str(), kernel_code.length()});

  cl::Program program(context, sources);
  if (program.build({default_device}) != CL_SUCCESS)
  {
    std::cout << "Error building: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device) << std::endl;
    exit(1);
  }

  // run the CPU code
  double CPUtime = timeAddVectorsCPU(n, k);

  // set up kernels and vectors for GPU code
  cl::CommandQueue queue(context, default_device);
  cl::Kernel add = cl::Kernel(program, "add");
  cl::Kernel add_looped_1 = cl::Kernel(program, "add_looped_1");
  cl::Kernel add_looped_2 = cl::Kernel(program, "add_looped_2");
  cl::Kernel add_single = cl::Kernel(program, "add_single");

  // construct vectors
  double A[n], B[n], C[n];
  for (auto i = 0; i < n; i++)
  {
    A[i] = i;
    B[i] = n - i - 1;
  }

  std::cout << "// attempt at warm-up..." << std::endl;
  warmup(context, queue, add, A, B, n);
  queue.finish();

  std::cout << "// VERSION 1 ==========================================" << std::endl;
  // allocate space
  cl::Buffer buffer_A(context, CL_MEM_READ_WRITE, sizeof(double) * n);
  cl::Buffer buffer_B(context, CL_MEM_READ_WRITE, sizeof(double) * n);
  cl::Buffer buffer_C(context, CL_MEM_READ_WRITE, sizeof(double) * n);

  // push write commands to queue
  queue.enqueueWriteBuffer(buffer_A, CL_TRUE, 0, sizeof(double) * n, A);
  queue.enqueueWriteBuffer(buffer_B, CL_TRUE, 0, sizeof(double) * n, B);

  // RUN ZE KERNEL
  add_looped_1.setArg(0, buffer_A);
  add_looped_1.setArg(1, buffer_B);
  add_looped_1.setArg(2, buffer_C);
  add_looped_1.setArg(3, n);
  add_looped_1.setArg(4, k);

  std::clock_t start_time = std::clock();

#ifdef _OPENMP
  std::cout << "Parallelized by OpenMP" << std::endl;
#pragma omp parallel for
#endif
  for (auto i = 0; i < 10; i++)
  {
    std::cout << i << std::endl;
    queue.enqueueNDRangeKernel(add_looped_1, cl::NullRange,    // kernel, offset
                               cl::NDRange(NUM_GLOBAL_WITEMS), // global number of work items
                               cl::NDRange(32));               // local number (per group)
  }

  double GPUtime1 = (std::clock() - start_time) / (double)CLOCKS_PER_SEC;

  // read result from GPU to here; including for the sake of timing
  queue.enqueueReadBuffer(buffer_C, CL_TRUE, 0, sizeof(double) * n, C);
  queue.finish();

  std::cout << "// VERSION 2 ==========================================" << std::endl;

  cl::Buffer buffer_A2(context, CL_MEM_READ_WRITE, sizeof(double) * n);
  cl::Buffer buffer_B2(context, CL_MEM_READ_WRITE, sizeof(double) * n);
  cl::Buffer buffer_C2(context, CL_MEM_READ_WRITE, sizeof(double) * n);
  queue.enqueueWriteBuffer(buffer_A2, CL_TRUE, 0, sizeof(double) * n, A);
  queue.enqueueWriteBuffer(buffer_B2, CL_TRUE, 0, sizeof(double) * n, B);

  add_looped_2.setArg(0, buffer_A2);
  add_looped_2.setArg(1, buffer_B2);
  add_looped_2.setArg(2, buffer_C2);
  add_looped_2.setArg(3, n);
  add_looped_2.setArg(4, k);

  start_time = std::clock();
  queue.enqueueNDRangeKernel(add_looped_2, cl::NullRange, cl::NDRange(NUM_GLOBAL_WITEMS), cl::NDRange(32));
  double GPUtime2 = (std::clock() - start_time) / (double)CLOCKS_PER_SEC;

  queue.enqueueReadBuffer(buffer_C2, CL_TRUE, 0, sizeof(int) * n, C);
  queue.finish();
  // let's compare!
  const int NUM_VERSIONS = 2;
  double GPUtimes[NUM_VERSIONS] = {GPUtime1, GPUtime2};
  if (verbose)
  {
    for (int i = 0; i < NUM_VERSIONS; i++)
      compareResults(CPUtime, GPUtimes[i], i + 1);
  }
  else
  {
    std::cout << CPUtime << ",";
    for (int i = 0; i < NUM_VERSIONS - 1; i++)
      std::cout << GPUtimes[i] << ",";
    std::cout << GPUtimes[NUM_VERSIONS - 1] << std::endl;
  }
  return 0;
}
