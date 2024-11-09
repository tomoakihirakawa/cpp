#include <CL/cl.hpp>
#include <chrono>
#include <iostream>
#include <vector>

int main() {
   const int arrayLength = 1024 * 1024 * 100;  // Size for vectorized operations
   std::vector<float> inA(arrayLength, 1.0f);
   std::vector<float> inB(arrayLength, 2.0f);
   std::vector<float> result(arrayLength, 0.0f);

   // === CPU calculation with OpenMP parallelization ===
   auto cpu_start = std::chrono::high_resolution_clock::now();
   float cpu_sum = 0.0f;

#pragma omp parallel for reduction(+ : cpu_sum)
   for (int i = 0; i < arrayLength; ++i) {
      cpu_sum += inA[i] * inB[i] * inB[i];
   }

   auto cpu_end = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> cpu_duration = cpu_end - cpu_start;

   // === OpenCL initialization and computation ===
   std::vector<cl::Platform> platforms;
   cl::Platform::get(&platforms);
   if (platforms.empty()) {
      std::cerr << "No OpenCL platforms found." << std::endl;
      return -1;
   }
   cl::Platform platform = platforms.front();

   std::vector<cl::Device> devices;
   platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);
   if (devices.empty()) {
      std::cerr << "No GPU devices found on the platform." << std::endl;
      return -1;
   }
   cl::Device device = devices.front();

   cl::Context context(device);
   cl::Program::Sources sources;
   std::string kernel_code =
       "__kernel void add_vectors(__global const float4* inA, __global const float4* inB, __global float4* result) {"
       "    int id = get_global_id(0);"
       "    float4 a = inA[id];"
       "    float4 b = inB[id];"
       "    result[id] = a * b * b;"
       "}";
   sources.push_back({kernel_code.c_str(), kernel_code.length()});

   cl::Program program(context, sources);
   if (program.build({device}) != CL_SUCCESS) {
      std::cerr << "Error building: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device) << std::endl;
      return -1;
   }

   cl::Buffer bufferA(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float) * arrayLength, inA.data());
   cl::Buffer bufferB(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float) * arrayLength, inB.data());
   cl::Buffer bufferResult(context, CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float) * arrayLength, result.data());

   cl::CommandQueue queue(context, device);
   cl::Kernel kernel(program, "add_vectors");

   kernel.setArg(0, bufferA);
   kernel.setArg(1, bufferB);
   kernel.setArg(2, bufferResult);

   cl::NDRange globalWorkSize(arrayLength / 4);  // Adjusted for vectorized operations
   cl::NDRange localWorkSize(128);               // Example workgroup size, tune this value based on performance tests

   auto gpu_start = std::chrono::high_resolution_clock::now();

   queue.enqueueNDRangeKernel(kernel, cl::NullRange, globalWorkSize, localWorkSize);
   queue.finish();

   queue.enqueueReadBuffer(bufferResult, CL_TRUE, 0, sizeof(float) * arrayLength, result.data());

   auto gpu_end = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> gpu_duration = gpu_end - gpu_start;

   // Collect and display the results
   float gpu_sum = 0.0f;
   for (int i = 0; i < arrayLength; ++i) {
      gpu_sum += result[i];
   }

   std::cout << "CPU Sum: " << cpu_sum << ", Time: " << cpu_duration.count() << " seconds" << std::endl;
   std::cout << "GPU Sum: " << gpu_sum << ", Time: " << gpu_duration.count() << " seconds" << std::endl;

   return 0;
}