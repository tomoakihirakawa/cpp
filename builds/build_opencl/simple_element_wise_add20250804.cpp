#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// OpenCL C++ラッパーを使うために、cl2.hppを準備してください
// https://github.com/KhronosGroup/OpenCL-CLHPP/blob/main/include/CL/cl2.hpp
// このファイルを cl.hpp として保存・インクルードするのが一般的です。
// #define CL_HPP_ENABLE_EXCEPTIONS
// #define CL_HPP_TARGET_OPENCL_VERSION 200  // macOSは1.2までですが、ヘッダーは新しいものでOK
#define CL_SILENCE_DEPRECATION
#include "opencl.hpp"

// ファイルを読み込むヘルパー関数
std::string read_file(const std::string& path) {
   std::ifstream file(path);
   if (!file.is_open()) {
      throw std::runtime_error("Kernel file not found: " + path);
   }
   return std::string(std::istreambuf_iterator<char>(file), (std::istreambuf_iterator<char>()));
}

int main() {

   const int num_elements = 10 * 1024 * 1024;  // 約1000万要素

   // --- CPUでの計算 ---
   std::cout << "1. Preparing data and running on CPU..." << std::endl;
   std::vector<float> a(num_elements);
   std::vector<float> b(num_elements);
   std::vector<float> c_cpu(num_elements);

   for (int i = 0; i < num_elements; ++i) {
      a[i] = static_cast<float>(i);
      b[i] = static_cast<float>(num_elements - i);
   }

   auto start_cpu = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
   for (int i = 0; i < num_elements; ++i) {
      float temp_c = 0.0f;
      for (int j = 0; j < 10000; j++) {  // ループを追加
         temp_c += a[i] + b[i];
      }
      c_cpu[i] = temp_c;
   }
   auto end_cpu = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double, std::milli> cpu_time = end_cpu - start_cpu;
   std::cout << "  CPU (OpenMP) time: " << cpu_time.count() << " ms" << std::endl;

   // --- GPUでの計算 ---
   std::cout << "\n2. Setting up OpenCL and running on GPU..." << std::endl;

   // プラットフォームとデバイスの選択
   std::vector<cl::Platform> platforms;
   cl::Platform::get(&platforms);
   cl::Platform platform = platforms.front();
   std::vector<cl::Device> devices;
   platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);
   cl::Device device = devices.front();

   std::cout << "  Using device: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;

   // コンテキストとコマンドキューの作成
   cl::Context context(device);
   cl::CommandQueue queue(context, device);

   // カーネルコードを読み込んでビルド
   std::string kernel_code = R"(
   __kernel void vector_add(__global const float *a, __global const float *b, __global float *c)
   {
    int i = get_global_id(0);
    
    // 計算結果を一時的に保持する変数
    float temp_c = 0.0f;

    // 同じ足し算を10000回繰り返す
    for (int j = 0; j < 10000; j++) {
        temp_c += a[i] + b[i];
    }

    c[i] = temp_c;
   })";

   cl::Program program(context, kernel_code);
   program.build("-cl-std=CL1.2");

   // GPUメモリ上のバッファを作成
   cl::Buffer buffer_a(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * num_elements, a.data());
   cl::Buffer buffer_b(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * num_elements, b.data());
   cl::Buffer buffer_c(context, CL_MEM_WRITE_ONLY, sizeof(float) * num_elements);

   // カーネルを作成し、引数を設定
   cl::Kernel kernel(program, "vector_add");
   kernel.setArg(0, buffer_a);
   kernel.setArg(1, buffer_b);
   kernel.setArg(2, buffer_c);

   // GPUでの計算時間を計測
   auto start_gpu = std::chrono::high_resolution_clock::now();

   // カーネルを実行 (グローバルサイズ = 要素数)
   // ローカルワークサイズを256に指定する例（コメントアウト）
   queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(num_elements), cl::NullRange);

   // 現在は自動選択を使用
   // queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(num_elements), cl::NullRange);
   queue.finish();  // 計算が終わるまで待つ

   auto end_gpu = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double, std::milli> gpu_time = end_gpu - start_gpu;
   std::cout << "  GPU time: " << gpu_time.count() << " ms" << std::endl;

   // 結果をGPUからCPUに読み戻す
   std::vector<float> c_gpu(num_elements);
   queue.enqueueReadBuffer(buffer_c, CL_TRUE, 0, sizeof(float) * num_elements, c_gpu.data());

   float total_c_gpu = 0.0f, total_c_cpu = 0.0f;
   for (int i = 0; i < num_elements; ++i) {
      total_c_gpu += c_gpu[i];
      total_c_cpu += c_cpu[i];
   }
   std::cout << "  Total of CPU results: " << total_c_cpu << std::endl;
   std::cout << "  Total of GPU results: " << total_c_gpu << std::endl;

   // --- 結果の検証 ---
   std::cout << "\n3. Verifying results..." << std::endl;
   bool success = true;
   for (int i = 0; i < num_elements; ++i) {
      if (std::abs(c_cpu[i] - c_gpu[i]) > 1e-5) {
         success = false;
         break;
      }
   }
   if (success) {
      std::cout << "  Success! CPU and GPU results match." << std::endl;
   } else {
      std::cout << "  Error! Results do not match." << std::endl;
   }

   std::cout << "\nSpeedup: " << cpu_time.count() / gpu_time.count() << "x" << std::endl;

   return 0;
}