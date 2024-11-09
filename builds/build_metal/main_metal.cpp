// clang++ -std=c++17 -Imetal-cpp -framework Metal -framework Foundation -o metal_example main_metal.cpp

#define NS_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>
#include <QuartzCore/QuartzCore.hpp>
#include <iostream>

int main() {
   // Metalデバイスの取得
   MTL::Device* device = MTL::CreateSystemDefaultDevice();

   // if (!device) {
   //    std::cerr << "Metal is not supported on this device." << std::endl;
   //    return -1;
   // } else {
   //    std::cout << "Metal is supported on this device." << std::endl;
   //    std::cout << "GPU Name: " << device->name()->cString(NS::UTF8StringEncoding) << std::endl;
   //    std::cout << "Max Threads per Threadgroup: " << device->maxThreadsPerThreadgroup().width << std::endl;
   // }

   // コマンドキューの作成
   MTL::CommandQueue* commandQueue = device->newCommandQueue();

   // if (!commandQueue) {
   //    std::cerr << "Failed to create command queue." << std::endl;
   //    return -1;
   // }

   // 簡単な計算用のバッファ作成
   const int dataSize = 1024;
   float input[dataSize];
   for (int i = 0; i < dataSize; ++i) {
      input[i] = static_cast<float>(i);
   }

   MTL::Buffer* inputBuffer = device->newBuffer(input, sizeof(input), MTL::ResourceStorageModeShared);
   MTL::Buffer* outputBuffer = device->newBuffer(sizeof(input), MTL::ResourceStorageModeShared);

   // シェーダーのロード
   NS::Error* error = nullptr;
   NS::String* filePath = NS::String::string("first_test.metallib", NS::UTF8StringEncoding);
   MTL::Library* library = device->newLibrary(filePath, &error);

   // if (!library) {
   //    std::cerr << "Failed to load library: " << error->localizedDescription()->cString(NS::UTF8StringEncoding) << std::endl;
   //    return -1;
   // }

   MTL::Function* function = library->newFunction(NS::String::string("add_arrays", NS::UTF8StringEncoding));
   MTL::ComputePipelineState* pipelineState = device->newComputePipelineState(function, &error);

   // if (!pipelineState) {
   //    std::cerr << "Failed to create pipeline state: " << error->localizedDescription()->cString(NS::UTF8StringEncoding) << std::endl;
   //    return -1;
   // }

   // コマンドバッファとコンピュートエンコーダの作成
   MTL::CommandBuffer* commandBuffer = commandQueue->commandBuffer();
   MTL::ComputeCommandEncoder* encoder = commandBuffer->computeCommandEncoder();

   // if (!encoder) {
   //    std::cerr << "Failed to create compute command encoder." << std::endl;
   //    return -1;
   // }

   // バッファの設定
   encoder->setBuffer(inputBuffer, 0, 0);
   encoder->setBuffer(outputBuffer, 0, 1);
   encoder->setComputePipelineState(pipelineState);

   // スレッドグループの設定
   MTL::Size gridSize = MTL::Size(dataSize, 1, 1);
   MTL::Size threadGroupSize = MTL::Size(16, 1, 1);
   encoder->dispatchThreads(gridSize, threadGroupSize);

   // エンコードと実行
   encoder->endEncoding();
   commandBuffer->commit();
   commandBuffer->waitUntilCompleted();

   // 結果の表示
   float* result = static_cast<float*>(outputBuffer->contents());
   std::cout << "First element in buffer: " << result[0] << std::endl;

   return 0;
}