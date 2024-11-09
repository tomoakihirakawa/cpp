// clang++ -std=c++17 -I/path/to/metal-cpp -framework Metal -framework Foundation -o metal_example main_metal.cpp

#define NS_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>
#include <iostream>

int main() {
    MTL::Device* device = MTL::CreateSystemDefaultDevice();
    
    if (!device) {
        std::cerr << "Metal is not supported on this device." << std::endl;
    } else {
        std::cout << "Metal is supported on this device." << std::endl;
    }

    return 0;
}


// #include <Metal/Metal.h>
// #include <Foundation/Foundation.h>
// #include <iostream>
// #include <vector>
// #include <chrono>

// int main() {
//     const int arrayLength = 1024 * 1024 * 100;
//     std::vector<float> inA(arrayLength, 1.0f);
//     std::vector<float> inB(arrayLength, 2.0f);
//     std::vector<float> result(arrayLength, 0.0f);

//     // === CPU calculation with OpenMP parallelization ===
//     auto cpu_start = std::chrono::high_resolution_clock::now();
//     float cpu_sum = 0.0f;

//     #pragma omp parallel for reduction(+ : cpu_sum)
//     for (int i = 0; i < arrayLength; ++i) {
//         cpu_sum += inA[i] * inB[i] * inB[i];
//     }

//     auto cpu_end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> cpu_duration = cpu_end - cpu_start;

//     // === Metal setup and computation ===
//     @autoreleasepool {
//         id<MTLDevice> device = MTLCreateSystemDefaultDevice();
//         if (!device) {
//             std::cerr << "Metal is not supported on this device." << std::endl;
//             return -1;
//         }

//         id<MTLCommandQueue> commandQueue = [device newCommandQueue];
//         NSError *error = nil;

//         NSString *shaderPath = [NSString stringWithUTF8String:"./add_vectors.metal"];
//         NSString *shaderSource = [NSString stringWithContentsOfFile:shaderPath encoding:NSUTF8StringEncoding error:&error];
//         if (!shaderSource) {
//             std::cerr << "Failed to load Metal shader source: " << error.localizedDescription.UTF8String << std::endl;
//             return -1;
//         }

//         id<MTLLibrary> library = [device newLibraryWithSource:shaderSource options:nil error:&error];
//         if (!library) {
//             std::cerr << "Failed to load Metal library: " << error.localizedDescription.UTF8String << std::endl;
//             return -1;
//         }

//         id<MTLFunction> function = [library newFunctionWithName:@"add_vectors"];
//         if (!function) {
//             std::cerr << "Failed to find Metal function." << std::endl;
//             return -1;
//         }

//         id<MTLComputePipelineState> pipelineState = [device newComputePipelineStateWithFunction:function error:&error];
//         if (!pipelineState) {
//             std::cerr << "Failed to create Metal pipeline state: " << error.localizedDescription.UTF8String << std::endl;
//             return -1;
//         }

//         id<MTLBuffer> bufferA = [device newBufferWithBytes:inA.data()
//                                                     length:sizeof(float) * arrayLength
//                                                    options:MTLResourceStorageModeShared];
//         id<MTLBuffer> bufferB = [device newBufferWithBytes:inB.data()
//                                                     length:sizeof(float) * arrayLength
//                                                    options:MTLResourceStorageModeShared];
//         id<MTLBuffer> bufferResult = [device newBufferWithBytes:result.data()
//                                                          length:sizeof(float) * arrayLength
//                                                         options:MTLResourceStorageModeShared];

//         id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
//         id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];
//         [computeEncoder setComputePipelineState:pipelineState];
//         [computeEncoder setBuffer:bufferA offset:0 atIndex:0];
//         [computeEncoder setBuffer:bufferB offset:0 atIndex:1];
//         [computeEncoder setBuffer:bufferResult offset:0 atIndex:2];

//         MTLSize gridSize = MTLSizeMake(arrayLength, 1, 1);
//         NSUInteger threadGroupSize = pipelineState.maxTotalThreadsPerThreadgroup;
//         if (threadGroupSize > 128) {
//             threadGroupSize = 128;
//         }
//         MTLSize threadgroupSize = MTLSizeMake(threadGroupSize, 1, 1);

//         auto gpu_start = std::chrono::high_resolution_clock::now();
//         [computeEncoder dispatchThreads:gridSize threadsPerThreadgroup:threadgroupSize];
//         [computeEncoder endEncoding];

//         [commandBuffer commit];
//         [commandBuffer waitUntilCompleted];

//         auto gpu_end = std::chrono::high_resolution_clock::now();
//         std::chrono::duration<double> gpu_duration = gpu_end - gpu_start;

//         memcpy(result.data(), [bufferResult contents], sizeof(float) * arrayLength);

//         float gpu_sum = 0.0f;
//         for (int i = 0; i < arrayLength; ++i) {
//             gpu_sum += result[i];
//         }

//         std::cout << "CPU Sum: " << cpu_sum << ", Time: " << cpu_duration.count() << " seconds" << std::endl;
//         std::cout << "GPU Sum: " << gpu_sum << ", Time: " << gpu_duration.count() << " seconds" << std::endl;
//     }

//     return 0;
// }