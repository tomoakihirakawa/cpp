#define NS_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#include <Metal/Metal.hpp>
#include <iostream>

int main() {
   // Initialize Metal device
   id<MTLDevice> device = MTLCreateSystemDefaultDevice();
   if (!device) {
      std::cerr << "Failed to create Metal device" << std::endl;
      return 1;
   }

   // Create Metal command queue
   id<MTLCommandQueue> commandQueue = [device newCommandQueue];
   if (!commandQueue) {
      std::cerr << "Failed to create command queue" << std::endl;
      return 1;
   }

   // Create Metal command buffer
   id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
   if (!commandBuffer) {
      std::cerr << "Failed to create command buffer" << std::endl;
      return 1;
   }

   // Create Metal render pass descriptor
   MTLRenderPassDescriptor *renderPassDescriptor = [MTLRenderPassDescriptor renderPassDescriptor];
   if (!renderPassDescriptor) {
      std::cerr << "Failed to create render pass descriptor" << std::endl;
      return 1;
   }

   // Set up render pass descriptor for clearing the screen
   renderPassDescriptor.colorAttachments[0].texture = nil;  // Assuming screen is the default drawable
   renderPassDescriptor.colorAttachments[0].loadAction = MTLLoadActionClear;
   renderPassDescriptor.colorAttachments[0].clearColor = MTLClearColorMake(0.0, 0.0, 1.0, 1.0);  // Clear to blue

   // Create Metal render command encoder
   id<MTLRenderCommandEncoder> renderEncoder = [commandBuffer renderCommandEncoderWithDescriptor:renderPassDescriptor];
   if (!renderEncoder) {
      std::cerr << "Failed to create render command encoder" << std::endl;
      return 1;
   }

   // End encoding commands
   [renderEncoder endEncoding];

   // Commit the command buffer for execution
   [commandBuffer commit];
   [commandBuffer waitUntilCompleted];

   std::cout << "Screen cleared to blue using Metal!" << std::endl;

   return 0;
}
