#include <metal_stdlib>
using namespace metal;

kernel void add_vectors(const device float* inA [[buffer(0)]],
                        const device float* inB [[buffer(1)]],
                        device float* result [[buffer(2)]],
                        uint id [[thread_position_in_grid]]) {
    result[id] = inA[id] * inB[id] * inB[id];
}