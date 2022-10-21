//
// multiply two arrays and store them in another array, for OpenCL
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                               Hiro KITAYAMA
//
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif //__APPLE__

#include <stdio.h>
#include <stdlib.h>

void verify(const int n, const float* a, const float *x, const float *y);

// get platform, device
int
getPlatFormDevideID(cl_device_type device_type,
    cl_platform_id* platformId, cl_device_id* deviceID)
{
    char message[1024];
    cl_uint numOfPlatforms;
    int rval = -1;

    // get list of platform
    cl_int status = clGetPlatformIDs(0, NULL, &numOfPlatforms);
    if (status != CL_SUCCESS || numOfPlatforms < 1)
    {
        fprintf(stderr, "clGetPlatformIDs function failed.\n");
        return -1;
    }

    cl_platform_id *platforms = (cl_platform_id*)malloc(sizeof(cl_platform_id)*numOfPlatforms);
    status = clGetPlatformIDs(numOfPlatforms, platforms, &numOfPlatforms);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clGetDeviceIDs function failed.\n");
        free(platforms);
        return -1;
    }

    for (unsigned plt = 0; plt < numOfPlatforms; plt++)
    {
        status = clGetPlatformInfo(platforms[plt], CL_PLATFORM_VERSION,
            sizeof(message), message, NULL);
        if (status != CL_SUCCESS)
        {
            fprintf(stderr, "clGetPlatformInfo function failed.\n");
            rval = -1;
            break;
        }
        fprintf(stdout, "platform: %s\n", message);

        cl_device_id deviceId[10];
        cl_uint numOfDevices;

        status = clGetDeviceIDs(platforms[plt], device_type,
            sizeof(deviceId) / sizeof(deviceId[0]), deviceId, &numOfDevices);
        if (status != CL_SUCCESS)
        {
            fprintf(stderr, "clGetDeviceIDs function failed.\n");
            rval = -1;
            break;
        }
        if (numOfDevices > 0)
        {
            clGetDeviceInfo(deviceId[0],
                CL_DEVICE_NAME, sizeof(message), message, NULL);
            fprintf(stdout, "device  : [%s]\n\n", message);

            *platformId = platforms[plt];
            *deviceID = deviceId[0];
            rval = 0;
            break;
        }
    }
    free(platforms);

    return rval;
}

#define N   4096

// main
int
main()
{
    cl_int status;
    cl_platform_id platformId;
    cl_device_id deviceID;

    float a[N], b[N], c[N];

    // initialize array
    for (int i = 0; i < N; i++)
    {
        a[i] = (float)(i + 1000);
        b[i] = (float)i / 10.f;
    }

    // get platform and device id
    if (getPlatFormDevideID(CL_DEVICE_TYPE_GPU, &platformId, &deviceID) < 0)
    {
        fprintf(stderr, "no opencl 2.x platform.");
        return -1;
    }

    // create Context
    cl_context context = clCreateContext(NULL, 1, &deviceID, NULL, NULL, NULL);

    // create Command Queue 2.x
    cl_command_queue queue = clCreateCommandQueueWithProperties(
        context, deviceID, NULL, NULL);

    // create program object
    static const char *src[] =
    {
        "__kernel void\n\
        mul(__global const float *a,\n\
            __global const float *b,\n\
            __global float *c)\n\
        {\n\
            int i=get_global_id(0);\n\
            c[i] = a[i] * b[i];\n\
        }\n"
    };
    cl_program prog = clCreateProgramWithSource(context,
        1, (const char**)&src, NULL, NULL);

    // build program
    const char* options = "-cl-std=CL2.0";
    status = clBuildProgram(prog, 1, &deviceID, options, NULL, NULL);

    // create kernel
    cl_kernel kernel = clCreateKernel(prog, "mul", NULL);

    // create memory object
    cl_mem mem_a = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
        sizeof(a), a, NULL);
    cl_mem mem_b = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
        sizeof(b), b, NULL);
    cl_mem mem_c = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(c), NULL, NULL);

    // set kernel parameters
    status = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&mem_a);
    status = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_b);
    status = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_c);

    // request execute kernel
    size_t globalSize[] = { sizeof(c) / sizeof(c[0]) };
    status = clEnqueueNDRangeKernel(queue, kernel, 1, NULL,
        globalSize, 0, 0, NULL, NULL);

    // get results
    status = clEnqueueReadBuffer(queue, mem_c, CL_TRUE, 0,
        sizeof(c), c, 0, NULL, NULL);

    // list results
    printf("(a * b = c)\n");
    for (int i = 0; i < 10; i++)
        printf("%f * %f = %f\n", a[i], b[i], c[i]);

    // flush queue
    status = clFlush(queue);

    // release resources
    clReleaseMemObject(mem_c);
    clReleaseMemObject(mem_b);
    clReleaseMemObject(mem_a);
    clReleaseKernel(kernel);
    clReleaseProgram(prog);
    clReleaseCommandQueue(queue);
    clReleaseContext(context);

    verify(N, a, b, c);

    return 0;
}
