#ifndef INCL_OPENCL
#define INCL_OPENCL

#include <vector>
#include <iomanip>
#include <string>
#include "cl.hpp"
struct OPENCL_kernel{
  std::vector<cl::Platform> all_platforms;
  cl::Platform platform;

  std::vector<cl::Device> all_devices;
  cl::Device device;

  std::string kernel_code;
  std::string func_name;
  cl::Context context;
  cl::Program::Sources sources;
  cl::Program program;
  cl::CommandQueue queue;/// vector is not needed??
  cl::Kernel kernel;
  std::vector<cl::Buffer> buffers;
  
  OPENCL_kernel(const std::string& kernel_codeIN, const std::string& func_nameIN):
    kernel_code(kernel_codeIN),
    func_name(func_nameIN)
  {
    /////////////////////////////////////
    // make platform object in a vector
    /////////////////////////////////////
    cl::Platform::get(&all_platforms);    
    for(const auto &pfm :all_platforms){
      std::cout << std::setw(40) << "CL_PLATFORM_EXTENSIONS: " << std::setw(30) << pfm.getInfo<CL_PLATFORM_EXTENSIONS>() << std::endl;
      std::cout << std::setw(40) << "CL_PLATFORM_NAME: "       << std::setw(30) << pfm.getInfo<CL_PLATFORM_NAME>() << std::endl;
      std::cout << std::setw(40) << "CL_PLATFORM_PROFILE: "    << std::setw(30) << pfm.getInfo<CL_PLATFORM_PROFILE>() << std::endl;
      std::cout << std::setw(40) << "CL_PLATFORM_VENDOR: "     << std::setw(30) << pfm.getInfo<CL_PLATFORM_VENDOR>() << std::endl;
      std::cout << std::setw(40) << "CL_PLATFORM_VERSIONN: "   << std::setw(30) << pfm.getInfo<CL_PLATFORM_VERSION>() << std::endl;    
    }
    platform = all_platforms[0]/* this is default platform */;
    /////////////////////////////////////
    // get all devices in vector of the default platform
    /////////////////////////////////////
    platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
    if(all_devices.size()==0){
      std::cout << std::setw(30) <<" No devices found. Check OpenCL installation!\n";
	exit(1);   
      }else
      for(const auto &dev :all_devices)
	{
	  std::cout << std::setw(40) << "CL_DEVICE_NAME: "                     << std::setw(30) << dev.getInfo<CL_DEVICE_NAME>() << std::endl;
	  std::cout << std::setw(40) << "CL_DEVICE_MAX_CLOCK_FREQUENCY: "      << std::setw(30) << dev.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << std::endl;
	  std::cout << std::setw(40) << "CL_DEVICE_MAX_WORK_ITEM_SIZES: "      << std::setw(30) << dev.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>() << std::endl; 
	  std::cout << std::setw(40) << "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS: " << std::setw(30) << dev.getInfo<CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>() << std::endl;
	  std::cout << std::setw(40) << "CL_DEVICE_MAX_WRITE_IMAGE_ARGS: "     << std::setw(30) << dev.getInfo<CL_DEVICE_MAX_WRITE_IMAGE_ARGS>() << std::endl;      
	}
    cl::Device device = all_devices[1];
    std::cout << "-> CL_DEVICE_NAME: " << device.getInfo<CL_DEVICE_NAME>() << " is chosen" << std::endl;
    ///////////////////////////////////////
    // construct program <- context and (source <- string), then buid
    ///////////////////////////////////////
    context = cl::Context({device});
    sources.push_back({kernel_code.c_str(), kernel_code.length()});
    program = cl::Program(context, sources);
    if(program.build({device}) != CL_SUCCESS){
      std::cout << "Error building: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device) << std::endl;
      exit(1);
    }else{
      std::cout << std::setw(40) << "CL_PROGRAM_KERNEL_NAMES: " << std::setw(30) << program.getInfo<CL_PROGRAM_KERNEL_NAMES>() << std::endl;
      std::cout << std::setw(40) << "CL_PROGRAM_SOURCE: "       << std::setw(30) << program.getInfo<CL_PROGRAM_SOURCE>() << std::endl;
    }
    ///////////////////////////////////////
    // set up kernels and vectors for GPU code
    ///////////////////////////////////////
    queue = cl::CommandQueue(context, device);
    
    // set up kernels and vectors for GPU code
    kernel = cl::Kernel(program, func_name.c_str());
    
    std::cout << std::setw(40) << "CL_KERNEL_ATTRIBUTES: "    << std::setw(30) << kernel.getInfo<CL_KERNEL_ATTRIBUTES>() << std::endl;
    std::cout << std::setw(40) << "CL_KERNEL_FUNCTION_NAME: " << std::setw(30) << kernel.getInfo<CL_KERNEL_FUNCTION_NAME>() << std::endl;
    
  };
  //========================================
  void enqueueNDRangeKernel(const int offset, const int work_item, const int num){
    queue.enqueueNDRangeKernel(kernel, cl::NDRange(offset), cl::NDRange(work_item), cl::NDRange(num));
  };
  //========================================  
  void enqueueNDRangeKernel(const int work_item, const int num){
    queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(work_item), cl::NDRange(num));
  };
  //========================================
  template<typename T>
  void setArg(const int i, const T v){
    kernel.setArg(i, v);
  };
  //========================================
  template<typename T>
  void makeBuffer(const int n, const std::vector<T>& vec){
    buffers.resize(n+1);     
    buffers[n]=cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(T)*(vec.size()));
    // queue.enqueueWriteBuffer(buffers[n], CL_TRUE, 0, sizeof(T)*(vec.size()), &vec[0]);
  };
  //========================================
  template<typename T>
  void enqueueWriteBuffer(const int n, const std::vector<T>& vec){
    // buffers.resize(n+1);
    // buffers[n]=cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(T)*(vec.size()));
    queue.enqueueWriteBuffer(buffers[n], CL_TRUE, 0, sizeof(T)*(vec.size()), &vec[0]);
  };
  //========================================
  template<typename T>
  void enqueueReadBuffer(const int n, std::vector<T>& vec){
    queue.enqueueReadBuffer(buffers[n], CL_TRUE, 0, sizeof(T)*(vec.size()), &vec[0]);
  };
};

#endif
