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

#include "scw.h"

#include "GNUPLOT3.h"
#include "fundamental.h"

//#include "cl.hpp"

// #include "interp2parm.h"
// #include "EarthScience.h"

// std::vector<std::vector<double>> base_spline_matrix(const std::vector<double>& h, const int s, const int K){
//   std::vector<std::vector<double>> ret(h.size(),std::vector<double>(s,0.));
//   std::vector<double> q = base_spline_node(h, K); 
//   for(auto i = 0; i<h.size(); i++)
//     for(auto j = 0; j<s; j++)
//       ret[i][j] = base_spline(h[i], q, j, K);
//   return ret;
// };
// std::vector<std::vector<double>> D_base_spline_matrix(const std::vector<double>& h, const int s, const int K, const int n){
//   std::vector<std::vector<double>> ret(h.size(),std::vector<double>(s,0.));
//   std::vector<double> q = base_spline_node(h, K); 
//   for(auto i = 0; i<h.size(); i++)
//     for(auto j = 0; j<s; j++)
//       ret[i][j] = D_base_spline(h[i], q, j, K, n);
//   return ret;
// };


// struct Interpolation{
// public:
//   int K;  
//   int s;
//   std::vector<double> samp2, q;
//   int s1, s2;
//   std::vector<double> q1, q2;
//   std::vector<std::vector<double>> samp3;
//   std::vector<std::vector<double>> inv_coefmat, inv_coefmatT;
//   /////////////////////////// 1 parm  ////////////////////////////
//   Interpolation(const std::vector<double>& samp2_IN,
// 		const std::vector<double>& samp2_x_IN,
// 		const int K_IN):
//     samp2(samp2_IN), s(samp2_IN.size()), K(K_IN), q(base_spline_node(samp2_x_IN, K_IN)),
//     inv_coefmat(Inverse(base_spline_matrix(samp2_x_IN, samp2_x_IN.size(), K_IN, base_spline_node(samp2_x_IN, K_IN)))){};
//   std::vector<std::vector<double>> shape(const std::vector<double>& a){return Dot(base_spline_matrix(a, s, K, q),inv_coefmat);};
//   std::vector<double> operator()(const std::vector<double>& a){return Dot(shape(a), samp2);};  
//   std::vector<double> shape(const double a){return Dot(base_spline_vector(a, s, K, q),inv_coefmat);};
//   double operator()(const double a){return Dot(shape(a), samp2);};
//   //============================================
//   std::vector<std::vector<double>> D_shape(const std::vector<double>& a, const int n){return Dot(D_base_spline_matrix(a, s, K, n, q),inv_coefmat);};  
//   std::vector<double> operator()(const std::vector<double>& a, const int n){return Dot(D_shape(a,n), samp2);};
//   std::vector<double> D_shape(const double a, const int n){return Dot(D_base_spline_vector(a, s, K, n, q),inv_coefmat);};  
//   double operator()(const double a, const int n){return Dot(D_shape(a, n), samp2);};
//   ////////////////////////// 2 parm  ///////////////////////////
//   Interpolation(const std::vector<std::vector<double>>& samp3_IN,
//   		const std::vector<std::vector<double>>& samp3_x_IN,
//   		const std::vector<std::vector<double>>& samp3_y_IN,
//   		const int K_IN):
//     samp3(samp3_IN),
//     s1(samp3_IN.size()),
//     s2(samp3_IN[0].size()),
//     K(K_IN),
//     q1(samp3_x_IN.size()),
//     q2(samp3_y_IN.size()),
//     inv_coefmat(samp3_IN.size(),std::vector<double>(samp3_IN.size(),0.)),
//     inv_coefmatT(samp3_IN.size(),std::vector<double>(samp3_IN.size(),0.))
//   {
//     double min_x = Min(samp3_x_IN);
//     double max_x = Max(samp3_x_IN);
//     double min_y = Min(samp3_y_IN);
//     double max_y = Max(samp3_y_IN);

//     int sx = samp3_x_IN.size();
//     int sy = samp3_y_IN.size();
    
//     q1 = base_spline_node(Subdivide(min_x, max_x, sx-1), K_IN);
//     q2 = base_spline_node(Subdivide(min_y, max_y, sy-1), K_IN);

//     inv_coefmatT = Inverse(Transpose(base_spline_matrix(samp3_x_IN, sx, K_IN, q1)));
//     inv_coefmat = Inverse(Transpose(base_spline_matrix(samp3_y_IN, sy, K_IN, q2)));

//     std::cout << inv_coefmat << std::endl;        
//     std::cout << inv_coefmatT << std::endl;
    
//   };
//   // std::vector<std::vector<double>> shape1(const std::vector<double>& a){
//   //   std::vector<std::vector<double>> coefmat_a = base_spline_matrix(a, s1, K, q1);
//   //   return Dot(coefmat_a ,inv_coefmat);
//   // };
//   // std::vector<std::vector<double>> shape2(const std::vector<double>& b){
//   //   std::vector<std::vector<double>> coefmat_bT = Transpose(base_spline_matrix(b, s2, K, q2));// s x b.size()
//   //   return Dot(inv_coefmatT ,coefmat_bT);
//   // };
//   std::vector<double> shape1(double a){
//     std::vector<double> coefmat_a = base_spline_vector(a, s1, K, q1);
//     return Dot(coefmat_a ,inv_coefmat);
//   };
//   std::vector<double> shape2(double b){
//     std::vector<double> coefmat_bT = base_spline_vector(b, s2, K, q2);
//     return Dot(inv_coefmatT ,coefmat_bT);
//   };
//   // std::vector<double> shape1(const double a){
//   //   std::vector<double> coefmat_a = base_spline_vector(a, s1, K);
//   //   return Dot(coefmat_a ,inv_coefmat);
//   // };
//   // std::vector<double> shape2(const double b){
//   //   std::vector<double> coefmat_bT = base_spline_vector(b, s2, K);
//   //   return Dot(inv_coefmatT ,coefmat_bT);
//   // };
//   //  std::vector<std::vector<double>> operator()(const std::vector<double>& a, const std::vector<double>& b){return Dot(shape1(a),Dot(samp3,shape2(b)));};  
//   double operator()(const double a, const double b){return Dot(shape1(a),Dot(shape2(b),samp3));};
//   //=============================================
//   // std::vector<std::vector<double>> D_shape1(const std::vector<double>& a, const int n){
//   //   std::vector<std::vector<double>> coefmat_a = D_base_spline_matrix(a, s1, K, n, q);
//   //   return Dot(coefmat_a ,inv_coefmat);
//   // };
//   // std::vector<std::vector<double>> D_shape2(const std::vector<double>& b, const int n){
//   //   std::vector<std::vector<double>> coefmat_bT = Transpose(D_base_spline_matrix(b, s2, K, n, q));// s x b.size()
//   //   return Dot(inv_coefmatT ,coefmat_bT);
//   // };    
//   // std::vector<double> D_shape1(const double a, const int n){
//   //   std::vector<double> coefmat_a = D_base_spline_vector(a, s1, K, n);
//   //   return Dot(coefmat_a ,inv_coefmat);
//   // };
//   // std::vector<double> D_shape2(const double b, const int n){
//   //   std::vector<double> coefmat_bT = D_base_spline_vector(b, s2, K, n);
//   //   return Dot(inv_coefmatT ,coefmat_bT);
//   // };
//   // std::vector<std::vector<double>> operator()(const std::vector<double>& a, const std::vector<double>& b, const std::vector<int>&n){
//   //   return Dot(D_shape1(a, n[0]),Dot(samp3,D_shape2(b, n[1])));};
//   // double operator()(const double a, const double b, const std::vector<int> &n){
//   //   return Dot(D_shape1(a,n[0]),Dot(samp3,D_shape2(b,n[1])));};  
// };
// //==================================================
// void Load(const std::string& filename, std::vector<std::vector<double>> &mat, const std::string &DIV){
//   std::ifstream strm(filename, std::ios::in);
//   if(!strm)
//     std::cout << "The file can not be opened" << std::endl;  
//   std::string str;
//   int s(0), e(0), len(0);
//   double v;
//   for(auto & MAT: mat)
//     {
//       s = 0;
//       std::getline(strm,str);
//       for(auto & VEC: MAT)
// 	{
// 	  if((e = str.find(DIV))!=std::string::npos)
// 	    {
// 	      v = (double) stod(str.substr(s,e));
// 	      VEC = v;
// 	      s = e;
// 	    }
// 	  else
// 	    {
// 	      v = (double) stod(str);
// 	      VEC = v;	
// 	    }
// 	}
//     }
//   strm.close();
// };

// //////////////////////////////////////////////////////////////////////////////////////////
// template<typename T>
// struct OPENCL_buffer{
// public:
//   cl::Context context;
//   std::vector<T> vec;
//   std::vector<T> vec_read;
//   std::vector<T> vec_write;  
//   cl::Buffer buffer;
//   bool b_write, b_read;
//   OPENCL_buffer(const cl::Context& contextIN, const std::vector<T>& vecIN, cl::CommandQueue &queue, bool b_writeIN, bool b_readIN):
//     context(contextIN), vec(vecIN), vec_write(vecIN), vec_read(vecIN), b_write(b_writeIN), b_read(b_readIN)
//   {
//     cl::Buffer buffer(context, CL_MEM_READ_WRITE, sizeof(T)*(vec.size()) );//make buffer space
//     if(b_write)
//       queue.enqueueWriteBuffer(buffer, CL_TRUE, 0, sizeof(T)*(vec_write.size()), &vec_write);
//     if(b_read)
//       queue.enqueueReadBuffer(buffer, CL_TRUE, 0, sizeof(T)*(vec_read.size()), &vec_read);
//   };
//   OPENCL_buffer(const cl::Context& contextIN, const std::vector<T>& vecIN):
//     context(contextIN), vec(vecIN), vec_write(vecIN), vec_read(vecIN), b_write(false), b_read(false)
//   {
//     cl::Buffer buffer(context, CL_MEM_READ_WRITE, sizeof(T)*(vec.size()) );//make buffer space
//   };
//   void set_buffer(const std::vector<T>& vecIN)
//   {
//     cl::Buffer buffer(context, CL_MEM_READ_WRITE, sizeof(T)*(vec.size()) );//make buffer space
//   };
//   void write(cl::CommandQueue &queue)
//   {
//     queue.enqueueWriteBuffer(buffer, CL_TRUE, 0, sizeof(T)*(vec_write.size()), &vec_write);
//   };
//   void read(cl::CommandQueue &queue)
//   {
//     queue.enqueueReadBuffer(buffer, CL_TRUE, 0, sizeof(T)*(vec_read.size()), &vec_read);
//   };
// };
// struct OPENCL{
// public:
//   std::string func_name;
//   std::vector<cl::Platform> all_platforms;
//   std::vector<cl::Device> all_devices;
//   cl::Platform platform;
//   cl::Device device;//chosen device
//   cl::Context context;
//   std::string kernel_code;
//   cl::Program::Sources sources;
//   cl::Program program;
//   cl::CommandQueue queue;/// vector is not needed??
//   std::map<std::string, cl::Kernel> kernels;
  
//   // arguments
//   std::vector<cl::Buffer> buff_vec;
//   std::vector<int> size_vec;
//   std::vector<bool> write_buff;

//   OPENCL(const std::string& kernel_codeIN): all_platforms(0), all_devices(0), kernel_code(kernel_codeIN)
//   {
//     /////////////////////////////////////
//     // make platform object in a vector
//     /////////////////////////////////////
//     cl::Platform::get(&all_platforms);    
//     for(const auto &pfm :all_platforms)
//       {
// 	std::cout << std::setw(40) << "CL_PLATFORM_EXTENSIONS: " << std::setw(30) << pfm.getInfo<CL_PLATFORM_EXTENSIONS>() << std::endl;
// 	std::cout << std::setw(40) << "CL_PLATFORM_NAME: "       << std::setw(30) << pfm.getInfo<CL_PLATFORM_NAME>() << std::endl;
// 	std::cout << std::setw(40) << "CL_PLATFORM_PROFILE: "    << std::setw(30) << pfm.getInfo<CL_PLATFORM_PROFILE>() << std::endl;
// 	std::cout << std::setw(40) << "CL_PLATFORM_VENDOR: "     << std::setw(30) << pfm.getInfo<CL_PLATFORM_VENDOR>() << std::endl;
// 	std::cout << std::setw(40) << "CL_PLATFORM_VERSIONN: "   << std::setw(30) << pfm.getInfo<CL_PLATFORM_VERSION>() << std::endl;    
//       }
//     platform = all_platforms[0]/* this is default platform */;
//     /////////////////////////////////////
//     // get all devices in vector of the default platform
//     /////////////////////////////////////
//     platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
//     // exception
//     if(all_devices.size()==0)
//       {
// 	std::cout << std::setw(30) <<" No devices found. Check OpenCL installation!\n";
// 	exit(1);   
//       }
//     else
//       for(const auto &dev :all_devices)
// 	{
// 	  std::cout << std::setw(40) << "CL_DEVICE_NAME: " << std::setw(30) << dev.getInfo<CL_DEVICE_NAME>() << std::endl;
// 	  std::cout << std::setw(40) << "CL_DEVICE_MAX_CLOCK_FREQUENCY: " << std::setw(30) << dev.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << std::endl;
// 	  std::cout << std::setw(40) << "CL_DEVICE_MAX_WORK_ITEM_SIZES: " << std::setw(30) << dev.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>() << std::endl; 
// 	  std::cout << std::setw(40) << "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS: " << std::setw(30) << dev.getInfo<CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>() << std::endl;
// 	  std::cout << std::setw(40) << "CL_DEVICE_MAX_WRITE_IMAGE_ARGS: " << std::setw(30) << dev.getInfo<CL_DEVICE_MAX_WRITE_IMAGE_ARGS>() << std::endl;      
// 	}
//     // use device[1] because that's a GPU; device[0] is the CPU
//     device = all_devices[1];
//     std::cout << "-> CL_DEVICE_NAME: "<< device.getInfo<CL_DEVICE_NAME>() << " is chosen" << std::endl;
//     ///////////////////////////////////////
//     // construct program <- context and (source <- string), then buid
//     ///////////////////////////////////////
//     context = {device};
//     sources.push_back({kernel_code.c_str(), kernel_code.length()});
//     cl::Program program(context, sources);
//     if(program.build({device}) != CL_SUCCESS)
//       {
// 	std::cout << "Error building: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device) << std::endl;
// 	exit(1);
//       }
//     ///////////////////////////////////////
//     // set up kernels and vectors for GPU code
//     ///////////////////////////////////////
//     cl::CommandQueue queue(context, device);
//   };

//   void ActivateKernel(const std::string& func_nameIN, const int num_of_args/* determine the size of buff_vec */)
//   {/* a function name which is written in input kernel_code */
//     // kernels.insert(
//     // 		   make_pair(
//     // 			     func_name,
//     // 			     cl::Kernel(program, func_name.c_str())
//     // 			     )
//     // 		   );
//   };

//   template<typename T>
//   void SetArg(const std::string& func_name,
// 	      const std::vector<T>& vecIN
// 	      const int n)
//   {
//     OPENCL_buffer buf(context, vecIN);
//     buff_vec[n] = buf;
    
//     // // exception
//     // if(!buff_vecIN.size()==size_vecIN.size() || !buff_vecIN.size()==write_vecIN.size()){
//     //   std::cout << "size different" << std::endl;
//     //   abort();
//     // };
    
//     // buff_vec.clear();
//     // buff_vec = buff_vecIN;

//     // size_vec.clear();
//     // size_vec = size_vecIN;

//     // write_vec.clear();
//     // write_vec = write_vecIN;

//     // int i(0);
//     // for(const auto& buff: buff_vec)
//     //   {
//     // 	kernels.at(func_name).setArg(i, buff);
//     // 	if(write_buff[i])
//     // 	  enqueueWriteBuffer(buff, CL_TRUE, 0, sizeof(double)*size_vec[i], A);
//     // 	i++;
//     //   }
//   };

  // void Execute(const std::string& func_name, const int num_workitem, const int num_local)
  // {// ex. OPENCL("add")
  //   queue.enqueueNDRangeKernel(kernels.at(func_name),
  // 			       cl::NullRange,  // kernel, offset
  // 			       cl::NDRange(num_workitem), // global number of work items
  // 			       cl::NDRange(num_local)); // local number (per group)
  // };
  
//   void Get(const int i, std::vector<double>& ret)
//   {
//     ret.resize(size_vec[i]);
//     queue.enqueueReadBuffer(buff_vec[i], CL_TRUE, 0, sizeof(double)*size_vec[i], &ret);
//   };

//   void Clear(){
//     queue.finish();
//   };
  
//   ~OPENCL(){
//     queue.finish();
//   };
 
// };
// //////////////////////////////////////////////////////////////////////////////////////////
///////////////
#define pi 3.1415926535897932385

using namespace std;

int main()
{



  std::vector<double> vec{1,2,3,4,5,6,7,8};
  std::cout << vec << std::endl;

  std::cout << vec << std::endl;
  
  //  std::string dir = "/Users/tomoaki/Desktop/sideE.dat";


  // std::string file = "/Users/tomoaki/Desktop/top.csv";
  // std::vector<std::vector<std::vector<double>>> vvv;
  // Load(file,vvv);
  // cout << vvv << endl;

  // GNUPLOT PLOT;
  // PLOT.SaveSplotData(vvv,{{"w","l"}});
  // PLOT.Splot();
  // std::cin.ignore();

  // cout << vvv << endl;
  // std::vector<std::vector<double>> mode = Load(dir,{"}","{","} ","{ "," ",",",", ",",  ",",   "});
  // cout << mode << endl;

  // std::vector<std::vector<int>> v{{1,2,3,4,4},{1,2,3,4,4},{1,2,3,433333,4}};
  // std::vector<std::vector<int>> w{{1,2,3,4,2,3,4,5,4},{1,2,3,4,2,3,4,5,4},{1,2,3,4,2,3,4,5,2,24,3,53,53,53,4}};

  // for(auto &f: w[0])
  //   cout << f << endl;

  // w.clear();
  // w=v;
  // cout << "srga" << endl;
  // for(auto &f: w[2])
  //   cout << f << endl;

  // if(false)
  //   {
  //     std::string dir = "/Users/tomoaki/Dropbox/research/scw_3d/workshop/result/50_6/2d0000EP00_1d0000EP02/1d0000000EM02/0d0000000EP00_2d0051028EP00/mode.dat";
  //     std::vector<std::vector<double>> mode = Load(dir,{",",", ",",  ",",   "});
  //     dir = "/Users/tomoaki/Dropbox/research/scw_3d/workshop/result/50_6/2d0000EP00_1d0000EP02/1d0000000EM02/0d0000000EP00_2d0051028EP00/data.dat";
  //     std::vector<std::vector<double>> data = Load(dir,{",",", ",",  ",",   "});
  //     std::vector<double> data0 = Flatten(data);
  
  //     scw obj(data0,mode);
  //     int max = data0[0];
  
  //     int M = 50, N = 50;
  //     double x, y, z;

  //     std::vector<std::vector<std::vector<double>>> mat;
  //     for(auto i = 0; i < M; i++)
  // 	{
  // 	  std::vector<std::vector<double>> row_vec;
  // 	  for(auto j = 0; j < N; j++)
  // 	    {     
  // 	      x = pi*(2.*i/(M-1.) - 1.);
  // 	      y = pi*(2.*j/(N-1.) - 1.);
  // 	      z = obj.PHI(x,y,0);
  // 	      row_vec.push_back({x,y,z});	  
  // 	    }
  // 	  std::cout << row_vec << std::endl;
  // 	  mat.push_back(row_vec);
  // 	}
  //     GNUPLOT p;
  //     p.Splot(mat,{{"w","l"}});

  //     std::cin.ignore();

  //   }



  
  // if(false)
  //   {
  //     std::vector<std::vector<double>> vec(100, std::vector<double>(3,0));
  
  //     int N = 30/*y*/;
  //     int M = 30/*x*/;  
  //     std::vector<std::vector<double>> x(N,std::vector<double>(M,0)), y(N,std::vector<double>(M,0)), z(N,std::vector<double>(M,0));
  //     std::vector<std::vector<double>> dzdy(z), dzdx(z);
  //     for(auto i = 0; i < M; i++){
  // 	for(auto j = 0; j < N; j++){
  // 	  double r = (double)i/M;	  
  // 	  double theta = M_PI * (2.*j/(M-1.) - 1.);
  // 	  x[j][i] = r*cos(theta);
  // 	  y[j][i] = r*sin(theta);
  // 	  z[j][i] = r;
  // 	}
  //     }
      
  //     ParametricInterpolation interpx(x,4);
  //     ParametricInterpolation interpy(y,4);
  //     ParametricInterpolation interpz(z,4);      

  //     GNUPLOT PLOT;
  //     PLOT.Set({{"hidden3d",""}});
  //     PLOT.Set({{"title","\'PLOT\'"}});
      
  //     // for(auto ii=0; ii<M; ii++)
  // 	// for(auto jj=0; jj<N; jj++)	
  //     //  	  {
  // 	    int div=100;
  // 	    std::vector<std::vector<std::vector<double>>> mat;
  // 	    for(auto j=0; j<div; j++)
  // 	      {
  // 		double eta = (double)(2.*j/(div-1.) - 1.);
  // 		std::vector<std::vector<double>> row_vector;      
  // 		for(auto i=0; i<div; i++)
  // 		  {
  // 		    double xi = (double)(2.*i/(div-1.) - 1.);
  // 		    //if(!((ii==3||ii==4)&&(jj==3||jj==4)))
  // 		    row_vector.push_back({interpx({xi,eta}),interpy({xi,eta}),interpz({xi,eta})});
  // 		  }
  // 		mat.push_back(row_vector);
  // 	      }
  // 	    PLOT.SaveSplotData(mat,{{"w","l"}});
  // 	    //}
  
  //     PLOT.Splot();
  //     std::cin.ignore();
      
  //   }


  
  // auto base_spline_ = [](const double h, const std::vector<double> &q, const int i, const int K)
  // 		      {return (q[i] <= h && h < q[i+1]) ? 1. : 0.;};
    

    //  std::cout << lambda(1) << std::endl;
  
  //////////////////////////
  // //subdivideで行列を作成
  // GNUPLOT PLOT;
  // int N = 100;
  // std::vector<std::vector<std::vector<double>>> data;  
  // for(auto i=0; i<N; i++)
  //   {
  //     std::vector<std::vector<double>> subdata;	  
  //     double x = 2.*(2.*(((double)i - N/2.)));
  //     for(auto j=0; j<N; j++)    
  // 	{
  // 	  double y = 2.*(2.*(((double)j - N/2.)));
  // 	  std::vector<double> xyz{x,y,0.};  
  // 	  std::vector<double> u1{0., 0.0068,0.};
  // 	  std::vector<double> u2{0., 0.0046,0.};
  // 	  std::vector<double> u3{0., 0.0068,0.};
  // 	  subdata.push_back(xyz + IntegrateStrikeDipOkada(u1,{x,y},{-16.1039, 37.7521},{35., 35.},{-15., 20.},13.,1.,1.)
  // 			    + IntegrateStrikeDipOkada(u2,{x,y},{-23.7321, 8.88264},{35., 35.},{15., 20.},13.,1.,1.)
  // 			    + IntegrateStrikeDipOkada(u3,{x,y},{-33.0555, -25.5369},{35., 35.},{15., 20.},13.,1.,1.));
	  
  // 	}
  //     data.push_back(subdata);      
  //   }
  // int counter = 1;
  // PLOT.Set({{"contour","base"}});
  // PLOT.Splot(data,{{"w","l"}});
  // PLOT.Clear();











  
  // std::vector<double> x(10), z(10);
  // std::vector<std::vector<double>> data2d;
  // for(auto i = 0; i<x.size(); i++){
  //   x[i] = (double)4.*(2.*i/(x.size()-1.) - 1.);
  //   z[i] = sin(x[i]);
  //   data2d.push_back({x[i],z[i]});
  // }

  // GNUPLOT PLOT;
  // PLOT.SaveData(data2d);

  // cout << x << endl;
  
  // Interpolation interp_zx(z, x, 3);

  // int N = 100;
  // std::vector<std::vector<double>> datavec;
  // for(int i=0; i<N; i++){
  //   double h = (double)3.*i/(N-1.) - 1.;
  //   datavec.push_back({h, interp_zx(h)});
  // }

  // PLOT.SaveData(datavec);
  // PLOT.Plot2D();
  /////////////////////////////////////////////

  // cout << QuotientMod(1,7) << endl;
  
  // std::vector<int> v{1,2,3};
  // std::vector<int> v1{1,2,1};  
  // cout << (v==v1) << endl;
  // cout << MemberQ(std::vector<std::vector<int>>{{0,0},{1,618},{1,618}},std::vector<int> {1,618}) << endl;
  // cout << MemberQ({{0,0},{1,618},{1,618}},{1,618}) << endl;
  // cout << MemberQ({{0,0},{1,618},{1,618}},{1}) << endl;
  // cout << MemberQ({1,618},1) << endl;      
  // abort();


  
  // ///////////////////////////////////////////////////////////////////////////////////////
  // // BiBspline is working properly
  // std::vector<std::vector<double>> topo=Import("/Users/tomoaki/Dropbox/research/tsunami/mathematica/data.asc");
  // std::vector<std::vector<double>> z=Take(topo,{0,(int)topo.size()-1,57},{0,(int)topo[0].size()-1,48});

  // cout << Take(topo,{0,(int)topo.size()-1,57},{0,(int)topo[0].size()-1,48})[0].size() << endl;  
  
  // //  std::vector<std::vector<double>> z(7,std::vector<double>(7,0.));
  // std::vector<std::vector<double>> x(z), y(z);
  
  // std::vector<std::vector<std::vector<double>>> data3d, data3d2;  
  // for(auto i = 0; i<x.size(); i++){
  //   std::vector<std::vector<double>> data2d, data2d2;    
  //   for(auto j = 0; j<x[i].size(); j++){    
  //     x[i][j] = (double)10.*(2.*i/(x.size()-1.) - 1.);
  //     y[i][j] = (double)10.*(2.*j/(y.size()-1.) - 1.);
  //     //      z[i][j] = sin(x[i][j] + y[i][j]);
  //     data2d.push_back({x[i][j],y[i][j]});
  //     data2d2.push_back({x[i][j],y[i][j],z[i][j]});
  //   }
  //   data3d.push_back(data2d);
  //   data3d2.push_back(data2d2);    
  // }
  
  // BiBspline13_3 method;
  
  // int N = 50;
  // std::vector<std::vector<std::vector<double>>> data;  
  // for(int i=0; i<N; i++){
  //   std::vector<std::vector<double>> data2d;    
  //   for(int j=0; j<N; j++){
  //     double h = (double)(2.*i/(N-1.) - 1.);
  //     double r = (double)(2.*j/(N-1.) - 1.);
  //     data2d.push_back({h, -r, method.interp(z,{r,h})});
  //   }
  //   data.push_back(data2d);
  // }
  
  // GNUPLOT PLOT;
  // PLOT.Set({{"contour","base"}});
  // PLOT.Splot(data,{{"w","l"}});

  // std::cin.ignore();
  
  //////////////////////////////////////////////////

  // std::vector<int> S{10,50,100,500,1000,5000};
  // for(auto & s: S)
  //   {
  //     std::vector<double> x(s), w(s);
  //     gauleg(0.,1,x,w);
  //     double sum(0.);
  //     for(auto i=0; i<x.size(); i++)
  // 	sum += log(x[i])*w[i];
  //     std::cout << s  << " " << std::setprecision(4) << sum + 1. << std::endl;
  //   }

  // double beta=2.;
  // double x_sing=0.;
  // std::vector<int> S{10,50,100,500,1000,5000};
  // for(auto & s: S)
  //   {
  //     std::vector<double> x(s), w(s);
  //     gauleg(InvSg(0,x_sing,beta),InvSg(1,x_sing,beta),x,w);      
  //     double sum(0.);
  //     for(auto i=0; i<x.size(); i++)
  // 	sum += log(Sg(x[i],x_sing,beta)) * w[i] * DSg(x[i],x_sing,beta);
  //     std::cout << s  << " " << std::setprecision(4) << sum + 1. << std::endl;
  //   }
  //////////////////////////////  
  // std::vector<double> Beta{1., 1.5, 2., 2.5, 3.5, 4.};

  // std::vector<double> alpha{1., 1.5, 2., 2.5, 3.5, 4.};

  // alpha *= Beta;
  
  // std::cout << alpha << std::endl;
  
  // for(auto &beta: Beta)
  // {
  //   std::cout << "-------------------" << std::endl;
  //   std::cout << "<beta = "<< beta << ">" << std::endl;
  //   double x_sing=0.;
  //   std::vector<int> S{10,20,30,40};
  //   for(auto & s: S)
  //     {
  // 	std::vector<double> x(s), w(s);
  // 	gauleg(0., InvSg(1,x_sing,beta),x,w);
  // 	double sum(0.);
  // 	for(auto i=0; i<x.size(); i++)
  // 	  {
  // 	    // sum += log(abs(Sg(x[i],x_sing,beta))) * w[i] * DSg(x[i],x_sing,beta);
  // 	    sum += 1/sqrt(abs(Sg(x[i],x_sing,beta))) * w[i] * DSg(x[i],x_sing,beta);
  // 	  }
  // 	std::cout << std::setw(5) << s << "   " << std::scientific << std::setprecision(2) << sum - 2.   << std::endl;
  //     }
  // }
  //////
  // std::vector<double> Beta{1., 1.2, 1.4, 1.6, 1.8, 2.};
  // for(auto &beta: Beta)
  // {
  //   std::cout << "-------------------" << std::endl;
  //   std::cout << "<beta = "<< beta << ">" << std::endl;
  //   double x_sing=0.;
  //   std::vector<int> S{10,50,100,500};
  //   for(auto & s: S)
  //     {
  // 	std::vector<double> x(s), w(s);
  // 	gauleg(InvSg(-1,x_sing,beta), InvSg(1.,x_sing,beta), x, w);      
  // 	double sum(0.);
  // 	for(auto i=0; i<x.size(); i++)
  // 	  {
  // 	    sum += log(Sg(x[i], x_sing, beta)) * w[i] * DSg(x[i], x_sing, beta);
  // 	    //sum += 1/sqrt(Sg(x[i],x_sing,beta)) * w[i] * DSg(x[i],x_sing,beta);
  // 	  }
  // 	std::cout << std::setw(5) << s << "   " << std::scientific << std::setprecision(2) << sum + 1. << std::endl;
  //     }    
  // }  
  ///////////////////////////////////////////////////////////////////////////////////////
  // std::vector<double> Beta{1., 1.2, 1.4, 1.6, 1.8, 2.};
  // for(auto &beta: Beta)
  // {
  //   std::cout << "-------------------" << std::endl;
  //   std::cout << "<beta = "<< beta << ">" << std::endl;
  //   double x_sing=0.;
  //   std::vector<int> S{10,50,100,500};
  //   for(auto & s: S)
  //     {
  // 	std::vector<double> x(s), w(s);
  // 	gauleg(-1E+1,log(1),x,w);
  // 	double sum(0.);
  // 	for(auto i=0; i<x.size(); i++)
  // 	  {
  // 	    sum += log(exp(x[i])) * w[i] * exp(x[i]);
  // 	    //sum += 1/sqrt(Sg(x[i],x_sing,beta)) * w[i] * DSg(x[i],x_sing,beta);
  // 	  }
  // 	std::cout << std::setw(5) << s << "   " << std::scientific << std::setprecision(2) << sum << std::endl;
  //     }    
  // }  
  //////////////////////////////
  // std::vector<double> Beta{1., 1.2, 1.4, 1.6, 1.8, 2.};
  // for(auto &beta: Beta)
  // {
  //   std::cout << "-------------------" << std::endl;
  //   std::cout << "<beta = "<< beta << ">" << std::endl;
  //   double x_sing=0.;
  //   std::vector<int> S{10,50,100,500};
  //   for(auto & s: S)
  //     {
  // 	std::vector<double> x(s), w(s);
  // 	gauleg(-1E+1,0.,x,w);      
  // 	double sum(0.);
  // 	for(auto i=0; i<x.size(); i++)
  // 	  {
  // 	    sum += log(abs(SgLog(x[i],x_sing))) * w[i] * DSgLog(x[i],x_sing);
  // 	    //sum += 1/sqrt(Sg(x[i],x_sing,beta)) * w[i] * DSg(x[i],x_sing,beta);
  // 	  }
  // 	std::cout << std::setw(5) << s << "   " << std::scientific << std::setprecision(2) << sum << std::endl;
  //     }    
  // }
  
  //////////////////////////////
  // std::vector<double> Hrange=   {1.2,1.26,1.32,1.38,1.44,1.5,1.56,1.62,1.68,1.74,1.8,1.86,1.92,1.98,2.04,2.1,2.16,2.22,2.28,2.34,2.4,2.46,2.52,2.58,2.64,2.7,2.76,2.82,2.88,2.94,3.,3.06,3.12,3.18,3.24,3.3,3.36,3.42,3.48,3.54,3.6,3.66,3.72,3.78,3.84,3.9,3.96,4.02,4.08,4.14,4.2};
  // std::vector<double> hLrange=  {0.035,0.0373,0.0396,0.0419,0.0442,0.0465,0.0488,0.0511,0.0534,0.0557,0.058,0.0603,0.0626,0.0649,0.0672,0.0695,0.0718,0.0741,0.0764,0.0787,0.081,0.0833,0.0856,0.0879,0.0902,0.0925,0.0948,0.0971,0.0994,0.1017,0.104,0.1063,0.1086,0.1109,0.1132,0.1155,0.1178,0.1201,0.1224,0.1247,0.127,0.1293,0.1316,0.1339,0.1362,0.1385,0.1408,0.1431,0.1454,0.1477,0.15};
  // std::vector<double> Lrange=  {0.7,0.738,0.776,0.814,0.852,0.89,0.928,0.966,1.004,1.042,1.08,1.118,1.156,1.194,1.232,1.27,1.308,1.346,1.384,1.422,1.46,1.498,1.536,1.574,1.612,1.65,1.688,1.726,1.764,1.802,1.84,1.878,1.916,1.954,1.992,2.03,2.068,2.106,2.144,2.182,2.22,2.258,2.296,2.334,2.372,2.41,2.448,2.486,2.524,2.562,2.6};
  // std::vector<double> Trange= {0.75,0.791,0.832,0.873,0.914,0.955,0.996,1.037,1.078,1.119,1.16,1.201,1.242,1.283,1.324,1.365,1.406,1.447,1.488,1.529,1.57,1.611,1.652,1.693,1.734,1.775,1.816,1.857,1.898,1.939,1.98,2.021,2.062,2.103,2.144,2.185,2.226,2.267,2.308,2.349,2.39,2.431,2.472,2.513,2.554,2.595,2.636,2.677,2.718,2.759,2.8};
    
  // std::vector<double> wx = {13.0744,59.6451,442.13,203.366,237.076,448.993,567.377,66.4711,142.221,228.549,127.73,715.386,231.932,615.922,34.787,133.242,158.234,102.14,450.35,267.827,670.213,23.2585,40.3444,69.2853,513.483,296.505,660.17,331.335,20.,38.2329,123.107,211.941,152.493,799.079,538.09,22.5763,100.322,413.352,562.681,115.419,280.271,283.225,6.36541,22.6323,85.9154,187.331,179.185,178.704,464.456,9.85771,23.5323,50.6554,96.2859,254.258,218.763,286.513,10.5451,10.5204,69.4681,253.317,560.217,280.2,160.703,6.97749,6.43085,24.6743,156.606,273.148,195.733,165.647,10.1162,5.83113,12.8169,62.5477,106.825,261.093,203.076};
  // std::vector<double> A = {0.223355,2.96027,0.480487,1.68838,2.53657,0.519238,0.603888,0.906279,0.205537,2.55722,2.12171,3.33323,0.511006,0.744374,2.73626,1.07587,2.76606,0.759013,2.46303,0.309348,2.58345,1.22965,0.243398,0.190388,0.491817,4.03177,0.655419,4.23352,1.53287,0.0629484,0.494069,0.275508,1.7916,0.896713,3.96509,3.99313,0.209686,0.396889,1.44635,0.711903,0.308246,3.35744,4.23385,1.87269,3.70568,0.215013,0.63793,1.33749,1.26718,2.2295,0.130964,0.0718069,0.133145,0.243027,1.60134,1.36898,0.0203072,2.87735,0.547377,0.271447,0.866459,2.94698,1.17654,0.643726,0.12053,0.0999192,0.110892,0.268783,0.306583,0.330919,0.151479,0.925811,1.38347,0.101317,0.390085,0.384652,1.37503};

  // std::vector<std::vector<double>> TH{{2.92571,0.203688},{2.92571,0.480345},{0.853333,2.18951},{1.07789,3.47918},{1.46286,3.04674},{0.853333,4.0194},{2.27556,4.78303},{0.758519,0.625976},{0.744727,0.850579},{0.999024,2.30256},{1.20471,2.96881},{0.8192,3.86987},{2.048,4.13772},{2.40941,4.61522},{0.999024,0.725535},{0.999024,1.45239},{0.975238,2.48128},{1.41241,2.79873},{0.787692,3.71403},{2.048,4.52566},{2.56,5.44434},{1.20471,0.802013},{1.20471,1.26906},{1.24121,2.01429},{1.57538,2.13626},{1.86182,3.41026},{2.15579,3.79502},{2.56,4.46143},{1.46286,0.602042},{1.46286,0.989318},{1.51704,1.74252},{1.46286,2.75097},{2.15579,2.73725},{2.27556,4.11973},{2.73067,4.4936},{1.70667,0.449244},{1.70667,1.1314},{1.70667,1.69904},{1.70667,2.15752},{2.27556,2.97169},{2.56,3.30059},{2.92571,4.72064},{1.86182,0.511566},{1.95048,1.13837},{1.95048,1.9504},{1.86182,2.53187},{2.048,3.55615},{2.92571,3.59937},{2.73067,5.14692},{2.15579,0.475986},{2.048,0.947512},{2.15579,1.70159},{2.27556,2.33178},{2.27556,3.10026},{2.92571,3.82964},{1.70667,4.06801},{2.40941,0.458722},{2.27556,1.06934},{2.40941,1.54574},{2.40941,2.23939},{2.40941,2.91526},{2.92571,3.72463},{1.86182,3.83707},{2.73067,0.473304},{2.73067,0.985823},{2.73067,1.68428},{2.73067,2.60111},{2.56,3.37264},{2.92571,3.42413},{1.95048,3.0777},{2.73067,0.308385},{2.73067,0.684938},{2.73067,1.34801},{2.92571,2.06167},{2.73067,3.20432},{2.92571,4.07058},{2.92571,3.22165}};
  // std::vector<std::vector<double>> LH{{2.4381,0.203688},{2.92571,0.480345},{0.761905,2.18951},{1.00738,3.47918},{1.43417,3.04674},{0.836601,4.0194},{2.20928,4.78303},{0.648306,0.625976},{0.636519,0.850579},{0.891986,2.30256},{1.12589,2.96881},{0.765607,3.86987},{1.95048,4.13772},{2.31674,4.61522},{0.900022,0.725535},{0.916536,1.45239},{0.88658,2.48128},{1.34516,2.79873},{0.722653,3.71403},{1.93208,4.52566},{2.56,5.44434},{1.11547,0.802013},{1.11547,1.26906},{1.18211,2.01429},{1.45869,2.13626},{1.80759,3.41026},{2.05313,3.79502},{2.4381,4.46143},{1.36716,0.602042},{1.36716,0.989318},{1.41779,1.74252},{1.40659,2.75097},{2.07287,2.73725},{2.23094,4.11973},{2.70363,4.4936},{1.64103,0.449244},{1.64103,1.1314},{1.65696,1.69904},{1.6254,2.15752},{2.12669,2.97169},{2.46154,3.30059},{2.86835,4.72064},{1.74002,0.511566},{1.8576,1.13837},{1.8576,1.9504},{1.75643,2.53187},{1.95048,3.55615},{2.8405,3.59937},{2.60063,5.14692},{2.11352,0.475986},{1.93208,0.947512},{2.07287,1.70159},{2.1672,2.33178},{2.20928,3.10026},{2.86835,3.82964},{1.68977,4.06801},{2.33923,0.458722},{2.1672,1.06934},{2.33923,1.54574},{2.29468,2.23939},{2.31674,2.91526},{2.81319,3.72463},{1.75643,3.83707},{2.5761,0.473304},{2.55202,0.985823},{2.60063,1.68428},{2.62564,2.60111},{2.41509,3.37264},{2.81319,3.42413},{1.87546,3.0777},{2.60063,0.308385},{2.73067,0.684938},{2.70363,1.34801},{2.78639,2.06167},{2.62564,3.20432},{2.76011,4.07058},{2.76011,3.22165}};
  // std::vector<std::vector<double>> hLH{{2.05078,0.203688},{1.70898,0.480345},{6.5625,2.18951},{4.96338,3.47918},{3.48633,3.04674},{5.97656,4.0194},{2.26318,4.78303},{7.7124,0.625976},{7.85522,0.850579},{5.60547,2.30256},{4.44092,2.96881},{6.53076,3.86987},{2.56348,4.13772},{2.1582,4.61522},{5.55542,0.725535},{5.45532,1.45239},{5.63965,2.48128},{3.71704,2.79873},{6.91895,3.71403},{2.58789,4.52566},{1.95313,5.44434},{4.48242,0.802013},{4.48242,1.26906},{4.22974,2.01429},{3.42773,2.13626},{2.76611,3.41026},{2.4353,3.79502},{2.05078,4.46143},{3.65723,0.602042},{3.65723,0.989318},{3.52661,1.74252},{3.55469,2.75097},{2.41211,2.73725},{2.24121,4.11973},{1.84937,4.4936},{3.04688,0.449244},{3.04688,1.1314},{3.01758,1.69904},{3.07617,2.15752},{2.35107,2.97169},{2.03125,3.30059},{1.74316,4.72064},{2.87354,0.511566},{2.69165,1.13837},{2.69165,1.9504},{2.84668,2.53187},{2.56348,3.55615},{1.76025,3.59937},{1.92261,5.14692},{2.36572,0.475986},{2.58789,0.947512},{2.41211,1.70159},{2.30713,2.33178},{2.26318,3.10026},{1.74316,3.82964},{2.95898,4.06801},{2.13745,0.458722},{2.30713,1.06934},{2.13745,1.54574},{2.17896,2.23939},{2.1582,2.91526},{1.77734,3.72463},{2.84668,3.83707},{1.94092,0.473304},{1.95923,0.985823},{1.92261,1.68428},{1.9043,2.60111},{2.07031,3.37264},{1.77734,3.42413},{2.66602,3.0777},{1.92261,0.308385},{1.83105,0.684938},{1.84937,1.34801},{1.79443,2.06167},{1.9043,3.20432},{1.81152,4.07058},{1.81152,3.22165}};
  
  // std::vector<double> x = Lrange;//change
  // std::vector<double> y = Hrange;//change

  // RBF_inversemultiquadric mult(1.);
  // RBF_interp interp(LH, A, mult, 0);
  
  // std::vector<std::vector<std::vector<double>>> data;  
  // for(int i=0; i<x.size(); i++){
  //   std::vector<std::vector<double>> data2d;    
  //   for(int j=0; j<y.size(); j++){
  //     data2d.push_back({x[i], y[j], interp({x[i], y[j]}) });
  //   }    
  //   data.push_back(data2d);
  // }

  // cout << data << endl;
  
  // GNUPLOT PLOT2;  
  // PLOT2.Set({{"contour","base"},{"title","\'PLOT\'"}});  
  // PLOT2.Set({{"hidden3d",""},{"title","\"interpolated\""}});
  // PLOT2.Splot(data,{{"w","l"}});

  // std::cin.ignore();
  
  ////////////////////////////////////////////////////////////////////////////////////////
  
  // int N = 5;
  // int M = 5;  
  // std::vector<std::vector<double>> x(N,std::vector<double>(M,0)), y(N,std::vector<double>(M,0)), z(N,std::vector<double>(M,0));
  // std::vector<std::vector<double>> dzdy(z), dzdx(z);
  // for(auto i = 0; i < M; i++){
  //   for(auto j = 0; j < N; j++){    
  //     // x[j][i] = pi*(2.*i/(M-1.) - 1.);
  //     // y[j][i] = pi*(2.*j/(N-1.) - 1.);

  //     x[j][i] = (i==0) ? -pi : pi*(2.*i/(M-1.) - 1.);
  //     y[j][i] = (i==0) ? 0. : pi*(2.*j/(N-1.) - 1.);
  //     z[j][i] = (i==0) ? 0. : cos(x[j][i])*cos(y[j][i]);
  //     dzdy[j][i] = -pi*cos(x[j][i])*sin(y[j][i]);
  //     dzdx[j][i] = -pi*sin(x[j][i])*cos(y[j][i]);      
  //   }
  // }

  // ParametricInterpolation interpx(x,3);
  // ParametricInterpolation interpy(y,3);
  // ParametricInterpolation interpz(z,3);
  
  // ParametricInterpolation interp_dzdx(dzdx,3);  
  // ParametricInterpolation interp_dzdy(dzdy,3);

  // GNUPLOT PLOT;
  // N = 100;
  // M = 100;  
  // PLOT.Set({{"contour","base"},{"title","\'PLOT\'"}});  
  // for(auto n=0; n<2; n++){
  //   std::vector<std::vector<std::vector<double>>> data;
  //   for(auto i=0; i<M; i++){
  //     std::vector<std::vector<double>> data2d;
  //     for(auto j=0; j<N; j++)
  // 	{
  // 	  double h = (double)(2.*i/(M-1.) - 1.);
  // 	  double r = (double)(2.*j/(N-1.) - 1.);
  // 	  if(n==1)
  // 	    data2d.push_back({interpx({r,h},{0,0}), interpy({r,h},{0,0}), interpz({r,h},{0,0})});
  // 	  else
  // 	    {}
  // 	    //	    data2d.push_back({h, r, interp_dzdy({r,h},{0,0})});
  // 	}
  //     data.push_back(data2d);
  //   }
  //   PLOT.SaveSplotData(data,{{"w","l"}});    
  // }  

  // PLOT.Splot();
  // //std::cin.ignore();
  
  ///////////////////////////////////////////////////////////////////////////////////////

  // std::vector<std::vector<double>> topo=Import("/Users/tomoaki/Dropbox/research/tsunami/mathematica/data.asc");
  // std::vector<std::vector<double>> z=Take(topo,{0,(int)topo.size()-1,40},{0,(int)topo[0].size()-1,40});
  // std::vector<std::vector<double>> x(z), y(z);
  
  // std::vector<std::vector<std::vector<double>>> data3d, data3d2;  
  // for(auto i = 0; i<x.size(); i++){
  //   std::vector<std::vector<double>> data2d, data2d2;    
  //   for(auto j = 0; j<x[i].size(); j++){    
  //     x[i][j] = (double)10.*(2.*i/(x.size()-1.) - 1.);
  //     y[i][j] = (double)10.*(2.*j/(y.size()-1.) - 1.);
  //     data2d.push_back({x[i][j],y[i][j]});
  //     data2d2.push_back({x[i][j],y[i][j],z[i][j]});
  //   }
  //   data3d.push_back(data2d);
  //   data3d2.push_back(data2d2);    
  // }

  
  // ParametricInterpolation interp(z,3);

  
  // int N = 50;
  // std::vector<std::vector<std::vector<double>>> data;  
  // for(int i=0; i<N; i++){
  //   std::vector<std::vector<double>> data2d;    
  //   for(int j=0; j<N; j++){
  //     double h = (double)(2.*i/(N-1.) - 1.);
  //     double r = (double)(2.*j/(N-1.) - 1.);
  //     data2d.push_back({h, -r, interp(r,h)});
  //   }
  //   data.push_back(data2d);
  // }
  
  // GNUPLOT PLOT;
  // PLOT.Set({{"contour","base"}});
  // PLOT.Splot(data,{{"w","l"}});

  // std::cin.ignore();
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  // std::vector<std::vector<double>> topo=Import("/Users/tomoaki/Dropbox/research/tsunami/mathematica/data.asc");
  // std::vector<std::vector<double>> z=Take(topo,{0,(int)topo.size()-1,30},{0,(int)topo[0].size()-1,30});
  // std::vector<std::vector<double>> x(z), y(z);

  // std::vector<std::vector<std::vector<double>>> data3d, data3d2;
  // for(auto i = 0; i<x.size(); i++){
  //   std::vector<std::vector<double>> data2d, data2d2;    
  //   for(auto j = 0; j<x[i].size(); j++){    
  //     x[i][j] = (double)10.*(2.*i/(x.size()-1.) - 1.);
  //     y[i][j] = (double)10.*(2.*j/(y.size()-1.) - 1.);
  //     data2d.push_back({x[i][j],y[i][j]});
  //     data2d2.push_back({x[i][j],y[i][j],z[i][j]});
  //   }
  //   data3d.push_back(data2d);
  //   data3d2.push_back(data2d2);    
  // }

  // RBF_multiquadric mult(0.);
  // RBF_interp interp(Flatten(data3d),Flatten(z),mult,0);
    
  // // GNUPLOT PLOT;
  // // PLOT.Set({{"title","\"given\""}});  
  // // PLOT.Set({{"contour","base"}});
  // // PLOT.SaveSplotData(data3d2,{{"w","l"}});   
  // // PLOT.Splot();
  
  // int N = 50;
  // std::vector<std::vector<std::vector<double>>> data;  
  // for(int i=0; i<N; i++){
  //   std::vector<std::vector<double>> data2d;    
  //   for(int j=0; j<N; j++){
  //     double h = (double)20.*(2.*i/(N-1.) - 1.);
  //     double r = (double)20.*(2.*j/(N-1.) - 1.);
  //     data2d.push_back({h, r, interp({h,r})});
  //   }    
  //   data.push_back(data2d);
  // }


  ///////////////////////////////////////////////////////////////////////////////////////
  
  // GNUPLOT PLOT2;  
  // PLOT.Set({{"hidden3d",""},{"title","\"interpolated\""}});
  // PLOT.Set({{"contour","base"}});
  // PLOT.Splot(data,{{"w","l"}});
  
  // int N = 100, M = 100;
  // std::vector<std::vector<double>> topo0=Import("/Users/tomoaki/Dropbox/research/tsunami/mathematica/data.asc");
  // std::vector<std::vector<double>> topo=Take(topo0,{0,670,10},{0,570,10});
  // std::vector<std::vector<double>> topo_x(topo), topo_y(topo);
  // for(auto i = 0; i<topo.size(); i++){
  //   for(auto j = 0; j<topo[0].size(); j++){
  //     topo_x[i][j] = (double)250.*(2.*j/(topo[0].size()-1.) - 1.);
  //     topo_y[i][j] = (double)250.*(2.*i/(topo.size()-1.) - 1.);
  //   }
  // }
  // ParametricInterpolation interp_z(topo, 3);
  // ParametricInterpolation interp_x(topo_x, 3);
  // ParametricInterpolation interp_y(topo_y, 3);
  // std::vector<std::vector<std::vector<double>>> datamat;
  // for(auto i=0; i<N; i++){
  //     std::vector<std::vector<double>> datavec;
  //     for(auto j=0; j<M; j++){
  // 	  double x = (double)2.*i/(N-1.) - 1.;
  // 	  double y = (double)2.*j/(M-1.) - 1.;
  // 	  datavec.push_back({interp_x(x,y), interp_y(x,y), interp_z(x,y)});
  // 	}
  //     datamat.push_back(datavec);
  //   }
  
  // GNUPLOT PLOT;
  // PLOT.Set({{"hidden3d",""},{"title","\"interpolated\""},{"contour","base"}});
  // PLOT.SaveSplotData(datamat,{{"w","l"}});
  // PLOT.Splot();

  // std::cin.ignore();
  
  // PLOT.MatrixPlot();
  // PLOT.Clear();
  // cout << Dot(
  // 	       Dot(coefmat_a ,coefmat),
  // 	       Dot(x, Dot(coefmatT,coefmat_b))
  // 	       ) << endl;
  
  // for(int i=0; i<N; i++)
  //   {
  //     interp.push_back(Dot(f,base_spline_interp(x, size, 3)));
  //     interp.push_back(Dot(f,base_spline_interp(y, size, 3)));
  //     interp.push_back(Dot(f,base_spline_interp(z, size, 3)));      
  //   }
  // cout << interp << endl;
  // cout << f << endl;  
  // GNUPLOT plot;  
  // plot.Plot1D(interp);
  
  // std::vector<std::vector<double>> mat0=Import("/Users/tomoaki/Dropbox/research/tsunami/mathematica/data.asc");
  // std::vector<std::vector<double>> mat=Take(mat0,{0,720,10},{0,600,10});    
  // GNUPLOT PLOT;
  // PLOT.MatrixPlot(mat,{{"w","image"}}); 
  // PLOT.Clear();
  
  // ifstream fin("/Users/tomoaki/Dropbox/research/tsunami/mathematica/data.asc");  
  // int max_x = 600;
  // int max_y = 720;
  // int i=0, counter=0;
  // std::vector<std::vector<double>> mat;
  // for(int y=0; y<max_y; y++)
  //   {
  //     std::vector<double> row(max_x);
  //     for(int x=0; x<max_x; x++)	
  // 	fin>>row[x];//>> stops when it encounters a space
  //     mat.push_back(row);
  //   }
  // fin.close();
  // GNUPLOT PLOT;
  // PLOT.SaveMatrixData(mat,{{"w","image"}});
  // PLOT.MatrixPlot(); 
  // PLOT.Clear();

  
  // cout << vec.size() << endl;

  // vector<vector<double>>  mat{{1,2,3},{1,2,3},{1,2,3}};
  // cout << Flatten(mat) << endl;

  // vector<vector<vector<double>>>  mat{{{1,2,3},{1,2,3},{1,2,3}},{{1,2,3},{1,2,1000}},{{1,2,3},{1,2,3},{1,2,3}}};
  // cout << Flatten(mat) << endl;
  
  
  // //subdivideで行列を作成
  // int M = 200;
  // for(int i=0; i<M; i++)
  //   {
  //     double t = (i - M/2.);
  //     double w = 0.1;
  //     std::vector<std::vector<std::vector<double>>> data;
  //     int N = 100;
  //     for(int i=0; i<N; i++)
  // 	{
  // 	  std::vector<std::vector<double>> subdata;	  
  // 	  double x = 2.*(((double)i - N/2.));
  // 	  for(int j=0; j<N; j++)    
  // 	    {
  // 	      double y = 2.*(((double)j - N/2.));
  // 	      std::vector<double> xyz{x,y,0.};

	      
  // 	      // std::vector<double> u1{0., 0.0068,0.};
  // 	      // std::vector<double> u2{0., 0.0046,0.};
  // 	      // std::vector<double> u3{0., 0.0068,0.};
  // 	      // u1 = u1*(1+tanh(w*t))/2.;
  // 	      // u2 = u2*(1+tanh(w*t))/2.;
  // 	      // u3 = u3*(1+tanh(w*t))/2.;
  // 	      // subdata.push_back(xyz + IntegrateStrikeDipOkada(u1,{x,y},{-16.1039, 37.7521},{35., 35.},{-15., 20.},13.,1.,1.)
  // 	      // 		        + IntegrateStrikeDipOkada(u2,{x,y},{-23.7321, 8.88264},{35., 35.},{15., 20.},13.,1.,1.)
  // 	      // 		        + IntegrateStrikeDipOkada(u3,{x,y},{-33.0555, -25.5369},{35., 35.},{15., 20.},13.,1.,1.));

  // 	      std::vector<double> u1{0.00054, 0.003};
  // 	      std::vector<double> u2{0., 0.0076};
  // 	      u1 = u1*(1+tanh(w*t))/2.;
  // 	      u2 = u2*(1+tanh(w*t))/2.;
  // 	      subdata.push_back(xyz + IntegrateStrikeDipOkada(u1,{x,y},{-19.4942, 4.4413},{60, 30},{-5., 25.},3.,1.,1.)
  // 	      			+ IntegrateStrikeDipOkada(u2,{x,y},{-34.7516, -32.1985},{40, 30},{22., 40.},2.,1.,1.));
	      
  // 	    }
  // 	  data.push_back(subdata);      
  // 	}
  //     //cout << data << endl;
  //     GNUPLOT PLOT;
  //     int counter = 1;
  //     PLOT.Set(PLOT.def_eps_set);
  //     PLOT.Set({{"output","\'~/Desktop/matuda/__" + to_string(i)+".eps\'"}});
  //     PLOT.Set({{"pm3d","at s"}});
  //     PLOT.Set({{"hidden3d",""}});
  //     PLOT.Set({{"zrange","[-0.005:0.006]"}});
  //     PLOT.Set({{"title","\'{/Helvetica-Oblique=24 t} = " + to_string(t) + "\'" }});      
  //     PLOT.Set({{"contour","base"}});
  //     PLOT.Set({{"contour","both"}});
  //     PLOT.Set({{"view","40,30"}});                   
  //     PLOT.Splot(data,{{"w","l"}});
  //     PLOT.Clear();
  //   };
  
  // data.clear();
  // for(int i=0; i<N; i++)
  //   {
  //     std::vector<std::vector<double>> subdata;	  
  //     double x = 3*(((double)i - N/2.));
  //     for(int j=0; j<N; j++)    
  // 	{
  // 	  double y = 3*(((double)j - N/2.));
  // 	  std::vector<double> xyz{x,y,0.};
  // 	  subdata.push_back(xyz + IntegrateStrikeDipOkada({0.00054, 0.003},{x,y},{-19.4942, 4.4413},{60, 30},{-5., 25.},3.,1.,1.)
  // 			    + IntegrateStrikeDipOkada({.0, 0.0076},{x,y},{-34.7516, -32.1985},{40, 30},{22., 40.},2.,1.,1.));
  // 	}
  //     data.push_back(subdata);      
  //   }
  // cout << data << endl;
  // GNUPLOT PLOT2;
  // PLOT2.Set({{"title","\'{/Helvetica-Oblique=24 t} = " + to_string(real_t) + "\'" }});      
  // PLOT2.Set(PLOT.def_eps_set);
  // PLOT2.Set({{"output","\'./Desktop/" + name + "_" + to_string(counter)+".eps\'"}});
  // PLOT2.Set({{"hidden3d",""}}); 
  // PLOT2.Splot(data,{{"w","l"}});
  // PLOT2.Clear();
  
  // for(int i=0; i<10000; i++)  
  //   cout << StrikeDip({1, 2, 0.}, std::vector<double>{1, 2, 0.} - std::vector<double>{2, 1, 0.}, {1, 2, 0.}, {1,2} , 1., 1., 1.) << endl;

  // int len = 13;
  // std::vector<std::vector<std::vector<double>>> samp3;
  // std::vector<std::vector<double>> sampx(len,std::vector<double>(len,0.));
  // std::vector<std::vector<double>> sampy(len,std::vector<double>(len,0.));
  // std::vector<std::vector<double>> sampz(len,std::vector<double>(len,0.));  
  // for(int i=0; i<len; i++)
  //   {
  //     std::vector<std::vector<double>> subdata;	  
  //     for(int j=0; j<len; j++)
  // 	{	

  // 	  sampx[i][j]=i;
  // 	  sampy[i][j]=j;
  // 	  sampz[i][j]=sin(i*pi/len)*sin(j*pi/len);
	  
  // 	  std::vector<double> X{(double)i,(double)j,sampz[i][j]};
  // 	  subdata.push_back(X);
  // 	};
  //     samp3.push_back(subdata);      
  //   };
  // GNUPLOT PLOTsamp;
  // PLOTsamp.Set({{"hidden3d",""}});
  // PLOTsamp.Splot(samp3,{{"w","l"}});

  
  // std::vector<std::vector<std::vector<double>>> top{{{-3.,-3.,0.,0.},{-2.5,-3.,0.,1.},{-2.,-3.,0.,2.},{-1.5,-3.,0.,3.},{-1.,-3.,0.,4.},{-0.5,-3.,0.,5.},{0.,-3.,0.,6.},{0.5,-3.,0.,7.},{1.,-3.,0.,8.},{1.5,-3.,0.,9.},{2.,-3.,0.,10.},{2.5,-3.,0.,11.},{3.,-3.,0.,12.}},{{-3.,-2.5,0.,13.},{-2.5,-2.5,1.44271*1E-8,14.},{-2.,-2.5,5.48715*1E-8,15.},{-1.5,-2.5,1.72318*1E-7,16.},{-1.,-2.5,4.20277*1E-7,17.},{-0.5,-2.5,7.44866*1E-7,18.},{0.,-2.5,9.07999*1E-7,19.},{0.5,-2.5,7.44866*1E-7,20.},{1.,-2.5,4.20277*1E-7,21.},{1.5,-2.5,1.72318*1E-7,22.},{2.,-2.5,5.48715*1E-8,23.},{2.5,-2.5,1.44271*1E-8,24.},{3.,-2.5,0.,25.}},{{-3.,-2.,0.,26.},{-2.5,-2.,5.48715*1E-8,27.},{-2.,-2.,2.44089*1E-7,28.},{-1.5,-2.,9.07999*1E-7,29.},{-1.,-2.,2.60965*1E-6,30.},{-0.5,-2.,5.24501*1E-6,31.},{0.,-2.,6.70925*1E-6,32.},{0.5,-2.,5.24501*1E-6,33.},{1.,-2.,2.60965*1E-6,34.},{1.5,-2.,9.07999*1E-7,35.},{2.,-2.,2.44089*1E-7,36.},{2.5,-2.,5.48715*1E-8,37.},{3.,-2.,0.,38.}},{{-3.,-1.5,0.,39.},{-2.5,-1.5,1.72318*1E-7,40.},{-2.,-1.5,9.07999*1E-7,41.},{-1.5,-1.5,4.12971*1E-6,42.},{-1.,-1.5,0.0000147668,43.},{-0.5,-1.5,0.0000358351,44.},{0.,-1.5,0.0000495747,45.},{0.5,-1.5,0.0000358351,46.},{1.,-1.5,0.0000147668,47.},{1.5,-1.5,4.12971*1E-6,48.},{2.,-1.5,9.07999*1E-7,49.},{2.5,-1.5,1.72318*1E-7,50.},{3.,-1.5,0.,51.}},{{-3.,-1.,0.,52.},{-2.5,-1.,4.20277*1E-7,53.},{-2.,-1.,2.60965*1E-6,54.},{-1.5,-1.,0.0000147668,55.},{-1.,-1.,0.0000698689,56.},{-0.5,-1.,0.000228428,57.},{0.,-1.,0.00036619,58.},{0.5,-1.,0.000228428,59.},{1.,-1.,0.0000698689,60.},{1.5,-1.,0.0000147668,61.},{2.,-1.,2.60965*1E-6,62.},{2.5,-1.,4.20277*1E-7,63.},{3.,-1.,0.,64.}},{{-3.,-0.5,0.,65.},{-2.5,-0.5,7.44866*1E-7,66.},{-2.,-0.5,5.24501*1E-6,67.},{-1.5,-0.5,0.0000358351,68.},{-1.,-0.5,0.000228428,69.},{-0.5,-0.5,0.001178,70.},{0.,-0.5,0.00265802,71.},{0.5,-0.5,0.001178,72.},{1.,-0.5,0.000228428,73.},{1.5,-0.5,0.0000358351,74.},{2.,-0.5,5.24501*1E-6,75.},{2.5,-0.5,7.44866*1E-7,76.},{3.,-0.5,0.,77.}},{{-3.,0.,0.,78.},{-2.5,0.,9.07999*1E-7,79.},{-2.,0.,6.70925*1E-6,80.},{-1.5,0.,0.0000495747,81.},{-1.,0.,0.00036619,82.},{-0.5,0.,0.00265802,83.},{0.,0.,0.01,84.},{0.5,0.,0.00265802,85.},{1.,0.,0.00036619,86.},{1.5,0.,0.0000495747,87.},{2.,0.,6.70925*1E-6,88.},{2.5,0.,9.07999*1E-7,89.},{3.,0.,0.,90.}},{{-3.,0.5,0.,91.},{-2.5,0.5,7.44866*1E-7,92.},{-2.,0.5,5.24501*1E-6,93.},{-1.5,0.5,0.0000358351,94.},{-1.,0.5,0.000228428,95.},{-0.5,0.5,0.001178,96.},{0.,0.5,0.00265802,97.},{0.5,0.5,0.001178,98.},{1.,0.5,0.000228428,99.},{1.5,0.5,0.0000358351,100.},{2.,0.5,5.24501*1E-6,101.},{2.5,0.5,7.44866*1E-7,102.},{3.,0.5,0.,103.}},{{-3.,1.,0.,104.},{-2.5,1.,4.20277*1E-7,105.},{-2.,1.,2.60965*1E-6,106.},{-1.5,1.,0.0000147668,107.},{-1.,1.,0.0000698689,108.},{-0.5,1.,0.000228428,109.},{0.,1.,0.00036619,110.},{0.5,1.,0.000228428,111.},{1.,1.,0.0000698689,112.},{1.5,1.,0.0000147668,113.},{2.,1.,2.60965*1E-6,114.},{2.5,1.,4.20277*1E-7,115.},{3.,1.,0.,116.}},{{-3.,1.5,0.,117.},{-2.5,1.5,1.72318*1E-7,118.},{-2.,1.5,9.07999*1E-7,119.},{-1.5,1.5,4.12971*1E-6,120.},{-1.,1.5,0.0000147668,121.},{-0.5,1.5,0.0000358351,122.},{0.,1.5,0.0000495747,123.},{0.5,1.5,0.0000358351,124.},{1.,1.5,0.0000147668,125.},{1.5,1.5,4.12971*1E-6,126.},{2.,1.5,9.07999*1E-7,127.},{2.5,1.5,1.72318*1E-7,128.},{3.,1.5,0.,129.}},{{-3.,2.,0.,130.},{-2.5,2.,5.48715*1E-8,131.},{-2.,2.,2.44089*1E-7,132.},{-1.5,2.,9.07999*1E-7,133.},{-1.,2.,2.60965*1E-6,134.},{-0.5,2.,5.24501*1E-6,135.},{0.,2.,6.70925*1E-6,136.},{0.5,2.,5.24501*1E-6,137.},{1.,2.,2.60965*1E-6,138.},{1.5,2.,9.07999*1E-7,139.},{2.,2.,2.44089*1E-7,140.},{2.5,2.,5.48715*1E-8,141.},{3.,2.,0.,142.}},{{-3.,2.5,0.,143.},{-2.5,2.5,1.44271*1E-8,144.},{-2.,2.5,5.48715*1E-8,145.},{-1.5,2.5,1.72318*1E-7,146.},{-1.,2.5,4.20277*1E-7,147.},{-0.5,2.5,7.44866*1E-7,148.},{0.,2.5,9.07999*1E-7,149.},{0.5,2.5,7.44866*1E-7,150.},{1.,2.5,4.20277*1E-7,151.},{1.5,2.5,1.72318*1E-7,152.},{2.,2.5,5.48715*1E-8,153.},{2.5,2.5,1.44271*1E-8,154.},{3.,2.5,0.,155.}},{{-3.,3.,0.,156.},{-2.5,3.,0.,157.},{-2.,3.,0.,158.},{-1.5,3.,0.,159.},{-1.,3.,0.,160.},{-0.5,3.,0.,161.},{0.,3.,0.,162.},{0.5,3.,0.,163.},{1.,3.,0.,164.},{1.5,3.,0.,165.},{2.,3.,0.,166.},{2.5,3.,0.,167.},{3.,3.,0.,168.}}};
  // std::vector<std::vector<double>> X(top.size(), std::vector<double>(top[0].size(),0));
  // std::vector<std::vector<double>> Y(top.size(), std::vector<double>(top[0].size(),0));
  // std::vector<std::vector<double>> Z(top.size(), std::vector<double>(top[0].size(),0));
  // BiBspline13_3 spline;
  // double h, r;
  // int N=100;
  // for(auto i = 0; i<top.size(); i++){
  //   for(auto j = 0; j<top[0].size(); j++){
  //     X[i][j] = top[i][j][0];
  //     Y[i][j] = top[i][j][1];
  //     Z[i][j] = top[i][j][2];	
  //   }
  // }

  // std::vector<std::vector<std::vector<double>>> data;
  // for(int i=0; i<N; i++)
  //   {
  //     std::vector<std::vector<double>> subdata;
  //     h = 2.*i/(N-1.) - 1.;
  //     for(int j=0; j<N; j++)
  // 	{
  // 	  r = 2.*j/(N-1.) - 1.;
  // 	  double x = spline.interp(X,{h,r});
  // 	  double y = spline.interp(Y,{h,r});
  // 	  double z = spline.D_interp(Z,{h,r},{0,0}); 	    
  // 	  subdata.push_back({x,y,z});
  // 	}
  //     data.push_back(subdata);
  //   }
  // cout << data << endl;
  // GNUPLOT PLOT;
  // // PLOT.Set({{"hidden3d",""}});
  // PLOT.Splot(data,{{"w","l"}});
  // std::cin.ignore();

  
  // std::vector<std::vector<std::vector<double>>> data;
  // double x, y, z;
  // int N=20; 
  // for(int i=0; i<N; i++)
  //   {
  //     std::vector<std::vector<double>> subdata;	  
  //     x = (double)i - N/2;
  //     for(int j=0; j<N; j++)    
  // 	{
  // 	  y = (double)j - N/2;	
  // 	  z = pow(cosh((x+y)/(2.*pi)),-2);
  // 	  std::vector<double> X{x,y,z};
  // 	  subdata.push_back(X);
  // 	}
  //     data.push_back(subdata);      
  //   }
  // cout << data << endl;
  // GNUPLOT PLOT;
  // PLOT.Set({{"hidden3d",""}});
  // PLOT.Splot(data,{{"w","l"}});


  
  //  cout << Red << Cross(vector<double>{1,2,3},vector<double>{1,2,2}) << reset<< endl;
  
  //  std::vector<double> samp0{2,4,2,9,4};
  //  std::vector<std::vector<double>> data0(100,std::vector<double>(2,0.));  
  //  double h0;
  //  Bspline5_3 spline0;
  //  for(int i=0; i<100; i++){
  //    h0 = 2.*i/99. - 1.;
  //    data0[i] = {h0,spline0.interp(samp0,h0)};
  //  }
  //  GNUPLOT plo0;
  //  plo0.Plot2D(data0);


  //  std::vector<double> v1{1,2};
  //  std::vector<double> v2{1,2,3};  

  //  cout << Red << v1 - v2 << reset << endl;
  
  //  std::vector<std::vector<double>> samp{{1,2},{2,4},{3,2},{4,9},{5,4}};
  //  std::vector<std::vector<double>> samp2{{1,2},{2,4},{3,2},{4,9},{5,4},{5,4}};

  //  cout << samp - 2.*samp2 << endl;
  
  //  std::vector<std::vector<double>> data(100,std::vector<double>(2,0.));  
  //  double h;
  //  Bspline9_3 spline;
  //  for(int i=0; i<100; i++){
  //    h = 2*i/99. - 1.;
  //    data[i]=spline.interp(samp,h);
  //  }
  //  GNUPLOT plo;
  //  plo.Plot2D(data);
  
  //  cout << Cross({1,2}) << endl;

  //  cout << Cross({1,2,0},{2,7,1}) << endl;
  //  cout << Cross({2,7,1},{1,2,0}) << endl;  
  

  // cout << "\033[1;31mbold red text\033[0m\n";

  //  cout << Dot({{1,2,3},{1,2,3}},{1,2,3}) << endl;
  //  cout << Dot({{1,2,3},{1,2,3},{1,2,3}},{1,2,3}) << endl;    

  //  cout << Transpose({{1,2},{3,4}}) << endl;
  //  cout << -Dot({1,2},{1,3}) << endl;  
  //  const std::string red("\033[0;31m");
  //  const std::string green("\033[0;32m");
  //  const std::string yellow("\033[0;33m");
  //  cout << red << "test\n";
  //  cout << green << "test\n";
  //  cout << "\033[1;33mbold red text\033[0m\n";
  //  cout << "\033[1;34mbold red text\033[0m\n";
  //  cout << "\033[1;35mbold red text\033[0m\n";  

  // cout << "\033[1;31mbold red text\033[0m\n";

  // vector<double> v{5,2};
  // vector<vector<double>> mat0{{5,3},{3,pi}};
  // vector<vector<double>> mat1{{2,1/2.},{5,4}};  

  // cout << Dot(mat1,Dot(mat0, v)) << endl;
  // cout << Dot(Dot(mat1,mat0),v) << endl;


  // cout << Subdivide(0,10,10) << endl;


 
 
  // FILE *f = fopen( "test" , "w");

  // vector<vector<double>> data(100,vector<double>(100,0.));

  // for(size_t i=0; i<100; i++)
  //   for(size_t j=0; j<100; j++)
  //     {
  //       double x = pi/99 * i;
  //       double y = pi/99 * j;
  //       data[i][j] = sin(pow(y,2) + x);
  //     }
  // GNUPLOT mat;
  // mat.Set({{"palette","model RGB rgbformulae 1,3,5"},{"xrange","[0:100]"},{"yrange","[0:100]"}});
  // mat.SaveMatrixData(data,{{"w","image"}});
  // mat.MatrixPlot(); 

 
  // vector<double> vec1{cos(0.),sin(0.)};
  // vector<double> vec2{cos(300/180.*pi),sin(300./180.*pi)};

  // cout << VectorAngle(vec1,vec2)/pi*180. << endl;

 
  // vec1 +={1,1};
  // cout << vec1 << endl;
  // vec1 +={1,2};
  // cout << vec1 << endl;

  //  cout << setprecision(10) << -Dot({{1,2},{1,2},{1,2},{1,2}},{1,2})/pi << endl;  

  // vector<vector<double>> v1=SubdivideExclude({{5., 0.},{5.,-.5},{1.,-.5},{1.,0.}},34);
  // vector<vector<double>> v2=SubdivideExclude({{5., 0.},{5.,-.5},{1.,-.5},{1.,0.}},34);

  // cout << v1.size() << endl;
  // cout << v2.size() << endl;

  

 
   
  // std::vector<double> vec{1, 2}, r{1,2,3};
  // std::vector< std::vector<double> > v{{1, 2}, {2, 3}};

  // cout << Dot({{2,3},{2,3},{2,3}},vec) << endl;

  
  // vector<bool> b{false,true,false};
  // bool a = false;

  // for(size_t i = 0; i<b.size(); i++)
  //   {
  //     a +=b[i];
  //     cout << accumulate(b.begin(), b.end(), 0) << endl; 
  // 	//      if(a)
  // 	//cout << a << endl;
  //   }
    
  // vector<vector<double>> tmp=(SubdivideExclude({{0.,-1.},{2.,-1},{2.5,0.}},38));
  // //  PLOT.Plot2D(Transpose(tmp),{{"w","lp"}});

  // cout << Cross({1,2,3},{1,4,3}) << endl;
  // cout << Norm({1,2,3}) << endl;  
  // cout << Dot({{1,3,3},{1,2,3},{1,2,3},{1,2,1},{1,2,3}},{{1,2,3},{1,2,3},{1,2,3}}) << endl;  
  // cout << Transpose(SubdivideExclude({{0.,0.},{0.,-1.},{6.,-1.},{6.,0.}},30)) << endl;

  // std::vector<std::vector<double>> data;
  // double x, y, z;
  // int N=100;
  // for(int j=0; j<N; j++)      
  //   for(int i=0; i<N; i++)
  //     {
  // 	x = (double)i - 50.;
  // 	y = (double)j - 50.;      
  // 	z = cos((x+y)/(2.*pi));
  // 	std::vector<double> X{x,y,z};
  // 	data.push_back(X);
  //     }
  // cout << data << endl;
  // GNUPLOT PLOT2;
  // PLOT2.Set({{"hidden3d",""}});
  // PLOT2.Plot3D(data,{{"w","l"}});
  
  // std::vector<std::vector<std::vector<double>>> data;
  // for(int i=0; i<N; i++)
  //   for(int j=0; j<N; j++)
  //     for(int k=0; k<N; k++)      
  // 	{
  // 	  x = 3.*i/N - 1.5;
  // 	  y = 3.*j/N - 1.5;
  // 	  z = 3.*k/N - 1.5;	
  // 	  std::vector<double> X{x,y,z};
  // 	  std::vector<double> dX{y, x - pow(x,3), z - pow(z,3)};	
  // 	  data.push_back({ X, dX/5. , {Norm(dX)} });
  // 	}

  // GNUPLOT test;
  // test.VectorPlot3D(data ,{{"lc","palette"}});

  // vector<double> vec{1,2,3};
  // vector<double>* pointer_vector;

  // pointer_vector = &vec;
  // cout << vec << endl;
  // cout << (*pointer_vector).size() << endl;

  // vector<double> a={1,2};
  // vector<double>* b=&a;    
  // test Test(b);

  // cout << *Test.v << endl;

  // a={1,2};
  // Test.v=&a;
  // cout << *Test.v << endl;
  // a={1,2};
  // cout << *Test.v << endl;

  // Test.input(a);
  // cout << *Test.v << endl;


  // vector<double> e={1,2.2};
  // Test.input(e);
  // cout << *Test.v << endl;
  // e.push_back(1000);
  // cout << *Test.v << endl;
  
  // map<string, string> mp;

  //  mp["lc"] = "palette";
  //  if(mp.find("lc")==mp.end())
  //    cout <<  << endl;


  // if(mp["lc"]=="pallete")
  //   cout << "fzg" << endl;
  // sub.Plot2D(Transpose(SubdivideExclude({{0,0},{0,1},{2,1},{2,0}},30)),{{"ls","3"},{"w","lp"}});
  
  // PLOT.Plot3D({{1,3,3},{1,23,3},{2,3,4}},{{"w","i"}});
  // cout << tmp[0].size() << endl;
  // PLOT.Plot2D({-Dot({{1,2},{2,4},{1,3},{2,2},{1,2}},{1,2}),-Dot({{9,2},{3,1},{1,0},{2,2},{3,2}},{2,2})});
  // PLOT.Plot2D({-Dot({{1,2},{2,4},{1,3},{2,2},{1,2}},{1,2}),-Dot({{9,2},{3,1},{1,0},{2,2},{3,2}},{1,2})});
  // PLOT.Plot2D({-Dot({{1,2},{2,4},{1,3},{2,2},{1,2}},{1,2}),-Dot({{9,2},{3,1},{1,0},{2,2},{3,2}},{2,3})},{{"ls","50"},{"w","l"}});
  // PLOT.Plot2D();

  // PLOT.Data.clear();
  // PLOT.SaveData(-Dot({{1,2},{1,2},{1,2},{1,2}},{1,2}));
  // PLOT.Plot1D(vector<int>{0});
  
  return 0;
}
