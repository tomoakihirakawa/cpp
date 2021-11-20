//==========================================================
// test.cpp
//==========================================================
// #include <iostream>
// #include <fstream>
// #include <string>
// #include <vector>
// #include <map>
// #include <algorithm>
// #include <utility>//swap
// #include <cmath>
// #include <limits>
// #include <iomanip>
// #include <numeric>
// #include <chrono>
// #include <typeinfo>
// #include <sys/stat.h>
// #include <functional> 

#include "fundamental.hpp"
#include "GNUPLOT.hpp"

struct ParametricInterpolationInsert{
public:
  int s;  
  std::vector<double> samp2;
  std::vector<double> q, xi, eta;
  //
  std::vector<std::vector<double>> samp3;
  std::vector<std::vector<double>> mat_xi, mat_eta;;  
  std::vector<double> q1, q2;
  std::vector<std::vector<double>> inv_coefmat, inv_coefmatT;
  int K, s1, s2;
  ///////////// 1 parm  ///////////
  ParametricInterpolationInsert(const std::vector<double>& samp2_IN, const int K_IN):
    samp2(samp2_IN),
    s(samp2_IN.size()),
    xi(Subdivide(-1,1,samp2_IN.size()-1)),
    q(OpenUniform(Subdivide(-1,1,samp2_IN.size()-1), K_IN)),
    K(K_IN),
    inv_coefmat(Inverse(base_spline_matrix(samp2_IN.size(), K_IN))),
    inv_coefmatT(Inverse(Transpose(base_spline_matrix(samp2_IN.size(), K_IN)))){};
  std::vector<double> N(const double a){return Dot(inv_coefmat, base_spline_vector(a, s, K, q));};
  double operator()(const double a){return Dot(N(a), samp2);};
  double basis(const double xi, const int n){return base_spline(xi, q, n, K);};
  std::vector<double> basis(const std::vector<double>& xi, const int n){return base_spline(xi, q, n, K);};
  std::vector<double> N(const double a, const int n){return Dot(inv_coefmat, D_base_spline_vector(a, s, K, n, q));};  
  double operator()(const double a, const int n){return Dot(N(a,n),D_base_spline_vector(a, s, K, n, q));};  
  void Insert(const int i){
    double u = (*(xi.begin()+i) + *(xi.begin()+i+1))/2.;
    xi.insert(xi.begin()+i+1,	u );
    samp2.insert( samp2.begin() + i + 1, (*this)(u) );
    q.insert(q.begin() + K + i, (*(q.begin()+K+i-1) + *(q.begin()+K+i))/2. );
    //
    s++;
    inv_coefmat = Inverse(base_spline_matrix(xi,s,K,q));
  };
  void Erase(const int i){
    xi.erase(xi.begin()+i+1);
    samp2.erase( samp2.begin() + i + 1);
    q.erase(q.begin() + K + i);
    //
    s--;
    inv_coefmat = Inverse(base_spline_matrix(xi,s,K,q));
  };
  /////////// 2 parm  ////////////
  /*
    s1 is the size of a column vector : points in y direction
    s2 is the size of a row vector : points in x direction
  */
  ParametricInterpolationInsert(const std::vector<std::vector<double>>& samp3_IN, const int K_IN):
    samp3(samp3_IN),
    s1(samp3_IN.size())/*y*/,
    s2(samp3_IN[0].size())/*x*/,
    q1(OpenUniform(samp3_IN.size(), K_IN)),
    q2(OpenUniform(samp3_IN[0].size(), K_IN)),
    K(K_IN),
    inv_coefmat(Inverse(base_spline_matrix(samp3_IN.size()/*y*/, K_IN))),
    inv_coefmatT(Inverse(Transpose(base_spline_matrix(samp3_IN[0].size()/*x*/, K_IN)))),
    mat_eta(samp3_IN.size(),std::vector<double>(samp3_IN.size(),0)),
    mat_xi(samp3_IN[0].size(),std::vector<double>(samp3_IN[0].size(),0)){

    eta=(Subdivide(-1,1,s1-1));
    xi=(Subdivide(-1,1,s2-1));
    
    for(auto row=0; row<s1; row++)    
      mat_eta[row]=(eta);
    mat_eta=Transpose(mat_eta);
    for(auto col=0; col<s2; col++)
      mat_xi[col]=(xi);
  };
  ParametricInterpolationInsert(){};
  void reset(const std::vector<std::vector<double>>& samp3_IN, const int K_IN){
    q1.clear();
    q2.clear();
    samp3.clear();
    inv_coefmat.clear();
    inv_coefmatT.clear();
    //insert
    samp3 = samp3_IN;
    s1 = samp3_IN.size()/*y*/;
    s2 = samp3_IN[0].size()/*x*/;
    q1 = OpenUniform(samp3_IN.size(), K_IN);
    q2 = OpenUniform(samp3_IN[0].size(), K_IN);
    K = K_IN;
    inv_coefmat = Inverse(base_spline_matrix(samp3_IN.size(), K_IN));
    inv_coefmatT = Inverse(Transpose(base_spline_matrix(samp3_IN[0].size(), K_IN)));
  };
  ////////////////////////////////////////////////////  
  void Insert3(const int i){
    double u = (*(xi.begin()+i) + *(xi.begin()+i+1))/2.;
    
    std::vector<double> smp;
    for(auto l = 0; l<mat_eta.size(); l++)
      smp.push_back((*this)(u,mat_eta[l][0]));

    mat_xi[0].insert(mat_xi[0].begin()+i+1, u );
    std::vector<double> tmp = mat_xi[0];
    mat_xi.clear();
    for(auto l=0; l<s2+1; l++)
      mat_xi.push_back(tmp);
      
    for(auto m = 0; m<samp3.size(); m++)
      samp3[m].insert( samp3[m].begin() + i + 1, smp[m] );
    
    q2.insert(q2.begin() + K + i, (*(q2.begin()+K+i-1) + *(q2.begin()+K+i))/2. );
    //
    s2++;

    // std::vector<std::vector<double>> h;
    // for(auto row=0; row<s2; row++)
    //   for(auto col=0; col<s1; col++)	
    // 	h[row][col]={};
      
    inv_coefmatT = Inverse(Transpose(base_spline_matrix(mat_xi,s2,K,q2)));
    inv_coefmat  = Inverse(base_spline_matrix(Transpose(mat_eta),s1,K,q1));
    
    //    inv_coefmatT = Inverse(base_spline_matrix(xi,s1,K,q1));
    //    std::cout << inv_coefmatT << std::endl;
    //    std::cout << q1 << std::endl;
    
  };
  ////////////////////////////////////////
  std::vector<double> Nx(const double b/*y*/){return Dot(inv_coefmatT,base_spline_vector(b, s2, K, q2)); };
  std::vector<double> Ny(const double a/*x*/){return Dot(base_spline_vector(a, s1, K, q1),inv_coefmat);};
  std::vector<std::vector<double>> N(const std::vector<double>& ba){return TensorProduct(Ny(ba[1]), Nx(ba[0]));};
  double operator()(const double b/*x*/, const double a/*y*/){return Dot(Ny(a),Dot(samp3,Nx(b)));};
  double operator()(const std::vector<double>& ba){return Dot(Ny(ba[1]),Dot(samp3,Nx(ba[0])));};
  ////////////////// DERIVATIVES /////////////////
  std::vector<double> DNx(const double b/*y*/, const int n){return Dot(inv_coefmatT,D_base_spline_vector(b, s2, K, n, q2));};
  std::vector<double> DNy(const double a/*x*/, const int n){return Dot(D_base_spline_vector(a, s1, K, n, q1),inv_coefmat);};
  std::vector<std::vector<double>> DN(const std::vector<double>& ba, const std::vector<int>& n){return TensorProduct(DNy(ba[1], n[1]), DNx(ba[0], n[0]));};
  double operator()(const double b/*x*/, const double a/*y*/, const std::vector<int>& n){return Dot(DNy(a, n[1]),Dot(samp3,DNx(b, n[0])));};
  double operator()(const std::vector<double>& ba, const std::vector<int>& n){return Dot(DNy(ba[1], n[1]),Dot(samp3,DNx(ba[0], n[0])));};
  ///////////////// specific shape function ////////////
  double Nx(const double b/*y*/, const int j){return Dot(inv_coefmatT[j],base_spline_vector(b, s2, K, q2));};  
  double Ny(const double a/*x*/, const int m){
    std::vector<double> By = base_spline_vector(a, s1, K, q1);
    double ret(0.);
    for(auto l=0; l < s1; l++)        
      ret += inv_coefmat[l][m] * By[l];
    return ret;
  };
  double N(const std::vector<double>& ba, const std::vector<int> &jm){return Ny(ba[1], jm[1]) * Nx(ba[0], jm[0]);};
  double N(const std::vector<double>& ba, const int mj){return Ny(ba[1], (int)(mj/s2)/*row*/) * Nx(ba[0], mj%s2/*col*/);};
  double basis(const std::vector<double>& ba, const std::vector<int> &jm){return samp3[jm[1]][jm[0]] * Ny(ba[1], jm[1]) * Nx(ba[0], jm[0]);};
};

struct ParametricInterpolation{
public:
  std::vector<std::vector<double>> samp3;
  std::vector<double> samp2;
  std::vector<double> q1, q2;  
  std::vector<std::vector<double>> inv_coefmat, inv_coefmatT;
  int s, K;
  int s1, s2;
  ///////////// 1 parm  ///////////
  ParametricInterpolation(const std::vector<double>& samp2_IN, const int K_IN):
    samp2(samp2_IN),
    s(samp2_IN.size()),
    K(K_IN),
    inv_coefmat(Inverse(base_spline_matrix(samp2_IN.size(), K_IN))),
    inv_coefmatT(Inverse(Transpose(base_spline_matrix(samp2_IN.size(), K_IN)))){};
  std::vector<double> N(const double a){return Dot(base_spline_vector(a, s, K),inv_coefmat);};
  double operator()(const double a){return Dot(N(a), samp2);};
  //===========================
  std::vector<double> N(const double a, const int n){return Dot(D_base_spline_vector(a, s, K, n),inv_coefmat);};  
  double operator()(const double a, const int n){return Dot(N(a, n), samp2);};
  /////////// 2 parm  ////////////
  /*
    s1 is the size of a column vector : points in y direction
    s2 is the size of a row vector : points in x direction
  */
  ParametricInterpolation(const std::vector<std::vector<double>>& samp3_IN, const int K_IN):
    samp3(samp3_IN),
    s1(samp3_IN.size())/*y*/,
    s2(samp3_IN[0].size())/*x*/,
    q1(OpenUniform(samp3_IN.size(), K_IN)),
    q2(OpenUniform(samp3_IN[0].size(), K_IN)),
    K(K_IN),
    inv_coefmat(Inverse(base_spline_matrix(samp3_IN.size()/*y*/, K_IN))),
    inv_coefmatT(Inverse(Transpose(base_spline_matrix(samp3_IN[0].size()/*x*/, K_IN)))){};
  ParametricInterpolation(){};
  void reset(const std::vector<std::vector<double>>& samp3_IN, const int K_IN){
    q1.clear();
    q2.clear();
    samp3.clear();
    inv_coefmat.clear();
    inv_coefmatT.clear();
    //insert
    samp3 = samp3_IN;
    s1 = samp3_IN.size()/*y*/;
    s2 = samp3_IN[0].size()/*x*/;
    q1 = OpenUniform(samp3_IN.size(), K_IN);
    q2 = OpenUniform(samp3_IN[0].size(), K_IN);
    K = K_IN;
    inv_coefmat = Inverse(base_spline_matrix(samp3_IN.size(), K_IN));
    inv_coefmatT = Inverse(Transpose(base_spline_matrix(samp3_IN[0].size(), K_IN)));
  };
  ////////////////////////////////////////////////////  
  std::vector<double> Nx(const double b/*y*/){
    return Dot(inv_coefmatT,base_spline_vector(b, s2, K, q2));
  };
  std::vector<double> Ny(const double a/*x*/){
    return Dot(base_spline_vector(a, s1, K, q1),inv_coefmat);    
  };
  std::vector<std::vector<double>> N(const std::vector<double>& ba){
    return TensorProduct(Ny(ba[1]), Nx(ba[0]));
  };
  double operator()(const double b/*x*/, const double a/*y*/){
    return Dot(Ny(a),Dot(samp3,Nx(b)));
  };
  double operator()(const std::vector<double>& ba){
    return Dot(Ny(ba[1]),Dot(samp3,Nx(ba[0])));
  };
  ////////////////// DERIVATIVES /////////////////
  std::vector<double> DNx(const double b/*y*/, const int n){
    return Dot(inv_coefmatT,D_base_spline_vector(b, s2, K, n, q2));
  };
  std::vector<double> DNy(const double a/*x*/, const int n){
    return Dot(D_base_spline_vector(a, s1, K, n, q1),inv_coefmat);
  };
  std::vector<std::vector<double>> DN(const std::vector<double>& ba, const std::vector<int>& n){
    return TensorProduct(DNy(ba[1], n[1]), DNx(ba[0], n[0]));
  };
  double operator()(const double b/*x*/, const double a/*y*/, const std::vector<int>& n){
    return Dot(DNy(a, n[1]),Dot(samp3,DNx(b, n[0])));
  };
  double operator()(const std::vector<double>& ba, const std::vector<int>& n){
    return Dot(DNy(ba[1], n[1]),Dot(samp3,DNx(ba[0], n[0])));
  };
  ///////////////// specific shape function ////////////
  double Nx(const double b/*y*/, const int j){
    return Dot(inv_coefmatT[j],base_spline_vector(b, s2, K, q2));
  };  
  double Ny(const double a/*x*/, const int m){
    std::vector<double> By = base_spline_vector(a, s1, K, q1);
    double ret(0.);
    for(auto l=0; l < s1; l++)        
      ret += inv_coefmat[l][m] * By[l];
    return ret;
  };
  double N(const std::vector<double>& ba, const std::vector<int> &jm){/*2D vector access*/
    return Ny(ba[1], jm[1]) * Nx(ba[0], jm[0]);
  };
  double N(const std::vector<double>& ba, const int mj){/*1D access*/
    return Ny(ba[1], (int)(mj/s2)/*row*/) * Nx(ba[0], mj%s2/*col*/);
  };
  double basis(const std::vector<double>& ba, const std::vector<int> &jm){
    return samp3[jm[1]][jm[0]] * Ny(ba[1], jm[1]) * Nx(ba[0], jm[0]);
  };
};

using namespace std;

int main(){
  
  if(true){
  
    int n=50;
    std::vector<double> vec(n,0);
    int i(0);
    for(auto& v:vec)
      {
	vec[i] = (sin((double)i/49.*M_PI*2.));
	i++;
      }
  
    ParametricInterpolationInsert intp(vec, 3);
  
    GNUPLOT plot;
    // std::vector<double> h(Subdivide(-1,1,1000-1));
    // // intp.Insert(3);
    // // intp.Insert(3);
    // // intp.Insert(3);
    // // intp.Insert(3);  
    // for(auto l=0; l<intp.s+3; l++)
    //   plot.SaveData(Transpose({h,intp.basis(h,l)}),{{"w","l"}});



    {
      n=25;
      std::vector<std::vector<double>> data;
      for(auto i=-n; i<n+1; i++)
	{
	  double h = (double)i/n;
	  data.push_back({h,intp(h)});
	};
      plot.SaveData(data);
    }

  
    plot.Plot2D();
    cin.ignore();
  
    {  
      // std::vector<std::vector<double>> vec
      //   {{6.,7.,4,7,3,0.,2.,1.,1.,1.},
      //    {1.,3.,2,7,8,1.,3.,5.,2.,1.},
      //    {1.,1.,2,1,1,0.,3.,3.,2.,3.},     
      //    {5.,3.,1,7,5,10.,6.,3.,5.,2.},
      //    {1.,7.,4,7,1,5.,5.,1.,5.,5.}};

      std::vector<std::vector<double>> vec(5,std::vector<double>(5,1.));
      // std::vector<double> vec2{1.,2.,3.,4.,5.,6.,7.};  
      ParametricInterpolationInsert intpy(vec, 3);
      //ParametricInterpolationInsert intpx(vec2, 2);  
      int n = 15;
      GNUPLOT plot;
      {
	std::vector<std::vector<std::vector<double>>> mat;
	for(auto j=-n; j<n+1; j++)
	  {
	    std::vector<std::vector<double>> data;
	    for(auto i=-n; i<n+1; i++)
	      {
		double xi = (double)i/n;
		double eta = (double)j/n;	  
		data.push_back({xi,eta,intpy(xi,eta)});
	      };
	    mat.push_back(data);
	  }
	//  plot.SaveSplotData(mat,{{"w","p"}});
	//  plot.SaveSplotData(mat,{{"w","l"}});    
      }

    
      // intpy.Insert3(1);
      // intpy.Insert3(3);
      // intpy.Erase(2);
      // intpy.Erase(3);
      // intpy.Erase(5);  
      // intpy.Erase(3);

      plot.Set({{"hidden3d",""}});
  
      {
	// for(auto row=0; row<intpy.s2; row++)
	//   for(auto col=0; col<intpy.s1; col++)
	{
	  std::vector<std::vector<std::vector<double>>> mat;
	  for(auto j=-n; j<n+1; j++)
	    {
	      std::vector<std::vector<double>> data;
	      for(auto i=-n; i<n+1; i++)
		{
		  double xi = (double)i/n;
		  double eta = (double)j/n;	  
		  data.push_back({xi,eta,intpy({xi,eta})});
		}
	      mat.push_back(data);	
	    };
	  //plot.SaveSplotData(mat,{{"w","p"}});
	  plot.SaveSplotData(mat,{{"w","l"}});    
	}
      }

      //    plot.Splot();

  
      // cin.ignore();

      // data.clear();

      // n=10;
      // for(auto i=-n; i<n+1; i++)
      //   {
      //     double h = (double)i/n;
      //     data.push_back({h,intpy(h)});
      //     cout << intpy(h) << endl;
      //   };

      // plot.SaveData(data,{{"w","p"}});  
      // plot.SaveData(data,{{"w","l"}});  

      // plot.Plot2D();
      cin.ignore();
    }
  }
  return 0;
};
