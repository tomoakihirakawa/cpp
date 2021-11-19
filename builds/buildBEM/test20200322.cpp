
#include "GNUPLOT.hpp"
#include "fundamental.hpp"


// ///////////////////////////////////////////////////////////
// struct StarInterpolation{
// public:
//   std::vector<double> sampC, sampCIN;
//   std::vector<std::vector<double>> samp3, Dsamp3;  
//   std::vector<double> samp2;
//   std::vector<double> q, xi, eta;  
//   std::vector<double> q1, q2;  
//   std::vector<std::vector<double>> InvMatB, InvMatBT;
//   int s, K;
//   int s1, s2;
//   ///////////// 1 parm  ///////////
//   StarInterpolation(const std::vector<double>& sampCIN,,
// 		    const std::vector<double>& DsampCIN,
// 		    const std::vector<std::vector<double>>& samp3_IN,
// 		    const std::vector<std::vector<double>>& Dsamp3_IN,
// 		    const int K_IN):
//     sampC(sampCIN),
//     DsampC(DsampCIN),    
//     samp3(samp3_IN),
//     Dsamp3(Dsamp3_IN),    
//     s1(samp3_IN.size())/*y*/,
//     s2(samp3_IN[0].size())/*x*/,
//     q1(OpenUniformKnots(samp3_IN.size(), K_IN)),
//     q2(OpenUniformKnots(samp3_IN[0].size(), K_IN)),
//     K(K_IN),
//     InvMatB(Inverse(Bspline_matrix(samp3_IN.size()/*y*/, K_IN))),
//     InvMatBT(Inverse(Transpose(Bspline_matrix(samp3_IN[0].size()/*x*/, K_IN)))){};
//   StarInterpolation(){};
//   ////////////////////////////////////////////////////  
//   std::vector<double> Nx(const double b/*y*/){
//     return Dot(InvMatBT,Bspline_vector(b, s2, K, q2));
//   };
//   std::vector<double> Ny(const double a/*x*/){
//     return Dot(Bspline_vector(a, s1, K, q1),InvMatB);    
//   };
//   std::vector<std::vector<double>> N(const std::vector<double>& ba){
//     return TensorProduct(Ny(ba[1]), Nx(ba[0]));
//   };
//   double operator()(const double b/*x*/, const double a/*y*/){
//     return Dot(Ny(a),Dot(samp3,Nx(b)));
//   };
//   double operator()(const std::vector<double>& ba){
//     return Dot(Ny(ba[1]),Dot(samp3,Nx(ba[0])));
//   };
//   ////////////////// DERIVATIVES /////////////////
//   std::vector<double> DNx(const double b/*y*/, const int n){
//     return Dot(InvMatBT,D_Bspline_vector(b, s2, K, n, q2));
//   };
//   std::vector<double> DNy(const double a/*x*/, const int n){
//     return Dot(D_Bspline_vector(a, s1, K, n, q1),InvMatB);
//   };
//   std::vector<std::vector<double>> DN(const std::vector<double>& ba, const std::vector<int>& n){
//     return TensorProduct(DNy(ba[1], n[1]), DNx(ba[0], n[0]));
//   };
//   double operator()(const double b/*x*/, const double a/*y*/, const std::vector<int>& n){
//     return Dot(DNy(a, n[1]),Dot(samp3,DNx(b, n[0])));
//   };
//   double operator()(const std::vector<double>& ba, const std::vector<int>& n){
//     return Dot(DNy(ba[1], n[1]),Dot(samp3,DNx(ba[0], n[0])));
//   };
//   ///////////////// specific shape function ////////////
//   double Nx(const double b/*y*/, const int j){
//     return Dot(InvMatBT[j],Bspline_vector(b, s2, K, q2));
//   };  
//   double Ny(const double a/*x*/, const int m){
//     std::vector<double> By = Bspline_vector(a, s1, K, q1);
//     double ret(0.);
//     for(auto l=0; l < s1; l++)        
//       ret += InvMatB[l][m] * By[l];
//     return ret;
//   };
//   double N(const std::vector<double>& ba, const std::vector<int> &jm){/*2D vector access*/
//     return Ny(ba[1], jm[1]) * Nx(ba[0], jm[0]);
//   };
//   double N(const std::vector<double>& ba, const int mj){/*1D access*/
//     return Ny(ba[1], (int)(mj/s2)/*row*/) * Nx(ba[0], mj%s2/*col*/);
//   };
//   // double Basis(const std::vector<double>& ba, const std::vector<int> &jm){
//   //   return samp3[jm[1]][jm[0]] * Nx(ba[0], jm[0]) * Ny(ba[1], jm[1]);
//   // };
//   double Basis(const std::vector<double>& ba, const std::vector<int> &jm){
//     return Bspline(ba[0], q2, jm[0], K) * Bspline(ba[1], q1, jm[1], K);
//   };  
//   double Shape(const std::vector<double>& ba, const std::vector<int> &jm){
//     return Nx(ba[0], jm[0]) * Ny(ba[1], jm[1]);
//   };  
// };


using namespace std;

int main(){

  vector<vector<double>> around={{-4.55816, 11.4788, 7.74909},
				 {-1.13035, 11.5006, 8.13536},
				 {-0.590454, 14.1412, 6.59582},
				 {-4.22943, 15.6584, 4.11027},
				 {-5.98386, 15.2974, 3.709},
				 {1.1849, -18.6772, -0.388354}};
  vector<vector<double>> Daround(around.size(),vector<double>(around[0].size(),0));

  vector<vector<double>> center(around.size(), {-5.70826, 13.6563, 6.26488});
  
  vector<vector<double>> Dcenter(around.size(), {0,0,0});  

  cout << center << endl;
  
  int K=2;

  std::vector<double> eta=Subdivide(-1,1,around.size()-1);
  std::vector<double> xi=Subdivide(-1,1,2-1);

  std::vector<std::vector<double>> s=Transpose({Transpose(center)[0],Transpose(around)[0]}); 

  //////////////////////////////////
  cout << xi << endl;
  vector<double> q_xi=OpenUniformKnots(xi,K);
  cout << q_xi << endl;  
  cout << Bspline(xi,q_xi,K) << endl;
  /////////////////////////////////
  cout << eta << endl;
  vector<double> q_eta=OpenUniformKnots(eta,K);
  cout << q_eta << endl;  
  cout << Bspline(eta,q_eta,K) << endl;
  ////////////////////////////////
  
  GNUPLOT mat;
  mat.SaveMatrixData(s,{{"w","image"}});
  mat.MatrixPlot(); 

  std::cin.ignore();
//  cout << TensorProductSet(eta,xi) << endl;
  
  // cout << Transpose(around) << endl;
  // cout << = Transpose(around)[0] << endl;  
  // cout << = Transpose(around)[1] << endl;
  // cout << = Transpose(around)[2] << endl;    

  
  // std::vector<int> vecs{1,2,3,4,1,9};

  // cout << Position(vecs,3) << endl;


  // std::vector<std::vector<int>> mat{{1, 3}, {9, 10}, {10, 3}, {414, 1}, {414, 9}};

  // cout << SortVectorChain(mat) << endl;
  // cout << Union(SortVectorChain(mat)) << endl;
  
  // std::vector<std::vector<double>> a;
  // std::vector<std::vector<double>> b;
  // std::vector<std::vector<double>> c;

  // Load(Directory(__FILE__)+"csv/a.csv",a);
  // Load(Directory(__FILE__)+"csv/b.csv",b);
  // Load(Directory(__FILE__)+"csv/c.csv",c);

  
  // {
  //   StarInterpolation intpX(a, 3);
  //   StarInterpolation intpY(b, 3);
  //   StarInterpolation intpZ(c, 3);
    
  //   GNUPLOT plot;
  //   {
  //     int n=10;
  //     for(auto nn=0; nn<intpX.s2; nn++)
  // 	for(auto mm=0; mm<intpX.s1; mm++){
  // 	  std::vector<std::vector<std::vector<double>>> data;
  // 	  for(auto i=-n; i<n+1; i++){
  // 	    std::vector<std::vector<double>> vec;
  // 	    for(auto j=-n; j<n+1; j++)	
  // 	      {
  // 		double h = (double)i/n;
  // 		double xi = (double)j/n;
  // 		vec.push_back({intpX({xi,h}),intpY({xi,h}),intpZ({xi,h})});
  // 	      };
  // 	    data.push_back(vec);
  // 	  }
  // 	  plot.SaveSplotData(data,{{"w","l"}});
  // 	}
  //   }
      
  //   plot.Splot();
  //   cin.ignore();
  // }
  
  return 0;

};

