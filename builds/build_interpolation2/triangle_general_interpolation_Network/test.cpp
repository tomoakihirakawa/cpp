#include "GNUPLOT.hpp"
#include "fundamental.hpp"

class interpolation{
public:
  std::vector<std::vector<double>> sample;
  interpolation(const std::vector<std::vector<double>>& sample_IN):sample(sample_IN){}; 
  virtual std::vector<double> operator()(const double&, const double&) const = 0;  
};


struct interp_linear : public interpolation{
  interp_linear(const std::vector<std::vector<double>>& sample_IN):interpolation(sample_IN){};
  std::vector<double> operator()(const double& t0, const double& t1) const override{
    return Dot({t0, t1, 1.-(t0+t1)},this->sample);
  };
};

struct interp_2nd : public interpolation{
  interp_2nd(const std::vector<std::vector<double>>& sample_IN):interpolation(sample_IN){};
  std::vector<double> operator()(const double& t0, const double& t1) const override{
    double t2 = 1-t0-t1;
    return Dot({t0*(2.*t0 - 1.)/*0*/,
		t1*(2.*t1 - 1.)/*1*/,
		t2*(2.*t2 - 1.)/*2*/,
		4.*t0*t1/*3*/,
		4.*t1*t2/*4*/,
		4.*t0*t2/*5*/},sample);            
  };
};  


#include "BEM_Network.hpp"


// ＊形状の空間微分や補間を瞬時に呼び出せるユーティリティ
// ＊変数の微分や補間を呼び出せるユーティリティ
// 積分
// 行列

// 時間発展を処理していくクラス？

int main(){
  NetworkObj obj("./obj/tank.obj");
  obj.rotate(2*M_PI, {1.,0.,0.});
  BEM_NetworkW water({1,1},{1.1,1.15,1/2.},.3);
  
  GNUPLOT plot,plotA,plotB,plotC;

  for(const auto& p:water.Points){
    Print(p->xyz);
    Print(p->getNormal());    
    plot.SaveVectorData({{p->xyz,p->getNormal()}});
  }

  // //  interpolation interpLoc(new interp_linear({{1,1,1},{0,1,1},{0,0,-1}}));
  // interp_2nd interpLoc({{1,1,1},{0,1,1},{0,0,-1},
  // 			{.5,.3,1},{0,.5,1},{0.5,0.2,0}});  
  // interp_linear interpRGB({{255.,0,0},{0,255.,0},{0,0,255.}});
  // interp_linear interpR({{255.,0,0},{0,0,0},{0,0,0}});
  // interp_linear interpG({{0,0,0},{0,255.,0},{0,0,0}});
  // interp_linear interpB({{0,0,0},{0,0,0},{0,0,255.}});  
  
  // int s = 60;
  // std::vector<double> t;    
  // for(auto i=0; i<s; i++)
  //   for(auto j=0; j<s-i; j++){
  //     plot.SaveData(  interpLoc(i/(s-1.),j/(s-1.)), {{"pt","7"},{"w","p"},{"lc", plot.rgb( interpRGB(i/(s-1.),j/(s-1.))) }});
  //     plotA.SaveData( interpLoc(i/(s-1.),j/(s-1.)), {{"pt","7"},{"w","p"},{"lc", plot.rgb( interpR(i/(s-1.),j/(s-1.))) }});
  //     plotB.SaveData( interpLoc(i/(s-1.),j/(s-1.)), {{"pt","7"},{"w","p"},{"lc", plot.rgb( interpG(i/(s-1.),j/(s-1.))) }});
  //     plotC.SaveData( interpLoc(i/(s-1.),j/(s-1.)), {{"pt","7"},{"w","p"},{"lc", plot.rgb( interpB(i/(s-1.),j/(s-1.))) }}); 
  //   }
  
  // plot.Plot3D();
  // std::cin.ignore();    
  // plotA.Plot3D();
  // std::cin.ignore();    
  // plotB.Plot3D();
  // std::cin.ignore();    

  plot.plot3d(); 
  std::cin.ignore();  
  
};
