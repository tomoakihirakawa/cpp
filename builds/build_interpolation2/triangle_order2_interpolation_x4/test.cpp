#include "GNUPLOT.hpp"
#include "fundamental.hpp"

struct interpolate3d_2nd{
  std::vector<std::vector<double>> sample;  
  interpolate3d_2nd(const std::vector<std::vector<double>>& sample_IN):sample(sample_IN){};  
  std::vector<double> operator()(const double& t0, const double& t1){
    double t2 = 1-t0-t1;
    //    0
    //  3   5
    // 1  4  2    
    return Dot({t0*(2.*t0 - 1.)/*0*/,
		t1*(2.*t1 - 1.)/*1*/,
		t2*(2.*t2 - 1.)/*2*/,
		4.*t0*t1/*3*/,
		4.*t1*t2/*4*/,
		4.*t0*t2/*5*/},sample);            
  };
};
  
int main(){

  GNUPLOT plot;
  int s = 50;

  std::vector<std::vector<double>> p = {{0,2,1},{-2,-1,-.5},{2,-1,0},
					{-1.,0,.1},{0,-1.5,0},{1,0,-.5}};

  double max = 227;
  std::vector<std::vector<double>> c = {{max,0,0},{0,max,0},{0,0,max},
					{max,max,0},{0,max,max},{max,0,max}};
  
  interpolate3d_2nd interpLoc(p);
  interpolate3d_2nd interpRGB(c);
  
  for(auto i=0; i<s; i++)
    for(auto j=0; j<s-i; j++)
      plot.SaveData(  interpLoc(i/(s-1.),j/(s-1.)), {{"ps","1"},{"pt","7"},{"w","p"},{"lc", plot.rgb( interpRGB(i/(s-1.),j/(s-1.))) }});
    
  plot.Plot3D();
  std::cin.ignore();    
  
};
