#include "GNUPLOT.hpp"
#include "fundamental.hpp"


struct interpolate3d{
  std::vector<std::vector<double>> sample;  
  interpolate3d(const std::vector<std::vector<double>>& sample_IN):sample(sample_IN){};  
  std::vector<double> operator()(const double& t0, const double& t1){
    return Dot({t0, t1, 1.-(t0+t1)},this->sample);
  };
};
  
int main(){

  GNUPLOT plot,plotA,plotB,plotC;

  interpolate3d interpLoc({{1,1,1},{0,1,1},{0,0,-1}});
  interpolate3d interpRGB({{255.,0,0},{0,255.,0},{0,0,255.}});
  interpolate3d interpR({{255.,0,0},{0,0,0},{0,0,0}});
  interpolate3d interpG({{0,0,0},{0,255.,0},{0,0,0}});
  interpolate3d interpB({{0,0,0},{0,0,0},{0,0,255.}});  
  
  int s = 60;
  std::vector<double> t;    
  for(auto i=0; i<s; i++)
    for(auto j=0; j<s-i; j++){
      plot.SaveData(  interpLoc(i/(s-1.),j/(s-1.)), {{"pt","7"},{"w","p"},{"lc", plot.rgb( interpRGB(i/(s-1.),j/(s-1.))) }});
      plotA.SaveData( interpLoc(i/(s-1.),j/(s-1.)), {{"pt","7"},{"w","p"},{"lc", plot.rgb( interpR(i/(s-1.),j/(s-1.))) }});
      plotB.SaveData( interpLoc(i/(s-1.),j/(s-1.)), {{"pt","7"},{"w","p"},{"lc", plot.rgb( interpG(i/(s-1.),j/(s-1.))) }});
      plotC.SaveData( interpLoc(i/(s-1.),j/(s-1.)), {{"pt","7"},{"w","p"},{"lc", plot.rgb( interpB(i/(s-1.),j/(s-1.))) }}); 
    }
  
  plot.Plot3D();
  std::cin.ignore();    
  plotA.Plot3D();
  std::cin.ignore();    
  plotB.Plot3D();
  std::cin.ignore();    
  plotC.Plot3D(); 
  std::cin.ignore();  
  
};
