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

  GNUPLOT plot;
  int s = 30;

  std::vector<std::vector<double>> p = {{0,2,1},{-2,-1,-.5},{2,-1,0},
					{-1.,0,.1},{0,-1.5,0},{1,0,-.5}};

  std::vector<std::vector<double>> c = {{255,0,0},{0,255,0},{0,0,255},
					{255,255,0},{0,255,255},{255,0,255}};
  
  {
    interpolate3d interpLoc({p[3],p[4],p[5]});
    interpolate3d interpRGB({c[3],c[4],c[5]});
    for(auto i=0; i<s; i++)
      for(auto j=0; j<s-i; j++)
	plot.SaveData(  interpLoc(i/(s-1.),j/(s-1.)), {{"ps","1"},{"pt","7"},{"w","p"},{"lc", plot.rgb( interpRGB(i/(s-1.),j/(s-1.))) }});
  }
  {
    interpolate3d interpLoc({p[0],p[3],p[5]});
    interpolate3d interpRGB({c[0],c[3],c[5]});
    for(auto i=0; i<s; i++)
      for(auto j=0; j<s-i; j++)
	plot.SaveData(  interpLoc(i/(s-1.),j/(s-1.)), {{"ps","1"},{"pt","7"},{"w","p"},{"lc", plot.rgb( interpRGB(i/(s-1.),j/(s-1.))) }});
  }
  {
    interpolate3d interpLoc({p[3],p[1],p[4]});
    interpolate3d interpRGB({c[3],c[1],c[4]});
    for(auto i=0; i<s; i++)
      for(auto j=0; j<s-i; j++)
	plot.SaveData(  interpLoc(i/(s-1.),j/(s-1.)), {{"ps","1"},{"pt","7"},{"w","p"},{"lc", plot.rgb( interpRGB(i/(s-1.),j/(s-1.))) }});
  }
  {
    interpolate3d interpLoc({p[5],p[4],p[2]});
    interpolate3d interpRGB({c[5],c[4],c[2]});
    for(auto i=0; i<s; i++)
      for(auto j=0; j<s-i; j++)
	plot.SaveData(  interpLoc(i/(s-1.),j/(s-1.)), {{"ps","1"},{"pt","7"},{"w","p"},{"lc", plot.rgb( interpRGB(i/(s-1.),j/(s-1.))) }});
  }
  
  plot.Plot3D();
  std::cin.ignore();    
  
};
