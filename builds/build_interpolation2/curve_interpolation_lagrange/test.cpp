#include "GNUPLOT.hpp"
#include "fundamental.hpp"

struct LagrangeInterpolation{
  std::vector<double> parameter;
  std::vector<double> sample;
  LagrangeInterpolation(const std::vector<double>& parameter_IN,
			const std::vector<double>& sample_IN):sample(sample_IN),parameter(parameter_IN){
  };
  double operator()(const double t){
    double ret=0, prod0=1, prod1=1;
    for(auto j=0; j<sample.size(); j++){
      prod0=1, prod1=1;    
      for(auto i=0; i<parameter.size(); i++){
	if(i!=j){
	  prod0 = prod0*(t - parameter[i]);
	  prod1 = prod1*(parameter[i] - parameter[j]);	  
	}
      }
      ret += prod0/prod1*sample[j];
    };
    return ret;
  };
};

int main(){
  GNUPLOT plot;

  std::vector<double> param = Subdivide(0.,1.,6);
  std::vector<double> x = {0,1,2,3,4,0,4.}, y = {0,2,2,4,3,4,1};  
  
  LagrangeInterpolation interp1(param,x), interp2(param,y);

  for(const auto& t:Subdivide(0.,1.,1000))
    plot.SaveData({interp1(t), interp2(t), 0.},{{"pt","1"},{"lc","1"},{"w","p"},{"lw","1"}});
    
  for(auto i=0; i<x.size(); i++)
    plot.SaveData({x[i], y[i], 0.},{{"pt","3"},{"lc","3"},{"w","p"},{"lw","4"}});   
  
  plot.Plot3D();
  std::cin.ignore();  
};
