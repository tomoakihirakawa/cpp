#include "GNUPLOT.hpp"
#include "fundamental.hpp"

struct interp_method{
  virtual ~interp_method(){};
  virtual std::vector<double> operator()(const double& t0,
					 const double& t1) const = 0;  
};

class interpolation{
private:
  interp_method *interp; 
public:
  interpolation(interp_method *interpIN) : interp(interpIN){};  
  ~interpolation(){delete this->interp;}
  void set_interp(interp_method *interpIN){
    delete this->interp;
    this->interp = interpIN;
  };

  std::vector<double> operator()(const double& t0, const double& t1){
    return (*(this->interp))(t0,t1);
  };
  
};


struct interp_linear : public interp_method{
  std::vector<std::vector<double>> sample;  
  interp_linear(const std::vector<std::vector<double>>& sample_IN):sample(sample_IN){};  
  std::vector<double> operator()(const double& t0, const double& t1) const override{
    return Dot({t0, t1, 1.-(t0+t1)},this->sample);
  };
};

struct interp_2nd : public interp_method{
  std::vector<std::vector<double>> sample;  
  interp_2nd(const std::vector<std::vector<double>>& sample_IN):sample(sample_IN){};  
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
  
int main(){

  GNUPLOT plot,plotA,plotB,plotC;

  //  interpolation interpLoc(new interp_linear({{1,1,1},{0,1,1},{0,0,-1}}));
  interpolation interpLoc(new interp_2nd({{1,1,1},{0,1,1},{0,0,-1},
					  {.5,.3,1},{0,.5,1},{0.5,0.2,0}}));  
  interpolation interpRGB(new interp_linear({{255.,0,0},{0,255.,0},{0,0,255.}}));
  interpolation interpR(new interp_linear({{255.,0,0},{0,0,0},{0,0,0}}));
  interpolation interpG(new interp_linear({{0,0,0},{0,255.,0},{0,0,0}}));
  interpolation interpB(new interp_linear({{0,0,0},{0,0,0},{0,0,255.}}));  
  
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
