
#include "GNUPLOT.hpp"
#include "fundamental.hpp"

using namespace std;

int main(){

  std::vector<std::vector<double>> a;
  std::vector<std::vector<double>> b;
  std::vector<std::vector<double>> c;

  Load(Directory(__FILE__)+"csv/a.csv",a);
  Load(Directory(__FILE__)+"csv/b.csv",b);
  Load(Directory(__FILE__)+"csv/c.csv",c);

  
  {
    ParametricInterpolation intpX(a, 3);
    ParametricInterpolation intpY(b, 3);
    ParametricInterpolation intpZ(c, 3);
    
    GNUPLOT plot;
    {
      int n=10;
      for(auto nn=0; nn<intpX.s2; nn++)
	for(auto mm=0; mm<intpX.s1; mm++){
	  std::vector<std::vector<std::vector<double>>> data;
	  for(auto i=-n; i<n+1; i++){
	    std::vector<std::vector<double>> vec;
	    for(auto j=-n; j<n+1; j++)	
	      {
		double h = (double)i/n;
		double xi = (double)j/n;
		vec.push_back({intpX({xi,h}),intpY({xi,h}),intpZ({xi,h})});
	      };
	    data.push_back(vec);
	  }
	  plot.SaveSplotData(data,{{"w","l"}});
	}
    }
      
    plot.Splot();
    cin.ignore();
  }
  
  return 0;

};

