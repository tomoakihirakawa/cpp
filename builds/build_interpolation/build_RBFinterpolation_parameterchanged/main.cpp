
#include "fundamental.hpp"
#include "InterpolationRBF.hpp"
#include <iostream>
#include <string>
#include <cmath>

#include "GNUPLOT.hpp"

#include <algorithm>

#include "Network.hpp"

double RBFscale(const VV_d& sample){
  V_d r;
  for(auto i=0; i<sample.size(); i++)
    for(auto j=i+1; j<sample.size(); j++)
      r.push_back(Norm(sample[i]-sample[j]));
  
  std::sort(r.begin(), r.end(), [](double lhs, double rhs) {return lhs < rhs;});

  auto s = 1+(int)((double)r.size()/3.);
  V_d v(s);
  for(auto i=0; i<s; i++)
    v[i] = r[i];
  
  return Mean(v);
  
};

double RBFscale(const VV_d& sample, const V_d& x){
  V_d r;
  for(auto i=0; i<sample.size(); i++)
    r.push_back( Norm( sample[i]-x ) );
  
  std::sort(r.begin(), r.end(), [](const auto& lhs, const auto& rhs) {return lhs < rhs;});

  auto s = 1+(int)((double)r.size()/3.);
  V_d v(s);
  for(auto i=0; i<s; i++)
    v[i] = r[i];
  
  return Mean(v);  
};

double func(double x, double y){
  double r = sqrt(x*x + y*y);
  return sin(r);  
};

double func2(double x, double y, double z){
  double r = sqrt(x*x + y*y);
  return sin(r);  
};

V_d dfunc(double x, double y){
  double r = sqrt(x*x + y*y);
  return {cos(r)*x/r,cos(r)*y/r,0};
};

using V_d = std::vector<double>;

int main(){
  
  NetworkX net({1,1}, {2,2,2}, 2, "no_name");
  
  VV_d xyz1({}), xyz2({});
  VV_d ql({});
  
  int ind = 10;
  double level=0;
  V_d axis={1.,0.,0.};
  
  auto o = net.Points[ind];
  auto normal = o->getFaces()[0]->getNormal();

  depth_searcher s1(3), s2(4);
  s1.set(o);
  s1.search();
  s2.set(o);
  s2.search();

  std::map<netPp,V_d> mp;
  for(auto& p:s2.getObjects()){
    xyz2.push_back(p->getX());
    mp[p]=p->getX();
  }
  for(auto& p:s1.getObjects()){
    xyz1.push_back(p->getX());
    mp[p]=p->getX();
  }

  for(auto& p:s1.getObjects()){
    xyz1.push_back(p->getX());
    mp[p]=p->getX();
  }
  
  InterpolationRBF interp(Transpose(VV_d{vx,vy}), vz, f, dfd);
  
  Print(xyz1,Green); 
  Print(xyz2,Green); 
  //サンプルデータを基に補間    
  GNUPLOT plt;
  plt.Set({{"key",""},{"contour","base"},{"zrange","[-1.5:1.5]"}});
  plt.SaveData(xyz1,{{"pt","7"},{"w","p"},{"lc","'blue'"},{"title","s1"}});
  plt.SaveData(xyz2,{{"pt","7"},{"w","p"},{"lc","'blue'"},{"title","s2"}});
  //plt.SaveVectorData(getVectorData(net.Faces),{{"lc","'magenta'"}});    
  plt.plot3d();
  std::cin.ignore();

}
