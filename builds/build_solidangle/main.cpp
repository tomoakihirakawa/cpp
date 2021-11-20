#include "fundamental.hpp"

using V_double = std::vector<double>;
using VV_double = std::vector<std::vector<double>>;
using VVV_double = std::vector<std::vector<std::vector<double>>>;

int main(){
  
  // V_d org = {1/2., 1/2., 0};
  // V_d p0 = {0, 0, 1};
  // V_d p1 = {0, .5, 1};
  // V_d p2 = {1, .1, 1};
  // V_d p3 = {1., 0., 1.};
  // VV_d ps = {p0, p1, p2, p3};
  
  // std::cout << SphericalVectorAngle(p0-org, p1-org, p2-org) << std::endl;
  // std::cout << SphericalVectorAngle(p3-org, p0-org, p1-org) << std::endl;
  // std::cout << SphericalVectorAngle(p2-org, p3-org, p0-org) << std::endl;
  // std::cout << SphericalVectorAngle(p1-org, p2-org, p3-org) << std::endl;  
  // std::cout << geometry::SolidAngle(org, {p0,p1,p2,p3}) << std::endl;


  V_d org = {0., 0., 0};

  for(auto i=-100; i<100 ;i++){
    V_d p0 = {-1+2.,-1., (double)i/100.};
    V_d p1 = { 1+2.,-1., (double)i/100.};
    V_d p2 = { 0.+2, 1., (double)i/100.};
    VV_d ps = {p0,p1,p2};
    std::cout << std::setprecision(20) << geometry::SolidAngle(org,ps) << std::endl;
  }
};
