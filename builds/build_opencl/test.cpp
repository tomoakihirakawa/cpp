#include "fundamental.hpp"

using V_double = std::vector<double>;
using VV_double = std::vector<std::vector<double>>;
using VVV_double = std::vector<std::vector<std::vector<double>>>;

int main(){

  VV_d tab = {{.51, 2, 1}, {.5, 4, 3}, {.5, 1., 1}};
  std::cout << std::setprecision(15) << geometry::SolidAngle({.5,.5,0.},tab[0],tab[1],tab[2])  << std::endl;;
  
  return 0;      
};
